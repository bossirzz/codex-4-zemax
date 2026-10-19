import os
import json
import math
import time
from pathlib import Path

import numpy as np
import zospy as zp
import zospy_zh

OUT = Path(__file__).with_name('zospy_test_results')
OUT.mkdir(exist_ok=True)
report_path = OUT / 'report.json'
report = json.loads(report_path.read_text(encoding='utf-8')) if report_path.exists() else {'zospy_version': zp.__version__, 'tests': []}
report.pop('fatal_error', None)
zos = None


def emit(stage, **data):
    print(json.dumps({'stage': stage, **data}, ensure_ascii=True, default=str), flush=True)


def save():
    (OUT / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding='utf-8')


def number(value):
    v = float(value)
    return v if math.isfinite(v) else str(v)


def snapshot(oss):
    rows = []
    for i in range(int(oss.LDE.NumberOfSurfaces)):
        emit('reading_surface', surface=i)
        row = oss.LDE.GetSurfaceAt(i)
        rows.append({'surface': i, 'comment': str(row.Comment), 'radius': number(row.Radius),
                     'thickness': number(row.Thickness), 'material': str(row.Material),
                     'semi_diameter': number(row.SemiDiameter), 'conic': number(row.Conic),
                     'is_stop': bool(row.IsStop), 'type': str(row.Type)})
    sd = oss.SystemData
    fields = [sd.Fields.GetField(i) for i in range(1, int(sd.Fields.NumberOfFields) + 1)]
    waves = [sd.Wavelengths.GetWavelength(i) for i in range(1, int(sd.Wavelengths.NumberOfWavelengths) + 1)]
    return {'file': str(oss.SystemFile), 'needs_save': bool(oss.NeedsSave), 'rows': rows,
            'fields': [{'x': float(f.X), 'y': float(f.Y), 'weight': float(f.Weight)} for f in fields],
            'wavelengths': [{'wavelength': float(w.Wavelength), 'weight': float(w.Weight)} for w in waves],
            'aperture_type': str(sd.Aperture.ApertureType), 'aperture_value': number(sd.Aperture.ApertureValue),
            'merit_operands': int(oss.MFE.NumberOfOperands)}


def test(name, operation):
    emit('start_test', name=name)
    start = time.monotonic()
    try:
        details = operation()
        item = {'name': name, 'status': 'passed', 'seconds': round(time.monotonic() - start, 2), 'details': details}
    except Exception as exc:
        item = {'name': name, 'status': 'failed', 'seconds': round(time.monotonic() - start, 2),
                'error_type': type(exc).__name__, 'error': str(exc)}
    report['tests'].append(item)
    save()
    emit('test_result', **item)
    return item


def wrapped_analysis(name, analysis, oss, extract):
    try:
        result = analysis.run(oss, oncomplete='Close')
    except Exception:
        if analysis.analysis is not None:
            analysis.analysis.Close()
        raise
    (OUT / f'{name}.json').write_text(result.to_json(), encoding='utf-8')
    if result.data is None:
        raise RuntimeError(f'Analysis returned no data; messages: {result.messages}')
    summary = extract(result.data)
    return {'messages': [str(m) for m in (result.messages or [])], **summary}


def dataframe_summary(df):
    if df.empty:
        raise RuntimeError('Analysis returned an empty data table')
    numeric = df.select_dtypes(include='number').to_numpy()
    finite_count = int(np.isfinite(numeric).sum())
    if numeric.size and finite_count == 0:
        raise RuntimeError('Analysis returned no finite numeric values')
    return {'shape': list(df.shape), 'columns': [str(c) for c in df.columns],
            'finite_numeric_values': finite_count, 'nonfinite_numeric_values': int(numeric.size - finite_count),
            'first_rows': json.loads(df.head(3).to_json(orient='split'))}


def native_analysis(name, kind, oss, settings=None, require_series=False):
    analysis = zp.analyses.new_analysis(oss, getattr(zp.constants.Analysis.AnalysisIDM, kind))
    try:
        for key, value in (settings or {}).items():
            setattr(analysis.Settings, key, value)
        analysis.ApplyAndWaitForCompletion()
        text = analysis.get_text_output(str(OUT / f'{name}.txt'), zos.get_txtfile_encoding())
        if not text.strip():
            raise RuntimeError('Native analysis returned empty text')
        lines = [s.strip() for s in text.splitlines() if s.strip()]
        series = [zp.utils.zputils.unpack_dataseries(analysis.Results.DataSeries[i])
                  for i in range(analysis.Results.NumberOfDataSeries)]
        if require_series and not series:
            raise RuntimeError('Native analysis returned no data series')
        for i, df in enumerate(series):
            df.to_csv(OUT / f'{name}_{i}.csv', encoding='utf-8-sig')
        return {'messages': [str(m) for m in analysis.messages], 'text_lines': len(lines),
                'text_excerpt': lines[:8] + lines[-5:], 'series_count': len(series),
                'series': [dataframe_summary(df) for df in series]}
    finally:
        analysis.Close()


try:
    emit('connect', mode='extension')
    zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    oss = zos.connect('extension', instance_number=0, message_logging=False)
    report['connection'] = {'license_status': str(zos.Application.LicenseStatus),
                            'api_license_valid': bool(zos.Application.IsValidLicenseForAPI),
                            'application_mode': str(zos.Application.Mode)}
    emit('connected', **report['connection'])
    before = next((t for t in report['tests'] if t['name'] == 'lens_and_system_read'), None)
    if before is None:
        before = test('lens_and_system_read', lambda: snapshot(oss))
    else:
        current = snapshot(oss)
        if current != before['details']:
            report.setdefault('previous_baselines', []).append(before['details'])
            before['details'] = current
            emit('new_session_baseline', explanation='Using the current lens state as the baseline for this run', details=current)
            save()
    if before['status'] != 'passed':
        raise RuntimeError('Cannot establish a baseline of lens data')

    if not any(t['name'] == 'single_ray_trace' and t['status'] == 'passed' for t in report['tests']):
        test('single_ray_trace', lambda: wrapped_analysis('single_ray_trace',
             zospy_zh.SingleRayTrace(px=0, py=0.5, wavelength=1), oss,
             lambda data: {'units': data.units, 'real_ray': dataframe_summary(data.real_ray_trace_data),
                           'paraxial_ray': dataframe_summary(data.paraxial_ray_trace_data)}))
    if not any(t['name'] == 'ray_fan' and t['status'] == 'passed' for t in report['tests']):
        test('ray_fan', lambda: wrapped_analysis('ray_fan',
             zospy_zh.RayFan(number_of_rays=10), oss,
             lambda data: dataframe_summary(data.to_dataframe())))
    test('single_ray_trace_native', lambda: native_analysis('single_ray_trace_native', 'RayTrace', oss,
         settings={'Hx': 0.0, 'Hy': 0.0, 'Px': 0.0, 'Py': 0.5}))
    test('ray_fan_native', lambda: native_analysis('ray_fan_native', 'RayFan', oss, require_series=True))

    def spot():
        analysis = zp.analyses.new_analysis(oss, zp.constants.Analysis.AnalysisIDM.StandardSpot)
        try:
            analysis.ApplyAndWaitForCompletion()
            text = analysis.get_text_output(str(OUT / 'standard_spot.txt'), zos.get_txtfile_encoding())
            if not text.strip():
                raise RuntimeError('Spot diagram returned empty text')
            lines = [s.strip() for s in text.splitlines() if s.strip()]
            metrics = [s for s in lines if any(k in s.lower() for k in ('rms', 'radius', 'airy', '半径', '均方根'))]
            return {'messages': [str(m) for m in analysis.messages], 'text_lines': len(lines),
                    'metrics': metrics[:12], 'text_excerpt': lines[-12:]}
        finally:
            analysis.Close()

    if not any(t['name'] == 'standard_spot' and t['status'] == 'passed' for t in report['tests']):
        test('standard_spot', spot)
    test('fft_mtf', lambda: wrapped_analysis('fft_mtf',
         zp.analyses.mtf.FFTMTF(sampling='64x64', maximum_frequency=100, show_diffraction_limit=True),
         oss, dataframe_summary))
    test('system_data_report', lambda: wrapped_analysis('system_data',
         zospy_zh.SystemData(), oss, lambda data: {'data_type': type(data).__name__}))
    if report['tests'][-1]['status'] != 'passed':
        test('system_data_report_native', lambda: native_analysis('system_data_native', 'SystemData', oss))
    after = test('lens_unchanged', lambda: {'unchanged': snapshot(oss) == before['details']})
    if after['status'] == 'passed' and not after['details']['unchanged']:
        after['status'] = 'failed'
        after['error'] = 'Lens snapshot changed during analysis'
    save()
    emit('completed', passed=sum(t['status'] == 'passed' for t in report['tests']), total=len(report['tests']))
except Exception as exc:
    report['fatal_error'] = {'type': type(exc).__name__, 'message': str(exc)}
    if zos is not None and zos.Application is not None:
        try:
            report['fatal_error']['license_status'] = str(zos.Application.LicenseStatus)
        except Exception:
            pass
    save()
    emit('fatal_error', **report['fatal_error'])
finally:
    if zos is not None and zos.Application is not None:
        emit('disconnect')
        zos.disconnect()
        emit('disconnected')
