import os
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import zospy as zp

OUT = Path(__file__).with_name('zospy_image_results')
OUT.mkdir(exist_ok=True)
EXTRA = '--extra' in sys.argv
report = json.loads((OUT / 'report.json').read_text(encoding='utf-8')) if EXTRA else {'zospy_version': zp.__version__, 'tests': []}
report.pop('fatal_error', None)
zos = None


def emit(stage, **data):
    print(json.dumps({'stage': stage, **data}, ensure_ascii=True, default=str), flush=True)


def save():
    (OUT / 'report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding='utf-8')


def safe_number(value):
    value = float(value)
    return value if math.isfinite(value) else str(value)


def snapshot(oss):
    rows = []
    for i in range(int(oss.LDE.NumberOfSurfaces)):
        row = oss.LDE.GetSurfaceAt(i)
        rows.append({'surface': i, 'radius': safe_number(row.Radius), 'thickness': safe_number(row.Thickness),
                     'material': str(row.Material), 'semi_diameter': safe_number(row.SemiDiameter),
                     'conic': safe_number(row.Conic), 'comment': str(row.Comment), 'is_stop': bool(row.IsStop)})
    return {'file': str(oss.SystemFile), 'rows': rows,
            'fields': int(oss.SystemData.Fields.NumberOfFields),
            'wavelengths': int(oss.SystemData.Wavelengths.NumberOfWavelengths)}


def test(name, operation):
    emit('start_test', name=name)
    start = time.monotonic()
    try:
        details = operation()
        record = {'name': name, 'status': 'passed', 'seconds': round(time.monotonic() - start, 2), 'details': details}
    except Exception as exc:
        record = {'name': name, 'status': 'failed', 'seconds': round(time.monotonic() - start, 2),
                  'error_type': type(exc).__name__, 'error': str(exc)}
    report['tests'].append(record)
    save()
    emit('test_result', **record)


def export_viewer(name, viewer, oss):
    path = OUT / f'{name}.png'
    try:
        result = viewer.run(oss, image_output_file=path, oncomplete='Close')
        if not path.is_file() or path.stat().st_size < 100:
            raise RuntimeError('No usable image file was exported')
        return {'file': str(path), 'bytes': path.stat().st_size, 'messages': [str(m) for m in (result.messages or [])]}
    except Exception:
        if viewer.analysis is not None:
            viewer.analysis.Close()
        raise


curves = {}
if EXTRA:
    for size in (64, 256, 512, 1024):
        path = OUT / f'mtf_{size}.json'
        if path.exists():
            curves[size] = zp.analyses.base.AnalysisResult.from_json(path.read_text(encoding='utf-8')).data


def mtf(size, oss):
    wrapper = zp.analyses.mtf.FFTMTF(sampling=f'{size}x{size}', maximum_frequency=100, show_diffraction_limit=True)
    try:
        result = wrapper.run(oss, oncomplete='Close')
    except Exception:
        if wrapper.analysis is not None:
            wrapper.analysis.Close()
        raise
    df = result.data
    if df is None or df.empty:
        raise RuntimeError('No MTF data returned')
    values = df.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError('MTF contains nonfinite data')
    df.to_csv(OUT / f'mtf_{size}.csv', encoding='utf-8-sig')
    (OUT / f'mtf_{size}.json').write_text(result.to_json(), encoding='utf-8')
    curves[size] = df
    details = {'sampling': size, 'points': len(df), 'curves': df.shape[1],
               'messages': [str(m) for m in (result.messages or [])],
               'range': [float(values.min()), float(values.max())],
               'zero_frequency': values[0].tolist()}
    previous = next((s for s in sorted(curves, reverse=True) if s < size), None)
    if previous is not None:
        old = curves[previous]
        if list(df.columns) == list(old.columns) and np.allclose(df.index, old.index):
            difference = np.abs(values - old.to_numpy(dtype=float))
            details['compared_with_sampling'] = previous
            details['max_absolute_difference'] = float(difference.max())
            details['mean_absolute_difference'] = float(difference.mean())
    return details


try:
    emit('connect', mode='extension')
    zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    oss = zos.connect('extension', instance_number=0, message_logging=False)
    report['connection'] = {'license_status': str(zos.Application.LicenseStatus),
                            'api_license_valid': bool(zos.Application.IsValidLicenseForAPI)}
    emit('connected', **report['connection'])
    before = snapshot(oss)
    if EXTRA:
        report['extra_before'] = before
    else:
        report['before'] = before
        report['needs_save_before'] = bool(oss.NeedsSave)
    save()
    if EXTRA:
        test('cross_section_correct_dimensions', lambda: export_viewer('cross_section_960x640',
             zp.analyses.systemviewers.CrossSection(number_of_rays=7, image_size=(640, 960)), oss))
    else:
        test('cross_section_image', lambda: export_viewer('cross_section',
             zp.analyses.systemviewers.CrossSection(number_of_rays=7, image_size=(960, 640)), oss))
        test('viewer3d_image', lambda: export_viewer('viewer3d',
             zp.analyses.systemviewers.Viewer3D(number_of_rays=3, image_size=(960, 640),
                 camera_viewpoint_angle_x=-15, camera_viewpoint_angle_y=30), oss))
    for size in ((2048, 4096) if EXTRA else (64, 256, 512, 1024)):
        test(f'fft_mtf_{size}', lambda size=size: mtf(size, oss))
    report['after'] = snapshot(oss)
    report['lens_parameters_unchanged'] = report['after'] == before
    report['needs_save_after'] = bool(oss.NeedsSave)
    save()
    emit('completed', passed=sum(t['status'] == 'passed' for t in report['tests']),
         total=len(report['tests']), lens_parameters_unchanged=report['lens_parameters_unchanged'])
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
