"""Regression checks against actual Chinese Zemax output; optional live checks."""
import os
import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
import warnings

import numpy as np
import zospy as zp
from zospy.analyses.raysandspots.single_ray_trace import SingleRayTraceResult, SingleRayTraceTransformer
from zospy.analyses.reports.system_data import SystemDataResult, SystemDataTransformer
import zospy_zh as zh

ROOT = Path(__file__).parent
FIXTURES = ROOT / 'zospy_test_results'
OUT = ROOT / 'zospy_zh_results'
OUT.mkdir(exist_ok=True)


def emit(obj):
    print(json.dumps(obj, ensure_ascii=True, default=str), flush=True)


def offline():
    ray_text = (FIXTURES / 'single_ray_trace_native.txt').read_text(encoding='utf-16')
    system_text = (FIXTURES / 'system_data_native.txt').read_text(encoding='utf-16')
    ray = zh.SingleRayTrace(px=0, py=0.5)
    ray.get_text_output = lambda: ray_text
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Header and row length mismatch.*')
        rays = ray.parse_output('single_ray_trace', SingleRayTraceTransformer, SingleRayTraceResult)
    assert rays.wavelength.value == 0.55 and rays.normalized_y_pupil_coord == 0.5
    # Compare every reported number, including infinities, for every real/paraxial row.
    section = None
    compared = 0
    for line in ray_text.splitlines():
        if line.strip() == '实际光线追迹数据:':
            section = rays.real_ray_trace_data
        elif line.strip() == '近轴光线追迹数据:':
            section = rays.paraxial_ray_trace_data
        elif section is not None and (line.strip().startswith('物面') or line.strip().split()[:1] in [['1'], ['2']]):
            tokens = line.split()
            surf = 'OBJ' if tokens[0] == '物面' else int(tokens[0])
            row = section[section['Surf'] == surf].iloc[0]
            for i, token in enumerate(tokens[1:], start=1):
                if token == '-':
                    assert math.isnan(row.iloc[i])
                else:
                    expected = math.inf if token == '无限' else float(token)
                    assert row.iloc[i] == expected, (surf, i, row.iloc[i], expected)
                compared += 1
    system = zh.SystemData()
    system.get_text_output = lambda: system_text
    data = system.parse_output('system_data', SystemDataTransformer, SystemDataResult)
    assert data.general_lens_data.surfaces == 2
    assert data.general_lens_data.stop == 1
    assert data.general_lens_data.system_aperture.value == 3
    assert data.general_lens_data.temperature == 20
    assert data.general_lens_data.effective_focal_length_air == 90.9067
    assert data.general_lens_data.object_space_na == 3e-10
    assert data.general_lens_data.ray_aiming is False
    assert data.wavelengths.wavelengths[0].value == 0.55
    assert data.abcd_matrix.A == 1.335 and data.abcd_matrix.D == 1.335
    series = []
    for i in range(2):
        with (FIXTURES / f'ray_fan_native_{i}.csv').open(encoding='utf-8-sig', newline='') as f:
            rows = list(csv.reader(f))
        nums = np.array([[float(v) for v in row] for row in rows[3:]])
        series.append(SimpleNamespace(Description=rows[0][1], SeriesLabels=rows[1][1:], XLabel=rows[2][0],
                                      XData=SimpleNamespace(Data=nums[:, 0]), YData=SimpleNamespace(Data=nums[:, 1:])))
    fan = zh.RayFan()
    fan._analysis = SimpleNamespace(Results=SimpleNamespace(NumberOfDataSeries=2, GetDataSeries=lambda i: series[i]))
    fan_data = fan.get_data_series()
    for actual, expected in zip([fan_data.tangential[0], fan_data.sagittal[0]], series):
        np.testing.assert_array_equal(actual.data.to_numpy(), expected.YData.Data)
        np.testing.assert_array_equal(actual.data.index.to_numpy(), expected.XData.Data)
        assert actual.field_number == 1 and actual.field_coordinate.value == (0, 0)
    # Unicode user paths/comments are not translated by the header adapter.
    original = '光线追迹数据\n文件: D:\\设计\\无限\\面.zos\n题目: 主波长：测试\n'
    normalized = zh.normalize_text(original, 'single_ray_trace')
    assert 'D:\\设计\\无限\\面.zos' in normalized and '主波长：测试' in normalized
    english = 'Ray Trace Data\nUnits: Millimeters\nCoordinates: Local\nWavelength: 0.55 µm\n'
    assert zh.normalize_text(english, 'single_ray_trace') == english
    assert zh.normalize_fan_description('Tangential fan, field number 1 = 0.0000 (deg)') == 'Tangential fan, field number 1 = 0.0000 (deg)'
    assert zh.normalize_fan_description('子午光扇，视场号2 = -1.2500, 2.5000 (度)').endswith('-1.2500, 2.5000 (度)')
    result = {'status': 'passed', 'analyses': ['SingleRayTrace', 'RayFan', 'SystemData'],
              'ray_numeric_cells_compared': compared, 'fan_values_compared': sum(s.YData.Data.size for s in series),
              'system_required_fields_validated': True, 'unicode_user_text_preserved': True}
    (OUT / 'offline_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    emit(result)


def snapshot(oss):
    return {'file': str(oss.SystemFile), 'needs_save': bool(oss.NeedsSave),
            'surfaces': [[str(getattr(oss.LDE.GetSurfaceAt(i), key))
                          for key in ['Radius', 'Thickness', 'Material', 'SemiDiameter', 'Conic', 'Comment']]
                         for i in range(oss.LDE.NumberOfSurfaces)],
            'fields': [[float(oss.SystemData.Fields.GetField(i).X), float(oss.SystemData.Fields.GetField(i).Y)]
                       for i in range(1, oss.SystemData.Fields.NumberOfFields + 1)],
            'wavelengths': [float(oss.SystemData.Wavelengths.GetWavelength(i).Wavelength)
                            for i in range(1, oss.SystemData.Wavelengths.NumberOfWavelengths + 1)]}


def live():
    zos = None
    report = {'status': 'running', 'tests': []}
    try:
        emit({'stage': 'connecting'})
        zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
        oss = zos.connect('extension', instance_number=0, message_logging=False)
        before = snapshot(oss)
        emit({'stage': 'connected', 'file': before['file']})
        for name, analysis in [('SingleRayTrace', zh.SingleRayTrace(py=0.5)),
                               ('RayFan', zh.RayFan(number_of_rays=10)), ('SystemData', zh.SystemData())]:
            try:
                result = analysis.run(oss, oncomplete='Close')
                assert result.data is not None
                (OUT / f'{name}.json').write_text(result.to_json(), encoding='utf-8')
                report['tests'].append({'name': name, 'status': 'passed'})
                emit(report['tests'][-1])
            finally:
                if analysis.analysis is not None:
                    analysis.analysis.Close()
        report['lens_unchanged'] = snapshot(oss) == before
        assert report['lens_unchanged']
        report['status'] = 'passed'
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = {'type': type(exc).__name__, 'message': str(exc)}
    finally:
        if zos is not None and zos.Application is not None:
            zos.disconnect()
        (OUT / 'live_report.json').write_text(json.dumps(report, ensure_ascii=True, indent=2), encoding='utf-8')
        emit(report)
    if report['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--live', action='store_true')
    args = parser.parse_args()
    live() if args.live else offline()
