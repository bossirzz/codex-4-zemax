"""Read-only stop-position sensitivity check on an independent Zemax copy."""
import os
import json
import math
import re
import sys
from pathlib import Path

import zospy as zp

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'zemax_import' / 'CN110703413B_example3' / 'stop_analysis'
OUT.mkdir(exist_ok=True)


def snapshot(oss):
    return {'file': str(oss.SystemFile), 'needs_save': bool(oss.NeedsSave),
            'configuration': int(oss.MCE.CurrentConfiguration),
            'surfaces': [[str(getattr(oss.LDE.GetSurfaceAt(i), key))
                          for key in ['Radius', 'Thickness', 'Material', 'SemiDiameter', 'Comment']]
                         for i in range(oss.LDE.NumberOfSurfaces)]}


def metric(text, label):
    match = re.search(r'^' + re.escape(label) + r'\s*:\s*([+\-\d.Ee]+)', text, re.MULTILINE)
    if match is None:
        raise ValueError(f'Missing system data metric: {label}')
    return float(match[1])


def offline_check():
    import numpy as np
    from import_patent_zospy import compile_prescription
    rows = compile_prescription()
    result = {'method': 'independent paraxial ABCD calculation using published Nd',
              'live_zemax_sweep_performed': False,
              'stop_surface': 7, 'configuration': 'WD500', 'total_gap_mm': 8.76, 'samples': []}
    for distance in [1.0, 4.38, 7.76]:
        matrix = np.eye(2)
        index = 1.0
        prefix = None
        for surface, row in enumerate(rows, 1):
            next_index = 1.0 if row['glass'] is None else row['glass']['nd']
            radius = math.inf if row['radius'] == 'Infinity' else row['radius']
            matrix = np.array([[1., 0.], [-(next_index - index) / radius, 1.]]) @ matrix
            index = next_index
            if surface == 7:
                prefix = matrix.copy()
            thickness = distance if surface == 6 else 8.76 - distance if surface == 7 else row['thickness']
            if surface < 16:
                matrix = np.array([[1., thickness / index], [0., 1.]]) @ matrix
        rear = matrix @ np.linalg.inv(prefix)
        efl = -1.0 / matrix[1, 0]
        epd = efl / 2.8
        item = {'L3R2_to_stop_mm': distance, 'stop_to_L4R1_mm': 8.76 - distance,
                'efl_mm': efl, 'paraxial_bfl_mm': -matrix[0, 0] / matrix[1, 0],
                'stop_radius_mm': abs(prefix[0, 0]) * epd / 2,
                'entrance_pupil_diameter_mm': epd,
                'entrance_pupil_position_mm': prefix[0, 1] / prefix[0, 0],
                'exit_pupil_position_relative_image_mm': -rear[0, 1] / rear[1, 1] - 14.65}
        result['samples'].append(item)
    assert max(s['efl_mm'] for s in result['samples']) - min(s['efl_mm'] for s in result['samples']) < 1e-9
    # Compare the midpoint prediction with the already saved native Zemax report.
    midpoint = result['samples'][1]
    result['midpoint_agrees_with_saved_zemax_report'] = (
        math.isclose(midpoint['efl_mm'], 29.86665, abs_tol=1e-4)
        and math.isclose(midpoint['stop_radius_mm'], 3.237858, abs_tol=1e-4))
    assert result['midpoint_agrees_with_saved_zemax_report']
    (OUT / 'paraxial_stop_position_check.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(result, indent=2, ensure_ascii=True), flush=True)


if '--live' not in sys.argv:
    offline_check()
    raise SystemExit(0)


zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
copy = None
try:
    primary = zos.connect('extension', instance_number=0, message_logging=False)
    before = snapshot(primary)
    expected = ROOT / 'zemax_import' / 'CN110703413B_example3' / 'CN110703413B_example3_modelglass.zos'
    if Path(str(primary.SystemFile)).resolve() != expected:
        raise RuntimeError('Current lens is no longer the patent model; refusing to run on a different lens.')
    copy = primary.copy_system()
    assert copy.SystemID != primary.SystemID
    copy.MCE.SetCurrentConfiguration(1)
    assert copy.LDE.GetSurfaceAt(6).Comment == 'L3R2'
    assert copy.LDE.GetSurfaceAt(7).IsStop
    assert copy.LDE.GetSurfaceAt(8).Comment == 'L4R1'
    result = {'primary_file': str(primary.SystemFile), 'stop_surface': 7,
              'before_stop_surface': 'L3R2', 'after_stop_surface': 'L4R1',
              'total_gap_mm': 8.76, 'configuration': 'WD500', 'samples': []}
    for distance in [1.0, 4.38, 7.76]:
        copy.LDE.GetSurfaceAt(6).Thickness = distance
        copy.LDE.GetSurfaceAt(7).Thickness = 8.76 - distance
        analysis = zp.analyses.new_analysis(copy, zp.constants.Analysis.AnalysisIDM.SystemData)
        try:
            analysis.ApplyAndWaitForCompletion()
            text = analysis.get_text_output(str(OUT / f'stop_{distance:.2f}mm_system_data.txt'),
                                            zos.get_txtfile_encoding())
            item = {'L3R2_to_stop_mm': distance, 'stop_to_L4R1_mm': 8.76 - distance,
                    'efl_mm': metric(text, 'Effective Focal Length'),
                    'paraxial_bfl_mm': metric(text, 'Back Focal Length'),
                    'stop_radius_mm': metric(text, 'Stop Radius'),
                    'entrance_pupil_diameter_mm': metric(text, 'Entrance Pupil Diameter'),
                    'exit_pupil_position_mm': metric(text, 'Exit Pupil Position')}
            result['samples'].append(item)
        finally:
            analysis.Close()
    result['primary_unchanged'] = snapshot(primary) == before
    assert result['primary_unchanged']
    assert max(s['efl_mm'] for s in result['samples']) - min(s['efl_mm'] for s in result['samples']) < 1e-6
    (OUT / 'stop_position_check.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(result, indent=2, ensure_ascii=True), flush=True)
finally:
    if copy is not None:
        copy.close(saveifneeded=False)
    if zos.Application is not None:
        zos.disconnect()
