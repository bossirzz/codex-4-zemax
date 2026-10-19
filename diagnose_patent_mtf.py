"""Read-only MTF diagnosis: saved ray data, paraxial focus, optional live settings."""
import os
import json
import math
import sys
from pathlib import Path

import numpy as np

from import_patent_zospy import compile_prescription


ROOT = Path(__file__).resolve().parent
MODEL = ROOT / 'zemax_import' / 'CN110703413B_example3'
OUT = MODEL / 'mtf_diagnosis'


def offline():
    matrix = np.eye(2)
    index = 1.0
    rows = compile_prescription()
    for surface, row in enumerate(rows, 1):
        next_index = 1.0 if row['glass'] is None else row['glass']['nd']
        radius = math.inf if row['radius'] == 'Infinity' else row['radius']
        matrix = np.array([[1., 0.], [-(next_index - index) / radius, 1.]]) @ matrix
        index = next_index
        if surface < len(rows):
            matrix = np.array([[1., row['thickness'] / index], [0., 1.]]) @ matrix
    finite = matrix @ np.array([[1., 500.], [0., 1.]])
    paraxial_focus = -finite[0, 1] / finite[1, 1]
    trace = json.loads((MODEL / 'ray_trace.json').read_text(encoding='utf-8'))['data']
    data = trace['real_ray_trace_data']
    image = dict(zip(data['columns'], data['data'][-1]))
    dy = image['Y-cosine']
    dz = image['Z-cosine']
    y = image['Y-coordinate']
    ray_focus_shift = -y / (dy / dz)
    return {
        'basis': 'Saved WD500 model and native ray trace; not a measurement of the current UI',
        'configuration': 'WD500', 'efl_mm': -1. / matrix[1, 0],
        'last_surface_to_current_image_mm': rows[-1]['thickness'],
        'paraxial_last_surface_to_finite_conjugate_focus_mm': paraxial_focus,
        'paraxial_image_shift_mm': paraxial_focus - rows[-1]['thickness'],
        'ray_wavelength_um': trace['wavelength']['value'],
        'ray_normalized_pupil_y': trace['normalized_y_pupil_coord'],
        'ray_y_at_current_image_mm': y,
        'ray_axis_intersection_shift_from_current_image_mm': ray_focus_shift,
        'ray_last_surface_to_axis_intersection_mm': rows[-1]['thickness'] + ray_focus_shift,
        'warning': 'Paraxial and single-ray focus are diagnostic estimates, not optimized polychromatic best focus.'
    }


def live():
    import zospy as zp
    zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    try:
        oss = zos.connect('extension', instance_number=0, message_logging=False)
        result = {'file': str(oss.SystemFile), 'configuration': int(oss.MCE.CurrentConfiguration),
                  'needs_save': bool(oss.NeedsSave), 'analyses_members':
                  sorted({str(m.Name) for m in oss.Analyses.GetType().GetMembers()})}
        windows = []
        for i in range(int(oss.Analyses.NumberOfAnalyses)):
            window = oss.Analyses.Get_AnalysisAtIndex(i)
            settings = window.GetSettings()
            item = {'index': i, 'settings_type': str(settings.GetType())}
            for key in ['SampleSize', 'PupilSampleSize', 'ImageSampleSize', 'MaximumFrequency',
                        'ImageDelta', 'ShowDiffractionLimit', 'Type', 'Surface', 'Field', 'Wavelength']:
                if hasattr(settings, key):
                    item[key] = str(getattr(settings, key))
            windows.append(item)
        result['open_analyses'] = windows
        return result
    finally:
        if zos.Application is not None:
            zos.disconnect()


if __name__ == '__main__':
    result = {'saved_model_diagnosis': offline()}
    if '--live' in sys.argv:
        try:
            result['live_settings'] = live()
        except Exception as exc:
            result['live_settings_error'] = f'{type(exc).__name__}: {exc}'
    OUT.mkdir(exist_ok=True)
    (OUT / 'mtf_diagnosis.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(result, indent=2, ensure_ascii=True), flush=True)
