"""Audit the published prescription against the patent's paraxial specifications."""
import itertools
import json
import math
from pathlib import Path

import numpy as np

from import_patent_zospy import compile_prescription


ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'zemax_import' / 'CN110703413B_example3' / 'mtf_diagnosis'
ROWS = compile_prescription()


def metrics(rows):
    matrix = np.eye(2)
    index = 1.
    for surface, row in enumerate(rows):
        following = 1. if row['glass'] is None else row['glass']['nd']
        radius = math.inf if row['radius'] == 'Infinity' else row['radius']
        matrix = np.array([[1., 0.], [-(following - index) / radius, 1.]]) @ matrix
        index = following
        if surface < len(rows) - 1:
            matrix = np.array([[1., row['thickness'] / index], [0., 1.]]) @ matrix
    finite = matrix @ np.array([[1., 500.], [0., 1.]])
    return {'efl_mm': -1. / matrix[1, 0], 'bfl_infinite_mm': -matrix[0, 0] / matrix[1, 0],
            'last_surface_to_paraxial_focus_WD500_mm': -finite[0, 1] / finite[1, 1]}


def lens_focal_length(r1, r2, thickness, nd):
    power = (nd - 1.) * (1. / r1 - 1. / r2 + (nd - 1.) * thickness / (nd * r1 * r2))
    return 1. / power


base = metrics(ROWS)
assert math.isclose(base['efl_mm'], 29.86665, abs_tol=1e-4)
result = {'method': 'Independent reduced-angle ABCD and thick-lens lensmaker calculations at Nd',
          'source': 'CN110703413B Table7 and Table8; no original or live model modified',
          'published_system_efl_mm': 35., 'published_lens7_efl_mm': -41.7,
          'published_lens8_efl_mm': 37.19, 'published_ratio': -1.12,
          'calculated_system': base,
          'calculated_lens7_efl_mm': lens_focal_length(114.24, 24.8, 1.11, 1.76),
          'calculated_lens8_efl_mm': lens_focal_length(32.62, 242.95, 4.065, 1.77),
          'current_ui': {'configuration': 1, 'image_distance_mm_displayed': 7.634,
                         'efl_mm_displayed': 29.8666,
                         'basis': 'UI text and screenshot captured during this diagnosis',
                         'image_distance_minus_paraxial_focus_mm': 7.634 - base['last_surface_to_paraxial_focus_WD500_mm'],
                         'warning': 'Near-axis focus is not polychromatic best focus; spherical and chromatic aberration may change the best plane.'}}
result['calculated_ratio'] = result['calculated_lens7_efl_mm'] / result['calculated_lens8_efl_mm']
glass_rows = [i for i, row in enumerate(ROWS) if row['glass'] is not None]
efls = []
for signs in itertools.product([-1, 1], repeat=len(glass_rows)):
    rows = json.loads(json.dumps(ROWS))
    for i, sign in zip(glass_rows, signs):
        rows[i]['glass']['nd'] += sign * .005
    efls.append(metrics(rows)['efl_mm'])
result['nd_rounding_corner_check'] = {'rounding_half_range': .005, 'corners': len(efls),
                                     'minimum_efl_mm': min(efls), 'maximum_efl_mm': max(efls),
                                     'scope': 'Nd rounding only; corner check, not a full continuous interval proof'}
OUT.mkdir(exist_ok=True)
(OUT / 'patent_consistency_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=True, indent=2), flush=True)
