"""Read-only reproduction diagnostics: no API calls or Zemax lens edits.

Independent reduced-angle ABCD and meridional Snell calculations. Hypothesis
changes are evidence probes, not corrections certified as the original design.
"""
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from import_patent_zospy import compile_prescription, CONFIGS

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / 'zemax_import' / 'CN110703413B_example3'
OUT = MODEL / 'reproduction_audit'


def paraxial(rows, wd=500):
    matrix = np.eye(2)
    n = 1.
    prefix = None
    for number, row in enumerate(rows):
        nn = 1. if row['glass'] is None else row['glass']['nd']
        radius = math.inf if row['radius'] == 'Infinity' else row['radius']
        matrix = np.array([[1., 0.], [-(nn - n) / radius, 1.]]) @ matrix
        n = nn
        if row.get('stop'):
            prefix = matrix.copy()
        if number < len(rows) - 1:
            matrix = np.array([[1., row['thickness'] / n], [0., 1.]]) @ matrix
    a, b, c, d = matrix.flat
    result = {'efl_mm': -1. / c, 'bfl_infinite_mm': -a / c,
              'matrix': matrix.tolist(), 'matrix_determinant': float(np.linalg.det(matrix))}
    result['paraxial_focus_from_last_vertex_mm'] = (-a / c if math.isinf(wd)
                                                   else -(a * wd + b) / (c * wd + d))
    q = rows[-1]['thickness']
    result['published_image_distance_mm'] = q
    result['published_image_minus_paraxial_focus_mm'] = q - result['paraxial_focus_from_last_vertex_mm']
    result['wd_required_for_published_image_mm'] = -(b + q * d) / (a + q * c)
    if prefix is not None:
        result['entrance_pupil_position_mm'] = prefix[0, 1] / prefix[0, 0]
    assert abs(result['matrix_determinant'] - 1.) < 1e-10
    return result


def focal(r1, r2, thickness, n):
    return 1. / ((n - 1.) * (1. / r1 - 1. / r2
                             + (n - 1.) * thickness / (n * r1 * r2)))


def trace(rows, wd, launch_tan):
    p = np.array([0., 0., -wd])
    direction = np.array([0., launch_tan, 1.])
    direction /= np.linalg.norm(direction)
    vertex = 0.
    n = 1.
    for number, row in enumerate(rows):
        nn = 1. if row['glass'] is None else row['glass']['nd']
        radius = math.inf if row['radius'] == 'Infinity' else row['radius']
        if math.isinf(radius):
            hit = p + (vertex - p[2]) / direction[2] * direction
            normal = np.array([0., 0., -1.])
        else:
            center = np.array([0., 0., vertex + radius])
            offset = p - center
            b = np.dot(direction, offset)
            discriminant = b * b - np.dot(offset, offset) + radius * radius
            if discriminant < 0:
                raise ValueError('Missed sphere')
            root = math.sqrt(discriminant)
            distance = -b - root if radius > 0 else -b + root
            if distance < -1e-6:
                raise ValueError('Backward intersection')
            hit = p + distance * direction
            normal = (hit - center) / abs(radius)
        if np.dot(normal, direction) > 0:
            normal = -normal
        cosine = -np.dot(normal, direction)
        eta = n / nn
        k = 1. - eta * eta * (1. - cosine * cosine)
        if k < 0:
            raise ValueError('Total internal reflection')
        direction = eta * direction + (eta * cosine - math.sqrt(k)) * normal
        direction /= np.linalg.norm(direction)
        p = hit
        n = nn
        if number < len(rows) - 1:
            vertex += row['thickness']
    slope = direction[1] / direction[2]
    return float(p[1] + (vertex - p[2]) * slope), float(slope)


def spot(rows, native_calibration, wd=500, fn=2.8, sample_count=600):
    metrics = paraxial(rows, wd)
    launch = native_calibration * metrics['efl_mm'] / (2. * fn * (wd + metrics['entrance_pupil_position_mm']))
    radii = (np.arange(sample_count) + .5) / sample_count
    rays = np.array([trace(rows, wd, launch * float(r)) for r in radii])
    weights = 2. * radii
    weights /= weights.sum()
    heights, slopes = rays.T
    best = -np.sum(weights * heights * slopes) / np.sum(weights * slopes * slopes)

    def rms(q):
        return float(np.sqrt(np.sum(weights * (heights + q * slopes) ** 2)) * 1000.)

    intersections = []
    for r in [.05, .25, .5, .75, 1.]:
        y, u = trace(rows, wd, launch * r)
        intersections.append({'normalized_pupil_radius': r, 'focus_from_last_vertex_mm': -y / u})
    return {'scope': 'Monochromatic on-axis geometric estimate at Nd, uniform paraxial entrance pupil, ignoring clipping; not polychromatic MTF',
            'image_space_f_number': fn, 'wd_mm': wd, 'rays': sample_count,
            'input_marginal_tan': launch,
            'geometric_best_focus_from_last_vertex_mm': float(best),
            'rms_at_geometric_best_focus_um': rms(best),
            'rms_at_source_14_65mm_um': rms(14.65),
            'rms_at_current_ui_7_634mm_um': rms(7.634),
            'axis_intersections': intersections,
            'axial_intersection_spread_nearaxis_to_marginal_mm': intersections[0]['focus_from_last_vertex_mm'] - intersections[-1]['focus_from_last_vertex_mm']}


def source_controls():
    # Manual transcription of the rendered PDF tables, retained for traceability.
    ex1 = [(28.85, 2.6, 1.77), (-183.07, .16, None),
           (12.5, 2.42, 1.59), (22.9, 1.09, None),
           (119.8, 1.85, 1.57), (9.08, 5.5, None),
           (-9.03, .89, 1.72), (55.08, 0., None),
           (55.08, 3.91, 1.8), (-13.13, .1, None),
           (51.15, 2.2, 1.77), (-29.08, 2.98, None),
           (25.06, 3.02, 1.78), (15.2, .88, None),
           (25.36, 1.6, 1.77), (79.6, 12.76, None)]
    ex2 = [(34.4, 2.04, 1.77), (116.99, .15, None),
           (8.72, 2.75, 1.59), (25.05, .16, None),
           (24.08, .95, 1.57), (7.08, 8.38, None),
           (-5.76, .65, 1.72), (-20.13, 0., None),
           (-20.13, 2.45, 1.8), (-8.25, .26, None),
           (225.58, 2.13, 1.77), (-21.59, 2.36, None),
           (41.65, 1.87, 1.78), (19.56, .77, None),
           (26.85, 1.59, 1.77), (184.26, 12.85, None)]
    results = []
    for example, faces, gap500, f7, f8 in [(1, ex1, 1.95, -56.97, 46.88), (2, ex2, 2.35, -48.54, 39.75)]:
        rows = [{'radius': r, 'thickness': t, 'glass': None if n is None else {'nd': n}} for r, t, n in faces]
        literal = paraxial(rows)
        rows[11]['thickness'] = gap500
        wd500 = paraxial(rows)
        results.append({'example': example, 'nominal_efl_mm': 25.,
                        'literal_table_metrics': literal, 'Table_WD500_gap_mm': gap500,
                        'with_WD500_gap_metrics': wd500,
                        'published_f7_mm': f7, 'computed_f7_mm': focal(faces[12][0], faces[13][0], faces[12][1], faces[12][2]),
                        'published_f8_mm': f8, 'computed_f8_mm': focal(faces[14][0], faces[15][0], faces[14][1], faces[14][2])})
    return results


def main():
    OUT.mkdir(exist_ok=True)
    original_lens = MODEL / 'CN110703413B_example3_modelglass.zos'
    before_hash = hashlib.sha256(original_lens.read_bytes()).hexdigest()
    rows = compile_prescription()
    base = paraxial(rows)
    assert math.isclose(base['efl_mm'], 29.86665, abs_tol=1e-4)
    saved = json.loads((MODEL / 'ray_trace.json').read_text(encoding='utf-8'))['data']['real_ray_trace_data']
    obj = dict(zip(saved['columns'], saved['data'][0]))
    image = dict(zip(saved['columns'], saved['data'][-1]))
    native_launch = 2. * obj['Y-cosine'] / obj['Z-cosine']
    calculated_launch = base['efl_mm'] / (5.6 * (500. + base['entrance_pupil_position_mm']))
    calibration = native_launch / calculated_launch
    y, u = trace(rows, 500., native_launch / 2.)
    native_error = y + 14.65 * u - image['Y-coordinate']
    assert abs(native_error) < 1e-4

    matching = json.loads((MODEL / 'glass_matching' / 'glass_candidates.json').read_text(encoding='utf-8'))
    selected = json.loads((MODEL / 'glass_matching' / 'selected_comparison_glasses.json').read_text(encoding='utf-8'))
    domestic = {lens: item['candidate'] for item in selected for lens in item['lenses']}
    nearest = {lens: group['nearest_ordinary_name'][0] for group in matching['groups'] for lens in group['lenses']}

    def apply_glasses(trial, mapping, wavelength='d'):
        used = []
        for row in trial:
            if row['glass'] is None:
                continue
            lens = 'L5' if row['comment'] == 'L4R2_L5R1_cemented' else row['comment'][:2]
            candidate = mapping[lens]
            row['glass']['nd'] = candidate['indices'][wavelength]
            used.append({'lens': lens, 'catalog': candidate['catalog'], 'name': candidate['name'],
                         'wave': wavelength, 'index': row['glass']['nd']})
        return used

    scenarios = []
    for label in ['published_table7', 'domestic_catalog_Nd', 'nearest_ordinary_catalog_Nd',
                  'stop_after_L3_1mm', 'stop_after_L3_7_76mm', 'L6R1_minus115_9_hypothesis',
                  'L6R1_plus1159_hypothesis', 'L8R2_minus242_95_hypothesis', 'L8_Nd_fitted_to_f8_only']:
        trial = copy.deepcopy(rows)
        changed = []
        if label == 'domestic_catalog_Nd':
            changed = apply_glasses(trial, domestic)
        elif label == 'nearest_ordinary_catalog_Nd':
            changed = apply_glasses(trial, nearest)
        elif label.startswith('stop_after'):
            gap = 1. if label.endswith('1mm') else 7.76
            trial[5]['thickness'] = gap
            trial[6]['thickness'] = 8.76 - gap
            changed = [{'parameter': 'L3R2_to_stop_mm', 'value': gap, 'total_L3_L4_gap_mm': 8.76}]
        elif label.startswith('L6R1'):
            value = -115.9 if 'minus' in label else 1159.
            trial[10]['radius'] = value
            changed = [{'parameter': 'L6R1', 'source': -1159, 'hypothesis': value}]
        elif label.startswith('L8R2'):
            trial[-1]['radius'] = -242.95
            changed = [{'parameter': 'L8R2', 'source': 242.95, 'hypothesis': -242.95}]
        elif label == 'L8_Nd_fitted_to_f8_only':
            value = brentq(lambda n: focal(32.62, 242.95, 4.065, n) - 37.19, 1.77, 2.5)
            trial[-2]['glass']['nd'] = value
            changed = [{'parameter': 'L8_Nd', 'source': 1.77, 'hypothesis': value,
                        'scope': 'Fits only Table8 f8; not evidence of the original glass'}]
        l8 = focal(trial[-2]['radius'], trial[-1]['radius'], trial[-2]['thickness'], trial[-2]['glass']['nd'])
        scenarios.append({'scenario': label, 'changes': changed, 'paraxial': paraxial(trial),
                          'L8_focal_mm': l8, 'onaxis': spot(trial, calibration)})

    configurations = []
    for config in CONFIGS:
        trial = copy.deepcopy(rows)
        trial[11]['thickness'] = config['gap_mm']
        wd = math.inf if config['object_mm'] == 'Infinity' else config['object_mm']
        configurations.append({'configuration': config, 'paraxial': paraxial(trial, wd),
                               'track_to_source_image_mm': sum(r['thickness'] for r in trial)})

    spectral = []
    for wave in ['F', 'd', 'C', '650nm', '850nm']:
        trial = copy.deepcopy(rows)
        used = apply_glasses(trial, domestic, wave)
        spectral.append({'wavelength': wave, 'glass_indices': used, 'paraxial': paraxial(trial)})

    bounds = []
    for signs in itertools.product([-1., 1.], repeat=4):
        r1, r2, t, n = [value + sign * half for value, sign, half in
                        zip([32.62, 242.95, 4.065, 1.77], signs, [.005, .005, .0005, .005])]
        bounds.append(focal(r1, r2, t, n))
    after_hash = hashlib.sha256(original_lens.read_bytes()).hexdigest()
    assert before_hash == after_hash
    result = {'status': 'Diagnostic complete; faithful patent reproduction is not established',
              'source': 'CN110703413B, source PDF Tables7/8/9 and Figures9/14/17',
              'methods': ['Reduced-angle ABCD at d-line', 'Thick-lens lensmaker formula',
                          'Independent on-axis Snell trace, calibrated and validated against saved native Zemax half-pupil ray'],
              'api_or_live_lens_mutations': False,
              'saved_lens_sha256_before': before_hash, 'saved_lens_sha256_after': after_hash,
              'native_ray_validation_absolute_error_mm': abs(native_error),
              'native_launch_over_paraxial_launch_calibration': calibration,
              'baseline': base, 'scenarios': scenarios, 'source_configurations': configurations,
              'catalog_spectral_paraxial_check': spectral,
              'L8_print_precision_corner_check': {'corners': len(bounds), 'min_focal_mm': min(bounds),
                                                  'max_focal_mm': max(bounds), 'published_focal_mm': 37.19,
                                                  'scope': 'Nearest-rounding endpoints for R1/R2/thickness/Nd; not a continuous interval proof'},
              'other_examples_cross_check': source_controls(),
              'figure_conditions': {'Figure9': {'example': 1, 'plot': 'Polychromatic diffraction MTF versus spatial frequency',
                                               'field_degrees_read_from_original_image': [0., 3., 6., 8., 10., 10.8],
                                               'max_frequency_cycles_per_mm': 200., 'nominal_f_number_from_Table2': 1.82},
                                    'Figure14': {'example': 3, 'plot': 'OTF modulus versus focus shift',
                                                 'field_degrees_read_from_original_image': [0., 6., 8., 10., 11., 13.5],
                                                 'focus_shift_plot_mm': [-.05, .05],
                                                 'frequency_and_wavelength_weights': 'Not present in the disclosed image; cannot infer an exact reproduction setting'},
                                    'current_import_fields': {'type': 'ParaxialImageHeight_mm', 'values': [0., 4.5, 6.3, 9.],
                                                              'basis': 'User Phi18 requirement; not copied from patent Figure14'},
                                    'initial_image_distance_mm': 14.65, 'current_UI_image_distance_mm': 7.634},
              'remaining_unknowns': ['Corrected original prescription matching Tables8 and9',
                                     'Exact glass catalog names or full spectral index/partial dispersion data',
                                     'Exact axial stop position and clear apertures',
                                     'WD reference plane', 'Spectral lines and weights; Figure14 evaluation frequency',
                                     'Full temperature and mechanical compensation data'],
              'limitation': 'Geometry and paraxial contradiction audit with single-color on-axis estimates. No newly computed Zemax polychromatic MTF or full-field verification.'}
    (OUT / 'reproduction_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'scenarios': [{'scenario': s['scenario'], 'efl': s['paraxial']['efl_mm'],
                                    'L8_focal': s['L8_focal_mm'],
                                    'best_geometric_focus': s['onaxis']['geometric_best_focus_from_last_vertex_mm'],
                                    'best_RMS_um': s['onaxis']['rms_at_geometric_best_focus_um']} for s in scenarios],
                      'native_ray_error_mm': abs(native_error),
                      'configurations': configurations,
                      'L8_rounding': result['L8_print_precision_corner_check'],
                      'other_examples': result['other_examples_cross_check'],
                      'lens_unchanged': before_hash == after_hash}, indent=2))


if __name__ == '__main__':
    main()
