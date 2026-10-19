"""Independent monochromatic Snell trace, validated against the saved Zemax ray."""
import json
import math
from pathlib import Path

import numpy as np

from import_patent_zospy import compile_prescription

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / 'zemax_import' / 'CN110703413B_example3'
ROWS = compile_prescription()
saved = json.loads((MODEL / 'ray_trace.json').read_text(encoding='utf-8'))['data']['real_ray_trace_data']
object_ray = dict(zip(saved['columns'], saved['data'][0]))
image_ray = dict(zip(saved['columns'], saved['data'][-1]))
tan_full = 2 * object_ray['Y-cosine'] / object_ray['Z-cosine']


def trace(pupil):
    p = np.array([0., 0., -500.])
    d = np.array([0., tan_full * pupil, 1.])
    d /= np.linalg.norm(d)
    vertex = 0.
    index = 1.
    for row in ROWS:
        next_index = 1. if row['glass'] is None else row['glass']['nd']
        radius = math.inf if row['radius'] == 'Infinity' else row['radius']
        if math.isinf(radius):
            hit = p + (vertex - p[2]) / d[2] * d
            normal = np.array([0., 0., -1.])
        else:
            center = np.array([0., 0., vertex + radius])
            offset = p - center
            b = np.dot(d, offset)
            discriminant = b * b - np.dot(offset, offset) + radius * radius
            if discriminant < 0:
                raise ValueError('Missed spherical surface')
            root = math.sqrt(discriminant)
            distance = -b - root if radius > 0 else -b + root
            if distance < -1e-6:
                raise ValueError('Non-forward intersection')
            hit = p + distance * d
            normal = (hit - center) / abs(radius)
        if np.dot(normal, d) > 0:
            normal = -normal
        cosine = -np.dot(normal, d)
        eta = index / next_index
        k = 1 - eta * eta * (1 - cosine * cosine)
        if k < 0:
            raise ValueError('Total internal reflection')
        d = eta * d + (eta * cosine - math.sqrt(k)) * normal
        d /= np.linalg.norm(d)
        p = hit
        index = next_index
        if row is not ROWS[-1]:
            vertex += row['thickness']
    slope = d[1] / d[2]
    # Express the emergent ray at the plane through the last surface vertex.
    y_at_vertex_plane = p[1] + (vertex - p[2]) * slope
    return y_at_vertex_plane, slope


half_a, half_b = trace(.5)
half_image_y = half_a + 14.65 * half_b
error = half_image_y - image_ray['Y-coordinate']
assert abs(error) < 1e-4, (half_image_y, image_ray['Y-coordinate'], error)
r = (np.arange(1000) + .5) / 1000
rays = np.array([trace(float(value)) for value in r])
weights = 2 * r
weights /= weights.sum()
a, b = rays.T
best = -np.sum(weights * a * b) / np.sum(weights * b * b)


def rms(distance):
    return float(np.sqrt(np.sum(weights * (a + distance * b) ** 2)))


result = {'method': 'Independent on-axis spherical refraction at published Nd; pupil launch validated against native Zemax',
          'configuration': 'WD500', 'native_half_pupil_image_y_mm': image_ray['Y-coordinate'],
          'independent_half_pupil_image_y_mm': half_image_y,
          'validation_absolute_error_mm': abs(error),
          'scope': 'Monochromatic on-axis estimate, ignoring clipping; not a polychromatic MTF or full design verification',
          'geometric_best_focus_last_surface_distance_mm': float(best),
          'rms_spot_radius_at_geometric_best_focus_um': rms(best) * 1000,
          'rms_spot_radius_at_current_ui_7_634mm_um': rms(7.634) * 1000,
          'rms_spot_radius_at_original_14_65mm_um': rms(14.65) * 1000,
          'ray_axis_intersections_last_surface_distance_mm': [
              {'normalized_pupil_radius': value, 'focus_distance_mm': -trace(value)[0] / trace(value)[1]}
              for value in [.05, .25, .5, .75, 1.]]}
path = MODEL / 'mtf_diagnosis' / 'onaxis_aberration_check.json'
path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=True, indent=2), flush=True)
