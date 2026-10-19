"""Independent dense sag/spacing check of an exported native candidate."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent


def run(tag):
    schema = json.loads((R / (tag + '_schema.json')).read_text())
    configs = json.loads((R / (tag + '_validation_configs.json')).read_text())
    aspheres = {q['surface']-1: q['coefficients'] for q in schema['aspheres']}
    caps = [float(schema['fixed_semi_diameters'].get(str(i+1), schema['fixed_stop_radius']))
            for i in range(15)]

    def sag(surface, radius):
        c = schema['curvatures'][surface]
        discriminant = 1-c*c*radius*radius
        if discriminant < 0:
            raise ValueError(f'Invalid spherical base at surface {surface+1}, radius {radius}')
        return c*radius*radius/(1+math.sqrt(discriminant)) + sum(
            coefficient*radius**power for coefficient, power in
            zip(aspheres.get(surface, [0., 0., 0.]), [4, 6, 8]))

    cases = []
    for cfg in configs:
        t = schema['thicknesses'].copy()
        for surface, key in [(5, 'iris_gap'), (10, 'gap'), (14, 'backfocus')]:
            t[surface] = cfg[key]
        for surface, value in cfg.get('extra_thicknesses', {}).items():
            t[int(surface)-1] = value
        gaps = []
        for s in range(14):
            glass = schema['media'][s] >= 0
            radius = caps[s] if s == 4 or glass else caps[6] if s == 5 else min(caps[s:s+2])
            radii = np.linspace(0., radius, 501)
            thickness = [t[s]+sag(s+1, y)-sag(s, y) for y in radii]
            index = int(np.argmin(thickness))
            required = .45 if glass else .30 if s in (4, 5) else .10
            gaps.append(dict(surface=s+1, material='glass' if glass else 'air',
                radius_checked_mm=radius, minimum_at_radius_mm=float(radii[index]),
                minimum_mm=thickness[index], required_mm=required,
                margin_mm=thickness[index]-required))
        matrix = np.eye(2)
        n_before = 1.
        for s, curvature in enumerate(schema['curvatures']):
            medium = schema['media'][s]
            n_after = schema['glasses'][medium]['indices']['d'] if medium >= 0 else 1.
            matrix = np.array([[1., 0.], [-(n_after-n_before)*curvature, 1.]]) @ matrix
            if s < 14:
                matrix = np.array([[1., t[s]/n_after], [0., 1.]]) @ matrix
            n_before = n_after
        cases.append(dict(wd=cfg['wd'], spacing_passed=all(q['margin_mm'] >= -1e-9 for q in gaps),
            worst=min(gaps, key=lambda q:q['margin_mm']), gaps=gaps,
            independent_paraxial_efl_mm=float(-1/matrix[1, 0]),
            independent_paraxial_bfl_mm=float(-matrix[0, 0]/matrix[1, 0]),
            group_positions_from_sensor_mm=[sum(t[start-1:]) for start in schema.get('focus_group_starts', [1, 7, 12])],
            geometric_back_clearance_mm=t[-1], ttl_mm=sum(t)))
    strokes = np.ptp(np.array([q['group_positions_from_sensor_mm'] for q in cases]), axis=0).tolist()
    result = dict(tag=tag, calculation='Independent geometry; not a native ray-trace or tolerance analysis',
        radial_samples_per_gap=501, maximum_lens_diameter_mm=2*max(caps),
        focus_group_starts=schema.get('focus_group_starts', [1, 7, 12]),
        group_strokes_mm=strokes, stroke_limit_mm=6.0,
        cases=cases, spacing_passed=all(q['spacing_passed'] for q in cases),
        length_passed=max(q['ttl_mm'] for q in cases) <= 60.,
        lens_diameter_passed=2*max(caps) <= 23.,
        geometric_back_clearance_passed=min(q['geometric_back_clearance_mm'] for q in cases) >= 8.,
        paraxial_bfl_passed=min(q['independent_paraxial_bfl_mm'] for q in cases) >= 8.,
        stroke_passed=max(strokes) <= 6.0+1e-6,
        complete_requirement_acceptance=False)
    (R / (tag + '_dense_geometry_audit.json')).write_text(json.dumps(result, indent=2), encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k != 'cases'}, indent=2))
    for q in cases:
        print(json.dumps({k:v for k,v in q.items() if k != 'gaps'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    run(parser.parse_args().tag)
