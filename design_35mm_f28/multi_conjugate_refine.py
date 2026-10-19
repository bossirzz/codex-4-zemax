"""Pure numerical multi-distance refinement; native Zemax is the acceptance authority.

The declared 500/200 mm MTF thresholds are hard numerical constraints with a
calibration margin. Other distances use engineering objectives, not newly
invented user acceptance criteria. No ZOS-API connection is made here.
"""
import argparse
import copy
import json
import sys
import time
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.optimize._numdiff import approx_derivative
from scipy.spatial import ConvexHull

# The legacy ray-engine import writes a shared glass JSON. This worker must
# keep its outputs isolated, so suppress only that known import side effect.
_write_text = Path.write_text


def _isolated_import_write(path, data, *args, **kwargs):
    if path.name == 'design_glasses.json':
        return len(data)
    return _write_text(path, data, *args, **kwargs)


with patch.object(Path, 'write_text', _isolated_import_write):
    import optical_engine as e
    import fast_otf_optimize as f
    from joint_design_model import JointModel, R

sys.stdout.reconfigure(encoding='utf8')
ADDED_WDS = (300., 800., 2000., 1e7)
FIELDS = np.array([0., 3., 6.3, 9.08])


def infinity_ray_geometry(c, t, wd, p, caps, stop):
    """Only the ray/pupil quantities needed for infinity's declared checks."""
    radii = np.zeros(15)
    heights, areas, clipped = [], [], []
    invalid = 0
    ray_caps = caps.copy()
    ray_caps[5] = np.inf
    radius = p['efl'] / 5.6
    for field in [0., 9.08]:
        field_areas, field_clipped = [], []
        oy = field / p['magnification']
        for n in e.N:
            rr, hh, v = e.trace_bundle_stop_aiming(c, t, n, wd, oy, p['ep'], f.PUPIL, radius, stop)
            invalid += int((~v).sum())
            radii = np.maximum(radii, hh.max(0))
            blocked = ~np.all(hh <= ray_caps + 1e-7, axis=1)
            field_clipped.append(int(np.sum(v & blocked)))
            good = v & ~blocked
            if good.sum() < 3:
                field_areas.append(0.)
            else:
                r = rr[good]
                angular = r[:, 2:4] * r[:, 11:12]
                field_areas.append(float(ConvexHull(angular).volume))
        chief = e.aim_chief(c, t, e.N[1], wd, oy, p['ep'])
        ray, _, _ = e.ray(c, t, e.N[1], wd, oy, 0., chief, p['ep'])
        heights.append(ray[1])
        areas.append(field_areas)
        clipped.append(field_clipped)
    areas = np.asarray(areas)
    diag = dict(relative_illumination_primary=(areas[:, 1] / max(1e-16, areas[0, 1])).tolist(),
                angular_pupil_areas=areas.tolist(), clipped_ray_counts=clipped,
                phase_centroid_colour_span_um=[])
    return None, radii, np.asarray(heights), invalid, diag


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    return value


class MultiModel:
    def __init__(self, source='candidate_v29', focus_file=None, focus_import_wds=None):
        self.joint = JointModel(source)
        self.source = source
        c, t, near, _, _ = self.joint.unpack(self.joint.x)
        native = json.loads((R / (source + '_native_acceptance.json')).read_text(encoding='utf8'))
        poses = {float(row['wd']): row for row in native}
        if focus_file:
            path = Path(focus_file)
            if not path.is_absolute():
                path = R / path
            extra = json.loads(path.read_text(encoding='utf8'))
            if isinstance(extra, dict):
                extra = extra.get('cases', extra.get('results', []))
            for row in extra:
                if (all(key in row for key in ['wd', 'gap', 'iris_gap', 'backfocus']) and
                        (focus_import_wds is None or float(row['wd']) in focus_import_wds)):
                    poses[float(row['wd'])] = row
        extra = []
        for wd in ADDED_WDS:
            if wd in poses:
                row = poses[wd]
                tt = t.copy()
                tt[-1], tt[5], tt[10] = row['backfocus'], row['iris_gap'], row['gap']
            else:
                # Interpolate in object vergence, then restore paraxial focus.
                fraction = (1 / wd - 1 / 500.) / (1 / 200. - 1 / 500.)
                tt = t.copy()
                tt[-1] += fraction * (near[-1] - t[-1])
                tt[5] += fraction * (near[5] - t[5])
                defocus = t[10] - e.focus_group_gap(c, t, 500.)
                tt[10] = e.focus_group_gap(c, tt, wd) + defocus
            extra.extend([tt[-1] - t[-1], tt[5] - t[5],
                          (tt[10] - e.focus_group_gap(c, tt, wd)) / .01])
        self.seed = np.r_[self.joint.x, extra]
        self.evaluations = 0
        self.cache_key = None
        self.cache_value = None
        self.optics_cache = {}

    def unpack(self, x):
        c, t, near, caps, _ = self.joint.unpack(x[:57])
        poses = [(500., t), (200., near)]
        for j, wd in enumerate(ADDED_WDS):
            delta = x[57 + 3*j:60 + 3*j]
            tt = t.copy()
            tt[-1] += delta[0]
            tt[5] += delta[1]
            tt[10] = e.focus_group_gap(c, tt, wd) + delta[2] * .01
            poses.append((wd, tt))
        return c, poses, caps

    def evaluate(self, x, args):
        key = np.asarray(x).tobytes()
        if key == self.cache_key:
            return self.cache_value
        c, poses, caps = self.unpack(x)
        constraints, labels, cases, positions = [], [], [], []
        engineering_residuals = []
        engineering_objective = 0.

        def add(name, values):
            values = np.atleast_1d(values).ravel()
            constraints.extend(values.tolist())
            labels.extend([name if len(values) == 1 else f'{name}[{i}]'
                           for i in range(len(values))])

        add('curvature_valid', 1 - c[:15]**2 * caps**2 - .01)
        add('max_lens_diameter', (11.5 - caps) * 10)
        pb = e.paraxial(c, poses[0][1], 500.)
        add('best_efl', (.045 - abs(pb['efl'] - 35.)) * 5)
        for wd, t in poses:
            p = e.paraxial(c, t, wd)
            marginal, _, marginal_ok = e.ray(c, t, e.N[1], wd, 0., 0., p['efl']/5.6, p['ep'])
            stop = abs(marginal[9])
            if wd == 200.:
                fc, fe = 150., 120.
            elif wd == 2000.:
                fc, fe = args.far_center_frequency, args.far_edge_frequency
            else:
                fc, fe = 225., 180.
            is_infinity = wd >= 1e7
            low_frequency = wd == 2000. and args.far_90_target is not None
            optical_key = (wd, fc, fe, low_frequency, c.tobytes(), t.tobytes(), caps.tobytes())
            optical_result = self.optics_cache.get(optical_key)
            if is_infinity:
                # Infinity has CRA/illumination requirements but no specified
                # high-frequency MTF target. Avoid spending six OTF fits here.
                if optical_result is None:
                    optical_result = infinity_ray_geometry(c, t, wd, p, caps, stop)
                _, radii, heights, invalid, diag = optical_result
                mtf = None
                field_values = np.array([0., 9.08])
            else:
                if optical_result is None:
                    optical_result = f.quantities(
                        c, t, wd=wd, fields=FIELDS, frequencies=(fc, fe),
                        fixed_stop_radius=stop, low_frequency=low_frequency, caps=caps,
                        return_diagnostics=True)
                mtf, radii, heights, invalid, diag = optical_result
                mtf90 = mtf[:, 2:].copy() if low_frequency else None
                mtf = mtf[:, :2]
                field_values = FIELDS
            self.optics_cache[optical_key] = optical_result
            # A full shape Jacobian changes the glass surfaces at every
            # perturbation. A bounded cache still reuses unaffected distances
            # during pose derivatives without accumulating a long-run archive.
            if len(self.optics_cache) > 768:
                for stale_key in list(self.optics_cache)[:256]:
                    del self.optics_cache[stale_key]
            prefix = f'wd{wd:g}'
            # blank_clearance subtracts .5 mm glass/.1 mm air/.3 mm iris.
            # The .05 addition restores the accepted .45 mm glass minimum.
            glass_adjust = np.where(e.MEDIA[:14] >= 0, .05, 0.)
            clearance = f.blank_clearance(c, t, caps) + glass_adjust
            material = e.MEDIA[:14] >= 0
            if wd == 500.:
                add('shared_material_clearance', 20 * clearance[material])
                safety_residual = 10 * np.maximum(.003 - clearance[material], 0.)
                engineering_residuals.append(safety_residual)
                engineering_objective += float(safety_residual @ safety_residual)
            add(prefix + '_air_clearance', 20 * clearance[~material])
            # Vignetting is modelled with the same physical lens radii; it is
            # allowed only through the explicit illumination constraint.
            add(prefix + '_ttl', (59.8 - np.sum(t)) * .5)
            add(prefix + '_sensor_clearance', (t[-1] - 8.) * .5)
            add(prefix + '_image_circle', (heights[-1] - 9.) * 20)
            add(prefix + '_invalid_rays', -.1 * (invalid + int(not marginal_ok)))
            ri = float(diag['relative_illumination_primary'][-1])
            # The angular-pupil proxy is biased versus native RELI. It is an
            # engineering penalty only, never a native illumination verdict.
            engineering_objective += .0001 * max(args.ri_floor - ri, 0.)**2
            engineering_residuals.append(np.array([.01 * max(args.ri_floor - ri, 0.)]))
            if wd in (500., 200.):
                add(prefix + '_mtf_center', 10 * (mtf[0] - args.primary_center))
                add(prefix + '_mtf_offaxis', 10 * (mtf[1:] - args.primary_edge))
                distortion_limit = .0029 if wd == 500. else .0058
                add(prefix + '_distortion', 20 * (field_values[1:] * distortion_limit -
                                                abs(heights[1:] - field_values[1:])))
                if args.colour_hard:
                    add(prefix + '_colour', .5 * (args.colour_limit -
                        np.asarray(diag['phase_centroid_colour_span_um'][1:])))
            elif not is_infinity:
                target_center = args.far_center if wd == 2000. else args.engineering_center
                target_edge = args.far_edge if wd == 2000. else args.engineering_edge
                desired = np.r_[np.full((1, 2), target_center),
                                np.full((3, 2), target_edge)]
                deficit = np.minimum(mtf - desired, 0.)
                weight = args.far_weight if wd == 2000. else args.mid_weight
                engineering_objective += weight * float(np.mean(deficit**2))
                engineering_residuals.append(np.sqrt(weight / deficit.size) * deficit.ravel())
                if low_frequency:
                    deficit90 = np.minimum(mtf90[1:] - args.far_90_target, 0.)
                    engineering_objective += args.far_90_weight * float(np.mean(deficit90**2))
                    engineering_residuals.append(np.sqrt(args.far_90_weight / deficit90.size) * deficit90.ravel())
            if not is_infinity:
                # The previous 2.4 um colour proxy was an engineering choice,
                # so keep it soft while retaining actual spans in all reports.
                colour = np.asarray(diag['phase_centroid_colour_span_um'][1:])
                engineering_objective += args.colour_weight * float(
                    np.mean(np.maximum(colour - args.colour_limit, 0.)**2))
                engineering_residuals.append(np.sqrt(args.colour_weight / len(colour)) *
                    np.maximum(colour - args.colour_limit, 0.))
            positions.append([t[5] + t[10] + t[-1], t[10] + t[-1], t[-1]])
            case = dict(wd=wd, gap=t[10], iris_gap=t[5], backfocus=t[-1],
                        fixed_stop_radius=stop, nominal_fno_target=2.8,
                        frequency_center=fc, frequency_edge=fe, predicted_mtf=mtf,
                        efl=p['efl'], ttl=np.sum(t), chief_image_heights=heights,
                        relative_illumination_approx=ri,
                        phase_centroid_colour_span_um=diag['phase_centroid_colour_span_um'],
                        clipped_ray_counts=diag['clipped_ray_counts'],
                        max_ray_radius=radii,
                        acceptance_basis='user thresholds' if wd in (500., 200.) else 'engineering check only')
            if is_infinity:
                oy = 9.08 / p['magnification']
                chief = e.aim_chief(c, t, e.N[1], wd, oy, p['ep'])
                ray, _, ok = e.ray(c, t, e.N[1], wd, oy, 0., chief, p['ep'])
                cra = abs(np.degrees(np.arctan(ray[3])))
                add('infinity_cra', args.cra_limit - cra)
                add('infinity_chief_valid', -float(not ok))
                case['predicted_cra_deg'] = cra
            elif low_frequency:
                case['engineering_90_mtf'] = mtf90
            cases.append(case)
        strokes = np.ptp(positions, axis=0)
        add('group_strokes', (args.stroke_limit - strokes) * .5)
        result = dict(constraints=np.asarray(constraints), labels=labels,
                      cases=cases, group_strokes=strokes,
                      engineering_objective=engineering_objective,
                      engineering_residuals=np.concatenate(engineering_residuals))
        self.cache_key, self.cache_value = key, result
        self.evaluations += 1
        return result

    def export(self, x, args, report):
        self.joint.export(x[:57], args.tag)
        # Retain every distance in an explicit validation file, while the
        # legacy focus_targets file remains compatible with two-distance tools.
        rows = []
        for case in report['cases']:
            row = {k: case[k] for k in ['wd', 'gap', 'iris_gap', 'backfocus',
                'fixed_stop_radius', 'nominal_fno_target', 'frequency_center', 'frequency_edge']}
            row['note'] = 'Pure numerical pose; requires native validation; ' + case['acceptance_basis']
            rows.append(row)
        (R / (args.tag + '_validation_configs.json')).write_text(
            json.dumps(serializable(rows), ensure_ascii=False, indent=2), encoding='utf8')
        np.save(R / (args.tag + '_multi_parameters.npy'), x)
        self.write_report(args.tag + '_prediction.json', report)

    @staticmethod
    def write_report(filename, report):
        report = copy.deepcopy(report)
        constraints = np.asarray(report.pop('constraints'))
        labels = report.pop('labels')
        bad = np.argsort(constraints)[:12]
        report['minimum_constraint'] = float(constraints.min())
        report['limiting_constraints'] = [dict(name=labels[i], value=float(constraints[i])) for i in bad]
        report['native_validation_required'] = True
        (R / filename).write_text(json.dumps(serializable(report), ensure_ascii=False, indent=2), encoding='utf8')


def run(args):
    if not args.tag.startswith('candidate_v30_mc'):
        raise ValueError('This isolated worker must write candidate_v30_mc* outputs.')
    model = MultiModel(args.source, args.focus_file, args.focus_import_wds)
    seed = model.seed
    if args.seed_file:
        path = Path(args.seed_file)
        seed = np.load(path if path.is_absolute() else R / path)
        if seed.shape != model.seed.shape:
            raise ValueError(f'Expected {len(model.seed)} parameters, got {seed.shape}.')
    baseline = model.evaluate(seed, args)
    model.write_report(args.tag + '_baseline.json', baseline)
    print('BASELINE', json.dumps(dict(min_constraint=float(baseline['constraints'].min()),
        engineering_objective=baseline['engineering_objective'], group_strokes=baseline['group_strokes'].tolist(),
        mtfs={str(q['wd']): serializable(q['predicted_mtf']) for q in baseline['cases']})), flush=True)
    if args.evaluate_only:
        model.export(seed, args, baseline)
        return
    # Open the shared shape only after diagnosing focus. Pose-only and low
    # dimensional stages preserve the accepted local design more readily.
    pose = np.r_[[19, 24, 28, 47, 48, 49], np.arange(57, 69)]
    if args.stage == 'pose':
        active = pose
    elif args.stage == 'shape':
        active = np.r_[np.arange(14), pose]
    elif args.stage == 'aspheric':
        active = np.r_[np.arange(14), [21, 23], np.arange(29, 41), pose]
    else:
        active = np.arange(69)
    delta = np.r_[np.full(14, .10), np.full(15, .65), np.full(18, 1.5),
                  [2., 1., 10., .6], np.full(6, 1.5), np.tile([3., 2., 20.], 4)]
    lo, hi = seed - delta, seed + delta
    lo[14:29] = np.maximum(lo[14:29], [1, .1, 1, .45, .3, .3, .45, 1, .1, 1, .15, 1, .1, .45, 8])
    lo[50], hi[50] = max(lo[50], 10.05), min(hi[50], 11.45)
    if args.stage == 'aspheric':
        lo[[21, 23]], hi[[21, 23]] = seed[[21, 23]] - .1, seed[[21, 23]] + .1
    scales = np.ones(69)
    if args.local_trust:
        scales[:14] = .0004
        scales[14:29] = .02
        scales[[21, 23]] = .10
        scales[24] = .005
        scales[29:47] = .01
        scales[47:50] = [.02, .02, .5]
        scales[50] = .02
        scales[51:57] = .01
        scales[57:69] = np.tile([.02, .02, .5], 4)
        lo = np.maximum(lo, seed - args.local_radius * scales)
        hi = np.minimum(hi, seed + args.local_radius * scales)
        alo, ahi = (lo[active] - seed[active]) / scales[active], (hi[active] - seed[active]) / scales[active]
        initial = np.zeros(len(active))
    else:
        alo, ahi = lo[active], hi[active]
        initial = seed[active]
    best_feasible = None
    best_feasible_objective = np.inf
    best_penalty = np.inf
    best_penalty_x = seed.copy()
    started = time.time()
    iterations = 0
    derivative_key = None
    derivative_value = None

    def unpack(active_x):
        x = seed.copy()
        x[active] = seed[active] + active_x * scales[active] if args.local_trust else active_x
        return x

    def evaluate(active_x):
        nonlocal best_feasible, best_feasible_objective, best_penalty, best_penalty_x
        x = unpack(active_x)
        result = model.evaluate(x, args)
        objective = result['engineering_objective'] + .00002 * float(np.mean((x - seed)**2))
        penalty = objective + 1000. * float(np.sum(np.minimum(result['constraints'], 0.)**2))
        if penalty < best_penalty:
            best_penalty, best_penalty_x = penalty, x.copy()
            np.save(R / (args.tag + '_best_penalty.npy'), x)
        if result['constraints'].min() >= -1e-5 and objective < best_feasible_objective:
            best_feasible, best_feasible_objective = x.copy(), objective
            np.save(R / (args.tag + '_best_feasible.npy'), x)
        return objective, result['constraints']

    def callback(active_x):
        nonlocal iterations
        iterations += 1
        objective, constraints = evaluate(active_x)
        print(json.dumps(dict(iteration=iterations, evaluations=model.evaluations,
            seconds=time.time() - started, objective=objective,
            minimum_constraint=float(constraints.min()),
            best_feasible_objective=None if best_feasible is None else best_feasible_objective)), flush=True)
        if best_feasible is not None:
            checkpoint_args = copy.copy(args)
            checkpoint_args.tag = args.tag + '_checkpoint'
            checkpoint = model.evaluate(best_feasible, args)
            model.export(best_feasible, checkpoint_args, checkpoint)
            print('FEASIBLE_CHECKPOINT', checkpoint_args.tag,
                  json.dumps({str(q['wd']): serializable(q['predicted_mtf'])
                              for q in checkpoint['cases'] if q['wd'] in (500., 200., 2000.)}), flush=True)

    def derivative(active_x):
        # SLSQP requests objective and constraint derivatives at the same
        # point. Differentiate them together so each expensive optical
        # perturbation is evaluated once rather than twice.
        nonlocal derivative_key, derivative_value
        key = np.asarray(active_x).tobytes()
        if key != derivative_key:
            def joined(y):
                objective, constraints = evaluate(y)
                return np.r_[objective, constraints]
            derivative_value = approx_derivative(joined, active_x, method='3-point',
                abs_step=(.005 if args.local_trust else 1e-5)*(1+abs(active_x)), bounds=(alo, ahi))
            derivative_key = key
        return derivative_value

    if args.solver == 'least_squares':
        last_log_call = 0

        def residual(active_x):
            nonlocal last_log_call
            evaluate(active_x)
            x = unpack(active_x)
            result = model.evaluate(x, args)
            con = result['constraints']
            if model.evaluations - last_log_call >= 150:
                callback(active_x)
                last_log_call = model.evaluations
            return np.r_[args.constraint_weight * np.minimum(con, 0.),
                         result['engineering_residuals'],
                         np.sqrt(.00002 / len(x)) * (x - seed)]

        opt = least_squares(residual, initial, bounds=(alo, ahi),
            jac=lambda y: approx_derivative(residual, y, method='3-point',
                abs_step=(.005 if args.local_trust else 1e-5)*(1+abs(y)), bounds=(alo, ahi)),
            x_scale='jac', max_nfev=args.steps, ftol=1e-9, xtol=1e-10,
            gtol=1e-7, verbose=1)
    else:
        opt = minimize(lambda y: evaluate(y)[0], initial, method='SLSQP',
            jac=lambda y: derivative(y)[0],
            constraints=[dict(type='ineq', fun=lambda y: evaluate(y)[1],
                jac=lambda y: derivative(y)[1:])],
            bounds=list(zip(alo, ahi)), callback=callback,
            options=dict(maxiter=args.steps, ftol=1e-8, disp=True))
    final = best_feasible if best_feasible is not None else best_penalty_x
    report = model.evaluate(final, args)
    report.update(solver_status=int(opt.status), solver_message=str(opt.message),
                  stage=args.stage, elapsed_seconds=time.time()-started,
                  numerical_constraints_feasible=bool(report['constraints'].min() >= -1e-5))
    model.export(final, args, report)
    print('EXPORTED', args.tag, 'PREDICTION_ONLY_REQUIRES_NATIVE_VALIDATION', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default='candidate_v29')
    parser.add_argument('--tag', default='candidate_v30_mc')
    parser.add_argument('--focus-file', help='Fixed-lens search poses to initialize 300/800/2000/infinity.')
    parser.add_argument('--focus-import-wds', type=float, nargs='+')
    parser.add_argument('--seed-file')
    parser.add_argument('--stage', choices=['pose', 'shape', 'aspheric', 'full'], default='shape')
    parser.add_argument('--steps', type=int, default=25)
    parser.add_argument('--solver', choices=['slsqp', 'least_squares'], default='slsqp')
    parser.add_argument('--constraint-weight', type=float, default=50.)
    parser.add_argument('--local-trust', action='store_true')
    parser.add_argument('--local-radius', type=float, default=1.)
    parser.add_argument('--evaluate-only', action='store_true')
    parser.add_argument('--primary-center', type=float, default=.307)
    parser.add_argument('--primary-edge', type=float, default=.208)
    parser.add_argument('--engineering-center', type=float, default=.28)
    parser.add_argument('--engineering-edge', type=float, default=.18)
    parser.add_argument('--far-center', type=float, default=.25)
    parser.add_argument('--far-edge', type=float, default=.10)
    parser.add_argument('--far-center-frequency', type=float, default=225.)
    parser.add_argument('--far-edge-frequency', type=float, default=180.)
    parser.add_argument('--far-weight', type=float, default=3.)
    parser.add_argument('--far-90-target', type=float)
    parser.add_argument('--far-90-weight', type=float, default=1.5)
    parser.add_argument('--mid-weight', type=float, default=1.25)
    parser.add_argument('--colour-limit', type=float, default=2.45,
                        help='PSF-centroid colour proxy in um; not a perceptual colour-fringe certificate.')
    parser.add_argument('--colour-hard', action='store_true')
    parser.add_argument('--colour-weight', type=float, default=.00001)
    parser.add_argument('--ri-floor', type=float, default=.78)
    parser.add_argument('--cra-limit', type=float, default=9.99)
    parser.add_argument('--stroke-limit', type=float, default=6.01)
    run(parser.parse_args())
