"""Native OpticStudio Contrast/DLS optimization on an isolated v31 copy.

Python prepares the declared configurations and monitors the tool. All search
steps and objective evaluations are performed by OpticStudio's optimizer.
Outputs remain experimental until separately accepted at high sampling.
"""
import os
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time

import zospy as zp

R = Path(__file__).resolve().parent
SOURCE = 'candidate_v31_thermal_probe'
TAG = 'candidate_v32_native_contrast'
sys.stdout.reconfigure(encoding='utf8')


def write(name, value):
    (R / (TAG + name + '.json')).write_text(json.dumps(value, indent=2), encoding='utf8')


def run(cycles=5, seconds=600, prepare_only=False, guarded=False):
    baseline = json.loads((R / (SOURCE + '_native_acceptance.json')).read_text())
    schema = json.loads((R / (SOURCE + '_schema.json')).read_text())
    old = json.loads((R / 'candidate_v29_native_acceptance.json').read_text())
    floating = json.loads((R / 'far_focus_diagnostic_v29_results.json').read_text())
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = opt = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(R / (SOURCE + '_review.zos'), saveifneeded=False)
        typ = zp.constants.Editors.MFE.MeritOperandType

        def val(kind, *args):
            pp = list(args) + [0.] * (8 - len(args))
            return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

        o.MCE.SetCurrentConfiguration(1)
        materials = []
        for surface, glass in zip([1, 3, 4, 7, 8, 10, 12, 14], schema['glasses']):
            actual = [val('INDX', surface, wave) for wave in [1, 2, 3]]
            expected = [glass['indices'][wave] for wave in ['F', 'd', 'C']]
            if max(abs(a - b) for a, b in zip(actual, expected)) > 1e-10:
                raise RuntimeError(f'Glass resolution differs at surface {surface}')
            materials.append(dict(surface=surface, name=str(o.LDE.GetSurfaceAt(surface).Material), indices=actual))

        m = o.MCE
        controllers = {(str(m.GetOperandAt(i).Type), int(m.GetOperandAt(i).Param1)): m.GetOperandAt(i)
                       for i in range(1, m.NumberOfOperands + 1)}
        best = next(q for q in baseline if q['wd'] == 500)
        near = next(q for q in baseline if q['wd'] == 200)
        cfgs = [dict(q, configuration=i) for i, q in enumerate(baseline, 1)]
        # IEEE infinity passed through this API produces invalid field rays.
        # The finite API value 1e10 is recognized as native Infinity and was
        # independently checked against the 10 km limit before this run.
        for cfg in cfgs:
            if cfg['wd'] >= 1e7:
                controllers[('THIC', 0)].GetOperandCell(cfg['configuration']).DoubleValue = 1e10 if guarded else 1e7
        q400 = dict(wd=400.)
        f400 = (1 / 400 - 1 / 500) / (1 / 200 - 1 / 500)
        for key in ['gap', 'iris_gap', 'backfocus', 'stop_radius']:
            q400[key] = best[key] + f400 * (near[key] - best[key])
        q300 = next(q for q in floating if q['wd'] == 300)
        q2000 = next(q for q in old if q['wd'] == 2000)
        for q in [q300, q400, q2000]:
            m.AddConfiguration(False)
            index = m.NumberOfConfigurations
            stop = q.get('stop_radius', q.get('fixed_stop_radius', best['stop_radius']))
            values = {('THIC', 0): q['wd'], ('THIC', 6): q['iris_gap'],
                      ('THIC', 11): q['gap'], ('THIC', 15): q['backfocus'], ('SDIA', 6): stop}
            for key, value in values.items():
                controllers[key].GetOperandCell(index).DoubleValue = float(value)
            m.SetCurrentConfiguration(index)
            for _ in range(3):
                cell = controllers[('SDIA', 6)].GetOperandCell(index)
                cell.DoubleValue *= val('ISFN') / 2.8
            cfgs.append(dict(wd=q['wd'], configuration=index, gap=q['gap'],
                             iris_gap=q['iris_gap'], backfocus=q['backfocus']))
        write('_contract', dict(source=SOURCE + '_review.zos', glasses=materials,
              configurations=[dict(configuration=q['configuration'], wd=q['wd']) for q in cfgs],
              optimizer='OpticStudio DampedLeastSquares', objective='Native Contrast MECS/MECT plus native constraints',
              main_user_thresholds='500mm:225/180 at .30/.20; 200mm:150/120 at .30/.20',
              other_distances='Engineering objectives; no new user acceptance thresholds',
              infinity_search_proxy_mm=None if guarded else 1e7,
              infinity_api_value=1e10 if guarded else None, true_infinity_acceptance_required=True,
              guard_profile='strong MTF and radial thickness guards' if guarded else 'initial contrast trial',
              temperature_cases_in_objective=False, high_sampling_acceptance_required=True))

        variables = []
        for surface in range(1, 16):
            if surface != 6:
                o.LDE.GetSurfaceAt(surface).RadiusCell.MakeSolveVariable()
                variables.append(dict(editor='LDE', surface=surface, parameter='Radius'))
        for a in schema['aspheres']:
            for par in [2, 3, 4]:
                o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                    getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{par}')).MakeSolveVariable()
                variables.append(dict(editor='LDE', surface=a['surface'], parameter=f'Par{par}'))
        for surface in [8, 10]:
            o.LDE.GetSurfaceAt(surface).ThicknessCell.MakeSolveVariable()
            variables.append(dict(editor='LDE', surface=surface, parameter='Thickness'))
        for cfg in cfgs:
            for kind, surface in [('THIC', 6), ('THIC', 11), ('THIC', 15), ('SDIA', 6)]:
                controllers[(kind, surface)].GetOperandCell(cfg['configuration']).MakeSolveVariable()
                variables.append(dict(editor='MCE', configuration=cfg['configuration'], type=kind, surface=surface))

        mf = o.MFE
        while mf.NumberOfOperands > 1:
            mf.RemoveOperandAt(mf.NumberOfOperands)
        first = True
        blocks = []
        for cfg in cfgs:
            wd, index = cfg['wd'], cfg['configuration']
            if wd >= 1e7:
                continue
            fc, fe = (150., 120.) if wd == 200 else (225., 180.)
            for centre_only, frequency in [(False, fe), (True, fc)]:
                if not first:
                    mf.AddOperand()
                w = mf.SEQOptimizationWizard2
                w.StartAt = 1 if first else mf.NumberOfOperands
                start = w.StartAt
                first = False
                w.Type = zp.constants.Wizards.OptimizationTypes.RMS
                w.Criterion = zp.constants.Wizards.CriterionTypes.Contrast
                w.Reference = zp.constants.Wizards.ReferenceTypes.Centroid
                w.UseAllConfigurations = False
                w.ConfigurationNumber = index
                w.FieldNumber = 1 if centre_only else 0
                w.UseAllFields = not centre_only
                w.UseGaussianQuadrature = True
                w.Rings = 3
                w.Arms = zp.constants.Wizards.PupilArmsCount.Arms_6
                w.SpatialFrequency = frequency
                w.IgnoreLateralColor = False
                w.DeleteVignetted = False
                w.OverallWeight = (2. if wd in (500, 200) else 1.) * (.6 if centre_only else 1.)
                w.UseGlassBoundaryValues = not centre_only
                w.GlassMin, w.GlassMax, w.GlassEdgeThickness = .45, 8., .45
                w.UseAirBoundaryValues = not centre_only
                w.AirMin, w.AirMax, w.AirEdgeThickness = .10, 25., .10
                w.UseMaximumDistortion = not centre_only
                w.MaxDistortionPct = .3 if wd in (400, 500) else .6
                w.Apply()
                blocks.append(dict(configuration=index, wd=wd, frequency=frequency,
                                   centre_only=centre_only, start=start, end=mf.NumberOfOperands))
                print('NATIVE_CONTRAST_BLOCK', json.dumps(blocks[-1]), flush=True)

        # Wizard defaults include surface zero. Exclude object and image
        # thicknesses so finite/infinite working distances are not air gaps.
        for i in range(1, mf.NumberOfOperands + 1):
            row = mf.GetOperandAt(i)
            if str(row.Type) in ['MNCA', 'MXCA', 'MNEA', 'MNCG', 'MXCG', 'MNEG']:
                row.GetCellAt(2).IntegerValue = 1
                row.GetCellAt(3).IntegerValue = 15

        explicit = []

        def operand(kind, params=(), target=0., weight=0.):
            row = mf.AddOperand()
            row.ChangeType(getattr(typ, kind))
            for i, value in enumerate(params, 2):
                cell = row.GetCellAt(i)
                if str(cell.DataType) == 'Integer':
                    cell.IntegerValue = int(value)
                else:
                    cell.DoubleValue = float(value)
            row.Target, row.Weight = float(target), float(weight)
            index = mf.NumberOfOperands
            explicit.append(dict(index=index, type=kind, parameters=list(params), target=target, weight=weight))
            return index

        def limits(raw, low=None, high=None, weight=None):
            if weight is None:
                weight = 200000. if guarded else 2000.
            if low is not None:
                operand('OPGT', [raw], low, weight)
            if high is not None:
                operand('OPLT', [raw], high, weight)

        position_rows = []
        for cfg in cfgs:
            wd, index = cfg['wd'], cfg['configuration']
            operand('CONF', [index])
            limits(operand('ISFN'), 2.7995, 2.8005)
            limits(operand('TTHI', [1, 15]), high=59.8)
            limits(operand('TTHI', [15, 15]), low=8.)
            limits(operand('TTHI', [6, 6]), low=.3)
            limits(operand('RELI', [20, 2, 4, 0]), low=.75)
            for h in [3., 6.3, 9.08]:
                raw = operand('REAY', [16, 2, 0., h / 9.08, 0., 0.])
                distortion = .003 if wd in (400, 500) else .006
                limits(raw, h * (1 - distortion), h * (1 + distortion))
                if h == 9.08:
                    limits(raw, low=9.)
            if wd == 500:
                raw = operand('EFFL', [0, 2], 35., 2.)
                limits(raw, 34.94, 35.06)
            if wd >= 1e7:
                limits(operand('RANG', [16, 2, 0., 1., 0., 0.]), high=math.radians(10.))
            if wd in (500, 200, 800) or (guarded and wd < 1e7):
                fc, fe = (150., 120.) if wd == 200 else (225., 180.)
                if guarded:
                    centre_floor, offaxis_floor, mtf_weight = (
                        (.312, .208, 30000.) if wd in (500, 200, 800) else
                        (.30, .15, 1000.) if wd in (300, 400) else (.25, .08, 1000.))
                else:
                    centre_floor, offaxis_floor, mtf_weight = .3, .2, 300.
                for field in range(1, 5):
                    for kind in ['MTFT', 'MTFS']:
                        raw = operand(kind, [4 if guarded else 3, 0, field, fc if field == 1 else fe, 1, 0])
                        limits(raw, low=centre_floor if field == 1 else offaxis_floor, weight=mtf_weight)
            if guarded:
                caps = [float(schema['fixed_semi_diameters'].get(str(i), schema['fixed_stop_radius']))
                        for i in range(1, 16)]
                for s in range(14):
                    glass = schema['media'][s] >= 0
                    radius = caps[s] if s == 4 or glass else caps[6] if s == 5 else min(caps[s:s+2])
                    required = .451 if glass else .30 if s in (4, 5) else .10
                    for fraction in [0., .25, .5, .75, 1.]:
                        raw = operand('STHI', [s+1, 0., radius*fraction, 0])
                        limits(raw, low=required)
            position_rows.append([operand('TTHI', [start, 15]) for start in [1, 7, 12]])
        # ZTHI directly bounds variation of a surface range over all
        # configurations, avoiding cached values across CONF blocks.
        for start in [1, 7, 12]:
            operand('ZTHI', [start, 15], 6.01, 200000. if guarded else 2000.)
        m.SetCurrentConfiguration(1)
        prepared = R / (TAG + '_prepared.zos')
        o.save_as(prepared)
        initial = float(mf.CalculateMeritFunction())
        audit = sorted([dict(index=i, type=str(mf.GetOperandAt(i).Type),
                      value=float(mf.GetOperandAt(i).Value), target=float(mf.GetOperandAt(i).Target),
                      weight=float(mf.GetOperandAt(i).Weight), contribution=float(mf.GetOperandAt(i).Contribution))
                      for i in range(1, mf.NumberOfOperands + 1)],
                      key=lambda q: q['contribution'] if math.isfinite(q['contribution']) else -1., reverse=True)
        write('_merit_audit', dict(initial_merit=initial, variables=variables, contrast_blocks=blocks,
              operand_counts=dict(Counter(str(mf.GetOperandAt(i).Type) for i in range(1, mf.NumberOfOperands + 1))),
              explicit_operands=explicit, largest_contributions=audit[:50]))
        print('NATIVE_PREPARED', initial, 'VARIABLES', len(variables), 'OPERANDS', mf.NumberOfOperands, flush=True)
        if not math.isfinite(initial) or initial >= 1e8:
            raise RuntimeError('Non-finite native merit function; optimization not started.')
        if prepare_only:
            return
        opt = o.Tools.OpenLocalOptimization()
        opt.Algorithm = zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
        opt.Cycles = getattr(zp.constants.Tools.Optimization.OptimizationCycles,
                             {1: 'Fixed_1_Cycle', 5: 'Fixed_5_Cycles', 10: 'Fixed_10_Cycles'}[cycles])
        opt.NumberOfCores = min(8, int(opt.MaxCores))
        start_time = datetime.now(timezone.utc).isoformat()
        info = dict(status='starting', started_utc=start_time, algorithm=str(opt.Algorithm),
                    variables=int(opt.Variables), targets=int(opt.Targets), cores=int(opt.NumberOfCores),
                    initial_merit=initial, source=SOURCE, full_requirements_passed=False)
        write('_progress', info)
        print('NATIVE_DLS_START', json.dumps(info), flush=True)
        opt.Run()
        clock = time.monotonic()
        last = clock - 10
        cancelled = False
        while opt.IsRunning:
            now = time.monotonic()
            if now - last >= 10:
                info.update(status='running', elapsed_seconds=now-clock,
                            progress=float(opt.Progress), current_merit=float(opt.CurrentMeritFunction))
                write('_progress', info)
                print('NATIVE_DLS_PROGRESS', json.dumps(info), flush=True)
                last = now
            if now - clock > seconds and not cancelled and opt.CanCancel:
                opt.Cancel()
                cancelled = True
            time.sleep(.5)
        info.update(status='cancelled_at_time_limit' if cancelled else 'finished',
                    elapsed_seconds=time.monotonic()-clock, succeeded=bool(opt.Succeeded),
                    current_merit=float(opt.CurrentMeritFunction), error=str(opt.ErrorMessage))
        opt.Close()
        opt = None
        write('_progress', info)
        print('NATIVE_DLS_DONE', json.dumps(info), flush=True)
        m.SetCurrentConfiguration(1)
        o.save_as(R / (TAG + '.zos'))
        output_schema = json.loads(json.dumps(schema))
        output_schema['curvatures'] = [0. if math.isinf(float(o.LDE.GetSurfaceAt(i).Radius)) else
                                       1 / float(o.LDE.GetSurfaceAt(i).Radius) for i in range(1, 16)]
        output_schema['thicknesses'] = [float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1, 16)]
        for a in output_schema['aspheres']:
            a['coefficients'] = [float(o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{j}')).DoubleValue) for j in [2, 3, 4]]
        output_schema['fixed_stop_radius'] = float(o.LDE.GetSurfaceAt(6).SemiDiameter)
        output_schema['native_optimizer_result_requires_independent_acceptance'] = True
        write('_schema', output_schema)
        poses = []
        for cfg in cfgs:
            m.SetCurrentConfiguration(cfg['configuration'])
            poses.append(dict(wd=cfg['wd'], gap=float(o.LDE.GetSurfaceAt(11).Thickness),
                iris_gap=float(o.LDE.GetSurfaceAt(6).Thickness), backfocus=float(o.LDE.GetSurfaceAt(15).Thickness),
                fixed_stop_radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter), nominal_fno_target=2.8,
                frequency_center=150. if cfg['wd']==200 else 225., frequency_edge=120. if cfg['wd']==200 else 180.))
        write('_validation_configs', poses)
        write('_run_summary', info)
    finally:
        if opt is not None:
            if opt.IsRunning and opt.CanCancel:
                opt.Cancel()
            opt.Close()
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('--cycles', type=int, choices=[1, 5, 10], default=5)
    a.add_argument('--seconds', type=float, default=600.)
    a.add_argument('--prepare-only', action='store_true')
    a.add_argument('--guarded', action='store_true')
    a.add_argument('--tag', default=TAG)
    q = a.parse_args()
    if not q.tag.startswith('candidate_v32_native_'):
        raise ValueError('Output tag must stay in the isolated v32 native candidate family.')
    TAG = q.tag
    run(q.cycles, q.seconds, q.prepare_only, q.guarded)
