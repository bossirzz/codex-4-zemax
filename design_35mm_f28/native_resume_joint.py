"""Warm-start a native joint DLS experiment from audited floating-group poses.

Python supplies the saved model, native merit weights and progress logging.
Every optical evaluation and every optimization step is performed by Zemax.
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
PARENT = 'candidate_v32_native_contrast_c'
POSE_SOURCE = 'candidate_v32_native_focus_refine'
sys.stdout.reconfigure(encoding='utf8')


def run(tag, cycles, seconds):
    if not tag.startswith('candidate_v33_native_'):
        raise ValueError('Use an isolated v33 native output tag')
    if (R / (tag + '.zos')).exists():
        raise FileExistsError('Preserve the existing experiment; choose a new tag')
    schema = json.loads((R / (PARENT + '_schema.json')).read_text())
    if schema != json.loads((R / (POSE_SOURCE + '_schema.json')).read_text()):
        raise RuntimeError('Floating-group seed has a different common prescription')
    poses = json.loads((R / (POSE_SOURCE + '_validation_configs.json')).read_text())
    parent_audit = json.loads((R / (PARENT + '_merit_audit.json')).read_text())
    explicit = {q['index']: q for q in parent_audit['explicit_operands']}

    def write(suffix, data):
        (R / (tag + suffix + '.json')).write_text(json.dumps(data, indent=2), encoding='utf8')

    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = opt = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(R / (PARENT + '.zos'), saveifneeded=False)
        typ = zp.constants.Editors.MFE.MeritOperandType

        def val(kind, *args):
            pp = list(args)+[0.]*(8-len(args))
            return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

        m = o.MCE
        if m.NumberOfConfigurations != len(poses):
            raise RuntimeError('Unexpected configuration count')
        controls = {(str(m.GetOperandAt(i).Type), int(m.GetOperandAt(i).Param1)): m.GetOperandAt(i)
                    for i in range(1, m.NumberOfOperands+1)}
        m.SetCurrentConfiguration(1)
        materials = []
        for surface, glass in zip([1, 3, 4, 7, 8, 10, 12, 14], schema['glasses']):
            measured = [val('INDX', surface, wave) for wave in [1, 2, 3]]
            expected = [glass['indices'][key] for key in ['F', 'd', 'C']]
            if max(abs(a-b) for a, b in zip(measured, expected)) > 1e-10:
                raise RuntimeError(f'Unresolved glass at surface {surface}')
            materials.append(dict(surface=surface, name=str(o.LDE.GetSurfaceAt(surface).Material),
                                  indices_F_d_C=measured))
        for configuration, q in enumerate(poses, 1):
            for key, value in {('THIC', 0): 1e10 if q['wd'] >= 1e7 else q['wd'],
                ('THIC', 6): q['iris_gap'], ('THIC', 11): q['gap'],
                ('THIC', 15): q['backfocus'], ('SDIA', 6): q['fixed_stop_radius']}.items():
                controls[key].GetOperandCell(configuration).DoubleValue = value
            m.SetCurrentConfiguration(configuration)
            for _ in range(3):
                cell = controls[('SDIA', 6)].GetOperandCell(configuration)
                cell.DoubleValue = float(o.LDE.GetSurfaceAt(6).SemiDiameter)*val('ISFN')/2.8

        # All floors outside the agreed 500/200 mm acceptance are engineering
        # objectives. Native penalty weights do not constitute hard acceptance.
        goals = {1: (.315, .212, 1e6), 2: (.315, .212, 1e6), 3: (.315, .212, 1e6),
                 5: (.310, .180, 10000.), 6: (.310, .205, 20000.), 7: (.300, .040, 2000.)}
        phase_scales = {1: .15, 2: .15, 3: .15, 5: 2., 6: 2., 7: .3}
        mf = o.MFE
        if mf.NumberOfOperands != 3108:
            raise RuntimeError('Parent native merit differs from the audited version')
        current_config = None
        mtf_guard_rows = set()
        changes = []
        counts = Counter()
        for i in range(1, mf.NumberOfOperands+1):
            row = mf.GetOperandAt(i)
            kind = str(row.Type)
            counts[kind] += 1
            if kind == 'CONF':
                current_config = int(row.GetCellAt(2).IntegerValue)
            elif kind in ['MECS', 'MECT']:
                factor = phase_scales[current_config]
                row.Weight = float(row.Weight)*factor
            elif kind in ['MTFT', 'MTFS']:
                row.GetCellAt(2).IntegerValue = 3
                cell = row.GetCellAt(4)
                field = int(cell.IntegerValue if str(cell.DataType) == 'Integer' else cell.DoubleValue)
                guard = mf.GetOperandAt(i+1)
                if str(guard.Type) != 'OPGT' or int(guard.GetCellAt(2).IntegerValue) != i:
                    raise RuntimeError('Unexpected native MTF guard layout')
                centre, offaxis, weight = goals[current_config]
                guard.Target = centre if field == 1 else offaxis
                guard.Weight = weight
                mtf_guard_rows.add(i+1)
                changes.append(dict(index=i, configuration=current_config, type=kind,
                    sampling=3, grid=1, field=field, floor=float(guard.Target), weight=weight))
            elif kind in ['OPGT', 'OPLT'] and i not in mtf_guard_rows:
                if float(row.Weight) < 200000.:
                    continue
                reference = explicit.get(int(row.GetCellAt(2).IntegerValue), {})
                raw_kind = reference.get('type')
                row.Weight = 1e6
                if raw_kind == 'STHI':
                    row.Weight = 1e7
                    surface = int(reference['parameters'][0])
                    if schema['media'][surface-1] >= 0:
                        row.Target = .452
                elif raw_kind == 'RANG' and kind == 'OPLT':
                    row.Target = math.radians(9.90)
                    row.Weight = 1e8
                elif raw_kind == 'REAY' and current_config == 2 and not math.isclose(float(row.Target), 9.):
                    height = reference['parameters'][3]*9.08
                    row.Target = height*(1-.0058 if kind == 'OPGT' else 1+.0058)
            elif kind == 'ZTHI':
                row.Target, row.Weight = 6., 1e6
        m.SetCurrentConfiguration(1)
        o.save_as(R / (tag + '_prepared.zos'))
        initial = float(mf.CalculateMeritFunction())
        if not math.isfinite(initial) or initial >= 1e8:
            raise RuntimeError('Invalid native merit; optimization was not started')
        write('_contract', dict(parent=PARENT, floating_group_seed=POSE_SOURCE,
            optimizer='OpticStudio native DampedLeastSquares', common_prescription_variables=40,
            floating_group_and_stop_variables=28, materials=materials, mtf_guards=changes,
            phase_scales=phase_scales, operand_counts=dict(counts), native_infinity_api_value=1e10,
            main_acceptance='500 mm:225/180 lp/mm at .30/.20; 200 mm:150/120 at .30/.20',
            other_distances='Engineering targets, not new user acceptance thresholds',
            temperature_in_objective=False, independent_acceptance_required=True))
        opt = o.Tools.OpenLocalOptimization()
        opt.Algorithm = zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
        opt.Cycles = getattr(zp.constants.Tools.Optimization.OptimizationCycles,
            {1: 'Fixed_1_Cycle', 5: 'Fixed_5_Cycles', 10: 'Fixed_10_Cycles'}[cycles])
        opt.NumberOfCores = min(8, int(opt.MaxCores))
        if int(opt.Variables) != 68:
            raise RuntimeError(f'Expected 68 native variables, got {opt.Variables}')
        info = dict(status='starting', started_utc=datetime.now(timezone.utc).isoformat(),
            source=PARENT, floating_group_seed=POSE_SOURCE, algorithm=str(opt.Algorithm),
            requested_cycles=cycles, variables=int(opt.Variables), targets=int(opt.Targets),
            cores=int(opt.NumberOfCores), initial_merit=initial, full_requirements_passed=False)
        write('_progress', info)
        print('NATIVE_JOINT_START', json.dumps(info), flush=True)
        opt.Run()
        clock = time.monotonic()
        last = clock-10
        cancelled = False
        while opt.IsRunning:
            now = time.monotonic()
            if now-last >= 10:
                info.update(status='running', elapsed_seconds=now-clock, current_merit=float(opt.CurrentMeritFunction))
                write('_progress', info)
                print('NATIVE_JOINT_PROGRESS', json.dumps(info), flush=True)
                last = now
            if now-clock > seconds and not cancelled and opt.CanCancel:
                opt.Cancel()
                cancelled = True
            time.sleep(.5)
        info.update(status='cancelled_at_time_limit' if cancelled else 'finished',
            succeeded=bool(opt.Succeeded), error=str(opt.ErrorMessage),
            elapsed_seconds=time.monotonic()-clock, current_merit=float(opt.CurrentMeritFunction))
        opt.Close()
        opt = None
        write('_progress', info)
        write('_run_summary', info)
        print('NATIVE_JOINT_DONE', json.dumps(info), flush=True)
        m.SetCurrentConfiguration(1)
        o.save_as(R / (tag + '.zos'))
        output = json.loads(json.dumps(schema))
        output['curvatures'] = [0. if math.isinf(float(o.LDE.GetSurfaceAt(i).Radius)) else
                               1/float(o.LDE.GetSurfaceAt(i).Radius) for i in range(1, 16)]
        output['thicknesses'] = [float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1, 16)]
        for a in output['aspheres']:
            a['coefficients'] = [float(o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{j}')).DoubleValue) for j in [2, 3, 4]]
        output['fixed_stop_radius'] = float(o.LDE.GetSurfaceAt(6).SemiDiameter)
        output['revision_reason'] = 'Native joint optimization from audited 400/2000 mm focus poses; independent acceptance required.'
        write('_schema', output)
        final_poses = []
        for configuration, pose in enumerate(poses, 1):
            m.SetCurrentConfiguration(configuration)
            final_poses.append(dict(pose, iris_gap=float(o.LDE.GetSurfaceAt(6).Thickness),
                gap=float(o.LDE.GetSurfaceAt(11).Thickness), backfocus=float(o.LDE.GetSurfaceAt(15).Thickness),
                fixed_stop_radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter)))
        write('_validation_configs', final_poses)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', default='candidate_v33_native_mainrange')
    parser.add_argument('--cycles', type=int, choices=[1, 5, 10], default=5)
    parser.add_argument('--seconds', type=float, default=600.)
    args = parser.parse_args()
    run(args.tag, args.cycles, args.seconds)
