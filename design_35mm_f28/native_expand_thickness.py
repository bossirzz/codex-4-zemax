"""Native DLS control experiment with additional axial thickness variables."""
import os
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time

import zospy as zp

R = Path(__file__).resolve().parent
PARENT = 'candidate_v33_native_mainrange'
ADDED = [1, 2, 3, 4, 5, 7, 9, 12, 13, 14]
sys.stdout.reconfigure(encoding='utf8')


def run(tag, cycles, seconds):
    if not tag.startswith('candidate_v34_native_') or (R / (tag + '.zos')).exists():
        raise ValueError('Choose an unused isolated v34 native tag')
    schema = json.loads((R / (PARENT + '_schema.json')).read_text())
    poses = json.loads((R / (PARENT + '_validation_configs.json')).read_text())

    def write(suffix, value):
        (R / (tag + suffix + '.json')).write_text(json.dumps(value, indent=2), encoding='utf8')

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
        m.SetCurrentConfiguration(1)
        materials = []
        for surface, glass in zip([1, 3, 4, 7, 8, 10, 12, 14], schema['glasses']):
            actual = [val('INDX', surface, w) for w in [1, 2, 3]]
            expected = [glass['indices'][key] for key in ['F', 'd', 'C']]
            if max(abs(a-b) for a, b in zip(actual, expected)) > 1e-10:
                raise RuntimeError(f'Glass resolution differs at surface {surface}')
            materials.append(dict(surface=surface, name=str(o.LDE.GetSurfaceAt(surface).Material), indices_F_d_C=actual))
        stop = next(m.GetOperandAt(i) for i in range(1, m.NumberOfOperands+1)
                    if str(m.GetOperandAt(i).Type) == 'SDIA' and int(m.GetOperandAt(i).Param1) == 6)
        for cfg in range(1, m.NumberOfConfigurations+1):
            m.SetCurrentConfiguration(cfg)
            for _ in range(3):
                stop.GetOperandCell(cfg).DoubleValue = float(o.LDE.GetSurfaceAt(6).SemiDiameter)*val('ISFN')/2.8
        for surface in ADDED:
            o.LDE.GetSurfaceAt(surface).ThicknessCell.MakeSolveVariable()
        mf = o.MFE
        changed = []
        for i in range(1, mf.NumberOfOperands+1):
            row = mf.GetOperandAt(i)
            if str(row.Type) == 'ISFN':
                for index in [i+1, i+2]:
                    guard = mf.GetOperandAt(index)
                    if str(guard.Type) not in ['OPGT', 'OPLT'] or int(guard.GetCellAt(2).IntegerValue) != i:
                        raise RuntimeError('Unexpected aperture guard layout')
                    guard.Weight = 1e8
                    changed.append(index)
        m.SetCurrentConfiguration(1)
        o.save_as(R / (tag + '_prepared.zos'))
        initial = float(mf.CalculateMeritFunction())
        if not math.isfinite(initial) or initial >= 1e8:
            raise RuntimeError('Invalid native merit; no optimization started')
        write('_contract', dict(parent=PARENT, optimizer='OpticStudio native DampedLeastSquares',
            added_LDE_thickness_variable_surfaces=ADDED, materials=materials,
            merit='Inherit audited main-range merit; aperture guard weight raised to 1e8',
            changed_aperture_guard_rows=changed, independent_acceptance_required=True,
            other_distances='Engineering goals, not new user acceptance thresholds',
            temperature_in_objective=False, seconds_limit=seconds))
        opt = o.Tools.OpenLocalOptimization()
        opt.Algorithm = zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
        opt.Cycles = getattr(zp.constants.Tools.Optimization.OptimizationCycles,
            {1: 'Fixed_1_Cycle', 5: 'Fixed_5_Cycles', 10: 'Fixed_10_Cycles'}[cycles])
        opt.NumberOfCores = min(8, int(opt.MaxCores))
        if int(opt.Variables) != 78:
            raise RuntimeError(f'Expected 78 native variables; got {opt.Variables}')
        info = dict(status='starting', source=PARENT, started_utc=datetime.now(timezone.utc).isoformat(),
            requested_cycles=cycles, time_limit_seconds=seconds, variables=int(opt.Variables),
            targets=int(opt.Targets), cores=int(opt.NumberOfCores), algorithm=str(opt.Algorithm),
            initial_merit=initial, full_requirements_passed=False)
        print('NATIVE_THICKNESS_START', json.dumps(info), flush=True)
        write('_progress', info)
        clock = time.monotonic()
        last = clock-10
        cancelled = False
        opt.Run()
        while opt.IsRunning:
            now = time.monotonic()
            if now-last >= 10:
                info.update(status='running', elapsed_seconds=now-clock, current_merit=float(opt.CurrentMeritFunction))
                write('_progress', info)
                print('NATIVE_THICKNESS_PROGRESS', json.dumps(info), flush=True)
                last = now
            if now-clock > seconds and not cancelled and opt.CanCancel:
                opt.Cancel()
                cancelled = True
            time.sleep(.5)
        info.update(status='cancelled_at_time_limit' if cancelled else 'finished',
            elapsed_seconds=time.monotonic()-clock, current_merit=float(opt.CurrentMeritFunction),
            succeeded=bool(opt.Succeeded), error=str(opt.ErrorMessage))
        opt.Close()
        opt = None
        write('_progress', info)
        write('_run_summary', info)
        print('NATIVE_THICKNESS_DONE', json.dumps(info), flush=True)
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
        output['revision_reason'] = 'Native DLS with 10 additional glass/air thickness variables; requires independent acceptance.'
        write('_schema', output)
        final_poses = []
        for cfg, pose in enumerate(poses, 1):
            m.SetCurrentConfiguration(cfg)
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
    parser.add_argument('--tag', default='candidate_v34_native_thickness')
    parser.add_argument('--cycles', type=int, choices=[1, 5, 10], default=5)
    parser.add_argument('--seconds', type=float, default=600.)
    args = parser.parse_args()
    run(args.tag, args.cycles, args.seconds)
