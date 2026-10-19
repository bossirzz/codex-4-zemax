"""Run native four-variable DLS focus refinements at weak conjugates."""
import os
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import time

import zospy as zp

R = Path(__file__).resolve().parent
SOURCE = 'candidate_v32_native_contrast_c'
TAG = 'candidate_v32_native_focus_refine'
sys.stdout.reconfigure(encoding='utf8')
schema = json.loads((R / (SOURCE + '_schema.json')).read_text())
reference = json.loads((R / (SOURCE + '_native_post_audit.json')).read_text())
poses = json.loads((R / (SOURCE + '_validation_configs.json')).read_text())
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = opt = None
results = []
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / (SOURCE + '.zos'), saveifneeded=False)
    m = o.MCE
    controllers = {(str(m.GetOperandAt(i).Type), int(m.GetOperandAt(i).Param1)): m.GetOperandAt(i)
                   for i in range(1, m.NumberOfOperands+1)}
    for s in range(1, 16):
        o.LDE.GetSurfaceAt(s).RadiusCell.MakeSolveFixed()
    for s in [8, 10]:
        o.LDE.GetSurfaceAt(s).ThicknessCell.MakeSolveFixed()
    for a in schema['aspheres']:
        for par in [2, 3, 4]:
            o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{par}')).MakeSolveFixed()
    typ = zp.constants.Editors.MFE.MeritOperandType

    def val(kind, *args):
        pp = list(args)+[0.]*(8-len(args))
        return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

    def mtf():
        return [val(k, 6, 0, f, 225. if f == 1 else 180., 1, 0)
                for f in range(1, 5) for k in ['MTFT', 'MTFS']]

    def score(v):
        return min(min(v[:2])-.3, min(v[2:])-.2)

    variable_keys = [('THIC', 6), ('THIC', 11), ('THIC', 15), ('SDIA', 6)]
    for configuration, wd in [(5, 300), (6, 400), (7, 2000)]:
        for key in variable_keys:
            for c in range(1, m.NumberOfConfigurations+1):
                controllers[key].GetOperandCell(c).MakeSolveFixed()
            controllers[key].GetOperandCell(configuration).MakeSolveVariable()
        m.SetCurrentConfiguration(configuration)
        original = [float(o.LDE.GetSurfaceAt(s).Thickness) if k == 'THIC' else
                    float(o.LDE.GetSurfaceAt(s).SemiDiameter) for k, s in variable_keys]
        mf = o.MFE
        while mf.NumberOfOperands > 1:
            mf.RemoveOperandAt(mf.NumberOfOperands)
        mf.GetOperandAt(1).ChangeType(typ.BLNK)
        mf.GetOperandAt(1).Weight = 0.
        for centre_only, frequency in [(False, 180.), (True, 225.)]:
            w = mf.SEQOptimizationWizard2
            w.StartAt = 1 if not centre_only else mf.NumberOfOperands+1
            w.Type = zp.constants.Wizards.OptimizationTypes.RMS
            w.Criterion = zp.constants.Wizards.CriterionTypes.Contrast
            w.Reference = zp.constants.Wizards.ReferenceTypes.Centroid
            w.UseAllConfigurations = False
            w.ConfigurationNumber = configuration
            w.FieldNumber = 1 if centre_only else 0
            w.UseAllFields = not centre_only
            w.UseGaussianQuadrature = True
            w.Rings = 3
            w.Arms = zp.constants.Wizards.PupilArmsCount.Arms_6
            w.SpatialFrequency = frequency
            w.OverallWeight = .05 if centre_only else .1
            w.UseGlassBoundaryValues = False
            w.UseAirBoundaryValues = False
            w.UseMaximumDistortion = False
            w.IgnoreLateralColor = False
            w.DeleteVignetted = False
            w.Apply()

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
            return mf.NumberOfOperands

        def limits(raw, low=None, high=None, weight=20000.):
            if low is not None:
                operand('OPGT', [raw], low, weight)
            if high is not None:
                operand('OPLT', [raw], high, weight)

        operand('CONF', [configuration])
        for field in range(1, 5):
            for kind in ['MTFT', 'MTFS']:
                raw = operand(kind, [3, 0, field, 225. if field == 1 else 180., 1, 0])
                limits(raw, .31 if field == 1 else .205, weight=1000.)
        limits(operand('ISFN'), 2.7997, 2.8003)
        limits(operand('TTHI', [1, 15]), high=59.9)
        limits(operand('TTHI', [6, 6]), low=.3)
        limits(operand('TTHI', [11, 11]), low=.1)
        limits(operand('TTHI', [15, 15]), low=8.)
        limits(operand('RELI', [20, 2, 4, 0]), low=.75)
        for height in [3., 6.3, 9.08]:
            raw = operand('REAY', [16, 2, 0, height/9.08, 0, 0])
            distortion = .003 if wd == 400 else .006
            limits(raw, height*(1-distortion), height*(1+distortion))
        for start in [1, 7, 12]:
            operand('ZTHI', [start, 15], 6., 100000.)
        initial = float(mf.CalculateMeritFunction())
        if not math.isfinite(initial) or initial >= 1e8:
            raise RuntimeError('Invalid focus merit function')
        opt = o.Tools.OpenLocalOptimization()
        opt.Algorithm = zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
        opt.Cycles = zp.constants.Tools.Optimization.OptimizationCycles.Fixed_10_Cycles
        opt.NumberOfCores = min(8, int(opt.MaxCores))
        if int(opt.Variables) != 4:
            raise RuntimeError(f'Expected four native focus variables, got {opt.Variables}')
        info = dict(wd=wd, configuration=configuration, started_utc=datetime.now(timezone.utc).isoformat(),
            optimizer='OpticStudio DLS', requested_cycles=10, variables=int(opt.Variables),
            targets=int(opt.Targets), initial_merit=initial, source=SOURCE)
        print('NATIVE_FOCUS_START', json.dumps(info), flush=True)
        clock = time.monotonic()
        opt.Run()
        cancelled = False
        while opt.IsRunning:
            if time.monotonic()-clock > 180 and not cancelled and opt.CanCancel:
                opt.Cancel()
                cancelled = True
            time.sleep(.5)
        info.update(elapsed_seconds=time.monotonic()-clock, final_merit=float(opt.CurrentMeritFunction),
                    succeeded=bool(opt.Succeeded), error=str(opt.ErrorMessage), time_limit_reached=cancelled)
        opt.Close()
        opt = None
        m.SetCurrentConfiguration(configuration)
        after = mtf()
        before = next(q['mtf_sampling6'] for q in reference if q['wd'] == wd)
        selected = score(after) > score(before)
        optimized_pose = [float(o.LDE.GetSurfaceAt(s).Thickness) if k == 'THIC' else
                          float(o.LDE.GetSurfaceAt(s).SemiDiameter) for k, s in variable_keys]
        if not selected:
            for key, value in zip(variable_keys, original):
                controllers[key].GetOperandCell(configuration).DoubleValue = value
        chosen = optimized_pose if selected else original
        for key, value in zip(['iris_gap', 'gap', 'backfocus', 'fixed_stop_radius'], chosen):
            poses[configuration-1][key] = value
        info.update(native_mtf_sampling6_before=before, native_mtf_sampling6_after=after,
                    selected=selected, original_pose=original, optimized_pose=optimized_pose,
                    engineering_goals_only=True, complete_requirements_passed=False)
        results.append(info)
        (R / (TAG + '_runs.json')).write_text(json.dumps(results, indent=2), encoding='utf8')
        print('NATIVE_FOCUS_DONE', json.dumps(info), flush=True)
    m.SetCurrentConfiguration(1)
    o.save_as(R / (TAG + '.zos'))
    (R / (TAG + '_schema.json')).write_text(json.dumps(schema, indent=2), encoding='utf8')
    (R / (TAG + '_validation_configs.json')).write_text(json.dumps(poses, indent=2), encoding='utf8')
finally:
    if opt is not None:
        if opt.IsRunning and opt.CanCancel:
            opt.Cancel()
        opt.Close()
    if o is not None:
        o.close(saveifneeded=False)
    if z.Application is not None:
        z.disconnect()
