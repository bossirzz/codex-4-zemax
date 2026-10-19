"""Native five-variable focus experiment adding independent L6 floating motion."""
import os
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
import zospy as zp

R = Path(__file__).resolve().parent
SOURCE = 'candidate_v34_native_thickness'
TAG = 'candidate_v35_native_fourgroup'
sys.stdout.reconfigure(encoding='utf8')
schema = json.loads((R / (SOURCE + '_schema.json')).read_text())
poses = json.loads((R / (SOURCE + '_validation_configs.json')).read_text())
reference = json.loads((R / (SOURCE + '_native_post_audit.json')).read_text())
parent_geometry = json.loads((R / (SOURCE + '_dense_geometry_audit.json')).read_text())
if (R / (TAG + '.zos')).exists():
    raise FileExistsError('Keep the existing four-group experiment')
if len(reference) != 7:
    raise RuntimeError('Parent must complete its independent native audit first')
if not all(next(q for q in reference if q['wd'] == wd)['primary_mtf_passed'] for wd in [500, 200]):
    raise RuntimeError('Parent must retain the primary MTF acceptance')
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = opt = None
results = []
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / (SOURCE + '.zos'), saveifneeded=False)
    m = o.MCE
    typ = zp.constants.Editors.MFE.MeritOperandType

    def val(kind, *args):
        pp = list(args)+[0.]*(8-len(args))
        return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

    m.SetCurrentConfiguration(1)
    for surface, glass in zip([1, 3, 4, 7, 8, 10, 12, 14], schema['glasses']):
        if max(abs(val('INDX', surface, w)-glass['indices'][key])
               for w, key in enumerate(['F', 'd', 'C'], 1)) > 1e-10:
            raise RuntimeError(f'Glass resolution differs at surface {surface}')
    common_gap = float(o.LDE.GetSurfaceAt(9).Thickness)
    for surface in range(1, 16):
        o.LDE.GetSurfaceAt(surface).RadiusCell.MakeSolveFixed()
    for surface in [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14]:
        o.LDE.GetSurfaceAt(surface).ThicknessCell.MakeSolveFixed()
    for a in schema['aspheres']:
        for par in [2, 3, 4]:
            o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{par}')).MakeSolveFixed()
    m.AddOperand()
    row = m.GetOperandAt(m.NumberOfOperands)
    row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.THIC)
    row.Param1 = 9
    for cfg in range(1, m.NumberOfConfigurations+1):
        row.GetOperandCell(cfg).DoubleValue = common_gap
        poses[cfg-1]['extra_thicknesses'] = {'9': common_gap}
    controls = {(str(m.GetOperandAt(i).Type), int(m.GetOperandAt(i).Param1)): m.GetOperandAt(i)
                for i in range(1, m.NumberOfOperands+1)}
    keys = [('THIC', 6), ('THIC', 11), ('THIC', 15), ('SDIA', 6), ('THIC', 9)]
    schema['focus_group_starts'] = [1, 7, 10, 12]
    schema['focus_architecture'] = 'Four floating groups: L1-L3+STOP / L4-L5 / L6 / L7-L8; mechanism not yet validated.'
    schema['revision_reason'] = 'Native focus experiment with independent L6 floating motion; optical common prescription inherited from v34.'

    def metric(values):
        return min(min(values[:2])-.3, min(values[2:])-.2)

    for cfg, wd in [(5, 300), (6, 400), (7, 2000)]:
        for key in keys:
            for c in range(1, m.NumberOfConfigurations+1):
                controls[key].GetOperandCell(c).MakeSolveFixed()
            controls[key].GetOperandCell(cfg).MakeSolveVariable()
        m.SetCurrentConfiguration(cfg)
        original = [float(o.LDE.GetSurfaceAt(s).Thickness) if kind == 'THIC' else
                    float(o.LDE.GetSurfaceAt(s).SemiDiameter) for kind, s in keys]
        mf = o.MFE
        while mf.NumberOfOperands > 1:
            mf.RemoveOperandAt(mf.NumberOfOperands)
        mf.GetOperandAt(1).ChangeType(typ.BLNK)
        mf.GetOperandAt(1).Weight = 0.
        for centre, frequency in [(False, 180.), (True, 225.)]:
            w = mf.SEQOptimizationWizard2
            w.StartAt = mf.NumberOfOperands+1 if centre else 1
            w.Type = zp.constants.Wizards.OptimizationTypes.RMS
            w.Criterion = zp.constants.Wizards.CriterionTypes.Contrast
            w.Reference = zp.constants.Wizards.ReferenceTypes.Centroid
            w.UseAllConfigurations = False
            w.ConfigurationNumber = cfg
            w.FieldNumber = 1 if centre else 0
            w.UseAllFields = not centre
            w.UseGaussianQuadrature = True
            w.Rings = 3
            w.Arms = zp.constants.Wizards.PupilArmsCount.Arms_6
            w.SpatialFrequency = frequency
            w.OverallWeight = .05 if centre else .1
            w.UseGlassBoundaryValues = False
            w.UseAirBoundaryValues = False
            w.UseMaximumDistortion = False
            w.IgnoreLateralColor = False
            w.DeleteVignetted = False
            w.Apply()

        def operand(kind, params=(), target=0., weight=0.):
            rr = mf.AddOperand()
            rr.ChangeType(getattr(typ, kind))
            for index, value in enumerate(params, 2):
                cell = rr.GetCellAt(index)
                if str(cell.DataType) == 'Integer':
                    cell.IntegerValue = int(value)
                else:
                    cell.DoubleValue = float(value)
            rr.Target, rr.Weight = float(target), float(weight)
            return mf.NumberOfOperands

        def limits(raw, low=None, high=None, weight=1e6):
            if low is not None:
                operand('OPGT', [raw], low, weight)
            if high is not None:
                operand('OPLT', [raw], high, weight)

        operand('CONF', [cfg])
        for field in range(1, 5):
            for kind in ['MTFT', 'MTFS']:
                limits(operand(kind, [3, 0, field, 225. if field == 1 else 180., 1, 0]),
                       .31 if field == 1 else .205, weight=1000.)
        operand('ISFN', target=2.8, weight=1e8)
        limits(operand('TTHI', [1, 15]), high=59.9)
        limits(operand('TTHI', [15, 15]), low=8.)
        caps = schema['fixed_semi_diameters']
        for surface, radius, required in [(6, caps['7'], .302),
                (9, min(caps['9'], caps['10']), .102),
                (11, min(caps['11'], caps['12']), .102)]:
            for fraction in [0., .25, .5, .75, 1.]:
                limits(operand('STHI', [surface, 0, radius*fraction, 0]), low=required)
        limits(operand('RELI', [20, 2, 4, 0]), low=.75)
        for height in [3., 6.3, 9.08]:
            delta = .003 if wd == 400 else .006
            limits(operand('REAY', [16, 2, 0, height/9.08, 0, 0]), height*(1-delta), height*(1+delta))
        for start in [1, 7, 10, 12]:
            operand('ZTHI', [start, 15], 6., 1e6)
        initial = float(mf.CalculateMeritFunction())
        if not math.isfinite(initial) or initial >= 1e8:
            raise RuntimeError('Invalid native four-group merit')
        opt = o.Tools.OpenLocalOptimization()
        opt.Algorithm = zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
        opt.Cycles = zp.constants.Tools.Optimization.OptimizationCycles.Fixed_10_Cycles
        opt.NumberOfCores = min(8, int(opt.MaxCores))
        if int(opt.Variables) != 5:
            raise RuntimeError(f'Expected five focus variables; got {opt.Variables}')
        info = dict(wd=wd, configuration=cfg, source=SOURCE, optimizer='OpticStudio DLS',
            requested_cycles=10, variables=int(opt.Variables), targets=int(opt.Targets), initial_merit=initial)
        print('NATIVE_FOURGROUP_START', json.dumps(info), flush=True)
        opt.Run()
        clock = time.monotonic()
        cancelled = False
        while opt.IsRunning:
            if time.monotonic()-clock > 180 and not cancelled and opt.CanCancel:
                opt.Cancel()
                cancelled = True
            time.sleep(.5)
        info.update(succeeded=bool(opt.Succeeded), error=str(opt.ErrorMessage),
            final_merit=float(opt.CurrentMeritFunction), elapsed_seconds=time.monotonic()-clock,
            time_limit_reached=cancelled)
        opt.Close()
        opt = None
        m.SetCurrentConfiguration(cfg)
        after = [val(kind, 6, 0, field, 225. if field == 1 else 180., 1, 0)
                 for field in range(1, 5) for kind in ['MTFT', 'MTFS']]
        before = next(q['mtf_sampling6'] for q in reference if q['wd'] == wd)
        optimized = [float(o.LDE.GetSurfaceAt(s).Thickness) if kind == 'THIC' else
                     float(o.LDE.GetSurfaceAt(s).SemiDiameter) for kind, s in keys]
        baseline_geometry = next(q for q in parent_geometry['cases'] if q['wd'] == wd)
        air_checks = []
        for surface, position in [(6, 0), (11, 1), (9, 4)]:
            gap = next(q for q in baseline_geometry['gaps'] if q['surface'] == surface)
            minimum = gap['minimum_mm']+optimized[position]-original[position]
            air_checks.append(dict(surface=surface, minimum_mm=minimum, required_mm=gap['required_mm']))
        actual_t = [float(o.LDE.GetSurfaceAt(s).Thickness) for s in range(1, 16)]
        matrix = np.eye(2)
        n_before = 1.
        for s, curvature in enumerate(schema['curvatures']):
            medium = schema['media'][s]
            n_after = schema['glasses'][medium]['indices']['d'] if medium >= 0 else 1.
            matrix = np.array([[1., 0.], [-(n_after-n_before)*curvature, 1.]]) @ matrix
            if s < 14:
                matrix = np.array([[1., actual_t[s]/n_after], [0., 1.]]) @ matrix
            n_before = n_after
        bfl = float(-matrix[0, 0]/matrix[1, 0])
        illumination = val('RELI', 64, 2, 4, 0)
        image = [val('REAY', 16, 2, 0, h/9.08, 0, 0) for h in [3., 6.3, 9.08]]
        distortion = max(abs(h/ref-1) for h, ref in zip(image, [3., 6.3, 9.08]))
        positions = []
        for c in range(1, m.NumberOfConfigurations+1):
            m.SetCurrentConfiguration(c)
            positions.append([val('TTHI', start, 15) for start in [1, 7, 10, 12]])
        strokes = np.ptp(np.array(positions), axis=0).tolist()
        m.SetCurrentConfiguration(cfg)
        feasible = all(q['minimum_mm'] >= q['required_mm']-1e-9 for q in air_checks)
        feasible = feasible and sum(actual_t) <= 60 and actual_t[-1] >= 8 and bfl >= 8
        feasible = feasible and max(strokes) <= 6.000001 and illumination >= .75 and image[-1] >= 9.
        feasible = feasible and distortion <= (.003 if wd == 400 else .006)
        selected = feasible and metric(after) > metric(before)+.001
        if not selected:
            for key, value in zip(keys, original):
                controls[key].GetOperandCell(cfg).DoubleValue = value
        chosen = optimized if selected else original
        for key, value in zip(['iris_gap', 'gap', 'backfocus', 'fixed_stop_radius'], chosen[:4]):
            poses[cfg-1][key] = value
        poses[cfg-1]['extra_thicknesses']['9'] = chosen[4]
        info.update(native_mtf_sampling6_before=before, native_mtf_sampling6_after=after,
            selected=selected, original_pose=original, optimized_pose=optimized,
            feasible_after=feasible, air_checks=air_checks, independent_paraxial_bfl_mm=bfl,
            group_strokes_mm=strokes, relative_illumination=illumination,
            max_distortion_percent=100*distortion,
            engineering_goals_only=True, complete_requirements_passed=False)
        results.append(info)
        (R / (TAG + '_runs.json')).write_text(json.dumps(results, indent=2), encoding='utf8')
        print('NATIVE_FOURGROUP_DONE', json.dumps(info), flush=True)
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
