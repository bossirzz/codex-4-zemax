"""Independently audit a native optimizer's saved configurations, without refocus."""
import os
import argparse
import json
import math
from pathlib import Path
import sys

import zospy as zp

R = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')


def run(tag, curves=False, configurations=None, source_suffix='', output_suffix='', verify_schema=False,
        extra_curve_configurations=None, configuration_order=None):
    declared = json.loads((R / (tag + '_validation_configs.json')).read_text())
    model_schema = json.loads((R / (tag + '_schema.json')).read_text())
    out = []
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(R / (tag + source_suffix + '.zos'), saveifneeded=False)
        output_tag = tag + output_suffix
        typ = zp.constants.Editors.MFE.MeritOperandType

        def val(kind, *args):
            pp = list(args) + [0.] * (8 - len(args))
            return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

        if verify_schema:
            schema = model_schema
            o.MCE.SetCurrentConfiguration(1)
            errors = []
            for i in range(1, 16):
                surface = o.LDE.GetSurfaceAt(i)
                radius = float(surface.Radius)
                curvature = 0. if math.isinf(radius) else 1/radius
                errors.extend([abs(curvature-schema['curvatures'][i-1]),
                               abs(float(surface.Thickness)-schema['thicknesses'][i-1])])
            for a in schema['aspheres']:
                for par, expected in zip([2, 3, 4], a['coefficients']):
                    actual = float(o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
                        getattr(zp.constants.Editors.LDE.SurfaceColumn, f'Par{par}')).DoubleValue)
                    errors.append(abs(actual-expected))
            for key, expected in schema['fixed_semi_diameters'].items():
                errors.append(abs(float(o.LDE.GetSurfaceAt(int(key)).SemiDiameter)-expected))
            for s, glass in zip([1, 3, 4, 7, 8, 10, 12, 14], schema['glasses']):
                errors.extend(abs(val('INDX', s, w)-glass['indices'][key])
                              for w, key in enumerate(['F', 'd', 'C'], 1))
            waves = [(float(o.SystemData.Wavelengths.GetWavelength(i).Wavelength),
                      float(o.SystemData.Wavelengths.GetWavelength(i).Weight)) for i in range(1, 6)]
            for actual, expected in zip(waves, [(.4861327, 1.), (.5875618, 1.), (.6562725, 1.), (.65, 0.), (.85, 0.)]):
                errors.extend(abs(a-b) for a, b in zip(actual, expected))
            for i, height in enumerate(schema['fields'], 1):
                errors.extend([abs(float(o.SystemData.Fields.GetField(i).X)),
                               abs(float(o.SystemData.Fields.GetField(i).Y)-height)])
            for index, cfg in enumerate(declared, 1):
                o.MCE.SetCurrentConfiguration(index)
                for s, key in [(6, 'iris_gap'), (11, 'gap'), (15, 'backfocus')]:
                    errors.append(abs(float(o.LDE.GetSurfaceAt(s).Thickness)-cfg[key]))
                errors.append(abs(float(o.LDE.GetSurfaceAt(6).SemiDiameter)-cfg['fixed_stop_radius']))
                for surface, value in cfg.get('extra_thicknesses', {}).items():
                    errors.append(abs(float(o.LDE.GetSurfaceAt(int(surface)).Thickness)-value))
            check = dict(native_prescription_verified=max(errors) <= 1e-10, maximum_difference=max(errors),
                         wavelengths=waves, geometry_source=tag + '_schema.json',
                         poses_source=tag + '_validation_configs.json')
            (R / (output_tag + '_native_schema_check.json')).write_text(json.dumps(check, indent=2), encoding='utf8')
            if not check['native_prescription_verified']:
                raise RuntimeError(f'Native saved model differs from declared prescription: {max(errors)}')
            print('NATIVE_PRESCRIPTION_CHECK', json.dumps(check), flush=True)

        order = configuration_order or list(range(1, len(declared)+1))
        if sorted(order) != list(range(1, len(declared)+1)):
            raise ValueError('Audit order must include each configuration exactly once')
        for index in order:
            cfg = declared[index-1]
            if configurations and index not in configurations:
                continue
            o.MCE.SetCurrentConfiguration(index)
            if cfg['wd'] >= 1e7:
                controls = [o.MCE.GetOperandAt(i) for i in range(1, o.MCE.NumberOfOperands + 1)
                            if str(o.MCE.GetOperandAt(i).Type) == 'THIC' and int(o.MCE.GetOperandAt(i).Param1) == 0]
                controls[0].GetOperandCell(index).MakeSolveFixed()
                # IEEE infinity through this API collapses all off-axis rays.
                # 1e10 is recognized by OpticStudio as native Infinity.
                controls[0].GetOperandCell(index).DoubleValue = 1e10
                o.MCE.SetCurrentConfiguration(1)
                o.MCE.SetCurrentConfiguration(index)
            fc, fe = cfg['frequency_center'], cfg['frequency_edge']
            mtf = [val(kind, 6, 0, field, fc if field == 1 else fe, 1, 0)
                   for field in range(1, 5) for kind in ['MTFT', 'MTFS']]
            if not all(math.isfinite(v) and 0 <= v <= 1 for v in mtf):
                raise RuntimeError(f'Invalid native MTF in configuration {index}')
            image = [val('REAY', 16, 2, 0, h / 9.08, 0, 0) for h in [3., 6.3, 9.08]]
            if not all(math.isfinite(h) and h > 0 for h in image):
                raise RuntimeError(f'Invalid off-axis image heights in configuration {index}')
            t = [float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1, 16)]
            rr = dict(wd=cfg['wd'], configuration=index, frequencies=[fc, fe],
                      mtf_sampling6=mtf, nominal_fno=val('ISFN'), working_fno=val('WFNO'),
                      efl=val('EFFL', 0, 2), ttl=val('TOTR', 0, 0),
                      gap=t[10], iris_gap=t[5], backfocus=t[-1],
                      stop_radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter),
                      cra_deg=abs(math.degrees(val('RANG', 16, 2, 0, 1, 0, 0))),
                      relative_illumination=val('RELI', 64, 2, 4, 0),
                      chief_image_heights=image, image_circle_diameter=2 * image[-1],
                      distortion_percent=[100 * (h / ref - 1) for h, ref in zip(image, [3., 6.3, 9.08])],
                      group_positions_from_sensor_mm=[sum(t[start-1:]) for start in model_schema.get('focus_group_starts', [1, 7, 12])],
                      focus_group_starts=model_schema.get('focus_group_starts', [1, 7, 12]),
                      primary_mtf_passed=(min(mtf[:2]) >= .3 and min(mtf[2:]) >= .2) if cfg['wd'] in (500, 200) else None,
                      true_infinity=bool(cfg['wd'] >= 1e7), refocus_applied=False,
                      actual_object_distance=float(o.LDE.GetSurfaceAt(0).Thickness),
                      native_object_is_infinite=math.isinf(float(o.LDE.GetSurfaceAt(0).Thickness)),
                      audit_valid_for_acceptance=True,
                      complete_requirement_acceptance=False)
            out.append(rr)
            (R / (output_tag + '_native_post_audit.json')).write_text(json.dumps(out, indent=2), encoding='utf8')
            print('NATIVE_POST_AUDIT', json.dumps(rr), flush=True)
            if index == 1:
                analysis = zp.analyses.new_analysis(o, zp.constants.Analysis.AnalysisIDM.SystemData)
                try:
                    analysis.ApplyAndWaitForCompletion()
                    analysis.get_text_output(str(R / (output_tag + '_native_system.txt')), z.get_txtfile_encoding())
                finally:
                    analysis.Close()
            if (curves and cfg['wd'] in (500, 200)) or (extra_curve_configurations and index in extra_curve_configurations):
                suffix = ('_review' if cfg['wd'] == 500 else '_near') if cfg['wd'] in (500, 200) else f'_WD{cfg["wd"]:g}'
                samplings = [512, 1024] if cfg['wd'] in (500, 200) else [1024]
                for sampling in samplings:
                    analysis = zp.analyses.mtf.FFTMTF(sampling=f'{sampling}x{sampling}',
                        maximum_frequency=300, field='All', wavelength='All')
                    result = analysis.run(o, oncomplete='Close')
                    result.data.to_csv(R / (output_tag + suffix + f'_mtf_{sampling}.csv'))
                    (R / (output_tag + suffix + f'_mtf_{sampling}.json')).write_text(result.to_json(), encoding='utf8')
                    print('NATIVE_CURVES', cfg['wd'], sampling,
                          'MESSAGES', [str(v) for v in result.messages or []], flush=True)
    finally:
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('tag')
    a.add_argument('--curves', action='store_true')
    a.add_argument('--configs', type=int, nargs='+')
    a.add_argument('--source-suffix', default='')
    a.add_argument('--output-suffix', default='')
    a.add_argument('--verify-schema', action='store_true')
    a.add_argument('--extra-curve-configs', type=int, nargs='+')
    a.add_argument('--order', type=int, nargs='+')
    q = a.parse_args()
    run(q.tag, q.curves, q.configs, q.source_suffix, q.output_suffix, q.verify_schema, q.extra_curve_configs, q.order)
