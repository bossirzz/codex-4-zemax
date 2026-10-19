"""Append independently verified focus configurations to a separate review file."""
import os
import argparse
import json
import math
import sys
from pathlib import Path

import zospy as zp

ROOT = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')


def run(tag, probe_file, output_tag, curves=False):
    if output_tag == tag:
        raise ValueError('An independent output name is required.')
    verified = json.loads((ROOT / probe_file).read_text(encoding='utf8'))
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(ROOT / (tag + '_review.zos'), saveifneeded=False)
        m = o.MCE
        controllers = {(str(m.GetOperandAt(i).Type), int(m.GetOperandAt(i).Param1)): m.GetOperandAt(i)
                       for i in range(1, m.NumberOfOperands + 1)}
        for cfg in verified:
            m.AddConfiguration(False)
            index = m.NumberOfConfigurations
            values = {('THIC', 0): cfg['wd'], ('THIC', 6): cfg['iris_gap'],
                      ('THIC', 11): cfg['gap'], ('THIC', 15): cfg['backfocus'],
                      ('SDIA', 6): cfg['stop_radius']}
            for key, value in values.items():
                if key not in controllers:
                    raise ValueError(f'Missing declared controller {key}')
                cell = controllers[key].GetOperandCell(index)
                cell.MakeSolveFixed()
                cell.DoubleValue = math.inf if key == ('THIC', 0) and value >= 1e7 else float(value)
            m.SetCurrentConfiguration(index)
            assert abs(float(o.LDE.GetSurfaceAt(11).Thickness) - cfg['gap']) < 1e-9
            assert abs(float(o.LDE.GetSurfaceAt(6).SemiDiameter) - cfg['stop_radius']) < 1e-9
            if curves:
                a = zp.analyses.mtf.FFTMTF(sampling='1024x1024', maximum_frequency=300,
                                         field='All', wavelength='All')
                result = a.run(o, oncomplete='Close')
                curve_prefix = output_tag + f'_wd{cfg["wd"]:g}_mtf_1024'
                result.data.to_csv(ROOT / (curve_prefix + '.csv'))
                (ROOT / (curve_prefix + '.json')).write_text(result.to_json(), encoding='utf8')
                print('EXTENDED_NATIVE_CURVE', curve_prefix, 'MESSAGES',
                      [str(x) for x in result.messages or []], flush=True)
        m.SetCurrentConfiguration(1)
        target = ROOT / (output_tag + '.zos')
        o.save_as(target)
        provenance = dict(source=tag + '_review.zos', appended_native_results=probe_file,
                          total_configurations=m.NumberOfConfigurations,
                          material_and_surface_shapes_unchanged=True,
                          complete_requirement_acceptance=False)
        (ROOT / (output_tag + '_provenance.json')).write_text(json.dumps(provenance, indent=2), encoding='utf8')
        print('EXTENDED_REVIEW_SAVED', str(target), json.dumps(provenance), flush=True)
    finally:
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('probe_file')
    parser.add_argument('output_tag')
    parser.add_argument('--curves', action='store_true')
    args = parser.parse_args()
    run(args.tag, args.probe_file, args.output_tag, args.curves)
