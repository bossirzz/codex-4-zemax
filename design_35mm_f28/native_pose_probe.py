"""Verify explicit focus poses in an isolated copy; never overwrite source data."""
import os
import argparse
import json
import math
import sys
from pathlib import Path

import zospy as zp

ROOT = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')


def run(tag, config_file, output_prefix, sampling=6):
    configs = json.loads((ROOT / config_file).read_text(encoding='utf8'))
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = None
    out = []
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        typ = zp.constants.Editors.MFE.MeritOperandType

        def val(kind, *args):
            pp = list(args) + [0.] * (8 - len(args))
            return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

        for cfg in configs:
            o.load(ROOT / (tag + '.zos'), saveifneeded=False)
            o.LDE.GetSurfaceAt(0).Thickness = math.inf if cfg['wd'] >= 1e7 else cfg['wd']
            for surface, key in [(6, 'iris_gap'), (11, 'gap'), (15, 'backfocus')]:
                o.LDE.GetSurfaceAt(surface).Thickness = cfg[key]
            iris = o.LDE.GetSurfaceAt(6)
            iris.SemiDiameter = cfg.get('fixed_stop_radius', float(iris.SemiDiameter))
            for _ in range(3):
                iris.SemiDiameter *= val('ISFN') / 2.8
            fc, fe = cfg['frequency_center'], cfg['frequency_edge']
            values = []
            for field in range(1, 5):
                for kind in ['MTFT', 'MTFS']:
                    values.append(val(kind, sampling, 0, field, fc if field == 1 else fe, 1, 0))
            heights = [val('REAY', 16, 2, 0, y / 9.08, 0, 0) for y in [3., 6.3, 9.08]]
            tt = [float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1, 16)]
            row = dict(wd=cfg['wd'], label=cfg.get('label', ''),
                       frequencies=[fc, fe], sampling=sampling, native_mtf=values,
                       gap=tt[10], iris_gap=tt[5], backfocus=tt[-1],
                       stop_radius=float(iris.SemiDiameter), nominal_fno=val('ISFN'),
                       efl=val('EFFL', 0, 2), ttl=val('TOTR', 0, 0),
                       relative_illumination=val('RELI', 64, 2, 4, 0),
                       cra_deg=abs(math.degrees(val('RANG', 16, 2, 0, 1, 0, 0))),
                       distortion_percent=[100 * (h / y - 1) for h, y in zip(heights, [3., 6.3, 9.08])],
                       image_circle_diameter=2 * heights[-1],
                       group_positions_from_sensor_mm=[sum(tt), sum(tt[6:]), sum(tt[11:])],
                       acceptance_basis=cfg.get('acceptance_basis', 'engineering diagnostic only'))
            out.append(row)
            (ROOT / (output_prefix + '.json')).write_text(json.dumps(out, indent=2), encoding='utf8')
            print('NATIVE_POSE', json.dumps(row), flush=True)
    finally:
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('config_file')
    parser.add_argument('output_prefix')
    parser.add_argument('--sampling', type=int, default=6)
    args = parser.parse_args()
    run(args.tag, args.config_file, args.output_prefix, args.sampling)
