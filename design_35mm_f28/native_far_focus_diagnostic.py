"""Read-only native through-focus diagnosis on an isolated OpticStudio copy."""
import os
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import zospy as zp

ROOT = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')


def run(tag, wd=2000., extent=.12, points=41):
    cfg = next(q for q in json.loads((ROOT / (tag + '_native_acceptance.json')).read_text())
               if q['wd'] == wd)
    out = dict(candidate=tag, wd=wd, units='mm', sampling=6,
               frequencies_lp_per_mm=[225., 180.], cases=[],
               purpose='Diagnosis only; far-distance thresholds have not been agreed.')
    target = ROOT / (tag + f'_native_through_focus_WD{wd:g}.json')
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(ROOT / (tag + '.zos'), saveifneeded=False)
        typ = zp.constants.Editors.MFE.MeritOperandType

        def val(kind, *args):
            pp = list(args) + [0.] * (8 - len(args))
            return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

        o.LDE.GetSurfaceAt(0).Thickness = wd
        for surface, key in [(6, 'iris_gap'), (11, 'gap'), (15, 'backfocus')]:
            o.LDE.GetSurfaceAt(surface).Thickness = cfg[key]
        iris = o.LDE.GetSurfaceAt(6)
        iris.SemiDiameter = cfg['stop_radius']
        for _ in range(3):
            iris.SemiDiameter *= val('ISFN') / 2.8

        for delta in np.linspace(-extent, extent, points):
            o.LDE.GetSurfaceAt(11).Thickness = cfg['gap'] + float(delta)
            row = dict(delta_gap_mm=float(delta), gap=float(o.LDE.GetSurfaceAt(11).Thickness))
            for wave, label in [(0, 'F_d_C'), (2, 'd_line')]:
                row[label] = [val(kind, 6, wave, field, 225. if field == 1 else 180., 1, 0)
                              for field in range(1, 5) for kind in ['MTFT', 'MTFS']]
            out['cases'].append(row)
            target.write_text(json.dumps(out, indent=2), encoding='utf8')
            print('NATIVE_FOCUS', json.dumps(row), flush=True)
        peaks = {}
        for label in ['F_d_C', 'd_line']:
            peaks[label] = [dict(field_index=i // 2 + 1, direction=['T', 'S'][i % 2],
                                **{k: max(out['cases'], key=lambda q: q[label][i])[k]
                                   for k in ['gap', 'delta_gap_mm']},
                                mtf=max(q[label][i] for q in out['cases'])) for i in range(8)]
        out['per_field_peaks'] = peaks
        out['best_shared_F_d_C'] = max(out['cases'], key=lambda q: min(q['F_d_C']))
        out['complete'] = True
        target.write_text(json.dumps(out, indent=2), encoding='utf8')
        print('NATIVE_FOCUS_COMPLETE', json.dumps(dict(peaks=peaks,
              best_shared=out['best_shared_F_d_C'])), flush=True)
    finally:
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('--wd', type=float, default=2000.)
    parser.add_argument('--extent', type=float, default=.12)
    parser.add_argument('--points', type=int, default=41)
    args = parser.parse_args()
    run(args.tag, args.wd, args.extent, args.points)
