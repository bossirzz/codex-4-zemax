"""Read native merit contributions and convergence on an isolated candidate."""
import os
import json
import math
from pathlib import Path
import sys

import zospy as zp

R = Path(__file__).resolve().parent
TAG = 'candidate_v32_native_contrast_b'
sys.stdout.reconfigure(encoding='utf8')
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = None
out = {}
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / (TAG + '.zos'), saveifneeded=False)
    typ = zp.constants.Editors.MFE.MeritOperandType

    def val(kind, *args):
        pp = list(args) + [0.] * (8-len(args))
        return float(o.MFE.GetOperandValue(getattr(typ, kind), int(pp[0]), int(pp[1]), *pp[2:]))

    out['merit'] = float(o.MFE.CalculateMeritFunction())
    rows = []
    cfg = None
    for i in range(1, o.MFE.NumberOfOperands+1):
        q = o.MFE.GetOperandAt(i)
        if str(q.Type) == 'CONF':
            cfg = int(q.GetCellAt(2).IntegerValue)
        if math.isfinite(float(q.Contribution)) and float(q.Contribution) > 0:
            rows.append(dict(index=i, configuration=cfg, type=str(q.Type), value=float(q.Value),
                             target=float(q.Target), weight=float(q.Weight), contribution=float(q.Contribution)))
    out['contributions'] = sorted(rows, key=lambda q:q['contribution'], reverse=True)[:50]
    print('MERIT_CONTRIBUTIONS', json.dumps(out), flush=True)
    out['sampling'] = []
    for cfg, freq in [(1, 180), (2, 120)]:
        o.MCE.SetCurrentConfiguration(cfg)
        for grid, sampling in [(1, 3), (1, 4), (0, 3), (0, 5), (0, 7)]:
            rr = dict(configuration=cfg, grid=grid, sampling=sampling,
                      mtf=[val(k, sampling, 0, 4, freq, grid, 0) for k in ['MTFT', 'MTFS']])
            out['sampling'].append(rr)
            print('SAMPLING', json.dumps(rr), flush=True)
    control = next(o.MCE.GetOperandAt(i) for i in range(1, o.MCE.NumberOfOperands+1)
                   if str(o.MCE.GetOperandAt(i).Type) == 'THIC' and int(o.MCE.GetOperandAt(i).Param1) == 0)
    out['infinity'] = []
    for distance in [1e7, 1e10, 1e15, math.inf]:
        control.GetOperandCell(4).DoubleValue = distance
        o.MCE.SetCurrentConfiguration(1)
        o.MCE.SetCurrentConfiguration(4)
        rr = dict(object_distance=float(o.LDE.GetSurfaceAt(0).Thickness),
                  image=val('REAY', 16, 2, 0, 1, 0, 0),
                  cra=math.degrees(val('RANG', 16, 2, 0, 1, 0, 0)),
                  fno=val('ISFN'))
        out['infinity'].append(rr)
        print('INFINITY', json.dumps(rr), flush=True)
    (R / (TAG + '_diagnosis.json')).write_text(json.dumps(out, indent=2), encoding='utf8')
finally:
    if o is not None:
        o.close(saveifneeded=False)
    if z.Application is not None:
        z.disconnect()
