"""Read installed native wizard/tool interfaces without changing the source lens."""
import os
import json
import sys
from pathlib import Path

import zospy as zp

R = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = opt = None
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / 'candidate_v31_thermal_probe_review.zos', saveifneeded=False)
    w = o.MFE.SEQOptimizationWizard
    names = [x for x in dir(w) if not x.startswith('_')]
    values = {}
    for name in names:
        try:
            value = getattr(w, name)
            if not callable(value):
                values[name] = str(value)
        except Exception as ex:
            values[name] = type(ex).__name__
    opt = o.Tools.OpenLocalOptimization()
    out = dict(wizard_names=names, wizard_defaults=values,
               wizard_enum_names=[x for x in dir(zp.constants.Wizards) if not x.startswith('_')],
               optimization_cycles=str(zp.constants.Tools.Optimization.OptimizationCycles),
               max_cores=int(opt.MaxCores), variables=int(opt.Variables), targets=int(opt.Targets))
    out['mfe_wizard_names'] = [x for x in dir(o.MFE) if 'Wizard' in x]
    out['legacy_data_names'] = [str(w.GetDataTypeAt(i)) for i in range(w.NumberOfDataTypes)]
    w2 = getattr(o.MFE, 'SEQOptimizationWizard2', None)
    if w2 is not None:
        out['wizard2_names'] = [x for x in dir(w2) if not x.startswith('_')]
        out['wizard2_defaults'] = {}
        for name in out['wizard2_names']:
            try:
                value = getattr(w2, name)
                if not callable(value):
                    out['wizard2_defaults'][name] = str(value)
            except Exception:
                pass
    (R / 'native_optimization_api_inspection.json').write_text(json.dumps(out, indent=2), encoding='utf8')
    print(json.dumps(out), flush=True)
finally:
    if opt is not None:
        opt.Close()
    if o is not None:
        o.close(saveifneeded=False)
    if z.Application is not None:
        z.disconnect()
