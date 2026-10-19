import os
import json
import math
from pathlib import Path
import sys
import numpy as np
import zospy as zp

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False)
 o=p.copy_system();o.load(ROOT/'candidate_v02.zos',saveifneeded=False)
 typ=zp.constants.Editors.MFE.MeritOperandType
 configs=json.loads((ROOT/'candidate_v02_focus_targets.json').read_text(encoding='utf8'))
 results=[]
 def value(kind,*params):
  pp=list(params)+[0.]*(8-len(params))
  return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 for cfg in configs:
  wd=cfg['wd'];o.LDE.GetSurfaceAt(0).Thickness=math.inf if wd>1e7 else wd
  gaprow=o.LDE.GetSurfaceAt(11)
  records=[]
  for d in np.linspace(-.032,.016,25):
   gaprow.Thickness=cfg['gap']+float(d)
   vals=[value('MTFT',4,0,1,cfg['frequency_center'],0.,0.),value('MTFT',4,0,4,cfg['frequency_edge'],0.,0.),value('MTFS',4,0,4,cfg['frequency_edge'],0.,0.)]
   records.append(dict(gap=float(gaprow.Thickness),values=vals,margin=float(min(vals[0]-.3,vals[1]-.2,vals[2]-.2))))
  best=max(records,key=lambda r:r['margin']);gaprow.Thickness=best['gap']
  best.update(wd=wd,efl=value('EFFL',0,2),ttl=sum(o.LDE.GetSurfaceAt(i).Thickness for i in range(1,16)),semi_diameters=[o.LDE.GetSurfaceAt(i).SemiDiameter for i in range(1,16)],frequencies=[cfg['frequency_center'],cfg['frequency_edge']])
  results.append(best)
  print(json.dumps(best),flush=True)
 ROOT.joinpath('candidate_v02_conjugate_audit.json').write_text(json.dumps(results,indent=2),encoding='utf8')
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
