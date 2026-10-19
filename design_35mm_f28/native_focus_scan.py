import os
import json
from pathlib import Path
import sys
import numpy as np
import zospy as zp

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False)
 o=p.copy_system();o.load(ROOT/'candidate_v01.zos',saveifneeded=False)
 row=o.LDE.GetSurfaceAt(o.LDE.NumberOfSurfaces-2);base=row.Thickness
 typ=zp.constants.Editors.MFE.MeritOperandType
 def mtf(field,freq,kind):
  return o.MFE.GetOperandValue(getattr(typ,kind),4,0,field,freq,0.,0.,0.,0.)
 records=[]
 for d in np.linspace(-.06,.06,61):
  row.Thickness=base+float(d)
  values=[mtf(1,225.,'MTFT'),mtf(4,180.,'MTFT'),mtf(4,180.,'MTFS')]
  records.append(dict(offset=float(d),values=values,margin=float(min(values[0]-.3,values[1]-.2,values[2]-.2))))
 records.sort(key=lambda r:r['margin'],reverse=True)
 ROOT.joinpath('native_focus_scan.json').write_text(json.dumps(records,indent=2),encoding='utf8')
 print(json.dumps(records[:8]),flush=True)
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
