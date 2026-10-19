import os
import json
from pathlib import Path
import sys
import numpy as np
import zospy as zp
sys.stdout.reconfigure(encoding='utf8')
ROOT=Path(__file__).resolve().parent
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False)
 o=p.copy_system();o.load(ROOT/'candidate_v08_refined.zos',saveifneeded=False)
 row=o.LDE.GetSurfaceAt(15);base=row.Thickness
 typ=zp.constants.Editors.MFE.MeritOperandType
 def mtf(field,freq,kind,samp=4):
  return o.MFE.GetOperandValue(getattr(typ,kind),samp,0,field,freq,0.,0.,0.,0.)
 records=[]
 for delta in np.linspace(-.020,.020,41):
  row.Thickness=base+float(delta)
  values=[mtf(1,225.,'MTFT'),mtf(4,180.,'MTFT'),mtf(4,180.,'MTFS')]
  record=dict(offset=float(delta),values=values,margin=float(min(values[0]-.3,values[1]-.2,values[2]-.2)))
  records.append(record)
 records.sort(key=lambda r:r['margin'],reverse=True)
 for r in records[:8]:
  row.Thickness=base+r['offset'];values=[mtf(1,225.,'MTFT',6),mtf(4,180.,'MTFT',6),mtf(4,180.,'MTFS',6)]
  r['high_sampling_values']=values;r['high_sampling_margin']=min(values[0]-.3,values[1]-.2,values[2]-.2)
 records.sort(key=lambda r:r.get('high_sampling_margin',-1),reverse=True)
 ROOT.joinpath('v08_focus_scan.json').write_text(json.dumps(records,indent=2),encoding='utf8')
 print('BEST',json.dumps(records[:8]),flush=True)
 best=records[0];row.Thickness=base+best['offset'];o.save_as(ROOT/'candidate_v09_focused.zos')
 schema=json.loads((ROOT/'candidate_v08_schema.json').read_text(encoding='utf8'));schema['thicknesses'][-1]=row.Thickness
 ROOT.joinpath('candidate_v09_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
