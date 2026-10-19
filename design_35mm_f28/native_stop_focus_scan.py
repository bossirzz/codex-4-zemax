"""Move the iris within a real existing air space and scan the common image plane."""
import os
import sys,json
from pathlib import Path
import numpy as np
import zospy as zp
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False);o=p.copy_system();o.load(ROOT/'candidate_v15_aperture.zos',saveifneeded=False)
 typ=zp.constants.Editors.MFE.MeritOperandType
 def val(kind,*args):
  pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 base=[o.LDE.GetSurfaceAt(s).Thickness for s in [5,6,15]]
 records=[]
 for shift in np.linspace(-.12,.72,15):
  o.LDE.GetSurfaceAt(5).Thickness=base[0]+float(shift);o.LDE.GetSurfaceAt(6).Thickness=base[1]-float(shift)
  local=[]
  for delta in np.linspace(-.012,.012,9):
   o.LDE.GetSurfaceAt(15).Thickness=base[2]+float(delta)
   values=[val('MTFT',4,0,1,225,0,0),val('MTFT',4,0,4,180,0,0),val('MTFS',4,0,4,180,0,0)]
   local.append(dict(stop_shift=float(shift),back_offset=float(delta),values=values,margin=min(values[0]-.3,values[1]-.2,values[2]-.2)))
  best=max(local,key=lambda r:r['margin']);records.extend(local);print('STEP_BEST',json.dumps(best),flush=True)
 best=max(records,key=lambda r:r['margin'])
 o.LDE.GetSurfaceAt(5).Thickness=base[0]+best['stop_shift'];o.LDE.GetSurfaceAt(6).Thickness=base[1]-best['stop_shift'];o.LDE.GetSurfaceAt(15).Thickness=base[2]+best['back_offset']
 best['high_sampling']=[val('MTFT',6,0,1,225,1,0),val('MTFT',6,0,4,180,1,0),val('MTFS',6,0,4,180,1,0)]
 best['ri']=val('RELI',64,2,4,0)
 o.save_as(ROOT/'candidate_v16_stop_focused.zos')
 schema=json.loads((ROOT/'candidate_v15_schema.json').read_text(encoding='utf8'))
 for s in [5,6,15]:schema['thicknesses'][s-1]=o.LDE.GetSurfaceAt(s).Thickness
 ROOT.joinpath('candidate_v16_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
 ROOT.joinpath('stop_focus_scan.json').write_text(json.dumps(dict(best=best,records=records),indent=2),encoding='utf8')
 print('BEST',json.dumps(best),flush=True)
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
