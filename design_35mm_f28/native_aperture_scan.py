"""Physical front clear-aperture scan with actual MTF and illumination."""
import os
import sys,json
from pathlib import Path
import numpy as np
import zospy as zp
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False);o=p.copy_system();o.load(ROOT/'candidate_v12_refined.zos',saveifneeded=False)
 typ=zp.constants.Editors.MFE.MeritOperandType
 def val(kind,*args):
  pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 records=[]
 for radius in np.linspace(6.6,9.0,13):
  for s in [1,2]:o.LDE.GetSurfaceAt(s).SemiDiameter=float(radius)
  vals=[val('MTFT',4,0,1,225,0,0),val('MTFT',4,0,4,180,0,0),val('MTFS',4,0,4,180,0,0)]
  record=dict(radius=float(radius),values=vals,ri=val('RELI',32,2,4,0),margin=min(vals[0]-.3,vals[1]-.2,vals[2]-.2))
  records.append(record);print(json.dumps(record),flush=True)
 best=max(records,key=lambda r:r['margin'] if r['ri']>=.75 else -1)
 for s in [1,2]:o.LDE.GetSurfaceAt(s).SemiDiameter=best['radius']
 best['high_sampling']=[val('MTFT',6,0,1,225,1,0),val('MTFT',6,0,4,180,1,0),val('MTFS',6,0,4,180,1,0)]
 o.save_as(ROOT/'candidate_v15_aperture.zos')
 schema=json.loads((ROOT/'candidate_v12_schema.json').read_text(encoding='utf8'))
 for s in [1,2]:schema['fixed_semi_diameters'][str(s)]=best['radius']
 ROOT.joinpath('candidate_v15_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
 ROOT.joinpath('physical_aperture_scan.json').write_text(json.dumps(dict(best=best,records=records),indent=2),encoding='utf8')
 print('BEST',json.dumps(best),flush=True)
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
