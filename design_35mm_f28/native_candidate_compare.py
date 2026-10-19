"""Compare the actual infinity chief angles before imposing finite-distance proxies."""
import os
import json,math,sys
from pathlib import Path
import numpy as np
import zospy as zp
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
def matrix(c,t,n):
 a,b,cc,d=1.,0.,0.,1.;ni=1.
 for i in range(15):
  power=(n[i]-ni)*c[i];cc,d=cc-power*a,d-power*b
  if i<14:a,b=a+t[i]/n[i]*cc,b+t[i]/n[i]*d
  ni=n[i]
 return a,b,cc,d
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('extension',instance_number=1,message_logging=False);o=p.copy_system()
 typ=zp.constants.Editors.MFE.MeritOperandType
 def val(kind,*args):
  pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 results=[]
 for tag,file in [('candidate_v08','candidate_v08_refined.zos'),('candidate_v10','candidate_v10_refined.zos'),('candidate_v11','candidate_v11_refined.zos')]:
  schema=json.loads((ROOT/(tag+'_schema.json')).read_text(encoding='utf8'))
  c=np.array(schema['curvatures']);t=np.array(schema['thicknesses']);n=np.array([1. if g<0 else schema['glasses'][g]['nd'] for g in schema['media']])
  o.load(ROOT/file,saveifneeded=False)
  for mode in ['WD500','infinity_whole_group','infinity_front_six']:
   tt=t.copy()
   if mode!='WD500':
    o.LDE.GetSurfaceAt(0).Thickness=math.inf
    if mode=='infinity_whole_group':
     a,b,cc,d=matrix(c,tt,n);tt[-1]=-a/cc
    else:
     tt[10]=0.;a,b,cc,d=matrix(c,tt,n);v0=a+t[-1]*cc
     tt[10]=1.;a,b,cc,d=matrix(c,tt,n);v1=a+t[-1]*cc;tt[10]=-v0/(v1-v0)
   o.LDE.GetSurfaceAt(11).Thickness=float(tt[10]);o.LDE.GetSurfaceAt(15).Thickness=float(tt[-1])
   results.append(dict(candidate=tag,mode=mode,efl=val('EFFL',0,2),gap=tt[10],backfocus=tt[-1],ttl=float(sum(tt)),
    chief_angle_deg=math.degrees(abs(val('RANG',16,2,0.,1.,0.,0.))),relative_illumination=val('RELI',32,2,4,0.),
    image_height=val('REAY',16,2,0.,1.,0.,0.)))
   print(json.dumps(results[-1]),flush=True)
 ROOT.joinpath('candidate_comparison.json').write_text(json.dumps(results,indent=2),encoding='utf8')
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
