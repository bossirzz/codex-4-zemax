"""Compare chief rays and MTF with exact stop aiming before acceptance."""
import os
import sys,json,math
from pathlib import Path
import zospy as zp
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('standalone',message_logging=False);o=p.copy_system();o.load(ROOT/'candidate_v17.zos',saveifneeded=False)
 import System
 aiming=o.SystemData.RayAiming;typ=zp.constants.Editors.MFE.MeritOperandType
 def val(kind,*args):
  pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 records=[]
 for mode in ['Off','Paraxial','Real']:
  aiming.RayAiming=System.Enum.Parse(aiming.RayAiming.GetType(),mode)
  rec=dict(mode=mode,mtf=[val('MTFT',6,0,1,225,1,0),val('MTFT',6,0,4,180,1,0),val('MTFS',6,0,4,180,1,0)],
   cra=math.degrees(val('RANG',16,2,0,1,0,0)),ri=val('RELI',64,2,4,0),image_height=val('REAY',16,2,0,1,0,0))
  records.append(rec);print(json.dumps(rec),flush=True)
 ROOT.joinpath('v17_ray_aiming_comparison.json').write_text(json.dumps(records,indent=2),encoding='utf8')
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
