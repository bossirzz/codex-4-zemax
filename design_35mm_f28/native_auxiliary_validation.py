import os
import argparse,json,sys
from pathlib import Path
import zospy as zp
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args()
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('standalone',message_logging=False);o=p.copy_system();o.load(R/(q.tag+'.zos'),saveifneeded=False)
 typ=zp.constants.Editors.MFE.MeritOperandType;ww=o.SystemData.Wavelengths
 def val(kind,*args):
  pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
 print('IMAGE_FNUM',val('ISFN'),'WORKING_FNUM',val('WFNO'),flush=True)
 profiles=json.loads((R/(q.tag+'_auxiliary_focus_profiles.json')).read_text(encoding='utf8'));out=[]
 for cfg in profiles:
  w=4 if cfg['wavelength_um']==.65 else 5
  for i in range(1,6):ww.GetWavelength(i).Weight=1. if i==w else 0.
  o.LDE.GetSurfaceAt(11).Thickness=cfg['gap'];o.LDE.GetSurfaceAt(15).Thickness=cfg['backfocus']
  rr=dict(wavelength_um=cfg['wavelength_um'],gap=cfg['gap'],backfocus=cfg['backfocus'],
   refocused_mtf=[val('MTFT',6,w,1,225,1,0),val('MTFT',6,w,4,180,1,0),val('MTFS',6,w,4,180,1,0)],
   relative_illumination=val('RELI',64,w,4,0),distortion_percent=100*(val('REAY',16,2,0,1,0,0)/9.08-1))
  out.append(rr);print(json.dumps(rr),flush=True)
  if cfg['wavelength_um']==.65 and rr['refocused_mtf'][0]>=.3 and min(rr['refocused_mtf'][1:])>=.2:
   target=R/(q.tag+'_650nm_refocused.zos');o.save_as(target);rr['saved_model']=str(target)
 (R/(q.tag+'_auxiliary_refocus_native.json')).write_text(json.dumps(out,indent=2),encoding='utf8')
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
