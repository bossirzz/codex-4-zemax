"""Native multi-conjugate audit, exact-frequency MTF, and saved configurations."""
import os
from datetime import datetime
import json
import math
from pathlib import Path
import sys
import numpy as np
import zospy as zp

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent
TAG='candidate_v04'

def run():
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
 try:
  p=z.connect('extension',instance_number=1,message_logging=False)
  o=p.copy_system();o.load(ROOT/(TAG+'.zos'),saveifneeded=False)
  typ=zp.constants.Editors.MFE.MeritOperandType
  configs=json.loads((ROOT/(TAG+'_focus_targets.json')).read_text(encoding='utf8'))
  results=[]
  def value(kind,*params):
   pp=list(params)+[0.]*(8-len(params))
   return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
  def mtfs(fc,fe,samp=4):
   return [value('MTFT',samp,0,1,fc,0.,0.),value('MTFT',samp,0,2,fe,0.,0.),value('MTFS',samp,0,2,fe,0.,0.),
    value('MTFT',samp,0,3,fe,0.,0.),value('MTFS',samp,0,3,fe,0.,0.),value('MTFT',samp,0,4,fe,0.,0.),value('MTFS',samp,0,4,fe,0.,0.)]
  for cfg in configs:
   wd=cfg['wd'];o.LDE.GetSurfaceAt(0).Thickness=math.inf if wd>1e7 else wd
   o.LDE.GetSurfaceAt(15).Thickness=cfg['backfocus']
   gaprow=o.LDE.GetSurfaceAt(11)
   records=[]
   for d in np.linspace(-.012,.012,49):
    gaprow.Thickness=cfg['gap']+float(d)
    vals=mtfs(cfg['frequency_center'],cfg['frequency_edge'])
    score=min(vals[0]-.3,min(vals[1:])-.2)
    records.append(dict(gap=float(gaprow.Thickness),values=vals,margin=float(score)))
   best=max(records,key=lambda r:r['margin']);gaprow.Thickness=best['gap']
   best.update(wd=wd,backfocus=cfg['backfocus'],efl=value('EFFL',0,2),
    ttl=sum(o.LDE.GetSurfaceAt(i).Thickness for i in range(1,16)),
    semi_diameters=[o.LDE.GetSurfaceAt(i).SemiDiameter for i in range(1,16)],
    frequencies=[cfg['frequency_center'],cfg['frequency_edge']],
    chief_angle_degrees=math.degrees(abs(value('RANG',16,2,0.,1.,0.,0.))),
    relative_illumination=value('RELI',64,2,4,0.),
    chief_image_heights=[value('REAY',16,2,0.,y/9.03,0.,0.) for y in [3.,6.3,9.03]])
   best['distortion_percent']=[100*(h/y-1) for h,y in zip(best['chief_image_heights'],[3.,6.3,9.03])]
   best['values_sampling6']=mtfs(cfg['frequency_center'],cfg['frequency_edge'],6)
   results.append(best)
   print(json.dumps(best),flush=True)
  ROOT.joinpath('native_full_audit.json').write_text(json.dumps(results,indent=2),encoding='utf8')
  # Best conjugate first. Only group air gaps and object distance vary.
  ordered=sorted(results,key=lambda r:(r['wd']!=500,r['wd']))
  m=o.MCE
  for _ in ordered[1:]:m.AddConfiguration(False)
  for number,surf,key in [(1,0,'wd'),(2,11,'gap'),(3,15,'backfocus')]:
   while m.NumberOfOperands<number:m.AddOperand()
   row=m.GetOperandAt(number);row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.THIC);row.Param1=surf
   for j,cfg in enumerate(ordered,1):row.GetOperandCell(j).DoubleValue=math.inf if key=='wd' and cfg[key]>1e7 else cfg[key]
  m.SetCurrentConfiguration(1)
  final=ROOT/'35mm_F28_Phi18_optimized.zos';o.save_as(final)
  for sample in [512,1024]:
   a=zp.analyses.mtf.FFTMTF(sampling=f'{sample}x{sample}',maximum_frequency=300)
   rr=a.run(o,oncomplete='Close')
   ROOT.joinpath(f'final_native_mtf_{sample}.json').write_text(rr.to_json(),encoding='utf8')
   rr.data.to_csv(ROOT/f'final_native_mtf_{sample}.csv')
  # Evaluate auxiliary wavelengths at the visible focal setting.
  aux={}
  for wave,label in [(4,'650nm'),(5,'850nm')]:
   aux[label]=[value('MTFT',6,wave,1,225.,0.,0.),value('MTFT',6,wave,4,180.,0.,0.),value('MTFS',6,wave,4,180.,0.,0.)]
  ROOT.joinpath('native_auxiliary_wavelengths.json').write_text(json.dumps(aux,indent=2),encoding='utf8')
  print('AUXILIARY',json.dumps(aux),flush=True)
  ROOT.joinpath('final_configurations.json').write_text(json.dumps(ordered,indent=2),encoding='utf8')
  o.save()
  print('FINAL_SAVED',str(final),flush=True)
 finally:
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':run()
