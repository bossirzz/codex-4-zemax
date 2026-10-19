"""Audit a declared real-stop-aimed candidate and its focus positions in Zemax.

Acceptance values come from the native engine. Candidate labels never imply
that untested thermal, tolerance, or stray-light requirements have passed.
"""
import os
import argparse,json,math,sys
from pathlib import Path
import numpy as np
import zospy as zp
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,focus_file=None,nominal_fno=None,all_config_curves=False,all_field_threshold=None,fine_focus_step=None):
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
 try:
  p=z.connect('standalone',message_logging=False);o=p.copy_system()
  o.load(ROOT/(tag+'.zos'),saveifneeded=False)
  typ=zp.constants.Editors.MFE.MeritOperandType
  def val(kind,*args):
   pp=list(args)+[0.]*(8-len(args))
   return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
  def mtf(fc,fe,samp=6):
   return [val('MTFT',samp,0,1,fc,1,0),val('MTFS',samp,0,1,fc,1,0)]+[
    val(k,samp,0,i,fe,1,0) for i in range(2,5) for k in ['MTFT','MTFS']]
  cfgs=json.loads((ROOT/focus_file).read_text(encoding='utf8')) if focus_file else [
   dict(wd=500.,gap=float(o.LDE.GetSurfaceAt(11).Thickness),backfocus=float(o.LDE.GetSurfaceAt(15).Thickness),frequency_center=225.,frequency_edge=180.)]
  fields=[float(o.SystemData.Fields.GetField(i).Y) for i in range(1,5)]
  base_iris_gap=float(o.LDE.GetSurfaceAt(6).Thickness)
  rows=[]
  for cfg in cfgs:
   wd=cfg['wd'];fc=cfg['frequency_center'];fe=cfg['frequency_edge']
   o.LDE.GetSurfaceAt(0).Thickness=math.inf if wd>=1e7 else wd
   o.LDE.GetSurfaceAt(15).Thickness=cfg['backfocus']
   o.LDE.GetSurfaceAt(6).Thickness=cfg.get('iris_gap',base_iris_gap)
   o.LDE.GetSurfaceAt(11).Thickness=cfg['gap']
   if 'fixed_stop_radius' in cfg:o.LDE.GetSurfaceAt(6).SemiDiameter=cfg['fixed_stop_radius']
   target_fno=cfg.get('nominal_fno_target',nominal_fno)
   if target_fno is not None:
    # Set the requested test aperture explicitly for each focus configuration.
    iris=o.LDE.GetSurfaceAt(6)
    for _ in range(3):
     measured=val('ISFN')
     if not math.isfinite(measured) or measured<=0:raise ValueError('Invalid nominal F number')
     iris.SemiDiameter*=measured/float(target_fno)
   gap=o.LDE.GetSurfaceAt(11);scan=[]
   for delta in np.linspace(-.004,.004,9):
    gap.Thickness=cfg['gap']+float(delta)
    vv=[val('MTFT',4,0,1,fc,1,0),val('MTFT',4,0,4,fe,1,0),val('MTFS',4,0,4,fe,1,0)]
    if all_field_threshold is not None and wd in (500,200,800):
     vv.extend(val(k,4,0,field,fe,1,0) for field in [2,3] for k in ['MTFT','MTFS'])
    merit=min(vv[0]-.3,min(vv[1:])-.2) if wd in (500,200) else min(vv)
    if all_field_threshold is not None and wd in (500,200,800):merit=min(vv[0]-.3,min(vv[1:])-all_field_threshold)
    scan.append(dict(gap=float(gap.Thickness),mtf=vv,merit=merit))
   best=max(scan,key=lambda q:q['merit']);fine_scan=[]
   if fine_focus_step is not None and wd in (500,200,800):
    anchor=best['gap']
    for delta in np.arange(-3,4)*fine_focus_step:
     gap.Thickness=anchor+float(delta)
     vv=[val('MTFT',4,0,1,fc,1,0)]+[val(k,4,0,field,fe,1,0) for field in range(2,5) for k in ['MTFT','MTFS']]
     merit=min(vv[0]-.3,min(vv[1:])-(all_field_threshold if all_field_threshold is not None else .2))
     fine_scan.append(dict(gap=float(gap.Thickness),mtf=vv,merit=merit))
    best=max([best,*fine_scan],key=lambda q:q['merit'])
   gap.Thickness=best['gap']
   image=[val('REAY',16,2,0,y/fields[-1],0,0) for y in fields[1:]]
   rr=dict(wd=wd,gap=best['gap'],backfocus=cfg['backfocus'],frequencies=[fc,fe],
    mtf_sampling6=mtf(fc,fe),efl=val('EFFL',0,2),ttl=val('TOTR',0,0),
    cra_deg=abs(math.degrees(val('RANG',16,2,0,1,0,0))),relative_illumination=val('RELI',64,2,4,0),
    chief_image_heights=image,distortion_percent=[100*(h/y-1) for h,y in zip(image,fields[1:])],
    image_circle_diameter=2*image[-1],front_group_position=best['gap']+cfg['backfocus'],rear_group_position=cfg['backfocus'],
    iris_gap=float(o.LDE.GetSurfaceAt(6).Thickness),nominal_fno=val('ISFN'),working_fno=val('WFNO'),
    group_positions=[float(o.LDE.GetSurfaceAt(6).Thickness)+best['gap']+cfg['backfocus'],best['gap']+cfg['backfocus'],cfg['backfocus']],
    stop_radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter),ray_aiming=str(o.SystemData.RayAiming.RayAiming),native_focus_scan=scan,native_fine_focus_scan=fine_scan)
   rr['nominal_fno_target']=target_fno
   rr['focus_merit_offaxis_threshold']=all_field_threshold
   rr['paraxial_back_focal_length']=val('BFLL') if hasattr(typ,'BFLL') else None
   rr['relative_illumination_F_d_C']=[val('RELI',64,w,4,0) for w in [1,2,3]]
   colour=[[val('REAY',16,w,0,y/fields[-1],0,0) for w in [1,2,3]] for y in fields[1:]]
   rr['chief_colour_heights_F_d_C']=colour
   rr['chief_colour_span_um']=[1000*(max(x)-min(x)) for x in colour]
   if wd==500:
    glass_indices=[dict(surface=surface,material=str(o.LDE.GetSurfaceAt(surface).Material),
       indices=[val('INDX',surface,w) for w in range(1,6)]) for surface in [1,3,4,7,8,10,12,14]]
    (ROOT/(tag+'_native_glass_indices.json')).write_text(json.dumps(glass_indices,indent=2),encoding='utf8')
   rows.append(rr);print(json.dumps(rr,ensure_ascii=False),flush=True)
   (ROOT/(tag+'_native_acceptance.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf8')
  ordered=sorted(rows,key=lambda q:(q['wd']!=500,q['wd']))
  m=o.MCE
  for _ in ordered[1:]:m.AddConfiguration(False)
  variable_iris=len({round(x['stop_radius'],10) for x in ordered})>1
  for number,surface,key in [(1,0,'wd'),(2,11,'gap'),(3,15,'backfocus'),(4,6,'iris_gap')]:
   while m.NumberOfOperands<number:m.AddOperand()
   row=m.GetOperandAt(number);row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.THIC);row.Param1=surface
   for j,cfg in enumerate(ordered,1):row.GetOperandCell(j).DoubleValue=math.inf if key=='wd' and cfg[key]>=1e7 else cfg[key]
  if variable_iris:
   while m.NumberOfOperands<5:m.AddOperand()
   row=m.GetOperandAt(5)
   row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.SDIA);row.Param1=6
   for j,cfg in enumerate(ordered,1):row.GetOperandCell(j).DoubleValue=cfg['stop_radius']
  m.SetCurrentConfiguration(1)
  target=ROOT/(tag+'_review.zos');o.save_as(target)
  for samp in [512,1024]:
   a=zp.analyses.mtf.FFTMTF(sampling=f'{samp}x{samp}',maximum_frequency=300,field='All',wavelength='All')
   result=a.run(o,oncomplete='Close')
   result.data.to_csv(ROOT/(tag+f'_review_mtf_{samp}.csv'))
   (ROOT/(tag+f'_review_mtf_{samp}.json')).write_text(result.to_json(),encoding='utf8')
   print('FFT_SAMPLE',samp,'MESSAGES',[str(x) for x in result.messages or []],flush=True)
  if all_config_curves:
   for index,cfg in enumerate(ordered,1):
    if cfg['wd']!=200:continue
    m.SetCurrentConfiguration(index)
    for samp in [512,1024]:
     a=zp.analyses.mtf.FFTMTF(sampling=f'{samp}x{samp}',maximum_frequency=300,field='All',wavelength='All')
     result=a.run(o,oncomplete='Close')
     result.data.to_csv(ROOT/(tag+f'_near_mtf_{samp}.csv'))
     (ROOT/(tag+f'_near_mtf_{samp}.json')).write_text(result.to_json(),encoding='utf8')
     print('NEAR_FFT_SAMPLE',samp,'MESSAGES',[str(x) for x in result.messages or []],flush=True)
   m.SetCurrentConfiguration(1)
  aux={};ww=o.SystemData.Wavelengths
  for wave,label in [(4,'650nm'),(5,'850nm')]:
   weights=[ww.GetWavelength(i).Weight for i in range(1,6)]
   for i in range(1,6):ww.GetWavelength(i).Weight=1. if i==wave else 0.
   aux[label]=dict(same_visible_focus=[val('MTFT',6,wave,1,225,1,0),val('MTFT',6,wave,4,180,1,0),val('MTFS',6,wave,4,180,1,0)])
   for i,weight in enumerate(weights,1):ww.GetWavelength(i).Weight=weight
  (ROOT/(tag+'_auxiliary.json')).write_text(json.dumps(aux,indent=2),encoding='utf8')
  o.save();print('REVIEW_SAVED',str(target),'AUX',json.dumps(aux),flush=True)
 finally:
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--focus-file');a.add_argument('--nominal-fno',type=float);a.add_argument('--all-config-curves',action='store_true');a.add_argument('--all-field-threshold',type=float);a.add_argument('--fine-focus-step',type=float);q=a.parse_args();run(q.tag,focus_file=q.focus_file,nominal_fno=q.nominal_fno,all_config_curves=q.all_config_curves,all_field_threshold=q.all_field_threshold,fine_focus_step=q.fine_focus_step)
