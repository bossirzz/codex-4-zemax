"""Verify declared compensation positions in the original native thermal models."""
import os
import argparse,json,math,sys
from pathlib import Path
import numpy as np
import zospy as zp
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,targets_file=None,temperatures=None,output_suffix='_thermal_compensated_native'):
 targets=json.loads((R/(targets_file or tag+'_thermal_focus_targets.json')).read_text(encoding='utf8'))
 if temperatures is not None:targets=[x for x in targets if x['temperature_c'] in temperatures]
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None;out=[]
 try:
  p=z.connect('standalone',message_logging=False);o=p.copy_system();typ=zp.constants.Editors.MFE.MeritOperandType
  def val(kind,*args):
   pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
  for wd in sorted({x['wd'] for x in targets},key=lambda wd:(wd!=500.,wd)):
   source=R/(tag+f'_thermal_base_WD{wd:g}.zos');o.load(source,saveifneeded=False);m=o.MCE
   def control(kind,surf,index,value):
    matches=[m.GetOperandAt(i) for i in range(1,m.NumberOfOperands+1)
       if str(m.GetOperandAt(i).Type)==kind and int(m.GetOperandAt(i).Param1)==surf]
    if not matches and kind=='SDIA':
     radius=float(o.LDE.GetSurfaceAt(surf).SemiDiameter);m.AddOperand();row=m.GetOperandAt(m.NumberOfOperands)
     row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.SDIA);row.Param1=surf
     for j in range(1,m.NumberOfConfigurations+1):row.GetOperandCell(j).DoubleValue=radius
     matches=[row]
    if len(matches)!=1:raise ValueError(f'Expected one {kind} controller for surface {surf}')
    cell=matches[0].GetOperandCell(index);cell.MakeSolveFixed();cell.DoubleValue=float(value)
   temp_map={}
   for i in range(1,m.NumberOfConfigurations+1):
    m.SetCurrentConfiguration(i);temp_map[round(float(o.SystemData.Environment.Temperature),6)]=i
   for cfg in [x for x in targets if x['wd']==wd]:
    index=temp_map[round(cfg['temperature_c'],6)];m.SetCurrentConfiguration(index)
    for surface,key in [(6,'iris_gap'),(11,'gap'),(15,'backfocus')]:control('THIC',surface,index,cfg[key])
    control('SDIA',6,index,cfg['stop_radius'])
    for _ in range(2):control('SDIA',6,index,float(o.LDE.GetSurfaceAt(6).SemiDiameter)*val('ISFN')/2.8)
    fc,fe=cfg['frequency_center'],cfg['frequency_edge'];scan=[]
    if wd<1e7:
     for delta in np.linspace(-.002,.002,9):
      control('THIC',11,index,cfg['gap']+float(delta))
      vv=[val(k,4,0,field,fc if field==1 else fe,1,0) for field in range(1,5) for k in ['MTFT','MTFS']]
      scan.append(dict(gap=float(o.LDE.GetSurfaceAt(11).Thickness),mtf=vv,margin=min(min(vv[:2])-.3,min(vv[2:])-.2)))
     best=max(scan,key=lambda q:q['margin']);control('THIC',11,index,best['gap'])
     mtf=[val(k,6,0,field,fc if field==1 else fe,1,0) for field in range(1,5) for k in ['MTFT','MTFS']]
    else:mtf=None
    heights=[val('REAY',16,2,0,y/9.08,0,0) for y in [3.,6.3,9.08]]
    tt=[float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1,16)]
    rr=dict(wd=wd,temperature_c=float(o.SystemData.Environment.Temperature),spacer_tce_ppm_per_k=cfg['spacer_tce_ppm_per_k'],
      frequencies=[fc,fe],native_mtf=mtf,nominal_fno=val('ISFN'),working_fno=val('WFNO'),efl=val('EFFL',0,2),
      ttl_mm=val('TOTR',0,0),backfocus_mm=tt[-1],iris_gap=tt[5],gap=tt[10],stop_radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter),
      cra_deg=abs(math.degrees(val('RANG',16,2,0,1,0,0))),relative_illumination=val('RELI',64,2,4,0),
      distortion_percent=[100*(h/y-1) for h,y in zip(heights,[3,6.3,9.08])],image_circle_mm=2*heights[-1],
      group_positions_from_sensor_mm=[sum(tt),sum(tt[6:]),sum(tt[11:])],native_focus_scan=scan,
      sampled_mtf_passed=None if mtf is None else min(mtf[:2])>=.3 and min(mtf[2:])>=.2)
    out.append(rr);print('NATIVE_THERMAL',json.dumps(rr),flush=True)
    (R/(tag+output_suffix+'.json')).write_text(json.dumps(dict(mounting_assumption='Native thin-shell spacers at 23 ppm/K',
      active_three_group_temperature_compensation=True,manufactured_assembly_validated=False,cases=out),ensure_ascii=False,indent=2),encoding='utf8')
   m.SetCurrentConfiguration(temp_map[20.])
   o.save_as(R/(tag+output_suffix+f'_WD{wd:g}.zos'))
 finally:
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--targets-file');a.add_argument('--temperatures',nargs='+',type=float);a.add_argument('--output-suffix',default='_thermal_compensated_native');q=a.parse_args();run(q.tag,q.targets_file,q.temperatures,q.output_suffix)
