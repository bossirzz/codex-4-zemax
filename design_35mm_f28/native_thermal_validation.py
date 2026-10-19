"""Native thermal screening with declared 23 ppm/K air spacers and refocusing.

Uses OpticStudio's Make Thermal tool for curvature, thickness and asphere
scaling. This assumed mounting model is not a manufactured assembly test.
"""
import os
import argparse,json,sys,math
from pathlib import Path
import numpy as np
import zospy as zp
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,spacer_tce=23.):
 base=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
 schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None;out=[]
 try:
  p=z.connect('standalone',message_logging=False);o=p.copy_system()
  typ=zp.constants.Editors.MFE.MeritOperandType
  def val(kind,*args):
   pp=list(args)+[0.]*(8-len(args))
   return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
  for wd in [500.,200.]:
   cfg=next(x for x in base if x['wd']==wd)
   o.load(R/(tag+'.zos'),saveifneeded=False)
   o.LDE.GetSurfaceAt(0).Thickness=wd
   for surface,key in [(6,'iris_gap'),(11,'gap'),(15,'backfocus')]:o.LDE.GetSurfaceAt(surface).Thickness=cfg[key]
   o.LDE.GetSurfaceAt(6).SemiDiameter=cfg['stop_radius']
   env=o.SystemData.Environment;env.AdjustIndexToEnvironment=True;env.Temperature=20.;env.Pressure=1.
   for i in range(1,16):
    row=o.LDE.GetSurfaceAt(i)
    if not row.Material:row.TCE=spacer_tce
   m=o.MCE
   if m.NumberOfConfigurations!=1:raise ValueError('Thermal source must be a single declared conjugate')
   m.RunTool_MakeThermal(1,2,-20.,60.,True)
   declared=[dict(row=i,type=str(m.GetOperandAt(i).Type),surface=int(m.GetOperandAt(i).Param1)) for i in range(1,m.NumberOfOperands+1)]
   print('THERMAL_TOOL',wd,'CONFIGS',m.NumberOfConfigurations,'ROWS',m.NumberOfOperands,flush=True)
   def control(kind,surface,index,value):
    matches=[m.GetOperandAt(i) for i in range(1,m.NumberOfOperands+1)
             if str(m.GetOperandAt(i).Type)==kind and int(m.GetOperandAt(i).Param1)==surface]
    if not matches and kind=='SDIA':
     # Make Thermal may omit a stop aperture which had a fixed solve.
     m.AddOperand();row=m.GetOperandAt(m.NumberOfOperands)
     row.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.SDIA);row.Param1=surface
     for config_index in range(1,m.NumberOfConfigurations+1):row.GetOperandCell(config_index).DoubleValue=cfg['stop_radius']
     matches=[row]
    if len(matches)!=1:raise ValueError(f'Missing unique thermal controller: {kind} {surface}; found {len(matches)}')
    cell=matches[0].GetOperandCell(index)
    cell.SetSolveData(cell.CreateSolveType(zp.constants.Editors.SolveType.Fixed))
    cell.DoubleValue=float(value)
   fc,fe=cfg['frequencies']
   def points(sampling=6):
    return [val(k,sampling,0,i,fc if i==1 else fe,1,0) for i in range(1,5) for k in ['MTFT','MTFS']]
   temps=[]
   for index in range(1,m.NumberOfConfigurations+1):
    m.SetCurrentConfiguration(index);temp=float(env.Temperature);temps.append(temp)
    geometry_before=[dict(surface=i,radius=str(o.LDE.GetSurfaceAt(i).Radius),
                         thickness=float(o.LDE.GetSurfaceAt(i).Thickness),
                         semi_diameter=float(o.LDE.GetSurfaceAt(i).SemiDiameter),
                         tce_ppm_per_k=(float(schema['glasses'][schema['media'][i-1]]['ED'][0])
                           if schema['media'][i-1]>=0 else float(o.LDE.GetSurfaceAt(i).TCE))) for i in range(1,16)]
    fixed=points();fixed_fno=val('ISFN');initial_gap=float(o.LDE.GetSurfaceAt(11).Thickness)
    # A selectable iris keeps the test aperture at the user's nominal F/2.8.
    radius=float(o.LDE.GetSurfaceAt(6).SemiDiameter)*fixed_fno/2.8
    control('SDIA',6,index,radius)
    scan=[]
    deltas=[0.] if abs(temp-20.)<1e-6 else np.linspace(-.06,.06,31)
    for delta in deltas:
     control('THIC',11,index,initial_gap+float(delta))
     vv=[val('MTFT',4,0,1,fc,1,0),val('MTFT',4,0,4,fe,1,0),val('MTFS',4,0,4,fe,1,0)]
     scan.append(dict(gap=float(o.LDE.GetSurfaceAt(11).Thickness),mtf=vv,margin=min(vv[0]-.30,min(vv[1:])-.20)))
    coarse=max(scan,key=lambda q:q['margin'])
    if abs(temp-20.)>=1e-6:
     for delta in np.linspace(-.002,.002,9):
      control('THIC',11,index,coarse['gap']+float(delta))
      vv=[val('MTFT',4,0,1,fc,1,0),val('MTFT',4,0,4,fe,1,0),val('MTFS',4,0,4,fe,1,0)]
      scan.append(dict(gap=float(o.LDE.GetSurfaceAt(11).Thickness),mtf=vv,margin=min(vv[0]-.30,min(vv[1:])-.20)))
    best=max(scan,key=lambda q:q['margin']);control('THIC',11,index,best['gap'])
    refocused=points()
    rr=dict(wd_mm=wd,temperature_c=temp,spacer_tce_ppm_per_k=spacer_tce,
     index_adjustment_enabled=bool(env.AdjustIndexToEnvironment),pressure_atm=float(env.Pressure),
     frequencies_lp_mm=[fc,fe],fixed_focus_mtf=fixed,refocused_mtf=refocused,
     group_gap_correction_mm=best['gap']-initial_gap,nominal_fno=val('ISFN'),working_fno=val('WFNO'),
     ttl_mm=val('TOTR',0,0),relative_illumination=val('RELI',64,2,4,0),
     sampled_target_passed=min(refocused[:2])>=.30 and min(refocused[-2:])>=.20,
     geometry_before_refocus=geometry_before,native_focus_scan=scan)
    out.append(rr);print(json.dumps(rr,ensure_ascii=False),flush=True)
    (R/(tag+'_thermal_native_results.json')).write_text(json.dumps(dict(
     mounting_assumption=f'All air spacers use {spacer_tce:g} ppm/K; native thin-shell thermal pickup model',
     source='Ansys OpticStudio Make Thermal, raw glass TD/ED and AdjustIndexToEnvironment enabled',
     full_assembly_validated=False,cases=out),ensure_ascii=False,indent=2),encoding='utf8')
   if not {-20.,20.,60.}.issubset(set(temps)):raise ValueError('Thermal tool did not generate all required temperatures')
   nominal_index=None
   for i in range(1,m.NumberOfConfigurations+1):
    m.SetCurrentConfiguration(i)
    if abs(env.Temperature-20.)<1e-6:nominal_index=i;break
   if nominal_index is None:raise ValueError('Missing nominal thermal configuration')
   m.SetCurrentConfiguration(nominal_index)
   o.save_as(R/(tag+f'_thermal_TCE{spacer_tce:g}_WD{wd:g}.zos'))
   (R/(tag+f'_thermal_WD{wd:g}_controllers.json')).write_text(json.dumps(declared,indent=2),encoding='utf8')
 finally:
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--spacer-tce',type=float,default=23.);q=a.parse_args();run(q.tag,q.spacer_tce)
