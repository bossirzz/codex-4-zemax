"""Export native, thermally scaled geometry and indices for joint focus fitting."""
import os
import argparse,json,math,sys
from pathlib import Path
import zospy as zp
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,spacer_tce=23.):
 source=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 configs=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None;out=[]
 try:
  p=z.connect('standalone',message_logging=False);o=p.copy_system();typ=zp.constants.Editors.MFE.MeritOperandType
  def val(kind,*args):
   pp=list(args)+[0.]*(8-len(args));return float(o.MFE.GetOperandValue(getattr(typ,kind),int(pp[0]),int(pp[1]),*pp[2:]))
  for wd in [500.,200.,1e7]:
   cfg=next(q for q in configs if q['wd']==wd)
   o.load(R/(tag+'.zos'),saveifneeded=False);o.LDE.GetSurfaceAt(0).Thickness=math.inf if wd>=1e7 else wd
   for surf,key in [(6,'iris_gap'),(11,'gap'),(15,'backfocus')]:o.LDE.GetSurfaceAt(surf).Thickness=cfg[key]
   o.LDE.GetSurfaceAt(6).SemiDiameter=cfg['stop_radius']
   env=o.SystemData.Environment;env.AdjustIndexToEnvironment=True;env.Temperature=20.;env.Pressure=1.
   for surf in range(1,16):
    row=o.LDE.GetSurfaceAt(surf)
    if not row.Material:row.TCE=spacer_tce
   o.MCE.RunTool_MakeThermal(1,2,-20.,60.,True)
   for index in range(1,o.MCE.NumberOfConfigurations+1):
    o.MCE.SetCurrentConfiguration(index)
    row=dict(wd=wd,temperature_c=float(env.Temperature),spacer_tce_ppm_per_k=spacer_tce,
     curvatures=[0. if math.isinf(float(o.LDE.GetSurfaceAt(i).Radius)) else 1./float(o.LDE.GetSurfaceAt(i).Radius) for i in range(1,16)],
     thicknesses=[float(o.LDE.GetSurfaceAt(i).Thickness) for i in range(1,16)],
     caps=[float(o.LDE.GetSurfaceAt(i).SemiDiameter) for i in range(1,16)],media=source['media'],
     aspheres=[dict(surface=a['surface'],coefficients=[float(o.LDE.GetSurfaceAt(a['surface']).GetSurfaceCell(
      getattr(zp.constants.Editors.LDE.SurfaceColumn,f'Par{j}')).DoubleValue) for j in [2,3,4]]) for a in source['aspheres']],
     indices=[[val('INDX',surf,w) for surf in range(1,16)] for w in [1,2,3]],
     nominal_fno=val('ISFN'),paraxial_back_focal_length=val('BFLL') if hasattr(typ,'BFLL') else None)
    if wd>=1e7:
     row.update(cra_deg=abs(math.degrees(val('RANG',16,2,0,1,0,0))),relative_illumination=val('RELI',64,2,4,0))
    out.append(row);print('THERMAL_SNAPSHOT',wd,row['temperature_c'],row['nominal_fno'],flush=True)
    (R/(tag+'_thermal_snapshots.json')).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
   o.MCE.SetCurrentConfiguration(1)
   o.save_as(R/(tag+f'_thermal_base_WD{wd:g}.zos'))
 finally:
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--spacer-tce',type=float,default=23.);q=a.parse_args();run(q.tag,q.spacer_tce)
