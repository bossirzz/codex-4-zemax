"""Direct native MTF optimization with fixed real glasses and shared blank apertures."""
import os
from pathlib import Path
import json
import sys
import time
import numpy as np
import zospy as zp

ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run():
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None;opt=None
 try:
  p=z.connect('extension',instance_number=1,message_logging=False);o=p.copy_system();o.load(ROOT/'candidate_v07_full_optimized.zos',saveifneeded=False)
  groups=[[1,2],[3,4,5],[7,8,9],[10,11],[12,13],[14,15]]
  caps={}
  for group in groups:
   cap=max(o.LDE.GetSurfaceAt(s).SemiDiameter for s in group)
   if cap>11.5:raise RuntimeError('Blank exceeds diameter limit')
   for s in group:
    row=o.LDE.GetSurfaceAt(s)
    sd=row.SemiDiameterCell.CreateSolveType(zp.constants.Editors.SolveType.Fixed);row.SemiDiameterCell.SetSolveData(sd)
    row.SemiDiameter=cap;caps[s]=cap
  variables=[]
  for s in [1,10,15]:
   for par in [2,3,4]:
    cell=o.LDE.GetSurfaceAt(s).GetSurfaceCell(getattr(zp.constants.Editors.LDE.SurfaceColumn,f'Par{par}'))
    cell.MakeSolveVariable();variables.append((s,f'Par{par}'))
  o.LDE.GetSurfaceAt(15).ThicknessCell.MakeSolveVariable();variables.append((15,'Thickness'))
  for s in range(1,16):
   row=o.LDE.GetSurfaceAt(s)
   if s!=6:row.RadiusCell.MakeSolveVariable();variables.append((s,'Radius'))
   row.ThicknessCell.MakeSolveVariable();variables.append((s,'Thickness'))
  mf=o.MFE
  while mf.NumberOfOperands>1:mf.RemoveOperandAt(mf.NumberOfOperands)
  operands=[]
  def operand(kind,params,target,weight):
   row=mf.GetOperandAt(1) if not operands else mf.AddOperand()
   row.ChangeType(getattr(zp.constants.Editors.MFE.MeritOperandType,kind))
   for i,v in enumerate(params,2):
    cell=row.GetCellAt(i)
    if str(cell.DataType)=='Integer':cell.IntegerValue=int(v)
    else:cell.DoubleValue=float(v)
   row.Target=target;row.Weight=weight
   operands.append(dict(type=kind,parameters=params,target=target,weight=weight))
   return len(operands)
  operand('MTFT',[4,0,1,225,0,0],.46,30)
  for field in [2,3,4]:
   for kind in ['MTFT','MTFS']:operand(kind,[4,0,field,180,0,0],.46,20)
  for field in [1,3,4]:
   for kind in ['MTFT','MTFS']:operand(kind,[4,0,field,90,0,0],.8,2.)
  for surf in [1,3,4,7,8,10,12,14]:
   operand('ETGT',[surf,0,1],.5,1000);operand('CTGT',[surf],.65,1000);operand('CTLT',[surf],8.,1000)
  for surf in [2,5,6,9,11,13]:operand('CTGT',[surf],.1,1000)
  for surf in [2,9,11,13]:operand('ETGT',[surf,0,1],.1,1000)
  operand('CTGT',[15],8.,1000)
  operand('EFFL',[0,2],35.,1000)
  tt=operand('TOTR',[1,16],0.,0.);operand('OPLT',[tt],59.5,1000)
  a=operand('TTHI',[5,5],0.,0.);b=operand('SAGY',[5],0.,0.);df=operand('DIFF',[a,b],0.,0.);operand('OPGT',[df],.3,1000)
  a=operand('TTHI',[6,6],0.,0.);b=operand('SAGY',[7],0.,0.);su=operand('SUMM',[a,b],0.,0.);operand('OPGT',[su],.3,1000)
  for h in [3.,6.3,9.08]:operand('REAY',[16,2,0.,h/9.08,0.,0.],h,30.)
  ca=operand('RANG',[16,2,0.,1.,0.,0.],0.,0.);operand('OPLT',[ca],np.deg2rad(9.2),1000)
  ri=operand('RELI',[20,2,4,0],0.,0.);operand('OPGT',[ri],.8,1000)
  mf.CalculateMeritFunction()
  print('INITIAL_OPERANDS',json.dumps([dict(index=i,type=str(mf.GetOperandAt(i).Type),value=mf.GetOperandAt(i).Value,target=mf.GetOperandAt(i).Target,contribution=mf.GetOperandAt(i).Contribution) for i in range(1,mf.NumberOfOperands+1)]),flush=True)
  print('INITIAL_MERIT',mf.CalculateMeritFunction(),flush=True)
  opt=o.Tools.OpenLocalOptimization();opt.Algorithm=zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
  opt.Cycles=zp.constants.Tools.Optimization.OptimizationCycles.Fixed_50_Cycles;opt.NumberOfCores=min(8,int(opt.MaxCores))
  opt.Run();stamp=time.monotonic()
  while opt.IsRunning:
   time.sleep(1)
   if time.monotonic()-stamp>10:
    print('OPT_PROGRESS',float(opt.Progress),float(opt.CurrentMeritFunction),flush=True);stamp=time.monotonic()
  print('OPT_DONE',opt.Succeeded,opt.InitialMeritFunction,opt.CurrentMeritFunction,str(opt.ErrorMessage),flush=True)
  opt.Close();opt=None
  out=ROOT/'candidate_v08_refined.zos';o.save_as(out)
  schema=json.loads((ROOT/'candidate_v07_schema.json').read_text(encoding='utf8'))
  for i in range(1,16):
   row=o.LDE.GetSurfaceAt(i);schema['curvatures'][i-1]=1/row.Radius if np.isfinite(row.Radius) else 0.;schema['thicknesses'][i-1]=row.Thickness
  for asp in schema['aspheres']:
   row=o.LDE.GetSurfaceAt(asp['surface']);asp['coefficients']=[row.GetSurfaceCell(getattr(zp.constants.Editors.LDE.SurfaceColumn,f'Par{par}')).DoubleValue for par in [2,3,4]]
  schema['fixed_semi_diameters']={str(k):v for k,v in caps.items()}
  ROOT.joinpath('candidate_v08_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
  ROOT.joinpath('native_refined_optimization_merit.json').write_text(json.dumps(dict(variables=variables,operands=operands,final_merit=mf.CalculateMeritFunction()),indent=2),encoding='utf8')
  a=zp.analyses.mtf.FFTMTF(sampling='1024x1024',maximum_frequency=300)
  r=a.run(o,oncomplete='Close');ROOT.joinpath('candidate_v08_native_mtf.json').write_text(r.to_json(),encoding='utf8');r.data.to_csv(ROOT/'candidate_v08_native_mtf.csv')
  for freq in [225,180]:print(r.data.loc[r.data.index.to_series().sub(freq).abs().nsmallest(1).index].to_string(),flush=True)
  print('MESSAGES',str(r.messages),flush=True)
 finally:
  if opt is not None:opt.Close()
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':run()
