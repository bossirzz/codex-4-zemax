"""Direct native MTF optimization with fixed real glasses and shared blank apertures."""
import os
from pathlib import Path
import json
import sys
import numpy as np
import zospy as zp

ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run():
 z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None;opt=None
 try:
  p=z.connect('extension',instance_number=1,message_logging=False);o=p.copy_system();o.load(ROOT/'candidate_v05.zos',saveifneeded=False)
  groups=[[1,2],[3,4,5],[7,8,9],[10,11],[12,13],[14,15]]
  caps={}
  for group in groups:
   cap=max(o.LDE.GetSurfaceAt(s).SemiDiameter for s in group)+.15
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
  operand('MTFT',[3,0,1,225,0,0],.46,3)
  for field in [2,3,4]:
   for kind in ['MTFT','MTFS']:operand(kind,[3,0,field,180,0,0],.46,2 if field==4 else 1)
  for field in [1,3,4]:
   for kind in ['MTFT','MTFS']:operand(kind,[3,0,field,90,0,0],.8,.2)
  for surf in [1,3,4,7,8,10,12,14]:operand('ETGT',[surf],.5,1000)
  operand('CTGT',[15],8.,1000)
  mf.CalculateMeritFunction()
  print('INITIAL_MERIT',mf.CalculateMeritFunction(),flush=True)
  opt=o.Tools.OpenLocalOptimization();opt.Algorithm=zp.constants.Tools.Optimization.OptimizationAlgorithm.DampedLeastSquares
  opt.Cycles=zp.constants.Tools.Optimization.OptimizationCycles.Fixed_50_Cycles;opt.NumberOfCores=min(8,int(opt.MaxCores))
  opt.RunAndWaitForCompletion();print('OPT_DONE',opt.Succeeded,opt.InitialMeritFunction,opt.CurrentMeritFunction,str(opt.ErrorMessage),flush=True)
  opt.Close();opt=None
  out=ROOT/'candidate_v06_native_mtf_optimized.zos';o.save_as(out)
  schema=json.loads((ROOT/'candidate_v05_schema.json').read_text(encoding='utf8'))
  for i in range(1,16):
   row=o.LDE.GetSurfaceAt(i);schema['curvatures'][i-1]=1/row.Radius if np.isfinite(row.Radius) else 0.;schema['thicknesses'][i-1]=row.Thickness
  for asp in schema['aspheres']:
   row=o.LDE.GetSurfaceAt(asp['surface']);asp['coefficients']=[row.GetSurfaceCell(getattr(zp.constants.Editors.LDE.SurfaceColumn,f'Par{par}')).DoubleValue for par in [2,3,4]]
  schema['fixed_semi_diameters']={str(k):v for k,v in caps.items()}
  ROOT.joinpath('candidate_v06_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
  ROOT.joinpath('native_optimization_merit.json').write_text(json.dumps(dict(variables=variables,operands=operands,final_merit=mf.CalculateMeritFunction()),indent=2),encoding='utf8')
  a=zp.analyses.mtf.FFTMTF(sampling='1024x1024',maximum_frequency=300)
  r=a.run(o,oncomplete='Close');ROOT.joinpath('candidate_v06_native_mtf.json').write_text(r.to_json(),encoding='utf8');r.data.to_csv(ROOT/'candidate_v06_native_mtf.csv')
  for freq in [225,180]:print(r.data.loc[r.data.index.to_series().sub(freq).abs().nsmallest(1).index].to_string(),flush=True)
  print('MESSAGES',str(r.messages),flush=True)
 finally:
  if opt is not None:opt.Close()
  if o is not None:o.close(saveifneeded=False)
  if z.Application is not None:z.disconnect()

if __name__=='__main__':run()
