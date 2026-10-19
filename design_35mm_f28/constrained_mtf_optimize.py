"""Find a feasible real-glass prescription using explicit nonlinear constraints."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
class Found(Exception):pass
def run():
 schema=json.loads((ROOT/'candidate_v12_schema.json').read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in schema['glasses']])
 scales=np.tile([1e-5,1e-7,1e-9],6)
 caps=np.array([schema['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 def unpack(xx):
  c=np.zeros(33);c[e.CURVED]=xx[:14]*.05;c[15:]=xx[29:]*scales;return c,xx[14:29].copy()
 x=np.load(ROOT/'full_mtf_checkpoint_before_constraints.npy')
 lastx=None;last=None;calls=0;success=None;best=1e99
 heights=np.array([0.,3.,6.3,9.08])
 def evaluate(xx):
  nonlocal lastx,last,calls,best,success
  if lastx is not None and np.array_equal(xx,lastx):return last
  c,t=unpack(xx);m,radii,chiefs,invalid=f.quantities(c,t);p=e.paraxial(c,t)
  clear=f.blank_clearance(c,t,caps)
  clear+=np.where(e.MEDIA[:14]>=0,.05,0.)
  target=np.array([[.34,.34],[.24,.27],[.24,.27],[.24,.24]])
  weights=np.array([[1,1],[.1,.1],[.1,.1],[2,2]])
  objective=float(np.sum(((m[:,:2]-target)*weights)**2))
  constraints=np.r_[clear,radii*0+caps-radii,
       .003*heights+1e-7-abs(chiefs-heights),
       [(.015-abs(p['efl']-35))*.1,(54.5-sum(t))*.1,-invalid*.01],
       1-c[:15]**2*caps*caps-.01,
       (m[0,:2]-.305)*2,(m[-1,:2]-.210)*2]
  score=objective+np.dot(np.minimum(constraints,0),np.minimum(constraints,0))*100
  calls+=1
  if score<best:best=score;np.save(ROOT/'constrained_best_parameters.npy',xx)
  if calls%100==0:print('STEP',calls,'SCORE',best,'MTF',m[:,:2].tolist(),'MIN',constraints.min(),'EFL',p['efl'],'TTL',sum(t),flush=True)
  if (constraints.min()>-2e-5 and m[0,:2].min()>.305 and m[-1,:2].min()>.210):
   success=xx.copy();raise Found()
  lastx=xx.copy();last=(objective,constraints);return last
 lo=np.r_[np.full(14,-3.),[1.,.1,1.,.65,.3,.3,.65,1.,.1,1.,.15,1.,.1,.65,8.],np.full(18,-20.)]
 hi=np.r_[np.full(14,3.),[8.,3.,8.,6.,10.,10.,6.,8.,3.,8.,12.,8.,8.,8.,25.],np.full(18,20.)]
 try:
  opt=minimize(lambda xx:evaluate(xx)[0],np.clip(x,lo+1e-8,hi-1e-8),method='SLSQP',jac='2-point',bounds=list(zip(lo,hi)),
     constraints=[dict(type='ineq',fun=lambda xx:evaluate(xx)[1])],
     options=dict(maxiter=120,ftol=1e-10,finite_diff_rel_step=2e-6,disp=True))
  final=opt.x;print('SOLVER',opt.success,opt.message,flush=True)
 except Found:
  final=success;print('FEASIBLE_FOUND',flush=True)
 c,t=unpack(final)
 schema['curvatures']=c[:15].tolist();schema['thicknesses']=t.tolist()
 schema['aspheres']=[dict(surface=s,coefficients=c[15+3*j:18+3*j].tolist()) for j,s in enumerate(f.SURFS)]
 ROOT.joinpath('candidate_v17_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
 np.save(ROOT/'constrained_final_parameters.npy',final)
 print('FINAL',f.quantities(c,t)[0].tolist(),flush=True)
 import verify_design as v
 v.N=e.N
 audit=v.diffraction(c,t,fields=heights,grid=192)
 ROOT.joinpath('candidate_v17_independent_mtf.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
 print('EXACT',[(d['MTF_225_tangential'],d['MTF_180_tangential'],d['MTF_180_sagittal']) for d in audit['fields']],flush=True)
if __name__=='__main__':run()
