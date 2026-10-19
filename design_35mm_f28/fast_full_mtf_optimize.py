"""MTF surrogate with lens curvatures, spacing and six aspheres as variables."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
class FeasibleFound(Exception):pass
def run():
 schema=json.loads((ROOT/'candidate_v12_schema.json').read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in schema['glasses']])
 c=np.r_[schema['curvatures'],np.zeros(18)];t=np.array(schema['thicknesses'])
 for a in schema['aspheres']:
  j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
 scales=np.tile([1e-5,1e-7,1e-9],6)
 old=np.load(ROOT/'otf_best_parameters.npy');c[15:]=old[:18]*scales;t[-1]=old[18]
 caps=np.array([schema['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 x=np.r_[c[e.CURVED]/.05,t,c[15:]/scales]
 def unpack(xx):
  cc=c.copy();cc[e.CURVED]=xx[:14]*.05;cc[15:]=xx[29:]*scales;return cc,xx[14:29].copy()
 best=np.inf;calls=0;successful=None
 def fun(xx):
  nonlocal best,calls,successful
  cc,tt=unpack(xx);m,radii,chiefs,invalid=f.quantities(cc,tt)
  clear=f.blank_clearance(cc,tt,caps);p=e.paraxial(cc,tt)
  target=np.array([[.40,.40,.72,.72],[.28,.30,.62,.62],[.28,.30,.62,.62],[.32,.32,.62,.62]])
  weights=np.array([[2.,2.,.15,.15],[.4,.4,.15,.15],[.4,.4,.15,.15],[2.,2.,.15,.15]])
  lower=np.array([[.335,.335],[.18,.18],[.18,.18],[.245,.245]])
  heights=np.array([0.,3.,6.3,9.08]);dist=abs(chiefs-heights)
  parts=[((m-target)*weights).ravel(),(np.maximum(0,lower-m[:,:2])*15).ravel(),
   np.maximum(0,-clear)*300,np.maximum(0,radii-caps+.02)*80,
   np.maximum(0,dist-.0027*heights)*200,np.array([invalid*10.,(p['efl']-35)*100.,max(0,sum(tt)-54.5)*50.])]
  rr=np.concatenate(parts);score=rr@rr;calls+=1
  if score<best:best=score;np.save(ROOT/'full_mtf_best_parameters.npy',xx)
  if calls%300==0:print('STEP',calls,best,'MTF',m.tolist(),'PHYS',clear.min(),'EFL',p['efl'],'TTL',sum(tt),flush=True)
  if (m[0,:2].min()>=.32 and m[-1,:2].min()>=.225 and m[1:3,:2].min()>=.165
      and clear.min()>=-.0001 and max(radii-caps)<0 and invalid==0
      and abs(p['efl']-35)<.015 and sum(tt)<54.55
      and np.all(dist<=.003*heights+1e-5)):
   successful=xx.copy();raise FeasibleFound()
  return rr
 lo=np.r_[np.full(14,-3.),[1.0,.1,1.,.65,.3,.3,.65,1.,.1,1.,.15,1.,.1,.65,8.],np.full(18,-20.)]
 hi=np.r_[np.full(14,3.),[8.,3.,8.,6.,10.,10.,6.,8.,3.,8.,12.,8.,8.,8.,25.],np.full(18,20.)]
 x=np.clip(x,lo+1e-7,hi-1e-7)
 print('INITIAL',f.quantities(*unpack(x))[0].tolist(),flush=True)
 try:
  opt=least_squares(fun,x,jac=lambda xx:approx_derivative(fun,xx,method='3-point',abs_step=5e-6*(1+abs(xx)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',max_nfev=150,ftol=1e-8,xtol=1e-9,gtol=1e-7)
  final=opt.x;print('SOLVER',opt.status,opt.cost,opt.nfev,flush=True)
 except FeasibleFound:
  final=successful;print('FEASIBLE_SURROGATE_FOUND',flush=True)
 cc,tt=unpack(final);np.save(ROOT/'full_mtf_optimized_parameters.npy',final)
 schema['curvatures']=cc[:15].tolist();schema['thicknesses']=tt.tolist()
 schema['aspheres']=[dict(surface=s,coefficients=cc[15+3*j:18+3*j].tolist()) for j,s in enumerate(f.SURFS)]
 (ROOT/'candidate_v14_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
 print('FINAL_SURROGATE',f.quantities(cc,tt)[0].tolist(),flush=True)
 import verify_design as v
 v.N=e.N
 audit=v.diffraction(cc,tt,fields=np.array([0.,3.,6.3,9.08]),grid=192)
 (ROOT/'candidate_v14_independent_mtf.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
 print('EXACT',[(d['MTF_225_tangential'],d['MTF_180_tangential'],d['MTF_180_sagittal']) for d in audit['fields']],flush=True)
if __name__=='__main__':run()
