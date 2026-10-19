"""Refine the actual acceptance margins before pursuing larger MTF targets."""
import argparse,json,sys,time
import numpy as np
from scipy.optimize import least_squares,minimize
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import JointModel,R
sys.stdout.reconfigure(encoding='utf8')

class Feasible(Exception):pass

def run(seed_file=None,steps=120,tag='candidate_v22',solver='least_squares',replacements=(),mid_best=.125,mid_near=.155,colour_limit_um=None):
 model=JointModel()
 if replacements:model.substitute([(int(k),v) for k,v in (x.split('=',1) for x in replacements)])
 seed=model.x if seed_file is None else np.load(R/seed_file)
 calls=0;best=float('inf');bestx=seed.copy();success=None;started=time.time()
 def evaluate(x,raw=False):
  nonlocal calls,best,bestx,success
  c,t,n,caps,stop=model.unpack(x)
  best_values=f.quantities(c,t,fixed_stop_radius=stop,low_frequency=False,return_diagnostics=colour_limit_um is not None)
  mb,rb,hb,ib=best_values[:4]
  # Use the same declared nominal F/2.8 at each conjugate, rather than
  # assuming that a constant iris diameter implies constant nominal F number.
  pn=e.paraxial(c,n,200.)
  marginal,_,ok=e.ray(c,n,e.N[1],200.,0.,0.,pn['efl']/5.6,pn['ep'])
  ns=abs(marginal[9])
  near_values=f.quantities(c,n,wd=200.,fields=np.array([0.,3.,6.3,9.08]),
   frequencies=(150.,120.),fixed_stop_radius=ns,low_frequency=False,return_diagnostics=colour_limit_um is not None)
  mn,rn,hn,inn=near_values[:4]
  pb=e.paraxial(c,t)
  glass=np.where(e.MEDIA[:14]>=0,.05,0.)
  # Units are scaled by their practical margins, not by their raw magnitude.
  constraints=np.r_[
   (f.blank_clearance(c,t,caps)+glass-.003)*20,
   (f.blank_clearance(c,n,caps)+glass-.003)*20,
   (caps-rb-.005)*10,(caps-rn-.005)*10,
   (np.array([3,6.3,9.08])*.0029-abs(hb[1:]-np.array([3,6.3,9.08])))*20,
   (np.array([3,6.3,9.08])*.0058-abs(hn[1:]-np.array([3,6.3,9.08])))*20,
   (0.045-abs(pb['efl']-35.))*5,
   (59.8-sum(t))*.5,(59.8-sum(n))*.5,
   1-c[:15]**2*caps*caps-.01,-.1*(ib+inn),
   10*(mb[0]-.307),10*(mb[-1]-.208),10*(mn[0]-.307),10*(mn[-1]-.208),
   5*(mb[1:-1]-mid_best).ravel(),5*(mn[1:-1]-mid_near).ravel()]
  ti=t.copy();ti[-1]+=5.;ti[10]=e.focus_group_gap(c,ti,1e7)+t[10]-e.focus_group_gap(c,t,500.)
  pi=e.paraxial(c,ti,1e7);ch=e.aim_chief(c,ti,e.N[1],1e7,9.08/pi['magnification'],pi['ep'])
  ray,_,ok=e.ray(c,ti,e.N[1],1e7,9.08/pi['magnification'],0.,ch,pi['ep'])
  cra=abs(np.rad2deg(np.arctan(ray[3])))
  constraints=np.r_[constraints,(9.96-cra),20*(f.blank_clearance(c,ti,caps)+glass-.003)]
  if colour_limit_um is not None:
   spans=best_values[4]['phase_centroid_colour_span_um'][1:]+near_values[4]['phase_centroid_colour_span_um'][1:]
   constraints=np.r_[constraints,.5*(colour_limit_um-np.array(spans))]
  residual=np.r_[np.minimum(constraints,0.),.0001*(x-seed)]
  score=float(residual@residual);calls+=1
  if score<best:
   best=score;bestx=x.copy();np.save(R/(tag+'_refine_best.npy'),x)
  if calls%100==0:
   print(json.dumps(dict(calls=calls,seconds=time.time()-started,best=best,mtf500=mb.tolist(),mtf200=mn.tolist(),
     minimum_constraint=float(constraints.min()),efl=pb['efl'],cra_infinity=cra,ttl200=float(sum(n))),ensure_ascii=False),flush=True)
  if constraints.min()>=-1e-4:
   success=x.copy();raise Feasible()
  return constraints if raw else residual
 delta=np.r_[np.full(14,.10),np.full(15,.75),np.full(18,2.),[2.,1.,8.,.6],np.full(6,2.)]
 lo=seed-delta;hi=seed+delta
 lo[14:29]=np.maximum(lo[14:29],[1,.1,1,.45,.3,.3,.45,1,.1,1,.15,1,.1,.45,8])
 lo[47:51]=[-1,-2,-5,10.05];hi[47:51]=[5,2.5,16,11.45]
 try:
  if solver=='slsqp':
   opt=minimize(lambda z:float(.0001*np.sum((z-seed)**2)),seed,method='SLSQP',
    jac=lambda z:.0002*(z-seed),bounds=list(zip(lo,hi)),
    constraints=[dict(type='ineq',fun=lambda z:evaluate(z,raw=True),
     jac=lambda z:approx_derivative(lambda zz:evaluate(zz,raw=True),z,method='3-point',
      abs_step=1e-5*(1+abs(z)),bounds=(lo,hi)))],options=dict(maxiter=steps,ftol=1e-10,disp=True))
  else:
   opt=least_squares(evaluate,seed,jac=lambda z:approx_derivative(evaluate,z,method='3-point',
    abs_step=1e-5*(1+abs(z)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',
    max_nfev=steps,ftol=1e-9,xtol=1e-10,gtol=1e-8,verbose=1)
  final=bestx;print('SOLVER_END',opt.status,opt.message,opt.nfev,flush=True)
 except Feasible:
  final=success;print('FEASIBLE_SURROGATE_ONLY',flush=True)
 model.export(final,tag)
 print('EXPORTED',tag,'REQUIRES_NATIVE_VALIDATION',flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--seed-file');a.add_argument('--steps',type=int,default=120);a.add_argument('--tag',default='candidate_v22');a.add_argument('--solver',choices=['least_squares','slsqp'],default='least_squares');a.add_argument('--replace',action='append',default=[]);a.add_argument('--mid-best',type=float,default=.125);a.add_argument('--mid-near',type=float,default=.155);a.add_argument('--colour-limit-um',type=float);q=a.parse_args();run(q.seed_file,q.steps,q.tag,q.solver,q.replace,q.mid_best,q.mid_near,q.colour_limit_um)
