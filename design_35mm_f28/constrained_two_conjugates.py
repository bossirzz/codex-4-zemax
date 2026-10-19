"""Joint best/near design with real catalogue dispersion and a common iris."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,least_squares
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
class Found(Exception):pass

def run(seed_file=None,extra_aspheres=False):
 tag='candidate_v21' if extra_aspheres else 'candidate_v20'
 prefix='two_conjugates_extra' if extra_aspheres else 'two_conjugates_high_quad'
 s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses']);scales=np.tile([1e-5,1e-7,1e-9],6)
 for a in s['aspheres']:
  j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
 cfg=next(x for x in json.loads((R/'candidate_v19_three_group_focus_targets.json').read_text(encoding='utf8')) if x['wd']==200.)
 near=t.copy();near[-1]=cfg['backfocus'];near[5]=cfg['iris_gap']
 x=np.r_[c[e.CURVED]/.05,t,c[15:]/scales,near[-1]-t[-1],near[5]-t[5],(cfg['gap']-e.focus_group_gap(c,near,200.))/.01,10.20]
 if seed_file is not None:
  x=np.load(R/seed_file)
  if x.shape!=(51,) and not(extra_aspheres and x.shape==(57,)):raise ValueError('Unexpected coupled prescription/focus parameter count')
 if extra_aspheres and len(x)==51:x=np.r_[x,np.zeros(6)]
 caps0=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 def unpack(x):
  cc=np.zeros(39 if extra_aspheres else 33);cc[e.CURVED]=x[:14]*.05;cc[15:33]=x[29:47]*scales
  if extra_aspheres:cc[33:]=x[51:57]*np.tile([1e-5,1e-7,1e-9],2)
  tt=x[14:29].copy();caps=caps0.copy();caps[9:11]=x[50]
  tn=tt.copy();tn[-1]+=x[47];tn[5]+=x[48];tn[10]=e.focus_group_gap(cc,tn,200.)+x[49]*.01
  pp=e.paraxial(cc,tt);rr,hh,ok=e.ray(cc,tt,e.N[1],500,0,0,pp['efl']/5.6,pp['ep'])
  return cc,tt,tn,caps,abs(rr[9])
 lastx=None;last=None;calls=0;best=np.inf;success=None
 def evaluate(x):
  nonlocal lastx,last,calls,best,success
  if lastx is not None and np.array_equal(x,lastx):return last
  cc,tt,tn,caps,stop=unpack(x)
  mb,rb,hb,ib=f.quantities(cc,tt,fields=np.array([0.,3.,6.3,9.08]),fixed_stop_radius=stop,low_frequency=False)
  mn,rn,hn,inn=f.quantities(cc,tn,fields=np.array([0.,6.3,9.08]),wd=200.,frequencies=(150.,120.),fixed_stop_radius=stop,low_frequency=False)
  pb=e.paraxial(cc,tt)
  glassmargin=np.where(e.MEDIA[:14]>=0,.05,0.)
  constraints=np.r_[f.blank_clearance(cc,tt,caps)+glassmargin,f.blank_clearance(cc,tn,caps)+glassmargin,
   caps-rb,caps-rn,.003*np.array([3,6.3,9.08])-abs(hb[1:]-np.array([3,6.3,9.08])),.0058*np.array([6.3,9.08])-abs(hn[1:]-np.array([6.3,9.08])),
   .1*(.02-abs(pb['efl']-35)),.1*(54.5-sum(tt)),.1*(59.8-sum(tn)),-.01*(ib+inn),
   1-cc[:15]**2*caps*caps-.01,2*(mb[0]-.305),2*(mb[-1]-.210),2*(mn[0]-.305),2*(mn[-1]-.215),2*(mb[1:-1]-.15).ravel(),2*(mn[1]-.15)]
  ti=tt.copy();ti[-1]+=5.;ti[10]=e.focus_group_gap(cc,ti,1e7)+tt[10]-e.focus_group_gap(cc,tt,500.)
  pi=e.paraxial(cc,ti,1e7);oy=9.08/pi['magnification'];ch=e.aim_chief(cc,ti,e.N[1],1e7,oy,pi['ep']);ray,hh,ok=e.ray(cc,ti,e.N[1],1e7,oy,0,ch,pi['ep'])
  cra=abs(np.rad2deg(np.arctan(ray[3])));constraints=np.r_[constraints,.1*(9.95-cra),f.blank_clearance(cc,ti,caps)+glassmargin]
  targetb=np.array([[.34,.34],[.22,.28],[.23,.28],[.25,.25]])
  optical_residual=np.r_[((mb-targetb)*np.array([[1,1],[.3,.3],[.3,.3],[2,2]])).ravel(),((mn-np.array([[.34,.34],[.23,.25],[.25,.25]]))*np.array([[1,1],[.5,.5],[2,2]])).ravel()]
  obj=float(np.dot(optical_residual,optical_residual))
  score=obj+100*np.dot(np.minimum(constraints,0),np.minimum(constraints,0));calls+=1
  if score<best:best=score;np.save(R/(prefix+'_best_parameters.npy'),x)
  if calls%100==0:print('CALL',calls,'BEST',best,'MTF500',mb.tolist(),'MTF200',mn.tolist(),'MIN',float(constraints.min()),'EFL',pb['efl'],'TTL200',sum(tn),'CRAinf',cra,flush=True)
  if constraints.min()>-1e-5:
   success=x.copy();raise Found()
  lastx=x.copy();last=(obj,constraints,np.r_[optical_residual,20*np.minimum(constraints,0)]);return last
 delta=np.r_[np.full(14,.20),np.full(15,.8),np.full(18,6.),np.array([4.,2.,10.,1.3])]
 if extra_aspheres:delta=np.r_[delta,np.full(6,6.)]
 lo=x-delta;hi=x+delta
 lo[14:29]=np.maximum(lo[14:29],np.array([1,.1,1,.65,.3,.3,.65,1,.1,1,.15,1,.1,.65,8]))
 lo[47:51]=[-1,-2,-5,10.1943];hi[47:51]=[5,2.5,16,11.49]
 jacx=None;jacval=None
 def jacobian(xx):
  nonlocal jacx,jacval
  if jacx is None or not np.array_equal(xx,jacx):
   jacx=xx.copy();jacval=approx_derivative(lambda z:np.r_[evaluate(z)[0],evaluate(z)[1]],xx,method='3-point',abs_step=1e-4*(1+abs(xx)),bounds=(lo,hi))
  return jacval
 try:
  initial=least_squares(lambda z:evaluate(z)[2],x,jac=lambda z:approx_derivative(lambda zz:evaluate(zz)[2],z,method='3-point',abs_step=1e-4*(1+abs(z)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',max_nfev=90,ftol=1e-8,xtol=1e-9,gtol=1e-8,verbose=1)
  x=initial.x;print('LEAST_SQUARES',initial.cost,initial.nfev,flush=True)
  opt=minimize(lambda z:evaluate(z)[0],x,method='SLSQP',jac=lambda z:jacobian(z)[0],bounds=list(zip(lo,hi)),constraints=[dict(type='ineq',fun=lambda z:evaluate(z)[1],jac=lambda z:jacobian(z)[1:])],options=dict(maxiter=100,ftol=1e-9,disp=True))
  final=opt.x;print('SOLVER',opt.success,opt.message,flush=True)
 except Found:final=success;print('FEASIBLE_FOUND',flush=True)
 cc,tt,tn,caps,stop=unpack(final)
 s['curvatures']=cc[:15].tolist();s['thicknesses']=tt.tolist();s['fixed_stop_radius']=float(stop);s['ray_aiming']='Real'
 surfaces=f.SURFS+([2,9] if extra_aspheres else [])
 s['aspheres']=[dict(surface=a,coefficients=cc[15+3*j:18+3*j].tolist()) for j,a in enumerate(surfaces)]
 for i in [9,10]:s['fixed_semi_diameters'][str(i+1)]=float(caps[i])
 s['focus_architecture']='Front three + STOP, middle three, rear two; fixed sensor and fixed physical iris'
 (R/(tag+'_schema.json')).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8');np.save(R/(prefix+'_final_parameters.npy'),final)
 cfgs=[dict(wd=wd,gap=float(tg[10]),iris_gap=float(tg[5]),backfocus=float(tg[-1]),frequency_center=fc,frequency_edge=fe,fixed_stop_radius=float(stop)) for wd,tg,fc,fe in [(500,tt,225,180),(200,tn,150,120)]]
 (R/(tag+'_focus_targets.json')).write_text(json.dumps(cfgs,indent=2),encoding='utf8');print('DONE',json.dumps(cfgs),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--seed-file');a.add_argument('--extra-aspheres',action='store_true');q=a.parse_args();run(q.seed_file,q.extra_aspheres)
