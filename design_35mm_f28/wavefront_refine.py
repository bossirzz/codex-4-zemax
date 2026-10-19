"""Smooth joint wavefront refinement of eight real-glass aspheric surfaces."""
import argparse,json,sys,time
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import JointModel,R
sys.stdout.reconfigure(encoding='utf8')

def run(tag='candidate_v25',steps=120,seed_file=None,replacements=()):
 model=JointModel()
 if replacements:model.substitute([(int(k),v) for k,v in (x.split('=',1) for x in replacements)])
 seed=model.x if seed_file is None else np.load(R/seed_file)
 pupil=np.vstack([e.pupil_sample(5,16),[(np.cos(q),np.sin(q)) for q in np.arange(32)*np.pi/16]])
 calls=0;best=float('inf');bestx=seed.copy();start=time.time()
 def residual(x):
  nonlocal calls,best,bestx
  c,t,n,caps,stop=model.unpack(x);parts=[];wfes=[];radii=np.zeros(15);chiefs=[];bad=0
  for wd,tt in [(500.,t),(200.,n)]:
   p=e.paraxial(c,tt,wd)
   sr=abs(p['stop_scale'])*p['efl']/5.6
   hh=[]
   for field in [0.,3.,6.3,9.08]:
    data=[]
    for ni in e.N:
     ray,h,valid=e.trace_bundle_stop_aiming(c,tt,ni,wd,field/p['magnification'],p['ep'],pupil,p['efl']/5.6,sr)
     radii=np.maximum(radii,h.max(0));data.append((ray,valid));bad+=int((~valid).sum())
    ref=data[1][0][data[1][1],:2].mean(0) if data[1][1].any() else np.zeros(2)
    for wave,(ray,valid) in zip(e.WAVES,data):
     vz=ray[:,11];vx=ray[:,2]*vz;vy=ray[:,3]*vz
     phase=ray[:,4]-(sum(tt)-ray[:,7])/np.maximum(vz,1e-8)+vx*(ref[0]-ray[:,5])+vy*(ref[1]-ray[:,6])+vz*(sum(tt)-ray[:,7])
     phase=(phase-(phase[valid].mean() if valid.any() else 0))/(wave*.001)
     weight=1.0 if field==0 else (1.5 if field>=9 else .8)
     if wd==200 and field==0:weight=.7
     parts.extend([np.where(valid,phase,0.)*weight/np.sqrt(len(pupil)),(~valid).astype(float)*10])
     wfes.append(float(np.std(phase[valid])) if valid.any() else 100.)
    oy=field/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep'])
    ray,_,ok=e.ray(c,tt,e.N[1],wd,oy,0.,ch,p['ep']);hh.append(ray[1])
   chiefs.append(np.array(hh))
  glass=np.where(e.MEDIA[:14]>=0,.05,0.)
  geometry=np.r_[100*(f.blank_clearance(c,t,caps)+glass-.005),100*(f.blank_clearance(c,n,caps)+glass-.005),
   50*(caps-radii-.01),50*(np.array([3,6.3,9.08])*.0028-abs(chiefs[0][1:]-np.array([3,6.3,9.08]))),
   50*(np.array([3,6.3,9.08])*.0057-abs(chiefs[1][1:]-np.array([3,6.3,9.08]))),
   20*(.04-abs(e.paraxial(c,t)['efl']-35.)),2*(59.8-sum(t)),2*(59.8-sum(n)),1-c[:15]**2*caps**2-.01]
  ti=t.copy();ti[-1]+=5.;ti[10]=e.focus_group_gap(c,ti,1e7)+t[10]-e.focus_group_gap(c,t,500.)
  pi=e.paraxial(c,ti,1e7);ch=e.aim_chief(c,ti,e.N[1],1e7,9.08/pi['magnification'],pi['ep']);ray,_,ok=e.ray(c,ti,e.N[1],1e7,9.08/pi['magnification'],0.,ch,pi['ep'])
  cra=abs(np.rad2deg(np.arctan(ray[3])))
  geometry=np.r_[geometry,3*(9.95-cra),100*(f.blank_clearance(c,ti,caps)+glass-.005)]
  parts.extend([np.minimum(geometry,0.),.001*(x-seed)])
  res=np.concatenate(parts);score=float(res@res);calls+=1
  if score<best:
   best=score;bestx=x.copy();np.save(R/(tag+'_wavefront_best.npy'),x)
  if calls%200==0:
   print(json.dumps(dict(calls=calls,seconds=time.time()-start,best=best,wfe_waves=wfes,
    efl=e.paraxial(c,t)['efl'],cra_infinity=cra,ttl200=float(sum(n)),min_geometry=float(geometry.min()))),flush=True)
  return res
 delta=np.r_[np.full(14,.30),np.full(15,1.2),np.full(18,8.),[3.,1.8,12.,1.2],np.full(6,8.)]
 lo=seed-delta;hi=seed+delta
 lo[14:29]=np.maximum(lo[14:29],[1,.1,1,.45,.3,.3,.45,1,.1,1,.15,1,.1,.45,8])
 lo[47:51]=[-1,-2,-5,10.05];hi[47:51]=[5,2.5,16,11.45]
 opt=least_squares(residual,seed,jac=lambda z:approx_derivative(residual,z,method='3-point',
  abs_step=1e-5*(1+abs(z)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',
  max_nfev=steps,ftol=1e-8,xtol=1e-9,gtol=1e-7)
 model.export(bestx,tag)
 c,t,n,caps,stop=model.unpack(bestx)
 print('DONE',opt.cost,opt.nfev,'MTF500',f.quantities(c,t,fixed_stop_radius=stop,low_frequency=False)[0].tolist(),
  'MTF200',f.quantities(c,n,wd=200,frequencies=(150,120),fixed_stop_radius=stop,low_frequency=False)[0].tolist(),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--tag',default='candidate_v25');a.add_argument('--steps',type=int,default=120);a.add_argument('--seed-file');a.add_argument('--replace',action='append',default=[]);q=a.parse_args();run(q.tag,q.steps,q.seed_file,q.replace)
