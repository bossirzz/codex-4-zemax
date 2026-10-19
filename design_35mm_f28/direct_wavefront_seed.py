"""Find a smooth wavefront seed before exact-MTF constrained refinement."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run():
 s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'));load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses']);scales=np.tile([1e-5,1e-7,1e-9],6)
 for a in s['aspheres']:
  j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
 cfg=next(x for x in json.loads((R/'candidate_v19_three_group_focus_targets.json').read_text(encoding='utf8')) if x['wd']==200)
 near=t.copy();near[-1]=cfg['backfocus'];near[5]=cfg['iris_gap']
 x=np.r_[c[e.CURVED]/.05,t,c[15:]/scales,near[-1]-t[-1],near[5]-t[5],(cfg['gap']-e.focus_group_gap(c,near,200))/.01,10.20]
 caps0=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 pupil=np.vstack([e.pupil_sample(5,16),[(np.cos(q),np.sin(q)) for q in np.arange(24)*np.pi/12]])
 def unpack(x):
  cc=np.zeros(33);cc[e.CURVED]=x[:14]*.05;cc[15:]=x[29:47]*scales;tt=x[14:29].copy();caps=caps0.copy();caps[9:11]=x[50]
  tn=tt.copy();tn[-1]+=x[47];tn[5]+=x[48];tn[10]=e.focus_group_gap(cc,tn,200)+x[49]*.01
  p=e.paraxial(cc,tt);r,hh,ok=e.ray(cc,tt,e.N[1],500,0,0,p['efl']/5.6,p['ep'])
  return cc,tt,tn,caps,abs(r[9])
 calls=0;best=np.inf
 def residual(x):
  nonlocal calls,best
  cc,tt,tn,caps,stop=unpack(x);pb=e.paraxial(cc,tt);pieces=[];maxr=np.zeros(15);wfes=[];chiefs=[]
  for wd,tg,fields in [(500,tt,[0,3,6.3,9.08]),(200,tn,[0,6.3,9.08])]:
   p=e.paraxial(cc,tg,wd);h=[]
   for field in fields:
    data=[]
    for n in e.N:
     r,hh,v=e.trace_bundle_stop_aiming(cc,tg,n,wd,field/p['magnification'],p['ep'],pupil,pb['efl']/5.6,stop)
     maxr=np.maximum(maxr,hh.max(0));data.append((r,v))
    ref=data[1][0][data[1][1],:2].mean(0) if data[1][1].any() else np.zeros(2)
    for wave,(r,v) in zip(e.WAVES,data):
     vz=r[:,11];vx=r[:,2]*vz;vy=r[:,3]*vz
     phase=r[:,4]-(sum(tg)-r[:,7])/np.maximum(vz,1e-8)+vx*(ref[0]-r[:,5])+vy*(ref[1]-r[:,6])+vz*(sum(tg)-r[:,7])
     phase=(phase-(phase[v].mean() if v.any() else 0))*v/(wave*.001)
     weight=(2. if field==0 else (1.5 if field>=9 else .7))
     pieces.extend([phase*weight/np.sqrt(len(pupil)),(~v).astype(float)*10])
     wfes.append(float(np.std(phase[v])) if v.any() else 1e3)
    oy=field/p['magnification'];ch=e.aim_chief(cc,tg,e.N[1],wd,oy,p['ep']);r,hh,v=e.ray(cc,tg,e.N[1],wd,oy,0,ch,p['ep']);h.append(r[1])
   chiefs.append(np.asarray(h))
  glassmargin=np.where(e.MEDIA[:14]>=0,.05,0.)
  con=np.r_[f.blank_clearance(cc,tt,caps)+glassmargin,f.blank_clearance(cc,tn,caps)+glassmargin,caps-maxr,
   .003*np.array([3,6.3,9.08])-abs(chiefs[0][1:]-np.array([3,6.3,9.08])),.0058*np.array([6.3,9.08])-abs(chiefs[1][1:]-np.array([6.3,9.08])),
   .1*(.02-abs(pb['efl']-35)),.1*(54.5-sum(tt)),.1*(59.8-sum(tn)),1-cc[:15]**2*caps*caps-.01]
  ti=tt.copy();ti[-1]+=5.;ti[10]=e.focus_group_gap(cc,ti,1e7)+tt[10]-e.focus_group_gap(cc,tt,500)
  pi=e.paraxial(cc,ti,1e7);oy=9.08/pi['magnification'];ch=e.aim_chief(cc,ti,e.N[1],1e7,oy,pi['ep']);ray,hh,ok=e.ray(cc,ti,e.N[1],1e7,oy,0,ch,pi['ep'])
  con=np.r_[con,.1*(9.95-abs(np.rad2deg(np.arctan(ray[3])))),f.blank_clearance(cc,ti,caps)+glassmargin]
  pieces.append(np.minimum(con,0)*40);res=np.concatenate(pieces);score=res@res;calls+=1
  if score<best:best=score;np.save(R/'wavefront_best_parameters.npy',x)
  if calls%200==0:print('CALL',calls,'BEST',best,'WFE',wfes,'EFL',pb['efl'],'TTL200',sum(tn),'MIN_GEOMETRY',float(con.min()),flush=True)
  return res
 delta=np.r_[np.full(14,.25),np.full(15,1.),np.full(18,6.),np.array([4,2,10,1.3])]
 lo=x-delta;hi=x+delta;lo[14:29]=np.maximum(lo[14:29],[1,.1,1,.65,.3,.3,.65,1,.1,1,.15,1,.1,.65,8])
 lo[47:]=[-1,-2,-5,10.1943];hi[47:]=[5,2.5,16,11.49]
 opt=least_squares(residual,x,jac=lambda z:approx_derivative(residual,z,method='3-point',abs_step=1e-4*(1+abs(z)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',max_nfev=100,ftol=1e-7,xtol=1e-8,gtol=1e-8)
 np.save(R/'wavefront_conditioned_seed.npy',opt.x)
 cc,tt,tn,caps,stop=unpack(opt.x)
 print('DONE',opt.cost,opt.nfev,'BEST_MTF',f.quantities(cc,tt,fixed_stop_radius=stop,low_frequency=False)[0].tolist(),'NEAR_MTF',f.quantities(cc,tn,fields=np.array([0,9.08]),wd=200,frequencies=(150,120),fixed_stop_radius=stop,low_frequency=False)[0].tolist(),flush=True)

if __name__=='__main__':run()
