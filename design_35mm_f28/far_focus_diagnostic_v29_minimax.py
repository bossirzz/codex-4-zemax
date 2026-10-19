"""Strict max-min follow-up of frozen v29 far-field floating-focus search."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import JointModel

R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
model=JointModel('candidate_v29');c,t,n,caps,stop=model.unpack(model.x)
records=json.loads((R/'far_focus_diagnostic_v29_native_focus_configs.json').read_text(encoding='utf8'))
positions=np.array([[q['iris_gap']+q['gap']+q['backfocus'],q['gap']+q['backfocus'],q['backfocus']] for q in records if q['wd']!=2000.])
wd=2000.;cache={};calls=0

def evaluate(z):
 global calls
 key=tuple(z[:3])
 if key in cache:return cache[key]
 tt=t.copy();tt[-1],tt[5]=z[:2];tt[10]=e.focus_group_gap(c,tt,wd)+z[2]
 p=e.paraxial(c,tt,wd);ray,_,ok=e.ray(c,tt,e.N[1],wd,0,0,p['efl']/5.6,p['ep']);iris=abs(ray[9])
 mtf,r,h,invalid=f.quantities(c,tt,wd=wd,frequencies=(225.,180.),fixed_stop_radius=iris,low_frequency=False,caps=caps)
 pos=np.array([tt[5]+tt[10]+tt[-1],tt[10]+tt[-1],tt[-1]])
 strokes=np.ptp(np.vstack([positions,pos]),axis=0)
 fields=np.array([3.,6.3,9.08]);dist=100*(h[1:]-fields)/fields
 con=np.r_[f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.),(60.-sum(tt))*.1,(6.01-strokes)*.1,(tt[-1]-8.)*.1,(h[-1]-9)*.1,(.6-np.abs(dist))*.1,-invalid*.1]
 d=dict(tt=tt,mtf=mtf,constraints=con,stop=iris,strokes=strokes,distortion=dist,efl=p['efl']);cache[key]=d;calls+=1
 if calls%100==0:print('PROGRESS',calls,float(mtf.min()),float(con.min()),z.tolist(),flush=True)
 return d

solutions=[]
for x in [[13.59612366193886,6.5842860339106775,.019027318844570495],[13.688271563853245,7.072966924542362,.03322016361187108]]:
 z=np.r_[x,evaluate(np.asarray(x))['mtf'].min()]
 opt=minimize(lambda zz:-zz[3]-.002*evaluate(zz)['mtf'].mean(),z,method='SLSQP',
  bounds=[(positions[:,-1].max()-6.01,positions[:,-1].min()+6.01),(.3,10),(-.25,.3),(0.,.6)],
  constraints=[dict(type='ineq',fun=lambda zz:np.r_[evaluate(zz)['constraints'],evaluate(zz)['mtf'].ravel()-zz[3]])],
  options=dict(maxiter=60,eps=1e-5,ftol=1e-9))
 d=evaluate(opt.x);solutions.append((opt.x,d));print('LOCAL',opt.success,opt.message,opt.x.tolist(),d['mtf'].tolist(),float(d['constraints'].min()),flush=True)
feasible=[(z,d) for z,d in solutions if d['constraints'].min()>=-1e-5]
z,d=max(feasible or solutions,key=lambda v:v[1]['mtf'].min());tt=d['tt']
row=dict(wd=2000.,gap=float(tt[10]),iris_gap=float(tt[5]),backfocus=float(tt[-1]),fixed_stop_radius=float(d['stop']),frequency_center=225.,frequency_edge=180.,nominal_fno_target=2.8,
 predicted_mtf_225_180=d['mtf'].tolist(),ttl_mm=float(tt.sum()),predicted_efl=d['efl'],distortion_percent=d['distortion'].tolist(),three_group_stroke_mm=d['strokes'].tolist(),minimum_constraint=float(d['constraints'].min()),note='Strict max-min numerical prediction for frozen optics. Physical pupil clipping retained. Native verification required.')
(R/'far_focus_diagnostic_v29_minimax_result.json').write_text(json.dumps(row,indent=2),encoding='utf8')
records=[q for q in records if q['wd']!=2000.]+[row]
(R/'far_focus_diagnostic_v29_minimax_configs.json').write_text(json.dumps(records,indent=2),encoding='utf8')
print('FINAL',json.dumps(row),flush=True)
