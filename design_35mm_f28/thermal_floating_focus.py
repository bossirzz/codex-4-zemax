"""Fit temperature compensation using native thermal geometry and indices."""
import argparse,json,sys
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R,SURFACES
sys.stdout.reconfigure(encoding='utf8')

def run(tag):
 cases=json.loads((R/(tag+'_thermal_snapshots.json')).read_text(encoding='utf8'))
 schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 native=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
 nominal=np.array(schema['thicknesses']);positions=[];out=[]
 for row in native:
  tt=nominal.copy();tt[5],tt[10],tt[-1]=row['iris_gap'],row['gap'],row['backfocus']
  positions.append(np.array([sum(tt),sum(tt[6:]),sum(tt[11:])]))
 for case in cases:
  wd=case['wd'];temp=case['temperature_c'];c=np.r_[case['curvatures'],np.zeros(24)]
  for a in case['aspheres']:
   j=15+3*SURFACES.index(a['surface']);c[j:j+3]=a['coefficients']
  t=np.array(case['thicknesses']);caps=np.array(case['caps']);e.N=np.array(case['indices'])
  offset=(t[10]-e.focus_group_gap(c,t,wd))/.01
  fc,fe=(225.,180.) if wd==500 else (150.,120.)
  cache={}
  def evaluate(x):
   key=tuple(x)
   if key in cache:return cache[key]
   tt=t.copy();tt[-1]+=x[0];tt[5]+=x[1];tt[10]=e.focus_group_gap(c,tt,wd)+x[2]*.01
   p=e.paraxial(c,tt,wd);stop=abs(p['stop_scale'])*p['efl']/5.6
   mtf,r,h,bad,diag=f.quantities(c,tt,wd=wd,frequencies=(fc,fe),fixed_stop_radius=stop,
     low_frequency=False,caps=caps,return_diagnostics=True)
   pos=np.array([sum(tt),sum(tt[6:]),sum(tt[11:])]);stroke=np.ptp(np.vstack([*positions,pos]),axis=0)
   # Allow temperature variation of the nominal manufacturing clearances;
   # the full material and air profiles must remain separated.
   clearance=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,.02)
   con=np.r_[clearance,.1*(60.-sum(tt)),.1*(6.01-stroke),.1*(tt[-1]-8),
    .1*(h[-1]-9),diag['relative_illumination_primary'][-1]-.80,-bad*.1]
   oy=9.08/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);ray,_,ok=e.ray(c,tt,e.N[1],wd,oy,0,ch,p['ep'])
   cra=abs(np.rad2deg(np.arctan(ray[3])))
   if wd>=1e7:
    con=np.r_[con,.1*(9.93-cra)];target=np.array([[.30,.30],[.20,.20],[.20,.20],[.20,.20]])
   else:
    con=np.r_[con,2*(mtf[0]-.308),2*(mtf[1:]-.210).ravel()]
    target=np.array([[.32,.32],[.215,.215],[.215,.215],[.215,.215]])
   objective=float(np.sum(np.maximum(target-mtf,0.)**2))
   score=objective+1000*np.dot(np.minimum(con,0),np.minimum(con,0))
   d=dict(score=score,objective=objective,constraints=con,t=tt,mtf=mtf,stop=stop,cra=cra,
     positions=pos,diagnostics=diag)
   cache[key]=d;return d
  bounds=[(-2.,2.),(-1.2,1.2),(-10.,20.)]
  if abs(temp-20.)<1e-6 or wd>=1e7:
   x=np.array([0.,0.,offset])
  else:
   seeds=[np.array([b,g,d]) for b in [-1.,0.,1.] for g in [-.6,0.,.6] for d in [offset-2.,offset,offset+2.]]
   seeds=[np.clip(x,[-2,-1.2,-10],[2,1.2,20]) for x in seeds]
   ranked=sorted(seeds,key=lambda x:evaluate(x)['score']);trials=[]
   for seed in ranked[:2]:
    a=minimize(lambda z:evaluate(z)['score'],seed,method='Nelder-Mead',bounds=bounds,
      options=dict(maxiter=120,xatol=2e-4,fatol=1e-8));trials.append(a.x)
    b=minimize(lambda z:evaluate(z)['objective'],a.x,method='SLSQP',bounds=bounds,
      constraints=[dict(type='ineq',fun=lambda z:evaluate(z)['constraints'])],
      options=dict(maxiter=50,eps=1e-4,ftol=1e-9));trials.append(b.x)
   feasible=[x for x in trials if evaluate(x)['constraints'].min()>=-1e-5]
   x=min(feasible or trials,key=lambda x:evaluate(x)['score'])
  d=evaluate(x);tt=d['t'];positions.append(d['positions'])
  record=dict(wd=wd,temperature_c=temp,spacer_tce_ppm_per_k=case['spacer_tce_ppm_per_k'],
   iris_gap=float(tt[5]),gap=float(tt[10]),backfocus=float(tt[-1]),stop_radius=float(d['stop']),
   frequency_center=fc,frequency_edge=fe,predicted_mtf=d['mtf'].tolist(),predicted_cra=d['cra'],
   predicted_illumination=d['diagnostics']['relative_illumination_primary'][-1],
   min_constraint=float(d['constraints'].min()),group_positions_from_sensor_mm=d['positions'].tolist(),
   focus_corrections_mm=[float(tt[-1]-t[-1]),float(tt[5]-t[5]),float(tt[10]-t[10])],
   status='PREDICTION_ONLY_REQUIRES_NATIVE_VALIDATION')
  out.append(record);print('THERMAL_FOCUS',json.dumps(record),flush=True)
  (R/(tag+'_thermal_focus_targets.json')).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args();run(q.tag)
