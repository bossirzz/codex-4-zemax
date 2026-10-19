"""Evaluate a third floating group without altering lens shapes or glasses."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,wds):
 s=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 surfaces=f.SURFS+[2,9]
 c=np.r_[s['curvatures'],np.zeros(24 if any(a['surface'] in [2,9] for a in s['aspheres']) else 18)];t=np.array(s['thicknesses'])
 for a in s['aspheres']:
  j=15+3*surfaces.index(a['surface']);c[j:j+3]=a['coefficients']
 caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 offset=t[10]-e.focus_group_gap(c,t,500.)
 baseline=np.array([t[5]+t[10]+t[14],t[10]+t[14],t[14]])
 records=[]
 for wd in wds:
  fc,fe=(150.,120.) if wd<=200 else (225.,180.)
  distlimit=.003 if 400<=wd<=500 else .006
  cache={}
  def evaluate(x):
   key=tuple(x)
   if key in cache:return cache[key]
   tt=t.copy();tt[-1]+=x[0];tt[5]+=x[1];tt[10]=e.focus_group_gap(c,tt,wd)+x[2]*.01
   m,r,chief,invalid=f.quantities(c,tt,fields=np.array([0.,9.08]),wd=wd,frequencies=(fc,fe),fixed_stop_radius=s['fixed_stop_radius'],low_frequency=False)
   p=e.paraxial(c,tt,wd);oy=9.08/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);ray,hh,ok=e.ray(c,tt,e.N[1],wd,oy,0,ch,p['ep'])
   cra=abs(np.rad2deg(np.arctan(ray[3])));positions=np.array([tt[5]+tt[10]+tt[-1],tt[10]+tt[-1],tt[-1]])
   clear=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
   con=np.r_[clear,.1*(59.8-sum(tt)),.1*(6.-abs(positions-baseline)),
    .1*(tt[-1]-8),distlimit*9.08-abs(chief[-1]-9.08),-invalid*.01]
   # Require clear rays where the near-distance image quality is evaluated.
   if wd<=200:con=np.r_[con,caps-r]
   if wd>=1e7:con=np.r_[con,.1*(9.97-cra)]
   target=np.array([[.34,.34],[.26,.26]])
   objective=float(np.sum(((m-target)*np.array([[1,1],[2,2]]))**2))
   score=objective+1000*np.dot(np.minimum(con,0),np.minimum(con,0))
   d=dict(score=score,objective=objective,constraints=con,tt=tt,m=m,r=r,chief=chief,cra=cra,positions=positions,p=p)
   cache[key]=d;return d
  if wd==500:x=np.array([0.,0.,offset/.01])
  else:
   seeds=[np.array([b,g,d]) for b in [-1.,.5,2.,3.5,5.] for g in [-1.5,-.75,0.,.75,1.5] for d in [0.,3.,6.,9.,12.]]
   ranked=sorted(seeds,key=lambda x:evaluate(x)['score'])
   print('GRID',wd,[(x.tolist(),evaluate(x)['score'],evaluate(x)['m'].tolist()) for x in ranked[:3]],flush=True)
   feasible=[x for x in ranked if evaluate(x)['constraints'].min()>=-1e-6]
   trials=[]
   for seed in ranked[:3]:
    opt=minimize(lambda x:evaluate(x)['score'],seed,method='Nelder-Mead',bounds=[(-1.,5.),(-2.,2.5),(-5.,16.)],options=dict(maxiter=100,xatol=2e-3,fatol=1e-7))
    trials.append(opt.x)
    strict=minimize(lambda x:evaluate(x)['objective'],opt.x,method='SLSQP',bounds=[(-1.,5.),(-2.,2.5),(-5.,16.)],constraints=[dict(type='ineq',fun=lambda x:evaluate(x)['constraints'])],options=dict(maxiter=45,eps=1e-3,ftol=1e-8))
    if evaluate(strict.x)['constraints'].min()>=-1e-6:feasible.append(strict.x)
   x=min(feasible or trials,key=lambda x:evaluate(x)['score'])
  d=evaluate(x);tt=d['tt'];row=dict(wd=wd,gap=float(tt[10]),gap_surface=11,iris_gap=float(tt[5]),backfocus=float(tt[-1]),frequency_center=fc,frequency_edge=fe,
   mtf=d['m'].tolist(),distortion_percent=float((d['chief'][-1]/9.08-1)*100),cra_deg=float(d['cra']),efl=d['p']['efl'],ttl=float(sum(tt)),
   min_constraint=float(d['constraints'].min()),group_positions=d['positions'].tolist(),max_ray_radii=d['r'].tolist(),fixed_stop_radius=s['fixed_stop_radius'],ray_aiming='Real')
  records.append(row);print('RESULT',json.dumps(row),flush=True)
  (R/(tag+'_three_group_focus_targets.json')).write_text(json.dumps(records,indent=2),encoding='utf8')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--wds',type=float,nargs='+',default=[500,200,300,400,450,800,2000,1e7]);q=a.parse_args();run(q.tag,q.wds)
