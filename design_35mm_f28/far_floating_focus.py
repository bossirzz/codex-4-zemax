"""Improve 2 m and infinity focus using declared three-group motion and clipping."""
import argparse,json,sys
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R,SURFACES
from optimize_best_conjugate import load_glasses
sys.stdout.reconfigure(encoding='utf8')

def run(tag):
 s=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 c=np.r_[s['curvatures'],np.zeros(24)];t=np.array(s['thicknesses'])
 for a in s['aspheres']:
  j=15+3*SURFACES.index(a['surface']);c[j:j+3]=a['coefficients']
 caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 records=json.loads((R/(tag+'_focus_targets.json')).read_text(encoding='utf8'))
 positions=[np.array([q['iris_gap']+q['gap']+q['backfocus'],q['gap']+q['backfocus'],q['backfocus']]) for q in records]
 offset=t[10]-e.focus_group_gap(c,t,500.)
 for wd in [1e7,2000.]:
  cache={}
  def evaluate(x):
   key=tuple(x)
   if key in cache:return cache[key]
   tt=t.copy();tt[-1]+=x[0];tt[5]+=x[1];tt[10]=e.focus_group_gap(c,tt,wd)+x[2]*.01
   p=e.paraxial(c,tt,wd);stop=abs(p['stop_scale'])*p['efl']/5.6
   mtf,r,h,invalid,diag=f.quantities(c,tt,wd=wd,frequencies=(150.,120.),fixed_stop_radius=stop,
    low_frequency=False,caps=caps,return_diagnostics=True)
   oy=9.08/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);ray,_,ok=e.ray(c,tt,e.N[1],wd,oy,0,ch,p['ep'])
   cra=abs(np.rad2deg(np.arctan(ray[3])));pos=np.array([tt[5]+tt[10]+tt[-1],tt[10]+tt[-1],tt[-1]])
   strokes=np.ptp(np.vstack([*positions,pos]),axis=0)
   con=np.r_[f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.),
    .1*(59.8-sum(tt)),.1*(6.-strokes),.1*(tt[-1]-8),.008*9.08-abs(h[-1]-9.08),
    .1*(h[-1]-9),diag['relative_illumination_primary'][-1]-.82,-invalid*.1]
   if wd>=1e7:con=np.r_[con,.1*(9.95-cra)]
   objective=float(np.sum((mtf-np.array([[.33,.33],[.22,.25],[.22,.25],[.24,.24]]))**2))
   score=objective+1000*np.dot(np.minimum(con,0),np.minimum(con,0))
   out=dict(score=score,objective=objective,constraints=con,tt=tt,mtf=mtf,cra=cra,positions=pos,diagnostics=diag,p=p)
   cache[key]=out;return out
  seeds=[np.array([b,g,d]) for b in [-2.,0.,2.5,5.] for g in [-1.5,0.,1.5,2.5] for d in [0.,5.,10.]]
  ranked=sorted(seeds,key=lambda x:evaluate(x)['score']);print('FAR_GRID',wd,[(x.tolist(),evaluate(x)['score']) for x in ranked[:3]],flush=True)
  solutions=[]
  for x in ranked[:3]:
   opt=minimize(lambda z:evaluate(z)['score'],x,method='Nelder-Mead',bounds=[(-2.5,5),(-2,2.5),(-10,20)],options=dict(maxiter=100,xatol=1e-3,fatol=1e-7))
   solutions.append(opt.x)
   strict=minimize(lambda z:evaluate(z)['objective'],opt.x,method='SLSQP',bounds=[(-2.5,5),(-2,2.5),(-10,20)],
    constraints=[dict(type='ineq',fun=lambda z:evaluate(z)['constraints'])],options=dict(maxiter=45,eps=1e-4,ftol=1e-8))
   solutions.append(strict.x)
  feasible=[x for x in solutions if evaluate(x)['constraints'].min()>=-1e-5]
  x=min(feasible or solutions,key=lambda x:evaluate(x)['score']);d=evaluate(x);tt=d['tt'];positions.append(d['positions'])
  row=dict(wd=wd,gap=float(tt[10]),iris_gap=float(tt[5]),backfocus=float(tt[-1]),frequency_center=150.,frequency_edge=120.,
   nominal_fno_target=2.8,fixed_stop_radius=s['fixed_stop_radius'],predicted_mtf=d['mtf'].tolist(),
   predicted_cra=d['cra'],predicted_illumination=d['diagnostics']['relative_illumination_primary'][-1],
   min_constraint=float(d['constraints'].min()),note='150/120 lp/mm at far distances is an engineering check, not a user-defined acceptance threshold')
  records.append(row);print('FAR_RESULT',json.dumps(row),flush=True)
  (R/(tag+'_floating_validation_configs.json')).write_text(json.dumps(records,indent=2),encoding='utf8')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args();run(q.tag)
