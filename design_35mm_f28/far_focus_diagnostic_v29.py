"""Frozen v29 optics: physical-pupil conjugate/floating-focus diagnostics.

Numerical predictions only. No ZOS API. Existing candidates are never changed.
"""
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import SURFACES
from optimize_best_conjugate import load_glasses

R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(wds=(2000.,300.,800.),nstarts=5,iterations=100):
 tag='candidate_v29';prefix='far_focus_diagnostic_v29'
 schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
 native=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in schema['glasses']])
 c=np.r_[schema['curvatures'],np.zeros(24)];t=np.array(schema['thicknesses'])
 for a in schema['aspheres']:
  j=15+3*SURFACES.index(a['surface']);c[j:j+3]=a['coefficients']
 caps=np.array([schema['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 best=next(q for q in native if q['wd']==500.)
 t[5],t[10],t[-1]=best['iris_gap'],best['gap'],best['backfocus']
 baseline_positions=np.array([[q['iris_gap']+q['gap']+q['backfocus'],q['gap']+q['backfocus'],q['backfocus']] for q in native])
 results=[];configs=[dict(wd=q['wd'],iris_gap=q['iris_gap'],gap=q['gap'],backfocus=q['backfocus'],frequency_center=150. if q['wd']==200. else 225.,frequency_edge=120. if q['wd']==200. else 180.,fixed_stop_radius=schema['fixed_stop_radius'],nominal_fno_target=2.8) for q in native]
 for wd in wds:
  started=time.time();cache={};calls=0
  q=next((q for q in native if q['wd']==wd),best)
  tt0=t.copy();tt0[5],tt0[10],tt0[-1]=q['iris_gap'],q['gap'],q['backfocus']
  if q['wd']!=wd:tt0[10]=e.focus_group_gap(c,tt0,wd)+(t[10]-e.focus_group_gap(c,t,500.))
  seed=np.array([tt0[-1],tt0[5],tt0[10]-e.focus_group_gap(c,tt0,wd)])
  def evaluate(x):
   nonlocal calls
   key=tuple(x)
   if key in cache:return cache[key]
   tt=t.copy();tt[-1],tt[5]=x[0],x[1];tt[10]=e.focus_group_gap(c,tt,wd)+x[2]
   p=e.paraxial(c,tt,wd);mr,_,ok=e.ray(c,tt,e.N[1],wd,0.,0.,p['efl']/5.6,p['ep']);stop=abs(mr[9])
   mtf,r,h,invalid,diag=f.quantities(c,tt,wd=wd,frequencies=(225.,180.),fixed_stop_radius=stop,low_frequency=False,caps=caps,return_diagnostics=True)
   pos=np.array([tt[5]+tt[10]+tt[-1],tt[10]+tt[-1],tt[-1]])
   strokes=np.ptp(np.vstack([baseline_positions,pos]),axis=0)
   clearance=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
   fields=np.array([3.,6.3,9.08]);dist=100.*(h[1:]-fields)/fields
   # Only best and 200 mm have specified distortion thresholds. Use 0.6%
   # elsewhere as an explicit engineering guard; do not invent acceptance.
   limit=.6 if wd!=500. else .3
   con=np.r_[clearance,(60.-sum(tt))*.1,(6.01-strokes)*.1,(tt[-1]-8.)*.1,
     (h[-1]-9.)*.1,(limit-np.abs(dist))*.1,-invalid*.1]
   softmin=-.018*logsumexp(-mtf.ravel()/.018)
   objective=-softmin-.025*float(mtf.mean())
   score=objective+2000.*float(np.square(np.minimum(con,0)).sum())
   value=dict(tt=tt,mtf=mtf,constraints=con,score=score,objective=objective,stop=stop,diagnostics=diag,distortion=dist,strokes=strokes,efl=p['efl'])
   cache[key]=value;calls+=1
   if calls%100==0:print(json.dumps(dict(wd=wd,calls=calls,seconds=time.time()-started,min_mtf=float(mtf.min()),min_constraint=float(con.min()),x=x.tolist())),flush=True)
   return value
  base=evaluate(seed)
  print('BASELINE',wd,json.dumps(dict(x=seed.tolist(),mtf=base['mtf'].tolist(),min_constraint=float(base['constraints'].min()))),flush=True)
  # Diagnose the best possible focus of each field without moving rear/iris.
  sweep=[]
  for defocus in np.linspace(seed[2]-.12,seed[2]+.12,41):
   x=seed.copy();x[2]=defocus;d=evaluate(x)
   mono=f.quantities(c,d['tt'],wd=wd,frequencies=(225.,180.),fixed_stop_radius=d['stop'],low_frequency=False,wavelengths=[e.WAVES[1]],caps=caps)[0]
   sweep.append(dict(defocus_mm=float(defocus),poly_mtf=d['mtf'].tolist(),mono_d_mtf=mono.tolist(),min_constraint=float(d['constraints'].min())))
  (R/(prefix+'_scan_'+str(int(wd))+'.json')).write_text(json.dumps(sweep,indent=2),encoding='utf8')
  feasible_sweep=[z for z in sweep if z['min_constraint']>=-1e-6]
  scanbest=max(feasible_sweep or sweep,key=lambda z:np.min(z['poly_mtf']))
  ss=seed.copy();ss[2]=scanbest['defocus_mm']
  lo=np.array([max(8.,baseline_positions[:,-1].max()-6.01),.3,-.25])
  hi=np.array([baseline_positions[:,-1].min()+6.01,10.,.30])
  seeds=[seed,ss]
  for back in np.linspace(lo[0]+.1,hi[0]-.1,4):
   for iris in [max(.4,seed[1]-2.),seed[1],min(9.8,seed[1]+2.)]:
    for defocus in [seed[2]-.04,seed[2],seed[2]+.04]:seeds.append(np.array([back,iris,defocus]))
  ranked=sorted(seeds,key=lambda x:evaluate(x)['score'])
  print('GRID_BEST',wd,[(z.tolist(),float(evaluate(z)['mtf'].min()),float(evaluate(z)['constraints'].min())) for z in ranked[:nstarts]],flush=True)
  solutions=[ss,seed]
  for i,x in enumerate(ranked[:nstarts]):
   opt=minimize(lambda z:evaluate(z)['score'],x,method='Nelder-Mead',bounds=list(zip(lo,hi)),options=dict(maxiter=iterations,xatol=1e-5,fatol=2e-7))
   solutions.append(opt.x)
   strict=minimize(lambda z:evaluate(z)['objective'],opt.x,method='SLSQP',bounds=list(zip(lo,hi)),constraints=[dict(type='ineq',fun=lambda z:evaluate(z)['constraints'])],options=dict(maxiter=40,eps=1e-5,ftol=1e-9))
   solutions.append(strict.x)
   print('LOCAL',wd,i,json.dumps(dict(x=strict.x.tolist(),mtf=evaluate(strict.x)['mtf'].tolist(),min_constraint=float(evaluate(strict.x)['constraints'].min()))),flush=True)
  feasible=[x for x in solutions if evaluate(x)['constraints'].min()>=-1e-5]
  x=min(feasible or solutions,key=lambda z:evaluate(z)['score']);d=evaluate(x);tt=d['tt']
  mono=f.quantities(c,tt,wd=wd,frequencies=(225.,180.),fixed_stop_radius=d['stop'],low_frequency=False,wavelengths=[e.WAVES[1]],caps=caps)[0]
  lower=f.quantities(c,tt,wd=wd,frequencies=(150.,120.),fixed_stop_radius=d['stop'],low_frequency=False,caps=caps)[0]
  row=dict(wd=wd,gap=float(tt[10]),iris_gap=float(tt[5]),backfocus=float(tt[-1]),frequency_center=225.,frequency_edge=180.,fixed_stop_radius=float(d['stop']),nominal_fno_target=2.8,
   predicted_mtf_225_180=d['mtf'].tolist(),predicted_mono_d_mtf_225_180=mono.tolist(),predicted_mtf_150_120=lower.tolist(),baseline_predicted_mtf_225_180=base['mtf'].tolist(),
   predicted_efl=d['efl'],ttl_mm=float(tt.sum()),distortion_percent=d['distortion'].tolist(),three_group_stroke_mm=d['strokes'].tolist(),minimum_constraint=float(d['constraints'].min()),
   predicted_relative_illumination_d=d['diagnostics']['relative_illumination_primary'][-1],note='Frozen optics numerical prediction; physical lens pupil clipping retained. Native verification required. 2m/300/800 thresholds have not been user specified; 150/120 is a supplemental engineering record only.')
  results.append(row)
  configs=[z for z in configs if z['wd']!=wd]+[row]
  (R/(prefix+'_results.json')).write_text(json.dumps(results,indent=2),encoding='utf8')
  (R/(prefix+'_native_focus_configs.json')).write_text(json.dumps(configs,indent=2),encoding='utf8')
  print('RESULT',json.dumps(row),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--wds',type=float,nargs='+',default=[2000.,300.,800.]);p.add_argument('--starts',type=int,default=5);p.add_argument('--iterations',type=int,default=100);a=p.parse_args();run(a.wds,a.starts,a.iterations)
