import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--waves',type=float,nargs='+',default=[.65,.85]);q=a.parse_args()
s=json.loads((R/(q.tag+'_schema.json')).read_text(encoding='utf8'));load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for asp in s['aspheres']:
 j=15+3*f.SURFS.index(asp['surface']);c[j:j+3]=asp['coefficients']
offset=t[10]-e.focus_group_gap(c,t,500.);caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
records=json.loads((R/(q.tag+'_auxiliary_focus_profiles.json')).read_text(encoding='utf8')) if (R/(q.tag+'_auxiliary_focus_profiles.json')).exists() else []
for wave in q.waves:
 def evaluate(x):
  tt=t.copy();tt[-1]+=x[0];tt[10]=e.focus_group_gap(c,tt,500.)+x[1]*.01
  m,r,h,iv=f.quantities(c,tt,fields=np.array([0,9.08]),fixed_stop_radius=s['fixed_stop_radius'],low_frequency=False,wavelengths=[wave])
  con=np.r_[f.blank_clearance(c,tt,caps),caps-r,.003*9.08-abs(h[-1]-9.08)]
  if wave==.65:con=np.r_[con,m[0]-.305,m[-1]-.208]
  score=float(np.sum(((m-np.array([[.36,.36],[.24,.24]]))*np.array([[1,1],[2,2]]))**2))+1000*np.dot(np.minimum(con,0),np.minimum(con,0))
  return score,tt,m,con
 initial=evaluate([0,offset/.01]);seeds=[np.array([b,g]) for b in [-.5,0,.5,1,1.5] for g in np.linspace(-15,20,15)]
 ranked=sorted(seeds,key=lambda x:evaluate(x)[0]);trials=[]
 for seed in ranked[:3]:
  opt=minimize(lambda x:evaluate(x)[0],seed,method='Nelder-Mead',bounds=[(-1,2),(-20,25)],options=dict(maxiter=70,xatol=1e-3,fatol=1e-8));trials.append(opt.x)
 x=min(trials,key=lambda x:evaluate(x)[0]);score,tt,m,con=evaluate(x)
 rr=dict(wavelength_um=wave,visible_focus_mtf=initial[2].tolist(),refocused_mtf=m.tolist(),gap=float(tt[10]),backfocus=float(tt[-1]),min_constraint=float(con.min()))
 records=[r for r in records if r['wavelength_um']!=wave]+[rr];print(json.dumps(rr),flush=True)
 (R/(q.tag+'_auxiliary_focus_profiles.json')).write_text(json.dumps(records,indent=2),encoding='utf8')
