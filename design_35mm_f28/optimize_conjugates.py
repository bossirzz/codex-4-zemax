"""Suggest floating front-six/rear-two focus positions, pending native checks."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
def run():
 s=json.loads((ROOT/'candidate_v17_schema.json').read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
 for a in s['aspheres']:
  j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
 offset=t[10]-e.focus_group_gap(c,t,500.);caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 output=[]
 for wd in [200.,300.,400.,450.,500.,800.,2000.,1e7,100.]:
  fc,fe=(150.,120.) if wd<=200 else (225.,180.)
  def trial(x):
   tt=t.copy();tt[-1]+=x[0];tt[10]=e.focus_group_gap(c,tt,wd)+x[1]*.01
   return tt
  def fun(x):
   tt=trial(x);m,radii,chiefs,invalid=f.quantities(c,tt,fields=np.array([0.,9.08]),wd=wd,frequencies=(fc,fe))
   desired=np.array([[.38,.38],[.28,.28]]) if wd<=200 else np.array([[.36,.36],[.26,.26]])
   parts=[((m[:,:2]-desired)*np.array([[1.,1.],[1.,1.]])).ravel()]
   if wd!=100:parts.append(np.array([max(0,sum(tt)-59.5)*2,max(0,.3-tt[10])*10,max(0,8-tt[-1])*10]))
   if wd>1e6:
    p=e.paraxial(c,tt,wd);oy=9.08/p['magnification'];chief=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);r,_,ok=e.ray(c,tt,e.N[1],wd,oy,0.,chief,p['ep'])
    parts.append(np.array([max(0,abs(np.rad2deg(np.arctan(r[3])))-9.7)*.5]))
   return np.concatenate(parts)
  if wd==500:x=np.array([0.,offset/.01])
  else:
   initial=[2.5 if wd<=200 else (3. if wd>1e6 else .5),offset/.01]
   opt=least_squares(fun,initial,bounds=([-1.,-10.],[6.,10.]),diff_step=1e-4,max_nfev=100,xtol=1e-8,ftol=1e-8,gtol=1e-7)
   x=opt.x
  tt=trial(x);m,radii,chiefs,invalid=f.quantities(c,tt,fields=np.array([0.,9.08]),wd=wd,frequencies=(fc,fe))
  p=e.paraxial(c,tt,wd);oy=9.08/p['magnification'];chief=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);r,_,ok=e.ray(c,tt,e.N[1],wd,oy,0.,chief,p['ep'])
  record=dict(wd=wd,gap=float(tt[10]),backfocus=float(tt[-1]),frequency_center=fc,frequency_edge=fe,efl=p['efl'],ttl=float(sum(tt)),
     mtf=m[:,:2].tolist(),distortion_percent=float((chiefs[-1]/9.08-1)*100),cra_deg=float(abs(np.rad2deg(np.arctan(r[3])))),
     max_ray_radii=radii.tolist(),front_group_position=float(tt[10]+tt[-1]),rear_group_position=float(tt[-1]))
  output.append(record);print(json.dumps(record),flush=True)
  ROOT.joinpath('candidate_v17_focus_targets.json').write_text(json.dumps(output,indent=2),encoding='utf8')
if __name__=='__main__':run()
