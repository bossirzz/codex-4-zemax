"""Search focus positions for one fixed optical prescription and physical iris."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def run(tag,gap_index=10,selected_wds=None,fno=2.8):
 s=json.loads((ROOT/(tag+'_schema.json')).read_text(encoding='utf8'))
 load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
 c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
 for a in s['aspheres']:
  j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
 caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 stop=s['fixed_stop_radius']*2.8/fno;offset=t[gap_index]-e.focus_group_gap(c,t,500.,gap_index)
 output=[];front=t[gap_index]+t[-1]
 for wd in (selected_wds or [500.,450.,400.,300.,800.,200.,2000.,1e7]):
  fc,fe=(150.,120.) if wd<=200 else (225.,180.)
  dist_limit=.003 if 400<=wd<=500 else .006
  cache={}
  def evaluate(x):
   key=tuple(x)
   if key in cache:return cache[key]
   tt=t.copy();tt[-1]+=x[0];tt[gap_index]=e.focus_group_gap(c,tt,wd,gap_index)+x[1]*.01
   try:
    m,r,chief,invalid=f.quantities(c,tt,fields=np.array([0.,9.08]),wd=wd,frequencies=(fc,fe),fixed_stop_radius=stop,low_frequency=False)
    p=e.paraxial(c,tt,wd);oy=9.08/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep']);ray,hh,ok=e.ray(c,tt,e.N[1],wd,oy,0,ch,p['ep'])
    cra=abs(np.rad2deg(np.arctan(ray[3])))
    clear=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
    con=np.r_[clear,.1*(59.8-sum(tt)),.1*(8-abs(tt[gap_index]+tt[-1]-front)),
       .1*(tt[-1]-8),dist_limit*9.08-abs(chief[-1]-9.08),-invalid*.01]
    if wd>=1e7:con=np.r_[con,.1*(9.97-cra)]
    desired=np.array([[.40,.40],[.30,.30]]) if wd<=200 else np.array([[.36,.36],[.26,.26]])
    score=float(np.sum((m-desired)**2))+500*np.dot(np.minimum(con,0),np.minimum(con,0))
    data=dict(score=score,constraints=con,tt=tt,mtf=m,p=p,radii=r,chief=chief,cra=cra)
   except (ValueError,ZeroDivisionError):data=dict(score=1e6)
   cache[key]=data;return data
  if wd==500:x=np.array([0.,offset/.01])
  else:
   seeds=[np.array([b,g]) for b in np.linspace(-1.,5.,7) for g in [-2.,0.,2.,4.,6.,8.,10.,12.]]
   ranked=sorted(seeds,key=lambda z:evaluate(z)['score'])
   attempts=[]
   for seed in ranked[:3]:
    opt=minimize(lambda z:evaluate(z)['score'],seed,method='Nelder-Mead',bounds=[(-1.,5.),(-5.,16.)],options=dict(maxiter=65,xatol=2e-3,fatol=1e-7))
    attempts.append(opt.x)
   x=min(attempts,key=lambda z:evaluate(z)['score'])
   # Distortion and material clearances are acceptance constraints, rather
   # than a cost that can be traded away for a little extra centre MTF.
   feasible=[z for z in ranked+attempts if evaluate(z).get('constraints',np.array([-1.])).min()>=-1e-6]
   for seed in attempts:
    opt=minimize(lambda z:evaluate(z)['score'],seed,method='SLSQP',bounds=[(-1.,5.),(-5.,16.)],
     constraints=[dict(type='ineq',fun=lambda z:evaluate(z)['constraints'])],
     options=dict(maxiter=30,ftol=1e-8,eps=1e-3))
    if evaluate(opt.x)['constraints'].min()>=-1e-6:feasible.append(opt.x)
   if feasible:x=min(feasible,key=lambda z:evaluate(z)['score'])
  d=evaluate(x);tt=d['tt']
  rr=dict(wd=wd,gap=float(tt[gap_index]),gap_surface=gap_index+1,backfocus=float(tt[-1]),frequency_center=fc,frequency_edge=fe,
   mtf=d['mtf'].tolist(),distortion_percent=float((d['chief'][-1]/9.08-1)*100),cra_deg=float(d['cra']),
   efl=d['p']['efl'],ttl=float(sum(tt)),min_constraint=float(d['constraints'].min()),
   front_group_position=float(tt[gap_index]+tt[-1]),rear_group_position=float(tt[-1]),max_ray_radii=d['radii'].tolist(),
   fixed_stop_radius=stop,nominal_fno=fno,ray_aiming='Real')
  output.append(rr);print(json.dumps(rr),flush=True)
  suffix='' if gap_index==10 else f'_group{gap_index+1}'
  if fno!=2.8:suffix+=f'_F{fno:g}'.replace('.','')
  (ROOT/(tag+suffix+'_real_focus_targets.json')).write_text(json.dumps(output,indent=2),encoding='utf8')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--gap-index',type=int,default=10);a.add_argument('--wds',type=float,nargs='+');a.add_argument('--fno',type=float,default=2.8);q=a.parse_args();run(q.tag,q.gap_index,q.wds,q.fno)
