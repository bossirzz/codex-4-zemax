"""Screen catalogue substitutions while preserving d-line paraxial powers."""
import json,sys
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
pairs=[(g['catalog'],g['name']) for g in s['glasses']]
load_glasses(pairs);original_nd=[g['nd'] for g in e.GLASSES]
x=np.load(R/'two_conjugates_high_quad_best_parameters.npy')
c=np.zeros(33);c[e.CURVED]=x[:14]*.05;c[15:]=x[29:47]*np.tile([1e-5,1e-7,1e-9],6);t=x[14:29].copy()
records=[];best=None
for swap3 in [False,True]:
 for swap7 in [False,True]:
  for swap8 in [False,True]:
   pp=pairs.copy()
   if swap3:pp[2]=('CDGM','H-ZLAF55D')
   if swap7:pp[6]=('CDGM','H-ZLAF55D')
   if swap8:pp[7]=('CDGM','H-QF8')
   load_glasses(pp);nd=[g['nd'] for g in e.GLASSES];cc=c.copy();tt=t.copy()
   if swap3:
    cc[3]*=(original_nd[2]-original_nd[1])/(nd[2]-nd[1]);cc[4]*=(original_nd[2]-1)/(nd[2]-1);tt[3]*=nd[2]/original_nd[2]
   if swap7:
    ratio=(original_nd[6]-1)/(nd[6]-1);cc[[11,12]]*=ratio;cc[24:27]*=ratio;tt[11]*=nd[6]/original_nd[6]
   if swap8:
    ratio=(original_nd[7]-1)/(nd[7]-1);cc[[13,14]]*=ratio;cc[18:21]*=ratio;tt[13]*=nd[7]/original_nd[7]
   p=e.paraxial(cc,tt);ray,hh,ok=e.ray(cc,tt,e.N[1],500,0,0,p['efl']/5.6,p['ep']);stop=abs(ray[9])
   tn=tt.copy();tn[-1]+=x[47];tn[5]+=x[48];tn[10]=e.focus_group_gap(cc,tn,200)+x[49]*.01
   mb,rb,hb,iv=f.quantities(cc,tt,fields=np.array([0,6.3,9.08]),fixed_stop_radius=stop,low_frequency=False)
   mn,rn,hn,iv=f.quantities(cc,tn,fields=np.array([0,6.3,9.08]),wd=200,frequencies=(150,120),fixed_stop_radius=stop,low_frequency=False)
   score=float(np.sum(np.maximum(0,np.array([[.305,.305],[.15,.15],[.210,.210]])-mb)**2)+np.sum(np.maximum(0,np.array([[.305,.305],[.15,.15],[.215,.215]])-mn)**2))
   rr=dict(glasses=pp,mtf500=mb.tolist(),mtf200=mn.tolist(),efl=p['efl'],score=score)
   records.append(rr);print(json.dumps(rr),flush=True)
   if best is None or score<best[0]:best=(score,cc.copy(),tt.copy(),stop,[dict(g) for g in e.GLASSES],rr)
score,cc,tt,stop,glasses,result=best
variant=dict(s);variant['glasses']=glasses;variant['curvatures']=cc[:15].tolist();variant['thicknesses']=tt.tolist();variant['fixed_stop_radius']=float(stop)
variant['fixed_semi_diameters']=dict(s['fixed_semi_diameters']);variant['fixed_semi_diameters']['10']=float(x[50]);variant['fixed_semi_diameters']['11']=float(x[50])
variant['aspheres']=[dict(surface=a,coefficients=cc[15+3*j:18+3*j].tolist()) for j,a in enumerate(f.SURFS)]
variant['focus_architecture']='Three floating groups, fixed sensor'
(R/'glass_variant_screen.json').write_text(json.dumps(records,indent=2),encoding='utf8')
(R/'glass_screen_seed_schema.json').write_text(json.dumps(variant,ensure_ascii=False,indent=2),encoding='utf8')
newx=x.copy();newx[:14]=cc[e.CURVED]/.05;newx[14:29]=tt;newx[29:47]=cc[15:]/np.tile([1e-5,1e-7,1e-9],6)
np.save(R/'glass_screen_seed_parameters.npy',newx)
print('BEST_VARIANT',json.dumps(result),flush=True)
