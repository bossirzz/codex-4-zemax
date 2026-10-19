import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v17_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
x=np.load(R/'real_stop_aiming_best_parameters.npy')
c=np.zeros(33);c[e.CURVED]=x[:14]*.05;c[15:]=x[29:]*np.tile([1e-5,1e-7,1e-9],6);t=x[14:29].copy()
caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
caps[:2]=9.4;caps[2:5]=7.1
m,r,chief,invalid=f.quantities(c,t)
print(json.dumps(dict(mtf=m.tolist(),radii=r.tolist(),margins=(caps-r).tolist(),clearances=f.blank_clearance(c,t,caps).tolist(),chief=chief.tolist(),invalid=invalid,paraxial=e.paraxial(c,t)),indent=2),flush=True)
for delta in [-2.,0.,2.,4.,5.5]:
 tt=t.copy();tt[-1]+=delta;tt[10]=e.focus_group_gap(c,tt,1e7)+(t[10]-e.focus_group_gap(c,t,500.))
 p=e.paraxial(c,tt,1e7);oy=9.08/p['magnification'];ch=e.aim_chief(c,tt,e.N[1],1e7,oy,p['ep']);ray,hh,ok=e.ray(c,tt,e.N[1],1e7,oy,0,ch,p['ep'])
 print('INF',delta,'gap',tt[10],'ttl',sum(tt),'CRA',abs(np.rad2deg(np.arctan(ray[3]))),'dist',(ray[1]/9.08-1)*100,flush=True)
s['curvatures']=c[:15].tolist();s['thicknesses']=t.tolist();s['ray_aiming']='Real'
s['aspheres']=[dict(surface=a,coefficients=c[15+3*j:18+3*j].tolist()) for j,a in enumerate(f.SURFS)]
for i in [0,1,2,3,4]:s['fixed_semi_diameters'][str(i+1)]=float(caps[i])
(R/'checkpoint_real_schema.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
