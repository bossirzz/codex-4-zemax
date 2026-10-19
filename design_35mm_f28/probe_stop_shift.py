import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v18_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for a in s['aspheres']:
 j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
for delta in [-.2,-.1,0.,.1,.2,.3]:
 tt=t.copy();tt[4]+=delta;tt[5]-=delta
 m,rr,cc,iv=f.quantities(c,tt,fields=np.array([0.,9.08]),low_frequency=False)
 ti=tt.copy();ti[-1]+=5.;ti[10]=e.focus_group_gap(c,ti,1e7)+tt[10]-e.focus_group_gap(c,tt,500.)
 pi=e.paraxial(c,ti,1e7);oy=9.08/pi['magnification'];ch=e.aim_chief(c,ti,e.N[1],1e7,oy,pi['ep']);ray,hh,ok=e.ray(c,ti,e.N[1],1e7,oy,0.,ch,pi['ep'])
 print(json.dumps(dict(stop_shift=delta,mtf=m.tolist(),cra_inf=float(abs(np.rad2deg(np.arctan(ray[3])))),distortion=float((cc[-1]/9.08-1)*100))),flush=True)
