import argparse,json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args()
s=json.loads((R/(q.tag+'_schema.json')).read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for asp in s['aspheres']:
 j=15+3*f.SURFS.index(asp['surface']);c[j:j+3]=asp['coefficients']
p=e.paraxial(c,t);r,hh,ok=e.ray(c,t,e.N[1],500.,0.,0.,p['efl']/5.6,p['ep'])
assert ok
s['fixed_stop_radius']=abs(float(r[9]));s['ray_aiming']='Real'
(R/(q.tag+'_schema.json')).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(dict(stop_radius=s['fixed_stop_radius'],paraxial=p),indent=2))
