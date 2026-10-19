"""Move the physical iris slightly while retaining positive full-blank space."""
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
t[4]-=.09;t[5]+=.09
for i in [3,4,5]:s['fixed_semi_diameters'][str(i)]=7.
caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
m,r,chief,invalid=f.quantities(c,t)
clear=f.blank_clearance(c,t,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
assert invalid==0 and min(caps-r)>=.025 and min(clear)>=-1e-6,(caps-r,clear)
p=e.paraxial(c,t);ray,hh,ok=e.ray(c,t,e.N[1],500,0,0,p['efl']/5.6,p['ep'])
s['fixed_stop_radius']=abs(float(ray[9]));s['thicknesses']=t.tolist();s['ray_aiming']='Real'
s['revision_reason']='Physical STOP moved 0.09 mm toward object; L2/L3 shared blank radius 7.0 mm clears all sampled rays and preserves iris clearance.'
(R/'candidate_v19_schema.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
report=dict(mtf=m.tolist(),clearance_margins=clear.tolist(),ray_clearance=(caps-r).tolist(),chief_image_heights=chief.tolist(),paraxial=p,stop_radius=s['fixed_stop_radius'])
(R/'candidate_v19_prebuild_geometry.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2),flush=True)
