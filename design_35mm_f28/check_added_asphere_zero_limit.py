import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent;s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']]);c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for a in s['aspheres']:
 j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
extended=np.r_[c,np.zeros(6)];p=e.paraxial(c,t);radius=p['efl']/5.6;max_error=0.
for n in e.N:
 for field in [0.,3.,6.3,9.08]:
  a,ha,va=e.trace_bundle_stop_aiming(c,t,n,500,field/p['magnification'],p['ep'],f.PUPIL,radius,s['fixed_stop_radius'])
  b,hb,vb=e.trace_bundle_stop_aiming(extended,t,n,500,field/p['magnification'],p['ep'],f.PUPIL,radius,s['fixed_stop_radius'])
  assert np.array_equal(va,vb);max_error=max(max_error,float(np.max(abs(a-b))),float(np.max(abs(ha-hb))))
caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
assert np.array_equal(f.blank_clearance(c,t,caps),f.blank_clearance(extended,t,caps))
assert max_error<1e-12
print('MAXIMUM_RAY_DIFFERENCE_WITH_ZERO_ADDED_ASPHERES',max_error)
