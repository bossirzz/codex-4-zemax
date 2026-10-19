import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
import verify_design as v
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']]);v.N=e.N
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for a in s['aspheres']:
 j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
cfg=next(x for x in json.loads((R/'candidate_v19_real_focus_targets.json').read_text(encoding='utf8')) if x['wd']==200)
t[10]=cfg['gap'];t[-1]=cfg['backfocus']
for g in [192,256]:
 caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 audit=v.diffraction(c,t,wd=200,fields=np.array([0,9.08]),grid=g,fixed_stop_radius=s['fixed_stop_radius'],caps=caps)
 (R/f'candidate_v19_near_scalar_{g}.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
 print('GRID',g,[(d['tangential'][150],d['tangential'][120],d['sagittal'][120]) for d in audit['fields']],flush=True)
