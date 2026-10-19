import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for a in s['aspheres']:
 j=15+3*f.SURFS.index(a['surface']);c[j:j+3]=a['coefficients']
cfg=next(x for x in json.loads((R/'candidate_v19_real_focus_targets.json').read_text(encoding='utf8')) if x['wd']==200)
t[10]=cfg['gap'];t[-1]=cfg['backfocus']
for degree in [10]:
 f.POWERS=np.array([(i,j) for total in range(degree+1) for i in range(total+1) for j in [total-i]])
 f.PUPIL=np.vstack([e.pupil_sample(14,40),[(np.cos(q),np.sin(q)) for q in np.arange(80)*np.pi/40]])
 for quad in [12,16,20,24]:
  f.TX,f.TW=np.polynomial.legendre.leggauss(quad);f.TX=(f.TX+1)/2;f.TW=f.TW/2
  m,r,h,iv=f.quantities(c,t,fields=np.array([0,9.08]),wd=200,frequencies=(150,120),fixed_stop_radius=s['fixed_stop_radius'],low_frequency=False)
  print('DEGREE',degree,'QUAD',quad,'MTF',m.tolist(),flush=True)
