import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v21_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(24)];t=np.array(s['thicknesses']);surfs=f.SURFS+[2,9]
for a in s['aspheres']:
 j=15+3*surfs.index(a['surface']);c[j:j+3]=a['coefficients']
caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
cfgs=json.loads((R/'candidate_v21_native_acceptance.json').read_text(encoding='utf8'))
for cfg in cfgs:
 tt=t.copy();tt[5]=cfg['iris_gap'];tt[10]=cfg['gap'];tt[-1]=cfg['backfocus'];wd=cfg['wd']
 clear=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
 vals=f.quantities(c,tt,fields=np.array([0.,3.,6.3,9.08]),wd=wd,frequencies=tuple(cfg['frequencies']),fixed_stop_radius=cfg['stop_radius'],low_frequency=False,return_complex=True)
 m,r,h,invalid,otf=vals
 out=dict(wd=wd,efl=e.paraxial(c,tt,wd)['efl'],mtf=m.tolist(),native_mtf=cfg['mtf_sampling6'],
  monochromatic_mtf=np.abs(otf).tolist(),edge_otf_real=otf[-1].real.tolist(),edge_otf_imag=otf[-1].imag.tolist(),
  geometric_clearances=clear.tolist(),ray_caps=caps.tolist(),max_ray_radii=r.tolist(),invalid=int(invalid))
 print(json.dumps(out),flush=True)
print('THICKNESSES',t.tolist(),'CURVATURES',c[:15].tolist(),flush=True)
