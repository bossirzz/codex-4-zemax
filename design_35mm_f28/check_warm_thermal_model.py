"""Ensure the thermal snapshot reproduces the independently audited warm model."""
import json,numpy as np
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R,SURFACES
rows=json.loads((R/'candidate_v27_nativecheck_native_acceptance.json').read_text())
cases=json.loads((R/'candidate_v27_nativecheck_thermal_snapshots.json').read_text())
for q in cases:
 if q['temperature_c']!=20 or q['wd']>=1e7:continue
 c=np.r_[q['curvatures'],np.zeros(24)];t=np.array(q['thicknesses']);e.N=np.array(q['indices'])
 for a in q['aspheres']:
  j=15+3*SURFACES.index(a['surface']);c[j:j+3]=a['coefficients']
 p=e.paraxial(c,t,q['wd']);stop=abs(p['stop_scale'])*p['efl']/5.6
 native=next(x for x in rows if x['wd']==q['wd'])
 vals=f.quantities(c,t,wd=q['wd'],frequencies=tuple(native['frequencies']),fixed_stop_radius=stop,
  low_frequency=False,caps=np.array(q['caps']),return_diagnostics=True)
 unclipped=f.quantities(c,t,wd=q['wd'],frequencies=tuple(native['frequencies']),fixed_stop_radius=stop,low_frequency=False)[0]
 assert np.max(abs(vals[0]-unclipped))<1e-10,'Unexpected clipping of the audited clear pupil'
 error=float(np.max(abs(vals[0].ravel()-np.array(native['mtf_sampling6']))))
 # The current corner has a known ~0.005 surrogate/native difference;
 # exact native operands remain the acceptance authority.
 assert error<.006,(q['wd'],error)
 print(json.dumps(dict(wd=q['wd'],maximum_native_mtf_difference=error,clipped_counts=vals[4]['clipped_ray_counts'])),flush=True)
