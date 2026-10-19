"""Check actual design MTF equivalence and timing against reference clipping."""
import json,time
import numpy as np
import fast_otf_optimize as f
from joint_design_model import JointModel
m=JointModel();c,t,n,caps,stop=m.unpack(m.x)
fast=f.clip_polygon
def evaluate():
 return [f.quantities(c,tt,wd=wd,fixed_stop_radius=stop,frequencies=ff,low_frequency=False)[0]
         for wd,tt,ff in [(500.,t,(225.,180.)),(200.,n,(150.,120.))]]
evaluate()
start=time.perf_counter();a=evaluate();new_seconds=time.perf_counter()-start
f.clip_polygon=f.clip_polygon_reference
start=time.perf_counter();b=evaluate();old_seconds=time.perf_counter()-start
f.clip_polygon=fast
difference=float(np.max(abs(np.array(a)-np.array(b))))
assert difference<1e-10
print(json.dumps(dict(maximum_mtf_difference=difference,compiled_seconds=new_seconds,reference_seconds=old_seconds)),flush=True)
