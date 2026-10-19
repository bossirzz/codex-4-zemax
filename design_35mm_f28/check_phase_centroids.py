"""Compare analytic PSF moments with the independent low-frequency OTF phase."""
import json,numpy as np
import fast_otf_optimize as f
from joint_design_model import JointModel,R
m=JointModel();x=np.load(R/'candidate_v27_parameters.npy');c,t,n,caps,stop=m.unpack(x)
for wd,tt in [(500.,t),(200.,n)]:
 values=f.quantities(c,tt,wd=wd,fixed_stop_radius=stop,low_frequency=False,return_diagnostics=True)
 analytic=np.array(values[4]['phase_centroids_y_mm'])
 otf=f.quantities(c,tt,wd=wd,frequencies=(.02,.02),fixed_stop_radius=stop,low_frequency=False,return_complex=True)[4]
 finite=-np.angle(otf[:,:,0])/(2*np.pi*.02)
 error=float(np.max(abs(analytic-finite))*1000)
 assert error<.01,(wd,error)
 print(json.dumps(dict(wd=wd,maximum_centroid_difference_um=error,colour_span_um=values[4]['phase_centroid_colour_span_um'])),flush=True)
