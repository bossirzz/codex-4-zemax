import json
import numpy as np
import fast_otf_optimize as f
from joint_design_model import JointModel
m=JointModel();x=m.x
indices=[0,2,4,6,8,9,11,12,24,28,29,35,41,47,48,49,51,54]
def metric(xx):
 c,t,n,caps,stop=m.unpack(xx)
 return f.quantities(c,n,fields=np.array([9.08]),wd=200,frequencies=(150,120),fixed_stop_radius=stop,low_frequency=False)[0][0]
for i in indices:
 values=[]
 for h0 in [1e-4,1e-5,1e-6]:
  h=h0*(1+abs(x[i]));xp=x.copy();xm=x.copy();xp[i]+=h;xm[i]-=h
  values.append(((metric(xp)-metric(xm))/(2*h)).tolist())
 print(json.dumps(dict(parameter=i,gradient_by_step=values)),flush=True)
