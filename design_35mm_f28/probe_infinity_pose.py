import json
import numpy as np
from scipy.optimize import least_squares
import optical_engine as e
from joint_design_model import JointModel
m=JointModel();c,t,n,caps,stop=m.unpack(m.x)
def evaluate(x):
 tt=t.copy();tt[-1]+=x[0];tt[5]+=x[1];tt[10]=e.focus_group_gap(c,tt,1e7)+t[10]-e.focus_group_gap(c,t,500.)
 p=e.paraxial(c,tt,1e7);ch=e.aim_chief(c,tt,e.N[1],1e7,9.08/p['magnification'],p['ep']);r,_,ok=e.ray(c,tt,e.N[1],1e7,9.08/p['magnification'],0.,ch,p['ep'])
 cra=abs(np.rad2deg(np.arctan(r[3])))
 return tt,p,cra
for seed in [[5.,0.],[0.,0.],[2.,1.],[-2.,-1.]]:
 def residual(x):
  tt,p,cra=evaluate(x)
  return [p['efl']-35.,cra-9.8,max(0,.15-tt[10])*10,max(0,8.-tt[-1])*10]
 opt=least_squares(residual,seed,bounds=([-4.,-3.],[6.,3.]),max_nfev=100,diff_step=1e-5)
 tt,p,cra=evaluate(opt.x)
 print(json.dumps(dict(seed=seed,solution=opt.x.tolist(),efl=p['efl'],cra=cra,gap=tt[10],iris_gap=tt[5],backfocus=tt[-1],ttl=float(sum(tt)))),flush=True)
