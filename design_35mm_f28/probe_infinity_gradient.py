import json
from pathlib import Path
import numpy as np
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v18_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
x=np.load(R/'infinity_constrained_best_parameters.npy')
caps=np.array([s['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
def unpack(x):
 c=np.zeros(33);c[e.CURVED]=x[:14]*.05;c[15:]=x[29:]*np.tile([1e-5,1e-7,1e-9],6);return c,x[14:29].copy()
def inspect(x):
 c,t=unpack(x);tt=t.copy();tt[-1]+=5.;tt[10]=e.focus_group_gap(c,tt,1e7)+t[10]-e.focus_group_gap(c,t,500.)
 pi=e.paraxial(c,tt,1e7);oy=9.08/pi['magnification'];ch=e.aim_chief(c,tt,e.N[1],1e7,oy,pi['ep'])
 r,hh,ok=e.ray(c,tt,e.N[1],1e7,oy,0.,ch,pi['ep'])
 return abs(np.rad2deg(np.arctan(r[3]))),tt,f.blank_clearance(c,tt,caps)
angle,tt,clear=inspect(x)
grad=approx_derivative(lambda z:inspect(z)[0],x,method='3-point',abs_step=1e-4*(1+abs(x))).ravel()
ix=np.argsort(abs(grad))[-8:]
print(json.dumps(dict(cra=angle,gap=float(tt[10]),clear=clear.tolist(),gradient_indices=ix.tolist(),gradient=grad[ix].tolist(),norm=float(np.linalg.norm(grad))),indent=2),flush=True)
