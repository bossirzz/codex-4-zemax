"""Export focus diagnostics; all curves are explicitly numerical predictions."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

r=Path(__file__).resolve().parent
a=json.loads((r/'far_focus_diagnostic_v29_scan_2000.json').read_text(encoding='utf8'))
x=np.array([z['defocus_mm'] for z in a])*1000
p=np.array([z['poly_mtf'] for z in a]);m=np.array([z['mono_d_mtf'] for z in a])
native=json.loads((r/'candidate_v29_native_acceptance.json').read_text(encoding='utf8'))
v=next(z for z in native if z['wd']==2000.)
# Middle sweep sample is the actual native baseline motor location.
origin=x[len(x)//2]
fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
for j,(ax,h) in enumerate(zip(axes.ravel(),[0.,3.,6.3,9.08])):
 for k,style in enumerate(['-','--']):
  if h==0 and k==1:continue
  ax.plot(x-origin,p[:,j,k],style,color='#1769aa',lw=2,label='F/d/C '+('T' if k==0 else 'S'))
  ax.plot(x-origin,m[:,j,k],style,color='#c36410',lw=1.5,label='d only '+('T' if k==0 else 'S'))
 ax.axvline(0,color='0.4',lw=.8);ax.axhline(.2,color='0.5',lw=.8,ls=':')
 ax.set_title(f'Image height {h:g} mm; {225 if h==0 else 180} lp/mm')
 ax.grid(alpha=.2);ax.set_ylim(0,.5);ax.set_xlim(-120,120)
 ax.set_xlabel('g11 motion from native 2 m configuration (micrometres)')
 ax.set_ylabel('MTF');ax.legend(fontsize=8,loc='upper right')
fig.suptitle('Frozen v29 optics: 2 m focus scan, nominal F/2.8, bare image plane\nNumerical predictions with physical lens clipping; native verification required',fontsize=14)
fig.tight_layout(rect=(0,0,1,.93))
fig.savefig(r/'far_focus_diagnostic_v29_focus_scan.png',dpi=220)
fig.savefig(r/'far_focus_diagnostic_v29_focus_scan.pdf')
print(r/'far_focus_diagnostic_v29_focus_scan.png')
