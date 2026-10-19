"""Half-micron diagnostic scan around the broad-search cold candidate."""
import json
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R,SURFACES

source=next(q for q in json.loads((R/'candidate_v29_thermal_snapshots.json').read_text())
            if q['wd']==500 and q['temperature_c']==-20)
target=json.loads((R/'candidate_v29_thermal_cold_broad_targets.json').read_text())[0]
c=np.r_[source['curvatures'],np.zeros(24)]
for a in source['aspheres']:
    j=15+3*SURFACES.index(a['surface'])
    c[j:j+3]=a['coefficients']
t=np.array(source['thicknesses'])
t[5],t[10],t[-1]=target['iris_gap'],target['gap'],target['backfocus']
caps=np.array(source['caps'])
e.N=np.array(source['indices'])
rows=[]
for delta in np.linspace(-.003,.003,13):
    tt=t.copy()
    tt[10]+=delta
    p=e.paraxial(c,tt,500)
    stop=abs(p['stop_scale'])*p['efl']/5.6
    mtf=f.quantities(c,tt,wd=500,frequencies=(225.,180.),fixed_stop_radius=stop,
        low_frequency=False,caps=caps)[0]
    rows.append(dict(delta_mm=float(delta),gap=float(tt[10]),predicted_mtf=mtf.tolist(),
        acceptance_margin=float(min(mtf[0].min()-.30,mtf[1:].min()-.20))))
(R/'candidate_v29_thermal_cold_broad_gap_scan.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(max(rows,key=lambda q:q['acceptance_margin'])))
