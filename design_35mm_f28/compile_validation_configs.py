import argparse,json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args()
s=json.loads((R/(q.tag+'_schema.json')).read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(18)];t=np.array(s['thicknesses'])
for asp in s['aspheres']:
 j=15+3*f.SURFS.index(asp['surface']);c[j:j+3]=asp['coefficients']
cfgs=json.loads((R/(q.tag+'_real_focus_targets.json')).read_text(encoding='utf8'))
cfgs=[dict(wd=500.,gap=float(t[10]),backfocus=float(t[-1]),frequency_center=225.,frequency_edge=180.,fixed_stop_radius=s['fixed_stop_radius'])]+[x for x in cfgs if x['wd']!=500.]
for cfg in cfgs:
 if cfg['wd']>=1e7:
  tt=t.copy();tt[-1]+=5.;tt[10]=e.focus_group_gap(c,tt,1e7)+t[10]-e.focus_group_gap(c,t,500.)
  cfg.update(gap=float(tt[10]),backfocus=float(tt[-1]))
 cfg.setdefault('fixed_stop_radius',s['fixed_stop_radius'])
(R/(q.tag+'_validation_configs.json')).write_text(json.dumps(cfgs,indent=2),encoding='utf8')
print([(x['wd'],x['gap'],x['backfocus']) for x in cfgs])
