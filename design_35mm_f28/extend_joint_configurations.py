"""Declare far configurations for an exported joint candidate, before native audit."""
import argparse,json
import numpy as np
import optical_engine as e
from joint_design_model import R,SURFACES
from optimize_best_conjugate import load_glasses
a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args()
s=json.loads((R/(q.tag+'_schema.json')).read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
c=np.r_[s['curvatures'],np.zeros(24)];t=np.array(s['thicknesses'])
for asp in s['aspheres']:
 j=15+3*SURFACES.index(asp['surface']);c[j:j+3]=asp['coefficients']
cfgs=json.loads((R/(q.tag+'_focus_targets.json')).read_text(encoding='utf8'))
offset=t[10]-e.focus_group_gap(c,t,500.)
for wd,back in [(2000.,0.),(1e7,5.)]:
 tt=t.copy();tt[-1]+=back;tt[10]=e.focus_group_gap(c,tt,wd)+offset
 cfgs.append(dict(wd=wd,gap=float(tt[10]),iris_gap=float(tt[5]),backfocus=float(tt[-1]),
  frequency_center=225.,frequency_edge=180.,fixed_stop_radius=s['fixed_stop_radius'],nominal_fno_target=2.8))
(R/(q.tag+'_validation_configs.json')).write_text(json.dumps(cfgs,indent=2),encoding='utf8')
print('DECLARED_CONFIGURATIONS',[x['wd'] for x in cfgs],flush=True)
