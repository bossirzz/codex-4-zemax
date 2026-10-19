"""Common packing and declared real-glass model for coupled-conjugate refinement."""
import copy,json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses

R=Path(__file__).resolve().parent
SURFACES=f.SURFS+[2,9]
SCALES=np.tile([1e-5,1e-7,1e-9],8)

class JointModel:
 def __init__(self,tag='candidate_v21'):
  self.schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
  load_glasses([(g['catalog'],g['name']) for g in self.schema['glasses']])
  c=np.r_[self.schema['curvatures'],np.zeros(24)];t=np.array(self.schema['thicknesses'])
  for a in self.schema['aspheres']:
   j=15+3*SURFACES.index(a['surface']);c[j:j+3]=a['coefficients']
  rows=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
  best=next(q for q in rows if q['wd']==500);near=next(q for q in rows if q['wd']==200)
  t[5],t[10],t[-1]=best['iris_gap'],best['gap'],best['backfocus']
  nt=t.copy();nt[5],nt[10],nt[-1]=near['iris_gap'],near['gap'],near['backfocus']
  self.x=np.r_[c[e.CURVED]/.05,t,c[15:33]/SCALES[:18],nt[-1]-t[-1],nt[5]-t[5],
               (nt[10]-e.focus_group_gap(c,nt,200))/.01,self.schema['fixed_semi_diameters']['10'],c[33:]/SCALES[18:]]
  self.caps=np.array([self.schema['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
 def substitute(self,changes):
  c,t,near,_,_=self.unpack(self.x)
  offset=t[10]-e.focus_group_gap(c,t,500.)
  pairs=[(g['catalog'],g['name']) for g in self.schema['glasses']]
  for lens,name in changes:
   if not 1<=lens<=8:raise ValueError('Lens index must be 1..8')
   pairs[lens-1]=('CDGM',name)
  load_glasses(pairs)
  self.schema['glasses']=copy.deepcopy(e.GLASSES)
  self.x[24]=e.focus_group_gap(c,t,500.)+offset
 def unpack(self,x):
  c=np.zeros(39);c[e.CURVED]=x[:14]*.05;c[15:33]=x[29:47]*SCALES[:18];c[33:]=x[51:57]*SCALES[18:]
  t=x[14:29].copy();caps=self.caps.copy();caps[9:11]=x[50]
  near=t.copy();near[-1]+=x[47];near[5]+=x[48];near[10]=e.focus_group_gap(c,near,200)+x[49]*.01
  p=e.paraxial(c,t);r,_,ok=e.ray(c,t,e.N[1],500.,0.,0.,p['efl']/5.6,p['ep'])
  if not ok:raise ValueError('Physical stop ray invalid')
  return c,t,near,caps,abs(r[9])
 def export(self,x,tag):
  c,t,near,caps,stop=self.unpack(x);s=copy.deepcopy(self.schema)
  s.update(curvatures=c[:15].tolist(),thicknesses=t.tolist(),fixed_stop_radius=float(stop),ray_aiming='Real')
  s['aspheres']=[dict(surface=a,coefficients=c[15+3*j:18+3*j].tolist()) for j,a in enumerate(SURFACES)]
  s['fixed_semi_diameters']={str(i+1):float(caps[i]) for i in range(15) if i!=5}
  s['focus_architecture']='Three floating groups; nominal F/2.8 explicitly checked by native ISFN at each conjugate'
  (R/(tag+'_schema.json')).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
  cfgs=[dict(wd=wd,gap=float(tt[10]),iris_gap=float(tt[5]),backfocus=float(tt[-1]),
    frequency_center=fc,frequency_edge=fe,fixed_stop_radius=float(stop),nominal_fno_target=2.8)
    for wd,tt,fc,fe in [(500.,t,225.,180.),(200.,near,150.,120.)]]
  (R/(tag+'_focus_targets.json')).write_text(json.dumps(cfgs,indent=2),encoding='utf8')
  np.save(R/(tag+'_parameters.npy'),x)

if __name__=='__main__':
 m=JointModel();np.save(R/'candidate_v21_native_seed.npy',m.x)
 print('PARAMETER_COUNT',len(m.x),'NEAR_FOCUS_OFFSETS',m.x[47:50].tolist(),flush=True)
