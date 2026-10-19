"""Export a saved joint-design iterate; native validation determines acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses

R = Path(__file__).resolve().parent

def run(checkpoint, tag):
    source = R / checkpoint
    x = np.load(source)
    if x.shape not in [(51,), (57,)] or not np.isfinite(x).all():
        raise ValueError('Invalid joint-design checkpoint')
    s = json.loads((R / 'candidate_v19_schema.json').read_text(encoding='utf8'))
    load_glasses([(g['catalog'], g['name']) for g in s['glasses']])
    c = np.zeros(39 if len(x) == 57 else 33)
    c[e.CURVED] = x[:14] * .05
    c[15:33] = x[29:47] * np.tile([1e-5, 1e-7, 1e-9], 6)
    if len(x) == 57:
        c[33:] = x[51:57] * np.tile([1e-5, 1e-7, 1e-9], 2)
    t = x[14:29].copy()
    near = t.copy()
    near[-1] += x[47]
    near[5] += x[48]
    near[10] = e.focus_group_gap(c, near, 200.) + x[49] * .01
    p = e.paraxial(c, t)
    ray, _, ok = e.ray(c, t, e.N[1], 500., 0., 0., p['efl'] / 5.6, p['ep'])
    if not ok:
        raise ValueError('Cannot determine physical iris radius')
    stop = float(abs(ray[9]))
    s.update(curvatures=c[:15].tolist(), thicknesses=t.tolist(),
             fixed_stop_radius=stop, ray_aiming='Real',
             focus_architecture='Three floating groups; test aperture explicitly F/2.8 in each configuration')
    surfaces = f.SURFS + ([2, 9] if len(x) == 57 else [])
    s['aspheres'] = [dict(surface=a, coefficients=c[15+3*j:18+3*j].tolist())
                     for j, a in enumerate(surfaces)]
    for i in [9, 10]:
        s['fixed_semi_diameters'][str(i+1)] = float(x[50])
    cfgs = []
    focus_offset = t[10] - e.focus_group_gap(c, t, 500.)
    for wd, tt in [(500., t), (200., near)]:
        cfgs.append(dict(wd=wd, gap=float(tt[10]), iris_gap=float(tt[5]),
                         backfocus=float(tt[-1]), frequency_center=225. if wd == 500 else 150.,
                         frequency_edge=180. if wd == 500 else 120.,
                         fixed_stop_radius=stop, nominal_fno_target=2.8))
    for wd, back_offset in [(2000., 0.), (1e7, 5.)]:
        tt = t.copy()
        tt[-1] += back_offset
        tt[10] = e.focus_group_gap(c, tt, wd) + focus_offset
        cfgs.append(dict(wd=wd, gap=float(tt[10]), iris_gap=float(tt[5]),
                         backfocus=float(tt[-1]), frequency_center=225., frequency_edge=180.,
                         fixed_stop_radius=stop, nominal_fno_target=2.8))
    (R / (tag+'_schema.json')).write_text(json.dumps(s, ensure_ascii=False, indent=2), encoding='utf8')
    (R / (tag+'_scheduled_focus_targets.json')).write_text(json.dumps(cfgs, indent=2), encoding='utf8')
    info = dict(checkpoint=str(source), parameter_count=len(x),
                sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                acceptance='NOT YET VALIDATED; saved iterate is not a qualified design')
    (R / (tag+'_checkpoint_provenance.json')).write_text(json.dumps(info, indent=2), encoding='utf8')
    print(json.dumps(info), flush=True)

if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('checkpoint')
    a.add_argument('tag')
    q = a.parse_args()
    run(q.checkpoint, q.tag)
