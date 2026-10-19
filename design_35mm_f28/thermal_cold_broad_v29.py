"""Broad-start, max-min search for the v29 500 mm cold snapshot.

Pure numerical search. It writes separate cold targets and preserves the
targets being consumed by native verification. Native acceptance is required.
"""
import argparse
import json
import sys

import numpy as np
from scipy.optimize import minimize

import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R, SURFACES

sys.stdout.reconfigure(encoding='utf8')


def run(maxiter=40):
    tag='candidate_v29'
    case=next(q for q in json.loads((R/(tag+'_thermal_snapshots.json')).read_text(encoding='utf8'))
              if q['wd']==500 and q['temperature_c']==-20)
    schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
    nominal=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
    held=json.loads((R/(tag+'_thermal_focus_targets.json')).read_text(encoding='utf8'))
    c=np.r_[case['curvatures'],np.zeros(24)]
    for a in case['aspheres']:
        j=15+3*SURFACES.index(a['surface'])
        c[j:j+3]=a['coefficients']
    t=np.array(case['thicknesses'])
    caps=np.array(case['caps'])
    e.N=np.array(case['indices'])
    positions=[]
    for q in nominal:
        tt=np.array(schema['thicknesses'])
        tt[5],tt[10],tt[-1]=q['iris_gap'],q['gap'],q['backfocus']
        positions.append(np.array([sum(tt),sum(tt[6:]),sum(tt[11:])]))
    for q in held:
        if q['wd']==500 and q['temperature_c']==-20:
            continue
        positions.append(np.array(q['group_positions_from_sensor_mm']))
    geo_names=['clearance_'+str(i+1) for i in range(14)]
    geo_names+=['ttl','stroke_front','stroke_middle','stroke_rear','backfocus',
                'image_circle','illumination','valid_rays','max_diameter',
                'distortion_1','distortion_2','distortion_3']
    cache={}
    center_cache={}
    best=None

    def pose(x):
        tt=t.copy()
        tt[-1]+=x[0]
        tt[5]+=x[1]
        tt[10]=e.focus_group_gap(c,tt,500)+x[2]*.01
        return tt

    def precheck(tt):
        pos=np.array([sum(tt),sum(tt[6:]),sum(tt[11:])])
        stroke=np.ptp(np.vstack([*positions,pos]),axis=0)
        clearance=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,.02)
        geo=np.r_[clearance,.1*(60-sum(tt)),.1*(6.01-stroke),.1*(tt[-1]-8)]
        return geo,pos,stroke

    def evaluate(x,center_only=False):
        local=center_cache if center_only else cache
        key=tuple(x)
        if key in local:
            return local[key]
        tt=pose(x)
        geo,pos,stroke=precheck(tt)
        if geo.min() < -1e-5:
            result=dict(margin=-1.,geometry=geo,mtf=np.zeros((1 if center_only else 4,2)),
                        t=tt,positions=pos,stroke=stroke,invalid=True)
            local[key]=result
            return result
        p=e.paraxial(c,tt,500)
        stop=abs(p['stop_scale'])*p['efl']/5.6
        fields=np.array([0.]) if center_only else np.array([0.,3.,6.3,9.08])
        mtf,r,h,bad,diag=f.quantities(c,tt,fields=fields,wd=500,frequencies=(225.,180.),
            fixed_stop_radius=stop,low_frequency=False,caps=caps,return_diagnostics=True)
        if center_only:
            result=dict(margin=float(mtf.min()-.30),geometry=geo,mtf=mtf,
                        t=tt,positions=pos,stroke=stroke,invalid=False)
        else:
            distortion=100*(h[1:]/np.array([3.,6.3,9.08])-1)
            geo=np.r_[geo,.1*(h[-1]-9),diag['relative_illumination_primary'][-1]-.80,
                      -bad*.1,.1*(11.5-max(caps)),.3-abs(distortion)]
            margins=np.r_[mtf[0]-.30,(mtf[1:]-.20).ravel()]
            result=dict(margin=float(margins.min()),margins=margins,geometry=geo,
                mtf=mtf,t=tt,positions=pos,stroke=stroke,stop=stop,diagnostics=diag,
                distortion=distortion,invalid=False,efl=p['efl'])
        local[key]=result
        return result

    def record(x,result):
        tt=result['t']
        return dict(wd=500.,temperature_c=-20.,spacer_tce_ppm_per_k=case['spacer_tce_ppm_per_k'],
            mounting_assumption='Native thin-shell spacers; spacer material is assumed',
            manufactured_assembly_validated=False,
            iris_gap=float(tt[5]),gap=float(tt[10]),backfocus=float(tt[-1]),
            stop_radius=float(result['stop']),frequency_center=225.,frequency_edge=180.,
            predicted_mtf=result['mtf'].tolist(),predicted_efl=float(result['efl']),
            predicted_illumination=float(result['diagnostics']['relative_illumination_primary'][-1]),
            predicted_distortion_percent=result['distortion'].tolist(),
            predicted_image_circle_mm=float(2*(1+result['distortion'][-1]/100)*9.08),
            predicted_colour_centroid_span_um=result['diagnostics']['phase_centroid_colour_span_um'],
            min_geometry_constraint=float(result['geometry'].min()),
            geometry_constraints=dict(zip(geo_names,map(float,result['geometry']))),
            predicted_mtf_acceptance_margin=result['margin'],
            predicted_group_strokes_mm=result['stroke'].tolist(),
            group_positions_from_sensor_mm=result['positions'].tolist(),
            focus_corrections_mm=[float(tt[-1]-t[-1]),float(tt[5]-t[5]),float(tt[10]-t[10])],
            parameters=list(map(float,x)),evaluated_full_poses=len(cache),
            evaluated_center_poses=len(center_cache),status='PREDICTION_ONLY_REQUIRES_NATIVE_VALIDATION')

    def consider(x,label):
        nonlocal best
        result=evaluate(x)
        if result['invalid'] or result['geometry'].min() < -1e-7:
            return
        if best is None or result['margin']>best[1]['margin']+1e-7:
            best=(np.array(x),result)
            q=record(x,result)
            (R/(tag+'_thermal_cold_broad_targets.json')).write_text(
                json.dumps([q],ensure_ascii=False,indent=2),encoding='utf8')
            print('IMPROVED',label,json.dumps(q),flush=True)

    old=next(q for q in held if q['wd']==500 and q['temperature_c']==-20)
    old_t=t.copy()
    old_t[-1]=old['backfocus']
    old_t[5]=old['iris_gap']
    old_x=np.array([old_t[-1]-t[-1],old_t[5]-t[5],
                    (old['gap']-e.focus_group_gap(c,old_t,500))/.01])
    consider(old_x,'original')

    centers=[]
    for db in [-1.,-.6,-.3,0.,.3,.6,1.]:
        for dg in [-1.,-.6,-.3,0.,.3,.6,1.]:
            poses=[np.array([db,dg,u]) for u in np.linspace(0.,12.,17)]
            xx=max(poses,key=lambda z:evaluate(z,True)['margin'])
            centers.append((evaluate(xx,True)['margin'],xx))
        print('CENTER_GRID_ROW',db,'best_center_margin',max(q[0] for q in centers),flush=True)
    centers=sorted(centers,key=lambda q:q[0],reverse=True)
    starting=[old_x]
    for margin,xx in centers[:14]:
        if margin<-.04:
            continue
        scan=[xx+np.array([0.,0.,du]) for du in np.linspace(-.7,.7,15)]
        scan=[z for z in scan if not evaluate(z)['invalid'] and evaluate(z)['geometry'].min()>=-1e-7]
        if not scan:
            continue
        yy=max(scan,key=lambda z:evaluate(z)['margin'])
        starting.append(yy)
        consider(yy,'fine_focus_grid')
    starting.sort(key=lambda z:evaluate(z)['margin'],reverse=True)
    distinct=[]
    for x in starting:
        if all(np.linalg.norm(x[:2]-y[:2])>.12 for y in distinct):
            distinct.append(x)
    print('LOCAL_STARTS',json.dumps([z.tolist() for z in distinct[:5]]),flush=True)
    for count,ss in enumerate(distinct[:5]):
        start=np.r_[ss,evaluate(ss)['margin']]
        def constraints(z):
            q=evaluate(z[:3])
            if q['invalid']:
                return np.r_[np.full(len(geo_names),-1.),np.full(8,-1.)]
            return np.r_[q['geometry'],q['margins']-z[3]]
        a=minimize(lambda z:-z[3],start,method='SLSQP',
            bounds=[(-2.,2.),(-1.2,1.2),(-10.,20.),(-.15,.05)],
            constraints=[dict(type='ineq',fun=constraints)],
            options=dict(maxiter=maxiter,eps=3e-5,ftol=1e-10))
        consider(a.x[:3],f'SLSQP_{count}')
        print('LOCAL_DONE',count,'success',a.success,'message',a.message,
              'margin',evaluate(a.x[:3])['margin'],flush=True)
    if best is None:
        raise RuntimeError('No physically feasible cold candidate found')
    q=record(*best)
    q['search_complete']=True
    (R/(tag+'_thermal_cold_broad_targets.json')).write_text(
        json.dumps([q],ensure_ascii=False,indent=2),encoding='utf8')
    print('COLD_BROAD_COMPLETE',json.dumps(q),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--maxiter',type=int,default=40)
    args=parser.parse_args()
    run(args.maxiter)
