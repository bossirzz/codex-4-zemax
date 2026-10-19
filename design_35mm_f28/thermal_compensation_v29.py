"""Numerical thermal compensation with explicit distortion and stroke checks.

This reads native thermal snapshots; it never connects to OpticStudio. All
results are predictions until native_apply_thermal_focus.py independently
checks them. The 23 ppm/K spacer model is an assumed mechanical assembly.
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


def run(tag, warm_only=False, maxiter=80):
    cases=json.loads((R/(tag+'_thermal_snapshots.json')).read_text(encoding='utf8'))
    schema=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
    native=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
    nominal=np.array(schema['thicknesses'])
    positions=[]
    out=[]
    for row in native:
        tt=nominal.copy()
        tt[5],tt[10],tt[-1]=row['iris_gap'],row['gap'],row['backfocus']
        positions.append(np.array([sum(tt),sum(tt[6:]),sum(tt[11:])]))
    prior_path=R/'candidate_v27_nativecheck_thermal_focus_targets.json'
    prior=json.loads(prior_path.read_text(encoding='utf8')) if prior_path.exists() else []
    cases=sorted(cases,key=lambda q:(q['wd'],abs(q['temperature_c']-20),q['temperature_c']))
    for case in cases:
        wd=case['wd']
        temp=case['temperature_c']
        if warm_only and temp!=20:
            continue
        c=np.r_[case['curvatures'],np.zeros(24)]
        for a in case['aspheres']:
            j=15+3*SURFACES.index(a['surface'])
            c[j:j+3]=a['coefficients']
        t=np.array(case['thicknesses'])
        caps=np.array(case['caps'])
        e.N=np.array(case['indices'])
        offset=(t[10]-e.focus_group_gap(c,t,wd))/.01
        fc,fe=(225.,180.) if wd==500 else (150.,120.)
        cache={}
        geometry_names=['clearance_'+str(i+1) for i in range(14)]
        geometry_names+=['ttl','stroke_front','stroke_middle','stroke_rear',
                         'backfocus','image_circle','illumination','valid_rays','max_diameter']
        if wd<1e7:
            geometry_names+=['distortion_'+str(i) for i in range(1,4)]
        else:
            geometry_names+=['cra']

        def evaluate(x):
            key=tuple(x)
            if key in cache:
                return cache[key]
            tt=t.copy()
            tt[-1]+=x[0]
            tt[5]+=x[1]
            tt[10]=e.focus_group_gap(c,tt,wd)+x[2]*.01
            p=e.paraxial(c,tt,wd)
            stop=abs(p['stop_scale'])*p['efl']/5.6
            mtf,r,h,bad,diag=f.quantities(c,tt,wd=wd,frequencies=(fc,fe),
                fixed_stop_radius=stop,low_frequency=False,caps=caps,return_diagnostics=True)
            pos=np.array([sum(tt),sum(tt[6:]),sum(tt[11:])])
            stroke=np.ptp(np.vstack([*positions,pos]),axis=0)
            clearance=f.blank_clearance(c,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,.02)
            geo=np.r_[clearance,.1*(60.-sum(tt)),.1*(6.01-stroke),
                .1*(tt[-1]-8),.1*(h[-1]-9),
                diag['relative_illumination_primary'][-1]-.80,-bad*.1,
                .1*(11.5-max(caps))]
            oy=9.08/p['magnification']
            ch=e.aim_chief(c,tt,e.N[1],wd,oy,p['ep'])
            ray,_,ok=e.ray(c,tt,e.N[1],wd,oy,0,ch,p['ep'])
            cra=abs(np.rad2deg(np.arctan(ray[3])))
            distortion=100*(h[1:]/np.array([3.,6.3,9.08])-1)
            if wd>=1e7:
                geo=np.r_[geo,.1*(10.-cra)]
                con=geo
                objective=0.
                mtf_margin=None
            else:
                limit=.3 if wd==500 else .6
                geo=np.r_[geo,limit-abs(distortion)]
                mtfcon=np.r_[2*(mtf[0]-.308),2*(mtf[1:]-.210).ravel()]
                con=np.r_[geo,mtfcon]
                target=np.array([[.32,.32],[.215,.215],[.215,.215],[.215,.215]])
                objective=float(np.sum(np.maximum(target-mtf,0.)**2))
                mtf_margin=min(float(mtf[0].min()-.30),float(mtf[1:].min()-.20))
            score=objective+1000*np.dot(np.minimum(con,0),np.minimum(con,0))
            d=dict(score=score,objective=objective,constraints=con,geometry=geo,t=tt,
                mtf=mtf,stop=stop,cra=cra,positions=pos,stroke=stroke,diagnostics=diag,
                distortion=distortion,mtf_margin=mtf_margin)
            cache[key]=d
            return d

        bounds=[(-2.,2.),(-1.2,1.2),(-10.,20.)]
        def seed(db,dg,dd):
            tt=t.copy()
            tt[-1]+=db
            tt[5]+=dg
            return np.clip(np.array([db,dg,(t[10]+dd-e.focus_group_gap(c,tt,wd))/.01]),
                           [-2,-1.2,-10],[2,1.2,20])
        baseline=seed(0,0,0)
        if abs(temp-20.)<1e-6 or wd>=1e7:
            x=baseline
        else:
            sign=(temp-20)/40
            seeds=[baseline,seed(.10*sign,0,0),seed(.08*sign,.06*sign,0),
                   seed(.05*sign,-.05*sign,0)]
            old=next((q for q in prior if q['wd']==wd and q['temperature_c']==temp),None)
            if old is not None:
                seeds.append(seed(*old['focus_corrections_mm']))
            ranked=sorted(seeds,key=lambda z:evaluate(z)['score'])
            trials=[baseline]
            for ss in ranked[:2]:
                a=minimize(lambda z:evaluate(z)['score'],ss,method='Nelder-Mead',bounds=bounds,
                    options=dict(maxiter=maxiter,xatol=1e-4,fatol=1e-9))
                trials.append(a.x)
                b=minimize(lambda z:evaluate(z)['objective'],a.x,method='SLSQP',bounds=bounds,
                    constraints=[dict(type='ineq',fun=lambda z:evaluate(z)['constraints'])],
                    options=dict(maxiter=maxiter,eps=3e-5,ftol=1e-10))
                trials.append(b.x)
            feasible=[z for z in trials if evaluate(z)['constraints'].min()>=-1e-5]
            x=min(feasible or trials,key=lambda z:evaluate(z)['score'])
        d=evaluate(x)
        tt=d['t']
        positions.append(d['positions'])
        record=dict(wd=wd,temperature_c=temp,
            spacer_tce_ppm_per_k=case['spacer_tce_ppm_per_k'],
            mounting_assumption='Native thin-shell spacers; spacer material is assumed',
            manufactured_assembly_validated=False,
            iris_gap=float(tt[5]),gap=float(tt[10]),backfocus=float(tt[-1]),
            stop_radius=float(d['stop']),frequency_center=fc,frequency_edge=fe,
            predicted_mtf=d['mtf'].tolist(),predicted_cra=d['cra'],
            predicted_illumination=float(d['diagnostics']['relative_illumination_primary'][-1]),
            predicted_distortion_percent=d['distortion'].tolist(),
            predicted_image_circle_mm=float(2*(d['distortion'][-1]/100+1)*9.08),
            predicted_colour_centroid_span_um=d['diagnostics']['phase_centroid_colour_span_um'],
            min_constraint=float(d['constraints'].min()),
            geometry_constraints=dict(zip(geometry_names,map(float,d['geometry']))),
            predicted_mtf_acceptance_margin=d['mtf_margin'],
            predicted_group_strokes_mm=d['stroke'].tolist(),
            group_positions_from_sensor_mm=d['positions'].tolist(),
            focus_corrections_mm=[float(tt[-1]-t[-1]),float(tt[5]-t[5]),float(tt[10]-t[10])],
            evaluated_poses=len(cache),status='PREDICTION_ONLY_REQUIRES_NATIVE_VALIDATION')
        if temp==20 and wd<1e7:
            row=next(q for q in native if q['wd']==wd)
            record['warm_native_mtf_max_difference']=float(np.max(abs(
                d['mtf'].ravel()-np.array(row['mtf_sampling6']))))
            record['warm_native_distortion_max_difference_percent']=float(np.max(abs(
                d['distortion']-np.array(row['distortion_percent']))))
        out.append(record)
        print('THERMAL_FOCUS',json.dumps(record),flush=True)
        suffix='_thermal_warm_model_check.json' if warm_only else '_thermal_focus_targets.json'
        (R/(tag+suffix)).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('--warm-only',action='store_true')
    parser.add_argument('--maxiter',type=int,default=80)
    args=parser.parse_args()
    run(args.tag,args.warm_only,args.maxiter)
