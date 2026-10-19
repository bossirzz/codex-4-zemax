"""Small local shape/focus refinement across warm conjugates and cold 500 mm.

Pure numerical work. Changed manufacturing geometry is mapped through native
Make Thermal scale factors, not substituted into an old thermal geometry.
The 23 ppm/K spacer assumption does not certify a manufactured assembly.
"""
import argparse
import copy
import json
import sys
import time

import numpy as np
from scipy.linalg import qr
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative

import optical_engine as e
import fast_otf_optimize as f
from joint_design_model import R, SURFACES, SCALES

sys.stdout.reconfigure(encoding='utf8')


class Feasible(Exception):
    pass


def run(steps=20, shape_count=6, seed_file=None, tag='candidate_v31_thermal_joint', near_floating=False):
    schema=json.loads((R/'candidate_v29_schema.json').read_text(encoding='utf8'))
    snaps=json.loads((R/'candidate_v29_thermal_snapshots.json').read_text(encoding='utf8'))
    warm=next(q for q in snaps if q['wd']==500 and q['temperature_c']==20)
    cold=next(q for q in snaps if q['wd']==500 and q['temperature_c']==-20)
    rows=json.loads((R/'candidate_v29_native_acceptance.json').read_text(encoding='utf8'))
    eight=json.loads((R/'candidate_v29_800_native.json').read_text(encoding='utf8'))[0]
    cold_target=json.loads((R/'candidate_v29_thermal_cold_broad_targets.json').read_text(encoding='utf8'))[0]
    thermal_targets=json.loads((R/'candidate_v29_thermal_focus_targets.json').read_text(encoding='utf8'))
    cw=np.r_[warm['curvatures'],np.zeros(24)]
    for a in warm['aspheres']:
        j=15+3*SURFACES.index(a['surface'])
        cw[j:j+3]=a['coefficients']
    tw=np.array(warm['thicknesses'])
    caps_w=np.array(warm['caps'])
    nw=np.array(warm['indices'])
    mappings={}
    for snap in snaps:
        ratio=np.ones(39)
        nz=cw[:15]!=0
        ratio[:15][nz]=np.array(snap['curvatures'])[nz]/cw[:15][nz]
        for a in snap['aspheres']:
            j=15+3*SURFACES.index(a['surface'])
            base=cw[j:j+3]
            fallback=ratio[a['surface']-1]**np.array([3,5,7])
            ratio[j:j+3]=np.divide(np.array(a['coefficients']),base,
                out=fallback.copy(),where=base!=0)
        mappings[(snap['wd'],snap['temperature_c'])]=dict(
            c=ratio,t=np.array(snap['thicknesses'])/tw,
            caps=np.array(snap['caps'])/caps_w,n=np.array(snap['indices']))
    map_cold=mappings[(500,-20)]
    original_cold=np.r_[cold['curvatures'],np.zeros(24)]
    for a in cold['aspheres']:
        j=15+3*SURFACES.index(a['surface'])
        original_cold[j:j+3]=a['coefficients']
    assert np.max(abs(cw*map_cold['c']-original_cold))<1e-15
    assert np.max(abs(tw*map_cold['t']-np.array(cold['thicknesses'])))<1e-12
    refs=[]
    for wd in [500,200,800]:
        row=eight if wd==800 else next(q for q in rows if q['wd']==wd)
        tt=tw.copy()
        tt[5],tt[10],tt[-1]=row['iris_gap'],row['gap'],row['backfocus']
        e.N=nw
        refs.append(dict(wd=wd,temp=20,pose=[tt[-1],tt[5]],
            offset=tt[10]-e.focus_group_gap(cw,tt,wd)))
    tc=tw*map_cold['t']
    tc[5],tc[10],tc[-1]=cold_target['iris_gap'],cold_target['gap'],cold_target['backfocus']
    e.N=map_cold['n']
    refs.append(dict(wd=500,temp=-20,pose=[tc[-1],tc[5]],
        offset=tc[10]-e.focus_group_gap(cw*map_cold['c'],tc,500)))
    inf=next(q for q in rows if q['wd']>=1e7)
    inf_t=tw.copy()
    inf_t[5],inf_t[10],inf_t[-1]=inf['iris_gap'],inf['gap'],inf['backfocus']
    e.N=nw
    pi=e.paraxial(cw,inf_t,1e7)
    si=abs(pi['stop_scale'])*pi['efl']/5.6
    inf_ri_base=f.quantities(cw,inf_t,fields=np.array([0.,9.08]),wd=1e7,
        frequencies=(1.,1.),fixed_stop_radius=si,low_frequency=False,
        caps=caps_w,return_diagnostics=True)[4]['relative_illumination_primary'][-1]
    seed=np.zeros(38) if seed_file is None else np.load(R/seed_file)
    cache={}
    calls=0
    started=time.time()
    best_score=np.inf
    best_x=seed.copy()
    best_data=None
    frozen=[q for q in thermal_targets if q['temperature_c']!=20 and
            not(q['wd']==500 and q['temperature_c']==-20)]
    far=next(q for q in rows if q['wd']==2000)

    def evaluate(x):
        nonlocal calls
        key=tuple(x)
        if key in cache:
            return cache[key]
        c=cw.copy()
        c[15:]+=x[:24]*SCALES
        t=tw.copy()
        t[7]+=x[24]*.01
        t[9]+=x[25]*.01
        constraints=[]
        metrics=[]
        poses=[]
        mtf_indices=[]
        try:
            for j,ref in enumerate(refs):
                mapping=mappings[(500,ref['temp'])]
                cc=c*mapping['c']
                tt=t*mapping['t']
                caps=caps_w*mapping['caps']
                e.N=mapping['n']
                dx=x[26+3*j:29+3*j]
                tt[-1]=ref['pose'][0]+dx[0]*.1
                tt[5]=ref['pose'][1]+dx[1]*.1
                tt[10]=e.focus_group_gap(cc,tt,ref['wd'])+ref['offset']+dx[2]*.001
                p=e.paraxial(cc,tt,ref['wd'])
                stop=abs(p['stop_scale'])*p['efl']/5.6
                frequencies=(150.,120.) if ref['wd']==200 else (225.,180.)
                mtf,r,h,bad,diag=f.quantities(cc,tt,wd=ref['wd'],frequencies=frequencies,
                    fixed_stop_radius=stop,low_frequency=False,caps=caps,return_diagnostics=True)
                dist=100*(h[1:]/np.array([3.,6.3,9.08])-1)
                limit=.595 if ref['wd']==200 else .295
                clearance=f.blank_clearance(cc,tt,caps)+np.where(e.MEDIA[:14]>=0,.05,0.)
                if ref['temp']!=20:
                    clearance+=np.where(e.MEDIA[:14]<0,.02,0.)
                    clearance-=.001
                else:
                    clearance-=.003
                geo=np.r_[20*clearance,.5*(59.8-sum(tt)),.5*(tt[-1]-8),
                    h[-1]-9,diag['relative_illumination_primary'][-1]-.80,
                    limit-abs(dist),-.1*bad]
                constraints.extend(geo.tolist())
                off=.209 if ref['temp']!=20 else (.210 if ref['wd']==800 else .206)
                start=len(constraints)
                constraints.extend(np.r_[10*(mtf[0]-.307),10*(mtf[1:]-off).ravel()].tolist())
                mtf_indices.extend(range(start,start+8))
                if ref['wd']==500 and ref['temp']==20:
                    constraints.append(5*(.045-abs(p['efl']-35)))
                pos=np.array([sum(tt),sum(tt[6:]),sum(tt[11:])])
                poses.append(pos)
                metrics.append(dict(wd=ref['wd'],temperature_c=ref['temp'],gap=float(tt[10]),
                    iris_gap=float(tt[5]),backfocus=float(tt[-1]),fixed_stop_radius=float(stop),
                    nominal_fno_target=2.8,frequency_center=frequencies[0],frequency_edge=frequencies[1],
                    predicted_mtf=mtf.tolist(),predicted_distortion_percent=dist.tolist(),
                    predicted_image_circle_mm=float(2*h[-1]),predicted_efl=p['efl'],ttl_mm=float(sum(tt)),
                    predicted_relative_illumination_primary=diag['relative_illumination_primary'][-1],
                    predicted_colour_centroid_span_um=diag['phase_centroid_colour_span_um'],
                    group_positions_from_sensor_mm=pos.tolist()))
            # Original infinity group gaps remain fixed; manufacturing thickness changes propagate.
            e.N=nw
            ti=t.copy()
            ti[5],ti[10],ti[-1]=inf['iris_gap'],inf['gap'],inf['backfocus']
            p=e.paraxial(c,ti,1e7)
            ch=e.aim_chief(c,ti,nw[1],1e7,9.08/p['magnification'],p['ep'])
            ray,_,ok=e.ray(c,ti,nw[1],1e7,9.08/p['magnification'],0,ch,p['ep'])
            cra=abs(np.rad2deg(np.arctan(ray[3])))
            stop=abs(p['stop_scale'])*p['efl']/5.6
            _,_,h,bad,diag=f.quantities(c,ti,fields=np.array([0.,9.08]),wd=1e7,
                frequencies=(1.,1.),fixed_stop_radius=stop,low_frequency=False,
                caps=caps_w,return_diagnostics=True)
            constraints.extend(np.r_[20*(f.blank_clearance(c,ti,caps_w)+np.where(e.MEDIA[:14]>=0,.05,0.)-.003),
                .5*(59.8-sum(ti)),.5*(ti[-1]-8),.5*(9.98-cra),h[-1]-9,
                diag['relative_illumination_primary'][-1]-(inf_ri_base-.015),-.1*bad].tolist())
            poses.append(np.array([sum(ti),sum(ti[6:]),sum(ti[11:])]))
            metrics.append(dict(wd=1e7,temperature_c=20,gap=inf['gap'],iris_gap=inf['iris_gap'],
                backfocus=inf['backfocus'],fixed_stop_radius=float(stop),nominal_fno_target=2.8,
                frequency_center=225.,frequency_edge=180.,predicted_cra_deg=cra,
                predicted_relative_illumination_primary=diag['relative_illumination_primary'][-1],
                illumination_constraint='Preserve the calibrated model baseline; native acceptance required',
                group_positions_from_sensor_mm=poses[-1].tolist()))
            for q in frozen:
                mp=mappings[(q['wd'],q['temperature_c'])]
                cc=c*mp['c']
                tt=t*mp['t']
                tt[5],tt[10],tt[-1]=q['iris_gap'],q['gap'],q['backfocus']
                cap=caps_w*mp['caps']
                constraints.extend((20*(f.blank_clearance(cc,tt,cap)+np.where(e.MEDIA[:14]>=0,.05,.02))).tolist())
                poses.append(np.array([sum(tt),sum(tt[6:]),sum(tt[11:])]))
            tf=t.copy()
            tf[5],tf[10],tf[-1]=far['iris_gap'],far['gap'],far['backfocus']
            poses.append(np.array([sum(tf),sum(tf[6:]),sum(tf[11:])]))
            strokes=np.ptp(np.array(poses),axis=0)
            constraints.extend((5*(6.01-strokes)).tolist())
            constraints.extend((1-c[:15]**2*caps_w**2-.01).tolist())
            constraints=np.array(constraints)
            data=dict(constraints=constraints,metrics=metrics,strokes=strokes,c=c,t=t,
                      mtf_indices=np.array(mtf_indices),calls=calls)
            calls+=1
            cache[key]=data
            return data
        finally:
            e.N=nw

    def export(x,data,feasible=False):
        s=copy.deepcopy(schema)
        c=data['c']
        t=data['t'].copy()
        best=data['metrics'][0]
        t[5],t[10],t[-1]=best['iris_gap'],best['gap'],best['backfocus']
        s.update(curvatures=c[:15].tolist(),thicknesses=t.tolist(),
            fixed_stop_radius=best['fixed_stop_radius'],ray_aiming='Real')
        s['aspheres']=[dict(surface=surf,coefficients=c[15+3*j:18+3*j].tolist())
                      for j,surf in enumerate(SURFACES)]
        s['thermal_design_assumption']='23 ppm/K thin-shell spacers; geometry mapped using native Make Thermal scales'
        s['thermal_design_requires_new_native_make_thermal']=True
        warm_configs=[{k:q[k] for k in ['wd','gap','iris_gap','backfocus','fixed_stop_radius',
            'nominal_fno_target','frequency_center','frequency_edge']} for q in data['metrics'] if q['temperature_c']==20]
        cold_ref=copy.deepcopy(next(q for q in data['metrics'] if q['temperature_c']==-20))
        cold_ref['spacer_tce_ppm_per_k']=23.
        cold_ref['stop_radius']=cold_ref['fixed_stop_radius']
        cold_ref['status']='PREDICTED_REFERENCE_REQUIRES_NEW_NATIVE_MAKE_THERMAL'
        prediction=dict(status='SURROGATE_FEASIBLE_REQUIRES_NATIVE_VALIDATION' if feasible else 'INCOMPLETE_CHECKPOINT',
            manufactured_assembly_validated=False,spacer_tce_ppm_per_k=23.,
            cases=data['metrics'],minimum_constraint=float(data['constraints'].min()),
            group_strokes_mm=data['strokes'].tolist(),selected_variables=selected_names,
            warm_native_800_calibration_error_max=baseline_800_error)
        (R/(tag+'_schema.json')).write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf8')
        (R/(tag+'_validation_configs.json')).write_text(json.dumps(warm_configs,indent=2),encoding='utf8')
        (R/(tag+'_focus_targets.json')).write_text(json.dumps(warm_configs,indent=2),encoding='utf8')
        (R/(tag+'_cold_reference_targets.json')).write_text(json.dumps([cold_ref],indent=2),encoding='utf8')
        (R/(tag+'_prediction.json')).write_text(json.dumps(prediction,ensure_ascii=False,indent=2),encoding='utf8')
        np.save(R/(tag+'_delta_parameters.npy'),x)

    initial=evaluate(seed)
    baseline_800_error=float(np.max(abs(np.array(initial['metrics'][2]['predicted_mtf']).ravel()-
                                      np.array(eight['native_mtf']))))
    print('BASELINE',json.dumps(dict(cold_mtf=initial['metrics'][3]['predicted_mtf'],
        minimum_constraint=float(initial['constraints'].min()),warm_800_max_native_error=baseline_800_error)),flush=True)
    # Screen 24 coefficients with the MTF Jacobian. Pivoting selects independent responses.
    scan_indices=np.r_[np.arange(24),[24,25,28,31,34,37]]
    def scanned(z):
        xx=seed.copy()
        xx[scan_indices]=z
        return evaluate(xx)['constraints']
    jac=approx_derivative(scanned,seed[scan_indices],method='3-point',
        abs_step=np.r_[np.full(24,2e-5),np.full(6,2e-4)])
    aj=jac[initial['mtf_indices'],:24].copy()
    aj[-8:]*=2.
    _,_,pivot=qr(aj,mode='economic',pivoting=True)
    selected=np.r_[pivot[:shape_count],[24,25,28,31,34,37],([29,30] if near_floating else [])].astype(int)
    names=[f'S{SURFACES[i//3]}_A{4+2*(i%3)}' for i in range(24)]
    names+=['L5_thickness','L6_thickness']
    names+=[f'{wd}_{temp}_{axis}' for wd,temp in [(500,20),(200,20),(800,20),(500,-20)]
            for axis in ['back','iris','g11']]
    selected_names=[names[i] for i in selected]
    print('SELECTED',json.dumps(selected_names),'seconds',time.time()-started,flush=True)
    (R/(tag+'_selected_jacobian.json')).write_text(json.dumps(dict(names=selected_names,
        pivots=pivot.tolist(),mtf_jacobian=aj.tolist()),indent=2),encoding='utf8')
    lo=np.full(38,-.025)
    hi=np.full(38,.025)
    lo[24:26]=-2.
    hi[24:26]=2.
    lo[26:]=-1.5
    hi[26:]=1.5
    lo[[28,31,34,37]]=-6.
    hi[[28,31,34,37]]=6.
    if near_floating:
        lo[[29,30]]=-.5
        hi[[29,30]]=.5
    def residual(z):
        nonlocal best_score,best_x,best_data
        xx=seed.copy()
        xx[selected]=z
        data=evaluate(xx)
        con=data['constraints']
        rr=np.r_[np.minimum(con,0.),.0001*(z-seed[selected])]
        score=float(rr@rr)
        if score<best_score:
            best_score=score
            best_x=xx.copy()
            best_data=data
            export(xx,data,float(con.min())>=-1e-4)
            print('CHECKPOINT',json.dumps(dict(calls=calls,seconds=time.time()-started,score=score,
                min_constraint=float(con.min()),cold_mtf=data['metrics'][3]['predicted_mtf'],
                warm_mtf_min=[float(np.array(q['predicted_mtf'])[1:].min()) for q in data['metrics'][:3]])),flush=True)
        if con.min()>=-1e-4:
            raise Feasible()
        return rr
    try:
        opt=least_squares(residual,seed[selected],bounds=(lo[selected],hi[selected]),
            jac=lambda z:approx_derivative(residual,z,method='3-point',
                abs_step=np.array([2e-5 if i<24 else 2e-4 for i in selected]),
                bounds=(lo[selected],hi[selected])),x_scale='jac',max_nfev=steps,
            ftol=1e-9,xtol=1e-10,gtol=1e-9,verbose=1)
        print('END',opt.message,opt.nfev,flush=True)
    except Feasible:
        print('FEASIBLE_SURROGATE_ONLY',flush=True)
    if best_data is not None:
        export(best_x,best_data,float(best_data['constraints'].min())>=-1e-4)
    print('COMPLETE',time.time()-started,'REQUIRES_NATIVE_VALIDATION',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--steps',type=int,default=20)
    parser.add_argument('--shape-count',type=int,default=6)
    parser.add_argument('--seed-file')
    parser.add_argument('--tag',default='candidate_v31_thermal_joint')
    parser.add_argument('--near-floating',action='store_true')
    args=parser.parse_args()
    run(args.steps,args.shape_count,args.seed_file,args.tag,args.near_floating)
