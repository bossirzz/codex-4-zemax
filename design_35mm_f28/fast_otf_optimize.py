"""Optimize a small set of MTF values with an angular-pupil polynomial surrogate.

The surrogate is checked against scalar diffraction and native OpticStudio;
it never provides the final acceptance result.
"""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative
from scipy.spatial import ConvexHull
from numba import njit
import optical_engine as e
from optimize_best_conjugate import load_glasses
ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
SURFS=[1,15,10,12,3,7]
POWERS=np.array([(i,j) for total in range(11) for i in range(total+1) for j in [total-i]])
PUPIL=np.vstack([e.pupil_sample(8,28),[(np.cos(q),np.sin(q)) for q in np.arange(80)*np.pi/40]])
GX,GW=np.polynomial.legendre.leggauss(20)
TX,TW=np.polynomial.legendre.leggauss(24)
TX=(TX+1)/2;TW=TW/2
CX,CW=np.polynomial.legendre.leggauss(6);CX=(CX+1)/2;CW=CW/2
def polynomial(xy):
    return xy[:,0,None]**POWERS[:,0]*xy[:,1,None]**POWERS[:,1]

@njit(cache=True)
def phase_difference(pa,pb,coeff,powers):
    """Evaluate the phase with Horner's rule, avoiding large power matrices."""
    degree=int(np.max(powers[:,0]+powers[:,1]));matrix=np.zeros((degree+1,degree+1))
    for j in range(len(coeff)):matrix[powers[j,0],powers[j,1]]=coeff[j]
    out=np.zeros(len(pa))
    for k in range(len(pa)):
        va=0.;vb=0.
        for i in range(degree,-1,-1):
            ra=0.;rb=0.
            for j in range(degree-i,-1,-1):
                ra=ra*pa[k,1]+matrix[i,j];rb=rb*pb[k,1]+matrix[i,j]
            va=va*pa[k,0]+ra;vb=vb*pb[k,0]+rb
        out[k]=va-vb
    return out

def clip_polygon_reference(subject,clip):
    result=subject
    for j,a in enumerate(clip):
        b=clip[(j+1)%len(clip)];edge=b-a
        if not len(result):break
        output=[];prev=result[-1];dp=edge[0]*(prev[1]-a[1])-edge[1]*(prev[0]-a[0])
        for point in result:
            dn=edge[0]*(point[1]-a[1])-edge[1]*(point[0]-a[0])
            if (dp>=0)!=(dn>=0):output.append(prev+(point-prev)*dp/(dp-dn))
            if dn>=0:output.append(point)
            prev=point;dp=dn
        result=np.array(output)
    return result

@njit(cache=True)
def clip_polygon(subject,clip):
    """Convex pupil intersection with the same clipping rule as the reference."""
    capacity=len(subject)+len(clip)+8
    source=np.empty((capacity,2));target=np.empty((capacity,2))
    count=len(subject)
    source[:count]=subject
    for j in range(len(clip)):
        if count==0:break
        ax=clip[j,0];ay=clip[j,1]
        bx=clip[(j+1)%len(clip),0];by=clip[(j+1)%len(clip),1]
        ex=bx-ax;ey=by-ay
        px=source[count-1,0];py=source[count-1,1]
        dp=ex*(py-ay)-ey*(px-ax)
        outcount=0
        for i in range(count):
            qx=source[i,0];qy=source[i,1]
            dn=ex*(qy-ay)-ey*(qx-ax)
            if (dp>=0)!=(dn>=0):
                fraction=dp/(dp-dn)
                target[outcount,0]=px+(qx-px)*fraction
                target[outcount,1]=py+(qy-py)*fraction
                outcount+=1
            if dn>=0:
                target[outcount,0]=qx;target[outcount,1]=qy;outcount+=1
            px=qx;py=qy;dp=dn
        source,target=target,source
        count=outcount
    return source[:count].copy()

def pair_otf(coeff,na,wave,freq,axis,hull=None):
    shift=wave*.001*freq/na[axis]
    if shift>=2:return 0j
    if hull is not None:
        delta=np.zeros(2);delta[axis]=shift/2
        poly=clip_polygon(hull-delta,hull+delta)
        if len(poly)<3:return 0j
        ctr=poly.mean(0);va=poly-ctr;vb=np.roll(poly,-1,axis=0)-ctr
        cross=va[:,0]*vb[:,1]-va[:,1]*vb[:,0]
        points=ctr+TX[None,:,None,None]*((1-TX)[None,None,:,None]*va[:,None,None,:]+TX[None,None,:,None]*vb[:,None,None,:])
        area=abs(np.sum(hull[:,0]*np.roll(hull[:,1],-1)-hull[:,1]*np.roll(hull[:,0],-1)))/2
        weights=(cross[:,None,None]*TX[None,:,None]*TW[None,:,None]*TW[None,None,:]/area).ravel()
        points=points.reshape(-1,2);pa=points+delta;pb=points-delta
        phase=phase_difference(pa,pb,coeff,POWERS)
        return np.sum(weights*np.exp(2j*np.pi*phase/(wave*.001)))
    half=1-shift/2
    xx=np.r_[(GX+1)*half/2,-(GX+1)*half/2]
    width=np.sqrt(1-(abs(xx)+shift/2)**2)
    yy=width[:,None]*GX[None,:]
    weights=np.r_[GW,GW][:,None]*half/2*width[:,None]*GW[None,:]/np.pi
    points=np.c_[np.broadcast_to(xx[:,None],yy.shape).ravel(),yy.ravel()]
    if axis==1:points=points[:,::-1]
    pa=points.copy();pb=points.copy();pa[:,axis]+=shift/2;pb[:,axis]-=shift/2
    phase=phase_difference(pa,pb,coeff,POWERS)
    return np.sum(weights.ravel()*np.exp(2j*np.pi*phase/(wave*.001)))

def phase_centroid_y(coeff,na,hull):
    """PSF first moment from the wavefront boundary integral, in mm.

    Green's theorem evaluates the pupil-average phase gradient. Six-point
    edge quadrature exactly integrates the degree-ten fitted polynomial.
    """
    nxt=np.roll(hull,-1,axis=0);edge=nxt-hull
    points=(hull[:,None,:]+CX[None,:,None]*edge[:,None,:]).reshape(-1,2)
    phase=phase_difference(points,np.zeros_like(points),coeff,POWERS).reshape(len(hull),len(CX))
    area=abs(np.sum(hull[:,0]*nxt[:,1]-hull[:,1]*nxt[:,0]))/2
    return float(np.sum(edge[:,0,None]*CW[None,:]*phase)/(area*na[1]))

def quantities(c,t,fields=np.array([0.,3.,6.3,9.08]),wd=500.,frequencies=(225.,180.),fixed_stop_radius=None,low_frequency=True,wavelengths=None,return_complex=False,caps=None,return_diagnostics=False):
    p=e.paraxial(c,t,wd);radius=p['efl']/5.6
    marginal,_,_=e.ray(c,t,e.N[1],wd,0.,0.,radius,p['ep'])
    stop_radius=abs(marginal[9]) if fixed_stop_radius is None else fixed_stop_radius
    wave_set=e.WAVES if wavelengths is None else np.asarray(wavelengths)
    refractive_indices=e.N if wavelengths is None else e.indices(wave_set)
    mtfs=[];complex_values=[];angular_areas=[];clipped_counts=[];centroid_values=[];maxr=np.zeros(15);chief_images=[];invalid=0
    for field in fields:
        data=[];field_clipped=[]
        for n in refractive_indices:
            rr,hh,v=e.trace_bundle_stop_aiming(c,t,n,wd,field/p['magnification'],p['ep'],PUPIL,radius,stop_radius)
            invalid+=int((~v).sum());maxr=np.maximum(maxr,hh.max(0))
            if caps is None:blocked=np.zeros(len(v),bool)
            else:
                ray_caps=np.asarray(caps,float).copy()
                # The physical iris is already imposed by stop aiming. Its
                # radius may be adjusted for F/2.8 in this configuration.
                ray_caps[5]=np.inf
                blocked=~np.all(hh<=ray_caps+1e-7,axis=1)
            field_clipped.append(int(np.sum(v&blocked)));v=v&~blocked;data.append((rr,v))
        primary_data=data[min(1,len(data)-1)]
        ref=primary_data[0][primary_data[1],:2].mean(0) if primary_data[1].any() else np.zeros(2)
        chief=e.aim_chief(c,t,e.N[1],wd,field/p['magnification'],p['ep'])
        rc,_,v=e.ray(c,t,e.N[1],wd,field/p['magnification'],0.,chief,p['ep']);chief_images.append(rc[1])
        wave_otf=[];field_areas=[];field_centroids=[]
        for wave,(r,v) in zip(wave_set,data):
            if v.sum()<70:wave_otf.append(np.zeros(4 if low_frequency else 2,complex));field_areas.append(0.);field_centroids.append(0.);continue
            r=r[v];vz=r[:,11];vx=r[:,2]*vz;vy=r[:,3]*vz;xy=np.c_[vx,vy]
            ctr=(xy.max(0)+xy.min(0))/2;na=(xy.max(0)-xy.min(0))/2
            q=(sum(t)-r[:,7])/vz
            phase=r[:,4]-q+vx*(ref[0]-r[:,5])+vy*(ref[1]-r[:,6])+vz*(sum(t)-r[:,7]);phase-=phase.mean()
            coeff=np.linalg.lstsq(polynomial((xy-ctr)/na),phase,rcond=1e-9)[0]
            fc=frequencies[0] if field==0 else frequencies[1]
            points=(xy-ctr)/na
            if return_diagnostics:field_areas.append(float(ConvexHull(xy).volume))
            hull=points[ConvexHull(points).vertices] if field>0 else None
            if return_diagnostics:field_centroids.append(phase_centroid_y(coeff,na,points[ConvexHull(points).vertices] if hull is None else hull))
            vv=[pair_otf(coeff,na,wave,fc,1,hull),pair_otf(coeff,na,wave,fc,0,hull)]
            if low_frequency:vv.extend([pair_otf(coeff,na,wave,90,1,hull),pair_otf(coeff,na,wave,90,0,hull)])
            wave_otf.append(vv)
        mtfs.append(abs(np.mean(wave_otf,axis=0)))
        if return_complex:complex_values.append(wave_otf)
        if return_diagnostics:angular_areas.append(field_areas);clipped_counts.append(field_clipped);centroid_values.append(field_centroids)
    if return_diagnostics:
        areas=np.array(angular_areas);primary=min(1,len(wave_set)-1)
        diagnostics=dict(relative_illumination_primary=(areas[:,primary]/max(1e-16,areas[0,primary])).tolist(),
                         angular_pupil_areas=areas.tolist(),clipped_ray_counts=clipped_counts,
                         phase_centroids_y_mm=centroid_values,phase_centroid_colour_span_um=(1000*np.ptp(np.array(centroid_values),axis=1)).tolist())
        return np.array(mtfs),maxr,np.array(chief_images),invalid,diagnostics
    if return_complex:return np.array(mtfs),maxr,np.array(chief_images),invalid,np.array(complex_values)
    return np.array(mtfs),maxr,np.array(chief_images),invalid

def blank_clearance(c,t,caps):
    minima=[]
    def sag(s,r):
        disc=1-c[s]**2*r*r
        value=c[s]*r*r/(1+np.sqrt(np.maximum(.001,disc)))
        j=15+3*SURFS.index(s+1) if s+1 in SURFS else (33 if s==1 else (36 if s==8 else -1))
        if j>=0 and len(c)>=j+3:value+=c[j]*r**4+c[j+1]*r**6+c[j+2]*r**8
        return value
    for s in range(14):
        radius=caps[s] if s in [4] or e.MEDIA[s]>=0 else (caps[6] if s==5 else min(caps[s:s+2]))
        r=np.linspace(0,radius,33)
        required=.5 if e.MEDIA[s]>=0 else (.3 if s in [4,5] else .1)
        minima.append(np.min(t[s]+sag(s+1,r)-sag(s,r))-required)
    return np.array(minima)

def run(tag='candidate_v12',steps=150):
    schema=json.loads((ROOT/(tag+'_schema.json')).read_text(encoding='utf8'))
    load_glasses([(g['catalog'],g['name']) for g in schema['glasses']])
    c=np.r_[schema['curvatures'],np.zeros(18)];t=np.array(schema['thicknesses'])
    for asp in schema['aspheres']:
        j=15+3*SURFS.index(asp['surface']);c[j:j+3]=asp['coefficients']
    caps=np.array([schema['fixed_semi_diameters'].get(str(i+1),10.) for i in range(15)])
    scales=np.tile([1e-5,1e-7,1e-9],6)
    x=np.r_[c[15:]/scales,t[-1]]
    calls=0;best=np.inf
    def unpack(xx):
        cc=c.copy();tt=t.copy();cc[15:]=xx[:18]*scales;tt[-1]=xx[18];return cc,tt
    def fun(xx):
        nonlocal best,calls
        cc,tt=unpack(xx);m,radii,chiefs,invalid=quantities(cc,tt)
        target=np.array([[.42,.42,.72,.72],[.26,.32,.60,.60],[.26,.30,.60,.60],[.32,.32,.60,.60]])
        weights=np.array([[2.,2.,.2,.2],[.35,.35,.2,.2],[.35,.35,.2,.2],[2.,2.,.2,.2]])
        lower=np.array([[.33,.33],[.16,.16],[.16,.16],[.245,.245]])
        parts=[((m-target)*weights).ravel(),(np.maximum(0,lower-m[:,:2])*15).ravel(),
               np.maximum(0,-blank_clearance(cc,tt,caps))*100,
               np.maximum(0,radii-caps+.02)*50,
               np.maximum(0,abs(chiefs-np.array([0,3,6.3,9.08]))-.0027*np.array([0,3,6.3,9.08]))*100,
               np.array([invalid*10.])]
        rr=np.concatenate(parts);score=rr@rr;calls+=1
        if score<best:best=score;np.save(ROOT/'otf_best_parameters.npy',xx)
        if calls%400==0:print('STEP',calls,best,m.tolist(),flush=True)
        return rr
    print('INITIAL_SURROGATE',quantities(c,t)[0].tolist(),flush=True)
    lo=np.r_[np.full(18,-20.),t[-1]-.08];hi=np.r_[np.full(18,20.),t[-1]+.08]
    opt=least_squares(fun,x,jac=lambda xx:approx_derivative(fun,xx,method='3-point',abs_step=2e-5*(1+abs(xx)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',max_nfev=steps,ftol=1e-9,xtol=1e-9,gtol=1e-7)
    cc,tt=unpack(opt.x)
    schema['thicknesses']=tt.tolist();schema['aspheres']=[dict(surface=s,coefficients=cc[15+3*j:18+3*j].tolist()) for j,s in enumerate(SURFS)]
    (ROOT/'candidate_v13_schema.json').write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
    print('DONE',opt.status,opt.cost,opt.nfev,quantities(cc,tt)[0].tolist(),flush=True)
    import verify_design as v
    v.N=e.N
    audit=v.diffraction(cc,tt,fields=np.array([0.,3.,6.3,9.08]),grid=192)
    (ROOT/'candidate_v13_independent_mtf.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
    print('EXACT',[(d['MTF_225_tangential'],d['MTF_180_tangential'],d['MTF_180_sagittal']) for d in audit['fields']],flush=True)
if __name__=='__main__':run()
