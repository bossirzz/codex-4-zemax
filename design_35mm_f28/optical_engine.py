"""Independent sequential spherical optics; lengths mm, wavelengths micrometres.

Coordinates use the usual OpticStudio radius sign. No licensing/API access.
This is a design aid, not a replacement for native optical validation.
"""
import json
import math
from pathlib import Path
import sys
import numpy as np
from numba import njit
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from match_patent_glasses import index, parse_agf, enrich, CATALOG_DIR

WAVES = np.array([.4861327, .5875618, .6562725])
MEDIA = np.array([0,-1,1,2,-1,-1,3,4,-1,5,-1,6,-1,7,-1])
CURVED = np.array([i for i in range(15) if i != 5])
RADII = np.array([27.3,223.5,11.7,256.,10.3,np.inf,-12.1,16.4,-16.4,79.8,-41.9,-2074.,-37.4,-22.,111.2])
THICKNESS = np.array([5.35,.10,2.6,2.33,3.60,1.35,.70,3.6,.10,5.1,1.74,5.1,4.18,3.5,8.73])

def select_glass():
    pairs=[('CDGM','H-LAF54'),('CDGM','H-FK61B'),('CDGM','H-LAF55'),
           ('CDGM','H-QF8'),('NHG','H-FK69'),('CDGM','H-LAF10LA'),
           ('CDGM','H-ZLAF53B'),('CDGM','H-ZF1')]
    profile=ROOT/'active_glass_profile.json'
    if profile.exists():pairs=[tuple(pair) for pair in json.loads(profile.read_text(encoding='utf8'))]
    allg={}
    for cat in {p[0] for p in pairs}:
        path=next(p for p in CATALOG_DIR.glob('*.agf') if p.stem.upper()==cat)
        for g in parse_agf(path)[0]:
            allg[(cat,g['name'])]=enrich(g)
    glasses=[allg[p] for p in pairs]
    ROOT.joinpath('design_glasses.json').write_text(json.dumps(glasses,ensure_ascii=False,indent=2),encoding='utf8')
    return glasses

GLASSES = select_glass()

def indices(waves=WAVES):
    result=np.ones((len(waves),15))
    for wi,w in enumerate(waves):
        for s,g in enumerate(MEDIA):
            if g>=0: result[wi,s]=index(GLASSES[g],float(w))
    return result

N = indices()

@njit(cache=True)
def matrix(c,t,n,end=15):
    a,b,cc,d=1.,0.,0.,1.
    ni=1.
    for s in range(end):
        nn=n[s]; power=(nn-ni)*c[s]
        cc,d=cc-power*a,d-power*b
        if s<end-1:
            a,b=a+t[s]/nn*cc,b+t[s]/nn*d
        ni=nn
    return np.array([a,b,cc,d])

def paraxial(c,t,wd=500.,n=None):
    if n is None:n=N[1]
    a,b,cc,d=matrix(c,t,n)
    ap,bp,_,_=matrix(c,t,n,6)
    q=-a/cc if np.isinf(wd) else -(a*wd+b)/(cc*wd+d)
    m=(a+q*cc) if np.isinf(wd) else 1./(cc*wd+d)
    return dict(efl=-1/cc,focus=q,magnification=m,ep=bp/ap,stop_scale=ap,ttl=sum(t),matrix=[a,b,cc,d])

@njit(cache=True)
def ray(c,t,n,wd,oy,px,py,ep):
    # Launch through a point in the paraxial entrance-pupil plane. Chief aiming
    # changes px/py; all other rays retain uniform entrance-pupil sampling.
    length=wd+ep
    vx=px/length;vy=(py-oy)/length;vz=1.
    norm=math.sqrt(vx*vx+vy*vy+1);vx/=norm;vy/=norm;vz/=norm
    z=-40.;x=px*(wd+z)/length;y=oy+(py-oy)*(wd+z)/length
    # OPL from the object is needed for wavefront validation.
    opl=(wd+z)/vz
    vertex=0.;ni=1.;maxrho=0.;stopx=0.;stopy=0.
    heights=np.zeros(15)
    for s in range(15):
        cv=c[s]
        if abs(cv)<1e-10:
            distance=(vertex-z)/vz
            nx,ny,nz=0.,0.,-1.
        else:
            r=1/cv;cz=vertex+r
            bx=x;by=y;bz=z-cz
            bb=bx*vx+by*vy+bz*vz
            disc=bb*bb-(bx*bx+by*by+bz*bz)+r*r
            if disc<0.: return np.zeros(12),heights,False
            distance=-bb-math.sqrt(disc) if r>0 else -bb+math.sqrt(disc)
            nx=(x+distance*vx)/abs(r);ny=(y+distance*vy)/abs(r);nz=(z+distance*vz-cz)/abs(r)
        j=-1
        if s==0:j=15
        elif s==14:j=18
        elif s==9:j=21
        elif s==11:j=24
        elif s==2:j=27
        elif s==6:j=30
        elif s==1:j=33
        elif s==8:j=36
        if j>=0 and len(c)>=j+3:
            a4,a6,a8=c[j],c[j+1],c[j+2]
            if a4!=0. or a6!=0. or a8!=0.:
                for _ in range(8):
                    ax=x+distance*vx;ay=y+distance*vy;ar2=ax*ax+ay*ay
                    dd=1-cv*cv*ar2
                    if dd<=0:return np.zeros(12),heights,False
                    root=math.sqrt(dd)
                    sag=cv*ar2/(1+root)+a4*ar2**2+a6*ar2**3+a8*ar2**4
                    grad=cv/root+4*a4*ar2+6*a6*ar2**2+8*a8*ar2**3
                    error=z+distance*vz-vertex-sag
                    derivative=vz-grad*(ax*vx+ay*vy)
                    if abs(derivative)<1e-8:return np.zeros(12),heights,False
                    delta=error/derivative;distance-=delta
                    if abs(delta)<1e-10:break
                ax=x+distance*vx;ay=y+distance*vy;ar2=ax*ax+ay*ay
                root=math.sqrt(max(1e-9,1-cv*cv*ar2))
                grad=cv/root+4*a4*ar2+6*a6*ar2**2+8*a8*ar2**3
                nnorm=math.sqrt(1+grad*grad*ar2)
                nx=grad*ax/nnorm;ny=grad*ay/nnorm;nz=-1/nnorm
        # Slightly negative distances can occur for touching cemented surfaces.
        if distance < -1e-5:return np.zeros(12),heights,False
        x+=distance*vx;y+=distance*vy;z+=distance*vz;opl+=distance*ni
        dot=nx*vx+ny*vy+nz*vz
        if dot>0.:nx=-nx;ny=-ny;nz=-nz;dot=-dot
        nn=n[s];eta=ni/nn;k=1-eta*eta*(1-dot*dot)
        if k<0.:return np.zeros(12),heights,False
        scale=eta*(-dot)-math.sqrt(k)
        vx=eta*vx+scale*nx;vy=eta*vy+scale*ny;vz=eta*vz+scale*nz
        # Analytical vector Snell preserves unit length to floating precision.
        rho=math.sqrt(x*x+y*y);heights[s]=rho
        maxrho=max(maxrho,rho)
        if s==5:stopx=x;stopy=y
        vertex+=t[s]
        ni=nn
    q=(vertex-z)/vz
    out=np.array([x+q*vx,y+q*vy,vx/vz,vy/vz,opl+q,
                  x,y,z,stopx,stopy,maxrho,vz])
    return out,heights,True

@njit(cache=True)
def aim_chief(c,t,n,wd,oy,ep):
    py=0.
    for _ in range(5):
        r,_,ok=ray(c,t,n,wd,oy,0.,py,ep)
        r2,_,ok2=ray(c,t,n,wd,oy,0.,py+.002,ep)
        if not(ok and ok2):return 0.
        derivative=(r2[9]-r[9])/.002
        if abs(derivative)<1e-6:return py
        step=r[9]/derivative
        py-=step
        if abs(step)<1e-8:break
    return py

@njit(cache=True)
def trace_bundle(c,t,n,wd,oy,ep,pupil,radius):
    chief=aim_chief(c,t,n,wd,oy,ep)
    out=np.zeros((len(pupil),12));heights=np.zeros((len(pupil),15));valid=np.ones(len(pupil),np.bool_)
    for i in range(len(pupil)):
        out[i],heights[i],valid[i]=ray(c,t,n,wd,oy,radius*pupil[i,0],chief+radius*pupil[i,1],ep)
    return out,heights,valid

@njit(cache=True)
def trace_bundle_stop_aiming(c,t,n,wd,oy,ep,pupil,radius,stop_radius):
    """Aim every sampled ray at the specified circular physical stop."""
    chief=aim_chief(c,t,n,wd,oy,ep)
    out=np.zeros((len(pupil),12));heights=np.zeros((len(pupil),15));valid=np.ones(len(pupil),np.bool_)
    for i in range(len(pupil)):
        px=radius*pupil[i,0];py=chief+radius*pupil[i,1]
        targetx=stop_radius*pupil[i,0];targety=stop_radius*pupil[i,1]
        for _ in range(6):
            r,h,ok=ray(c,t,n,wd,oy,px,py,ep)
            if not ok:break
            ex=r[8]-targetx;ey=r[9]-targety
            if abs(ex)+abs(ey)<1e-9:break
            rx,_,okx=ray(c,t,n,wd,oy,px+.001,py,ep)
            ry,_,oky=ray(c,t,n,wd,oy,px,py+.001,ep)
            if not(okx and oky):break
            jxx=(rx[8]-r[8])/.001;jyx=(rx[9]-r[9])/.001
            jxy=(ry[8]-r[8])/.001;jyy=(ry[9]-r[9])/.001
            det=jxx*jyy-jxy*jyx
            if abs(det)<1e-9:break
            px-=(jyy*ex-jxy*ey)/det;py-=(-jyx*ex+jxx*ey)/det
        out[i],heights[i],valid[i]=ray(c,t,n,wd,oy,px,py,ep)
        valid[i]=valid[i] and abs(out[i,8]-targetx)+abs(out[i,9]-targety)<1e-6
    return out,heights,valid

def pupil_sample(rings=4,arms=12):
    # Equal-area radial quadrature; circle centre is separately traced.
    rr=np.sqrt((np.arange(rings)+.5)/rings)
    return np.array([(r*np.cos(a),r*np.sin(a)) for r in rr for a in np.arange(arms)*2*np.pi/arms])

PUPIL=np.vstack([pupil_sample(),np.array([(np.cos(a),np.sin(a)) for a in np.arange(12)*2*np.pi/12])])
FIELDS=np.array([0.,3.,6.3,9.])

def physical_clearances(c,t,ray_radii,margin=.20):
    """Shared blank diameters for singlets and cemented groups, not ray-only caps."""
    groups=[[0,1],[2,3,4],[5],[6,7,8],[9,10],[11,12],[13,14]]
    blank=np.zeros(15)
    for group in groups:blank[group]=max(ray_radii[group])+margin
    def sag(s,r):
        disc=1-c[s]**2*r*r
        value=c[s]*r*r/(1+np.sqrt(max(.001,disc)))
        if len(c)>15 and (s in [0,14] or (s==9 and len(c)>21) or (s==11 and len(c)>24) or (s==2 and len(c)>27) or (s==6 and len(c)>30)):
            j=15 if s==0 else (18 if s==14 else (21 if s==9 else (24 if s==11 else (27 if s==2 else 30))))
            value+=c[j]*r**4+c[j+1]*r**6+c[j+2]*r**8
        return value
    clear=[]
    for s in range(14):
        if MEDIA[s]>=0:
            r=blank[s];needed=.60
        elif s==4:
            r=blank[4];needed=.30
        elif s==5:
            r=blank[6];needed=.30
        else:
            r=min(blank[s],blank[s+1]);needed=.10
        # Check the entire radial profile; edge-only checks miss minima between
        # the vertex and rim of a polynomial aspheric surface.
        sep=min(t[s]+sag(s+1,rr)-sag(s,rr) for rr in np.linspace(0,r,17))
        clear.append(sep-needed)
    return np.array(clear),blank

def unpack(x):
    c=np.zeros(24 if len(x)>39 else (21 if len(x)>29 else 15));c[CURVED]=x[:14]*.05
    if len(x)>29:c[15:21]=x[29:35]*np.tile([1e-5,1e-7,1e-9],2)
    if len(x)>39:c[21:24]=x[39:42]*np.array([1e-5,1e-7,1e-9])
    return c,x[14:29].copy()

def pack(c,t):
    xx=np.r_[c[CURVED]/.05,t]
    return np.r_[xx,c[15:21]/np.tile([1e-5,1e-7,1e-9],2)] if len(c)>15 else xx

def focus_group_gap(c,t,wd,gap_index=10):
    trial=t.copy();trial[gap_index]=0.
    a,b,cc,d=matrix(c,trial,N[1]);f0=a*wd+b+t[-1]*(cc*wd+d)
    trial[gap_index]=1.
    a,b,cc,d=matrix(c,trial,N[1]);f1=a*wd+b+t[-1]*(cc*wd+d)
    return -f0/(f1-f0)

def analyze(c,t,wd=500.,waves=None,pupil=PUPIL,fields=FIELDS):
    ns=N if waves is None else indices(waves)
    p=paraxial(c,t,wd,ns[min(1,len(ns)-1)])
    oy=fields/p['magnification']
    bundles=[]
    for wi,n in enumerate(ns):
        for f,y in zip(fields,oy):
            r,h,v=trace_bundle(c,t,n,wd,y,p['ep'],pupil,p['efl']/5.6)
            bundles.append((wi,f,r,h,v))
    return p,bundles

def residual(x,focus_weight=0.,stage=1):
    c,t=unpack(x);p=paraxial(c,t)
    if p['efl']<=0 or not(np.isfinite(p['focus'])):return np.full(2367,1e3)
    p,bundles=analyze(c,t)
    pieces=[]
    maxdiam=np.zeros(15)
    centroids=np.zeros((3,4,2));chiefs=np.zeros((3,4,2))
    for wi,f,r,h,v in bundles:
        fi=np.argmin(abs(FIELDS-f));centroids[wi,fi]=np.mean(r[v,:2],axis=0) if np.any(v) else 0.
    for wi,f,r,h,v in bundles:
        fi=np.argmin(abs(FIELDS-f))
        maxdiam=np.maximum(maxdiam,np.max(h,axis=0))
        centroid=centroids[wi,fi]
        centroids[wi,fi]=centroid
        fieldweight=.25
        pieces.append(((r[:,:2]-centroid)*v[:,None]*fieldweight*1000/np.sqrt(len(r))).ravel())
        pieces.append(np.where(v,0.,100.))
        # Wavefront at a common reference point includes lateral colour.
        ref=centroids[1,fi];vz=r[:,11];vx=r[:,2]*vz;vy=r[:,3]*vz
        q=(sum(t)-r[:,7])/np.maximum(vz,1e-9)
        phase=r[:,4]-q+vx*(ref[0]-r[:,5])+vy*(ref[1]-r[:,6])+vz*(sum(t)-r[:,7])
        phase=np.where(v,phase-np.mean(phase[v]) if np.any(v) else 0.,0.)
        weight=2.0 if f>=9. else 1.5
        pieces.append(phase*50000.*weight/np.sqrt(len(r)))
    for fi in range(4):
        pieces.append(((centroids[:,fi]-centroids[1,fi])*1000).ravel())
    pieces.append(np.array([(p['efl']-35)*100.,(p['focus']-t[-1])*focus_weight]))
    # Distortion relative to paraxial image height at the same conjugate.
    dist=centroids[1,:,1]-FIELDS
    pieces.append(np.maximum(0.,abs(dist)-.0028*FIELDS)*2000)
    # Infinity chief-ray angle at the full 9 mm image height.
    ti=t.copy()
    if len(x)>37:ti[-1]+=x[38]*4./3.
    ti[10]=focus_group_gap(c,ti,1e8)
    pi=paraxial(c,ti,np.inf)
    ci=aim_chief(c,ti,N[1],1e7,-9e7/pi['efl'],pi['ep'])
    rayi,_,oki=ray(c,ti,N[1],1e7,-9e7/pi['efl'],0.,ci,pi['ep'])
    cra=abs(np.arctan(rayi[3]))
    pieces.append(np.array([max(0.,cra-np.deg2rad(9.5))*1200.,max(0.,sum(t)-58.5)*20]))
    pieces.append(np.maximum(0.,maxdiam-11.15)*100)
    # Lens edge thickness / inter-element air clearance at full used aperture.
    sagittae=np.zeros(15)
    for s in range(15):
        rho=maxdiam[s]+.15
        disc=1-c[s]**2*rho**2
        sagittae[s]=c[s]*rho*rho/(1+np.sqrt(max(.001,disc)))
        if len(c)>15 and (s in [0,14] or (s==9 and len(c)>21)):
            j=15 if s==0 else (18 if s==14 else 21)
            sagittae[s]+=c[j]*rho**4+c[j+1]*rho**6+c[j+2]*rho**8
    clearance=[]
    for s in range(14):
        needed=.65 if MEDIA[s]>=0 else .05
        # Do not constrain the stop-to-lens gaps as if the stop were glass.
        if s in [4,5]:needed=.02
        clearance.append(max(0.,needed-(t[s]+sagittae[s+1]-sagittae[s]))*3000)
    pieces.append(np.array(clearance))
    g200=focus_group_gap(c,t,200.);ginf=focus_group_gap(c,t,1e8)
    max_ttl=sum(t)+max(g200,ginf)-t[10]
    pieces.append(np.array([max(0.,max_ttl-59.5)*100.,max(0.,.15-min(g200,ginf))*100.]))
    if len(c)>15:
        positions=[]
        # Front six and rear two form separately moving focusing groups; sensor
        # stays fixed. Rear motion changes the physical air gap to the sensor.
        for wi_config,(wd,offset,weight) in enumerate([(200.,x[35],.6),(2000.,x[36],.6)]):
            tt=t.copy()
            if len(x)>37:tt[-1]+=x[37+wi_config]
            tt[10]=focus_group_gap(c,tt,wd)+offset
            positions.append(tt[10]+tt[-1])
            pp,bb=analyze(c,tt,wd,fields=np.array([0.,6.3,9.]))
            ccs=np.zeros((3,3,2))
            for wi,f,r,h,v in bb:ccs[wi,np.argmin(abs(np.array([0.,6.3,9.])-f))]=np.mean(r[v,:2],axis=0) if np.any(v) else 0.
            maxd=np.zeros(15)
            for wi,f,r,h,v in bb:
                fi=np.argmin(abs(np.array([0.,6.3,9.])-f));ref=ccs[1,fi]
                vz=r[:,11];vx=r[:,2]*vz;vy=r[:,3]*vz
                q=(sum(tt)-r[:,7])/np.maximum(vz,1e-9)
                phase=r[:,4]-q+vx*(ref[0]-r[:,5])+vy*(ref[1]-r[:,6])+vz*(sum(tt)-r[:,7])
                phase=np.where(v,phase-np.mean(phase[v]) if np.any(v) else 0.,0.)
                pieces.append(phase*50000.*weight/np.sqrt(len(r)))
                pieces.append(np.where(v,0.,100.))
                maxd=np.maximum(maxd,h.max(0))
            pieces.append(np.maximum(0.,maxd-11.15)*100.)
            maxdiam=np.maximum(maxdiam,maxd)
            distnear=ccs[1,:,1]-np.array([0.,6.3,9.])
            pieces.append(np.maximum(0.,abs(distnear)-.0058*np.array([0.,6.3,9.]))*2000)
            # The plane of the stop must lie outside the preceding glass even at
            # the full pupil boundary. Centre-only rays cannot prove clearance.
            rstop=maxd[4]+.15
            sagstop=c[4]*rstop*rstop/(1+np.sqrt(max(.001,1-c[4]**2*rstop**2)))
            pieces.append(np.array([max(0.,sagstop+.05-t[4])*3000.]))
            pieces.append(np.array([max(0.,sum(tt)-59.5)*1000.,max(0.,8.-tt[-1])*1000.,max(0.,.05-tt[10])*1000.]))
        travel=abs(positions[0]-positions[1])
        pieces.append(np.array([max(0.,travel-6.2)*300.]))
    geo,blank=physical_clearances(c,t,maxdiam)
    pieces.append(np.maximum(0.,-geo)*5000.)
    return np.concatenate(pieces)

def metrics(c,t):
    p,bundles=analyze(c,t,pupil=pupil_sample(8,24),fields=FIELDS)
    result=dict(paraxial=p,fields=[])
    maxr=np.zeros(15)
    for fi,f in enumerate(FIELDS):
        rays=[b[2] for b in bundles if b[1]==f]
        rr=np.concatenate([r[:,:2] for r in rays]);avg=rr.mean(0)
        for b in bundles:
            if b[1]==f:maxr=np.maximum(maxr,b[3].max(0))
        result['fields'].append(dict(image_height=float(f),rms_um=float(np.sqrt(np.mean(np.sum((rr-avg)**2,axis=1)))*1000),
            centroid=avg.tolist(),distortion_percent=float((avg[1]/f-1)*100) if f else 0,
            valid_rays=int(sum(np.sum(b[4]) for b in bundles if b[1]==f))))
    result['ray_semidiameters']=maxr.tolist()
    return result

def optimize(max_nfev=400,continue_file=None,aspheric=False,floating=False,third=False):
    c=1/RADII;t=THICKNESS.copy()
    scale=35/paraxial(c,t)['efl'];c/=scale;t*=scale
    t[-1]=paraxial(c,t)['focus']
    x=pack(c,t)
    if continue_file:x=np.load(continue_file)
    if aspheric and len(x)==29:x=np.r_[x,np.zeros(8)]
    if floating and len(x)==37:x=np.r_[x,np.zeros(2)]
    if third and len(x)==39:x=np.r_[x,np.zeros(3)]
    lo=np.r_[np.full(14,-3.),[1.4,.05,1.4,.65,.3,.3,.65,1.4,.05,1.4,.25,1.2,.1,.65,8.]]
    hi=np.r_[np.full(14,3.),[7.,2.,6.,4.,8.,8.,4.,7.,2.,8.,8.,8.,8.,5.,22.]]
    if len(x)>29:
        lo=np.r_[lo,np.full(6,-5.),-.2,-.2];hi=np.r_[hi,np.full(6,5.),.2,.2]
    if len(x)>37:
        lo=np.r_[lo,-5.,-5.];hi=np.r_[hi,5.,5.]
    if len(x)>39:
        lo=np.r_[lo,-5.,-5.,-5.];hi=np.r_[hi,5.,5.,5.]
    x=np.clip(x,lo+1e-6,hi-1e-6)
    counter=0
    best=np.inf
    def fun(xx):
        nonlocal counter,best
        r=residual(xx);score=np.dot(r,r);counter+=1
        if score<best:
            best=score;np.save(ROOT/'best_parameters.npy',xx)
        if counter%2000==0:
            c1,t1=unpack(xx);m=metrics(c1,t1)
            print(json.dumps(dict(evaluation=counter,best=best,efl=m['paraxial']['efl'],ttl=sum(t1),rms=[f['rms_um'] for f in m['fields']]),ensure_ascii=False),flush=True)
        return r
    print(json.dumps(dict(initial=metrics(*unpack(x))),ensure_ascii=False),flush=True)
    from scipy.optimize._numdiff import approx_derivative
    def jac(xx):
        return approx_derivative(fun,xx,method='3-point',abs_step=1e-6*(1+abs(xx)),bounds=(lo,hi))
    opt=least_squares(fun,x,bounds=(lo,hi),max_nfev=max_nfev,jac=jac,x_scale='jac',ftol=2e-9,xtol=2e-10,gtol=2e-6)
    c,t=unpack(opt.x)
    np.save(ROOT/'optimized_parameters.npy',opt.x)
    result=dict(status=int(opt.status),message=opt.message,nfev=opt.nfev,cost=opt.cost,curvatures=c.tolist(),thicknesses=t.tolist(),metrics=metrics(c,t))
    ROOT.joinpath('optimized_design.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return opt

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--steps',type=int,default=300);ap.add_argument('--continue-file');ap.add_argument('--aspheric',action='store_true');ap.add_argument('--floating',action='store_true');ap.add_argument('--third',action='store_true')
    args=ap.parse_args();optimize(args.steps,args.continue_file,args.aspheric,args.floating,args.third)
