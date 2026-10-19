"""Best-conjugate redesign, keeping real catalog glasses and full blank clearance."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative
import optical_engine as e
from match_patent_glasses import parse_agf, enrich, CATALOG_DIR

ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')

def load_glasses(pairs):
    db={}
    for cat in {p[0] for p in pairs}:
        for g in parse_agf(CATALOG_DIR/(cat+'.agf'))[0]:db[cat,g['name']]=enrich(g)
    e.GLASSES=[db[tuple(p)] for p in pairs]
    e.N=e.indices()

def unpack(x):
    return e.unpack(np.r_[x[:35],np.zeros(4),x[35:38]])

def residual(x):
    c,t=unpack(x)
    p,bundles=e.analyze(c,t,fields=e.FIELDS)
    parts=[];radii=np.zeros(15);centroids=np.zeros((3,len(e.FIELDS),2))
    for wi,f,r,h,v in bundles:
        fi=np.argmin(abs(e.FIELDS-f));centroids[wi,fi]=np.mean(r[v,:2],axis=0) if v.any() else 0
    for wi,f,r,h,v in bundles:
        fi=np.argmin(abs(e.FIELDS-f));ref=centroids[1,fi]
        radii=np.maximum(radii,h.max(0))
        vz=r[:,11];vx=r[:,2]*vz;vy=r[:,3]*vz
        q=(sum(t)-r[:,7])/np.maximum(vz,1e-9)
        phase=r[:,4]-q+vx*(ref[0]-r[:,5])+vy*(ref[1]-r[:,6])+vz*(sum(t)-r[:,7])
        phase=np.where(v,phase-np.mean(phase[v]) if v.any() else 0,0)
        parts.append(phase*60000/np.sqrt(len(r)))
        parts.append(((r[:,:2]-ref)*v[:,None]*1000*.10/np.sqrt(len(r))).ravel())
        parts.append(np.where(v,0.,1000.))
    parts.append(np.array([(p['efl']-35)*300.]))
    parts.append(((centroids-centroids[1])*500).ravel())
    parts.append(np.maximum(0,abs(centroids[1,:,1]-e.FIELDS)-.0027*e.FIELDS)*5000)
    # Full mechanical profiles, including the iris plane, not individual ray caps.
    clear,blank=e.physical_clearances(c,t,radii,margin=.20)
    parts.append(np.maximum(0,-clear)*5000)
    parts.append(np.maximum(0,blank-11.45)*3000)
    ti=t.copy();ti[-1]=e.paraxial(c,t,np.inf)['focus']
    pi=e.paraxial(c,ti,1e7)
    chief=e.aim_chief(c,ti,e.N[1],1e7,9.05/pi['magnification'],pi['ep'])
    rr,hh,valid=e.ray(c,ti,e.N[1],1e7,9.05/pi['magnification'],0.,chief,pi['ep'])
    parts.append(np.array([max(0,abs(np.arctan(rr[3]))-np.deg2rad(9.2))*1500,
                           max(0,sum(t)-57)*1000]))
    return np.concatenate(parts)

def main(seed='gauss',steps=300):
    e.FIELDS=np.array([0.,3.,6.3,9.05])
    if seed=='gauss':
        pairs=[('SCHOTT',n) for n in ['N-SK2','N-SK16','F5','F5','N-SK16','N-SK16','N-LAK34','N-SF5']]
        load_glasses(pairs)
        scale=35/99.50067458127944
        c=np.array([.018466113683,.0065564346019,.027815928521,0,.044903609561,0,-.038933179444,0,-.027041482725,.0050912003495,-.01489257612,0,0,0,0])/scale
        t=np.array([8.74665785,.5,14,3.77696589,14.2530593,12.4281291,3.77696589,10.8339285,.5,6.85817491,8,8,.5,6,36.814537905])*scale
        c[11:15]=[.025,-.008,-.015,.015]
        ratio=35/e.paraxial(c,t)['efl'];c/=ratio;t*=ratio
        t[-1]=e.paraxial(c,t)['focus']
        c=np.r_[c,np.zeros(9)]
        x=np.r_[e.pack(c,t),c[21:24]/[1e-5,1e-7,1e-9]]
    else:
        s=json.loads((ROOT/(seed+'_schema.json')).read_text(encoding='utf8'))
        pairs=[(g['catalog'],g['name']) for g in s['glasses']];load_glasses(pairs)
        c=np.r_[s['curvatures'],np.zeros(9)];t=np.array(s['thicknesses'])
        for a in s.get('aspheres',[]):
            j={1:15,15:18,10:21}[a['surface']];c[j:j+3]=a['coefficients']
        if seed=='candidate_v04':
            # Move the iris within the existing air space without moving the
            # glass surfaces. Repair small material/air clearances separately.
            m=e.metrics(c,t)
            clear,_=e.physical_clearances(c,t,np.array(m['ray_semidiameters']))
            shift=max(0.,-clear[4]+.10);t[4]+=shift;t[5]-=shift
            for si in range(14):
                if si not in [4,5]:t[si]+=max(0.,-clear[si]+.03)
            t[-1]=e.paraxial(c,t)['focus']
        x=np.r_[e.pack(c,t),c[21:24]/[1e-5,1e-7,1e-9]]
    lo=np.r_[np.full(14,-3.),[1.0,.1,1.,.65,.3,.3,.65,1.,.1,1.,.25,1.,.1,.65,8.],np.full(9,-30.)]
    hi=np.r_[np.full(14,3.),[8.,3.,8.,6.,10.,10.,6.,8.,3.,8.,12.,8.,8.,8.,25.],np.full(9,30.)]
    x=np.clip(x,lo+1e-8,hi-1e-8);tag='gauss_best' if seed=='gauss' else seed+'_best'
    best=np.inf;calls=0
    def fun(xx):
        nonlocal best,calls
        rr=residual(xx);score=rr@rr;calls+=1
        if score<best:best=score;np.save(ROOT/(tag+'_best.npy'),xx)
        if calls%2000==0:
            cc,tt=unpack(xx);m=e.metrics(cc,tt)
            print(json.dumps(dict(calls=calls,score=best,efl=m['paraxial']['efl'],ttl=sum(tt),rms=[f['rms_um'] for f in m['fields']])),flush=True)
        return rr
    print('INITIAL',json.dumps(e.metrics(*unpack(x))),flush=True)
    opt=least_squares(fun,x,jac=lambda xx:approx_derivative(fun,xx,method='3-point',abs_step=1e-6*(1+abs(xx)),bounds=(lo,hi)),bounds=(lo,hi),x_scale='jac',max_nfev=steps,ftol=1e-9,xtol=1e-10,gtol=1e-6)
    c,t=unpack(opt.x);np.save(ROOT/(tag+'_parameters.npy'),opt.x)
    schema=dict(wd=500.,curvatures=c[:15].tolist(),thicknesses=t.tolist(),media=e.MEDIA.tolist(),stop_index=5,fields=e.FIELDS.tolist(),glasses=e.GLASSES,
        aspheres=[dict(surface=surf,coefficients=c[j:j+3].tolist()) for surf,j in [(1,15),(15,18),(10,21)]])
    (ROOT/(tag+'_schema.json')).write_text(json.dumps(schema,ensure_ascii=False,indent=2),encoding='utf8')
    print('DONE',opt.status,opt.cost,opt.nfev,json.dumps(e.metrics(c,t)),flush=True)
    import verify_design as v
    v.N=e.N;v.GLASSES=e.GLASSES
    audit=v.diffraction(c,t,fields=e.FIELDS,grid=192)
    (ROOT/(tag+'_independent_mtf.json')).write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    print('MTF',json.dumps([{k:d[k] for k in ['image_height_mm','MTF_225_tangential','MTF_180_tangential','MTF_180_sagittal']} for d in audit['fields']]),flush=True)

if __name__=='__main__':main(*(sys.argv[1:] or ['gauss']),steps=300)
