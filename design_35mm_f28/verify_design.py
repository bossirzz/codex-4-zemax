"""Scalar pupil diffraction estimates and independent design audit.

Angular pupil coordinates retain image-space NA and pupil clipping. A common
image reference retains lateral colour when averaging complex OTFs. No detector
MTF is applied. Native OpticStudio verification remains a separate requirement.
"""
import json
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull
from optical_engine import *

def autocorrelation(pupil):
    size=pupil.shape[0]*2
    padded=np.zeros((size,size),complex)
    padded[:pupil.shape[0],:pupil.shape[1]]=pupil
    corr=np.fft.fftshift(np.fft.ifft2(abs(np.fft.fft2(padded))**2))
    return corr/corr[size//2,size//2]

def ideal_check():
    size=512;u=np.linspace(-1.05,1.05,size);xx,yy=np.meshgrid(u,u)
    otf=autocorrelation((xx*xx+yy*yy<=1).astype(complex))
    shifts=np.arange(-size,size)*(u[1]-u[0]);freq=np.linspace(0,.98,100)
    computed=np.interp(2*freq,shifts,otf[size].real)
    expected=2/np.pi*(np.arccos(freq)-freq*np.sqrt(1-freq*freq))
    return float(np.max(abs(computed-expected)))

def diffraction(c,t,wd=500.,waves=WAVES,fields=FIELDS,grid=192,fixed_stop_radius=None,caps=None):
    p=paraxial(c,t,wd)
    radius=p['efl']/5.6
    rr,_,ok=ray(c,t,N[1],wd,0.,0.,radius,p['ep'])
    stop_radius=abs(rr[9]) if fixed_stop_radius is None else fixed_stop_radius
    # A wider input window includes pupil aberration and field-dependent shape.
    sample=np.linspace(-1.32,1.32,129);xp,yp=np.meshgrid(sample,sample)
    inputs=np.c_[xp.ravel(),yp.ravel()]
    frequencies=np.arange(0,301,1.)
    result=dict(wd_mm=wd,waves_um=np.asarray(waves).tolist(),grid=grid,
        stop_semidiameter_mm=float(stop_radius),method='Scalar angular-pupil OTF; equal spectral weights; independent estimate, pending native OpticStudio verification',fields=[])
    angular_areas=[]
    for field in fields:
        bundles=[];reference=[]
        for n in indices(np.asarray(waves)):
            rays,heights,valid=trace_bundle(c,t,n,wd,float(field/p['magnification']),p['ep'],inputs,radius)
            valid &= rays[:,8]**2+rays[:,9]**2<=stop_radius**2
            valid &= np.all(heights<=(11.35 if caps is None else caps),axis=1)
            rv=rays[valid];hv=heights[valid]
            if len(rv)<100:raise RuntimeError('Insufficient valid pupil samples')
            bundles.append(rv)
            reference.append(np.mean(rv[:,:2],axis=0))
        ref=np.mean(reference,axis=0)
        curves=[];wfes=[];areas=[];maxr=np.zeros(15)
        for wave,rv in zip(waves,bundles):
            vz=rv[:,11];vx=rv[:,2]*vz;vy=rv[:,3]*vz
            image_z=float(np.sum(t))
            opl_last=rv[:,4]-(image_z-rv[:,7])/vz
            # Eikonal of each exiting plane wave, evaluated at common image Q.
            phase=opl_last+vx*(ref[0]-rv[:,5])+vy*(ref[1]-rv[:,6])+vz*(image_z-rv[:,7])
            phase-=np.mean(phase)
            xy=np.c_[vx,vy]
            ctr=(xy.max(0)+xy.min(0))/2
            extent=float(np.max(xy.max(0)-xy.min(0))/2*1.03)
            u=np.linspace(-extent,extent,grid)
            ux,uy=np.meshgrid(u+ctr[0],u+ctr[1])
            wavefront=LinearNDInterpolator(xy,phase,fill_value=np.nan)(ux,uy)
            mask=np.isfinite(wavefront)
            pupil=np.zeros((grid,grid),complex)
            pupil[mask]=np.exp(2j*np.pi*wavefront[mask]/(wave*.001))
            otf=autocorrelation(pupil)
            shifts=np.arange(-grid,grid)*(u[1]-u[0])
            wave_mm=wave*.001
            # Sagittal = x, tangential = y for a y-field rotational system.
            sag=np.interp(frequencies*wave_mm,shifts,otf[grid,:].real)+1j*np.interp(frequencies*wave_mm,shifts,otf[grid,:].imag)
            tan=np.interp(frequencies*wave_mm,shifts,otf[:,grid].real)+1j*np.interp(frequencies*wave_mm,shifts,otf[:,grid].imag)
            curves.append((sag,tan))
            # Remove a fitted tilt for reporting WFE, not for polychromatic OTF.
            fit=np.linalg.lstsq(np.c_[np.ones(len(xy)),xy],phase,rcond=None)[0]
            wfes.append(float(np.std(phase-np.c_[np.ones(len(xy)),xy]@fit)/(wave*.001)))
            areas.append(float(np.sum(mask)*(u[1]-u[0])**2))
        sag=abs(np.mean([cr[0] for cr in curves],axis=0));tan=abs(np.mean([cr[1] for cr in curves],axis=0))
        angular_areas.append(np.mean(areas))
        result['fields'].append(dict(image_height_mm=float(field),reference_image=ref.tolist(),wfe_rms_waves=wfes,
            MTF_225_sagittal=float(sag[225]),MTF_225_tangential=float(tan[225]),
            MTF_180_sagittal=float(sag[180]),MTF_180_tangential=float(tan[180]),
            sagittal=sag.tolist(),tangential=tan.tolist()))
    for i,f in enumerate(result['fields']):f['relative_illumination_estimate']=float(angular_areas[i]/angular_areas[0])
    return result

def audit(filename='optimized_parameters.npy',grid=192):
    c,t=unpack(np.load(ROOT/filename))
    result=dict(ideal_circular_pupil_MTF_max_error=ideal_check(),design=metrics(c,t),visible=diffraction(c,t,grid=grid))
    ROOT.joinpath('independent_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(ideal_check=result['ideal_circular_pupil_MTF_max_error'],fields=[{k:v for k,v in f.items() if k not in ['sagittal','tangential']} for f in result['visible']['fields']]),ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':audit()
