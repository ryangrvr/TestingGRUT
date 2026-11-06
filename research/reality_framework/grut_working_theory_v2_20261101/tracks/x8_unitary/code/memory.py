#!/usr/bin/env python3
"""Leading two-time vacuum emission for archived modes on a prescribed H55 trace."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

from transport import initial_data, read_columns


def radiated(eta,amplitude,k,mass,width,window,pad):
    """Positive-frequency massless two-body spectrum for each chiral component."""
    dt=eta[1]-eta[0]
    n=pad*len(eta)
    freq=2*np.pi*np.fft.rfftfreq(n,dt)
    out=np.zeros((len(k),2,2))
    truncated=np.zeros_like(out)
    for j,sign in enumerate((1,-1)):
        signal=window[:,None,None]*amplitude[:,:,j,:]
        freq_amp=np.fft.ifft(signal,n=n,axis=0)[:len(freq)]*(n*dt)
        density=np.sum(abs(freq_amp)**2,axis=2)*(width/(2*np.pi*mass))
        for i,q in enumerate(abs(k)):
            ix=int(np.searchsorted(freq,q,side='right'))
            spec=np.r_[np.interp(q,freq,density[:,i]),density[ix:,i]]
            omega=np.r_[q,freq[ix:]]
            p=(omega-sign*k[i])*spec
            for a,z in enumerate((p,p*omega)):
                out[i,j,a]=trapezoid(z,omega)
                keep=omega<freq[-1]/2
                truncated[i,j,a]=trapezoid(z[keep],omega[keep])
    return out,truncated


def local(eta,k,mass,width,a,delta,omega,occ,window):
    out=np.zeros((len(k),2,2))
    for j,sign in enumerate((1,-1)):
        frequency=omega-sign*delta
        rate=width/mass*(a[:,None]*mass)**2/(2*omega)
        weight=rate*occ*window[:,None]**2*(frequency>abs(k))
        out[:,j,0]=trapezoid(weight,eta,axis=0)
        out[:,j,1]=trapezoid(weight*frequency,eta,axis=0)
    return out


def stationary():
    n=4096
    eta=np.arange(n)*100/n
    k=np.array([1.2,-1.2])
    mass=np.full(2,2.)
    width=np.full(2,0.001)
    d=0.4
    omega=np.hypot(mass,k+d)
    z=(k+d)/omega
    positive=np.stack([np.sqrt((1+z)/2),np.sqrt((1-z)/2)],axis=1)
    negative=np.stack([-positive[:,1],positive[:,0]],axis=1)
    phase=np.exp(1j*d*eta[:,None,None]*np.array([1.,-1.]))
    g=np.sin(np.pi*eta/100)**4
    good=phase*positive[None,:,:]*np.exp(-1j*eta[:,None,None]*omega[:,None])
    bad=phase*negative[None,:,:]*np.exp(1j*eta[:,None,None]*omega[:,None])
    measured,_=radiated(eta,good[...,None],k,mass,width,g,16)
    vacuum,_=radiated(eta,bad[...,None],k,mass,width,g,16)
    expected=local(eta,k,mass,width,np.ones(n),np.full((n,1),d),
                   np.broadcast_to(omega,(n,2)),np.ones((n,2)),g)
    sums=expected.sum(axis=1)
    return dict(probability_relative_error=float(np.max(abs(measured.sum(axis=1)[:,0]/sums[:,0]-1))),
                energy_relative_error=float(np.max(abs(measured.sum(axis=1)[:,1]/sums[:,1]-1))),
                negative_branch_probability=float(np.max(vacuum.sum(axis=1)[:,0])))


def select(model,terminal,state,widths):
    result=[]
    q=state['q']
    a=terminal['a']
    for i,(parent,mass,f,_) in enumerate(model['species']):
        k=np.array([-1,1])[:,None]*q/a
        E=np.hypot(k,mass)
        r=state['bloch'][i]
        occupation=(1+(mass*r[...,0]+k*r[...,2])/E)/2
        cumulative=np.cumsum(state['w']*np.sum(E*occupation,axis=0))
        for percentile in (0.25,0.75):
            j=int(np.searchsorted(cumulative,percentile*cumulative[-1]))
            for h in (-1,1):
                result.append(dict(species=i+1,parent=parent,mass=mass,f=f,width=widths[i],
                                   percentile=percentile,helicity=h,q=float(q[j]),
                                   bloch=r[(h+1)//2,j].tolist()))
    return result


def replay(history,modes,n):
    t=history['t']
    a=history['a']
    a_of_t=CubicHermiteSpline(t,a,history['H']*a)
    eta_map=np.r_[0.,np.cumsum((1/a[:-1]+1/a[1:])*np.diff(t)/2)]
    t_of_eta=PchipInterpolator(eta_map,t)
    eta=np.linspace(0.,eta_map[-1],n,endpoint=False)
    time=t_of_eta(eta)
    midtime=t_of_eta(eta+0.5*(eta[1]-eta[0]))
    scale=a_of_t(time)
    scale_mid=a_of_t(midtime)
    parent=np.array([r['parent'] for r in modes])
    f=np.array([r['f'] for r in modes])
    m=np.array([r['mass'] for r in modes])
    width=np.array([r['width'] for r in modes])
    k=np.array([r['q']*r['helicity'] for r in modes])
    x,y,z=np.array([r['bloch'] for r in modes]).T
    cov=np.zeros((len(modes),2,2),complex)
    cov[:,0,0]=(1+z)/2
    cov[:,1,1]=(1-z)/2
    cov[:,0,1]=(x-1j*y)/2
    cov[:,1,0]=(x+1j*y)/2
    vals,vecs=np.linalg.eigh(cov)
    correction=float(np.max(abs(vals-np.clip(vals,0,1))))
    if correction>1e-10:
        raise ValueError('Archive exceeds its Bloch ball tolerance')
    psi=np.empty((n,*cov.shape),complex)
    psi[0]=vecs*np.sqrt(np.clip(vals,0,1))[:,None,:]

    def scalar(at):
        phi=np.empty((len(at),2))
        vel=np.empty_like(phi)
        for i,(field,velocity) in enumerate((('sigma','dsigma'),('chi','dchi'))):
            spline=CubicHermiteSpline(t,history[field],history[velocity])
            phi[:,i]=spline(at)
            vel[:,i]=spline(at,1)
        return phi[:,parent]/f,vel[:,parent]/f

    theta,theta_dot=scalar(time)
    _,mid_dot=scalar(midtime)
    X=scale_mid[:,None]*m
    Z=k+scale_mid[:,None]*mid_dot/2
    W=np.hypot(X,Z)
    dt=eta[1]-eta[0]
    cc=np.cos(W*dt)[...,None]
    sx=(-1j*X*np.sin(W*dt)/W)[...,None]
    sz=(-1j*Z*np.sin(W*dt)/W)[...,None]
    for i in range(n-1):
        left,right=psi[i,:,0,:],psi[i,:,1,:]
        psi[i+1,:,0,:]=(cc[i]+sz[i])*left+sx[i]*right
        psi[i+1,:,1,:]=sx[i]*left+(cc[i]-sz[i])*right
    delta=scale[:,None]*theta_dot/2
    X=scale[:,None]*m
    Z=k+delta
    omega=np.hypot(X,Z)
    spin=np.stack([np.sqrt((1+Z/omega)/2),np.sqrt((1-Z/omega)/2)],axis=2)
    occupied=np.sum(abs(np.einsum('tmj,tmja->tma',spin,psi))**2,axis=2)
    phase=np.exp(0.5j*theta[:,:,None]*np.array([1.,-1.]))
    amplitude=phase[:,:,:,None]*psi
    return eta,time,scale,k,m,width,delta,omega,occupied,amplitude,psi,theta,correction


def frozen(eta,i,a,k,m,w,d,omega,psi,theta,window,pad):
    """Same exposure and central state, with the Hamiltonian and phase speed frozen."""
    center=psi[i]
    freq=omega[i]
    p=np.stack([np.sqrt((1+(k+d[i])/freq)/2),np.sqrt((1-(k+d[i])/freq)/2)],axis=1)
    projected=np.einsum('mi,mia->ma',p,center)
    positive=p[:,:,None]*projected[:,None,:]
    negative=center-positive
    tau=eta-eta[i]
    carrier=np.exp(-1j*tau[:,None]*freq)[:,:,None,None]
    phase=np.exp(0.5j*(theta[i]+2*tau[:,None]*d[i])[:,:,None]*np.array([1.,-1.]))
    amp=phase[...,None]*(carrier*positive+carrier.conj()*negative)
    two_time,_=radiated(eta,amp,k,m,w,window,pad)
    occ=np.broadcast_to(np.sum(abs(projected)**2,axis=1),(len(eta),len(k)))
    ref=local(eta,k,m,w,np.full(len(eta),a[i]),d[i],freq,occ,window)
    return two_time.sum(axis=1)[:,1],ref.sum(axis=1)[:,1]


def evaluate(history,modes,n,pad,include_controls=False):
    eta,t,a,k,m,w,d,omega,occ,amp,psi,theta,correction=replay(history,modes,n)
    offsets=(0.02,0.065)
    result=[]
    for start in offsets:
        duration=0.04
        x=(t-t[0]-start)/duration
        envelope=np.sin(np.pi*np.clip(x,0,1))**4
        two_time,cut=radiated(eta,amp,k,m,w,envelope,pad)
        markov=local(eta,k,m,w,a,d,omega,occ,envelope)
        if include_controls:
            center=int(np.argmin(abs(t-(t[0]+start+duration/2))))
            frozen_energy,frozen_local=frozen(eta,center,a,k,m,w,d,omega,psi,theta,envelope,pad)
            alternative=np.sin(np.pi*np.clip(x,0,1))**2
            switched,_=radiated(eta,amp,k,m,w,alternative,pad)
            switched_local=local(eta,k,m,w,a,d,omega,occ,alternative)
        rows=[]
        for i,info in enumerate(modes):
            e=float(np.sum(two_time[i,:,1]))
            p=float(np.sum(two_time[i,:,0]))
            ref=float(np.sum(markov[i,:,1]))
            rows.append(dict(mode=i,exact_probability=p,exact_energy=e,local_energy=ref,
                             exact_to_local_energy=e/ref,
                             uv_relative_error=abs(float(np.sum(cut[i,:,1]))/e-1),
                             exact_energy_by_branch=two_time[i,:,1].tolist(),
                             local_energy_by_branch=markov[i,:,1].tolist()))
            if include_controls:
                rows[-1].update(frozen_energy=float(frozen_energy[i]),frozen_local_energy=float(frozen_local[i]),
                                alternative_envelope_energy_ratio=float(np.sum(switched[i,:,1])/
                                                                        np.sum(switched_local[i,:,1])))
        result.append(dict(start=start,duration=duration,modes=rows))
    return result,correction


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path,help='X7 fine coherent history')
    p.add_argument('output',type=Path)
    p.add_argument('--coarse-history',type=Path)
    args=p.parse_args()
    history=read_columns(args.history)
    model,terminal,_,state,sha,widths,_=initial_data()
    if len(history['t'])!=916 or abs(history['t'][0]-terminal['t'])>1e-12 or max(history['rho_sm'])!=0:
        raise ValueError('Expected the X7 fine coherent history')
    modes=select(model,terminal,state,widths)
    time_coarse,_=evaluate(history,modes,8192,32)
    frequency_coarse,_=evaluate(history,modes,16384,16)
    fine,change=evaluate(history,modes,16384,32,include_controls=True)

    def compare(control,key):
        return max(abs(a[key]/b[key]-1) for x,y in zip(control,fine,strict=True)
                   for a,b in zip(x['modes'],y['modes'],strict=True))

    rows=[row for window in fine for row in window['modes']]
    background=None
    if args.coarse_history:
        lower,_=evaluate(read_columns(args.coarse_history),modes,16384,32)
        background=compare(lower,'exact_energy')
    result=dict(
        scope='two-time order-Yukawa-squared vacuum emission with smooth switching on a prescribed homogeneous trace',
        checkpoint_sha256=sha,history_sha256=hashlib.sha256(args.history.read_bytes()).hexdigest(),
        stationary=stationary(),mode_selection=modes,windows=fine,
        time_coarse_resolution=dict(conformal_samples=8192,fft_padding=32),
        frequency_coarse_resolution=dict(conformal_samples=16384,fft_padding=16),
        fine_resolution=dict(conformal_samples=16384,fft_padding=32),
        conformal_timestep_refinement_max_relative_error=compare(time_coarse,'exact_energy'),
        spectral_spacing_refinement_max_relative_error=compare(frequency_coarse,'exact_energy'),
        background_timestep_refinement_max_relative_error=background,
        local_quadrature_refinement_max_relative_error=compare(time_coarse,'local_energy'),
        uv_refinement_max_relative_error=max(row['uv_relative_error'] for row in rows),
        archive_positivity_roundoff_adjustment=change,
        maximum_exact_probability=max(row['exact_probability'] for row in rows),
        memory_to_local_energy_ratio_range=[min(row['exact_to_local_energy'] for row in rows),
                                            max(row['exact_to_local_energy'] for row in rows)],
        maximum_selected_sample_energy_difference=max(abs(row['exact_to_local_energy']-1) for row in rows),
        samples_differing_by_over_10_percent=sum(abs(row['exact_to_local_energy']-1)>0.1 for row in rows),
        scientific_acceptance=False,
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Two-time versus local energy ratio',result['memory_to_local_energy_ratio_range'],flush=True)


if __name__=='__main__':
    main()
