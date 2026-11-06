#!/usr/bin/env python3
"""Actual early H55 axial seed and a continuously coupled N1 momentum probe."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp, trapezoid
from scipy.interpolate import CubicHermiteSpline

from coupled import Configuration, CoupledSlab, Geometry, State, adjoint
from h55_seed import archive, axial, reconstruct
from run import read_columns

BEGIN,END=-0.02,0.10


def hessian(P,sigma,chi,s):
    g2=P['g']**2
    gs2=P['g_s']**2
    u=g2*sigma*sigma-P['lam']*P['v']**2
    log=np.log(max(abs(u),1e-50)/P['mu2'])
    loop=((2*g2*sigma)**2*log+2*g2*u*(log-1))/(32*np.pi**2)
    return np.array([
        [g2*chi**2-P['m2']+loop,2*g2*sigma*chi,0.],
        [2*g2*sigma*chi,P['lam']*(3*chi**2-P['v']**2)+g2*sigma**2+gs2*s**2,2*gs2*chi*s],
        [0.,2*gs2*chi*s,gs2*chi**2],
    ])


def make_source(run,fields,history):
    use=(history['t']<=END+0.003)
    data={k:history[k][use] for k in ('t','a','H','sigma','dsigma','chi','dchi','s','ds')}
    t=data['t']
    if abs(t[0]-BEGIN)>1e-12 or t[-1]<END:
        raise ValueError('Need the original X9 early history')
    scale=CubicHermiteSpline(t,data['a'],data['a']*data['H'])
    mean=[CubicHermiteSpline(t,data[k],data[d]) for k,d in
          (('sigma','dsigma'),('chi','dchi'),('s','ds'))]
    k=2*np.pi/run['L']
    keys=('sigma','chi','spectator_a','dsigma','dchi','dspectator_a')
    initial=np.array([np.fft.fft(axial(fields[key])[0])[1]/run['N'] for key in keys])

    def rhs(t,y):
        a=scale(t)
        H=scale(t,1)/a
        V=hessian(run['P'],*(f(t) for f in mean))
        return np.r_[y[3:],-(V+(k/a)**2*np.eye(3)) @ y[:3]-3*H*y[3:]]

    mode=solve_ivp(rhs,(BEGIN,END),initial,method='DOP853',rtol=1e-10,atol=1e-11,
                   max_step=0.001,dense_output=True)
    if not mode.success:
        raise RuntimeError(mode.message)
    sha=hashlib.sha256(np.column_stack(list(data.values())).astype('<f8').tobytes()).hexdigest()
    return scale,mean,mode.sol,sha


def prepare(slab,a,adot,phi,velocity):
    h=slab.full_hamiltonian(a,phi)
    e,u=np.linalg.eigh(h)
    occ=slab.dim//2
    if np.any(e[:,occ-1]>=0) or np.any(e[:,occ]<=0):
        raise ValueError('Ungapped initial H55 probe')
    empty,occupied=u[:,:,occ:],u[:,:,:occ]
    eps=1e-6
    tangent=(slab.full_hamiltonian(a+eps*adot,phi+eps*velocity)-
             slab.full_hamiltonian(a-eps*adot,phi-eps*velocity))/(2*eps)
    gap=e[:,occ:,None]-e[:,None,:occ]
    X=-1j*(adjoint(empty) @ tangent @ occupied)/gap**2
    vals,axes=np.linalg.eigh(np.eye(occ)+adjoint(X) @ X)
    factor=(axes*vals[:,None,:]**-0.5) @ adjoint(axes)
    frame=(occupied+empty @ X) @ factor
    return frame,float(np.max(np.linalg.norm(X,axis=(-1,-2))))


def correlations(slab,frame):
    C=(frame @ adjoint(frame))[0]
    ster=slab.sterile_indices
    bath=np.setdiff1d(np.arange(slab.dim),ster)
    q=np.repeat(np.arange(slab.n),2*slab.block)
    cross=C[ster[:,None],bath[None,:]]
    off=q[ster,None]!=q[None,bath]
    return dict(sterile_bath_norm=float(np.linalg.norm(cross)),
                spatial_sterile_bath_norm=float(np.linalg.norm(cross[off])),
                spatial_covariance_norm=float(np.linalg.norm(C[q[:,None]!=q[None,:]])))


def probe(run,source,dt,spatial):
    scale,means,mode,_=source
    mean=means[0]
    card=run['channels'][0]
    mass,phase,width=card['M'],card['f'],card['Grest']
    cfg=Configuration(longitudinal=1,transverse_bins=1,transverse_cut=mass,points=25,
                      scalar_modes=1,mass=mass,phase_scale=phase,length=run['L'],expand=False)
    slab=CoupledSlab(Geometry(cfg,bath_bins=16,bath_cutoff=3*mass,gamma=width/(2*np.pi*mass)))
    x=2*np.pi*slab.box.x/run['L']
    factor=np.exp(1j*x)

    def values(time):
        perturbation=mode(time)
        f,v=(perturbation[0],perturbation[3]) if spatial else (0.,0.)
        return mean(time)+2*(f*factor).real,mean(time,1)+2*(v*factor).real

    a=float(scale(BEGIN))
    phi,v=values(BEGIN)
    frame,graph=prepare(slab,a,float(scale(BEGIN,1)),phi,v)
    start=correlations(slab,frame)
    steps=round((END-BEGIN)/dt)
    times=np.linspace(BEGIN,END,steps+1)
    dt=(END-BEGIN)/steps
    energy,work,metric=[],[],[]
    max_phase=0.
    for i,t in enumerate(times):
        a=float(scale(t))
        phi,v=values(t)
        obs=slab.observe(State(a,phi,0*phi,frame),source=True)
        z,mx,kin,bath,cross,plus,minus,count,mx_total,_,force=obs
        K=float(np.sum(slab.weight*(kin+bath+cross)))
        m=slab.mass(a)
        energy.append(0.5*(K/a+m*mx_total))
        work.append(m*np.mean(force*v))
        metric.append(float(scale(t,1))/a/2*(-K/a+slab.mass_slope(a)*mx_total))
        max_phase=max(max_phase,float(np.max(abs(phi-mean(t)))/phase))
        if i<steps:
            middle=(times[i]+times[i+1])/2
            mid_a=float(scale(middle))
            mid_phi,_=values(middle)
            frame=slab.quantum_free(mid_a,frame,dt/2)
            frame,_=slab.local(mid_a,mid_phi,frame,dt,False)
            frame=slab.quantum_free(mid_a,frame,dt/2)
    E=np.array(energy)
    scalar=cumulative_trapezoid(work,times,initial=0.)
    redshift=cumulative_trapezoid(metric,times,initial=0.)
    eta=float(trapezoid(1/scale(times),times))
    return dict(dt=dt,spatial=spatial,sterile_mass=mass,phase_scale=phase,rest_width=width,
                transverse_momentum=mass/2,bath_cutoff=3*mass,dimension=slab.dim,
                initial_graph_norm=graph,initial=start,terminal=correlations(slab,frame),
                conformal_duration=eta,minimum_discrete_recurrence=slab.recurrence,
                phase_perturbation_peak=max_phase,integrated_scalar_work=float(scalar[-1]),
                integrated_redshift_work=float(redshift[-1]),
                terminal_raw_quantum_energy=float(E[-1]),
                max_relative_energy_work_error=float(np.max(abs(E-E[0]-scalar-redshift))/abs(E[0])),
                orthonormality_error=float(np.max(abs(adjoint(frame) @ frame-np.eye(frame.shape[-1])))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    run,_,sha=archive()
    source=make_source(run,reconstruct(run),read_columns(args.history))
    jobs={'spatial':(4e-5,True),'spatial_half':(2e-5,True),'homogeneous':(4e-5,False)}
    with ThreadPoolExecutor(max_workers=3) as pool:
        cases=dict(zip(jobs,pool.map(lambda job:probe(run,source,*job),jobs.values()),strict=True))
    result=dict(scope='early N1 spatial probe on the reconstructed axial seed and a prescribed linearized history',
                interval=[BEGIN,END],diagnostic_sha256=sha,source_trace_sha256=source[-1],
                cases=cases,scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Physical H55 probe completed.',flush=True)


if __name__=='__main__':
    main()
