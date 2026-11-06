#!/usr/bin/env python3
"""Independent finite-operator, continuum-response and evolution checks."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, expm_frechet
from scipy.special import roots_legendre

from unitary import Configuration, State, UnitaryMass, UnitarySlab, adjoint
from paired import PairedSlab
from regulator import continuum_periodic, response_parts
from slab import Slab


def continuum_cylinder(cfg,n,transverse,a,harmonic,frequency_squared,paired=False):
    """Closed radial integral for a sharp cylinder, with the contact vertex paired."""
    x,w=roots_legendre(64)
    x=(1+x)/2
    q=2*np.pi/(a*cfg.length)*np.arange(-n,n+1)
    k=2*np.pi*harmonic/(a*cfg.length)
    base=cfg.mass**2+q*q
    p2=(transverse/a)**2
    shifted=base[None,:-harmonic]+x[:,None]*(2*k*q[:-harmonic]+k*k) \
        -(x*(1-x))[:,None]*frequency_squared

    def radial(d):
        return 1/np.sqrt(d)-1/np.sqrt(d+p2)

    pair=np.sum(w[:,None]*radial(shifted))/2
    local=np.sum(radial(base[:-harmonic])) if paired else np.sum(radial(base))
    return float(cfg.mass**2*(pair-local)/(8*np.pi*a*cfg.length*cfg.phase_scale**2))


def spectral(cfg,a,harmonic,frequency_squared):
    box=UnitaryMass(cfg,None)
    wave=np.cos(2*np.pi*harmonic*box.x/cfg.length)/cfg.phase_scale
    T=box.matrix(wave)
    first=np.kron(box.chiral(1j*cfg.mass*T),np.eye(2))
    second=np.kron(box.chiral(-cfg.mass*T @ T),np.eye(2))
    e,v=np.linalg.eigh(box.full_hamiltonian(a,0*wave))
    occupied,empty=v[:,:,:box.d],v[:,:,box.d:]
    contact=np.trace(adjoint(occupied) @ second @ occupied,axis1=-2,axis2=-1).real
    transition=abs(adjoint(empty) @ first @ occupied)**2
    gap=e[:,box.d:,None]-e[:,None,:box.d]
    virtual=np.sum(2*transition*gap/(gap**2-frequency_squared),axis=(1,2))
    answer=box.weight @ (contact-virtual)/a**3
    return float(answer)


def static():
    cfg=Configuration()
    mass=np.array([cfg.mass])
    coefficient=np.ones(1)
    checks=[]
    response=[]
    for a in (1.,1.2):
        for r in (1,2,3):
            for freq2 in (-2.6**2,0.,0.6**2):
                edge,interior,local=response_parts(cfg,mass,coefficient,a,r,frequency_squared=freq2)
                s=(r/a)**2-freq2
                direct=spectral(cfg,a,r,freq2)
                want=s*interior
                checks.append(dict(a=a,harmonic=r,frequency_squared=freq2,
                                   paired_kernel=float(want),matrix_kernel=direct,
                                   old_contact_excess=float(edge),absolute_error=float(abs(direct-want))))
                reference=continuum_periodic(cfg,mass,coefficient,a,r,modes=12000,nodes=80,
                                             frequency_squared=freq2)
                exact=continuum_periodic(cfg,mass,coefficient,a,r,modes=24000,nodes=96,
                                         frequency_squared=freq2)
                ladder=[]
                for n in (12,24,48,96,192,384,768,1536):
                    finite=continuum_cylinder(cfg,n,4*n/3,a,r,freq2)
                    paired=continuum_cylinder(cfg,n,4*n/3,a,r,freq2,paired=True)
                    ladder.append(dict(longitudinal=n,transverse=4*n/3,matched=finite,
                                       relative_error=abs(finite/exact-1),
                                       paired=paired,paired_relative_error=abs(paired/exact-1)))
                error0,error1=(x['relative_error'] for x in ladder[-2:])
                paired0,paired1=(x['paired_relative_error'] for x in ladder[-2:])
                response.append(dict(a=a,harmonic=r,frequency_squared=freq2,unbounded=exact,
                                     reference_refinement_error=abs(reference/exact-1),
                                     ladder=ladder,asymptotic_order=math.log2(error0/error1),
                                     paired_asymptotic_order=math.log2(paired0/paired1)))
    nonlinear=[]
    for n in (6,12,18,24,32):
        config=replace(cfg,longitudinal=n,transverse_cut=4*n/3,points=4*n+17)
        box=UnitarySlab(config,12)
        x=box.spatial.x
        phi=0.55*np.sin(x)+0.25*np.cos(2*x)
        h=box.spatial.full_hamiltonian(1.,phi/cfg.phase_scale)
        h0=box.spatial.full_hamiltonian(1.,0*phi)
        energy=np.linalg.eigvalsh(h)[:,:box.spatial.d]-np.linalg.eigvalsh(h0)[:,:box.spatial.d]
        rho=0.5*box.weight @ energy.sum(axis=1)
        A,*_=box.matching(1.)
        local=(1-A)*np.mean(box.spatial.derivative(phi)**2)/2
        paired=PairedSlab(config,12)
        coeff,*_=paired.paired_coefficients(1.)
        grad=box.spatial.derivative(phi)
        grad_modes=np.fft.fft(grad)/len(grad)
        correction=0.5*np.sum((1-coeff)*abs(grad_modes)**2)
        nonlinear.append(dict(longitudinal=n,raw=float(rho),matched=float(rho-local),
                              paired=float(rho-correction)))
    return dict(contact_checks=checks,periodic_response=response,finite_amplitude=nonlinear)


def small(paired=False,**changes):
    cfg=Configuration(longitudinal=2,transverse_cut=3.,transverse_bins=3,points=17,scalar_modes=2)
    return (PairedSlab if paired else UnitarySlab)(replace(cfg,**changes),None)


def vacuum_state(box,phi,v):
    sector=box.spatial
    _,basis=np.linalg.eigh(sector.full_hamiltonian(1.,phi/box.cfg.phase_scale))
    occupied=basis[:,:,:sector.d].reshape(len(sector.perp),sector.d,2,sector.d)
    A,*_=box.matching(1.)
    return State(1.,phi,A*v,(occupied,))


def homogeneous():
    box=small()
    x=box.spatial.x
    reference=Slab(box.cfg)
    phi,v=np.full_like(x,0.23),np.full_like(x,0.85)
    baseline=reference.initial(phi,v,adiabatic=False)
    state=State(1.,phi,baseline.pi,(baseline.occupied.copy(),))
    obs=box.observe(state)
    for _ in range(160):
        baseline,_=reference.step(baseline,1/160)
        state,obs,_=box.step(state,obs,1/160)
    target=reference.measure(baseline,work=True)
    measured=box.measure(state,obs,work=True)
    l=reference.blocks(baseline.occupied)
    r=box.spatial.blocks(state.orbitals[0])
    return dict(covariance_max_abs_error=float(np.max(abs(l @ adjoint(l)-r @ adjoint(r)))),
                field_max_abs_error=float(np.max(abs(state.phi-baseline.phi))),
                observable_max_abs_error=float(max(abs(measured[k]-target[k]) for k in
                                                   ('rho','P','rho_F','P_F','Jv','Jv_matched'))))


def spatial_ode():
    box=small(expand=False)
    sector=box.spatial
    x=sector.x
    state=vacuum_state(box,0.18*np.sin(x)+0.07*np.cos(2*x),0.47*np.cos(x))
    n=len(x)
    d=sector.d
    cfg=box.cfg
    A,*_=box.matching(1.)
    free=sector.full_hamiltonian(1.,0*x)
    basis=sector.plane

    def pack(s):
        return np.r_[s.phi,s.pi,s.orbitals[0].ravel()]

    y0=pack(state)

    def rhs(_,y):
        phi,pi=y[:n].real,y[n:2*n].real
        o=y[2*n:].reshape(state.orbitals[0].shape)
        T=adjoint(basis) @ ((phi/cfg.phase_scale)[:,None]*basis)
        h=free+np.kron(sector.chiral(cfg.mass*(expm(1j*T)-np.eye(sector.n))),np.eye(2))
        C=np.einsum('lism,ljsm,l->ij',o,o.conj(),sector.weight,optimize=True)
        pullback=expm_frechet(1j*T,C[sector.n:,:sector.n],compute_expm=False)
        force=(1j*np.diag(basis @ pullback @ adjoint(basis))).real*n*cfg.mass/cfg.phase_scale
        laplace=np.fft.ifft(-sector.scalar_k**2*np.fft.fft(phi)).real
        dy=-1j*h @ o.reshape(len(sector.perp),2*d,d)
        return np.r_[pi/A,A*laplace-cfg.scalar_mass**2*phi-sector.project_scalar(force),dy.ravel()]

    duration=0.5
    solution=solve_ivp(rhs,(0,duration),y0,method='DOP853',rtol=2e-13,atol=1e-14)
    if not solution.success:
        raise RuntimeError(solution.message)
    comparisons=[]
    for steps in (40,80,160):
        current=state
        for _ in range(steps):
            current=box.matter(current,duration/steps)
        comparisons.append(dict(steps=steps,max_abs_error=float(np.max(abs(pack(current)-solution.y[:,-1])))))
    for _ in range(steps):
        current=box.matter(current,-duration/steps)
    inverse=float(np.max(abs(pack(current)-y0)))
    return dict(comparisons=comparisons,reversibility_max_abs_error=inverse)


def action(paired=False):
    box=small(paired=paired)
    x=box.spatial.x
    state=replace(box.initial(0.11*np.sin(x)+0.07*np.cos(2*x),0.42*np.cos(x)),a=1.11)
    obs=box.observe(state)
    eps=1e-4
    direction=0.3*np.cos(x)+0.2*np.sin(2*x)
    plus,minus=(replace(state,phi=state.phi+s*eps*direction) for s in (1,-1))
    ep,em=(box.measure(s,box.observe(s))['rho_F'] for s in (plus,minus))
    force=abs(state.a**3*(ep-em)/(2*eps)-np.mean(obs.force*direction))
    d=box.measure(state,obs)
    aa=[replace(state,a=state.a*math.exp(s*eps)) for s in (1,-1)]
    energies=[s.a**3*box.measure(s,obs)['rho'] for s in aa]
    pressure=abs((energies[0]-energies[1])/(2*eps)/(-3*state.a**3*d['P'])-1)
    phase=state.phi/box.cfg.phase_scale
    sector=box.spatial
    new=sector.residual(phase)+sector.chiral(box.cfg.mass*np.eye(sector.n))
    target=box.cfg.mass**2*np.eye(sector.d)
    norm=float(np.max(abs(new @ new-target)))
    return dict(force_abs_error=float(force),pressure_relative_error=float(pressure),mass_squared_error=norm)


def vacuum():
    data={}
    for label,mass,expand in (('static',0.7,False),('massless',0.,True)):
        box=small(mass=mass,expand=expand)
        x=box.spatial.x
        state=box.initial(0.12*np.sin(x)*expand,0.3*np.cos(x)*expand)
        obs=box.observe(state)
        for _ in range(100):
            state,obs,_=box.step(state,obs,.01)
        d=box.measure(state,obs,work=True)
        data[label]={k:float(d[k]) for k in ('rho_F','P_F','Jv','Jv_matched')}
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    result=dict(static=static(),homogeneous=homogeneous(),spatial_ode=spatial_ode(),
                action=action(),paired_action=action(paired=True),vacuum=vacuum())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Unitary source controls complete.',flush=True)


if __name__=='__main__':
    main()
