#!/usr/bin/env python3
"""Independent helicity and spatial ODEs, and finite-action matching controls."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from slab import Configuration, SIGMA_X, SIGMA_Z, Slab, State, adjoint


def small_box(**kwargs):
    cfg = Configuration(longitudinal=2,transverse_cut=3.0,transverse_bins=3,points=25,scalar_modes=2)
    return Slab(replace(cfg,**kwargs))


def homogeneous():
    box = small_box()
    cfg = box.cfg
    q = np.sqrt(box.q2)[...,None]*np.array([-1,1])
    w = np.broadcast_to(box.weight[:,None,None],q.shape)
    M, f = cfg.mass,cfg.phase_scale
    initial_phi, initial_speed = 0.23,0.85
    initial = box.initial(np.full(cfg.points,initial_phi),np.full(cfg.points,initial_speed),adiabatic=False)
    duration = 1.0

    def kinetics(a):
        k = q/a
        E = np.hypot(k,M)
        zt = M*M*np.sum(w/E**3)/(8*f*f*a**3)
        zg = M*M*np.sum(w*k*k/E**5)/(8*a**3)
        return k,E,1-zt,zg

    def axes(a, phi):
        k,E,A,zg = kinetics(a)
        b = np.stack([np.full_like(k,M*math.cos(phi/f)),
                      np.full_like(k,-M*math.sin(phi/f)),k],axis=-1)
        return b,E,A,zg

    def observe(y):
        a,phi,pi = y[:3]
        r = y[3:].reshape(*q.shape,3)
        b,E,A,zg = axes(a,phi)
        v = pi/(a**3*A)
        occ = (1+np.sum(b*r,axis=-1)/E)/2
        rho_f = float(np.sum(w*E*occ)/a**3)
        pressure = float(np.sum(w*b[...,2]*(r[...,2]+b[...,2]/E))/(6*a**3))
        rho = cfg.vacuum_density+0.5*A*v*v+0.5*cfg.scalar_mass**2*phi*phi+rho_f
        H = math.sqrt(rho/(3*cfg.planck**2+zg/2))
        return a,phi,r,b,v,occ,rho_f,pressure,H

    def rhs(_,y):
        a,phi,r,b,v,_,_,_,H = observe(y)
        J = -M*np.sum(w*(math.sin(phi/f)*r[...,0]+math.cos(phi/f)*r[...,1]))/(2*f)
        return np.r_[a*H,v,-a**3*cfg.scalar_mass**2*phi-J,(2*np.cross(b,r)).ravel()]

    b0,E0,A0,_ = axes(1.0,initial_phi)
    y0 = np.r_[1.0,initial_phi,A0*initial_speed,(-b0/E0[...,None]).ravel()]
    ode = solve_ivp(rhs,(0,duration),y0,rtol=3e-13,atol=1e-14,method='DOP853')
    assert ode.success
    a,phi,_,b,v,f_ref,rho_f,p_ref,_ = observe(ode.y[:,-1])
    comparisons = []
    for n in (100,200,400):
        state=initial
        for _ in range(n):
            state,_ = box.step(state,duration/n)
        moments=box.measure(state,work=True)
        speed = np.mean(state.pi)/(state.a**3*moments['A'])
        blocks = box.blocks(state.occupied)
        covariance = blocks @ adjoint(blocks)
        f_slab = np.zeros_like(q)
        for ell,perp in enumerate(box.perp):
            for j,k in enumerate(box.k):
                _, helicities = np.linalg.eigh(k*SIGMA_Z+perp*SIGMA_X)
                phi_s = np.mean(state.phi)
                for h in (0,1):
                    projector = np.kron(np.eye(2),helicities[:,h,None])
                    c = adjoint(projector) @ covariance[ell,j] @ projector
                    kq=q[ell,j,h]/state.a
                    phase=phi_s/cfg.phase_scale
                    term = np.array([[kq,M*np.exp(1j*phase)],[M*np.exp(-1j*phase),-kq]])
                    pplus = (np.eye(2)+term/np.hypot(kq,M))/2
                    f_slab[ell,j,h] = np.trace(pplus @ c).real
        relative = dict(a=abs(state.a/a-1), phi=abs(np.mean(state.phi)/phi-1),
                        speed=abs(speed/v-1), rho_F=abs(moments['rho_F']/rho_f-1),
                        pressure_F=abs(moments['P_F']/p_ref-1))
        comparisons.append(dict(steps=n, relative_errors=relative,
                                mode_occupation_max_absolute=float(max(abs(f_slab-f_ref).ravel())),
                                field_spread=float(np.std(state.phi))))
    assert max(comparisons[-1]['relative_errors'].values()) < 5e-4
    assert comparisons[-1]['mode_occupation_max_absolute'] < 1e-5
    assert comparisons[-1]['mode_occupation_max_absolute'] < comparisons[0]['mode_occupation_max_absolute']/10
    assert comparisons[-1]['field_spread'] < 1e-10
    return comparisons


def spatial_ode():
    box = small_box(expand=False)
    cfg = box.cfg
    duration = 0.6
    initial = box.initial(0.1*np.sin(box.x),0.4*np.cos(box.x),adiabatic=False)
    a = initial.a
    points = cfg.points
    A,B,*_ = box.matching(a)
    original = np.r_[initial.phi,initial.pi,initial.occupied.ravel()]

    def rhs(_, y):
        phi,pi = y[:points].real,y[points:2*points].real
        occupied = y[2*points:].reshape(initial.occupied.shape)
        phase = phi/cfg.phase_scale
        laplace = np.fft.ifft(-box.scalar_k**2*np.fft.fft(phi)).real
        gradient = box.project_scalar(box.mass_derivative(phase,box.cov(occupied)))
        h = box.full_hamiltonian(a,phase)
        dpsi = -1j*h @ occupied.reshape(len(box.perp),2*box.d,box.d)
        return np.r_[pi/(a**3*A),a**3*(B*laplace/a**2-cfg.scalar_mass**2*phi)-gradient,dpsi.ravel()]

    ode = solve_ivp(rhs,(0,duration),original,rtol=3e-13,atol=1e-14,method='DOP853')
    assert ode.success
    values = []
    for n in (60,120,240):
        state=initial
        for _ in range(n):
            state=box.matter(state,duration/n)
        final=np.r_[state.phi,state.pi,state.occupied.ravel()]
        values.append(dict(steps=n,max_absolute_error=float(max(abs(final-ode.y[:,-1])))))
    assert values[-1]['max_absolute_error'] < 1e-5
    assert values[-1]['max_absolute_error'] < values[0]['max_absolute_error']/10
    for _ in range(n):
        state=box.matter(state,-duration/n)
    inverse=max(abs(np.r_[state.phi,state.pi,state.occupied.ravel()]-original))
    assert inverse<1e-10
    return dict(ode=values,inversion_max_absolute=float(inverse))


def static_and_massless():
    box = small_box(expand=False)
    zero = np.zeros(box.cfg.points)
    state = box.initial(zero,zero,adiabatic=False)
    for _ in range(100):
        state=box.matter(state,0.01)
    still=box.measure(state,work=True)
    assert abs(still['rho_F'])<1e-11 and abs(still['Jv'])<1e-11
    zero_mass=small_box(mass=0.0)
    state=zero_mass.initial(np.sin(zero_mass.x),np.cos(zero_mass.x))
    for _ in range(100):
        state,_=zero_mass.step(state,0.01)
    decoupled=zero_mass.measure(state,work=True)
    assert abs(decoupled['rho_F'])<1e-11 and decoupled['Jv']==0
    return dict(static_vacuum=still,massless_wave=decoupled)


def action_matching():
    box=small_box()
    cfg=box.cfg
    a=1.19
    wave=np.cos(2*math.pi*box.x/cfg.length)
    eps=0.002

    def vacuum(phi):
        energies=np.linalg.eigvalsh(box.full_hamiltonian(a,phi/cfg.phase_scale))
        return 0.5*(box.weight @ energies[:,:box.d].sum(axis=1))/a**3

    response=(vacuum(eps*wave)+vacuum(-eps*wave)-2*vacuum(0*wave))/eps**2
    zx=box.spatial_match(a)
    expected=zx*np.mean((box.derivative(wave)/a)**2)
    error=abs(response/expected-1)
    assert error<1e-5
    assert abs(vacuum(eps*np.ones(cfg.points))-vacuum(0*wave))<1e-11

    state=box.initial(0.18*np.sin(box.x),0.52*wave,a=a,adiabatic=True)
    s=box.measure(state,work=True)
    delta=1e-5
    u=box.project_scalar(0.3*wave+0.2*np.sin(2*box.x))
    de=(box.measure(replace(state,phi=state.phi+delta*u))['rho_F']
        -box.measure(replace(state,phi=state.phi-delta*u))['rho_F'])*a**3/(2*delta)
    work=np.mean(u*box.mass_derivative(state.phi/cfg.phase_scale,box.cov(state.occupied)))
    force_error=abs(de-work)
    assert force_error<1e-8

    plus=a*math.exp(delta)
    minus=a*math.exp(-delta)
    derivative=(plus**3*box.measure(replace(state,a=plus))['rho']
                -minus**3*box.measure(replace(state,a=minus))['rho'])/(2*delta)
    pressure_error=abs(derivative/(-3*a**3*s['P'])-1)
    assert pressure_error<1e-8
    return dict(spatial_subtraction=zx,static_curvature_relative_error=float(error),
                field_variation_absolute_error=float(force_error),pressure_variation_relative_error=float(pressure_error))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    result=dict(helicity=homogeneous(),spatial=spatial_ode(),vacuum=static_and_massless(),matching=action_matching())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Independent helicity and spatial ODEs, vacuum and matching controls passed.',flush=True)


if __name__=='__main__':
    main()
