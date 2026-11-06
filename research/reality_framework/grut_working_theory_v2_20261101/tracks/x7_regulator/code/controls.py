#!/usr/bin/env python3
"""Independent signed helicity and action checks for the regulated source."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from evolution import Configuration, MASSES, RegulatedSlab, State
from slab import SIGMA_X, SIGMA_Z, adjoint


def helicity():
    cfg = replace(Configuration(),longitudinal=2,transverse_cut=3.,transverse_bins=3,
                  points=25,scalar_modes=2)
    slab = RegulatedSlab(cfg,MASSES)
    box = slab.spatial
    q = np.sqrt(box.q2)[None,:,:,None]*np.array([-1,1])
    m = slab.masses[:,None,None,None]
    w = slab.signed[:,None,None,None]*box.weight[None,:,None,None]
    phase_scale = cfg.phase_scale

    def axes(a,phi):
        k = np.broadcast_to(q/a,(len(slab.masses),*q.shape[1:]))
        e = np.hypot(m,k)
        b = np.stack([np.broadcast_to(m*math.cos(phi/phase_scale),e.shape),
                      np.broadcast_to(-m*math.sin(phi/phase_scale),e.shape),k],axis=-1)
        z = np.sum(w*m*m/e**3)/(8*a**3*phase_scale**2)
        zg = np.sum(w*m*m*k*k/e**5)/(8*a**3)
        return k,e,b,1-z,zg

    def measure(y):
        a,phi,pi = y[:3]
        k,e,b,A,zg = axes(a,phi)
        r = y[3:].reshape(*e.shape,3)
        f = (1+np.sum(b*r,axis=-1)/e)/2
        v = pi/(a**3*A)
        e_F = np.sum(w*e*f)/a**3
        p_F = np.sum(w*k*(r[...,2]+k/e))/(6*a**3)
        density = cfg.vacuum_density+(A*v*v+cfg.scalar_mass**2*phi*phi)/2+e_F
        H = math.sqrt(density/(3*cfg.planck**2+zg/2))
        return f,e_F,p_F,r,b,v,H

    def rhs(_,y):
        a,phi,_ = y[:3]
        _,_,_,r,b,v,H = measure(y)
        J = -np.sum(w*m*(math.sin(phi/phase_scale)*r[...,0]+math.cos(phi/phase_scale)*r[...,1]))/(2*phase_scale)
        return np.r_[a*H,v,-a**3*cfg.scalar_mass**2*phi-J,(2*np.cross(b,r)).ravel()]

    phi0,vel0 = 0.23,0.85
    _,e,b,A,_ = axes(1.,phi0)
    y0 = np.r_[1.,phi0,A*vel0,(-b/e[...,None]).ravel()]
    solution = solve_ivp(rhs,(0,1),y0,method='DOP853',rtol=3e-13,atol=1e-14)
    if not solution.success:
        raise RuntimeError(solution.message)
    y = solution.y[:,-1]
    f_ref,e_ref,p_ref,_,_,v_ref,_ = measure(y)
    field = np.full(cfg.points,phi0)
    occupied = []
    for sector in slab.sectors:
        _,u = np.linalg.eigh(sector.full_hamiltonian(1.,field/phase_scale))
        occupied.append(u[:,:,:sector.d].reshape(len(sector.perp),sector.d,2,sector.d))
    initial = State(1.,field,np.full(cfg.points,A*vel0),tuple(occupied))
    results = []
    for steps in (100,200,400):
        state = initial
        obs = slab.observe(state)
        for _ in range(steps):
            state,obs,_ = slab.step(state,obs,1/steps)
        a,phi = state.a,np.mean(state.phi)
        d = slab.measure(state,obs)
        speed = np.mean(state.pi)/(a**3*d['A'])
        k,e,b,*_ = axes(a,phi)
        f = np.zeros_like(f_ref)
        for i,(sector,orbits) in enumerate(zip(slab.sectors,state.orbitals,strict=True)):
            blocks = sector.blocks(orbits)
            covariance = blocks @ adjoint(blocks)
            for l,p in enumerate(sector.perp):
                for j,qx in enumerate(sector.k):
                    _,vectors = np.linalg.eigh(qx*SIGMA_Z+p*SIGMA_X)
                    for h in (0,1):
                        u = np.kron(np.eye(2),vectors[:,h,None])
                        c = adjoint(u) @ covariance[l,j] @ u
                        bx,by,bz = b[i,l,j,h]
                        hamiltonian = np.array([[bz,bx-1j*by],[bx+1j*by,-bz]])
                        projector = 0.5*(np.eye(2)+hamiltonian/e[i,l,j,h])
                        f[i,l,j,h] = np.trace(projector @ c).real
        results.append(dict(steps=steps,max_occupation_error=float(np.max(abs(f-f_ref))),
                            field_spread=float(np.std(state.phi)),
                            relative_errors={key:float(abs(got/want-1)) for key,got,want in
                                             [('a',a,y[0]),('phi',phi,y[1]),('velocity',speed,v_ref),
                                              ('rho_F',d['rho_F'],e_ref),('P_F',d['P_F'],p_ref)]}))
    return results


def action():
    cfg = replace(Configuration(),longitudinal=2,transverse_cut=3.,transverse_bins=3,
                  points=25,scalar_modes=2,expand=False)
    slab = RegulatedSlab(cfg,MASSES)
    x = slab.spatial.x
    initial = slab.initial(0.11*np.sin(x),0.42*np.cos(x))
    a = 1.11
    state = replace(initial,a=a)
    direction = 0.3*np.cos(x)+0.2*np.sin(2*x)
    obs = slab.observe(state)
    eps = 1e-4
    fields = [replace(state,phi=state.phi+sgn*eps*direction) for sgn in (-1,1)]
    density = [slab.measure(s,slab.observe(s))['rho_F'] for s in fields]
    force = a**3*(density[1]-density[0])/(2*eps)
    force_error = abs(force-np.mean(obs.force*direction))
    d = slab.measure(state,obs)
    scales = [replace(state,a=a*math.exp(sgn*eps)) for sgn in (-1,1)]
    energies = [s.a**3*slab.measure(s,obs)['rho'] for s in scales]
    pressure_error = abs((energies[1]-energies[0])/(2*eps)/(-3*a**3*d['P'])-1)
    for dt in (0.01,-0.01):
        for _ in range(30):
            state = slab.matter(state,dt)
    inverse = max(np.max(abs(state.phi-initial.phi)),np.max(abs(state.pi-initial.pi)),
                  *(np.max(abs(x-y)) for x,y in zip(state.orbitals,initial.orbitals,strict=True)))
    return dict(force_error=float(force_error),pressure_relative_error=float(pressure_error),
                reversibility_max_error=float(inverse))


def vacuum():
    result={}
    for name,mass,expanding in (('static',0.7,False),('massless',0.,True)):
        cfg=replace(Configuration(),longitudinal=2,transverse_cut=3.,transverse_bins=3,
                    points=25,scalar_modes=2,mass=mass,expand=expanding)
        slab=RegulatedSlab(cfg,MASSES)
        x=slab.spatial.x
        state=slab.initial(0.13*np.sin(x)*expanding,0.65*np.cos(x)*expanding)
        obs=slab.observe(state)
        for _ in range(100):
            state,obs,_=slab.step(state,obs,0.01)
        d=slab.measure(state,obs,work=True)
        result[name]={k:float(d[k]) for k in ('rho_F','P_F','Jv','Jv_matched')}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    result=dict(helicity=helicity(),action=action(),vacuum=vacuum())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Signed helicity, reciprocity, reversibility and vacuum controls completed.')


if __name__=='__main__':
    main()
