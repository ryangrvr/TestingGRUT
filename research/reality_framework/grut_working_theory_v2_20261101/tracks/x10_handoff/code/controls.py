#!/usr/bin/env python3
"""Dense, stress and static controls for the coupled slab."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from coupled import Configuration, CoupledSlab, Geometry, State, adjoint
from handoff import evolve, ledger


def model(bins=2,**changes):
    cfg=Configuration(longitudinal=1,transverse_cut=1.6,transverse_bins=1,
                      points=17,scalar_modes=2,expand=False)
    return CoupledSlab(Geometry(replace(cfg,**changes),bath_bins=bins))


def homogeneous():
    slab=model()
    x=slab.box.x
    shift=0.37
    scale=slab.cfg.phase_scale
    T=0.5
    frame=slab.vacuum()
    z=np.broadcast_to(np.array([1.,-1.])[None,None,:,None],
                      (slab.n,2,2,slab.nb)).ravel()
    h=slab.full_hamiltonian(1.,0*x)[0]
    rot=h+shift*np.diag(z)
    want=(np.exp(1j*shift*T*z)[:,None]*(expm(-1j*rot*T) @ frame[0]))[None,:,:]
    cases=[]
    for n in (50,100,200):
        dt=T/n
        start=frame.copy()
        for j in range(n):
            phi=np.full_like(x,2*scale*shift*(j+0.5)*dt)
            start=slab.quantum_free(1.,start,dt/2)
            start,_=slab.local(1.,phi,start,dt,False)
            start=slab.quantum_free(1.,start,dt/2)
        cases.append(dict(steps=n,max_state_error=float(np.max(abs(start-want)))))
    phase=.72
    U=np.exp(0.5j*phase*z)
    target=(U[:,None]*h)*U.conj()[None,:]
    got=slab.full_hamiltonian(1.,np.full_like(x,phase*scale))[0]
    return dict(steps=cases,axial_rotation_error=float(np.max(abs(got-target))))


def spatial_ode():
    slab=model(bins=1)
    state=slab.prepare(.01)
    dim=len(state.phi)
    A,*_=slab.matching(1.)

    def pack(q):
        return np.r_[q.phi,q.pi,q.frame.ravel()]

    y0=pack(state)

    def rhs(_,y):
        phi,pi=y[:dim].real,y[dim:2*dim].real
        frame=y[2*dim:].reshape(state.frame.shape)
        inst=State(1.,phi,pi,frame)
        force=slab.mass(1.)*slab.observe(inst,source=True)[-1]
        laplace=slab.box.derivative(slab.box.derivative(phi))
        dframe=-1j*slab.full_hamiltonian(1.,phi) @ frame
        return np.r_[pi/A,A*laplace-slab.cfg.scalar_mass**2*phi-force,dframe.ravel()]

    T=0.4
    ref=solve_ivp(rhs,(0,T),y0,rtol=2e-13,atol=1e-14,method='DOP853')
    if not ref.success:
        raise RuntimeError(ref.message)
    rows=[]
    for n in (40,80,160):
        current=state
        obs=slab.observe(current,source=True)
        for _ in range(n):
            current,obs,_=slab.step(current,obs,T/n)
        rows.append(dict(steps=n,max_error=float(np.max(abs(pack(current)-ref.y[:,-1])))))
    for _ in range(n):
        current,obs,_=slab.step(current,obs,-T/n)
    return dict(comparisons=rows,reversal_error=float(np.max(abs(pack(current)-y0))))


def action():
    slab=model()
    state=replace(slab.prepare(.01),a=1.08)
    obs=slab.observe(state,source=True)
    current=slab.measure(state,obs,source=True)
    eps=1e-4
    delta=0.20*np.sin(slab.box.x)+0.08*np.cos(2*slab.box.x)
    values=[]
    for sign in (-1,1):
        changed=replace(state,phi=state.phi+sign*eps*delta)
        values.append(state.a**3*slab.measure(changed)['rho_q'])
    field_error=abs((values[1]-values[0])/(2*eps)-slab.mass(state.a)*np.mean(obs[-1]*delta))
    energies=[]
    for sign in (-1,1):
        changed=replace(state,a=state.a*np.exp(sign*eps))
        energies.append(changed.a**3*slab.measure(changed,obs)['rho'])
    pressure_error=abs((energies[1]-energies[0])/(2*eps)/(-3*state.a**3*current['P'])-1)
    sample=(1.031,1.072,1.117,1.163,1.204,1.251)
    residual=[]
    for a in sample:
        energy,slope,aux=slab.vacuum_row(a)
        residual.append(dict(a=a,energy=float(abs(energy-slab.vac_energy(np.log(a)))),
                             loga_derivative=float(abs(slope-slab.vac_energy(np.log(a),1))),
                             subtraction=float(np.max(abs(aux-slab.vac_aux(np.log(a)))))))
    return dict(field_variation_error=float(field_error),pressure_relative_error=float(pressure_error),
                vacuum_table=residual)


def vacuum():
    data={}
    for name,slab in (('static',model()),('massless',model(mass=0.,expand=True))):
        if name=='static':
            zero=np.zeros_like(slab.box.x)
            initial=State(1.,zero,zero,slab.vacuum())
        else:
            initial=slab.prepare(.01)
        trace,_=evolve(slab,initial,0.,.01,end=.5)
        data[name]=dict(ledger=ledger(trace),
                        max_quantum_energy=float(np.max(abs(trace['rho_q']))),
                        max_quantum_pressure=float(np.max(abs(trace['P_q']))),
                        max_quantum_work=float(np.max(abs(trace['work']))))
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    result=dict(homogeneous=homogeneous(),spatial_ode=spatial_ode(),
                action=action(),vacuum=vacuum(),scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Independent coupled controls completed.',flush=True)


if __name__=='__main__':
    main()
