#!/usr/bin/env python3
"""Reciprocal N1 bath feedback on one physical H55 spatial Fourier mode."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from coupled import Configuration, CoupledSlab, Geometry, State
from h55_probe import BEGIN, END, hessian, make_source, prepare, correlations
from h55_seed import archive, reconstruct
from run import read_columns

CASES={'reference':(4e-5,16,True),'half_step':(2e-5,16,True),
       'bath_spacing':(4e-5,24,True),'no_bath':(4e-5,16,False)}


def scalar_matrix(P,a,means,k,t):
    return hessian(P,*(f(t) for f in means))+(k/a)**2*np.eye(3)


def oscillator(z,p,a,A,dt):
    squared,u=np.linalg.eigh(A)
    freq=np.sqrt(squared.astype(complex))
    c=np.cos(freq*dt).real
    s=(dt*np.sinc(freq*dt/np.pi)).real
    x=u.T @ z
    v=u.T @ (p/a**3)
    return u @ (c*x+s*v),a**3*(u @ (c*v-s*squared*x))


def quantum(slab,a,phi,frame):
    data=slab.observe(State(a,phi,0*phi,frame),source=True)
    _,_,kin,bath,cross,_,_,_,mx,_,force=data
    K=float(np.sum(slab.weight*(kin+bath+cross)))
    m=slab.mass(a)
    return (K/a+m*mx)/2,(-K/a+slab.mass_slope(a)*mx)/2,m*force


def run(run_card,source,name):
    dt,bins,coupled=CASES[name]
    scale,means,mode,_=source
    card=run_card['channels'][0]
    M,f,width=card['M'],card['f'],card['Grest']
    cfg=Configuration(longitudinal=1,transverse_bins=1,transverse_cut=M,points=25,
                      scalar_modes=1,mass=M,phase_scale=f,length=run_card['L'],expand=False)
    gamma=width/(2*np.pi*M) if coupled else 0.
    slab=CoupledSlab(Geometry(cfg,bath_bins=bins,bath_cutoff=3*M,gamma=gamma))
    k=2*np.pi/cfg.length
    plane=np.exp(1j*k*slab.box.x)
    a=float(scale(BEGIN))
    y=mode(BEGIN)
    z,p=y[:3],a**3*y[3:]
    zero=np.ones_like(slab.box.x)

    def fields(z,t):
        return means[0](t)*zero+2*(z[0]*plane).real

    frame,_=prepare(slab,a,float(scale(BEGIN,1)),fields(z,BEGIN),
                    means[0](BEGIN,1)*zero+2*(y[3]*plane).real)
    homogeneous,_=prepare(slab,a,float(scale(BEGIN,1)),means[0](BEGIN)*zero,
                          means[0](BEGIN,1)*zero)
    count=round((END-BEGIN)/dt)
    times=np.linspace(BEGIN,END,count+1)
    dt=(END-BEGIN)/count
    energy=[]
    external=[]
    axial_work=[]
    delta_q=[]
    quantum_external=[]
    for i,t in enumerate(times):
        a=float(scale(t))
        H=float(scale(t,1))/a
        A=scalar_matrix(run_card['P'],a,means,k,t)
        eps=1e-6
        dA=(scalar_matrix(run_card['P'],float(scale(t+eps)),means,k,t+eps)-
            scalar_matrix(run_card['P'],float(scale(t-eps)),means,k,t-eps))/(2*eps)
        kinetic=float(np.vdot(p,p).real)/a**3
        potential=a**3*float(np.vdot(z,A @ z).real)
        eq,dq,j=quantum(slab,a,fields(z,t),frame)
        e0,d0,j0=quantum(slab,a,means[0](t)*zero,homogeneous)
        axis_external=3*H*(potential-kinetic)+a**3*float(np.vdot(z,dA @ z).real)
        q_external=H*(dq-d0)+float(means[0](t,1))*np.mean(j-j0)
        j1=float(np.mean(j*2*(p[0]*plane).real))/a**3
        energy.append(kinetic+potential+eq-e0)
        delta_q.append(eq-e0)
        external.append(axis_external+q_external)
        quantum_external.append(q_external)
        axial_work.append(j1)
        if i==count:
            continue
        mid=(times[i]+times[i+1])/2
        a0=float(scale(mid))
        A0=scalar_matrix(run_card['P'],a0,means,k,mid)
        z,p=oscillator(z,p,a0,A0,dt/2)
        frame=slab.quantum_free(a0,frame,dt/2)
        frame,kick=slab.local(a0,fields(z,mid),frame,dt,True)
        frame=slab.quantum_free(a0,frame,dt/2)
        p[0]+=np.mean(kick*plane.conj())
        z,p=oscillator(z,p,a0,A0,dt/2)
        homogeneous=slab.quantum_free(a0,homogeneous,dt/2)
        homogeneous,_=slab.local(a0,means[0](mid)*zero,homogeneous,dt,False)
        homogeneous=slab.quantum_free(a0,homogeneous,dt/2)
    work=cumulative_trapezoid(axial_work,times,initial=0)
    ext=cumulative_trapezoid(external,times,initial=0)
    qext=cumulative_trapezoid(quantum_external,times,initial=0)
    energy=np.array(energy)
    q=np.array(delta_q)
    ref=max(abs(energy[0]),1)
    expected=mode(END)[:3]
    return name,dict(dt=dt,bath_bins=bins,bath_coupled=coupled,
                     initial_axial_and_quantum_energy=float(energy[0]),
                     max_total_energy_work_error=float(np.max(abs(energy-energy[0]-ext))/ref),
                     max_quantum_energy_work_error=float(np.max(abs(q-q[0]-qext-work))/ref),
                     integrated_axial_work=float(work[-1]),terminal_quantum_energy_difference=float(q[-1]),
                     sigma_mode_terminal=[float(z[0].real),float(z[0].imag)],
                     sigma_difference_from_no_feedback=float(abs(z[0]-expected[0])/abs(expected[0])),
                     terminal_correlations=correlations(slab,frame))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('history',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    card,_,sha=archive()
    source=make_source(card,reconstruct(card),read_columns(args.history))
    with ThreadPoolExecutor(max_workers=4) as pool:
        cases=dict(pool.map(lambda name:run(card,source,name),CASES))
    out=dict(scope='physical H55 early axial mode and one N1 momentum shell, prescribed background',
             diagnostic_sha256=sha,source_trace_sha256=source[-1],interval=[BEGIN,END],
             cases=cases,scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('Physical H55 axial reciprocal probe completed.',flush=True)


if __name__=='__main__':
    main()
