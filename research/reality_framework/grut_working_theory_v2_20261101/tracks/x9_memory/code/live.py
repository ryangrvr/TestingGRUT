#!/usr/bin/env python3
"""Conserving live scalar control, dense reference, and stationary phase averaging."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from bath import build, engine, evolve, mass_matched, prepare


def observables(star,m,c):
    n=len(star.frequency)
    prob=np.sum(abs(c)**2,axis=1)
    sterile=c[:2]
    bath=c[2:].reshape(2,n,-1)
    dx=float(2*np.vdot(sterile[0],sterile[1]).real)
    z=prob[0]-prob[1]
    sterile_energy=m*dx+star.k*z
    bprob=prob[2:].reshape(2,n)
    bath_energy=float(np.sum(bprob*star.frequency))
    mix=2*float(np.sum(sterile.conj()*np.einsum('bi,bic->bc',star.coupling,bath)).real)
    return dict(fermion_energy=float(sterile_energy),bath_energy=bath_energy,interaction=mix,
                charge=float(z+np.sum(bprob[0]-bprob[1])),sterile_count=float(prob[:2].sum()),
                bath_positive_count=float(np.sum(bprob[:,star.frequency>0])),
                bath_negative_count=float(np.sum(bprob[:,star.frequency<0])))


def discrete_step(star,m,c,phi,p,dt,inertia,f,curvature,charge):
    """One constant-speed rotating solve with an energy-consistent scalar gradient."""
    velocity=p/inertia
    z=2j/dt
    for it in range(20):
        delta=velocity/(2*f)
        next_phi=phi+velocity*dt
        gradient=curvature*(phi+next_phi)/2
        new=2*z*star.resolve(m,z,c,delta)-c
        updated=observables(star,m,new)
        new_p=p+(updated['charge']-charge)/(2*f)-dt*gradient
        midpoint=(p+new_p)/(2*inertia)
        if abs(midpoint-velocity)<2e-13*(1+abs(midpoint)):
            return new,next_phi,new_p,updated,gradient,it+1
        velocity=midpoint
    raise RuntimeError('Scalar midpoint did not converge')


def scalar_case(star,m,c,dt_goal,curvature):
    inertia,f=6.,2.5
    phi,p=0.36,2.4*inertia
    T=8.
    steps=int(np.ceil(T/dt_goal))
    dt=T/steps
    values=[]
    charge_integral=0
    max_iterations=0
    obs=observables(star,m,c)
    pzero=p-obs['charge']/(2*f)
    energy0=None
    worst_norm=0.
    for i in range(steps+1):
        scalar=p*p/(2*inertia)+curvature*phi*phi/2
        energy=scalar+obs['fermion_energy']+obs['bath_energy']+obs['interaction']
        if energy0 is None:
            energy0=energy
        residual=(energy-energy0)/max(1.,abs(energy0))
        charge=p-obs['charge']/(2*f)+charge_integral-pzero
        values.append(dict(t=i*dt,phi=phi,velocity=p/inertia,scalar_energy=scalar,
                           energy=energy,residual=residual,charge_residual=charge,**obs))
        if i<steps:
            c,phi,p,obs,gradient,it=discrete_step(star,m,c,phi,p,dt,inertia,f,curvature,obs['charge'])
            charge_integral+=dt*gradient
            max_iterations=max(it,max_iterations)
        if i%100==0 or i==steps:
            worst_norm=max(worst_norm,float(np.max(abs(c.conj().T @ c-np.eye(c.shape[1])))))
    result=dict(dt=dt,initial=values[0],terminal=values[-1],
                max_energy_error=max(abs(r['residual']) for r in values),
                max_charge_error=max(abs(r['charge_residual']) for r in values),
                max_pauli_projector_error=worst_norm,iterations=max_iterations)
    return result


def full_feedback():
    M,G,k,cutoff=1.4,0.08,1.1,10.
    gamma=G/(2*np.pi*M)
    m=mass_matched(M,gamma,cutoff)
    star,spacing=build(k,gamma,cutoff,0.36)
    energy,vecs=np.linalg.eigh(star.matrix(m))
    filled=vecs[:,energy<0]
    pair,prep=prepare(star,m,M,width=0.65)
    exact=vecs[:,energy>0]
    positive=exact @ (exact.T @ pair.T)
    prep['dense_positive_error']=float(np.max(abs(positive-pair.T)))
    frame=[]
    for i in (0,1):
        frame.append(np.column_stack([filled,pair[i]]))
    runs={}
    for tag,dt in (('base',0.04),('half',0.02)):
        data=[scalar_case(star,m,c,dt,curvature=0.65) for c in frame]
        keys=('scalar_energy','fermion_energy','bath_energy','interaction',
              'sterile_count','bath_positive_count','bath_negative_count','phi','velocity')
        difference={part:{k:data[0][part][k]-data[1][part][k] for k in keys}
                    for part in ('initial','terminal')}
        runs[tag]=dict(cases=data,differences=difference)
    return dict(regulator=cutoff,frequency_spacing=spacing,dimension=len(energy),
                vacuum_width=G,initial_preparation=prep,runs=runs)


def dense_and_born(work):
    lib=engine(work)
    M,k,delta=1.3,1.1,0.43
    G=0.03
    cutoff=8.
    m=mass_matched(M,G/(2*np.pi*M),cutoff)
    star,_=build(k,G/(2*np.pi*M),cutoff,0.33)
    pair,prep=prepare(star,m,M,width=0.7)
    T=3.0
    exact=(expm(-1j*T*star.matrix(m,delta)) @ pair.T).T
    splits=[]
    for steps in (300,600,1200):
        eta=np.linspace(0,T,steps+1)
        path=np.tile([m,delta,0.,M],(steps+1,1))
        _,d,end,_=evolve(lib,star,path,eta,pair,steps)
        splits.append(dict(steps=steps,error=float(np.max(abs(end-exact))),
                           norm_error=float(np.max(abs(d['norm']-1))),
                           energy_work_error=float(np.max(abs(d['total_energy'][-1]-
                                                              d['total_energy'][0]-d['scalar_work'][-1])))))
    eta=np.linspace(0,T,121)
    path=np.tile([m,delta,0.,M],(len(eta),1))
    _,d,end,_=evolve(lib,star,path,eta,pair,len(eta)-1,pade=True)
    pade=dict(steps=120,error=float(np.max(abs(end-exact))),
              max_energy_work_error=float(np.max(abs(d['total_energy'][-1]-d['total_energy'][0]-
                                                     d['scalar_work'][-1]))))
    weak,_=build(k,1e-4*G/(2*np.pi*M),cutoff,0.33)
    e=np.hypot(m,k+delta)
    spin=np.sqrt((1+np.array([1.,-1.])*(k+delta)/e)/2)
    start=np.zeros((2,2+2*len(weak.frequency)),complex)
    start[:,:2]=spin
    eta=np.linspace(0,T,4801)
    path=np.tile([m,delta,0.,M],(len(eta),1))
    _,_,end,_=evolve(lib,weak,path,eta,start,len(eta)-1)
    branch=weak.frequency[None,:]+np.array([delta,-delta])[:,None]-e
    prediction=weak.coupling**2*spin[:,None]**2*(T*np.sinc(branch*T/(2*np.pi)))**2
    measured=abs(end[0,2:].reshape(prediction.shape))**2
    return dict(split=splits,pade=pade,leading_born_relative_error=float(np.sum(abs(measured-prediction))/prediction.sum()),
                first_order_bath_probability=float(prediction.sum()),projector=prep)


def stationary_control(work):
    lib=engine(work)
    k,M,G=0.7,2.,0.15
    gamma=G/(2*np.pi*M)
    cutoff=15.
    m=mass_matched(M,gamma,cutoff)
    outcome=[]
    for spacing in (0.025,0.0125):
        star,actual=build(k,gamma,cutoff,spacing)
        pair,prep=prepare(star,m,M,width=0.4)
        eta=np.linspace(0,120,12001)
        path=np.tile([m,0.,0.,M],(len(eta),1))
        keep,d,_,_=evolve(lib,star,path,eta,pair,25)
        late=eta[keep]>=100
        early=eta[keep]<=30
        ster=abs(d['sterile_free_count'][:,0]-d['sterile_free_count'][:,1])
        daughter=abs(d['bath_chiral_count'][:,0]-d['bath_chiral_count'][:,1])
        total=abs(d['bath_positive'][:,0]-d['bath_positive'][:,1])
        outcome.append(dict(
            frequency_spacing=actual,recurrence_time=2*np.pi/actual,
            max_early_sterile_count_difference=float(max(ster[early])),
            max_late_sterile_count_difference=float(max(ster[late])),
            max_late_total_positive_bath_count_difference=float(max(total[late])),
            mean_late_chiral_bath_count_difference=float(np.mean(daughter[late])),
            max_static_energy_error=float(np.max(abs(d['total_energy']-d['total_energy'][0]))),
            max_norm_error=float(np.max(abs(d['norm']-1))),**prep,
        ))
    # Exact finite-star dephasing of conjugate states for real observables.
    star,_=build(k,gamma,cutoff,0.25)
    h=star.matrix(m)
    _,u=np.linalg.eigh(h)
    pair,_=prepare(star,m,M,width=0.4)
    populations=abs(u.T @ pair.T)**2
    return dict(runs=outcome,exact_stationary_spectral_weight_difference=float(
        np.max(abs(populations[:,0]-populations[:,1]))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    result=dict(scope='finite conserving scalar feedback and stationary controls for the quadratic vacuum bath',
                microscopic=dense_and_born(work),stationary=stationary_control(work),
                live=full_feedback(),scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Live scalar and stationary controls completed.',flush=True)


if __name__=='__main__':
    main()
