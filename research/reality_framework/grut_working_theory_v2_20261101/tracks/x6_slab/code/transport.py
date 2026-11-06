#!/usr/bin/env python3
"""Coherent overlap and a pressure-certified particle transition in a control slab."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from slab import Configuration, Slab, State, adjoint

PULSE_DURATION = 3.0
PULSE_FREQUENCY = 2.7
REST_RATE = 0.015
DEPHASING = 0.6
TOLERANCE = 1e-6
TAIL = np.linspace(0,80,1601)
TRACK = Path(__file__).resolve().parents[1]
X5 = TRACK.parent/'x5_coupled'


def pulse_profile(x,time):
    if time<=0 or time>=PULSE_DURATION:
        zero=np.zeros_like(x)
        return zero,zero,zero
    k=math.pi/PULSE_DURATION
    sn=math.sin(k*time)
    amplitude=sn**4
    derivative=4*k*sn**3*math.cos(k*time)
    angle=x-PULSE_FREQUENCY*time
    phase=amplitude*np.cos(angle)
    slope=derivative*np.cos(angle)+amplitude*PULSE_FREQUENCY*np.sin(angle)
    spatial=-amplitude*np.sin(angle)
    return phase,slope,spatial


def pulse(box,steps):
    times=np.linspace(0,PULSE_DURATION,steps+1)
    zero=np.zeros(box.cfg.points)
    state=box.initial(zero,zero,adiabatic=False)
    occupied=state.occupied
    dt=times[1]-times[0]
    values=[]
    for i,t in enumerate(times):
        phase,tdot,xdot=pulse_profile(box.x,t)
        vertex=box.mass_derivative(phase,box.cov(occupied))
        phi=phase*box.cfg.phase_scale
        moments=box.measure(State(1.0,phi,zero,occupied))
        power=np.mean(vertex*tdot*box.cfg.phase_scale)
        momentum_rate=-np.mean(vertex*xdot*box.cfg.phase_scale)
        values.append((moments['rho_F'],moments['fermion_momentum'],power,momentum_rate))
        if i<steps:
            midpoint,_,_=pulse_profile(box.x,t+dt/2)
            occupied=box.free_step(1.0,occupied,dt/2)
            occupied,_=box.local(midpoint,occupied,dt,False)
            occupied=box.free_step(1.0,occupied,dt/2)
    data=np.array(values)
    end_e,end_p=data[-1,:2]
    work=np.trapezoid(data[:,2],times)
    impulse=np.trapezoid(data[:,3],times)
    assert end_e>0 and abs(end_p)>0.01*end_e
    result=dict(steps=steps,energy=float(end_e),momentum=float(end_p),injected_work=float(work),
                injected_momentum=float(impulse),work_relative_error=float(abs(end_e-data[0,0]-work)/end_e),
                momentum_relative_error=float(abs(end_p-data[0,1]-impulse)/end_e),
                phase_and_derivatives_zero_after_pulse=True)
    return occupied,result


def daughter_pressure(box,energy):
    """Isotropic two-body decay in the rest frame, boosted before averaging."""
    p=np.sqrt(box.q2)
    beta=p/energy
    nodes,weights=leggauss(96)
    fraction=0.5*((beta[...,None]+nodes)**2/(1+beta[...,None]*nodes) @ weights)
    transverse=(1-fraction)/2
    x=energy*(transverse+(fraction-transverse)*(box.k[None,:]/p)**2)
    assert np.min(x)>0 and np.max(x/energy)<=1
    return x,(energy-x)/2


class FreeTail:
    def __init__(self,box,occupied):
        self.box=box
        self.energy=box.energy(1)
        self.rate=REST_RATE*box.cfg.mass/self.energy
        self.weight=box.weight[:,None]/2
        self.eigen,self.vectors=np.linalg.eigh(box.free_blocks(1))
        states=adjoint(self.vectors) @ box.blocks(occupied)
        states=states.reshape(len(box.perp),4*box.n,box.d)
        self.vacuum=np.tile([1,1,0,0],box.n)
        self.delta=states @ adjoint(states)-np.diag(self.vacuum)
        diagonal=np.diagonal(self.delta,axis1=-2,axis2=-1).real.reshape(len(box.perp),box.n,4)
        self.particles=diagonal[...,2:].sum(axis=-1)
        self.holes=-diagonal[...,:2].sum(axis=-1)
        self.energy_scale=float(np.sum(self.weight*self.energy*(self.particles+self.holes)))
        self.number_scale=float(np.sum(self.weight*(self.particles+self.holes)))
        self.parent_x=box.k[None,:]**2/self.energy
        self.parent_p=box.perp[:,None]**2/(2*self.energy)
        self.daughter_x,self.daughter_p=daughter_pressure(box,self.energy)
        self.stress=(adjoint(self.vectors) @ box.kx @ self.vectors,
                     adjoint(self.vectors) @ (box.kt/2) @ self.vectors)
        diagonal_blocks=self.blocks(self.delta)
        self.coefficients=np.stack([
            2*self.weight*np.einsum('...ij,...ji->...',m[...,:2,2:],diagonal_blocks[...,2:,:2])
            for m in self.stress])

    def blocks(self,covariance):
        n=self.box.n
        c=covariance.reshape(len(self.box.perp),n,4,n,4)
        return np.stack([c[:,j,:,j,:] for j in range(n)],axis=1)

    def covariance(self,t,dephase):
        rates=np.repeat(self.rate,4,axis=1)
        energies=self.eigen.reshape(len(self.box.perp),-1)
        decay=np.exp(-0.5*rates*t-1j*energies*t)
        opposite=self.vacuum[:,None]!=self.vacuum[None,:]
        mix=np.exp(-dephase*t*opposite)
        return np.diag(self.vacuum)+decay[:,:,None]*self.delta*decay[:,None,:].conj()*mix

    def gas(self,elapsed,start=0.0,daughters=None):
        t=np.atleast_1d(elapsed)
        initial=np.exp(-self.rate*start)
        later=np.exp(-t[:,None,None]*self.rate)
        number=initial*(self.particles+self.holes)
        difference=initial*(self.particles-self.holes)
        population=later*number
        lost=(1-later)*number

        def integrate(f):
            return np.sum(self.weight*f,axis=(-2,-1))

        result=dict(
            number=integrate(population),
            energy_F=integrate(population*self.energy),
            momentum_F=integrate(later*difference*self.box.k[None,:]),
            P_F_x=integrate(population*self.parent_x),
            P_F_perp=integrate(population*self.parent_p),
            energy_R=integrate(lost*self.energy),
            momentum_R=integrate((1-later)*difference*self.box.k[None,:]),
            P_R_x=integrate(lost*self.daughter_x),
            P_R_perp=integrate(lost*self.daughter_p),
        )
        if daughters is not None:
            for key,value in daughters.items():
                result[key]+=value
        result['energy_total']=result['energy_F']+result['energy_R']
        result['momentum_total']=result['momentum_F']+result['momentum_R']
        result['P_x']=result['P_F_x']+result['P_R_x']
        result['P_perp']=result['P_F_perp']+result['P_R_perp']
        return result

    def coherence_pressure(self,times,dephase):
        t=np.atleast_1d(times)
        envelope=np.exp(-t[:,None,None]*(self.rate+dephase))
        phase=np.exp(-2j*t[:,None,None]*self.energy)
        p=np.sum(self.coefficients[:,None]*envelope*phase,axis=(-2,-1)).real.T
        bound=np.sum(abs(self.coefficients[:,None])*envelope,axis=(-2,-1)).T
        return p,bound

    def full_check(self,time,dephase):
        c=self.covariance(time,dephase)
        eig=np.linalg.eigvalsh(c)
        assert min(eig.ravel())>-1e-10 and max(eig.ravel())<1+1e-10
        reduced=self.blocks(c)-np.diag([1,1,0,0])
        direct=np.array([np.sum(self.weight*np.einsum('...ij,...ji->...',m,reduced).real)
                         for m in self.stress])
        gas=self.gas(time)
        pressure,bound=self.coherence_pressure(time,dephase)
        expected=pressure[0]+np.array([gas['P_F_x'][0],gas['P_F_perp'][0]])
        error=float(max(abs(direct-expected))/self.energy_scale)
        assert error<1e-11
        # Only spatially averaged observables were matched in this tail.
        c4=c.reshape(len(self.box.perp),self.box.n,4,self.box.n,4)
        q_nonlocal=~np.eye(self.box.n,dtype=bool)
        same_band=(np.array([1,1,0,0])[:,None]==np.array([1,1,0,0])[None,:])
        mask=q_nonlocal[:,None,:,None]*same_band[None,:,None,:]
        nonlocal_norm=float(np.linalg.norm(c4*mask)/np.linalg.norm(self.delta))
        return dict(t=float(time),covariance_eigenvalue_range=[float(min(eig.ravel())),float(max(eig.ravel()))],
                    pressure_reference_relative_error=error,
                    remaining_same_band_spatial_coherence_relative_norm=nonlocal_norm)


def assess(tail,dephase):
    scale=tail.energy_scale
    gas=tail.gas(TAIL)
    pressure,bound=tail.coherence_pressure(TAIL,dephase)
    pressure_norm=np.max(bound,axis=1)/scale
    surviving=gas['energy_F']/scale
    assert np.all(abs(pressure)<=bound+1e-15)
    energy_residual=float(max(abs(gas['energy_total']-scale))/scale)
    mom_residual=float(max(abs(gas['momentum_total']-gas['momentum_total'][0]))/scale)
    candidates=np.flatnonzero((pressure_norm<TOLERANCE)&(surviving>0.5))
    earliest=int(candidates[0]) if len(candidates) else None
    boundaries=[0] if earliest is None else [0,earliest,earliest+50,earliest+100]
    full=dict(gas)
    full['P_x']=gas['P_x']+pressure[:,0]
    full['P_perp']=gas['P_perp']+pressure[:,1]
    trials=[]
    accepted_series=[]
    for idx in boundaries:
        if idx>=len(TAIL):
            raise ValueError('Boundary exceeds downstream interval')
        time=TAIL[idx]
        carried={k:gas[k][idx] for k in ('energy_R','momentum_R','P_R_x','P_R_perp')}
        predicted=tail.gas(TAIL[idx:]-time,start=time,daughters=carried)
        errors={k:float(max(abs(v-full[k][idx:]))/(tail.number_scale if k=='number' else scale))
                for k,v in predicted.items() if k not in ('P_F_x','P_F_perp')}
        adequate=(pressure_norm[idx]<TOLERANCE and surviving[idx]>0.5
                  and max(errors.values())<TOLERANCE)
        trials.append(dict(boundary=float(time),pressure_bound_over_initial_energy=float(pressure_norm[idx]),
                           surviving_parent_energy_fraction=float(surviving[idx]),
                           maximum_downstream_relative_errors=errors,
                           final_daughter_energy=float(predicted['energy_R'][-1]),
                           final_daughter_momentum=float(predicted['momentum_R'][-1]),
                           accepted=bool(adequate)))
        if adequate:
            accepted_series.append((idx,predicted))
    if earliest is not None:
        assert all(x['accepted'] for x in trials[1:])
    points=sorted({0,200,len(TAIL)-1,*boundaries})
    checks=[tail.full_check(TAIL[idx],dephase) for idx in points]
    half_time=brentq(lambda t: tail.gas(t)['energy_F'][0]/scale-0.5,0,10000,xtol=1e-9)
    _,half_envelope=tail.coherence_pressure(half_time,dephase)
    boundary_difference=None
    if accepted_series:
        base_idx,base=accepted_series[0]
        boundary_difference=max(float(max(abs(x[key]-base[key][idx-base_idx:]))
                                      /(tail.number_scale if key=='number' else scale))
                                for idx,x in accepted_series[1:] for key in x)
    return dict(
        elastic_dephasing_rate=dephase, max_energy_relative_error=energy_residual,
        max_momentum_relative_error=mom_residual,
        boundary=None if earliest is None else float(TAIL[earliest]),
        trials=trials, full_covariance_controls=checks,
        max_difference_between_accepted_boundaries=boundary_difference,
        parent_energy_half_time=float(half_time),
        pressure_envelope_at_parent_half_time=float(max(half_envelope[0])/scale),
        min_pressure_bound_while_half_energy_remains=float(min(pressure_norm[surviving>0.5])),
        final_daughter_energy=float(gas['energy_R'][-1]),
        final_daughter_momentum=float(gas['momentum_R'][-1]),
        final_daughter_pressure_x=float(gas['P_R_x'][-1]),
        final_daughter_pressure_perp=float(gas['P_R_perp'][-1]),
    )


def inherited_h55():
    audit=json.loads((X5/'results/diagnosis.json').read_text())
    fine=next(row for row in audit['handoffs'] if row['case']=='eight_bins')
    energy=fine['replacing_subtracted_stress_with_dressed_energy_gap']
    pressure=fine['total_pressure_change_if_dephased_over_rho']
    work=fine['total_coherent_fermion_power_over_rho']
    source=json.loads((X5/'results/h55_refined.json').read_text())
    hubble=next(r['terminal']['H'] for r in source['cases'] if r['case']=='eight_bins')
    derivative_jump=-work-3*hubble*pressure
    gate=abs(energy)<TOLERANCE and abs(pressure)<TOLERANCE and abs(work)<TOLERANCE
    return dict(time=fine['t'],dressed_energy_gap=energy,free_dephasing_pressure_gap=pressure,
                coherent_work_over_rho=work,normalized_fermion_energy_derivative_jump=derivative_jump,
                terminal_refinement_pass=audit['numerical_terminal_pass'],
                tested_boundary_pass=gate,applied_control_tail_to_H55=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    box=Slab(replace(Configuration(),expand=False,subtract=False))
    _,coarse=pulse(box,600)
    occupied,fine=pulse(box,1200)
    tail=FreeTail(box,occupied)
    assert abs(tail.energy_scale/fine['energy']-1)<1e-10
    assert fine['work_relative_error']<2e-5 and fine['momentum_relative_error']<2e-5
    elastic=assess(tail,DEPHASING)
    no_elastic=assess(tail,0.0)
    h55=inherited_h55()
    assert not h55['tested_boundary_pass']
    assert elastic['boundary'] is not None and no_elastic['boundary'] is None
    result=dict(scope='fixed-volume controlled preparation with coherent overlap and specified test interactions',
                pulse=[coarse,fine],parent_rest_decay_rate=REST_RATE,
                tolerance=TOLERANCE,controlled_elastic=elastic,controlled_no_elastic=no_elastic,
                inherited_H55_interface=h55,local_transport_adequacy_claimed=False,scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Earliest controlled boundary',elastic['boundary'],'without dephasing',no_elastic['boundary'],
          'H55 boundary',h55['tested_boundary_pass'],flush=True)


if __name__=='__main__':
    main()
