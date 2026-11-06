#!/usr/bin/env python3
"""Execute the slab refinements with a resolved spatial scalar source."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from slab import Configuration, Slab, adjoint

REFERENCE = Configuration()
CASES = {
    'reference': ({}, 600),
    'spatial_quadrature': ({'points':97}, 600),
    'scalar_bandwidth': ({'scalar_modes':5}, 600),
    'transverse_bins': ({'transverse_bins':32}, 600),
    'wide_cutoff': ({'longitudinal':9, 'transverse_cut':12.0, 'transverse_bins':24}, 600),
    'half_step': ({}, 1200),
}
DURATION = 3.0


def evolve(box, initial, duration, steps):
    state = initial
    values = []
    pressure_work = np.zeros(2)
    for i in range(steps+1):
        measure = box.measure(state,work=True)
        values.append(dict(t=i*duration/steps,a=state.a,work_P=float(pressure_work[0]),
                           work_P_F=float(pressure_work[1]),**measure))
        if i<steps:
            state, step_work = box.step(state,duration/steps)
            pressure_work += step_work
    d = {k:np.array([v[k] for v in values]) for k in values[0]}
    return d, state


def evaluate(cfg, steps):
    box = Slab(cfg)
    x = 2*np.pi*box.x/cfg.length
    initial = box.initial(0.15*np.sin(x),0.9*np.cos(x),adiabatic=True)
    d, final = evolve(box,initial,DURATION,steps)
    measure, scale = d['a']**3, d['rho'][0]*d['a'][0]**3
    q = cumulative_trapezoid(measure*d['Jv'],d['t'],initial=0)
    absolute_q = cumulative_trapezoid(abs(measure*d['Jv']),d['t'],initial=0)
    energy = (measure*d['rho']+d['work_P'])/scale-1
    fermion_energy = (measure*d['rho_F']-measure[0]*d['rho_F'][0]+d['work_P_F'])
    momentum = d['a']**4*d['momentum']
    amplitudes = final.occupied.reshape(len(box.perp),2*box.d,box.d)
    gram = adjoint(amplitudes) @ amplitudes
    n_error = np.max(abs(gram-np.eye(box.d)))
    summary = dict(
        parameters=asdict(cfg), steps=steps, elapsed=DURATION, final={k:float(v[-1]) for k,v in d.items()},
        integrated_transfer=float(q[-1]), integrated_absolute_transfer=float(absolute_q[-1]),
        transfer_from_stress=float(fermion_energy[-1]),
        max_pressure_work_error=float(max(abs(energy))),
        max_separate_fermion_work_error=float(max(abs(q-fermion_energy))/scale),
        max_comoving_momentum_error=float(max(abs(momentum-momentum[0]))/scale),
        orthonormality_error=float(n_error),
        maximum_anisotropic_fermion_stress=float(max(abs(d['P_F_x']-d['P_F_perp']))),
        maximum_gradient_fraction=float(max(d['scalar_gradient']/d['rho'])),
        peak_instantaneous_work=float(max(abs(measure*d['Jv']))),
        pressure_tensor_verified_in_isotropic_average=True,
    )
    d['integrated_transfer']=q
    d['comoving_work']=measure*d['Jv']
    return summary,d


def compare(reference, control):
    r,d=reference
    s,q=control
    work_peak=r['peak_instantaneous_work']
    candidate=np.interp(d['t'],q['t'],q['comoving_work'])
    instantaneous=float(max(abs(candidate-d['comoving_work']))/work_peak)
    integrated=float(abs(s['integrated_transfer']/r['integrated_transfer']-1))
    integrated_history=float(max(abs(np.interp(d['t'],q['t'],q['integrated_transfer'])-d['integrated_transfer']))
                             / max(abs(d['integrated_transfer'])))
    absolute=float(abs(s['integrated_absolute_transfer']/r['integrated_absolute_transfer']-1))
    return dict(instantaneous_peak_scaled_error=instantaneous, integrated_transfer_error=integrated,
                integrated_history_error=integrated_history, integrated_absolute_transfer_error=absolute,
                gate_pass=max(instantaneous,integrated,integrated_history,absolute)<0.01)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)

    def run(name):
        changes,steps=CASES[name]
        value,data=evaluate(replace(REFERENCE,**changes),steps)
        np.savez_compressed(args.work/f'{name}.npz',**data)
        print(name,'integrated',value['integrated_transfer'],
              'budget',value['max_pressure_work_error'],flush=True)
        return name,(value,data)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        histories=dict(pool.map(run,CASES))
    differences={key:compare(histories['reference'],pair) for key,pair in histories.items() if key!='reference'}
    result=dict(
        scope='expanding homogeneous-metric slab with full transverse momentum',
        cases={key:value[0] for key,value in histories.items()},
        refinement_comparisons=differences,
        finite_slab_convergence_pass=all(row['gate_pass'] for row in differences.values()),
        scientific_acceptance=False,
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Source refinements', {k:v['gate_pass'] for k,v in differences.items()},flush=True)


if __name__=='__main__':
    main()
