#!/usr/bin/env python3
"""Collect the H55 runs and compare source and continuation errors separately."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from production import ARCHIVE, CASES, DT, END, TRACKS, h55_model, read_columns, sha
from run import compile_engine, run_engine

MOMENTS=('N_fraction','N_count_comoving',*(f'{key}_{i}' for key in ('rho','n') for i in range(1,5)))


def json_read(path):
    return json.loads(path.read_text())


def independent_homogeneous(folder, reference):
    folder.mkdir(parents=True,exist_ok=True)
    meta=json_read(ARCHIVE)
    model=h55_model(meta,(360,45000,DT,True,False,True))
    model.update(t1=END,nsteps=math.ceil((END-model['t0'])/DT))
    model['stride']=round(.0004/((END-model['t0'])/model['nsteps']))
    data,_,_=run_engine(compile_engine(folder),folder,'x5',model)
    if not np.array_equal(data['t'],reference['t']):
        raise ValueError('Homogeneous control times differ')
    columns=('a','rho','sigma','chi',*(f'{key}_{i}' for key in ('rho','n') for i in range(1,5)))
    errors={k:float(max(abs(data[k]-reference[k]))/max(abs(reference[k]))) for k in columns}
    return dict(comparison='direct X5 homogeneous split',max_error=max(errors.values()),history_errors=errors)


def quantum_measures(data):
    t,a,rho=(data[k] for k in ('t','a','rho'))
    energy=sum(data[f'rho_{j}'] for j in range(1,5))
    pressure=sum(data[f'p_{j}'] for j in range(1,5))
    work=sum(data[f'work_{j}'] for j in range(1,5))/rho
    # The pressure integral uses the recorded pressure and volume. It
    # resolves the accumulated transfer even if output cadence aliases J.
    net=energy[-1]*a[-1]**3-energy[0]*a[0]**3+trapezoid(pressure,a**3)
    t50=np.flatnonzero(energy/rho>=.5)
    return dict(work=work,first_sample_over_half=float(t[t50[0]]) if len(t50) else None,
                integrated_quantum_stress_transfer=float(net))


def compare(cases, data, quantum, a, b):
    ta,wa=data[a]['t'],quantum[a]['work']
    tb,wb=data[b]['t'],quantum[b]['work']
    i=(ta>=tb[0])&(ta<=tb[-1])
    return dict(
        baseline=a,variation=b,
        terminal_relative_changes={k:cases[b]['terminal'][k]/cases[a]['terminal'][k]-1 for k in MOMENTS},
        signed_work_trace_difference_over_baseline_peak=float(max(abs(np.interp(ta[i],tb,wb)-wa[i]))/max(abs(wa[i]))),
        terminal_signed_work_difference_over_baseline_magnitude=float((wb[-1]-wa[-1])/abs(wa[-1])),
        integrated_transfer_relative_change=quantum[b]['integrated_quantum_stress_transfer']/
                                             quantum[a]['integrated_quantum_stress_transfer']-1)


def production(root):
    work=root/'production'
    cases={name:json_read(work/name/'summary.json') for name in CASES}
    for row in cases.values():
        row.setdefault('completed','terminal' in row)
        if row['completed'] and 'max_reciprocal_work_ledger_error' not in row:
            row['sampled_quadrature_scalar_work']=row.pop('integrated_comoving_scalar_work')
    traces={name:read_columns(work/name/'history.txt') for name in cases if 'terminal' in cases[name]}
    q={name:quantum_measures(d) for name,d in traces.items()}
    pairs=dict(timestep=('reference','half_step'),radial=('reference','radial'),
               cutoff=('half_step','cutoff_half_step'),spatial=('reference','spatial'),
               homogeneous=('reference','homogeneous'),late_spatial=('late_homogeneous','late'))
    differences={name:compare(cases,traces,q,*p) for name,p in pairs.items()}
    reference=traces['reference']
    n=cases['reference']['grid']
    L=json_read(ARCHIVE)['run']['L']
    t=np.r_[.18,reference['t'][(reference['t']>.18)&(reference['t']<.35)],.35]
    travel=float(trapezoid(1/np.interp(t,reference['t'],reference['a']),t))
    metrics={k:{name:value for name,value in row.items() if name!='work'} for k,row in q.items()}
    sources=('x2_strip/results/lower_H_lattice.json','x2_strip/results/lower_H_diag.txt.gz',
             'lf_fermion/code/lfcommon.py','x5_coupled/code/coupled.c','x5_coupled/code/run.py',
             'x10_handoff/code/h55_seed.py','e2e_kinetic/code/ek_modes.c',
             'x4_closure/code/kinetic_engine.py','x3_actual/results/kinetic/actual_native_6400.in.gz',
             'mf_ideal/code/rpmfi.c','e2e_relic/code/relic.py','e2e_card/CANDIDATE_M2.json')
    return dict(scope='H55 shared FRW scalar lattice with site-local coherent sterile production',
                cases=cases,comparisons=differences,work_from_quantum_stress=metrics,
                light_travel_comoving_Hi_inverse=travel,travel_in_archived_scalar_cells=travel/(L/64),
                travel_in_reference_cells=travel/(L/n),travel_interval_Hi_inverse=[.18,.35],
                independent_homogeneous=independent_homogeneous(root/'control',traces['homogeneous']),
                input_sha256={s:sha(TRACKS/s) for s in sources},
                spatially_controlled_H55_production=False)


def continuations(root):
    cases={stage:{p.parent.name:json_read(p) for p in sorted((root/stage).glob('*/summary.json'))}
           for stage in ('kinetic','abundance','baryon')}
    if not all(cases.values()):
        raise ValueError('Missing kinetic, abundance or baryon archives')
    a=cases['abundance']
    out=dict(cases=cases,accepted_abundance=False,accepted_GRUT_law=False,strict_candidate_ledger='0/135')
    if {'at035','at055','late','at055_half'}<=set(a):
        ref=a['at055']['conditional_Omega_h2']
        out['boundary_relative_abundance_difference']={k:a[k]['conditional_Omega_h2']/ref-1
                                                       for k in ('at035','late')}
        out['kinetic_step_relative_abundance_difference']=a['at055_half']['conditional_Omega_h2']/ref-1
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    for name,value in (('production',production(args.work)),('continuation',continuations(args.work))):
        (args.output/f'{name}.json').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
        print(f'{name} archived',flush=True)


if __name__=='__main__':
    main()
