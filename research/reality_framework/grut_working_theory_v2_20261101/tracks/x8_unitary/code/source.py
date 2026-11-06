#!/usr/bin/env python3
"""Run the single-mass unitary slab at independent source resolutions."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

from unitary import Configuration, UnitarySlab, adjoint
from evolution import compare
from paired import PairedSlab
from preparation import ThirdOrderSlab

DURATION=3.
# Longitudinal cutoff, transverse cutoff, Gaussian order, collocation, scalar band, steps.
CASES={
    'low':(9,12,12,49,3,1200),
    'reference':(12,16,12,65,3,1200),
    'wide':(16,22,12,81,3,1200),
    'longitudinal':(16,16,12,81,3,1200),
    'transverse':(12,22,12,65,3,1200),
    'radial':(12,16,8,65,3,1200),
    'half_step':(12,16,12,65,3,2400),
    'collocation':(12,16,12,129,3,1200),
    'bandwidth':(12,16,12,65,5,1200),
}
PAIRED_CASES={
    'paired_low':(9,12,12,49,3,2400),
    'paired_reference':(12,16,12,65,3,2400),
    'paired_wide':(16,22,12,81,3,2400),
    'paired_longitudinal':(16,16,12,81,3,2400),
    'paired_transverse':(12,22,12,65,3,2400),
    'paired_radial':(12,16,8,65,3,2400),
    'paired_half_step':(12,16,12,65,3,4800),
    'paired_collocation':(12,16,12,129,3,2400),
    'paired_bandwidth':(12,16,12,65,5,2400),
}

PREPARATION_CASES={
    'uv_reference':(12,16,12,65,3,2400,(3,5)),
    'uv_wide':(16,22,12,81,3,2400,(3,5)),
    'uv_later':(12,16,12,65,3,2400,(4,6)),
    'uv_wide_later':(16,22,12,81,3,2400,(4,6)),
    'uv_wide_half':(16,22,12,81,3,4800,(3,5)),
}


def evaluate(cfg,order,steps,duration=DURATION,paired=False,transition=None):
    if transition is not None:
        box=ThirdOrderSlab(cfg,order,transition)
    else:
        box=(PairedSlab if paired else UnitarySlab)(cfg,order)
    x=box.spatial.x
    state=box.initial(0.15*np.sin(x),0.9*np.cos(x))
    observation=box.observe(state)
    records=[]
    pressure=np.zeros(3)
    dt=duration/steps
    for i in range(steps+1):
        result=box.measure(state,observation,work=True)
        result.pop('acceleration')
        records.append(dict(t=i*dt,work_P=pressure[0],work_P_F=pressure[1],
                            work_P_F_matched=pressure[2],**result))
        if i<steps:
            state,observation,change=box.step(state,observation,dt)
            pressure+=change
    d={k:np.array([r[k] for r in records]) for k in records[0]}
    volume=d['a']**3
    scale=volume[0]*d['rho'][0]
    momentum=d['a']**4*d['momentum']
    orbitals=state.orbitals[0].reshape(len(box.weight),2*box.spatial.d,box.spatial.d)
    summary=dict(
        configuration=asdict(cfg),uv_gap_transition=transition,radial_order=order,radial_nodes=len(box.weight),dt=dt,elapsed=duration,
        max_pressure_work_error=float(np.max(abs((volume*d['rho']+d['work_P'])/scale-1))),
        max_comoving_momentum_error=float(np.max(abs(momentum-momentum[0]))/scale),
        orthonormality_error=float(np.max(abs(adjoint(orbitals) @ orbitals-np.eye(box.spatial.d)))),
        max_fermion_anisotropic_stress=float(np.max(abs(d['P_F_x']-d['P_F_perp']))),
        max_scalar_gradient_fraction=float(np.max(d['scalar_gradient']/d['rho'])),
        min_scalar_inertia=float(np.min(d['A'])),
    )
    for suffix in ('','_matched'):
        instantaneous=volume*d['Jv'+suffix]
        work=cumulative_trapezoid(instantaneous,d['t'],initial=0.)
        stress=volume*d['rho_F'+suffix]
        stress=stress-stress[0]+d['work_P_F'+suffix]
        d['comoving_work'+suffix]=instantaneous
        d['integrated_transfer'+suffix]=work
        summary.update({
            'integrated_transfer'+suffix:float(work[-1]),
            'integrated_absolute_transfer'+suffix:float(trapezoid(abs(instantaneous),d['t'])),
            'transfer_from_stress'+suffix:float(stress[-1]),
            'max_fermion_work_error'+suffix:float(np.max(abs(work-stress))/scale),
            'peak_instantaneous_work'+suffix:float(np.max(abs(instantaneous))),
        })
    return summary,d


def run_case(task):
    name,work=task
    paired=name in PAIRED_CASES
    settings=PREPARATION_CASES[name] if name in PREPARATION_CASES else (PAIRED_CASES if paired else CASES)[name]
    n,qt,order,points,band,steps,*extra=settings
    transition=extra[0] if extra else None
    cfg=replace(Configuration(),longitudinal=n,transverse_cut=qt,
                points=points,scalar_modes=band)
    summary,history=evaluate(cfg,order,steps,paired=paired,transition=transition)
    (work/f'{name}.json').write_text(json.dumps(summary,allow_nan=False)+'\n')
    np.savez_compressed(work/f'{name}.npz',**history)
    print(name,'matched work',summary['integrated_transfer_matched'],
          'budget',summary['max_pressure_work_error'],flush=True)


def gather(work,names):
    series={}
    for name in names:
        with np.load(work/f'{name}.npz') as f:
            series[name]=(json.loads((work/f'{name}.json').read_text()),dict(f))
    comparisons={}
    if 'uv_reference' in names:
        pairs={'wide':('uv_reference','uv_wide'),'later':('uv_reference','uv_later'),
               'wide_later':('uv_wide','uv_wide_later'),'wide_half':('uv_wide','uv_wide_half'),
               'cutoff_later':('uv_later','uv_wide_later')}
    else:
        reference='paired_reference' if 'paired_reference' in names else 'reference'
        pairs={name:(reference,name) for name in series if name!=reference}
    for name,(a,b) in pairs.items():
        comparisons[name]=compare(series[a],series[b])
    return dict(scope='synthetic one-axis expanding slab, one physical mass',
                cases={key:value[0] for key,value in series.items()},
                matched_comparisons=comparisons,scientific_acceptance=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--only',nargs='+',choices=[*CASES,*PAIRED_CASES,*PREPARATION_CASES])
    p.add_argument('--paired',action='store_true')
    p.add_argument('--preparation',action='store_true')
    p.add_argument('--gather',action='store_true')
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    names=args.only or (PREPARATION_CASES if args.preparation else PAIRED_CASES if args.paired else CASES)
    if not args.gather:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(run_case,[(name,args.work) for name in names]))
        if args.only:
            return
    result=gather(args.work,names)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
