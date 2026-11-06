#!/usr/bin/env python3
"""Evolve retained quantum bath directions with passive channel occupations."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from handoff import (BOUNDARIES, CASES, OBSERVABLES, adjoint, boundary_spread, comparison, evolve,
                     ledger, make_box, measure)
from retention import CASES_TESTED, adaptive

HOLDOUT="bath_holdout"
METHODS=("energy","energy_plus_one","energy_plus_two")
CASES[HOLDOUT]=(1,2,13,8.5,0.01)
HYBRID_CASES=(*CASES_TESTED,HOLDOUT)


def cut_model(slab,projected,basis,layout):
    n=slab.geometry.bath_bins
    f=np.empty_like(slab.bath_frequency)
    g=np.empty_like(slab.g)
    transform=basis.copy()
    omitted=0.
    for ell in range(slab.l):
        for i in range(slab.n):
            for h in (0,1):
                for b in (0,1):
                    begin=(2*i+h)*slab.block+b*slab.nb+1
                    for sign in (0,1):
                        section=slice(begin+sign*n,begin+(sign+1)*n)
                        leg=slice(sign*n,(sign+1)*n)
                        r=layout[ell,i,h,b,sign]
                        U=transform[ell,section,section]
                        raw=slab.frequency[ell,i,h,leg]
                        H=(U.T*raw) @ U
                        e,A=np.linalg.eigh(H[:r,:r])
                        U[:,:r]=U[:,:r] @ A
                        f[ell,i,h,b,leg]=np.r_[e,np.diag(H)[r:]]
                        v=U.T @ slab.g[ell,i,h,b,leg]
                        v[r:]=0.
                        g[ell,i,h,b,leg]=v
                        omitted=max(omitted,float(np.linalg.norm(H[:r,r:])))
    other=copy(slab)
    other.g=g
    other.bath_frequency=f
    other.kinetic_blocks=slab.kinetic_blocks.copy()
    for b in (0,1):
        off=b*slab.nb
        bath=np.arange(off+1,off+slab.nb)
        other.kinetic_blocks[...,bath,bath]=f[...,b,:]
        other.kinetic_blocks[...,off,bath]=g[...,b,:]
        other.kinetic_blocks[...,bath,off]=g[...,b,:]
    rotated=replace(projected,frame=adjoint(transform) @ projected.frame)
    return other,rotated,omitted


def run(args):
    name,folder=args
    slab=make_box(CASES[name])
    dt=CASES[name][-1]
    ref,checkpoints=evolve(slab,slab.prepare(dt),0.,dt,checkpoints=BOUNDARIES)
    tests={}
    for t,state in sorted(checkpoints.items()):
        row={}
        for extra,label in enumerate(METHODS):
            projected,info,basis,layout=adaptive(slab,state,extra,return_layout=True)
            reduced,rotated,strength=cut_model(slab,projected,basis,layout)
            old=measure(slab,projected,slab.observe(projected,source=True))
            new=measure(reduced,rotated,reduced.observe(rotated,source=True))
            trace,_=evolve(reduced,rotated,t,dt)
            row[label]=dict(**info,**comparison(ref,trace),**ledger(trace),
                            omitted_coupling_norm=strength,
                            model_switch_difference={k:float(new[k]-old[k])
                                                     for k in ('rho','P','work','momentum',*OBSERVABLES)},
                            terminal={k:float(trace[k][-1]) for k in OBSERVABLES})
        tests[str(t)]=row
        print(name,t,{k:max(v['terminal_relative_errors'].values()) for k,v in row.items()},flush=True)
    spreads={key:boundary_spread(tests,key) for key in METHODS}
    result=dict(boundaries=tests,reference_terminal={k:float(ref[k][-1]) for k in OBSERVABLES},
                reference_ledger=ledger(ref),boundary_spreads=spreads,
                passes_all_boundaries={key:max(spreads[key].values())<0.01 and all(
                    r[key]['jump_pass'] and r[key]['terminal_pass'] and r[key]['energy_pressure_error']<1e-5
                    for r in tests.values()) for key in METHODS})
    (folder/f'{name}.json').write_text(json.dumps(result,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run,[(name,args.work) for name in HYBRID_CASES]))
    result=dict(scope='coherent retained bath directions with their boundary response to passive occupations removed',
                cases={k:json.loads((args.work/f'{k}.json').read_text()) for k in HYBRID_CASES},
                scientific_acceptance=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
