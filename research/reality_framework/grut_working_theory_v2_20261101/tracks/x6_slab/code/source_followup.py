#!/usr/bin/env python3
"""Test finer cutoff and preparation, and diagnose the sharp-cylinder edge."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from slab import Configuration, SIGMA_Z, Slab, adjoint
from source import REFERENCE, compare, evaluate

EXTRA = {
    'wider_first': ({'longitudinal':12,'transverse_cut':16.0,'transverse_bins':32,'points':65},1200,1),
    'transverse_64_first': ({'transverse_bins':64},600,1),
    'reference_second': ({},1200,2),
    'wide_second': ({'longitudinal':9,'transverse_cut':12.0,'transverse_bins':24},1200,2),
}


def edge(cfg,a=1.0):
    box=Slab(cfg)
    e=box.energy(a)
    dk=box.k[1]-box.k[0]
    pref=cfg.mass**2/(cfg.phase_scale**2*a*dk**2)
    end=pref*(box.weight @ (1/e[:,0]+1/e[:,-1]))/2
    lo,hi=e[:,:-1],e[:,1:]
    interior=pref*dk**2/a**2*(box.weight @ (1/(2*lo*hi*(lo+hi))).sum(axis=1))
    combined=box.spatial_match(a)
    assert abs((end+interior)/combined-1)<5e-13
    return dict(longitudinal=cfg.longitudinal,transverse_cut=cfg.transverse_cut,
                edge=float(end),interior=float(interior),total=float(combined),
                edge_fraction=float(end/combined))


def rotating_vacuum():
    cfg=Configuration(longitudinal=2,transverse_cut=3,transverse_bins=3,points=25,
                      scalar_modes=2,expand=False,scalar_mass=0)
    box=Slab(cfg)
    gamma=np.kron(SIGMA_Z,np.eye(2))
    controls=[]
    for speed in (0.10,0.05,0.025):
        target=box.free_blocks(1)+speed/(2*cfg.phase_scale)*gamma
        _,vectors=np.linalg.eigh(target)
        c_target=vectors[...,:,:2] @ adjoint(vectors[...,:,:2])
        errors=[]
        for order in (1,2):
            state=box.initial(np.zeros(cfg.points),np.full(cfg.points,speed),adiabatic=order)
            blocks=box.blocks(state.occupied)
            diff=blocks @ adjoint(blocks)-c_target
            errors.append(float(np.max(abs(diff))))
        controls.append(dict(speed=speed,first_error=errors[0],second_error=errors[1]))
    assert all(row['second_error']<row['first_error'] for row in controls)
    assert controls[-1]['second_error']<controls[0]['second_error']/30
    return controls


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prior_work',type=Path)
    parser.add_argument('prior_results',type=Path)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    original=json.loads(args.prior_results.read_text())['cases']
    inherited={}
    for case in ('wide_half_step','half_step','transverse_bins'):
        with np.load(args.prior_work/f'{case}.npz') as data:
            inherited[case]=(original[case],dict(data))

    def run(name):
        changes,steps,order=EXTRA[name]
        summary,data=evaluate(replace(REFERENCE,**changes),steps,order)
        np.savez_compressed(args.work/f'{name}.npz',**data)
        print(name,summary['integrated_transfer_matched'],flush=True)
        return name,(summary,data)

    vacuum=rotating_vacuum()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        cases=dict(pool.map(run,EXTRA))
    all_cases={**inherited,**cases}
    pairs=(('wide_half_step','wider_first'),('transverse_bins','transverse_64_first'),
           ('reference_second','wide_second'),('half_step','reference_second'),
           ('wide_half_step','wide_second'))
    comparisons={f'{a}__{b}':dict(raw=compare(all_cases[a],all_cases[b]),
                                   matched=compare(all_cases[a],all_cases[b],'_matched'))
                 for a,b in pairs}
    edges=[edge(REFERENCE),edge(replace(REFERENCE,longitudinal=9,transverse_cut=12.0,transverse_bins=24)),
           edge(replace(REFERENCE,**EXTRA['wider_first'][0]))]
    payload=dict(
        cases={k:v[0] for k,v in cases.items()},comparisons=comparisons,
        rotating_vacuum=vacuum,spatial_subtraction_decomposition=edges,
        scientific_acceptance=False,
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    for key,pair in comparisons.items():
        print(key,'matched instantaneous',pair['matched']['instantaneous_peak_scaled_error'],flush=True)


if __name__=='__main__':
    main()
