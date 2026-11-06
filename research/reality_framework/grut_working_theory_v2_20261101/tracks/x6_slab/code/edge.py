#!/usr/bin/env python3
"""Separate the missing-neighbor cutoff response from the static interior response."""
import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np

from slab import Configuration, Slab


def harmonic(box, order, a, finite_difference):
    cfg=box.cfg
    k=order*2*math.pi/cfg.length
    e=box.energy(a)
    lo,hi=e[:,:-order],e[:,order:]
    dot=cfg.mass**2+(box.perp[:,None]**2+box.k[:-order]*box.k[order:])/a**2
    raw=(1/e).sum(axis=1)-((1+dot/(lo*hi))/(lo+hi)).sum(axis=1)
    weight=cfg.mass**2*box.weight/(cfg.phase_scale**2*a*k**2)
    total=weight @ raw
    edge=weight @ (0.5*(1/e[:,:order]+1/e[:,-order:]).sum(axis=1))
    interior=weight @ ((k/a)**2/(2*lo*hi*(lo+hi))).sum(axis=1)
    assert abs((edge+interior)/total-1)<5e-13
    error=None
    if finite_difference:
        eps=0.002
        wave=np.cos(k*box.x)

        def vacuum(amplitude):
            h=box.full_hamiltonian(a,amplitude*wave/cfg.phase_scale)
            energies=np.linalg.eigvalsh(h)
            return 0.5*(box.weight @ energies[:,:box.d].sum(axis=1))/a**3

        curvature=(vacuum(eps)+vacuum(-eps)-2*vacuum(0))/eps**2
        measured=curvature*2*a*a/k**2
        error=float(abs(measured/total-1))
        assert error<2e-5
    return dict(harmonic=order,coefficient=float(total),edge=float(edge),interior=float(interior),
                edge_fraction=float(edge/total),static_curvature_relative_error=error)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    tests=[]
    for n in (6,9,12,24,48,96):
        cfg=replace(Configuration(),longitudinal=n,transverse_cut=4*n/3,transverse_bins=round(8*n/3),
                    points=max(49,4*n+9))
        box=Slab(cfg)
        rows=[harmonic(box,m,1.0,n in (6,9)) for m in (1,2,3)]
        tests.append(dict(longitudinal_cutoff=n,transverse_cutoff=cfg.transverse_cut,harmonics=rows))
    last=tests[-1]['harmonics']
    ratios=[r['edge']/last[0]['edge'] for r in last]
    assert abs(ratios[1]-0.5)<0.02 and abs(ratios[2]-1/3)<0.02
    result=dict(scope='finite Galerkin static susceptibility, independent of time evolution',
                cylinders=tests,largest_cylinder_edge_ratios=ratios)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Largest-cylinder edge ratios',ratios,flush=True)


if __name__=='__main__':
    main()
