#!/usr/bin/env python3
"""Exact domain-wall constraint and finite-system formation consequence."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from calculate import BASE, hamiltonian, state, evolve

HERE=Path(__file__).resolve().parent


def main():
    target=HERE/'algebra_results.json'
    if target.exists():raise SystemExit('Refusing to overwrite algebra record')
    out={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runs':[]}
    for L in [3,4,5]:
        H0,z=hamiltonian(L,{**BASE,'epsilon':0.})
        D=sum((1-z[L+i]*z[L+(i+1)%L])/2 for i in range(L))
        comm=H0*(D[None,:]-D[:,None])
        psi=state(L,'formation')
        zero=(D==0)
        initial=float(np.sum(abs(psi[zero])**2))
        H,_=hamiltonian(L,BASE)
        E,V=eigh(H)
        ts=np.linspace(0,40,201)
        trajectory=evolve(E,V,psi,ts)
        probability=np.sum(abs(trajectory[zero])**2,axis=0)
        out['runs'].append(dict(L=L,hard_constraint_commutator_max=float(np.max(abs(comm))),
            initial_ordered_sector_probability=initial,expected_initial_probability=2.**(1-L),
            soft_constraint_times=ts.tolist(),soft_constraint_ordered_probability=probability.tolist(),
            maximum_ordered_probability=float(probability.max()),
            final_ordered_probability=float(probability[-1])))
    out['checks']={'exact_domain_wall_conservation_at_epsilon0':all(x['hard_constraint_commutator_max']==0 for x in out['runs']),
                   'initial_ordered_sector_weight':all(abs(x['initial_ordered_sector_probability']-x['expected_initial_probability'])<1e-14 for x in out['runs'])}
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks':out['checks'],'ordered_probabilities':[{k:x[k] for k in ['L','initial_ordered_sector_probability','maximum_ordered_probability','final_ordered_probability']} for x in out['runs']]}))


if __name__=='__main__':main()
