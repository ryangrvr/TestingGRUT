#!/usr/bin/env python3
"""Stoichiometric invariants, not a simulated cosmological population.

Comoving counts: cores Q=0..12, free X+, free X-, neutral scalars.
Core energies and all reaction kernels remain supplied data, not inferred here.
"""
from fractions import Fraction as F
import json
from moments import RES, dump


def main():
    labels=[f'core_{q}' for q in range(13)]+['Xplus','Xminus','neutral']
    charge=list(range(13))+[1,-1,0]
    reactions=[]
    for q in range(2,13):
        v=[0]*16;v[q]-=1;v[q-2]+=1;v[13]+=1;v[14]-=1
        reactions.append({'name':f'forward_{q}','nu':v})
        reactions.append({'name':f'reverse_{q-2}','nu':[-n for n in v]})
    v=[0]*16;v[13]=v[14]=1;v[15]=-1
    reactions.append({'name':'neutral_to_pair','nu':v})
    for r in reactions:r['charge_change']=sum(a*b for a,b in zip(charge,r['nu']))
    assert all(r['charge_change']==0 for r in reactions)
    # Exact rest-energy example: event energy is fixed using a supplied mass table.
    # Values below are synthetic unit tests, NOT the M2 soliton spectrum.
    core={q:F(4,5)*q+F(1,20)*q*q for q in range(13)}
    cases=[]
    for q in range(2,13):
        ein=F(8);eout=ein+core[q]-core[q-2]
        mismatch=core[q]+ein-core[q-2]-eout
        assert mismatch==0
        cases.append({'Q':q,'incoming':float(ein),'outgoing':float(eout),
                      'core_energy_lost':float(core[q]-core[q-2]),'energy_error':float(mismatch)})
    # Same net charge admits different object populations. No abundance follows.
    examples=[{'core_charge':F(0),'free_plus':F(1),'free_minus':F(0)},
              {'core_charge':F(1),'free_plus':F(0),'free_minus':F(0)},
              {'core_charge':F(3,2),'free_plus':F(0),'free_minus':F(1,2)}]
    for row in examples:assert row['core_charge']+row['free_plus']-row['free_minus']==1
    dump(RES/'population_invariants.json',dict(labels=labels,charge_vector=charge,reactions=reactions,
        exact_energy_examples=cases,net_charge_one_examples=[{k:float(v) for k,v in x.items()} for x in examples],
        scope='23 charge-neutral reaction vectors; synthetic energy identities only',
        cosmological_population_computed=False,formation_kernel_supplied=False,
        missing=['formation spectrum','core lifetime/transport','forward and reverse kernels',
                 'neutral reservoir decay spectrum','coherent-to-particle matching']))
    print('23 reaction charge invariants and 11 exact event-energy balances pass; no population inferred.')


if __name__=='__main__':main()
