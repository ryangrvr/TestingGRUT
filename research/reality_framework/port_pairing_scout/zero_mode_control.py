#!/usr/bin/env python3
"""Exact N-fermion zero-mode degeneracy and low-energy thermal counting.

For c disconnected cycles, Dirac incidence rank = N-c, leaving 2c zero
one-body modes. Filling the N lowest energies at fixed N fermions gives
binom(2c,c) degenerate zero-mode fillings. This counts ONLY the fermionic
GROUND MANIFOLD on each geometry, not the full excited fermion Gibbs trace.
"""
import json
import math
from pathlib import Path
import numpy as np

def delta(n,t=1.):
    return 4*t*n/math.pi-2*t/math.tan(math.pi/(2*n))

def ground_manifold_connected(n,beta):
    w=np.array([0.]+[math.exp(-beta*delta(l)) for l in range(1,n+1)])
    # Q[k,c] is the coefficient of z^k u^c in
    # exp[(u/2)*sum_l w_l z^l/l].
    Q=np.zeros((n+1,n+1));Q[0,0]=1
    for k in range(1,n+1):
        for c in range(1,k+1):
            Q[k,c]=sum(.5*w[l]*Q[k-l,c-1]
                       for l in range(1,k+1))/k
    Z=sum(math.comb(2*c,c)*Q[n,c] for c in range(1,n+1))
    return w[n]/(n*Z)

def sanity_check(n):
    """Independent direct one-body Dirac zero-mode count for cycle graph."""
    # Each component c contributes one Laplacian zero and TWO Dirac zeros.
    groups=((n,),tuple([1]*n),tuple([2]*(n//2)+([1] if n%2 else [])))
    result=[]
    for p in groups:
        B=np.zeros((n,n),float)
        col=0; offset=0
        for length in p:
            for j in range(length):
                u=offset+j;v=offset+(j+1)%length
                B[u,col]+=1;B[v,col]-=1
                col+=1
            offset+=length
        singular=np.linalg.svd(B,compute_uv=False)
        rank=sum(singular>1e-10)
        count_zero=2*(n-rank)
        c=len(p)
        assert count_zero==2*c
        result.append(dict(cycles=list(p),zero_onebody_modes=count_zero,
                           ground_fermion_fillings=math.comb(2*c,c)))
    return result

def main():
    checks={str(n):sanity_check(n) for n in (3,4,5,6)}
    rows=[]
    for n in (4,8,16,32,64,128):
        rows.append(dict(N=n,probabilities={
            str(beta):ground_manifold_connected(n,beta)
            for beta in (0.,1.,5.,20.)}))
    assert abs(rows[0]['probabilities']['0.0']-.14328358208955222)<1e-12
    assert abs(rows[-1]['probabilities']['1.0']-.0033198923656691287)<1e-12
    assert rows[-1]['probabilities']['5.0']<.05
    result=dict(
        assumptions="Massless incidence fermions, N fixed fermions, conditional ground manifold only",
        exact_degeneracy="binomial(2c,c) zero-mode fillings for c cycles",
        first_lowest_modes=checks,counts=rows,
        caveat="NOT full thermal trace, no nonzero-energy fermion excitations and no coherent quantum graph reconnection")
    Path(__file__).with_name("zero_mode_results.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print("SCOUT_RG02_ZERO_MODE_DEGENERACY_PASS")

if __name__=="__main__":main()
