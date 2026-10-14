#!/usr/bin/env python3
"""Independent CP12 polynomial all-time graph-assembly exclusion controls.

Stdlib only. Checks combinatorial parent loop-packing bound against exact
N=6 graph enumeration and evaluates the electric Chernoff / min-entropy terms.
Does NOT diagonalize CP12 or prove an experimental spacetime result.
"""
from itertools import combinations
from math import comb, exp, factorial, lgamma, log, pi, sqrt

def logbinom(n,k):
    return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)

def loop_A(n,q,R):
    m=comb(n,2)
    return sum((n*q/m)**r/(2*r) for r in range(3,R+1))

def exact_count():
    n=6;q=6;R=3
    edges=list(combinations(range(n),2))
    for k in (1,2):
        successes=0
        for subset in combinations(edges,q):
            G=set(subset)
            tri=[set(v) for v in combinations(range(n),3)
                 if all(pair in G for pair in combinations(v,2))]
            if (k==1 and len(tri)>0) or (k==2 and any(a.isdisjoint(b) for a,b in combinations(tri,2))):
                successes+=1
        exact=successes/comb(comb(n,2),q)
        bound=min(1.,loop_A(n,q,R)**k/factorial(k))
        assert exact<=bound+1e-12
        print('exact packing',n,q,R,k,exact,'parent_upper',bound)

def theta_small(s):
    assert 0<s<1
    # Poisson summation, truncation below double-precision error for these s.
    return sqrt(pi/s)*(1+2*sum(exp(-pi*pi*j*j/s) for j in range(1,8)))

def terms(n,nu,alpha,delta,R=4,eps=.25,hop=1.,gamma=.1):
    assert 0<delta<alpha/nu and nu*n<comb(n,2)
    q=nu*n; k=int(alpha*n+0.999999)
    log_fraction=min(0.,k*log(loop_A(n,q,R))-lgamma(k+1))
    log_rotor_upper=2*n/eps+2*q*log(theta_small(n**(-delta)))
    loglow=log_fraction+logbinom(n+q,q)+log_rotor_upper
    high=(2*hop+gamma)*q/(n**(1+delta))
    return loglow/(n*log(n)), high

def main():
    exact_count()
    for nu in (2,3):
        alpha=.25; delta=alpha/(2*nu)
        for n in (100,1000,1000000,10**12,10**18):
            lo,hi=terms(n,nu,alpha,delta)
            print('nu',nu,'N',n,'delta',delta,'loglow/(NlogN)',round(lo,6),'high_tail_upper',round(hi,6))
        assert abs((nu*delta)-alpha/2)<1e-12
    print('CP12_PREPARATION_MINENTROPY_FORMULA_CHECK_PASS')

if __name__=='__main__': main()
