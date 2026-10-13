#!/usr/bin/env python3
"""Exact graph-diagonal statistics of CP12's coherent matching trial.

Conditional on k locked matching edges, the remaining Q-k graph edges are
uniform among M-k. The script independently enumerates tiny cases and checks
hypergeometric first and second moments. Does NOT identify the ground state.
"""
from itertools import combinations
from math import comb,exp,lgamma

def avoid(M,A,R):
    if M-A<R:return 0.0
    return exp(lgamma(M-A+1)-lgamma(M-A-R+1)-lgamma(M+1)+lgamma(M-R+1))

def moments(n,k):
    assert 0<=k<n/2
    M=comb(n,2);Q=n;R=Q-k;pool=M-k
    u=n-2*k
    p1=avoid(pool,n-1,R)
    p2=avoid(pool,2*n-3,R)
    mean=u*p1
    variance=u*p1*(1-p1)+u*(u-1)*(p2-p1*p1)
    assert variance>=-0.01
    return dict(N=n,k=k,guaranteed_unlocked_vertices=u,
                isolated_vertex_probability=p1,
                expected_isolates_lower_bound=mean,
                fraction_isolated_lower_mean=mean/n,
                variance=max(0.0,variance),
                chebyshev_probability_below_half_mean_upper=min(1.,max(0.,variance)*4/(mean*mean)))

def exhaustive(n,k):
    edges=list(combinations(range(n),2))
    forced=set(edges[:k]);endpoints={v for e in forced for v in e}
    U=list(set(range(n))-endpoints)[:n-2*k]
    rest=[e for e in edges if e not in forced]
    total=first=second=0
    for chosen in combinations(rest,n-k):
        G=set(chosen)|forced
        I=sum(all(v not in e for e in G) for v in U)
        total+=1;first+=I;second+=I*(I-1)
    p=moments(n,k)
    assert abs(first/total-p["expected_isolates_lower_bound"])<1e-10
    want=len(U)*(len(U)-1)*avoid(len(rest),2*n-3,n-k)
    assert abs(second/total-want)<1e-10

if __name__=="__main__":
    for n,k in ((6,1),(7,1),(7,2)):exhaustive(n,k)
    for n in (70,140,700,1400,7000):print(moments(n,n//7))
    print("asymptotic_density_lower",5/7*exp(-12/7))
    print("CP12_COHERENT_TRIAL_ISOLATES_EXACT_PASS")
