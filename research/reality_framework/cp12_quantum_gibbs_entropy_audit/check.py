#!/usr/bin/env python3
"""CP12 regulated quantum Gibbs entropy bound: analytical finite-N verifier.

The inequality assumes the complete frozen CP12 Hamiltonian; this script
checks its combinatorial/dimension evaluation and limiting coefficient.
It does not solve the infinite-cutoff theory or execute many-body evolution.
"""
import math, json
from pathlib import Path
KAPPA=.5; HOP=1.; GAMMA=.1

def lc(n,k):
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def case(N,nu,L,beta=1.):
    M=N*(N-1)//2; Q=nu*N
    assert 0<Q<M and L>=1
    # Upper bound for all labelings of one unlabelled lattice graph.
    logS=math.lgamma(N+1)
    logU=(logS-lc(M,Q)+lc(N+Q,Q)+2*Q*math.log(2*L+1)
          +beta*(2*abs(HOP)*Q+GAMMA*Q*(M-Q)/M))
    return dict(N=N,nu=nu,L=L,beta=beta,
        log_probability_bound=logU,
        nontrivial=logU<0,
        rate=logU/(N*math.log(N)))

def main():
    for N,Q in ((3,2),(4,3),(5,5),(6,6)):
        for L in (1,2):
            # Independent small combinatorial count of active fermion choices
            from itertools import combinations
            count=sum(1 for _ in combinations(range(N+Q),Q))
            assert count==math.comb(N+Q,Q)
            assert count*(2*L+1)**(2*Q)>=count
    rows=[case(N,nu,L) for nu in (1,2,3) for L in (1,2,3)
          for N in (100,1000,10000,100000,1000000)]
    for nu in (2,3):
        sub=[r for r in rows if r['nu']==nu and r['L']==1]
        assert sub[-1]['log_probability_bound']<0
        assert all(sub[i+1]['rate']<sub[i]['rate'] for i in range(len(sub)-1))
    data=dict(status="CP12_QUANTUM_GIBBS_ENTROPY_ANALYTICAL_CHECK_PASS",
      bound="(|S|/C(M,Q))*C(N+Q,Q)*(2L+1)^(2Q)*exp(beta*(2|t|Q+Gamma*Q(M-Q)/M))",
      scope="complete frozen quantum H, fixed Q, finite electric cutoff, finite beta and restricted |S|",
      rows=rows)
    Path(__file__).with_name("results.json").write_text(json.dumps(data,indent=2)+"\n")
    print(data['status'])
    for nu in (1,2,3):
        print("nu",nu, [(r['N'],round(r['rate'],4)) for r in rows if r['nu']==nu and r['L']==1])
if __name__=="__main__":main()
