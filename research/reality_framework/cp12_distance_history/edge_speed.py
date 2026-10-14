"""Exact graph-pool cross-sector norm and a geometry-change speed bound.

Finite incidence matrices certify the counting identity in selected cases.
Large-N rows only evaluate proved formulas, not quantum time evolutions.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,sqrt,sin,pi
import json
import numpy as np
from scipy.sparse import coo_matrix

P=Path(__file__).resolve().parent

def incidence(M,Q,L):
    assert 1<=L and 1<=Q<=M-L
    outside=list(range(L,M));left=list(combinations(outside,Q))
    right=[(e,c) for e in range(L) for c in combinations(outside,Q-1)]
    ids={(e,c):j for j,(e,c) in enumerate(right)}
    rr=[];cc=[]
    for j,c in enumerate(left):
        for e in range(L):
            for f in c:rr.append(ids[e,tuple(x for x in c if x!=f)]);cc.append(j)
    C=coo_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc)),shape=(len(right),len(left))).tocsr()
    out_degree=Q*L;in_degree=M-L-Q+1
    assert np.all(C@np.ones(len(left),np.int64)==in_degree)
    assert np.all(C.T@np.ones(len(right),np.int64)==out_degree)
    assert len(left)*out_degree==len(right)*in_degree
    return dict(M=M,Q=Q,L=L,zero_pool_graphs=len(left),one_pool_graphs=len(right),cross_entries=C.nnz,
                source_degree=out_degree,target_degree=in_degree,exact_cross_norm_square=str(F(out_degree*in_degree,100*M*M)))

def main():
    controls=[]
    for n in [4,5,6,7]:
        M=n*(n-1)//2;q=n-1
        controls.append(incidence(M,q,1));controls.append(incidence(M,q,n-1))
    controls.append(incidence(21,14,6))
    scaling=[]
    for n in [8,16,32,64,128,256,512,1024]:
        M=n*(n-1)//2;q=n-1;L=n-1
        one=F(q*(M-q),100*M*M);pool=F(q*L*(M-L-q+1),100*M*M)
        scaling.append(dict(N=n,Q=q,single_edge_cross_norm_square=str(one),star_pool_cross_norm_square=str(pool),
                            single_edge_N_scaled=float(n*one),star_pool_norm=float(pool)**.5))
    history=[]
    for f in sorted(P.glob('N*_Q*_B*_W*_T*.json')):
        a=json.loads(f.read_text());n=a['N'];q=a['Q'];M=n*(n-1)//2
        g=sqrt(q*(M-q))/(10*M);tau=a['lag'];bound=sin(min(g*tau,pi/2))**2
        classes=a['classes'];Qm=np.array(a['geometry_transitions']);idx=classes.index(1)
        tests=[]
        for i,r in enumerate(classes):
            pr=a['past'][i]['probability']
            changed=Qm[i].sum()-Qm[i,idx] if r==1 else Qm[i,idx]
            assert changed<=pr*bound+1e-11
            tests.append(dict(initial_distance=r,measured_weight=float(changed),conditional_speed_upper=bound,weighted_upper=pr*bound))
        history.append(dict(file=f.name,tests=tests))
    out=dict(status='CP12_POOL_SPEED_EXACT_COUNTS_PASS',incidence_controls=controls,asymptotic_formula_evaluations=scaling,history_speed_checks=history,
             theorem='g_L^2 = Gamma^2 Q L (M-L-Q+1)/M^2 for a nonempty zero-pool subspace; one-edge L=1. Uniform neutral strata attain the Schur upper bound, including in the unbounded physical Hilbert space.',
             limit='Q~nu*N, L~lambda*N: g_L -> Gamma*sqrt(2*nu*lambda); L=1: g_1^2~2*nu*Gamma^2/N.',
             limit_scope='Operator capacity, not actual long-time response, a uniform Taylor extrapolation, or a no-go for state-dependent emergent locality.')
    (P/'EDGE_SPEED.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],finite_controls=len(controls),history_records=len(history))),flush=True)

if __name__=='__main__':main()
