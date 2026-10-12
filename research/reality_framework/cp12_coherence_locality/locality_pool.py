"""Exact rewiring norms for a calibrated sparse neighborhood pool.

F is supplied independently; |F|=S. No spatial embedding is inferred.
k counts present bonds outside F. Johnson-graph symmetry gives exact
neutral norms, which also equal the complete physical operator norms.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
import itertools
import json
import math
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, LinearOperator, expm_multiply

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parent/"cp12_remainder_and_mixing"))
from variational_exclusion import regular_lower, frac_dec, sqrt_upper
from near_regular_exclusion import conservative_trial_upper
sys.path.insert(0, str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def quotient(n, q, S, K):
    M = n*(n-1)//2
    L = M-S
    assert q <= S <= M and 0 <= K <= min(q,L)
    diag = [(q-k)*(S-q+k)+k*(L-k) for k in range(K+1)]
    c2 = [(q-k)*(L-k)*(k+1)*(S-q+k+1) for k in range(K)]
    boundary2 = (q-K)*(L-K)*(K+1)*(S-q+K+1) if K < min(q,L) else 0
    return diag, c2, boundary2


def neutral_control(n, q, S):
    M = n*(n-1)//2
    masks = [sum(1 << e for e in es) for es in itertools.combinations(range(M),q)]
    index = {mask:i for i,mask in enumerate(masks)}
    pool = (1 << S)-1
    long_count = np.array([(mask & ~pool).bit_count() for mask in masks])
    rr,cc = [],[]
    for i,mask in enumerate(masks):
        for e in range(M):
            if not mask >> e & 1: continue
            for f in range(M):
                if mask >> f & 1: continue
                rr.append(i); cc.append(index[mask ^ (1 << e) ^ (1 << f)])
    adjacency = coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(masks),len(masks))).tocsr()
    rows=[]
    for K in range(min(q,M-S)+1):
        ids=np.flatnonzero(long_count <= K)
        outside=np.flatnonzero(long_count > K)
        A=adjacency[ids][:,ids].toarray()
        B=adjacency[ids][:,outside].toarray()
        normA=float(max(abs(np.linalg.eigvalsh(A))))
        normB=float(np.linalg.norm(B,2)) if B.size else 0.
        diag,c2,b2=quotient(n,q,S,K)
        radial=np.diag(diag)+np.diag(np.sqrt(c2),1)+np.diag(np.sqrt(c2),-1)
        radial_norm=float(np.linalg.eigvalsh(radial)[-1])
        assert abs(normA-radial_norm)<1e-9 and abs(normB**2-b2)<1e-8
        rows.append(dict(K=K,sector_graphs=len(ids),internal_adjacency_norm=normA,
                         quotient_norm=radial_norm,crossing_adjacency_norm_square=normB**2,
                         exact_boundary_square=b2))
    return dict(N=n,Q=q,S=S,graph_count=len(masks),rows=rows)


def physical_control(n=4,q=3,Lcut=3):
    # A cycle supplies a declared pool and genuine coherent initial geometry.
    pairs=list(itertools.combinations(range(n),2))
    pool=sum(1 << pairs.index(tuple(sorted((i,(i+1)%n)))) for i in range(n))
    S=pool.bit_count(); M=len(pairs); outside=M-S
    H,states,_,_=build(n,q,Lcut)
    masks=np.array([s[0] for s in states])
    counts=np.array([int(m & ~pool).bit_count() for m in masks])
    coo=H.tocoo();changing=masks[coo.row] != masks[coo.col]
    V=coo_matrix((coo.data[changing],(coo.row[changing],coo.col[changing])),shape=H.shape).tocsr()
    norms=[]
    for K in range(min(q,outside)+1):
        ids=np.flatnonzero(counts<=K); other=np.flatnonzero(counts>K)
        A=V[ids][:,ids];B=V[ids][:,other]
        actualA=float(max(abs(eigsh(A,k=2,which='BE',return_eigenvectors=False,tol=1e-12))))
        if B.shape[1]:
            op=LinearOperator((len(ids),len(ids)),matvec=lambda x:B @ (B.T @ x),dtype=float)
            actualB=float(eigsh(op,k=1,which='LA',return_eigenvectors=False,tol=1e-12)[0])
        else:actualB=0.
        diag,c2,b2=quotient(n,q,S,K)
        radial=np.diag(diag)+np.diag(np.sqrt(c2),1)+np.diag(np.sqrt(c2),-1)
        expectedA=.1/M*float(np.linalg.eigvalsh(radial)[-1])
        expectedB=.01/M**2*b2
        assert abs(actualA-expectedA)<1e-10 and abs(actualB-expectedB)<1e-10
        norms.append(dict(K=K,actual_internal_norm=actualA,expected_internal_norm=expectedA,
                          actual_crossing_norm_square=actualB,expected_crossing_norm_square=expectedB))
    sea=sum(1 << (n+e) for e in range(M))
    start=np.array([i for i,(mask,bits,E) in enumerate(states)
                    if mask & ~pool == 0 and bits==sea and not any(E)])
    psi=np.zeros(len(states));psi[start]=1/math.sqrt(len(start))
    g2=.01*q*outside*S/M**2
    early2=.01/M**2*q*outside*(S-q+1)
    Hpsi=H @ psi
    mean_energy=float(psi @ Hpsi)
    variance=float(Hpsi @ Hpsi)-mean_energy**2
    predicted_mean=.5*q+.1/M*q*outside
    predicted_variance=2*q+early2
    assert abs(mean_energy-predicted_mean)<1e-10 and abs(variance-predicted_variance)<1e-10
    time_rows=[]
    for T in [.001,.01,.1,.4,1.,2.,5.]:
        evolved=expm_multiply(-1j*T*H,psi)
        mean=float(counts @ abs(evolved)**2)
        assert mean <= g2*T*T+1e-10
        if T <= .001: assert abs(mean/T**2-early2)<1e-6
        time_rows.append(dict(T=T,mean_long_bonds=mean,proven_mean_bound=g2*T*T))
    return dict(N=n,Q=q,S=S,L=Lcut,dimension=len(states),initial_neutral_graphs=len(start),
                norm_controls=norms,early_time_mean_coefficient=early2,
                dynamical_coefficient_bound=g2,time_controls=time_rows,
                prepared_mean_energy=mean_energy,predicted_mean_energy=predicted_mean,
                prepared_energy_variance=variance,predicted_energy_variance=predicted_variance,
                scope="Complete matter/gauge/graph dynamics. Prepared pool is supplied; density stability is not a light-cone theorem.")


def certificate(n, fraction):
    q=2*n;S=4*n;M=n*(n-1)//2;L=M-S
    K=(fraction*n).numerator//(fraction*n).denominator
    diag,c2,b2=quotient(n,q,S,K)
    maxdiag=F(max(diag),10*M)
    maxc2=F(max(c2,default=0),100*M*M)
    upper_trial=conservative_trial_upper(n)
    with localcontext() as ctx:
        ctx.prec=80;ctx.rounding=ROUND_CEILING
        r_upper=frac_dec(maxdiag,ROUND_CEILING)+2*sqrt_upper(maxc2)
        # sqrt_upper(0) gives a tiny positive outward padding, also safe.
        with localcontext() as low:
            low.prec=80;low.rounding=ROUND_FLOOR
            gap=Decimal(n)*regular_lower(n)-r_upper-upper_trial
            gap2=gap*gap if gap>0 else Decimal(0)
        probability=min(Decimal(1),frac_dec(F(b2,100*M*M),ROUND_CEILING)/gap2) if gap>0 else Decimal(1)
        # All k sectors, not only this local subspace, control the mean.
        _,all_c2,_=quotient(n,q,S,min(q,L))
        cross_upper=2*sqrt_upper(F(max(all_c2),100*M*M))
        short_upper=frac_dec(F(S*S,40*M),ROUND_CEILING)
        with localcontext() as low:
            low.prec=80;low.rounding=ROUND_FLOOR
            record_cost_lower=Decimal(n)*regular_lower(n)-upper_trial
            mean_long_lower=max(Decimal(0),(record_cost_lower-short_upper-cross_upper)/Decimal('.1'))
            normalized_coherence_lower=record_cost_lower/frac_dec(F(q*(M-q),10*M),ROUND_CEILING)
    return dict(N=n,Q=q,S=S,K=K,threshold_long_density=str(fraction),
                internal_norm_upper=str(r_upper),local_sector_gap_lower=str(gap),
                local_sector_ground_probability_upper=str(probability),
                total_graph_record_energy_lower=str(record_cost_lower),
                normalized_rewiring_coherence_lower=str(normalized_coherence_lower),
                mean_long_bonds_lower=str(mean_long_lower),
                scope="Outward full-H bounds; F is an independently specified pool, not an emergent frame.")


if __name__ == "__main__":
    neutral=[neutral_control(4,3,4),neutral_control(5,5,7),neutral_control(6,12,12)]
    physical=physical_control()
    certs=[certificate(n,f) for n in [128,512,2048,8192,32768]
           for f in [F(1,10),F(1,4),F(1,2)]]
    out=dict(status="CALIBRATED_POOL_AND_COHERENCE_CHECKS_PASS",neutral_norm_controls=neutral,
             full_physical_norm_and_dynamics=physical,certificates=certs)
    (P/"locality_pool_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
