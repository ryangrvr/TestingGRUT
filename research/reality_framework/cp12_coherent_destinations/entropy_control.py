"""Independent finite-matrix controls of the unregulated entropy argument.

The proof is regulator-free. These finite-cutoff matrices test its identities,
not its thermodynamic conclusion. Reference beta is a mathematical parameter
in the nonequilibrium bound, not a claim that the evolving state is Gibbs.
"""
from pathlib import Path
import hashlib
import itertools
import json
import math
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.special import logsumexp

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
assert hashlib.sha256((P/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
from frozen_model import build


def binary(p,q):
    value=0.
    if p:value+=p*math.log(p/q)
    if p<1:value+=(1-p)*math.log((1-p)/(1-q))
    return value


def finite(n,Q,L=1):
    H,states,lookup,_=build(n,Q,L)
    H=H.toarray();M=n*(n-1)//2;G=math.comb(M,Q);beta=1.
    A=np.array([Q/2+2*(bits&((1 << n)-1)).bit_count()+sum(x*x for x in E)/8 for g,bits,E in states])
    B=H-np.diag(A)
    b_theory=2*Q+2*.1*Q*(M-Q)/M
    b_actual=float(np.max(abs(np.linalg.eigvalsh(B))))
    assert b_actual<=b_theory+1e-10
    w,V=eigh(H)
    logZA=float(logsumexp(-beta*A));logZH=float(logsumexp(-beta*w))
    assert abs(logZH-logZA)<=beta*b_theory+1e-10
    thermal=np.exp(-beta*w-logZH)
    diagonal=abs(V)**2@thermal
    E=float(thermal@w);EB=E-float(diagonal@A)
    D=-beta*EB+logZA-logZH
    assert -1e-12<=D<=2*beta*b_theory+1e-10
    ref=np.exp(-beta*A-logZA)
    host=min(g for g,_,_ in states);target=np.array([g==host for g,_,_ in states],float)
    p=float(diagonal@target);q=float(ref@target)
    assert binary(p,q)<=D+1e-10
    assert p*math.log(1/q)-math.log(2)<=D+1e-10
    sea=sum(1 << (n+e) for e in range(M));zero=(0,)*(2*M)
    neutral=[lookup[g,sea,zero] for g in sorted(set(s[0] for s in states))]
    W0=np.eye(len(states))[:,neutral]
    pure=W0.sum(axis=1)/math.sqrt(G)
    Emix=float(np.trace(W0.T@H@W0)/G);Epure=float(pure@H@pure)
    assert abs(Emix-(Q/2+.1*Q*(M-Q)/M))<1e-12
    assert abs(Epure-Q/2)<1e-12
    Flog=math.lgamma(n+Q+1)-math.lgamma(Q+1)-math.lgamma(n+1)+2*Q*math.log(1+math.sqrt(8*math.pi))
    assert logZA<=math.log(G)-beta*Q/2+Flog+1e-10
    h_bound=2*Q
    exact_mixed_bound=beta*(Emix+h_bound)+logZA-math.log(G)
    general_mixed_bound=beta*(Emix-Q/2+h_bound)+Flog
    rows=[]
    for T in [0.,.4,2.,8.]:
        phase=np.exp(-1j*T*w)
        mixed=V@(phase[:,None]*(V.T@W0))
        coherent=V@(phase*(V.T@pure))
        pmix=np.mean(abs(mixed)**2,axis=1);ppure=abs(coherent)**2
        Amix=float(pmix@A);Apure=float(ppure@A)
        Dmix=beta*Amix+logZA-math.log(G);Dpure=beta*Apure+logZA
        assert Dmix<=exact_mixed_bound+1e-10<=general_mixed_bound+1e-10
        assert binary(float(pmix@target),q)<=Dmix+1e-10
        assert binary(float(ppure@target),q)<=Dpure+1e-10
        assert abs(np.trace(mixed.conj().T@H@mixed).real/G-Emix)<1e-10
        assert abs(np.vdot(coherent,H@coherent).real-Epure)<1e-10
        rows.append(dict(T=T,mixed_relative_entropy=Dmix,pure_relative_entropy=Dpure,
                         mixed_graph_probability=float(pmix@target),pure_graph_probability=float(ppure@target),
                         mixed_data_processing_gap=Dmix-binary(float(pmix@target),q)))
        field=A-Q/2-2*np.array([(bits&((1 << n)-1)).bit_count() for g,bits,E in states])
        for K in [.05,.25,1.]:
            high=field>K*n;low_rank=int(np.count_nonzero((target==1)&~high))
            measured_high=float(pmix@high)
            assert measured_high<=float(pmix@field)/(K*n)+1e-10
            assert float(pmix@target)<=low_rank/G+measured_high+1e-10
            assert float(pmix@field)<=Emix-Q/2+h_bound+1e-10
    assert abs(rows[0]['pure_relative_entropy']-rows[0]['mixed_relative_entropy']-math.log(G))<1e-10
    return dict(N=n,Q=Q,L=L,dimension=len(states),graph_count=G,
                B_norm=b_actual,B_norm_upper=b_theory,thermal_relative_entropy=D,
                reference_graph_probability=q,thermal_graph_probability=p,
                partition_log_difference=logZH-logZA,
                mixed_entropy=math.log(G),pure_entropy=0.,
                initial_mixed_energy=Emix,initial_pure_energy=Epure,
                exact_all_time_mixed_entropy_upper=exact_mixed_bound,
                general_all_time_mixed_entropy_upper=general_mixed_bound,rows=rows)


def noncommuting_countercontrol():
    A=np.diag([0.,20.]);B=np.array([[0.,.2],[.2,0.]])
    w,V=eigh(A+B);weights=np.exp(-w-logsumexp(-w));rho=(V*weights)@V.T
    q=float(np.exp(-20.)/(1+np.exp(-20.)));p=float(rho[1,1])
    incorrect=float(np.exp(2*np.linalg.norm(B,2))*q)
    assert p>incorrect*1000
    D=float(np.trace(rho@A))+float(logsumexp(-np.diag(A)))-float(-np.sum(weights*np.log(weights)))
    assert binary(p,q)<=D+1e-12
    return dict(actual_probability=p,reference_probability=q,invalid_multiplicative_bound=incorrect,
                relative_entropy=D,meaning='Exponential is not operator monotone; the direct proof uses data processing, not an operator Gibbs domination claim.')


def asymptotic(n,nu):
    Q=nu*n;M=n*(n-1)//2;beta=1.
    logG=math.lgamma(M+1)-math.lgamma(Q+1)-math.lgamma(M-Q+1)
    logS=math.lgamma(n+1)
    logF=math.lgamma(n+Q+1)-math.lgamma(Q+1)-math.lgamma(n+1)+2*Q*math.log(1+math.sqrt(8*math.pi))
    denominator=logG-logS-logF
    b=2*Q+.2*Q*(M-Q)/M
    thermal=1. if denominator<=0 else min(1.,(2*beta*b+math.log(2))/denominator)
    # Uniform neutral mixture, isolated unitary dynamics, conserved exact energy.
    mixedD=beta*(.1*Q*(M-Q)/M+b)+logF
    mixed=1. if denominator<=0 else min(1.,(mixedD+math.log(2))/denominator)
    return dict(N=n,nu=nu,reference_rarity_log_lower=denominator,
                canonical_spatial_probability_upper=thermal,
                all_time_neutral_mixture_spatial_probability_upper=mixed,
                pure_minus_mixed_initial_relative_entropy=logG)


def gauge_sector(n,Q,L):
    from invariant_evolution import Invariant
    sub=Invariant(n,Q,L);H=sub.matrix().toarray()
    sea=sum(1 << (n+e) for e in range(sub.M))
    neutral=[i for i,(g,bits,E) in enumerate(sub.states) if bits==sea and not any(E)]
    K=len(neutral);g_bound=.1*Q*(sub.M-Q)/sub.M
    assert K>=math.comb(sub.M,Q)/math.factorial(n)
    assert min(H[i,i] for i in neutral)>=Q/2-1e-12
    assert max(H[i,i] for i in neutral)<=Q/2+g_bound+1e-12
    A=np.array([Q/2+2*(bits&((1 << n)-1)).bit_count()+sum(x*x for x in E)/8
                for g,bits,E in sub.states])
    w,V=eigh(H);logZ=float(logsumexp(-w));logZA=float(logsumexp(-A))
    assert logZ>=math.log(K)-Q/2-g_bound-1e-12
    assert logZA>=math.log(K)-Q/2-1e-12
    weights=np.exp(-w-logZ);p=abs(V)**2@weights
    B=H-np.diag(A);D=-float(np.trace((V*weights)@V.T@B))+logZA-logZ
    b=2*Q+2*g_bound
    assert 0<=D+1e-12<=2*b+1e-10
    host=sub.states[neutral[0]][0]
    read=np.array([g==host for g,_,_ in sub.states],float)
    q=float(read@np.exp(-A-logZA));prob=float(read@p)
    assert binary(prob,q)<=D+1e-10
    return dict(N=n,Q=Q,L=L,invariant_dimension=len(sub.states),
                neutral_graph_orbits=K,neutral_diagonal_min=min(H[i,i] for i in neutral),
                neutral_diagonal_max=max(H[i,i] for i in neutral),thermal_relative_entropy=D,
                one_orbit_reference_probability=q,one_orbit_Gibbs_probability=prob,
                scope='A separately declared signed invariant-sector Gibbs ensemble; not the full labelled Gibbs state.')


if __name__=='__main__':
    out=dict(status='CP12_DIRECT_ENTROPY_AND_UNITARY_CONTROL_PASS',
             finite_controls=[finite(3,2),finite(4,2),finite(4,3)],
             noncommuting_adversarial_control=noncommuting_countercontrol(),
             conditional_gauge_sector_controls=[gauge_sector(4,3,3),gauge_sector(5,4,2)],
             asymptotic_formula_controls=[asymptotic(n,nu) for nu in [2,3]
                                           for n in [1000,10000,1000000,10**9]],
             scope='Finite-matrix identity checks and analytic bound evaluation, not a thermodynamic simulation.')
    (P/'entropy_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],finite_dimensions=[r['dimension'] for r in out['finite_controls']]),indent=2))
