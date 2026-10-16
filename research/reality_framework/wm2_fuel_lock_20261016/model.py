"""WM2-S2: reversible head-coupled energy storage, no supplied adjacency."""
from __future__ import annotations
from itertools import permutations
import math
import numpy as np
from scipy import sparse


def words(n):
    return [(0,)+tail for m in range(1,n+1)
        for tail in permutations(range(1,n),m-1)]


def spent(w):return sum(1<<v for v in w)


def mismatch(n,w,q):
    # Root is bound/spent; every other bit is retained in the component.
    mask=spent(w)
    return (mask^q).bit_count()


def diagonal(n,w,q,U=16.,penalty=64.,delta=0.):
    f=n-len(w);nq=q.bit_count();m=mismatch(n,w,q)
    return U*(f+nq-n)+penalty*m+delta*(nq-1)


def neighbors(n,w,q,p,eps=.25):
    """Return target, symbolic move kind, coefficient without 1/sqrt(N)."""
    if len(w)<n:
        u=w[-1]
        for v in range(n):
            if v in w:continue
            if not(q>>v&1):yield (w+(v,),q|(1<<v),p),'stored_append',1.
            if p is not None and p in (u,v):
                yield (w+(v,),q,p),'mismatch_append',eps
    if len(w)>1:
        v=w[-1];u=w[-2]
        if q>>v&1:yield (w[:-1],q^(1<<v),p),'stored_pop',1.
        if p is not None and p in (u,v):
            yield (w[:-1],q,p),'mismatch_pop',eps
    if p is not None and p in w:
        i=w.index(p)
        if i:yield (w,q,w[i-1]),'matter',1.
        if i+1<len(w):yield (w,q,w[i+1]),'matter',1.


def matrices(n,U=16.,penalty=64.,g=1.,eps=.25,J=.25,delta=0.,probe=True):
    # The whole invariant root component is enumerated, including all q.
    # In probe vacuum mismatch coupling vanishes, so only consistent q is reachable.
    states=[(w,q,p) for w in words(n)
        for q in (range(1,1<<n,2) if probe else [spent(w)])
        for p in (range(n) if probe else [None])]
    index={s:i for i,s in enumerate(states)}
    rr=[];cc=[];vv=[]
    for col,(w,q,p) in enumerate(states):
        rr.append(col);cc.append(col);vv.append(diagonal(n,w,q,U,penalty,delta))
        for target,kind,amp in neighbors(n,w,q,p,eps):
            value=-J if kind=='matter' else -(g if kind.startswith('stored') else amp)/math.sqrt(n)
            rr.append(index[target]);cc.append(col);vv.append(value)
    h=sparse.coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr();h.eliminate_zeros()
    pi=np.array([i for i,(w,q,p) in enumerate(states) if q==spent(w)],dtype=int)
    qi=np.array([i for i,(w,q,p) in enumerate(states) if q!=spent(w)],dtype=int)
    return dict(H=h,states=states,index=index,P=pi,Q=qi)


def constants(U=16.,penalty=64.,g=1.,eps=.25,J=.25):
    t=2*g+2*eps+2*J;aP=2*g+2*J
    a=penalty-U-aP;gap=a-t;assert gap>0
    ratio=t/a
    return dict(U=U,Lambda=penalty,g=g,epsilon=eps,J=J,
        diagonal_mismatch_min=penalty-U,
        all_offdiagonal_norm_upper=t,P_block_norm_upper=aP,
        protected_spectral_gap_lower=gap,PQ_norm_upper=eps,
        spectral_projector_distance_upper=eps/gap,
        Sigma_C0=eps*eps/gap,Sigma_ratio=ratio,
        Sigma_distance_power='d',
        Sigma_exponential_C=eps*eps/gap,
        xi=1/math.log(1/ratio),near_packet_lower=J-eps*eps/gap)


def radial(n,g=1.,delta=0.):
    # Order m=1,...,N. This is exact in probe vacuum.
    b=-g*np.sqrt(np.arange(n-1,0,-1)/n)
    return np.diag(delta*np.arange(n))+np.diag(b,1)+np.diag(b,-1)


def formation_constants(n,g=1.):
    harmonic=sum(1/k for k in range(1,n))
    spacing=math.pi*g/math.sqrt(n*(n+.5))
    T=32*harmonic/spacing
    L=max(2,math.ceil(n/2));rank=n-L+1
    return dict(N=n,threshold=L,diagonal_ensemble_macroscopic=rank/n,
        diagonal_ensemble_full=1/n,spacing_lower=spacing,
        observation_window=T,average_error_upper=.125,
        finite_window_macroscopic_lower=rank/n-.125,
        uniform_375_lower=.375,
        time_fraction_probability_at_least_quarter_lower=1/6)


def detuning_bound(n,delta,g=1.):
    L=max(2,math.ceil(n/2));sigma=g*math.sqrt((n-1)/n)
    b=g*math.sqrt((n-L+1)/n);den=abs(delta)*(L-1)-2*g
    return dict(N=n,delta=delta,threshold=L,
        all_time_macroscopic_upper=1. if den<=0 else min(1.,((sigma+b)/den)**2))


def register_encode(n,w,q,p):
    return (frozenset(zip(w[:-1],w[1:])),w[-1],frozenset(set(range(n))-set(w)),q,p)


def register_neighbors(n,c,eps=.25):
    bonds,head,free,q,p=c
    for u in range(n):
        for v in range(n):
            if u==v:continue
            if head==u and v in free and (u,v) not in bonds:
                new=(bonds|{(u,v)},v,free-{v})
                if not(q>>v&1):yield (*new,q|(1<<v),p),'stored_append',1.
                if p is not None and p in (u,v):yield (*new,q,p),'mismatch_append',eps
            if head==v and v not in free and (u,v) in bonds:
                new=(bonds-{(u,v)},u,free|{v})
                if q>>v&1:yield (*new,q^(1<<v),p),'stored_pop',1.
                if p is not None and p in (u,v):yield (*new,q,p),'mismatch_pop',eps
    if p is not None:
        for u,v in bonds:
            if p==u:yield (bonds,head,free,q,v),'matter',1.
            elif p==v:yield (bonds,head,free,q,u),'matter',1.
