"""WM2-S1: all-label end assembly, no supplied spatial neighbor relation.

Tuples encode the exactly invariant path component of finite-support
register operations. They are not a primitive lattice or an imposed
ordering of event labels. All permutations and all partial prefixes occur.
"""
from __future__ import annotations
from collections import deque
from itertools import permutations
import math
import numpy as np
from scipy import sparse


def words(n):
    return [w for m in range(1,n+1)
        for tail in permutations(range(1,n),m-1) for w in [(0,)+tail]]


def neighbors(n,w,p,virtual=False):
    """Unweighted primitive edges: pop/push and current-bond probe hop."""
    if len(w)>1:yield (w[:-1],p),'assembly'
    if len(w)<n-(1 if virtual else 0):
        used=set(w)
        for u in range(n):
            if u not in used:yield (w+(u,),p),'assembly'
    if p in w:
        i=w.index(p)
        if i:yield (w,w[i-1]),'matter'
        if i+1<len(w):yield (w,w[i+1]),'matter'


def basis(n,probe=True):
    ww=words(n)
    states=[(w,p) for w in ww for p in (range(n) if probe else [None])]
    return states,{s:i for i,s in enumerate(states)}


def matrices(n,U=16.,eps=1.,J=.25,probe=True):
    states,index=basis(n,probe); rr=[];cc=[];geo=[];matter=[]
    for col,(w,p) in enumerate(states):
        if len(w)<n:
            for u in range(n):
                if u not in w:
                    row=index[w+(u,),p]
                    rr.extend([row,col]);cc.extend([col,row]);geo.extend([-eps/math.sqrt(n)]*2);matter.extend([0.,0.])
        if probe:
            if p in w:
                i=w.index(p)
                if i+1<len(w):
                    row=index[w,w[i+1]]
                    rr.extend([row,col]);cc.extend([col,row]);geo.extend([0.,0.]);matter.extend([-J,-J])
    shape=(len(states),len(states))
    hg=sparse.coo_matrix((geo,(rr,cc)),shape=shape).tocsr();hg.eliminate_zeros()
    hm=sparse.coo_matrix((matter,(rr,cc)),shape=shape).tocsr();hm.eliminate_zeros()
    defects=np.array([n-len(w) for w,p in states])
    h=hg+hm+sparse.diags(U*defects,format='csr')
    pi=np.flatnonzero(defects==0); qi=np.flatnonzero(defects>0)
    return dict(H=h,geometry=hg,matter=hm,states=states,index=index,P=pi,Q=qi,defects=defects)


def constants(n,U=16.,eps=1.,J=.25,Emax=None):
    Emax=2*J if Emax is None else Emax
    t=2*eps+2*J; a=U-Emax
    gap=U-2*eps-4*J
    assert a>t and gap>0
    return dict(N=n,U=U,epsilon=eps,J=J,E_max=Emax,
        geometry_kinetic_norm_upper=2*eps,all_offdiagonal_norm_upper=t,
        PQ_norm=eps/math.sqrt(n),Q_energy_lower=U-t,
        protected_spectral_gap_lower=gap,
        spectral_projector_distance_upper=eps/(math.sqrt(n)*gap),
        bare_defect_probability_all_times_upper=4*eps*eps/(n*gap*gap),
        Sigma_prefactor=eps*eps/(n*(a-t)),Sigma_ratio=t/a,
        xi=1/math.log(a/t),nearest_coupling_lower=J,
        minimum_virtual_word='at least d_op+2',
        source_phase_variance=.25)


def formation_bound(n,U=16.,eps=1.):
    E0=U*(n-1);sigma=eps*math.sqrt((n-1)/n);b=eps/math.sqrt(n)
    threshold=max(2,math.ceil(n/2))
    bl=eps*math.sqrt((n-threshold+1)/n)
    return dict(N=n,U=U,epsilon=eps,initial_auxiliary_energy=E0,
        initial_energy_standard_deviation=sigma,
        extensive_length_threshold=threshold,
        extensive_assembly_probability_upper_all_times=min(1.,((sigma+bl)/(U*(threshold-1)-2*eps))**2),
        full_assembly_probability_upper_all_times=min(1.,((sigma+b)/E0)**2))


def config_distance(n,source,virtual=False):
    distance={source:0};queue=deque([source])
    while queue:
        c=queue.popleft()
        for target,kind in neighbors(n,*c,virtual=virtual):
            if virtual and len(target[0])==n:continue
            if target not in distance:
                distance[target]=distance[c]+1;queue.append(target)
    return distance


def register_encode(n,w,p=None):
    """Finite-qubit representation used for an independent primitive check."""
    bonds=frozenset(zip(w[:-1],w[1:]));free=frozenset(set(range(n))-set(w))
    return (bonds,w[-1],free,p)


def register_neighbors(n,c):
    bonds,head,free,p=c
    # X_uv uses only b_uv,h_u,h_v,f_v. No whole-graph degree test.
    for u in range(n):
        for v in range(n):
            if u==v:continue
            if head==u and v in free and (u,v) not in bonds:
                yield (bonds|{(u,v)},v,free-{v},p),'assembly'
            if head==v and v not in free and (u,v) in bonds:
                yield (bonds-{(u,v)},u,free|{v},p),'assembly'
    if p is not None:
        for u,v in bonds:
            if p==u:yield (bonds,head,free,v),'matter'
            elif p==v:yield (bonds,head,free,u),'matter'
