"""Hermitian, permutation-covariant split of the *frozen* CP-12 hopping.

chi_e is the diagonal indicator that edge e lies on a directed through-flow
cycle. K_in = sum_e {chi_e,K_e}/2; V_out = K-K_in. These are diagnostic
comparison operators. The physical evolution always retains the original H.
"""
from __future__ import annotations
from array import array
from collections import deque
from itertools import combinations, permutations
from math import comb, lcm, sqrt
from pathlib import Path
import json, os, pickle, sys, time
import numpy as np
from scipy.sparse import coo_matrix

HERE=Path(__file__).resolve().parent
RF=HERE.parent
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from fast_symmetry import FastInvariant
from exact_destination import transitions, gauss


def flow(n,state):
    g,bits,E=state
    reach=[1<<i for i in range(n)]; arcs=[]
    for e,(a,b) in enumerate(combinations(range(n),2)):
        x,y=E[2*e:2*e+2]
        if g>>e&1 and x*y<0:
            u,v=(a,b) if x>0 else (b,a)
            arcs.append((e,u,v));reach[u]|=1<<v
    for k in range(n):
        for u in range(n):
            if reach[u]>>k&1:reach[u]|=reach[k]
    mask=sum(1<<e for e,u,v in arcs if reach[v]>>u&1)
    return mask,reach


def hopping_edge(n,s,t):
    if s[0]!=t[0] or s[1]==t[1]:return None
    changed=(s[1]^t[1])>>n
    assert changed.bit_count()==1
    return changed.bit_length()-1


def outside_weight(n,s,t,e,info=None):
    """1-(chi_e(s)+chi_e(t))/2, without redoing a full closure at t.

    A single physical hop changes one half-leg by one unit. A through-flow
    cannot reverse in that step. A return path for the changed edge is thus
    the same before and after the hop, including when that edge first appears.
    """
    mask,reach=flow(n,s) if info is None else info
    a,b=list(combinations(range(n),2))[e]
    x,y=t[2][2*e:2*e+2]
    target=False
    if x*y<0:
        u,v=(a,b) if x>0 else (b,a)
        target=bool(reach[v]>>u&1)
    return 1.-.5*(bool(mask>>e&1)+target)


class Triplets:
    def __init__(self):self.r=array('i');self.c=array('i');self.v=array('d')
    def add(self,r,c,v):
        if v:self.r.append(r);self.c.append(c);self.v.append(v)
    def csr(self,shape):
        return coo_matrix((np.frombuffer(self.v,dtype=np.float64),
            (np.frombuffer(self.r,dtype=np.int32),np.frombuffer(self.c,dtype=np.int32))),shape=shape).tocsr()


def quotient_model(n,q,budget,source=False):
    """Complete reachable budget, all outgoing rows, exact signed quotient."""
    cache=Path(os.environ['GRUT_SHORTCUT_CACHE']) if os.environ.get('GRUT_SHORTCUT_CACHE') else None
    name=f'split_n{n}_q{q}_b{budget}_s{int(source)}_v1.pkl'
    if cache and (cache/name).exists():
        with (cache/name).open('rb') as f:return pickle.load(f)
    started=time.monotonic();sub=FastInvariant.__new__(FastInvariant)
    sub.n=n;sub.Q=q;sub.es=list(combinations(range(n),2));sub.M=len(sub.es)
    group=[(0,)+p for p in permutations(range(1,n))] if source else permutations(range(n))
    sub.group=[sub.permutation(p) for p in group]
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)};sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for a,b in enumerate(p['vertex']):inv[b]=a
        sub.inverse.append(gid[tuple(inv)])
    sub.gcache={};sub.cache={}
    masks=[sum(1<<e for e in es) for es in combinations(range(sub.M),q)]
    reps=sorted({sub.graph_info(g)[0] for g in masks})
    sea=sum(1<<(n+e) for e in range(sub.M));zero=(0,)*(2*sub.M)
    states=[];index={};sizes=[];queue=deque()
    def add(s):
        key,sign,size=sub.canonical(s)
        if sign and key not in index:
            index[key]=len(states);states.append(key);sizes.append(size);queue.append(len(states)-1)
        return key,sign,size
    for g in reps:add((g,sea,zero))
    sc=lcm(8,10*sub.M);hh=Triplets();vv=Triplets();bb=Triplets();bv=Triplets();outer={}
    counts={0:0,.5:0,1:0};last_report=started
    while queue:
        col=queue.popleft();s=states[col];gauss(n,s);info=flow(n,s)
        for target,c in transitions(n,s,sc):
            e=hopping_edge(n,s,target)
            weight=outside_weight(n,s,target,e,info) if e is not None else 0.
            if e is not None:counts[weight]+=1
            inside=sum(x*x for x in target[2])<=budget
            key,sign,size=add(target) if inside else sub.canonical(target)
            if not sign:continue
            value=c/sc*sign*sqrt(sizes[col]/size)
            if inside:
                row=index[key];hh.add(row,col,value);vv.add(row,col,weight*value)
            else:
                row=outer.setdefault(key,len(outer));bb.add(row,col,value);bv.add(row,col,weight*value)
        if time.monotonic()-last_report>15:
            print(json.dumps(dict(stage='split_closure',N=n,Q=q,budget=budget,processed=col,
                states=len(states),queued=len(queue),seconds=time.monotonic()-started)),flush=True)
            last_report=time.monotonic()
    dim=len(states);shape=(dim,dim);bs=(len(outer),dim)
    H=hh.csr(shape);V=vv.csr(shape);B=bb.csr(bs);BV=bv.csr(bs)
    checks=dict(H_hermiticity=float(np.linalg.norm((H-H.T).data)),
        V_hermiticity=float(np.linalg.norm((V-V.T).data)))
    assert max(checks.values())<1e-8
    psi=np.asarray([sqrt(sizes[i]/comb(sub.M,q)) if bits==sea and not any(E) else 0.
        for i,(g,bits,E) in enumerate(states)],complex)
    assert abs(np.linalg.norm(psi)-1)<1e-12
    checks['neutral_mean_energy']=float(np.vdot(psi,H@psi).real)
    checks['neutral_variance']=float(np.linalg.norm(H@psi-q/2*psi)**2)
    checks['neutral_core_generator_residual']=float(np.linalg.norm((H-V)@psi-q/2*psi))
    checks['neutral_outside_hopping_square']=float(np.linalg.norm(V@psi)**2+np.linalg.norm(BV@psi)**2)
    assert abs(checks['neutral_mean_energy']-q/2)<1e-10
    assert abs(checks['neutral_variance']-2*q)<1e-9
    assert checks['neutral_core_generator_residual']<1e-10
    sub.cache.clear()
    result=dict(sub=sub,H=H,V=V,B=B,BV=BV,states=states,sizes=np.asarray(sizes),psi=psi,
        checks=checks,construction_seconds=time.monotonic()-started,
        physical_hop_weight_counts=counts)
    if cache:
        cache.mkdir(parents=True,exist_ok=True)
        with (cache/name).open('wb') as f:pickle.dump(result,f)
    print(json.dumps(dict(stage='split_matrix',N=n,Q=q,budget=budget,dimension=dim,
        outside=B.shape[0],nnz=H.nnz,seconds=time.monotonic()-started,checks=checks)),flush=True)
    return result


def physical_model(n,q,budget):
    """Independent full labelled basis for a recorded-core experiment."""
    from frozen_model import graph_basis
    cache=Path(os.environ['GRUT_SHORTCUT_CACHE']) if os.environ.get('GRUT_SHORTCUT_CACHE') else None
    name=f'physical_n{n}_q{q}_b{budget}_v1.pkl'
    if cache and (cache/name).exists():
        with (cache/name).open('rb') as f:return pickle.load(f)
    started=time.monotonic();M=n*(n-1)//2;sc=lcm(8,10*M)
    gs=[sum(1<<e for e in c) for c in combinations(range(M),q)]
    states=[s for g in gs for s in graph_basis(n,g,int(sqrt(budget))) if sum(x*x for x in s[2])<=budget]
    index={s:i for i,s in enumerate(states)};outer={};hh=Triplets();bb=Triplets()
    edge_hops=[Triplets() for _ in range(M)];edge_outer=[Triplets() for _ in range(M)]
    for col,s in enumerate(states):
        gauss(n,s)
        for t,a in transitions(n,s,sc):
            e=hopping_edge(n,s,t)
            if t in index:
                row=index[t];hh.add(row,col,a/sc)
                if e is not None:edge_hops[e].add(row,col,a/sc)
            else:
                gauss(n,t);row=outer.setdefault(t,len(outer));bb.add(row,col,a/sc)
                if e is not None:edge_outer[e].add(row,col,a/sc)
    shape=(len(states),len(states));bs=(len(outer),len(states))
    H=hh.csr(shape);B=bb.csr(bs)
    Ks=[x.csr(shape) for x in edge_hops];BKs=[x.csr(bs) for x in edge_outer]
    assert np.linalg.norm((H-H.T).data)<1e-11
    assert all(np.linalg.norm((k-k.T).data)<1e-12 for k in Ks)
    sea=sum(1<<(n+e) for e in range(M));psi=np.array([1/sqrt(comb(M,q)) if b==sea and not any(E) else 0.
        for g,b,E in states],complex)
    assert abs(np.linalg.norm(psi)-1)<1e-12
    result=dict(H=H,B=B,Ks=Ks,BKs=BKs,states=states,psi=psi,construction_seconds=time.monotonic()-started)
    if cache:
        cache.mkdir(parents=True,exist_ok=True)
        with (cache/name).open('wb') as f:pickle.dump(result,f)
    print(json.dumps(dict(stage='physical_matrix',N=n,Q=q,budget=budget,dimension=len(states),
        outside=B.shape[0],seconds=time.monotonic()-started)),flush=True)
    return result
