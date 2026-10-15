"""Complete electric-budget closure with integer symmetry acceleration.

Identical transition law and initial ray to cp12_unmasked_formation/run.py.
Array-backed sparse construction avoids Python tuple storage per matrix entry.
"""
from __future__ import annotations
from array import array
from collections import deque
from itertools import combinations,permutations
from math import comb,sqrt,lcm
from pathlib import Path
import hashlib,json,os,pickle,sys,time
import numpy as np
from scipy.sparse import coo_matrix
from fast_symmetry import FastInvariant

P=Path(__file__).resolve().parent
RF=P.parent
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,gauss


def model(n,q,budget,source=False):
    cache=Path(os.environ['GRUT_DENSITY_CACHE']) if os.environ.get('GRUT_DENSITY_CACHE') else None
    cache_name=f'n{n}_q{q}_b{budget}_s{int(source)}_fast_v1.pkl'
    if cache and (cache/cache_name).exists():
        with (cache/cache_name).open('rb') as f:return pickle.load(f)
    started=time.monotonic()
    sub=FastInvariant.__new__(FastInvariant)
    sub.n=n;sub.Q=q;sub.es=list(combinations(range(n),2));sub.M=len(sub.es)
    group=[(0,)+p for p in permutations(range(1,n))] if source else permutations(range(n))
    sub.group=[sub.permutation(p) for p in group]
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)}
    sub.inverse=[]
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
    sc=lcm(8,10*sub.M)
    rr=array('i');cc=array('i');vv=array('d')
    br=array('i');bc=array('i');bv=array('d');outer={}
    last_report=started
    while queue:
        col=queue.popleft();s=states[col];gauss(n,s)
        for target,c in transitions(n,s,sc):
            inside=sum(x*x for x in target[2])<=budget
            key,sign,size=add(target) if inside else sub.canonical(target)
            if not sign:continue
            value=c/sc*sign*sqrt(sizes[col]/size)
            if inside:
                rr.append(index[key]);cc.append(col);vv.append(value)
            else:
                if key not in outer:outer[key]=len(outer)
                br.append(outer[key]);bc.append(col);bv.append(value)
        if time.monotonic()-last_report>15:
            print(json.dumps(dict(stage='closure',N=n,Q=q,budget=budget,processed=col,states=len(states),queued=len(queue),seconds=time.monotonic()-started)),flush=True)
            last_report=time.monotonic()
    dim=len(states)
    H=coo_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),shape=(dim,dim)).tocsr()
    B=coo_matrix((np.frombuffer(bv,dtype=np.float64),(np.frombuffer(br,dtype=np.int32),np.frombuffer(bc,dtype=np.int32))),shape=(len(outer),dim)).tocsr()
    assert np.linalg.norm((H-H.T).data)<1e-9
    psi=np.asarray([sqrt(sizes[i]/comb(sub.M,q)) if bits==sea and not any(E) else 0. for i,(g,bits,E) in enumerate(states)],complex)
    assert abs(np.linalg.norm(psi)-1)<1e-12
    assert abs(np.vdot(psi,H@psi).real-q/2)<1e-10
    assert abs(np.linalg.norm(H@psi)**2-q*q/4-2*q)<1e-9
    sub.cache.clear()
    result=(sub,H,B,states,np.asarray(sizes),psi)
    if cache:
        cache.mkdir(parents=True,exist_ok=True)
        with (cache/cache_name).open('wb') as f:pickle.dump(result,f)
    print(json.dumps(dict(stage='matrix',N=n,Q=q,budget=budget,dimension=dim,outside=B.shape[0],nnz=H.nnz,graphs=len(reps),seconds=time.monotonic()-started)),flush=True)
    return result
