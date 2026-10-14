"""Frozen CP12, exact signed quotient fixing an ordered source/receiver pair.

Only the symmetry subgroup differs from the parent builder. A total electric
square budget is a numerical projection. All outgoing transitions are kept.
"""
from pathlib import Path
from collections import deque
from itertools import combinations,permutations
from math import comb,lcm,sqrt
import hashlib,json,os,pickle,sys,time
import numpy as np
from scipy.sparse import coo_matrix

P=Path(__file__).resolve().parent
G=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((G/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
sys.path.insert(0,str(G))
from invariant_evolution import Invariant
from exact_destination import transitions,gauss


def model(n,q,budget):
    cache=Path(os.environ['CP12_FORMATION_CACHE']) if os.environ.get('CP12_FORMATION_CACHE') else None
    name=f'ordered_n{n}_q{q}_b{budget}_v1.pkl'
    if cache and (cache/name).exists():
        with (cache/name).open('rb') as f:return pickle.load(f)
    sub=Invariant.__new__(Invariant);sub.n=n;sub.Q=q
    sub.es=list(combinations(range(n),2));sub.M=len(sub.es)
    sub.group=[sub.permutation((0,1)+p) for p in permutations(range(2,n))]
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)};sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for a,b in enumerate(p['vertex']):inv[b]=a
        sub.inverse.append(gid[tuple(inv)])
    sub.gcache={};sub.cache={}
    masks=[sum(1<<e for e in c) for c in combinations(range(sub.M),q)]
    reps=sorted({sub.graph_info(g)[0] for g in masks})
    sea=sum(1<<(n+e) for e in range(sub.M));zero=(0,)*(2*sub.M)
    states=[];sizes=[];index={};todo=deque()
    def add(s):
        key,sign,size=sub.canonical(s)
        if sign and key not in index:
            index[key]=len(states);states.append(key);sizes.append(size);todo.append(key)
        return key,sign,size
    for g in reps:add((g,sea,zero))
    sc=lcm(8,10*sub.M);columns=[];outer={};br=[];bc=[];bv=[]
    while todo:
        s=todo.popleft();col=index[s];gauss(n,s);entries=[]
        for target,c in transitions(n,s,sc):
            inside=sum(x*x for x in target[2])<=budget
            key,sign,size=add(target) if inside else sub.canonical(target)
            if not sign:continue
            value=c/sc*sign*sqrt(sizes[col]/size)
            if inside:entries.append((index[key],value))
            else:
                if key not in outer:outer[key]=len(outer)
                br.append(outer[key]);bc.append(col);bv.append(value)
        columns.append(entries)
        if col and col%10000==0:print(json.dumps(dict(stage='closure',N=n,Q=q,budget=budget,processed=col,queued=len(todo))),flush=True)
    rr=[];cc=[];vv=[]
    for j,col in enumerate(columns):
        for i,v in col:rr.append(i);cc.append(j);vv.append(v)
    H=coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr()
    B=coo_matrix((bv,(br,bc)),shape=(len(outer),len(states))).tocsr()
    assert np.linalg.norm((H-H.T).data)<1e-9
    psi=np.array([sqrt(sizes[i]/comb(sub.M,q)) if bits==sea and not any(E) else 0 for i,(g,bits,E) in enumerate(states)],complex)
    assert abs(np.linalg.norm(psi)-1)<1e-12
    assert abs(np.vdot(psi,H@psi).real-q/2)<1e-10
    assert abs(np.linalg.norm(H@psi)**2-q*q/4-2*q)<1e-9
    out=(sub,H,B,states,np.array(sizes),psi)
    print(json.dumps(dict(stage='matrix',N=n,Q=q,budget=budget,dimension=len(states),outside=B.shape[0],nnz=H.nnz,graph_representatives=len(reps))),flush=True)
    if cache:
        cache.mkdir(parents=True,exist_ok=True);sub.cache.clear()
        with (cache/name).open('wb') as f:pickle.dump(out,f)
    return out


def distance(n,states,es):
    values={}
    for g,_,_ in states:
        if g in values:continue
        adj=[set() for _ in range(n)]
        for e,(a,b) in enumerate(es):
            if g>>e&1:adj[a].add(b);adj[b].add(a)
        dist=[n]*n;dist[0]=0;todo=[0]
        for a in todo:
            for b in adj[a]:
                if dist[b]==n:dist[b]=dist[a]+1;todo.append(b)
        values[g]=dist[1]
    return np.array([values[g] for g,_,_ in states],int)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--q',type=int);ap.add_argument('--budget',type=int,required=True);a=ap.parse_args();model(a.n,a.q or a.n-1,a.budget)
