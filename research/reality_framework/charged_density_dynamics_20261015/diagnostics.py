"""Dimension-neutral physical circulation, including charged edge holes."""
from __future__ import annotations
from itertools import combinations
import numpy as np


def cycle_edges(n,state,balanced_only=False):
    g,bits,E=state
    arcs=[];reach=[1<<v for v in range(n)]
    for e,(a,b) in enumerate(combinations(range(n),2)):
        if not (g>>e&1):continue
        ea,eb=E[2*e:2*e+2]
        if ea*eb>=0 or (balanced_only and ea!=-eb):continue
        u,v=(a,b) if ea>0 else (b,a)
        arcs.append((e,u,v));reach[u]|=1<<v
    for k in range(n):
        for u in range(n):
            if reach[u]>>k&1:reach[u]|=reach[k]
    return [(e,u,v) for e,u,v in arcs if reach[v]>>u&1]


def state_diagnostics(n,q,state):
    g,bits,E=state
    nv=(bits&((1<<n)-1)).bit_count()
    circ=cycle_edges(n,state);balanced=cycle_edges(n,state,True)
    assert len(balanced)<=q-nv
    assert sum(abs(x) for x in E)>=2*len(circ)+nv
    adjacency=[set() for _ in range(n)]
    for e,u,v in circ:adjacency[u].add(v);adjacency[v].add(u)
    unseen={v for v in range(n) if adjacency[v]};covered=len(unseen);sizes=[]
    while unseen:
        root=min(unseen);seen={root};todo=[root]
        for a in todo:
            for b in adjacency[a]-seen:seen.add(b);todo.append(b)
        unseen-=seen;sizes.append(len(seen))
    largest=max(sizes,default=0)
    degree=[0]*n
    for e,(a,b) in enumerate(combinations(range(n),2)):
        if g>>e&1:degree[a]+=1;degree[b]+=1
    return (len(circ),nv,len(balanced),covered,largest,
            len(circ)-covered+len(sizes),sum((d-2*q/n)**2 for d in degree)/n)


def arrays(n,q,states):
    return np.asarray([state_diagnostics(n,q,s) for s in states],float).T


def event(diag,q,n,a=.25,b=.25):
    return (diag[0]>=int(np.ceil(a*q-1e-12))) & (diag[1]>=int(np.ceil(b*n-1e-12)))
