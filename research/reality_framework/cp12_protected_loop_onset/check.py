#!/usr/bin/env python3
"""Gauge-compatible three-event reachability check (not a quantum amplitude calculation)."""
from collections import deque
from itertools import combinations
N=3
edges=list(combinations(range(N),2))
start=((0,)*N,(1,)*len(edges),(0,)*(2*len(edges)))

def neighbors(s):
 a,b,E=s
 for ei,(u,v) in enumerate(edges):
  for leg,vertex in enumerate((u,v)):
   j=2*ei+leg
   if b[ei] and not a[vertex]:
    aa=list(a);bb=list(b);EE=list(E)
    aa[vertex]=1;bb[ei]=0;EE[j]+=1
    yield tuple(aa),tuple(bb),tuple(EE)
   if a[vertex] and not b[ei]:
    aa=list(a);bb=list(b);EE=list(E)
    aa[vertex]=0;bb[ei]=1;EE[j]-=1
    yield tuple(aa),tuple(bb),tuple(EE)

def protected(s):
 a,b,E=s
 return [i for i in range(len(edges)) if E[2*i] or E[2*i+1]]

def fluxonly(s):
 a,b,E=s
 return not any(a) and all(b) and len(protected(s))==3

d={start:0};q=deque([start]);min_core=min_flux=None
while q:
 s=q.popleft();dep=d[s]
 if len(protected(s))==3 and min_core is None:min_core=dep
 if fluxonly(s) and min_flux is None:min_flux=dep
 if dep==6:continue
 for nxt in neighbors(s):
  if nxt not in d:d[nxt]=dep+1;q.append(nxt)
print('3-cycle first protected core at hop order:',min_core)
print('3-cycle first flux-only core at hop order:',min_flux)
print('states reached up to six hops:',len(d))
assert min_core==3 and min_flux==6 and len(d)==44
print('PASS')
