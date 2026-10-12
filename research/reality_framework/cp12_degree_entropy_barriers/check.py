"""Exact labelled simple-graph counting. Classical degree-only comparator."""
from pathlib import Path
from itertools import combinations
from collections import Counter
from math import comb,factorial,exp,sqrt,pi
import json
P=Path(__file__).resolve().parent
R=[1,0,0]
for n in range(3,501):
 R.append(sum(comb(n-1,l-1)*(factorial(l-1)//2)*R[n-l] for l in range(3,n+1)))
rows=[]
for n in range(4,8):
 edges=list(combinations(range(n),2));lookup={e:i for i,e in enumerate(edges)};regular=[];hist=Counter()
 for occ in combinations(range(len(edges)),n):
  deg=[0]*n
  for e in occ:
   a,b=edges[e];deg[a]+=1;deg[b]+=1
  W=sum(d*(d-1)//2 for d in deg);hist[W]+=1
  if all(d==2 for d in deg):regular.append(sum(1<<e for e in occ))
 assert len(regular)==R[n]
 predecessors=Counter()
 for mask in regular:
  count=0
  for i,(u,v) in enumerate(edges):
   if not mask>>i&1:continue
   for a,b in [(u,v),(v,u)]:
    for c in range(n):
     if c==a:continue
     j=lookup[tuple(sorted((a,c)))]
     if mask>>j&1:continue
     new=(mask^(1<<i))|(1<<j);predecessors[new]+=1;count+=1
  assert count==2*n*(n-3)
 assert max(predecessors.values())<=3
 probs=[]
 for x in [0.,1.,5.,10.]:
  prob=R[n]/sum(count*exp(-x*(W-n)) for W,count in hist.items())
  bound=1/(1+2*n*(n-3)*exp(-x)/3)
  assert prob<=bound+1e-12
  probs.append(dict(beta_lambda=x,regular_probability=prob,bound=bound))
 rows.append(dict(N=n,regular_count=R[n],max_regular_predecessors=max(predecessors.values()),probabilities=probs))
ring=[]
for n in [6,10,20,50,100,500]:
 p=(factorial(n-1)//2)/R[n];asym=exp(.75)*sqrt(pi)/(2*sqrt(n))
 ring.append(dict(N=n,conditional_ring_probability=p,asymptotic=asym,ratio=p/asym))
assert abs(ring[-1]['ratio']-1)<.003
out=dict(status='DEGREE_AND_COMPONENT_ENTROPY_CONTROLS_PASS',enumeration=rows,conditional_ring=ring,scope='Classical labelled degree-only model; not full CP12 Gibbs state.')
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
