"""Bare physical charge-pair energy versus incidence distance, with independent BFS."""
from pathlib import Path
from collections import deque
import itertools,json
from model import graph_basis,ELECTRIC,MASS,COST

def distances(n,mask,start):
 edges=list(itertools.combinations(range(n),2));A=[[] for _ in range(n+len(edges))]
 for e,(i,j) in enumerate(edges):
  if mask>>e&1:
   for v in [i,j]:A[n+e].append(v);A[v].append(n+e)
 D={start:0};queue=deque([start])
 while queue:
  i=queue.popleft()
  for j in A[i]:
   if j not in D:D[j]=D[i]+1;queue.append(j)
 return D
rows=[];count=0
for n,q in [(3,2),(3,3),(4,3),(4,4)]:
 M=n*(n-1)//2;checked=0
 for es in itertools.combinations(range(M),q):
  mask=sum(1<<e for e in es);minimum={}
  for _,bits,E in graph_basis(n,mask,1):
   occupied=[i for i in range(n) if bits>>i&1]
   if len(occupied)!=1:continue
   holes=[e for e in range(M) if not(bits>>(n+e)&1)]
   assert len(holes)==1
   key=(holes[0],occupied[0]);cost=sum(x*x for x in E)
   minimum[key]=min(minimum.get(key,10**9),cost)
  for e in es:
   D=distances(n,mask,n+e)
   for i in range(n):
    if i in D:
     assert minimum[e,i]==D[i];checked+=1
    else:assert (e,i) not in minimum
 count+=checked;rows.append(dict(N=n,Q=q,connected_charge_pairs_checked=checked))
result=dict(total_pairs=count,controls=rows,formula='E_bare - kappa Q = 2m + epsilon_E d_incidence/2',qualification='Minimum diagonal energy at fixed charge occupations; not the interacting screened potential or a spacetime metric.')
Path(__file__).with_name('distance_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print('CP12_CHARGE_DISTANCE_PASS')
