#!/usr/bin/env python3
"""Independent CP12 fixed-graph first-rewire audit. Stdlib only."""
import itertools, json, math
from pathlib import Path
from collections import Counter
GAMMA=.1
def gauss(n,mask,bits,E,edges):
 qv=[-((bits>>i)&1) for i in range(n)];qe=[]
 for e,(i,j) in enumerate(edges):
  a,b=E[2*e:2*e+2];qv[i]+=a;qv[j]+=b
  qe.append(-a-b-((bits>>(n+e))&1)+1)
  if not ((mask>>e)&1):
   assert a==b==0 and ((bits>>(n+e))&1)==1
 return qv,qe
def run(n,radius):
 assert n>=2*radius+1 and radius>=2
 edges=list(itertools.combinations(range(n),2))
 mid={e:i for i,e in enumerate(edges)}
 present={mid[tuple(sorted((i,(i+1)%n)))] for i in range(n)}
 M=len(edges);Q=n;mask=sum(1<<e for e in present)
 bits=sum(1<<(n+e) for e in range(M));E=[0]*(2*M)
 assert gauss(n,mask,bits,E,edges)==([0]*n,[0]*M)
 absent=set(range(M))-present
 def dist(e):
  u,v=edges[e];diff=abs(u-v)
  return min(diff,n-diff)
 nearby={f for f in absent if 2<=dist(f)<=radius}
 assert len(nearby)==n*(radius-1)
 assert len(absent)==n*(n-3)//2
 stats=[]
 for locked in (0,1):
  e_locked=min(present) if locked else None
  EE=E.copy();bb=bits
  if locked:
   source_vertex=edges[e_locked][0]
   bb^=1<<(n+e_locked);bb|=1<<source_vertex;EE[2*e_locked]=1
  assert gauss(n,mask,bb,EE,edges)==([0]*n,[0]*M)
  removable=[e for e in present if e!=e_locked]
  counts=Counter();seen=set()
  for e in removable:
   for f in absent:
    mask2=mask^(1<<e)^(1<<f)
    if n<=12 or (e==removable[0] and f==min(absent)):
     assert gauss(n,mask2,bb,EE,edges)==([0]*n,[0]*M)
    counts[f]+=1
    if n<=12:seen.add(mask2)
  if n<=12:assert len(seen)==len(removable)*len(absent)
  assert set(counts)==absent and set(counts.values())=={len(removable)}
  frac=sum(counts[f] for f in nearby)/sum(counts.values())
  assert math.isclose(frac,2*(radius-1)/(n-3),abs_tol=1e-14)
  stats.append(dict(locked_connections=locked,removable=len(removable),
    distinct_one_swap_states=len(removable)*len(absent),
    destination_multiplicity=min(counts.values()),
    per_absent_edge_probability_T2_coefficient=len(removable)*(GAMMA/M)**2,
    total_one_swap_probability_T2_coefficient=len(removable)*len(absent)*(GAMMA/M)**2,
    fraction_of_first_swaps_with_destination_at_distance_at_most_R=frac))
 return dict(N=n,M=M,Q=Q,R=radius,total_absent=len(absent),nearby_absent=len(nearby),
             local_destination_fraction=len(nearby)/len(absent),checks=stats)
def main():
 data=[run(n,r) for n,r in [(7,2),(9,3),(12,3),(25,3),(101,3)]]
 out=dict(status="CP12_DESTINATION_UNIFORMITY_EXACT_COUNT_PASS",gamma=GAMMA,
  note="Fixed-N early-time T^2 coefficients; one definite physical graph state, not arbitrary coherent graph preparation.",
  data=data)
 Path(__file__).with_name("results.json").write_text(json.dumps(out,indent=2)+"\n")
 print(json.dumps([dict(N=x["N"],R=x["R"],local_fraction=x["local_destination_fraction"],
   neutral_departure=x["checks"][0]["total_one_swap_probability_T2_coefficient"],
   charged_departure=x["checks"][1]["total_one_swap_probability_T2_coefficient"]) for x in data],indent=2))
 print(out["status"])
if __name__=="__main__":main()
