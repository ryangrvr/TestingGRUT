"""Independent electric-field enumeration and signed-permutation audit."""
from pathlib import Path
import itertools,json,math
import numpy as np
from model import graph_basis,build,gauss,wilson_vector

def brute_basis(n,mask,L):
 edges=list(itertools.combinations(range(n),2));present=[e for e in range(len(edges)) if mask>>e&1];M=len(edges);out=set()
 # All electric assignments; Gauss law independently determines allowed occupations.
 for vals in itertools.product(range(-L,L+1),repeat=2*len(present)):
  vertex=np.zeros(n,int);edgeocc=np.ones(M,int);E=[0]*(2*M)
  for a,e in enumerate(present):
   i,j=edges[e];x,y=vals[2*a:2*a+2];vertex[i]+=x;vertex[j]+=y;edgeocc[e]=1-x-y;E[2*e]=x;E[2*e+1]=y
  if min(vertex,default=0)<0 or max(vertex,default=0)>1 or min(edgeocc,default=0)<0 or max(edgeocc,default=0)>1:continue
  bits=sum(int(v)<<i for i,v in enumerate(np.r_[vertex,edgeocc]));out.add((mask,bits,tuple(E)))
 return out
rows=[]
for n,mask,L in [(3,3,1),(3,7,1),(3,7,2),(4,45,1)]:
 a=set(graph_basis(n,mask,L));b=brute_basis(n,mask,L);assert a==b
 rows.append(dict(N=n,mask=mask,cutoff=L,states=len(a)))

sym=[]
for n,q,L in [(3,2,2),(3,3,2),(4,3,2),(4,4,2)]:
 H,states,lookup,degree=build(n,q,L);edges=list(itertools.combinations(range(n),2));ei={e:i for i,e in enumerate(edges)};M=len(edges)
 perm=list(range(n));perm[0],perm[1]=perm[1],perm[0];mode=list(perm);sg=[1]*n;edge_map=[];leg_map=[]
 for e,(i,j) in enumerate(edges):
  a,b=perm[i],perm[j];f=ei[tuple(sorted((a,b)))];edge_map.append(f);mode.append(n+f);sg.append(1 if a<b else -1)
  leg_map.extend([2*f+(0 if a<b else 1),2*f+(1 if a<b else 0)])
 mapping=[];sign=[]
 for mask,bits,E in states:
  newmask=sum(1<<edge_map[e] for e in range(M) if mask>>e&1)
  occ=[i for i in range(n+M) if bits>>i&1];mapped=[mode[i] for i in occ]
  parity=sum(mapped[i]>mapped[j] for i in range(len(mapped)) for j in range(i+1,len(mapped)))
  newbits=sum(1<<i for i in mapped);newE=[0]*(2*M)
  for i,j in enumerate(leg_map):newE[j]=E[i]
  mapping.append(lookup[newmask,newbits,tuple(newE)]);sign.append((-1)**parity*math.prod(sg[i] for i in occ))
 sign=np.array(sign);difference=H[mapping,:][:,mapping].multiply(sign[:,None]).multiply(sign[None,:])-H
 err=float(np.linalg.norm(difference.data));assert err<1e-12
 # Constraints vanish exactly, and the graph Laplacian row degree matches reversible moves.
 assert all(not np.any(gauss(n,s)) for s in states)
 assert np.all(degree<=q*(M-q))
 sym.append(dict(N=n,Q=q,cutoff=L,dimension=len(states),permutation_error=err))
result=dict(independent_electric_enumeration=rows,signed_permutation_and_constraints=sym)
Path(__file__).with_name('independent_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));print('CP12_INDEPENDENT_CONSTRAINT_CHECKS_PASS')
