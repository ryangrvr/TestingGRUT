"""Complete labelled d-regular graph ensembles and their exact neutral-sector
regular-to-nonregular rewiring norm. Pure graph calculation controls the FULL
physical off-block operator because locked matter/field fibers are submatrices.
"""
from pathlib import Path
import itertools,math,json
from collections import Counter
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh,LinearOperator
P=Path(__file__).resolve().parent
GAMMA=.1

def regular_masks(n,d):
 if d>n-1-d:
  full=(1<<(n*(n-1)//2))-1
  return [full^x for x in regular_masks(n,n-1-d)]
 pairs=list(itertools.combinations(range(n),2));edgeid={e:i for i,e in enumerate(pairs)};degrees=[d]*n;out=[]
 def rec(i,mask):
  if i==n:
   if not any(degrees):out.append(mask)
   return
  possible=[j for j in range(i+1,n) if degrees[j]>0]
  k=degrees[i]
  if k<0 or k>len(possible):return
  for dest in itertools.combinations(possible,k):
   degrees[i]=0;newmask=mask
   for j in dest:degrees[j]-=1;newmask|=1<<edgeid[i,j]
   if all(0<=degrees[j]<=n-i-2 for j in range(i+1,n)):rec(i+1,newmask)
   for j in dest:degrees[j]+=1
   degrees[i]=k
 rec(0,0)
 assert len(out)==len(set(out))
 return out

def run(n,d):
 assert n*d%2==0
 q=n*d//2;M=n*(n-1)//2;regs=regular_masks(n,d);rset=set(regs);children={};rr=[];cc=[];parents=Counter()
 for i,mask in enumerate(regs):
  occupied=[e for e in range(M) if mask>>e&1];absent=[e for e in range(M) if not mask>>e&1]
  for e in occupied:
   for f in absent:
    new=mask^(1<<e)^(1<<f);assert new not in rset
    if new not in children:children[new]=len(children)
    rr.append(i);cc.append(children[new]);parents[new]+=1
 bound=GAMMA/M*math.sqrt((d+1)*q*(M-q))
 B=coo_matrix((np.full(len(rr),GAMMA/M),(rr,cc)),shape=(len(regs),len(children))).tocsr()
 assert max(parents.values())<=d+1
 if len(regs)==1:actual=float(np.linalg.norm(B.data))
 else:
  op=LinearOperator((len(regs),len(regs)),matvec=lambda x:B@(B.T@x),dtype=float)
  eigen=float(eigsh(op,k=1,which='LA',tol=1e-12,v0=np.ones(len(regs)),return_eigenvectors=False)[0]);actual=math.sqrt(eigen)
 assert actual<=bound+1e-11
 out=dict(N=n,d=d,Q=q,M=M,regular_graph_count=len(regs),adjacent_nonregular_count=len(children),transitions=len(rr),max_nonregular_predecessors=max(parents.values()),exact_neutral_coupling_norm=actual,universal_physical_bound=bound,bound_per_vertex=bound/n,old_full_offdiag_norm=GAMMA*q*(M-q)/M,scope='Exact norm of neutral graph matrix, equal to norm of the full physical regular/off-regular coupling. No sector ground-energy ordering.')
 (P/f'regular_bound_N{n}_d{d}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 for n,d in [(4,2),(5,2),(6,2),(7,2),(8,2),(6,4),(7,4)]:run(n,d)
