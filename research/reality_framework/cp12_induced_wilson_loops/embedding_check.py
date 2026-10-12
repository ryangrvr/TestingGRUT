"""Hostile embedding control: extra edges, including triangular shortcuts."""
from pathlib import Path
from fractions import Fraction as F
import itertools,sys,json
import numpy as np
from scipy.sparse import diags
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,wilson_vector
from loop_paths import coefficient
rows=[]
for label,n,pairs in [('square_diagonal',4,[(0,1),(1,2),(2,3),(0,3),(0,2)]),('square_pendant',5,[(0,1),(1,2),(2,3),(0,3),(3,4)])]:
 edges=list(itertools.combinations(range(n),2));M=len(edges);q=len(pairs);a=F(M-q,10*M)
 active=[edges.index(tuple(sorted(e))) for e in pairs];mask=sum(1<<e for e in active)
 H,states,lookup,_=build(n,q,2,masks=[mask],rewire=False);V=H-diags(H.diagonal())
 sea=sum(1<<(n+e) for e in range(M));vac=(mask,sea,(0,)*(2*M));i=lookup[vac];j=lookup[wilson_vector(n,[0,1,2,3],vac)]
 charged=np.array([k for k,(g,b,E) in enumerate(states) if b&((1<<n)-1)])
 rem=np.array([sum((b>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for g,b,E in states])
 gap=H.diagonal()+float(a)*rem;gap-=gap[i]
 inv=-1/gap[charged];Vc=V[charged][:,charged];u=V[charged,i].toarray().ravel();w=V[j,charged].toarray().ravel()
 x=inv*u
 for _ in range(6):x=inv*(Vc@x)
 measured=float(w@x);expected=coefficient(4,a)['float_coefficient']
 assert abs(measured-expected)<1e-12
 rows.append(dict(graph=label,dimension=len(states),a=str(a),measured=measured,expected=expected))
(P/'embedding_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
