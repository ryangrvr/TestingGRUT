"""Fourth-order degree coefficient with the actual graph-compression diagonal."""
from pathlib import Path
import itertools,sys,json
import numpy as np
from scipy.sparse import diags
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,GAMMA,MASS,ELECTRIC
rows=[]
for label,n,pairs in [('path4',4,[(0,1),(1,2),(2,3)]),('star4',4,[(0,1),(0,2),(0,3)]),('square',4,[(0,1),(1,2),(2,3),(0,3)]),('triangle_isolate',4,[(0,1),(1,2),(0,2)])]:
 edges=list(itertools.combinations(range(n),2));M=len(edges);q=len(pairs);active=[edges.index(tuple(sorted(e))) for e in pairs];mask=sum(1<<e for e in active)
 H,states,lookup,_=build(n,q,2,masks=[mask],rewire=False);V=H-diags(H.diagonal());a=GAMMA*(M-q)/M
 rem=np.array([sum((bits>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for g,bits,E in states])
 sea=sum(1<<(n+e) for e in range(M));i=lookup[mask,sea,(0,)*(2*M)];v=np.zeros(len(states));v[i]=1
 D=H.diagonal()+a*rem;D-=D[i];R=np.zeros_like(D);R[D>1e-12]=1/D[D>1e-12]
 w=V@v;E2=-float(w@(R*w));E4=-float(w@(R*(V@(R*(V@(R*w))))))+float(w@(R*R*w))*float(w@(R*w))
 deg=[sum(x in e for e in pairs) for x in range(n)];W=sum(d*(d-1)//2 for d in deg);delta=2*MASS+ELECTRIC/2-a
 assert abs(E2+2*q/delta)<1e-11 and abs(E4-(4*q+2*W)/delta**3)<1e-11
 rows.append(dict(graph=label,a=a,delta=delta,E2=E2,E4=E4,degree_coefficient=2/delta**3))
(P/'compressed_degree_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
