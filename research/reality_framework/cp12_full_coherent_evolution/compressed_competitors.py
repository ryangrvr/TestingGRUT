"""Compare two N=6,Q=6 graph blocks, retaining the true rewiring diagonal."""
from pathlib import Path
import sys,itertools,json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,GAMMA
rows=[];n=6;q=6;edges=list(itertools.combinations(range(n),2));M=len(edges);a=GAMMA*(M-q)/M
for L in [3,4]:
 for label,pairs in [('ring6',[(i,(i+1)%6) for i in range(6)]),('two_triangles',[(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)])]:
  mask=sum(1<<edges.index(tuple(sorted(e))) for e in pairs)
  H,states,_,_=build(n,q,L,masks=[mask],rewire=False)
  active=[e for e in range(M) if mask>>e&1]
  neutral=np.array([sum((bits>>(n+e)&1) and E[2*e]==0 and E[2*e+1]==0 for e in active) for g,bits,E in states])
  H=H+diags(a*neutral);E,U=eigsh(H,k=1,which='SA',tol=1e-11);v=U[:,0]
  row=dict(label=label,N=n,Q=q,L=L,a=a,dimension=len(states),compressed_energy=float(E[0]),
           per_vertex=float(E[0]/n),neutral_edges=float(abs(v)**2@neutral),
           residual=float(np.linalg.norm(H@v-E[0]*v)))
  rows.append(row);print(json.dumps(row),flush=True)
Path(__file__).with_name('compressed_competitors_results.json').write_text(json.dumps(rows,indent=2)+'\n')
