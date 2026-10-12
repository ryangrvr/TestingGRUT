"""Independent explicit orbit isometry checks against the full physical Hamiltonian."""
from pathlib import Path
import json,math
import numpy as np
from scipy.sparse import coo_matrix,eye
from model import build
from symmetric_sector import reduced_build
rows=[]
for n,q,L in [(3,2,1),(3,3,2),(4,4,1)]:
 H,states,lookup,degree=build(n,q,L)
 R,basis,sizes,meta,S,local=reduced_build(n,q,L)
 rr=[];cc=[];vv=[]
 for i,(mask,bits,E) in enumerate(states):
  rep,p=S.canonical[mask];b,f,sg=S.transform(p,bits,E);a,ss=local[rep][b,f]
  if a>=0:rr.append(i);cc.append(a);vv.append(sg*ss/math.sqrt(sizes[a]))
 W=coo_matrix((vv,(rr,cc)),shape=(len(states),len(basis))).tocsr()
 orth=np.linalg.norm((W.T@W-eye(len(basis))).data)
 matrix=np.linalg.norm((W.T@H@W-R).data)
 leakage=np.linalg.norm((H@W-W@R).data)
 assert max(orth,matrix,leakage)<1e-10
 rows.append(dict(N=n,Q=q,L=L,full=len(states),reduced=len(basis),orthogonality_error=orth,projected_matrix_error=matrix,sector_leakage=leakage))
Path(__file__).with_name('reduction_validation.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
