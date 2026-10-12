"""Independent full-basis charged-sector elimination, not a path-subset model."""
from pathlib import Path
from fractions import Fraction as F
import itertools,sys,json
import numpy as np
from scipy.sparse import diags,eye
from scipy.sparse.linalg import spsolve
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,wilson_vector
from loop_paths import coefficient
rows=[]
for n in [3,4]:
 edges=list(itertools.combinations(range(n),2));M=len(edges)
 mask=sum(1<<edges.index(tuple(sorted((i,(i+1)%n)))) for i in range(n))
 for L in [1,2,3]:
  H,states,lookup,_=build(n,n,L,masks=[mask],rewire=False)
  D=H.diagonal();V=H-diags(D)
  sea=sum(1<<(n+e) for e in range(M));vac=(mask,sea,(0,)*(2*M));flux=wilson_vector(n,list(range(n)),vac)
  i,j=lookup[vac],lookup[flux]
  charged=np.array([k for k,(g,b,E) in enumerate(states) if b&((1<<n)-1)])
  active=[e for e in range(M) if mask>>e&1]
  removable=np.array([sum((b>>(n+e)&1) and E[2*e]==0 and E[2*e+1]==0 for e in active) for g,b,E in states])
  Vc=V[charged][:,charged];u=V[charged,i].toarray().ravel();w=V[j,charged].toarray().ravel()
  for a in [F(0),F(3,50)]:
   gap=(D+float(a)*removable)-(D[i]+float(a)*n)
   inv=-1/gap[charged]
   vec=inv*u
   for _ in range(2*n-2):vec=inv*(Vc@vec)
   exact_matrix=float(w@vec)
   expected=coefficient(n,a)['float_coefficient']
   assert abs(exact_matrix-expected)<1e-12
   for t in [.02,.04,.08,.16]:
    A=diags(-gap[charged])-t*Vc
    x=spsolve(A,t*u);amp=float(t*w@x)
    residual=float(np.linalg.norm(A@x-t*u))
    row=dict(length=n,L=L,a=str(a),hopping=t,amplitude=amp,scaled_amplitude=amp/t**(2*n),leading_coefficient=expected,relative_correction=amp/t**(2*n)/expected-1,residual=residual,matrix_series_coefficient=exact_matrix)
    assert residual<1e-10;rows.append(row)
# Full charged-sector series and path enumeration agree; convergence at small t.
assert max(abs(r['relative_correction']) for r in rows if r['hopping']==.02)<.01
out=dict(status='INDUCED_WILSON_LOOP_CONTROLS_PASS',rows=rows)
(P/'feshbach_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'rows':len(rows),'max_small_hopping_relative_correction':max(abs(r['relative_correction']) for r in rows if r['hopping']==.02)},indent=2))
