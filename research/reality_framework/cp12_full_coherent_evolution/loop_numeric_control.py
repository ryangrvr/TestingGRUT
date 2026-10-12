"""Independent sparse-eigenvalue comparison against the exact loop series.
Small-hopping diagnostic, not the frozen-coupling full graph Hamiltonian.
"""
from pathlib import Path
import itertools,json,sys
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build
from fractions import Fraction
series=json.loads((P/'loop_series_results.json').read_text())
energies={};residuals={}
for n in [3,6]:
 edges=list(itertools.combinations(range(n),2))
 mask=sum(1<<edges.index(tuple(sorted((i,(i+1)%n)))) for i in range(n))
 H,*_=build(n,n,3,masks=[mask],rewire=False)
 D=diags(H.diagonal());V=H-D
 for t in [.12,.16,.20,.24]:
  A=D+t*V
  e,u=eigsh(A,k=1,which='SA',tol=1e-14,v0=np.ones(A.shape[0]))
  energies[n,t]=float(e[0]);residuals[n,t]=float(np.linalg.norm(A@u[:,0]-e[0]*u[:,0]))
rows=[]
for t in [.12,.16,.20,.24]:
 delta=energies[6,t]-2*energies[3,t]
 pred=sum(float(Fraction(r['ring6_minus_two_triangles']))*t**r['order'] for r in series)
 assert delta>0
 assert max(residuals[n,t] for n in [3,6])<1e-11
 rows.append(dict(hopping=t,energy_difference=delta,series_through_14=pred,residual=max(residuals[n,t] for n in [3,6])))
(P/'loop_numeric_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
