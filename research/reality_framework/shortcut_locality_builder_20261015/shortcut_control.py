"""Independent one-particle transport control. NOT a CP-12 simulation."""
import json
import numpy as np
from scipy.linalg import expm

def transport(n, chords, t, source=0):
    core=np.zeros((n,n))
    for i in range(n):
        j=(i+1)%n
        core[i,j]=core[j,i]=1.
    full=core.copy()
    for i,j in chords:
        full[i,j]=full[j,i]=1.
    initial=np.eye(n,dtype=complex)[:,source]
    return core,full,expm(-1j*t*core)@initial,expm(-1j*t*full)@initial

def run():
    n=10; chords=[(0,5),(2,7)]
    rows=[]
    for t in [0.1,0.25,0.5,1.0]:
        core,full,a,b=transport(n,chords,t)
        v=full-core
        grid=np.linspace(0,t,1001)
        integrand=np.array([np.linalg.norm(v@expm(-1j*s*core)[:,0]) for s in grid])
        integral=np.trapezoid(integrand,grid)
        rows.append(dict(time=t,core_opposite_prob=float(abs(a[5])**2),
          full_opposite_prob=float(abs(b[5])**2),
          state_distance=float(np.linalg.norm(a-b)),
          duhamel_integral=float(integral),
          inequality_holds=bool(np.linalg.norm(a-b)<=integral+1e-8)))
    result=dict(status="PASS" if all(x["inequality_holds"] for x in rows) else "FAIL",
      scope="10-site one-particle cycle plus two chords; illustrative control only; not CP-12",
      n=n,chords=chords,rows=rows)
    print(json.dumps(result,indent=2))
if __name__=="__main__":
    run()
