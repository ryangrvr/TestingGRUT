import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
import qball
import numpy as np
from scipy.integrate import solve_bvp
out=Path(__file__).parent
O=.8
trials=[(2.3,3.,1.),(2.3,4.,1.),(2.3,5.,1.),(2.,4.,1.5),(1.9,5.,1.),(2.3,6.,1.),(2.3,8.,1.)]
ans={'Omega_fixed':O,'a':1,'b':.1,'predeclared_numerical_guesses':trials,'attempts':[]}
(out/'results_preparation_search.json').write_text(json.dumps(ans,indent=2))
r=np.linspace(0,30,601)
def fun(r,y): return qball.fun(r,y,np.array([O]))
def bc(ya,yb): return np.array([ya[1],ya[3],yb[0],yb[2]-1])
for A,R,w in trials:
 t=np.tanh((r-R)/w);dt=(1-t*t)/w
 F=A*.5*(1-t); S=.5*(1+t)
 y=np.array([F,-A*.5*dt,S,.5*dt]);y[1,0]=y[3,0]=0
 start=time.perf_counter()
 sol=solve_bvp(fun,bc,r,y,S=np.diag([0.,-2.,0.,-2.]),tol=1e-5,max_nodes=30000)
 rr=np.linspace(0,30,10001); yy=sol.sol(rr)
 a={'A':A,'R':R,'w':w,'status':sol.status,'max_residual':np.max(sol.rms_residuals),'F0':yy[0,0],'S0':yy[2,0],'minF':yy[0].min(),'minS':yy[2].min(),'seconds':time.perf_counter()-start}
 ans['attempts'].append(a)
 print(a,flush=True)
 (out/'results_preparation_search.json').write_text(json.dumps(ans,indent=2,default=qball.serial))
 if sol.success and yy[0,0]>.1 and yy[0].min()>-1e-6 and yy[2].min()>0:
  np.savez(out/'results_preparation_seed.npz',r=rr,y=yy,Omega=O)
  break
