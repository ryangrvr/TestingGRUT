import importlib.util,json
from pathlib import Path
import numpy as np
from scipy.sparse import kron,eye
from scipy.sparse.linalg import eigsh,expm_multiply
p=Path(__file__).resolve();root=p.parents[2]
src=(root/'constructive_physics/attempt09/check.py').read_text()
ns={'__file__':str(root/'constructive_physics/attempt09/check.py')};exec(compile(src.split("out={'cases':[]}")[0],ns['__file__'],'exec'),ns)
c=ns['c'];build=ns['build'];out=[]
for n in [2,3,4]:
 H,states,pairs,index,U,edges,degree=build(n);m=1<<n
 vals,vec=np.linalg.eigh(H.toarray()) if n==2 else eigsh(H,k=12,which='SA',tol=1e-10,v0=np.ones(H.shape[0]));order=np.argsort(vals);vals=vals[order];vec=vec[:,order];g=vec[:,0]
 A=kron(eye(len(states)),c.pauli(n,1,'Z'),format='csr');B=kron(eye(len(states)),c.pauli(n,0,'X'),format='csr')
 rows=[]
 for t in [.1,.4,1.]:
  w=expm_multiply(-1j*t*H,g);wb=expm_multiply(-1j*t*H,B@g)
  R=float((-2*np.vdot(w,A@wb).imag))
  eps=1e-5;vp=expm_multiply(-1j*t*H,np.cos(eps)*g+1j*np.sin(eps)*(B@g));vm=expm_multiply(-1j*t*H,np.cos(eps)*g-1j*np.sin(eps)*(B@g))
  derivative=float((np.vdot(vp,A@vp).real-np.vdot(vm,A@vm).real)/(2*eps));assert abs(derivative-R)<1e-8
  # Same-site noise from the full exponential, not a 12-mode truncation.
  wa=expm_multiply(-1j*t*H,A@g);corr=np.vdot(w,A@wa);mean=float(np.vdot(g,A@g).real)
  rows.append(dict(t=t,cross_response=R,kick_derivative=derivative,detector_connected_noise=float(corr.real-mean**2)))
 residual=max(np.linalg.norm(H@vec[:,i]-vals[i]*vec[:,i]) for i in range(12))
 assert residual<1e-7
 gp=c.probabilities(g,n)
 out.append(dict(N=n,dimension=H.shape[0],lowest_12_energies=vals.tolist(),gap=float(vals[1]-vals[0]),max_eigen_residual=float(residual),ground_geometry=c.geometry(gp,states,pairs,n),responses=rows))
 print('completed N',n,flush=True)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_LOW_ENERGY_RESPONSE_PASS');print(json.dumps(out,indent=2))
