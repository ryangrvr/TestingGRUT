import json
from pathlib import Path
import numpy as np
from scipy.sparse import block_diag
from scipy.sparse.linalg import eigsh,expm_multiply
p=Path(__file__).resolve();f=p.parents[2]/'constructive_physics/attempt09/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split("out={'cases':[]}")[0],str(f),'exec'),ns);c=ns['c'];build=ns['build'];rows=[]
for n in [2,3,4]:
 H,states,pairs,index,U,edges,degree=build(n);m=1<<n
 source=[];sink=[];sizes=[];disjoint_source=[];disjoint_sink=[]
 for mask in states:
  C=np.zeros((n,n),bool)
  for k,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>k&1)
  lower=np.flatnonzero(~C.any(axis=0));upper=np.flatnonzero(~C.any(axis=1))
  # Minimal and maximal roles can overlap. Do not conceal that overlap.
  source.append(sum(c.pauli(n,int(i),'X') for i in lower)/len(lower))
  sink.append(sum(c.pauli(n,int(i),'Z') for i in upper)/len(upper))
  sizes.append((len(lower),len(upper),len(set(lower)&set(upper))))
  isolated=set(lower)&set(upper);lo=[int(i) for i in lower if i not in isolated];up=[int(i) for i in upper if i not in isolated]
  zero=c.pauli(n,0,'X')*0
  disjoint_source.append(sum((c.pauli(n,i,'X') for i in lo),zero)/max(1,len(lo)))
  disjoint_sink.append(sum((c.pauli(n,i,'Z') for i in up),zero)/max(1,len(up)))
 B=block_diag(source,format='csr');A=block_diag(sink,format='csr')
 vals,v=eigsh(H,k=1,which='SA',tol=1e-11,v0=np.ones(H.shape[0]));g=v[:,0];gp=c.probabilities(g,n)
 overlap=float(sum(w*o/(a*b) for w,(a,b,o) in zip(gp,sizes)))
 tests=[]
 for t in [0.,.1,.4,1.]:
  w=expm_multiply(-1j*t*H,g);wb=expm_multiply(-1j*t*H,B@g);R=float(-2*np.vdot(w,A@wb).imag)
  eps=1e-5;vp=expm_multiply(-1j*t*H,expm_multiply(1j*eps*B,g));vm=expm_multiply(-1j*t*H,expm_multiply(-1j*eps*B,g))
  derivative=float((np.vdot(vp,A@vp).real-np.vdot(vm,A@vm).real)/(2*eps));assert abs(R-derivative)<1e-8
  tests.append(dict(t=t,response=R,impulse_derivative=derivative))
 BD=block_diag(disjoint_source,format='csr');AD=block_diag(disjoint_sink,format='csr')
 assert np.linalg.norm((AD@BD-BD@AD).data)<1e-12
 clean=[]
 for t in [0.,.1,.4,1.]:
  w=expm_multiply(-1j*t*H,g);wb=expm_multiply(-1j*t*H,BD@g)
  R=float(-2*np.vdot(w,AD@wb).imag);eps=1e-5
  vp=expm_multiply(-1j*t*H,expm_multiply(1j*eps*BD,g));vm=expm_multiply(-1j*t*H,expm_multiply(-1j*eps*BD,g))
  derivative=float((np.vdot(vp,AD@vp).real-np.vdot(vm,AD@vm).real)/(2*eps));assert abs(R-derivative)<1e-8
  clean.append(dict(t=t,response=R,impulse_derivative=derivative))
 rows.append(dict(disjoint_role_response=clean,N=n,role_overlap_weight=overlap,operator_commutator_norm=float(np.linalg.norm((A@B-B@A).data)),response=tests))
 print(n,tests,flush=True)
p.with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n');print('CP09_RELATIONAL_ROLE_PROBE_PASS')
