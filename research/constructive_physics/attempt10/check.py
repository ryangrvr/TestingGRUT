from pathlib import Path
import numpy as np,json
from scipy.sparse import kron,eye,diags,block_diag
from scipy.sparse.linalg import eigsh,expm_multiply
p=Path(__file__).resolve();f=p.parents[2]/'reality_framework/cp09_five_event_geometry/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split('out=[]')[0],str(f),'exec'),ns)
rows=[]
for n in [2,3]:
 H,states,pairs=ns['build'](n);m=1<<n;energies,v=eigsh(H,k=1,which='SA',v0=np.ones(H.shape[0]),tol=1e-11);ground=v[:,0]
 laps=[];roles=[]
 for mask in states:
  C=np.zeros((n,n),int)
  for k,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>k&1)
  Adj=C+C.T;L=np.diag(Adj.sum(axis=1))-Adj;laps.append(kron(eye(m),L));roles.extend(np.tile((C.sum(axis=1)-C.sum(axis=0))/(n-1),m))
 HR=.25*block_diag(laps,format='csr');Htotal=kron(H,eye(n),format='csr')+HR
 diag=[]
 for mask in states:
  for bits in range(m):
   diag.extend([.25*(1-(1-2*((bits>>(n-1-i))&1))) for i in range(n)])
 interaction=diags(diag);Htotal+=interaction
 vals,_=eigsh(Htotal,k=1,which='SA',v0=np.ones(Htotal.shape[0]),tol=1e-10);assert vals[0]>-1e-9
 # Permutation swaps event 0/1 on ALL subsystems.
 perm=list(range(n));perm[0],perm[1]=perm[1],perm[0];lookup={mask:a for a,mask in enumerate(states)};pi={e:k for k,e in enumerate(pairs)};mapping=[]
 for mask in states:
  new=sum(1<<pi[(perm[i],perm[j])] for k,(i,j) in enumerate(pairs) if mask>>k&1)
  for bits in range(m):
   newbits=sum(((bits>>(n-1-i))&1)<<(n-1-perm[i]) for i in range(n))
   mapping.extend([((lookup[new]*m+newbits)*n+perm[i]) for i in range(n)])
 err=Htotal[mapping,:][:,mapping]-Htotal;assert np.linalg.norm(err.data)<1e-12
 initial=np.kron(ground,np.ones(n)/np.sqrt(n));E=float(np.vdot(initial,Htotal@initial).real);output=[]
 role=np.asarray(roles)
 for t in [0.,.4,1.]:
  w=expm_multiply(-1j*t*Htotal,initial);assert abs(np.linalg.norm(w)-1)<1e-10;assert abs(np.vdot(w,Htotal@w).real-E)<1e-10
  matrix=w.reshape(len(states)*m,n);rhoR=matrix.conj().T@matrix;purity=float(np.trace(rhoR@rhoR).real)
  joint=np.sum(abs(w.reshape(len(states),m,n))**2,axis=1);pc=joint.sum(axis=1);pr=joint.sum(axis=0)
  MI=float(sum(joint[a,i]*np.log(joint[a,i]/(pc[a]*pr[i])) for a in range(len(states)) for i in range(n) if joint[a,i]>1e-15))
  output.append(dict(t=t,reference_purity=purity,reference_geometry_mutual_information=MI,relational_orientation=float(np.sum(abs(w)**2*role)),energy=float(np.vdot(w,Htotal@w).real)))
 back=expm_multiply(1j*Htotal,expm_multiply(-1j*Htotal,initial));assert np.linalg.norm(back-initial)<1e-10
 rows.append(dict(N=n,dimension=Htotal.shape[0],ground_energy=float(vals[0]),permutation_error=float(np.linalg.norm(err.data)),initial_energy=E,controls=output));print(rows[-1],flush=True)
p.with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n');print('CP10_JOINT_DYNAMICAL_REFERENCE_PASS')
