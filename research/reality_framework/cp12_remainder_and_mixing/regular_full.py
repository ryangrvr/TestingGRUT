"""Full CP12 versus exact-regular and nonregular sectors.
Graph-diagonal gauge-invariant density response from a coherent regular
preparation. No physical distance, clock, or large-N ground phase is inferred.
"""
from pathlib import Path
import sys,json,itertools,argparse,math,gc
import numpy as np
from collections import defaultdict
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import LinearOperator,eigsh,expm_multiply
from regular_graph_bound import regular_masks
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,GAMMA

def lowest(H):
 if H.shape[0]<3:
  E,V=np.linalg.eigh(H.toarray());return float(E[0]),V[:,0]
 E,V=eigsh(H,k=1,which='SA',tol=1e-12,v0=np.random.default_rng(120127).normal(size=H.shape[0]));return float(E[0]),V[:,0]

def run(n,d,L,response=True):
 q=n*d//2;M=n*(n-1)//2;H,states,lookup,_=build(n,q,L);del lookup
 masks=np.array([s[0] for s in states]);rset=set(regular_masks(n,d));regular=np.array([i for i,m in enumerate(masks) if m in rset]);other=np.array([i for i,m in enumerate(masks) if m not in rset]);A=H[regular][:,regular];C=H[other][:,other];B=H[regular][:,other]
 coo=A.tocoo();assert not np.any(masks[regular][coo.row]!=masks[regular][coo.col])
 b=GAMMA/M*math.sqrt((d+1)*q*(M-q))
 op=LinearOperator((len(regular),len(regular)),matvec=lambda x:B@(B.T@x),dtype=float)
 actual=math.sqrt(float(eigsh(op,k=1,which='LA',tol=1e-12,v0=np.ones(len(regular)),return_eigenvectors=False)[0]));assert actual<=b+1e-10
 E0,psi=lowest(H);Er,phir=lowest(A);Ec,phic=lowest(C);lower=min(Er,Ec)-actual
 assert lower<=E0+1e-9 and E0<=min(Er,Ec)+1e-9
 p_regular=float(np.linalg.norm(psi[regular])**2)
 fiber_amplitudes=defaultdict(lambda:[0.,0])
 for amplitude,(_,bits,E) in zip(psi,states):
  item=fiber_amplitudes[bits,E];item[0]+=amplitude;item[1]+=1
 kernel_weight=float(sum(abs(x)**2/count for x,count in fiber_amplitudes.values()))
 out=dict(N=n,d=d,Q=q,L=L,dimension=len(states),regular_dimension=len(regular),nonregular_dimension=len(other),full_energy=E0,regular_sector_energy=Er,nonregular_sector_energy=Ec,regular_ground_weight=p_regular,rewiring_kernel_ground_weight=kernel_weight,coupling_norm=actual,universal_coupling_bound=b,two_sector_energy_lower=lower,two_sector_energy_upper=min(Er,Ec),residual=float(np.linalg.norm(H@psi-E0*psi)),cutoff_boundary=float(sum(abs(psi[i])**2 for i,s in enumerate(states) if max(map(abs,s[2]))==L)),scope='Full regulated finite-host Hamiltonian. No thermodynamic phase or locality selection.')
 if n<=4:
  out['best_nonregular_fixed_block_energy']=min(lowest(H[ids][:,ids])[0] for mask in set(masks.tolist())-rset for ids in [np.flatnonzero(masks==mask)])
 if response:
  # Equal-amplitude coherent sum of physical matter/gauge ground states in
  # every regular graph, not a dephased mixture and not the full ground state.
  prep=np.zeros(len(states))
  for mask in sorted(rset):
   ids=np.flatnonzero(masks==mask);_,v=lowest(H[ids][:,ids]);prep[ids]=v/math.sqrt(len(rset))
  Ad=np.array([(s[1]>>1)&1 for s in states],float);Bd=np.array([s[1]&1 for s in states],float)
  coo=H.tocoo();cross=(np.isin(coo.row,regular)!=np.isin(coo.col,regular))
  H0=H-coo_matrix((coo.data[cross],(coo.row[cross],coo.col[cross])),shape=H.shape).tocsr()
  initial=np.column_stack([prep,Bd*prep]);rows=[]
  for t in [.1,.4,1.,2.]:
   U=expm_multiply(-1j*t*H,initial);V=expm_multiply(-1j*t*H0,initial)
   R=-2*float(np.imag(np.vdot(U[:,0],Ad*U[:,1])));R0=-2*float(np.imag(np.vdot(V[:,0],Ad*V[:,1])));error=abs(R-R0);bound=4*b*b*t*t
   departure=float(np.linalg.norm(U[other,0])**2)
   assert error<=bound+1e-10 and departure<=b*b*t*t+1e-10
   rows.append(dict(T=t,full_response=R,regular_block_response=R0,response_error=error,proven_error_bound=bound,departure=departure,departure_bound=b*b*t*t))
  out['coherent_density_response']=rows
 (P/f'regular_full_N{n}_d{d}_L{L}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--d',type=int,default=2);p.add_argument('--L',type=int,default=2);p.add_argument('--no-response',action='store_true');a=p.parse_args();run(a.n,a.d,a.L,not a.no_response)
