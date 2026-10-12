"""Distance-resolved response of the exact compressed prepared-ring Hamiltonian."""
from pathlib import Path
import sys,itertools,json,time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,hop,GAMMA
rows=[]
for n,L in [(6,2),(6,3),(6,4),(8,2),(8,3),(8,4)]:
 start=time.time();edges=list(itertools.combinations(range(n),2));M=len(edges);eid={e:i for i,e in enumerate(edges)};mask=sum(1<<eid[tuple(sorted((i,(i+1)%n)))] for i in range(n));q=n
 H,states,lookup,_=build(n,q,L,masks=[mask],rewire=False)
 removable=np.array([sum((s[0]>>e&1) and (s[1]>>(n+e)&1) and not(s[2][2*e] or s[2][2*e+1]) for e in range(M)) for s in states])
 H=H+diags(GAMMA*(M-q)/M*removable)
 sea=sum(1<<(n+e) for e in range(M));F=[0]*(2*M);v=lookup[mask,sea,tuple(F)];bits,sg=hop(sea,0,n+eid[0,1]);F[2*eid[0,1]]=1;r=lookup[mask,bits,tuple(F)]
 psi=np.zeros(len(states),complex);psi[v]=1/np.sqrt(2);psi[r]=1j*sg/np.sqrt(2);B=np.array([(s[1]>>0)&1 for s in states]);A=np.array([[(s[1]>>j)&1 for s in states] for j in range(1,n//2+1)]);boundary=np.array([max(abs(e) for e in s[2])==L for s in states]);w=np.column_stack([psi,B*psi]);last=0;values=[]
 for T in [.4,1.,2.,4.]:
  w=expm_multiply(-1j*(T-last)*H,w);last=T;response=(-2*np.imag((w[:,0].conj()*w[:,1])@A.T)).tolist();values.append(dict(T=T,response_by_vertex_distance=response,cutoff_boundary=float(abs(w[:,0])**2@boundary)))
 row=dict(N=n,L=L,dimension=len(states),seconds=time.time()-start,rows=values);rows.append(row);print(n,L,len(states),values,flush=True)
Path(__file__).with_name('fixed_propagation_results.json').write_text(json.dumps(rows,indent=2)+'\n')
