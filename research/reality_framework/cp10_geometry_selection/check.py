from pathlib import Path
import numpy as np,json
from scipy.sparse import kron,eye,diags,block_diag
from scipy.sparse.linalg import eigsh
p=Path(__file__).resolve();f=p.parents[1]/'cp09_five_event_geometry/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split('out=[]')[0],str(f),'exec'),ns);out=[]
for n in [2,3,4]:
 H,states,pairs=ns['build'](n);m=1<<n;laps=[];heights=[]
 for mask in states:
  C=np.zeros((n,n),int)
  for k,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>k&1)
  Adj=C+C.T;laps.append(kron(eye(m),np.diag(Adj.sum(axis=1))-Adj))
  h=np.ones(n,int)
  for _ in range(n):h=np.array([1+max(h[np.flatnonzero(C[:,j])],default=0) for j in range(n)])
  heights.append(int(h.max()))
 diag=np.tile(np.array([.5*((bits>>(n-1-i))&1) for bits in range(m) for i in range(n)]),len(states))
 H10=kron(H,eye(n),format='csr')+.25*block_diag(laps,format='csr')+diags(diag)
 E09,V09=eigsh(H,k=1,which='SA',tol=1e-10,v0=np.ones(H.shape[0]));E,V=eigsh(H10,k=5,which='SA',tol=1e-10,v0=np.ones(H10.shape[0]));ix=np.argsort(E);E=E[ix];V=V[:,ix]
 p09=np.sum(abs(V09[:,0].reshape(-1,m))**2,axis=1);p10=np.sum(abs(V[:,0].reshape(-1,m,n))**2,axis=(1,2))
 dist=lambda prob:{str(k):float(sum(w for w,h in zip(prob,heights) if h==k)) for k in range(1,n+1)}
 residual=float(np.linalg.norm(H10@V[:,0]-E[0]*V[:,0]));assert residual<1e-7
 row=dict(N=n,dimension=H10.shape[0],ground_energy=float(E[0]),gap=float(E[1]-E[0]),ground_residual=residual,cp09_height=dist(p09),cp10_height=dist(p10),geometry_total_variation=float(sum(abs(p10-p09))/2));out.append(row);print(row,flush=True)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP10_GROUND_GEOMETRY_CONTROL_PASS')
