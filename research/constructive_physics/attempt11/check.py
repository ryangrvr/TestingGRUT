"""CP11: exact conditional enumeration and joint canonical-fermion diagonalization."""
from pathlib import Path
import itertools,json,math,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh

MASS=HOP=1.; COST=.5; GAMMA=.1

def adjacency(n,mask):
 A=np.zeros((n,n))
 for e,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>e&1:A[i,j]=A[j,i]=1
 return A

def stats(A):
 n=len(A);seen=set();sizes=[]
 for i in range(n):
  if i in seen:continue
  stack=[i];seen.add(i);size=0
  while stack:
   j=stack.pop();size+=1
   for v in np.flatnonzero(A[j]):
    v=int(v)
    if v not in seen:seen.add(v);stack.append(v)
  sizes.append(size)
 return sorted(sizes,reverse=True)

def energy(A):
 L=np.diag(A.sum(axis=1))-A
 lam=eigvalsh(L);assert lam.min()>-1e-9
 return float(COST*A.sum()/2-np.sum(np.sqrt(MASS*MASS+HOP*HOP*lam)-MASS))

def conditional(n):
 M=n*(n-1)//2;best=1e9;minima=[];bestconnected=1e9
 for mask in range(1<<M):
  A=adjacency(n,mask);E=energy(A)
  assert E>=-n/4-1e-10
  if len(stats(A))==1:bestconnected=min(bestconnected,E)
  if E<best-1e-10:best=E;minima=[mask]
  elif abs(E-best)<1e-10:minima.append(mask)
 groups={}
 for mask in minima:
  A=adjacency(n,mask);key=str((stats(A),sorted(A.sum(axis=1).astype(int).tolist())))
  groups[key]=groups.get(key,0)+1
 return dict(N=n,graphs=1<<M,minimum=best,minimizer_count=len(minima),minimizer_components_and_degrees=groups,best_connected=bestconnected)

def hop(bits,dest,src):
 if not(bits>>src&1) or bits>>dest&1:return None
 sign=(-1)**((bits&((1<<src)-1)).bit_count());b=bits^(1<<src)
 sign*=(-1)**((b&((1<<dest)-1)).bit_count())
 return b|(1<<dest),sign

def build_joint(n):
 edges=list(itertools.combinations(range(n),2));M=len(edges);F=n+M
 occ=[sum(1<<i for i in ss) for ss in itertools.combinations(range(F),M)]
 lookup={b:i for i,b in enumerate(occ)};d=len(occ);D=(1<<M)*d
 rr=[];cc=[];vv=[]
 for graph in range(1<<M):
  for a,bits in enumerate(occ):
   index=graph*d+a;nv=(bits&((1<<n)-1)).bit_count()
   rr.append(index);cc.append(index);vv.append(GAMMA*M/n+COST*graph.bit_count()+2*MASS*nv)
   for e,(i,j) in enumerate(edges):
    rr.append((graph^(1<<e))*d+a);cc.append(index);vv.append(-GAMMA/n)
    if graph>>e&1:
     for vertex,orient in [(i,1),(j,-1)]:
      for dest,src in [(vertex,n+e),(n+e,vertex)]:
       v=hop(bits,dest,src)
       if v is not None:
        new,sign=v;rr.append(graph*d+lookup[new]);cc.append(index);vv.append(HOP*orient*sign)
 H=coo_matrix((vv,(rr,cc)),shape=(D,D)).tocsr()
 assert np.linalg.norm((H-H.T).data)<1e-12
 return H,occ,edges

def joint(n):
 H,occ,edges=build_joint(n);M=len(edges);d=len(occ);D=H.shape[0]
 values,vec=eigsh(H,k=4,which='SA',tol=1e-10,v0=np.random.default_rng(1100+n).normal(size=D))
 order=np.argsort(values);values=values[order];vec=vec[:,order]
 residual=max(np.linalg.norm(H@vec[:,i]-values[i]*vec[:,i]) for i in range(4))
 assert residual<1e-7
 assert values[0]>=-n/4-1e-10
 assert values[0]<=GAMMA*(n-1)/2+1e-10
 p=np.sum(abs(vec[:,0].reshape(-1,d))**2,axis=1)
 comps={};edgesdist={};meanE=0
 for mask,w in enumerate(p):
  comp=stats(adjacency(n,mask));key=str(comp);comps[key]=comps.get(key,0)+float(w)
  e=mask.bit_count();edgesdist[e]=edgesdist.get(e,0)+float(w);meanE+=e*w
 # one-particle spectrum independently validates filled-sector formula for every graph
 for mask in range(1<<M):
  B=np.zeros((n,M))
  for e,(i,j) in enumerate(edges):
   if mask>>e&1:B[i,e]=1;B[j,e]=-1
  h=np.block([[MASS*np.eye(n),HOP*B],[HOP*B.T,-MASS*np.eye(M)]])
  E=np.sum(eigvalsh(h)[:M])+MASS*M+COST*mask.bit_count()
  assert abs(E-energy(adjacency(n,mask)))<1e-10
 return dict(N=n,dimension=D,nonzeros=H.nnz,energies=values.tolist(),gap=float(values[1]-values[0]),residual=float(residual),component_probability=comps,edge_probability=edgesdist,mean_edges=float(meanE))

if __name__=='__main__':
 results={'conditional':[],'joint':[]}
 for n in range(2,7):
  row=conditional(n);results['conditional'].append(row);print('conditional',row,flush=True)
 for n in range(2,5):
  row=joint(n);results['joint'].append(row);print('joint',row,flush=True)
 Path(__file__).with_name('results.json').write_text(json.dumps(results,indent=2)+'\n')
 print('CP11_CONDITIONAL_AND_JOINT_CHECKS_PASS')
