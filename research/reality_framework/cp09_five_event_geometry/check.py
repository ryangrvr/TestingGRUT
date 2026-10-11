from pathlib import Path
import itertools,json,math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
p=Path(__file__).resolve()
def orders(n):
 pairs=list(itertools.permutations(range(n),2));pi={e:k for k,e in enumerate(pairs)};und=list(itertools.combinations(range(n),2));states=[]
 for digits in itertools.product(range(3),repeat=len(und)):
  C=np.zeros((n,n),bool);mask=0
  for z,(i,j) in zip(digits,und):
   if z:
    a,b=(i,j) if z==1 else (j,i);C[a,b]=1;mask|=1<<pi[a,b]
  if np.all((C.astype(int)@C.astype(int)>0)<=C):states.append(mask)
 return sorted(states),pairs
def build(n):
 states,pairs=orders(n);lookup={x:i for i,x in enumerate(states)};M=1<<n;D=len(states)*M
 rows=[];cols=[];data=[];degree=np.zeros(len(states));v=np.arange(M)
 for a,mask in enumerate(states):
  for k,(i,j) in enumerate(pairs):
   new=mask|(1<<k)
   if new==mask or new not in lookup:continue
   b=lookup[new];degree[a]+=1;degree[b]+=1
   # U=cos I-i sin Z_i X_j. Since i != j, Z_i commutes with X_j.
   flip=v^(1<<(n-1-j));z=1-2*((v>>(n-1-i))&1)
   for dest,amp in [(v,np.full(M,-math.cos(math.pi/4),complex)),(flip,1j*math.sin(math.pi/4)*z)]:
    rr=b*M+dest;cc=a*M+v
    rows.extend([rr,cc]);cols.extend([cc,rr]);data.extend([amp,amp.conj()])
 for a in range(len(states)):
  rr=a*M+v;rows.append(rr);cols.append(rr);data.append(np.full(M,degree[a]+.25*n,complex))
  for i in range(n):rows.append(a*M+(v^(1<<(n-1-i))));cols.append(rr);data.append(np.full(M,-.25,complex))
 H=coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(D,D)).tocsr()
 return H,states,pairs
out=[]
for n,expected in [(3,19),(4,219),(5,4231)]:
 H,states,pairs=build(n);assert len(states)==expected
 values,vectors=eigsh(H,k=5,which='SA',tol=1e-9,v0=np.ones(H.shape[0]));ix=np.argsort(values);values=values[ix];vectors=vectors[:,ix];g=vectors[:,0]
 residual=max(np.linalg.norm(H@vectors[:,k]-values[k]*vectors[:,k]) for k in range(5));assert residual<1e-6
 if n==3:assert abs(values[0]-.3109792266442)<1e-8
 if n==4:assert abs(values[0]-.6037576453941)<1e-8
 prob=np.sum(abs(g.reshape(-1,1<<n))**2,axis=1);hist={};counts={};relations=0
 for w,mask in zip(prob,states):
  C=np.zeros((n,n),bool)
  for k,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>k&1)
  h=np.ones(n,int)
  for _ in range(n):h=np.array([1+max(h[np.flatnonzero(C[:,j])],default=0) for j in range(n)])
  height=int(h.max());counts[height]=counts.get(height,0)+1;hist[height]=hist.get(height,0)+float(w);relations+=w*mask.bit_count()
 row=dict(N=n,posets=len(states),dimension=H.shape[0],nonzeros=H.nnz,lowest_energies=values.tolist(),gap=float(values[1]-values[0]),eigen_residual=float(residual),height_distribution=hist,uniform_height_distribution={k:v/len(states) for k,v in counts.items()},mean_height=sum(k*w for k,w in hist.items()),mean_relations=float(relations));out.append(row);print(json.dumps(row),flush=True)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_FIVE_EVENT_GEOMETRY_PASS')
