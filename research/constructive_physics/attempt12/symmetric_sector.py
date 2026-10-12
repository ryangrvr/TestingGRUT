"""Vacuum-relative permutation sector; validate against full CP12 before extrapolating.
This is a symmetry-sector calculation, not proof of the full-system ground for untested N.
"""
from pathlib import Path
import itertools,math,json,time,argparse
from collections import defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from model import graph_basis,hop,MASS,HOP,COST,ELECTRIC,GAMMA

class Symmetry:
 def __init__(self,n,q):
  self.n=n;self.edges=list(itertools.combinations(range(n),2));self.M=len(self.edges);ei={e:i for i,e in enumerate(self.edges)}
  self.perms=list(itertools.permutations(range(n)))
  self.modes=[];self.orients=[];self.legs=[];self.vacuum=[];self.graphmap=[]
  self.masks=np.array([sum(1<<e for e in es) for es in itertools.combinations(range(self.M),q)],dtype=np.int64)
  self.maskindex={int(m):a for a,m in enumerate(self.masks)}
  maps=np.zeros((len(self.perms),len(self.masks)),dtype=np.int64)
  for a,p in enumerate(self.perms):
   modes=list(p);orient=[1]*n;legs=[]
   for e,(i,j) in enumerate(self.edges):
    u,v=p[i],p[j];f=ei[tuple(sorted((u,v)))];modes.append(n+f);orient.append(1 if u<v else -1)
    legs.extend([2*f+(0 if u<v else 1),2*f+(1 if u<v else 0)])
    maps[a]|=((self.masks>>e)&1)<<f
   parity=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
   self.modes.append(modes);self.orients.append(orient);self.legs.append(legs);self.vacuum.append((-1)**(parity*(n-1)))
  smallest=maps.min(axis=0);which=maps.argmin(axis=0)
  self.canonical={int(m):(int(smallest[i]),int(which[i])) for i,m in enumerate(self.masks)}
  self.reps=sorted(set(smallest.tolist()))
  self.aut={m:np.flatnonzero(maps[:,self.maskindex[m]]==m).tolist() for m in self.reps}
  self.bitscache={}
 def transform(self,p,bits,E):
  key=(p,bits)
  if key in self.bitscache:newbits,sign=self.bitscache[key]
  else:
   mapped=[self.modes[p][i] for i in range(self.n+self.M) if bits>>i&1]
   parity=0;seen=0
   for v in mapped:parity+=(seen>>(v+1)).bit_count();seen|=1<<v
   sign=self.vacuum[p]*(-1)**parity*math.prod(self.orients[p][i] for i in range(self.n+self.M) if bits>>i&1)
   newbits=seen;self.bitscache[key]=(newbits,sign)
  newE=[0]*(2*self.M)
  for i,j in enumerate(self.legs[p]):newE[j]=E[i]
  return newbits,tuple(newE),sign

def reduced_build(n,q,L):
 start=time.time();S=Symmetry(n,q);M=S.M;local={};basis=[];sizes=[];counts=[];annihilated=0
 for mask in S.reps:
  states=[(bits,E) for _,bits,E in graph_basis(n,mask,L)];allset=set(states);table={};aut=S.aut[mask];multiplicity=math.factorial(n)//len(aut)
  for bits,E in states:
   if (bits,E) in table:continue
   orbit={};zero=False
   for p in aut:
    b,f,sg=S.transform(p,bits,E);assert (b,f) in allset
    if (b,f) in orbit and orbit[b,f]!=sg:zero=True
    orbit[b,f]=sg
   if zero:
    for state in orbit:table[state]=(-1,0)
    annihilated+=len(orbit)*multiplicity
   else:
    a=len(basis);basis.append((mask,bits,E));sizes.append(len(orbit)*multiplicity)
    for state,sg in orbit.items():table[state]=(a,sg)
  local[mask]=table;counts.append(dict(mask=mask,physical_states=len(states),automorphisms=len(aut),labelled_graphs=multiplicity))
 print('basis',dict(N=n,Q=q,L=L,orbits=len(basis),unlabelled_graphs=len(S.reps),seconds=time.time()-start),flush=True)
 rows=[];cols=[];data=[];sizes=np.array(sizes,float)
 def add(a,newmask,bits,E,coefficient):
  rep,p=S.canonical[newmask];b,f,sign=S.transform(p,bits,E);dest,sg=local[rep][b,f]
  if dest<0:return
  rows.append(dest);cols.append(a);data.append(coefficient*sign*sg*np.sqrt(sizes[a]/sizes[dest]))
 for a,(mask,bits,E) in enumerate(basis):
  diagonal=COST*q+2*MASS*(bits&((1<<n)-1)).bit_count()+ELECTRIC/2*sum(v*v for v in E)
  present=[e for e in range(M) if mask>>e&1];deg=0
  for e in present:
   i,j=S.edges[e]
   for leg,vertex,orient in [(0,i,1),(1,j,-1)]:
    for dest,src,shift in [(vertex,n+e,1),(n+e,vertex,-1)]:
     h=hop(bits,dest,src)
     if h is None:continue
     newbits,sign=h;F=list(E);F[2*e+leg]+=shift
     if abs(F[2*e+leg])>L:continue
     add(a,mask,newbits,tuple(F),HOP*orient*sign)
   if E[2*e] or E[2*e+1] or not(bits>>(n+e)&1):continue
   for f in range(M):
    if mask>>f&1:continue
    deg+=1;add(a,mask^(1<<e)^(1<<f),bits,E,-GAMMA/M)
  diagonal+=GAMMA/M*deg;rows.append(a);cols.append(a);data.append(diagonal)
 H=coo_matrix((data,(rows,cols)),shape=(len(basis),len(basis))).tocsr();err=float(np.linalg.norm((H-H.T).data));assert err<1e-10,err
 return H,basis,sizes,dict(graph_counts=counts,projected_out_states=annihilated,full_basis_count=int(sum(sizes))+annihilated,hermiticity_error=err,seconds=time.time()-start),S,local

def run(n,q,L):
 start=time.time();H,basis,sizes,meta,S,local=reduced_build(n,q,L)
 E,V=eigsh(H,k=3,which='SA',tol=1e-10,v0=np.random.default_rng(12000+n+q).normal(size=H.shape[0]));ix=np.argsort(E);E=E[ix];V=V[:,ix]
 residual=max(np.linalg.norm(H@V[:,i]-E[i]*V[:,i]) for i in range(3));assert residual<1e-7
 p=abs(V[:,0])**2;graphs=defaultdict(float);boundary=0
 for w,(mask,bits,fields) in zip(p,basis):graphs[mask]+=float(w);boundary+=w*(max(abs(v) for v in fields)==L)
 result=dict(N=n,Q=q,cutoff=L,sector='vacuum-relative simultaneous permutation invariant',dimension=H.shape[0],energies=E.tolist(),residual=residual,cutoff_boundary_probability=float(boundary),unlabelled_graph_probability=dict(graphs),metadata=meta,seconds=time.time()-start)
 f=Path(__file__).with_name(f'symmetric_N{n}_Q{q}_L{L}.json');f.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['metadata','unlabelled_graph_probability']}),flush=True)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--q',type=int,default=4);p.add_argument('--L',type=int,default=2);a=p.parse_args();run(a.n,a.q,a.L)
