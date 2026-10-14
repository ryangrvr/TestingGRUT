"""CP12 exact physical U(1) basis. Electric cutoff is a regulator, not the law."""
import itertools
from array import array
import numpy as np
from scipy.sparse import coo_matrix

MASS=HOP=1.;COST=.5;ELECTRIC=.25;GAMMA=.1

def setup(n,mask):
 all_edges=list(itertools.combinations(range(n),2));present=[e for e in range(len(all_edges)) if mask>>e&1]
 B=np.zeros((n,len(present)),dtype=int)
 parent=list(range(n));tree=[];chords=[]
 def root(i):
  while parent[i]!=i:i=parent[i]
  return i
 for local,e in enumerate(present):
  i,j=all_edges[e];B[i,local]=1;B[j,local]=-1
  a,b=root(i),root(j)
  if a!=b:parent[a]=b;tree.append(local)
  else:chords.append(local)
 comps={}
 for i in range(n):comps.setdefault(root(i),[]).append(i)
 # For any divergence rhs, solve tree flows by peeling leaves; root residual checks neutrality.
 adjacency=[[] for _ in range(n)]
 for local in tree:
  i,j=all_edges[present[local]];adjacency[i].append((j,local));adjacency[j].append((i,local))
 elimination=[]
 for comp in comps.values():
  origin=comp[0];visited={origin};stack=[origin];order=[];parents={}
  while stack:
   i=stack.pop();order.append(i)
   for j,e in adjacency[i]:
    if j not in visited:visited.add(j);stack.append(j);parents[j]=(i,e)
  for leaf in reversed(order[1:]):
   neighbor,e=parents[leaf];elimination.append((leaf,neighbor,e,int(B[leaf,e])))
 return all_edges,present,B,tree,chords,list(comps.values()),elimination

def solve_tree(rhs,q,elimination):
 r=rhs.copy();f=np.zeros(q,dtype=int)
 for leaf,neighbor,e,sign in elimination:
  f[e]=r[leaf]*sign;r[neighbor]+=r[leaf];r[leaf]=0
 return f,r

def graph_basis(n,mask,L):
 all_edges,present,B,tree,chords,comps,elim=setup(n,mask);M=len(all_edges);q=len(present)
 inactive=sum(1<<(n+e) for e in range(M) if not(mask>>e&1))
 modes=list(range(n))+[n+e for e in present]
 free=np.array(list(itertools.product(range(-L,L+1),repeat=len(chords))),int).reshape(-1,len(chords)) if chords else np.zeros((1,0),int)
 cycle=np.zeros((q,len(chords)),int)
 for a,e in enumerate(chords):
  f,_=solve_tree(-B[:,e],q,elim);f[e]=1;cycle[:,a]=f
 offsets=free@cycle.T
 for occupied in itertools.combinations(modes,q):
  bits=inactive+sum(1<<i for i in occupied)
  nv=np.array([(bits>>i)&1 for i in range(n)],int)
  holes=np.array([1-((bits>>(n+e))&1) for e in present],int)
  rhs=nv.copy()
  for a,e in enumerate(present):rhs[all_edges[e][1]]-=holes[a]
  f,residual=solve_tree(rhs,q,elim)
  if np.any(residual):continue
  flows=offsets+f;tails=holes-flows
  keep=np.all((abs(flows)<=L)&(abs(tails)<=L),axis=1)
  for left,right in zip(flows[keep],tails[keep]):
   E=[0]*(2*M)
   for a,e in enumerate(present):E[2*e]=int(left[a]);E[2*e+1]=int(right[a])
   yield (mask,bits,tuple(E))

def hop(bits,dest,src):
 if not(bits>>src&1) or bits>>dest&1:return None
 b=bits^(1<<src);sign=(-1)**((bits&((1<<src)-1)).bit_count()+(b&((1<<dest)-1)).bit_count())
 return b|(1<<dest),sign

def gauss(n,state):
 mask,bits,E=state;edges=list(itertools.combinations(range(n),2));M=len(edges)
 Gv=np.array([-((bits>>i)&1) for i in range(n)],int);Ge=[]
 for e,(i,j) in enumerate(edges):
  Gv[i]+=E[2*e];Gv[j]+=E[2*e+1];Ge.append(-E[2*e]-E[2*e+1]-((bits>>(n+e))&1)+1)
  if not(mask>>e&1):assert E[2*e]==E[2*e+1]==0 and bits>>(n+e)&1
 assert bits.bit_count()==M
 return np.r_[Gv,Ge]

def build(n,q,L,masks=None,rewire=True):
 edges=list(itertools.combinations(range(n),2));M=len(edges)
 if masks is None:masks=[sum(1<<e for e in es) for es in itertools.combinations(range(M),q)]
 states=[s for mask in masks for s in graph_basis(n,mask,L)];lookup={s:a for a,s in enumerate(states)}
 assert len(lookup)==len(states)
 rr=array('i');cc=array('i');vv=array('d');degrees=np.zeros(len(states),int)
 for a,(mask,bits,E) in enumerate(states):
  assert not np.any(gauss(n,(mask,bits,E)))
  present=[e for e in range(M) if mask>>e&1]
  diagonal=COST*q+2*MASS*(bits&((1<<n)-1)).bit_count()+ELECTRIC/2*sum(x*x for x in E)
  for e in present:
   i,j=edges[e]
   for leg,vertex,orient in [(0,i,1),(1,j,-1)]:
    for dest,src,shift in [(vertex,n+e,1),(n+e,vertex,-1)]:
     h=hop(bits,dest,src)
     if h is None:continue
     newbits,sign=h;Ef=list(E);Ef[2*e+leg]+=shift
     if abs(Ef[2*e+leg])>L:continue
     key=(mask,newbits,tuple(Ef));b=lookup.get(key)
     assert b is not None,('gauge hop missing',key)
     rr.append(b);cc.append(a);vv.append(HOP*orient*sign)
  if rewire:
   for e in present:
    if E[2*e] or E[2*e+1] or not(bits>>(n+e)&1):continue
    for f in range(M):
     if mask>>f&1:continue
     key=(mask^(1<<e)^(1<<f),bits,E);b=lookup.get(key)
     if b is not None:
      degrees[a]+=1;rr.append(b);cc.append(a);vv.append(-GAMMA/M)
   diagonal+=GAMMA/M*degrees[a]
  rr.append(a);cc.append(a);vv.append(diagonal)
 H=coo_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),shape=(len(states),len(states))).tocsr()
 assert np.linalg.norm((H-H.T).data)<1e-12
 return H,states,lookup,degrees

def wilson_vector(n,cycle,state):
 mask,bits,E=state;edges=list(itertools.combinations(range(n),2));edgeid={e:i for i,e in enumerate(edges)};Ef=list(E)
 for i,j in zip(cycle,cycle[1:]+cycle[:1]):
  e=edgeid[tuple(sorted((i,j)))]
  if not(mask>>e&1):return None
  # U_ie U_je^dagger, ordered around the loop: divergence remains zero.
  if i<j:Ef[2*e]+=1;Ef[2*e+1]-=1
  else:Ef[2*e]-=1;Ef[2*e+1]+=1
 return mask,bits,tuple(Ef)
