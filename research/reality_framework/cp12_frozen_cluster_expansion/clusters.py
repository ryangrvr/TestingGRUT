"""Edge linked-cluster ground-state expansion at frozen hopping, not t-series.
Bipartite max-degree-four clusters; the atlas supplies orders through six edges,
and exhaustive edge extensions generate higher orders. a=.1 is the inherited
fixed-Q thermodynamic graph-compression coefficient.
"""
from pathlib import Path
import sys,itertools,json,argparse,time
import networkx as nx
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,COST,HOP

def catalog(k):
 base=sorted([nx.convert_node_labels_to_integers(g) for g in nx.graph_atlas_g() if len(g)>1 and nx.is_connected(g) and nx.is_bipartite(g) and g.number_of_edges()<=min(k,6) and max(dict(g.degree()).values())<=4],key=lambda g:(g.number_of_edges(),len(g),nx.to_graph6_bytes(g,header=False)))
 for q in range(7,k+1):
  buckets={};new=[]
  for old in [g for g in base if g.number_of_edges()==q-1]:
   n=len(old)
   for i in range(n):
    if old.degree(i)>=4:continue
    for j in range(i+1,n+1):
     if j<n and (old.has_edge(i,j) or old.degree(j)>=4):continue
     h=old.copy();h.add_edge(i,j)
     if not nx.is_bipartite(h):continue
     code=(len(h),tuple(sorted(dict(h.degree()).values())),nx.weisfeiler_lehman_graph_hash(h))
     found=buckets.setdefault(code,[])
     if any(nx.is_isomorphic(h,g) for g in found):continue
     found.append(h);new.append(h)
  base.extend(new)
 return base

def key(g):return nx.to_graph6_bytes(g,header=False).decode().strip()

def subcounts(g,graphs):
 counts={};edges=list(g.edges())
 for q in range(1,len(edges)):
  for subset in itertools.combinations(edges,q):
   h=nx.Graph();h.add_edges_from(subset)
   if not nx.is_connected(h):continue
   for c in graphs:
    if c.number_of_edges()==q and len(c)==len(h) and nx.is_isomorphic(c,h):
     counts[key(c)]=counts.get(key(c),0)+1;break
   else:raise AssertionError('missing subcluster')
 return counts

def rooted_count(g,host,root):
 # Injective non-induced embeddings with cluster vertex0 fixed to host root.
 order=list(nx.bfs_tree(g,0));assigned={0:root};used={root}
 def rec(i):
  if i==len(order):return 1
  u=order[i];neighbors=[assigned[v] for v in g[u] if v in assigned]
  assert neighbors
  choices=set(host[neighbors[0]])
  for v in neighbors[1:]:choices.intersection_update(host[v])
  answer=0
  for v in choices-used:
   assigned[u]=v;used.add(v);answer+=rec(i+1);used.remove(v);del assigned[u]
  return answer
 return rec(1)

def run(L,maxedges,a,embedded_only=False):
 graphs=catalog(maxedges);rows=[];weights={}
 square=nx.convert_node_labels_to_integers(nx.grid_2d_graph(18,18,periodic=True))
 thin=nx.convert_node_labels_to_integers(nx.grid_2d_graph(9,36,periodic=True))
 # Root has4 branches, descendants3; radius maxedges suffices.
 tree=nx.Graph();tree.add_node(0);front=[0];nextid=1
 for depth in range(maxedges):
  new=[]
  for u in front:
   for _ in range(4 if depth==0 else 3):tree.add_edge(u,nextid);new.append(nextid);nextid+=1
  front=new
 for g in graphs:
  aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(g,g).isomorphisms_iter())
  es=rooted_count(g,square,0)/aut;et=rooted_count(g,tree,0)/aut;eth=rooted_count(g,thin,0)/aut
  assert es==eth
  if embedded_only and es==0 and et==0:continue
  n=len(g);q=g.number_of_edges();edges=list(itertools.combinations(range(n),2));M=len(edges)
  active=[edges.index(tuple(sorted(e))) for e in g.edges()];mask=sum(1<<e for e in active)
  H,states,_,_=build(n,q,L,masks=[mask],rewire=False)
  rem=np.array([sum((bits>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for _,bits,E in states])
  H=H+diags(a*rem)
  e,u=eigsh(H,k=1,which='SA',tol=2e-12,v0=np.ones(len(states)));v=u[:,0];energy=float(e[0]-(COST+a)*q)
  counts=subcounts(g,graphs);w=energy-sum(count*weights[k] for k,count in counts.items());weights[key(g)]=w
  aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(g,g).isomorphisms_iter())
  es=rooted_count(g,square,0)/aut;et=rooted_count(g,tree,0)/aut;eth=rooted_count(g,thin,0)/aut
  assert es==eth
  row=dict(key=key(g),vertices=n,edges=q,pairs=list(g.edges()),dimension=len(states),relative_energy=energy,weight=w,subclusters=counts,square_embeddings_per_vertex=es,thin_embeddings_per_vertex=eth,tree_embeddings_per_vertex=et,residual=float(np.linalg.norm(H@v-e[0]*v)),cutoff_boundary=float(sum(abs(v[i])**2 for i,s in enumerate(states) if max(map(abs,s[2]))==L)))
  assert row['residual']<1e-9;rows.append(row);print(json.dumps(row),flush=True)
 sums=[]
 for order in range(1,maxedges+1):
  sums.append(dict(order=order,square_relative_per_vertex=sum(r['weight']*r['square_embeddings_per_vertex'] for r in rows if r['edges']<=order),tree_relative_per_vertex=sum(r['weight']*r['tree_embeddings_per_vertex'] for r in rows if r['edges']<=order)))
 out=dict(L=L,a=a,hopping=HOP,maxedges=maxedges,rows=rows,partial_sums=sums,scope='Exact finite-cluster ground energies, approximate thermodynamic sum without a remainder bound. No off-diagonal graph rewiring.')
 (P/(f'clusters_L{L}.json' if maxedges==6 else f'clusters_k{maxedges}_L{L}.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(sums,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=3);p.add_argument('--maxedges',type=int,default=6);p.add_argument('--a',type=float,default=.1);p.add_argument('--embedded-only',action='store_true');a=p.parse_args();run(a.L,a.maxedges,a.a,a.embedded_only)
