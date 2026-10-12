"""Independently diagonalize complete periodic hosts and compare the entire
finite-host edge-cluster sum. These exact finite remainders are NOT bounds
on the square lattice's infinite-host remainder.
"""
from pathlib import Path
import sys,json,itertools
import networkx as nx
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build

def record(g,rows):
 for r in rows:
  h=nx.Graph();h.add_nodes_from(range(r['vertices']));h.add_edges_from(r['pairs'])
  if len(h)==len(g) and h.number_of_edges()==g.number_of_edges() and nx.is_isomorphic(g,h):return r
 raise AssertionError('host cluster not present')

def run():
 output=[]
 for L in [3,4]:
  rows=json.loads((P.parent/'cp12_frozen_cluster_expansion'/f'clusters_k8_L{L}.json').read_text())['rows']
  for n in [4,6,8]:
   edges=list(itertools.combinations(range(n),2));M=len(edges);active=[edges.index(tuple(sorted((i,(i+1)%n)))) for i in range(n)];mask=sum(1<<e for e in active)
   H,states,lookup,_=build(n,n,L,masks=[mask],rewire=False);del lookup
   rem=np.array([sum((bits>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for _,bits,E in states])
   H=H+diags(.1*rem);E,V=eigsh(H,k=1,which='SA',tol=1e-13,v0=np.random.default_rng(120128+n).normal(size=len(states)));exact=float(E[0]/n-.6)
   r=record(nx.cycle_graph(n),rows);assert abs(exact-r['relative_energy']/n)<1e-11
   partial=[]
   for k in range(1,n):
    p1=record(nx.path_graph(k+1),rows)['relative_energy'];p0=record(nx.path_graph(k),rows)['relative_energy'] if k>1 else 0.
    approx=p1-p0;partial.append(dict(edge_order=k,partial_relative_per_vertex=approx,exact_remainder=exact-approx))
   # Only path subclusters plus the complete cycle occur. At the final order
   # their weights reconstruct the direct periodic host eigenenergy.
   reconstruction=partial[-1]['partial_relative_per_vertex']+r['weight']/n
   assert abs(reconstruction-exact)<1e-11
   out=dict(N=n,L=L,a=.1,hopping=1,dimension=len(states),exact_relative_per_vertex=exact,reconstructed_per_vertex=reconstruction,residual=float(np.linalg.norm(H@V[:,0]-E[0]*V[:,0])),partials=partial,scope='Complete periodic graph block; off-diagonal changes of host geometry excluded.')
   output.append(out);print(json.dumps({k:v for k,v in out.items() if k!='partials'}),flush=True)
 (P/'finite_host_tail_results.json').write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':run()
