"""Exact inclusion-exclusion expression for the tenth-edge difference.
Only connected components matter because the graph-block ground energy is additive.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools,json,networkx as nx
P=Path(__file__).resolve().parent
sources=json.loads((P/'tenth_boundary_results.json').read_text())['different']
registry=[];coeff={}
def identify(g):
 for i,h in enumerate(registry):
  if len(h)==len(g) and h.number_of_edges()==g.number_of_edges() and nx.is_isomorphic(g,h):return i
 registry.append(nx.convert_node_labels_to_integers(g));return len(registry)-1
for r in sources:
 edges=[tuple(e) for e in r['pairs']];q=len(edges);d=F(r['density_difference']).limit_denominator()
 for mask in range(1<<q):
  h=nx.Graph();h.add_edges_from(e for j,e in enumerate(edges) if mask>>j&1)
  c=d*((-1)**(q-mask.bit_count()))
  for comp in nx.connected_components(h):
   sub=h.subgraph(comp).copy();i=identify(sub);coeff[i]=coeff.get(i,F(0))+c
out=[]
for i,c in coeff.items():
 if c:
  g=registry[i];out.append(dict(coefficient=str(c),N=len(g),Q=g.number_of_edges(),pairs=list(g.edges())))
assert sum(F(r['coefficient'])*r['Q'] for r in out)==0
(P/'tenth_energy_expression.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
