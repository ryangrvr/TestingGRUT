"""Enumerate every relevant cluster topology through nine edges.
Thin width9 has no odd cycle shorter than9; the sole non-bipartite cluster
with <=9 edges that can embed is C9. All bipartite clusters are generated.
"""
from pathlib import Path
import json,networkx as nx
from collections import Counter
from clusters import catalog,rooted_count,key
P=Path(__file__).resolve().parent
A=nx.convert_node_labels_to_integers(nx.grid_2d_graph(27,27,periodic=True));B=nx.convert_node_labels_to_integers(nx.grid_2d_graph(9,81,periodic=True))
G=catalog(9)+[nx.cycle_graph(9)];different=[]
for g in G:
 x=rooted_count(g,A,0);y=rooted_count(g,B,0)
 if x!=y:
  aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(g,g).isomorphisms_iter())
  different.append(dict(key=key(g),vertices=len(g),edges=g.number_of_edges(),degrees=sorted(dict(g.degree()).values()),balanced_rooted=x,thin_rooted=y,density_difference=(y-x)/aut))
assert len(different)==2
assert sorted(r['density_difference'] for r in different)==[-1.,1/9]
assert all(r['edges']==9 for r in different)
out=dict(status='NINE_EDGE_EMBEDDING_BOUNDARY_PASS',counts=dict(Counter(g.number_of_edges() for g in G)),different=different)
(P/'embedding_boundary_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
