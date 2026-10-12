"""Locate the first correction after the width-nine linked-cluster difference."""
from pathlib import Path
import networkx as nx,json
from clusters import catalog,rooted_count,key
P=Path(__file__).resolve().parent
A=nx.convert_node_labels_to_integers(nx.grid_2d_graph(27,27,periodic=True));B=nx.convert_node_labels_to_integers(nx.grid_2d_graph(9,81,periodic=True))
leaf=nx.cycle_graph(9);leaf.add_edge(0,9)
G=[g for g in catalog(10) if g.number_of_edges()==10]+[leaf];rows=[]
for g in G:
 x=rooted_count(g,A,0);y=rooted_count(g,B,0)
 if x!=y:
  aut=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(g,g).isomorphisms_iter())
  rows.append(dict(key=key(g),vertices=len(g),edges=g.number_of_edges(),pairs=list(g.edges()),degrees=sorted(dict(g.degree()).values()),balanced_rooted=x,thin_rooted=y,density_difference=(y-x)/aut))
(P/'tenth_boundary_results.json').write_text(json.dumps(dict(catalog_count=len(G),different=rows),indent=2)+'\n');print(json.dumps(rows,indent=2))
