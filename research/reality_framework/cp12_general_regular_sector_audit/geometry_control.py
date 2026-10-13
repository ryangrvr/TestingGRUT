#!/usr/bin/env python3
"""Optional graph metric adversary. Requires networkx (tested 3.6.1)."""
import json
from pathlib import Path
import networkx as nx

L=7
N=L**3
D=6

def balls(g,origin=0):
    d=nx.single_source_shortest_path_length(g,origin)
    return [sum(x<=r for x in d.values()) for r in range(7)]

def sample(g,label):
    assert g.number_of_nodes()==N
    assert g.number_of_edges()==D*N//2
    assert all(x==D for _,x in g.degree())
    assert nx.is_connected(g)
    return dict(label=label,N=N,Q=g.number_of_edges(),k=D,
                W=N*D*(D-1)//2,connected=True,
                ball_sizes_radius_0_to_6=balls(g),diameter=nx.diameter(g))

lattice=nx.convert_node_labels_to_integers(nx.grid_graph(dim=[range(L)]*3,periodic=True))
rows=[sample(lattice,"periodic_cubic_7x7x7")]
assert rows[0]['ball_sizes_radius_0_to_6']==[1,7,25,63,123,195,263]
assert rows[0]['diameter']==9
for seed in range(1,6):
    rows.append(sample(nx.random_regular_graph(D,N,seed=seed),f"random_regular_seed_{seed}"))
Path(__file__).with_name("geometry_results.json").write_text(json.dumps(rows,indent=2)+"\n")
print(json.dumps(rows,indent=2))
print("CP12_SPATIAL_METRIC_ADVERSARY_PASS")
