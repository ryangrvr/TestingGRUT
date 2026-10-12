"""Full frozen CP12 ground, graph statistics, and physical record back-action.

The controlled record is an ideal measurement isometry, not an autonomous
apparatus. Finite electric cutoffs are explicit. CP12 itself is unchanged.
"""
from pathlib import Path
import argparse
from collections import defaultdict
import itertools
import json
import math
import sys

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def run(n, q, L):
    H, states, _, _ = build(n, q, L)
    energy, vectors = eigsh(H, k=1, which="SA", tol=1e-12,
                           v0=np.random.default_rng(121003).normal(size=H.shape[0]))
    E0, psi = float(energy[0]), vectors[:, 0]
    masks = np.array([s[0] for s in states])
    M = n*(n-1)//2
    coo = H.tocoo()
    cross = masks[coo.row] != masks[coo.col]
    V = coo_matrix((coo.data[cross], (coo.row[cross], coo.col[cross])),
                   shape=H.shape).tocsr()
    coherence_energy = -float(psi @ (V @ psi))
    normV = .1*q*(M-q)/M
    probabilities = defaultdict(float)
    for mask, p in zip(masks, psi**2):
        probabilities[int(mask)] += float(p)
    pairs = list(itertools.combinations(range(n), 2))
    graph_rows = []
    averaged = defaultdict(float)
    for mask, p in sorted(probabilities.items()):
        graph = nx.Graph()
        graph.add_nodes_from(range(n))
        graph.add_edges_from(pair for e, pair in enumerate(pairs) if mask >> e & 1)
        deg = [graph.degree(i) for i in range(n)]
        components = sorted((len(c) for c in nx.connected_components(graph)), reverse=True)
        row = dict(mask=mask, probability=p, degree=deg,
                   components=components, isolated=deg.count(0),
                   connected=len(components)==1,
                   triangles=sum(nx.triangles(graph).values())//3)
        graph_rows.append(row)
        averaged["connected_probability"] += p*row["connected"]
        averaged["isolated_vertices"] += p*row["isolated"]
        averaged["largest_component_fraction"] += p*components[0]/n
        averaged["triangle_count"] += p*row["triangles"]
        averaged["degree_variance"] += p*sum((d-2*q/n)**2 for d in deg)/n
    assert abs(sum(probabilities.values())-1) < 1e-10

    # Each rewiring changes two edge bits. Sum of the energy costs of
    # separately recording each edge is exactly twice the full graph cost.
    contributions = coo.data[cross]*psi[coo.row[cross]]*psi[coo.col[cross]]
    changes = masks[coo.row[cross]] ^ masks[coo.col[cross]]
    assert all(int(x).bit_count()==2 for x in changes)
    edge_costs = [-float(np.sum(contributions[(changes >> e & 1).astype(bool)]))
                  for e in range(M)]
    assert min(edge_costs) >= -1e-10
    assert abs(sum(edge_costs)-2*coherence_energy) < 1e-10

    # Independently build the actual system x reference pure state for
    # small hosts. Ancilla bit 0 is |0>; present-edge bit is c|0>+s|1>.
    record_rows = []
    if n <= 4:
        for theta in [.1, .3, math.pi/4, math.pi/2]:
            c, s = math.cos(theta), math.sin(theta)
            records = {}
            for mask in probabilities:
                v = np.ones(1)
                for e in range(M):
                    v = np.kron(v, [c, s] if mask >> e & 1 else [1., 0.])
                records[mask] = v
            joint = psi[:, None]*np.array([records[int(mask)] for mask in masks])
            assert abs(np.linalg.norm(joint)-1) < 1e-10
            actual = float(np.sum(joint*(H @ joint))) - E0
            predicted = s*s*coherence_energy
            assert abs(actual-predicted) < 1e-10
            record_rows.append(dict(theta=theta, equal_prior_edge_distinguishability=abs(s),
                                    actual_system_energy_increase=actual,
                                    predicted_energy_increase=predicted))
    out = dict(N=n,Q=q,L=L,dimension=len(states),graph_count=len(probabilities),
               full_energy=E0,residual=float(np.linalg.norm(H @ psi-E0*psi)),
               electric_boundary_weight=float(sum(p for p,s in zip(psi**2,states)
                                                  if max(map(abs,s[2]))==L)),
               offdiagonal_norm=normV,coherence_energy=coherence_energy,
               normalized_rewiring_coherence=coherence_energy/normV if normV else 0.,
               completely_graph_dephased_energy=E0+coherence_energy,
               edge_record_costs=edge_costs,record_controls=record_rows,
               geometry_statistics=dict(averaged),graph_probabilities=graph_rows,
               scope="Complete finite regulated Hamiltonian. Record statistics are an ideal isometry; no apparatus energy or spatial phase is inferred.")
    (P/f"full_ground_N{n}_Q{q}_L{L}.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="graph_probabilities"},indent=2),flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=4)
    ap.add_argument("--q",type=int,default=3)
    ap.add_argument("--L",type=int,default=3)
    args = ap.parse_args()
    run(args.n,args.q,args.L)
