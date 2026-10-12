#!/usr/bin/env python3
"""Independent even-width tube comparator for frozen CP-12 cluster dataset.

Uses only the exact finite-cluster energies from owner's *pinned*
clusters_k8_L3.json and clusters_k8_L4.json (frozen t=1, a=.1).
Enumerates all non-induced injective rooted embeddings, and computes graph
automorphisms independently with the standard library. No new Hamiltonian
diagonalization, full graph mixing, or cluster tail is claimed.

Run: python research/reality_framework/cp12_even_width_cluster_control/check.py
"""
import json,math
from pathlib import Path
from collections import deque

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/"cp12_frozen_cluster_expansion"

def host(w,h):
    assert w>3 and h>3
    return [tuple(sorted({y*w+(x+1)%w, y*w+(x+w-1)%w,
                          ((y+1)%h)*w+x, ((y+h-1)%h)*w+x}))
            for y in range(h) for x in range(w)]

def adjacency(row):
    n=row["vertices"]
    out=[set() for _ in range(n)]
    for a,b in row["pairs"]:
        out[a].add(b)
        out[b].add(a)
    return out

def rooted_embeddings(row, host_neighbors):
    graph=adjacency(row);n=len(graph)
    order=[0];seen={0}
    for u in order:
        for v in sorted(graph[u]):
            if v not in seen:
                seen.add(v);order.append(v)
    assert len(order)==n
    mapping={0:0};used={0}
    def rec(i):
        if i==n:return 1
        u=order[i]
        previous=[mapping[v] for v in graph[u] if v in mapping]
        choices=set(host_neighbors[previous[0]])
        for v in previous[1:]:
            choices.intersection_update(host_neighbors[v])
        result=0
        for c in choices-used:
            mapping[u]=c;used.add(c)
            result+=rec(i+1)
            used.remove(c);del mapping[u]
        return result
    return rec(1)

def aut_count(row):
    graph=adjacency(row);n=len(graph)
    order=sorted(range(n),key=lambda i:(-len(graph[i]),i))
    mapping={};used=set()
    def rec(i):
        if i==n:return 1
        u=order[i];total=0
        for v in range(n):
            if v in used or len(graph[v])!=len(graph[u]):continue
            if any((x in graph[u])!=(mapping[x] in graph[v])
                   for x in mapping):continue
            used.add(v);mapping[u]=v
            total+=rec(i+1)
            used.remove(v);del mapping[u]
        return total
    return rec(0)

def paths_and_cycle(rows):
    out={}
    for row in rows:
        adj=adjacency(row);degrees=sorted(map(len,adj))
        e=row["edges"];v=row["vertices"]
        if e==8 and v==8 and degrees==[2]*8:out["C8"]=row
        if e==7 and v==8 and degrees==[1,1]+[2]*6:out["P8"]=row
        if e==8 and v==9 and degrees==[1,1]+[2]*7:out["P9"]=row
    assert set(out)=={"C8","P8","P9"}
    return out

def main():
    square=host(18,18)
    even_tube=host(8,36)
    distant_tube=host(10,36)
    all_results=[]
    for cutoff in (3,4):
        rows=json.loads((SOURCE/f"clusters_k8_L{cutoff}.json").read_text())["rows"]
        deviations=[]
        for row in rows:
            aut=aut_count(row)
            sq=rooted_embeddings(row,square)
            tube=rooted_embeddings(row,even_tube)
            wide=rooted_embeddings(row,distant_tube)
            assert sq/aut==row["square_embeddings_per_vertex"],row["key"]
            assert wide==sq,(row["key"],"width10 should first differ at order >=10")
            if tube!=sq:
                deviations.append(dict(
                    key=row["key"],edges=row["edges"],vertices=row["vertices"],
                    aut=aut,rooted_square=sq,rooted_width8=tube,
                    embedding_density_difference=(tube-sq)/aut,
                    cluster_weight=row["weight"],
                    energy=row["relative_energy"]))
        assert len(deviations)==2, deviations
        deviations.sort(key=lambda r:r["vertices"])
        assert [(x["edges"],x["vertices"],x["embedding_density_difference"])
           for x in deviations]==[(8,8,1/8),(8,9,-1)],deviations
        target=paths_and_cycle(rows)
        direct=(target["C8"]["relative_energy"]/8
                -target["P9"]["relative_energy"]
                +target["P8"]["relative_energy"])
        weighted=sum(x["embedding_density_difference"]*x["cluster_weight"]
                     for x in deviations)
        assert direct>0
        assert abs(direct-weighted)<1e-13
        all_results.append(dict(
            cutoff_L=cutoff,clusters_tested=len(rows),
            host_comparison="square 18x18 versus even width8 x length36",
            control="even width10 x length36 shows no differences <=8 edges",
            contributing_clusters=deviations,
            E_cycle8=target["C8"]["relative_energy"],
            E_path8=target["P8"]["relative_energy"],
            E_path9=target["P9"]["relative_energy"],
            exact_embedding_weight_method=weighted,
            independent_energy_telescoping_method=direct,
            formula="delta_e8=w(C8)/8-w(P9)=E(C8)/8-[E(P9)-E(P8)]",
            delta_e8_sign="positive: extended square lower AT EDGE ORDER 8 only",
            no_thermodynamic_remainder_bound=True))
    assert abs(all_results[0]["independent_energy_telescoping_method"]
              -all_results[1]["independent_energy_telescoping_method"])<2e-14
    n=200;M=n*(n-1)//2;Q=2*n;gamma=.1
    # Exact Schur/row-sum upper bound for graph-changing off-diagonal
    # CP12 operator on the entire fixed-Q constrained joint Hilbert space:
    # any neutral state has no more than Q(M-Q) outgoing links, and
    # each transition matrix element has magnitude gamma/M.
    offdiag_norm_bound=gamma*Q*(M-Q)/M
    conclusion=dict(
      model="Frozen CP12 t=1 graph-block expansion, no full graph mixing",
      geometry="width8 bipartite tube versus even square-lattice host",
      results=all_results,
      comparison_width9_at_9_edges_L4=1.0414431046e-9,
      comparison_width9_through_10_edges_L4=-5.85430687e-9,
      graph_mixing_operator_schur_bound={
        "N":n,"Q":Q,"M":M,"Gamma":gamma,
        "upper_full_offdiag_operator_norm":offdiag_norm_bound,
        "upper_full_offdiag_norm_per_vertex":offdiag_norm_bound/n,
        "general_formula":"Gamma Q(M-Q)/M, fixed-Q full graph-changing block",
        "caveat":"Vastly loose maximum-norm control; not actual mixing energy or geometry comparison."},
      caution=[
       "The width9 tube is non-bipartite; the width8 tube remains bipartite.",
       "A positive first boundary coefficient is NOT a final thermodynamic energy ordering.",
       "All source finite-block energies from the owner's pinned 27010b8 dataset; not independently diagonalized here.",
       "Finite electric regulator agreement is not an edge-cluster remainder certificate.",
       "Quantum off-diagonal geometry mixing is not evaluated in linked graph-block energies."])
    (HERE/"results.json").write_text(json.dumps(conclusion,indent=2)+"\n")
    print(json.dumps({"checks":"PASS","differences":all_results},indent=2))
    print("CP12_EVEN_WIDTH_CLUSTER_CONTROL_PASS")
if __name__=="__main__":main()
