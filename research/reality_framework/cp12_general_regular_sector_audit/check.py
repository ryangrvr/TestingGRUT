#!/usr/bin/env python3
"""Independent fixed-degree regular-sector adversarial audit.

Exact finite simple-graph switching enumeration and analytic formula checks.
This is not a solution to the full CP12 quantum Gibbs state.
"""
from collections import Counter
from itertools import combinations
from math import comb, exp, log, factorial, lgamma, sqrt
import json
from pathlib import Path

GAMMA = 0.1

def finite_graph_switch_check(n, k):
    assert 2 <= k <= n-2 and (n*k)%2 == 0
    vertices = range(n)
    edges = list(combinations(vertices,2))
    M = len(edges)
    Q = n*k//2
    regular = 0
    successors = Counter()
    enumerated = 0
    for graph in combinations(edges,Q):
        enumerated += 1
        degrees = [0]*n
        for u,v in graph:
            degrees[u]+=1;degrees[v]+=1
        if any(d!=k for d in degrees):continue
        regular += 1
        es = set(graph)
        for u,v in graph:
            for a,b in ((u,v),(v,u)):
                for c in vertices:
                    if c==a or tuple(sorted((a,c))) in es:continue
                    assert c!=b
                    modified=frozenset(es-{(u,v)}|{tuple(sorted((a,c)))})
                    successors[modified] += 1
    assert regular > 0
    assert enumerated == comb(M,Q)
    forward=k*n*(n-1-k)
    assert sum(successors.values()) == regular*forward
    max_reverse = max(successors.values())
    assert max_reverse <= k+1
    for output in successors:
        degree=[0]*n
        for u,v in output:degree[u]+=1;degree[v]+=1
        assert sorted(degree)==[k-1]+[k]*(n-2)+[k+1]
    # Z>=Z_reg + Z_one_defect and every one-defect raises W by one.
    lower_ratio=forward/(k+1)
    exact_count_ratio=len(successors)/regular
    assert exact_count_ratio >= lower_ratio - 1e-12
    M2=Q*(M-Q)
    analytic_cut_bound=GAMMA/M*sqrt((k+1)*M2)
    return dict(n=n,k=k,Q=Q,all_graphs=enumerated,regular_graphs=regular,
                distinct_defective_successors=len(successors),
                total_switches=sum(successors.values()),
                max_reverse_multiplicity=max_reverse,
                exact_defect_count_ratio=exact_count_ratio,
                theorem_defect_ratio_lower_bound=lower_ratio,
                theoretical_sector_cut_norm_upper=analytic_cut_bound,
                p_regular_upper_at_a4=1/(1+lower_ratio*exp(-4)))

def spatial_adversary(d,L):
    assert L>=4
    n=L**d; k=2*d;Q=d*n
    assert n%(k+1)==0
    # A d-dimensional periodic cubic lattice and disjoint (k+1)-cliques:
    # both are simple k-regular graphs with Q=kN/2 and W=N*k*(k-1)/2.
    return dict(dimension=d,side=L,vertices=n,degree=k,connections=Q,
        degree_wedge_energy_index=n*comb(k,2),
        disconnected_cliques=n//(k+1),clique_size=k+1)

def log_regular_graph_asymptotic(n,k):
    assert n*k % 2 == 0
    # McKay/Wormald fixed-k asymptotic count:
    # (kn)!/((kn/2)! 2^(kn/2) (k!)^n) exp(-(k^2-1)/4) (1+o(1)).
    Q=n*k//2
    return (lgamma(n*k+1)-lgamma(Q+1)-Q*log(2)
            -n*lgamma(k+1)-(k*k-1)/4)

def main():
    cases=[(6,2),(6,3),(6,4),(7,2),(7,4)]
    finite=[finite_graph_switch_check(n,k) for n,k in cases]
    spatial=[spatial_adversary(d,L) for d,L in [(1,6),(2,5),(3,7)]]
    asym=[]
    for k in (4,6):
        for n in (100,1000,10000):
            regular_log=log_regular_graph_asymptotic(n,k)
            log_upper=lgamma(n+1)-regular_log
            asym.append(dict(n=n,k=k,log_upper_fraction_single_spatial_isomorphism_class=log_upper,
                             log_upper_over_nlogn=log_upper/(n*log(n)),
                             expected_leading_coefficient=1-k/2))
    out=dict(title="general fixed degree regularity and quantum cut adversaries",
             scope="classical degree-only Gibbs and graph-cut bound under CP12 neutral edge-swap law",
             theorem="P(k-regular) <= (1+k*N*(N-1-k)/(k+1)*exp(-beta*lambda))^-1",
             sector_coupling="b_Nk <= Gamma/M*sqrt[(k+1)*Q*(M-Q)] (fixed k: Gamma*sqrt[k(k+1)/N]*(1+o(1)))",
             dynamical_leakage="Pr(outside regular sector,t) <= min(1, b_Nk^2*t^2) for any initial state wholly in regular sector",
             comparator_examples=spatial,finite_checks=finite,
             asymptotic_count_diagnostics=asym,
             warnings=["No full quantum Gibbs theorem", "Fixed-k asymptotic enumeration cited, not proved by this script", "Graph-sector propagation bound does not select which regular graph", "Full CP12 may correlate matter/gauge fields and reshape geometry weights"])
    Path(__file__).with_name("results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"finite_switching": [{k:x[k] for k in ("n","k","regular_graphs","distinct_defective_successors","max_reverse_multiplicity")} for x in finite],"spatial_competitors":spatial},indent=2))
    print("CP12_GENERAL_REGULAR_SECTOR_CHECK_PASS")

if __name__=='__main__':main()
