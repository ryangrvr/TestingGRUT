"""Independent commutator controls for the CONDITIONAL local Hamiltonian.

The hosts and their incidence metric are supplied. Even fermionic probes
are gauge invariant. Retain the actual graph-compressed rewiring diagonal.
Onsite electric terms are unbounded in the law; hopping interactions have
uniform norm <=1. The analytical chain estimate is cutoff independent.
"""
from pathlib import Path
from fractions import Fraction
import itertools
import json
import math
import sys

import networkx as nx
import numpy as np
from scipy.sparse import diags

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def chain_tail(x, distance):
    # Positive-term summation avoids catastrophic subtraction at small x.
    if x==0:return 0.
    term=x**distance/math.factorial(distance)
    total=term;n=distance
    while True:
        n+=1;term*=x/n;total+=term
        if term < max(total*1e-16,1e-300):break
    return total


def run(n,edges,L,name):
    pairs=list(itertools.combinations(range(n),2));M=len(pairs);q=len(edges)
    active=[pairs.index(tuple(sorted(e))) for e in edges]
    mask=sum(1 << e for e in active)
    H,states,_,_=build(n,q,L,masks=[mask],rewire=False)
    a=Fraction(M-q,10*M)
    rem=np.array([sum((bits >> (n+e) & 1) and E[2*e]==E[2*e+1]==0
                     for e in active) for _,bits,E in states])
    H=H+diags(float(a)*rem)
    bits=np.array([s[1] for s in states]);coo=H.tocoo()
    e01=pairs.index((0,1))
    select=np.array([int(bits[r]^bits[c])==((1 << 0)|(1 << (n+e01)))
                     for r,c in zip(coo.row,coo.col)])
    from scipy.sparse import coo_matrix
    B=coo_matrix((coo.data[select],(coo.row[select],coo.col[select])),shape=H.shape).toarray()
    assert abs(np.linalg.norm(B,2)-1)<1e-10
    graph=nx.Graph();graph.add_nodes_from(range(n))
    for e in active:
        i,j=pairs[e];graph.add_edges_from([(i,n+e),(j,n+e)])
    D=max(dict(graph.degree()).values())
    source_support={0,n+e01}
    dense=H.toarray();E,W=np.linalg.eigh(dense)
    denominator=math.lcm(8,2,a.denominator)
    scaled=np.rint(dense*denominator).astype(np.int64)
    assert np.max(abs(scaled/denominator-dense))<1e-12
    Bint=np.rint(B).astype(np.int64)
    rows=[];zero_orders=[]
    for detector in range(1,n):
        A=np.diag((bits >> detector & 1).astype(float))
        distance=min(nx.shortest_path_length(graph,x,detector) for x in source_support)
        # Exact integer nested-commutator zeros independently test support
        # growth before the first chain capable of reaching the source.
        nested=A.astype(np.int64)
        for order in range(distance):
            comm=nested @ Bint-Bint @ nested
            assert not np.any(comm)
            zero_orders.append(dict(detector=detector,distance=distance,order=order))
            nested=scaled @ nested-nested @ scaled
        for T in [.01,.05,.1,.4,1.]:
            U=(W*np.exp(-1j*T*E)) @ W.T
            At=U.conj().T @ A @ U
            norm=float(np.linalg.norm(At @ B-B @ At,2))
            bound=2*D/(D+1)*chain_tail(2*(D+1)*T,distance)
            assert norm <= bound+1e-11
            rows.append(dict(detector=detector,T=T,incidence_distance=distance,
                             actual_commutator_norm=norm,chain_bound=min(2.,bound)))
    return dict(name=name,N=n,Q=q,L=L,dimension=len(states),a=str(a),
                maximum_incidence_degree=D,source_norm=1,integer_nested_zero_controls=zero_orders,
                rows=rows,scope="Conditional gauge-invariant fermionic propagation on a supplied host; no relativistic cone or geometry selection.")


if __name__=="__main__":
    out=dict(status="CONDITIONAL_LIGHT_CONE_CONTROL_PASS",controls=[
        run(4,[(0,1),(1,2),(2,3)],3,"path"),
        run(4,[(0,1),(1,2),(2,3),(3,0)],3,"cycle_L3"),
        run(4,[(0,1),(1,2),(2,3),(3,0)],4,"cycle_L4")])
    (P/"light_cone_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
