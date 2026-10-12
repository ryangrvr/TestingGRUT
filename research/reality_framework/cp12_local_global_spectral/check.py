#!/usr/bin/env python3
"""Independent CP12 local-geometry nonselection comparator.

No CP12 quantum dynamics is imported. Two connected graphs match in
degree, triangle count, 4-cycle count, and EVERY rooted radius-3 ball,
yet have different macroscopic Laplacian scaling and diffusion dimension.
Includes a separate sampled random 4-regular graph control.

Requires numpy scipy networkx. Run: OPENBLAS_NUM_THREADS=1 python check.py
"""
import json, math
from pathlib import Path
import networkx as nx
import numpy as np
from scipy.sparse import eye,csr_matrix
from scipy.sparse.linalg import eigsh

WIDTH,RADIUS=9,3

def torus(a,b):
    return nx.convert_node_labels_to_integers(
        nx.grid_2d_graph(a,b,periodic=True))

def rooted_ball(G,r):
    ids=nx.single_source_shortest_path_length(G,0,cutoff=r)
    return G.subgraph(ids).copy()

def count_4cycles(G):
    A=csr_matrix(nx.to_scipy_sparse_array(G,format="csr",dtype=np.int64))
    C=(A@A).tocoo()
    # Opposite pairs count every 4-cycle twice.
    x=sum(int(v)*(int(v)-1)//2 for i,j,v in zip(
        C.row,C.col,C.data) if i<j)
    assert x%2==0
    return x//2

def graph_features(G):
    n=len(G)
    assert G.number_of_edges()==2*n
    assert nx.is_connected(G)
    assert all(d==4 for _,d in G.degree())
    W=sum(d*(d-1)//2 for _,d in G.degree())
    assert W==6*n
    return dict(N=n,Q=2*n,W=W,connected=True,
        triangles=sum(nx.triangles(G).values())//3,
        squares=count_4cycles(G),
        radius3_vertices=len(rooted_ball(G,RADIUS)))

def exact_gap(a,b):
    return 4*min(math.sin(math.pi/a)**2,
                 math.sin(math.pi/b)**2)

def heat_return(a,b,time):
    def one(m):
        lam=4*np.sin(np.pi*np.arange(m)/m)**2
        return float(np.exp(-time*lam).mean())
    return one(a)*one(b)

def exact_spectral_dimension(a,b,time):
    # d_s(T)=-2 d log P(T)/d log T =2T <lambda>_heat.
    def mean_lambda(m):
        lam=4*np.sin(np.pi*np.arange(m)/m)**2
        weight=np.exp(-time*lam)
        return float(weight@lam/weight.sum())
    return 2*time*(mean_lambda(a)+mean_lambda(b))

def random_gap(G):
    A=csr_matrix(nx.to_scipy_sparse_array(G,format="csr",dtype=float))
    Lap=4*eye(len(G),format="csr")-A
    # Request enough eigenmodes to handle degeneracies in the lower shell.
    E,V=eigsh(Lap,k=10,which="SM",tol=1e-10,maxiter=100000)
    order=np.argsort(E);E=E[order];V=V[:,order]
    res=max(np.linalg.norm(Lap@V[:,i]-E[i]*V[:,i])
            for i in range(len(E)))
    assert abs(E[0])<1e-8 and res<1e-7
    return float(E[1])

def main():
    exact=[];random=[]
    for k in (2,3,6):
        L=WIDTH*k;N=L*L
        sq=(L,L);thin=(WIDTH,N//WIDTH)
        A,B=torus(*sq),torus(*thin)
        fa,fb=graph_features(A),graph_features(B)
        assert fa==fb
        assert nx.is_isomorphic(rooted_ball(A,RADIUS),
                               rooted_ball(B,RADIUS))
        assert fa["triangles"]==0 and fa["squares"]==N
        gap_a,gap_b=exact_gap(*sq),exact_gap(*thin)
        dims=[dict(time=t,square=exact_spectral_dimension(*sq,t),
                    thin=exact_spectral_dimension(*thin,t))
              for t in (1.,5.,10.,20.,40.)]
        exact.append(dict(N=N,square_dimensions=sq,
            thin_dimensions=thin,invariants=fa,square_gap=gap_a,
            thin_gap=gap_b,gap_ratio=gap_a/gap_b,diffusion=dims,
            heat_return_t20_square=heat_return(*sq,20.),
            heat_return_t20_thin=heat_return(*thin,20.)))
        if k in (2,3):
            for seed in range(3):
                G=nx.random_regular_graph(4,N,seed=11222026+seed)
                f=graph_features(G)
                random.append(dict(N=N,seed=seed,gap=random_gap(G),
                    triangles=f["triangles"],squares=f["squares"],
                    radius3_vertices=f["radius3_vertices"]))
    assert exact[1]["gap_ratio"]>8
    assert exact[2]["gap_ratio"]>30
    assert exact[2]["diffusion"][3]["square"]>1.6
    assert .95<exact[2]["diffusion"][3]["thin"]<1.10
    assert min(v["gap"] for v in random if v["N"]==729)>.35
    output=dict(
        scope="Pure graph and Laplacian comparator, NOT full CP12 quantum thermal state",
        theorem="For any fixed radius r, sufficiently wide fixed strips and growing square tori have isomorphic rooted r-balls and identical sum-of-local-r-ball energies, yet differ in macroscopic spectral gap and diffusion.",
        exact_gaps="lambda2(C_a Cartesian C_b)=4 min(sin(pi/a)^2,sin(pi/b)^2)",
        asymptotics="Square gap ~4 pi^2/N. Fixed-width w strip gap ~4 pi^2 w^2/N^2.",
        data=exact,random_regular_controls=random,
        limitations=[
            "The graph geometries are mathematical hostile controls, not CP12-generated states.",
            "Full CP12 includes nonlocal graph rewiring, quantum fields and higher-order topology-sensitive dynamics; local-energy degeneracy is NOT full-CP12 degeneracy.",
            "Equal radius-r neighborhoods do not force equal global quantum spectra.",
            "The thin graph is a finite 2D torus that becomes effectively 1D at long scales with its width fixed.",
            "Random regular graph samples are finite demonstrations, not asymptotic proofs."
        ])
    Path(__file__).with_name("results.json").write_text(
        json.dumps(output,indent=2)+"\n")
    print(json.dumps({"exact":[dict(N=v["N"],square_gap=v["square_gap"],
         thin_gap=v["thin_gap"],spectral_dimension_t20=v["diffusion"][3])
         for v in exact],
         "random":random},indent=2))
    print("CP12_LOCAL_GLOBAL_SPECTRAL_NON_SELECTION_PASS")

if __name__=="__main__":main()
