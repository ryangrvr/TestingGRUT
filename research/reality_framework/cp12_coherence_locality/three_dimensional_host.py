"""Supplied 3D diamond-type host and macroscopic prepared-response window.

This generates a BENCHMARK host; it does not derive dimension or select
geometry. Large rows evaluate exact analytical bounds, not full-H solves.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import math

import networkx as nx
import numpy as np
from scipy.sparse.linalg import eigsh

P=Path(__file__).resolve().parent


def host(L):
    graph=nx.Graph()
    def index(x,y,z,s):return 2*((x % L)*L*L+(y % L)*L+(z % L))+s
    for x in range(L):
        for y in range(L):
            for z in range(L):
                a=index(x,y,z,0)
                for dx,dy,dz in [(0,0,0),(-1,0,0),(0,-1,0),(0,0,-1)]:
                    graph.add_edge(a,index(x+dx,y+dy,z+dz,1))
    assert graph.number_of_nodes()==2*L**3 and graph.number_of_edges()==4*L**3
    assert nx.is_connected(graph) and set(dict(graph.degree()).values())=={4}
    return graph,index


def spectrum(L):
    k=2*np.pi*np.arange(L)/L
    x,y,z=np.meshgrid(k,k,k,indexing="ij")
    f=1+np.exp(1j*x)+np.exp(1j*y)+np.exp(1j*z)
    lam=np.r_[4-abs(f).ravel(),4+abs(f).ravel()]
    lam=np.maximum(lam,0)
    assert np.count_nonzero(lam<1e-12)==1
    return np.sort(lam)


def run():
    finite=[]
    for L in [3,4,6,8]:
        graph,index=host(L);N=len(graph);r=math.isqrt(L)
        extra=[(index(2*j,0,0,0),index(2*j+1,0,0,0)) for j in range(r)]
        expanded=graph.copy();expanded.add_edges_from(extra)
        assert expanded.number_of_edges()==2*N+r
        assert max(dict(expanded.degree()).values())==5
        lap=nx.laplacian_matrix(graph).astype(float)
        eig=np.sort(eigsh(lap,k=3,which="SA",return_eigenvectors=False,tol=1e-11))
        bloch=spectrum(L)
        assert abs(eig[1]-bloch[1])<1e-9
        finite.append(dict(L=L,N=N,host_edges=2*N,pool_edges=2*N+r,
                           pool_maximum_degree=5,sparse_gap=float(eig[1]),
                           independently_computed_bloch_gap=float(bloch[1])))
    diffusion=[]
    for L in [16,32,64]:
        eig=spectrum(L)
        for s in [5.,20.,50.]:
            w=np.exp(-s*eig)
            ds=2*s*float(eig @ w/w.sum())
            diffusion.append(dict(L=L,N=2*L**3,diffusion_time=s,spectral_dimension=ds))
    bounds=[]
    for L in [4,8,16,32,64,128,256]:
        N=2*L**3;M=N*(N-1)//2;Q=2*N
        for r in [0,1,math.isqrt(L)]:
            S=Q+r;T=L
            b2=F(Q*(M-S)*(r+1),100*M*M)
            vin=F(Q*r,10*M)
            error=4*b2*T*T+4*vin*T
            mean_energy=F(Q,2)+F(Q*(M-S),10*M)
            variance=2*Q+b2
            bounds.append(dict(L=L,N=N,Q=Q,S=S,vacancies=r,observation_time=T,
                               exact_response_error_bound=str(error),response_error_bound=float(error),
                               prepared_mean_energy=str(mean_energy),energy_per_vertex=float(mean_energy/N),
                               prepared_energy_variance=str(variance),
                               relative_energy_standard_deviation=math.sqrt(float(variance))/float(mean_energy)))
    out=dict(status="SUPPLIED_3D_MACROSCOPIC_WINDOW_CHECK_PASS",host_controls=finite,
             supplied_host_diffusion=diffusion,analytical_response_and_preparation_bounds=bounds,
             scope="A supplied 3D bounded-degree benchmark. Macroscopic bounds follow from full CP12; no large-H simulation, dimensional selection, clock derivation, or gravity.")
    (P/"three_dimensional_host_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":run()
