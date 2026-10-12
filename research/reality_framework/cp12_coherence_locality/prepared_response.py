"""Gauge-invariant local source and detector in a coherent prepared pool.

Compare full CP12, the pool compression, and independent graph blocks.
All keep the ACTUAL rewiring diagonal. The source is a frozen hopping
operator at one vertex; no geometry-role oracle is used.
"""
from pathlib import Path
import itertools
import json
import math
import sys

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def run(host_edges, name, n=4, q=3, L=3):
    pairs=list(itertools.combinations(range(n),2));M=len(pairs)
    pool=sum(1 << pairs.index(tuple(sorted(e))) for e in host_edges)
    S=pool.bit_count();r=S-q;outside=M-S
    H,states,_,_=build(n,q,L)
    masks=np.array([s[0] for s in states]);bits=np.array([s[1] for s in states])
    allowed=np.flatnonzero((masks & ~pool)==0)
    coo=H.tocoo(); graph_change=masks[coo.row]!=masks[coo.col]
    V=coo_matrix((coo.data[graph_change],(coo.row[graph_change],coo.col[graph_change])),shape=H.shape).tocsr()
    Hblock=H-V
    HP=H[allowed][:,allowed]
    source_degree=sum(0 in e for e in host_edges)
    select=[]
    for row,col in zip(coo.row,coo.col):
        diff=int(bits[row]^bits[col]);edge_diff=diff >> n
        # Drive the ENTIRE frozen hopping term at vertex0. The apparatus
        # does not query F or select occupied neighbors; n_e in H activates
        # the existing bonds. Only the scalar source strength is calibrated.
        select.append(not graph_change[len(select)] and bool(diff & 1) and
                      edge_diff.bit_count()==1)
    select=np.array(select,dtype=bool)
    B=coo_matrix((coo.data[select]/math.sqrt(source_degree),
                  (coo.row[select],coo.col[select])),shape=H.shape).tocsr()
    assert np.linalg.norm((B-B.T).data)<1e-12
    normB=float(max(abs(eigsh(B,k=2,which="BE",return_eigenvectors=False,tol=1e-12))))
    BP=B[allowed][:,allowed]
    prepared_norm=float(max(abs(eigsh(BP,k=2,which="BE",return_eigenvectors=False,tol=1e-12))))
    assert prepared_norm <= 1+1e-10
    assert not np.any(masks[B.tocoo().row]!=masks[B.tocoo().col])
    sea=sum(1 << (n+e) for e in range(M))
    initial=[i for i,(mask,b,E) in enumerate(states) if mask & ~pool == 0 and b==sea and not any(E)]
    psi=np.zeros(len(states));psi[initial]=1/math.sqrt(len(initial))
    rhs=np.column_stack([psi,B @ psi])
    b2=.01/M**2*q*outside*(r+1)
    inside_norm=.1/M*q*r
    incidence=nx.Graph();incidence.add_nodes_from(range(n))
    source_support={0}
    for i,j in host_edges:
        e=pairs.index(tuple(sorted((i,j))));center=n+e
        incidence.add_edges_from([(i,center),(j,center)])
        if i==0 or j==0:source_support.add(center)
    rows=[];geometry_rows=[]
    for T in [.01,.05,.1,.4,1.,2.]:
        full=expm_multiply(-1j*T*H,rhs)
        projected=expm_multiply(-1j*T*HP,rhs[allowed])
        fixed=expm_multiply(-1j*T*Hblock,rhs)
        for edge,pair in enumerate(pairs):
            A=(masks >> edge & 1).astype(float)
            RF=-2*float(np.imag(np.vdot(full[:,0],A*full[:,1])))
            RB=-2*float(np.imag(np.vdot(fixed[:,0],A*fixed[:,1])))
            assert abs(RB)<1e-11
            assert abs(RF)<=4*b2*T*T+4*inside_norm*T+1e-10
            geometry_rows.append(dict(T=T,edge=list(pair),full_connectivity_response=RF,
                                      exact_block_connectivity_response=RB,
                                      proven_response_bound=4*b2*T*T+4*inside_norm*T))
        for detector in range(1,n):
            A=(bits >> detector & 1).astype(float)
            RF=-2*float(np.imag(np.vdot(full[:,0],A*full[:,1])))
            RP=-2*float(np.imag(np.vdot(projected[:,0],A[allowed]*projected[:,1])))
            RB=-2*float(np.imag(np.vdot(fixed[:,0],A*fixed[:,1])))
            assert abs(RF-RP) <= 4*b2*T*T+1e-10
            assert abs(RP-RB) <= 4*inside_norm*T+1e-10
            assert abs(RF-RB) <= 4*b2*T*T+4*inside_norm*T+1e-10
            item=dict(T=T,detector=detector,
                      independently_declared_incidence_distance=min(nx.shortest_path_length(incidence,x,detector) for x in source_support),
                      full_response=RF,pool_compressed_response=RP,graph_block_response=RB,
                      full_pool_error=abs(RF-RP),pool_error_bound=4*b2*T*T,
                      full_fixed_error=abs(RF-RB),total_error_bound=4*b2*T*T+4*inside_norm*T)
            if T==.4:
                eps=1e-4
                positive=expm_multiply(-1j*T*H,expm_multiply(1j*eps*B,psi))
                negative=expm_multiply(-1j*T*H,expm_multiply(-1j*eps*B,psi))
                impulse=float((A @ abs(positive)**2-A @ abs(negative)**2)/(2*eps))
                assert abs(impulse-RF)<1e-7
                item.update(independent_impulse_derivative=impulse,impulse_error=abs(impulse-RF))
            rows.append(item)
    out=dict(name=name,N=n,Q=q,S=S,L=L,dimension=len(states),pool_dimension=len(allowed),
             host_edges=host_edges,initial_coherent_graphs=len(initial),source_global_operator_norm=normB,
             source_prepared_sector_norm=prepared_norm,
             source_rule="All frozen hopping terms at vertex0, regardless of edge or pool membership; normalization alone is calibrated.",
             crossing_norm_square=b2,internal_offdiagonal_norm=inside_norm,rows=rows,
             connectivity_response_controls=geometry_rows,
             scope="Prepared and independently calibrated host. Gauge-invariant local source/detectors; no emergent dimension, physical clock, or Lorentz cone is claimed.")
    (P/f"prepared_response_{name}.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
    return out


if __name__ == "__main__":
    run([(0,1),(1,2),(2,3),(3,0)],"cycle_pool")
    run([(0,1),(0,2),(0,3),(1,2)],"branched_pool")
