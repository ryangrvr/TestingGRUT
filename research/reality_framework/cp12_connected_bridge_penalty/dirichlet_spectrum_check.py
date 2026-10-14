#!/usr/bin/env python3
"""Connected-sector quantum graph Laplacian Dirichlet check (neutral fiber).

Checks all connected Q=N graphs at N=5 and N=6 with SciPy/networkx.
Other physical matter/electric fibers follow the analytic edge-deletion proof.
"""
from itertools import combinations
import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh

for n in (5,6):
    edge_pairs=list(combinations(range(n),2));M=len(edge_pairs);Q=n
    masks=[];bridgecounts=[];exitcounts=[]
    for occupied in combinations(range(M),Q):
        G=nx.Graph();G.add_nodes_from(range(n));G.add_edges_from(edge_pairs[e] for e in occupied)
        if not nx.is_connected(G):continue
        masks.append(sum(1<<e for e in occupied))
        bridges=list(nx.bridges(G));bridgecounts.append(len(bridges));ex=0
        for e in bridges:
            H=G.copy();H.remove_edge(*e)
            s=len(next(iter(nx.connected_components(H))))
            ex+=M-Q+1-s*(n-s)
        exitcounts.append(ex)
    idx={mask:i for i,mask in enumerate(masks)}
    rr=[];cc=[];vv=[]
    for i,mask in enumerate(masks):
        used=[e for e in range(M) if mask>>e&1]
        absent=[e for e in range(M) if not mask>>e&1]
        for e in used:
            for f in absent:
                j=idx.get(mask^(1<<e)^(1<<f))
                if j is not None:rr.append(i);cc.append(j);vv.append(-.1/M)
    R=coo_matrix((vv,(rr,cc)),shape=(len(masks),len(masks))).tocsr()
    R=R+diags(np.full(len(masks),.1*Q*(M-Q)/M))
    boundary=diags(.1/M*np.asarray(exitcounts))
    w=.1*max(0,M-Q+1-(n*n//4))/M
    bridge=diags(w*np.asarray(bridgecounts))
    for label,A in [('interior_positive_laplacian',R-boundary),('bridge_penalty_inequality',R-bridge)]:
        eig=float(eigsh(A,k=1,which='SA',tol=1e-10,return_eigenvectors=False)[0])
        assert eig>-1e-8
        print('N=',n,'connected=',len(masks),label,'lambda_min=',eig)
print('CP12_CONNECTED_DIRICHLET_NUMERICAL_SPECTRUM_PASS')
