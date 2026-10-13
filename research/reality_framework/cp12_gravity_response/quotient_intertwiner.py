"""Independent full-builder check of the signed symmetry compression."""
from pathlib import Path
import json
import sys
import math
import numpy as np
from scipy.sparse import coo_matrix, eye

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build
from marked_pair_quotient import Quotient


def run(n,q,L):
    quotient=Quotient(n,q,L);HQ,C,_=quotient.matrices()
    full,states,_,_=build(n,q,L)
    assert len(states)==quotient.full_dimension
    rows=[];cols=[];data=[]
    for row,state in enumerate(states):
        key,sign,size=quotient.canonical(state)
        if sign:rows.append(row);cols.append(quotient.index[key]);data.append(sign/math.sqrt(size))
    T=coo_matrix((data,(rows,cols)),shape=(len(states),len(quotient.states))).tocsr()
    orthogonal=T.T @ T-eye(T.shape[1],format='csr')
    residual=full @ T-T @ HQ
    assert np.linalg.norm(orthogonal.data)<1e-10 and np.linalg.norm(residual.data)<1e-9
    out=dict(N=n,Q=q,L=L,full_dimension=full.shape[0],quotient_dimension=HQ.shape[0],
             isometry_residual=float(np.linalg.norm(orthogonal.data)),
             full_H_intertwining_residual=float(np.linalg.norm(residual.data)),
             exact_integer_checks=quotient.integer_controls)
    print(json.dumps(out),flush=True)
    return out


if __name__=="__main__":
    out=dict(status="CP12_SIGNED_QUOTIENT_INTERTWINER_PASS",controls=[run(4,3,3),run(5,4,2)],
             scope="Complete frozen builder compared to an independently generated signed-orbit Hamiltonian. No fixed-graph or hopping expansion.")
    (P/"quotient_intertwiner_results.json").write_text(json.dumps(out,indent=2)+"\n")
