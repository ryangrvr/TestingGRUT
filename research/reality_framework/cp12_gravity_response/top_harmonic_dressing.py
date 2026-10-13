"""Exact graph-harmonic locking/unlocking and complete small-law controls.

This is a projection analysis of unchanged CP-12. No new constraint or
Hamiltonian is imposed. The projected space is not claimed invariant.
"""
from pathlib import Path
from collections import defaultdict
import itertools
import json
import math
import sys

import numpy as np
from scipy.linalg import null_space, svdvals
from scipy.sparse.linalg import eigsh

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import build


def top_projector(K,r):
    subsets=list(itertools.combinations(range(K),r)); index={s:i for i,s in enumerate(subsets)}
    if r==0:return np.ones((1,1)),subsets,index
    lower=list(itertools.combinations(range(K),r-1)); li={s:i for i,s in enumerate(lower)}
    down=np.zeros((len(lower),len(subsets)))
    for j,s in enumerate(subsets):
        for a in s:down[li[tuple(x for x in s if x!=a)],j]=1
    frame=null_space(down)
    assert frame.shape[1]==math.comb(K,r)-math.comb(K,r-1)
    return frame @ frame.T,subsets,index


def nested_controls():
    rows=[]
    for K,r in [(4,1),(5,1),(5,2),(7,2),(9,2)]:
        old,subsets,index=top_projector(K,r)
        new,newsets,newindex=top_projector(K+1,r+1)
        lift=np.zeros((len(newsets),len(subsets)))
        for j,s in enumerate(subsets):lift[newindex[s+(K,)],j]=1
        leakage=(np.eye(len(newsets))-new) @ lift @ old
        denominator=K+1-2*r
        # All singular values on the old top space agree, including r=0
        # special case outside these controls.
        operator=old @ lift.T @ (np.eye(len(newsets))-new) @ lift @ old
        error=float(np.linalg.norm(operator-old/denominator))
        assert error<1e-11
        locking_residual=float(np.linalg.norm((np.eye(len(subsets))-old) @ lift.T @ new))
        assert locking_residual<1e-11
        rows.append(dict(old_K=K,old_r=r,new_K=K+1,new_r=r+1,
                         exact_unlock_leakage_norm_square=f"1/{denominator}",
                         numerical_unlock_leakage_norm_square=float(svdvals(leakage)[0]**2),
                         projector_identity_residual=error,
                         top_preserving_locking_residual=locking_residual,
                         reverse_top_transport_norm_square=float(svdvals(old @ lift.T @ new)[0]**2),
                         expected_top_transport_norm_square=1-1/denominator))
    return rows


def full_control():
    n=4;q=3;L=3;M=n*(n-1)//2
    H,states,_,_=build(n,q,L); masks=np.array([s[0] for s in states]);coo=H.tocoo()
    V=np.zeros(H.shape)
    for i,j,v in zip(coo.row,coo.col,coo.data):
        if masks[i]!=masks[j]:V[i,j]=v
    projection=np.zeros(H.shape);r_values=np.zeros(len(states));fibers=defaultdict(list)
    for i,(_,bits,E) in enumerate(states):fibers[(bits,E)].append(i)
    for (bits,E),ids in fibers.items():
        locked=[e for e in range(M) if not(bits >> (n+e) & 1) or E[2*e] or E[2*e+1]]
        eligible=[e for e in range(M) if e not in locked];K=len(eligible);r=q-len(locked)
        block,subsets,index=top_projector(K,r)
        order=[index[tuple(eligible.index(e) for e in eligible if states[i][0] >> e & 1)] for i in ids]
        projection[np.ix_(ids,ids)]=block[np.ix_(order,order)];r_values[ids]=r
    assert np.linalg.norm(projection @ projection-projection)<1e-11
    Videntity=projection @ V @ projection-.1/M*projection*r_values[None,:]
    assert np.linalg.norm(Videntity)<1e-11
    hopping_leakage=(np.eye(len(states))-projection) @ (H @ projection)
    actual=float(svdvals(hopping_leakage)[0]);bound=math.sqrt(2*M*q/(M-2*q+1))
    assert actual<=bound+1e-10
    energy,ground=eigsh(H,k=1,which='SA',tol=1e-12);g=ground[:,0]
    sea=sum(1 << (n+e) for e in range(M));neutral=[i for i,s in enumerate(states) if s[1]==sea and not any(s[2])][0]
    exact_neutral_top_weight=1-q/(M-q+1)
    assert abs(projection[neutral,neutral]-exact_neutral_top_weight)<1e-12
    dressed=projection[:,neutral]/math.sqrt(projection[neutral,neutral])
    before=float(dressed @ (H @ dressed))
    after=0.
    for mask in np.unique(masks):
        component=dressed*(masks==mask)
        after+=float(component @ (H @ component))
    predicted_record_shift=-.1/M*float(r_values @ dressed**2)
    assert abs(after-before-predicted_record_shift)<1e-11
    # Does this particular nonlocal projection manufacture a geometry/source
    # commutator? Yes at finite N, but its norm has a vanishing sparse-N bound.
    bits=np.array([s[1] for s in states]); source=np.zeros(H.shape)
    for i,j,v in zip(coo.row,coo.col,coo.data):
        diff=int(bits[i]^bits[j])
        if masks[i]==masks[j] and diff & 1 and (diff >> n).bit_count()==1:source[i,j]=v/math.sqrt(n-1)
    projected_source=projection @ source @ projection
    source_leak=float(svdvals((np.eye(len(states))-projection) @ source @ projection)[0])
    commutators=[]
    pairs=list(itertools.combinations(range(n),2))
    for e in [pairs.index((0,1)),pairs.index((1,2))]:
        A=(masks >> e & 1).astype(float)
        assert np.linalg.norm(A[:,None]*source-source*A[None,:])<1e-12
        projected_A=(projection*A[None,:]) @ projection
        commutator=projected_A @ projected_source-projected_source @ projected_A
        norm=float(svdvals(commutator)[0])
        readout_leak=float(svdvals((np.eye(len(states))-projection) @ (A[:,None]*projection))[0])
        assert norm<=2*readout_leak*source_leak+1e-11
        assert readout_leak**2<=1/(M-2*q+2)+1e-11
        commutators.append(dict(edge=list(pairs[e]),bare_commutator_norm=0.,
             projected_commutator_norm=norm,exact_complement_product_bound=2*readout_leak*source_leak,
             readout_leakage_norm_square=readout_leak**2))
    return dict(N=n,Q=q,L=L,full_dimension=len(states),physical_fibers=len(fibers),
           top_projection_rank=int(round(np.trace(projection))),
           exact_V_top_identity_residual=float(np.linalg.norm(Videntity)),
           actual_full_H_crossing_norm=actual,proven_crossing_norm_upper_bound=bound,
           ground_energy=float(energy[0]),ground_top_weight=float(g @ projection @ g),
           definite_neutral_graph_top_weight=float(projection[neutral,neutral]),
           exact_definite_neutral_graph_top_weight=exact_neutral_top_weight,
           projected_preparation_record_control=dict(before=before,after_graph_record=after,
               system_energy_change=after-before,exact_change=predicted_record_shift),
           projected_source_leakage_norm=source_leak,projected_geometry_commutator_controls=commutators,
           global_projected_commutator_bound=2*math.sqrt(n-1)/math.sqrt((M-2*q+2)*(M-2*q+1)))


if __name__=="__main__":
    out=dict(status="CP12_TOP_HARMONIC_DRESSING_PASS",nested_fiber_controls=nested_controls(),
             full_control=full_control(),
             scope="Exact frozen-law projection and transition identities. Top harmonics are not an invariant physical sector, a new law, or gravitational dressing.")
    (P/"top_harmonic_dressing_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out['full_control'],indent=2))
