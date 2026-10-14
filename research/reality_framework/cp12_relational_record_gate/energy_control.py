"""A physical mean-energy-matched reference, not a classical flux shuffle.

Energy-dephase the INITIAL ray in the already used N4 retained Hamiltonian.
No additional dynamical evolution or physical change to H is introduced.
The resulting positive density operator is supported on physical Gauss states.
It preserves all *projected* energy probabilities and exactly the full-law
mean energy; higher full-law moments can involve the outgoing B matrix.
"""
from pathlib import Path
import json, time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from analyze import P, model, resources, setup


def run():
    start=time.monotonic();n=q=4;budget=256
    sub,H,B,states,sizes,psi=model(n,q,budget)
    r,k,nv,e2=resources(n,q,sub,states)
    S=setup(n,q,states,sizes,r,nv,e2)
    values,V=eigh(H.toarray())
    coefficients=V.T.conj()@psi
    # Retain coherence inside numerically degenerate energy blocks.
    blocks=[];begin=0
    for j in range(1,len(values)+1):
        if j==len(values) or values[j]-values[begin]>1e-10:
            blocks.append((begin,j));begin=j
    columns=[];eigen_spreads=[]
    for a,b in blocks:
        columns.append(V[:,a:b]@coefficients[a:b])
        eigen_spreads.append(float(values[b-1]-values[a]))
    X=np.column_stack(columns)
    p=np.sum(abs(X)**2,axis=1)
    weights=np.sum(abs(X)**2,axis=0)
    energy=float(np.sum(X.conj()*(H@X)).real)
    second=float(np.linalg.norm(H@X)**2+np.linalg.norm(B@X)**2)
    projected_second=float(np.linalg.norm(H@X)**2)
    variance=second-energy*energy
    norm=float(weights.sum())
    assert abs(norm-1)<1e-11 and abs(energy-q/2)<1e-10
    rho=X@X.conj().T
    comm=H@rho-(H@rho.conj().T).conj().T
    assert np.linalg.norm(comm)<1e-9
    hist=np.bincount(S['classes'],weights=p,minlength=len(S['catalog']['keys']))
    # The maximum-entropy reference is restricted to this already retained
    # symmetry sector. It is a genuine physical mean-energy comparator,
    # not a claim about the full unbounded model's thermal state.
    def gibbs(beta):
        w=np.exp(-beta*(values-values[0]));return w/w.sum()
    hi=1.
    while gibbs(hi)@values>q/2:hi*=2
    beta=brentq(lambda x:gibbs(x)@values-q/2,0.,hi,xtol=1e-13)
    gw=gibbs(beta);gp=abs(V)**2@gw
    BV=B@V
    genergy=float(gw@values)
    gsecond=float(gw@(values**2+np.sum(abs(BV)**2,axis=0)))
    maximum_entropy=dict(inverse_temperature=beta,full_mean_energy=genergy,
        full_energy_variance=gsecond-genergy**2,
        entropy_bits=float(-sum(w*np.log2(w) for w in gw if w>0)),
        means=dict(zip(S['names'],map(float,gp@S['actual']))),
        removable_fraction=float(gp@r/q),flux_protected_fraction=float(gp@k/q),
        resources=dict(vertex_fermions=float(gp@nv),electric_square=float(gp@e2)),
        scope='Maximum entropy at the declared mean energy within the existing finite retained signed-invariant sector. Positive physical state, coherent in the configuration basis. Not the full unbounded-H Gibbs state, not full-law stationary, and not resource-, variance- or purity-matched to the evolved pure state.')
    assert abs(genergy-q/2)<1e-11
    out=dict(N=n,Q=q,budget=budget,dimension=len(states),energy_blocks=len(blocks),
        maximum_degenerate_block_spread=max(eigen_spreads),
        initial_and_reference_full_mean_energy=q/2,
        numerical_reference_full_mean_energy=energy,
        reference_full_energy_variance=variance,initial_full_energy_variance=2*q,
        reference_projected_energy_second_moment=projected_second,
        full_second_moment_outgoing_correction=float(np.linalg.norm(B@X)**2),
        reference_mixed_state_entropy_bits=float(-sum(w*np.log2(w) for w in weights if w>0)),
        projected_stationarity_residual=float(np.linalg.norm(comm)),
        means=dict(zip(S['names'],map(float,p@S['actual']))),
        core_shape_distribution=hist.tolist(),
        removable_fraction=float(p@r/q),flux_protected_fraction=float(p@k/q),
        resources=dict(vertex_fermions=float(p@nv),electric_square=float(p@e2)),
        maximum_entropy_reference=maximum_entropy,
        scope='Positive physical comparator supported on the existing Gauss-law basis. Exact mean-energy equality follows because its projected energy probabilities equal those of the initial ray. Projected spectral dephasing is not asserted to be stationary under the unbounded full H, to preserve its complete energy distribution, or to be an autonomous state-selection/preparation mechanism. It is mixed and does not match each instantaneous resource distribution.',
        seconds=time.monotonic()-start)
    (P/'MEAN_ENERGY_CONTROL.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)
    return out


if __name__=='__main__':run()
