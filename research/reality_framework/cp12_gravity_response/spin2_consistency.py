"""Independent conventional spin-2/energy-constraint target regressions.

These test the GR benchmark, NOT CP-12 predictions. They explicitly state
which structure a successful microscopic limit would need to generate.
"""
from pathlib import Path
import json
import math

import numpy as np
from scipy.linalg import null_space

P=Path(__file__).resolve().parent


def tensor_basis():
    out=[]
    for i in range(3):
        A=np.zeros((3,3));A[i,i]=1;out.append(A)
    for i,j in [(1,2),(0,2),(0,1)]:
        A=np.zeros((3,3));A[i,j]=A[j,i]=1/math.sqrt(2);out.append(A)
    return np.array(out)


def constraints(k):
    basis=tensor_basis();k=np.asarray(k,float);k2=k @ k
    # canonical phase vector (h_1,...,h_6,pi_1,...,pi_6)
    C=np.zeros((4,12))
    C[0,:6]=np.einsum('i,aij,j->a',k,basis,k)-k2*np.trace(basis,axis1=1,axis2=2)
    C[1:,6:]=np.einsum('aij,j->ia',basis,k)
    J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])
    assert np.linalg.matrix_rank(C)==4 and np.linalg.norm(C @ J @ C.T)<1e-12
    gauge=J @ C.T
    assert np.linalg.matrix_rank(gauge)==4
    Pk=np.eye(3)-np.outer(k,k)/k2
    TT=np.empty((6,6))
    for a,A in enumerate(basis):
        projected=Pk @ A @ Pk-.5*Pk*np.trace(Pk @ A)
        TT[:,a]=np.einsum('aij,ij->a',basis,projected)
    assert np.linalg.norm(TT @ TT-TT)<1e-12 and np.linalg.matrix_rank(TT,tol=1e-10)==2
    # Transverse traceless canonical h/pi automatically satisfy all constraints.
    tt_frame=np.linalg.eigh(TT)[1][:,-2:]
    physical=np.block([[tt_frame,np.zeros((6,2))],[np.zeros((6,2)),tt_frame]])
    assert np.linalg.norm(C @ physical)<1e-12
    assert np.linalg.norm(physical.T @ gauge)<1e-12
    # In TT variables H=pi^2+(k^2/4)h^2: two healthy omega^2=k^2 modes.
    kinetic=np.block([[np.eye(2)*k2/2,np.zeros((2,2))],[np.zeros((2,2)),2*np.eye(2)]])
    Jp=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]])
    generator=Jp @ kinetic
    assert np.linalg.norm(generator @ generator+k2*np.eye(4))<1e-12
    return dict(momentum=k.tolist(),constraint_rank=4,first_class_bracket_norm=float(np.linalg.norm(C @ J @ C.T)),
                unconstrained_phase_dimension=12,physical_phase_dimension=12-2*4,tt_projector_rank=2,
                expected_frequency=float(math.sqrt(k2)),dynamical_eigenvalues_imaginary=np.sort(np.linalg.eigvals(generator).imag).tolist())


def minkowski_dot(a,b):return a[0]*b[0]-a[1:] @ b[1:]


def soft_ward(theta):
    m1=1.;m2=2.;p=.7
    p1=np.array([math.sqrt(m1*m1+p*p),0,0,p]);p2=np.array([math.sqrt(m2*m2+p*p),0,0,-p])
    k=np.array([p*math.sin(theta),0,p*math.cos(theta)])
    p3=np.r_[p1[0],k];p4=np.r_[p2[0],-k]
    particles=np.array([p1,p2,p3,p4]);eta=np.array([-1,-1,1,1])
    assert np.linalg.norm(eta @ particles)<1e-14
    q=np.array([.01,.006,0,.008]);xi=np.array([.3,-.7,.2,.5])
    longitudinal=np.outer(q,xi)+np.outer(xi,q)
    # Lowering epsilon indices implements p^mu p^nu epsilon_mu,nu.
    metric=np.diag([1,-1,-1,-1]);lower=metric @ longitudinal @ metric
    rows=[]
    for g1,g2 in [(1.,1.),(1.,1.1),(1.,2.)]:
        couplings=np.array([g1,g2,g1,g2])
        gauge_change=sum(sign*g*(pvec @ lower @ pvec)/minkowski_dot(pvec,q)
                         for sign,g,pvec in zip(eta,couplings,particles))
        independently=2*minkowski_dot(np.sum((eta*couplings)[:,None]*particles,axis=0),xi)
        assert abs(gauge_change-independently)<1e-12
        if g1==g2:assert abs(gauge_change)<1e-12
        else:assert abs(gauge_change)>1e-4
        rows.append(dict(species_couplings=[g1,g2],gauge_polarization_change=float(gauge_change),
                         conserved_momentum_ward_expression=float(independently)))
    return dict(scattering_angle=theta,controls=rows)


if __name__=="__main__":
    out=dict(status="CONVENTIONAL_SPIN2_TARGET_CHECKS_PASS",
             canonical_constraint_controls=[constraints(k) for k in [[0,0,1],[1,2,3],[.3,-.7,1.1]]],
             universal_soft_coupling_controls=[soft_ward(t) for t in [.3,.7,1.2]],
             scope="Supplied conventional linearized GR and Weinberg soft-pole structure. Not emergent constraints, Lorentz symmetry, or universal coupling derived from CP-12.")
    (P/"spin2_consistency_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(status=out['status'],physical_modes=2),indent=2))
