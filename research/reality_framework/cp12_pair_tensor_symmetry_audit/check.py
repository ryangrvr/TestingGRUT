#!/usr/bin/env python3
"""Exact finite-N symmetry audit of CP-12 pair-multiplet interpretation.

This script studies the S_N action on unordered event pairs and a
SEPARATELY DECLARED uniform fixed-Q annealed graph-Laplacian comparator.
It does NOT re-diagonalize CP-12 nor assume its spectra are mean graph
Laplacian spectra. numpy required; no GRUT imports. No new law.
"""
import itertools, json, math
from pathlib import Path
import numpy as np

def analyze(n, rng):
    pairs=list(itertools.combinations(range(n),2))
    m=len(pairs); lookup={e:i for i,e in enumerate(pairs)}
    B=np.zeros((n,m))
    for ix,(i,j) in enumerate(pairs):
        B[i,ix]=B[j,ix]=1.
    share=B.T@B-2*np.eye(m)
    disjoint=np.ones((m,m))-np.eye(m)-share
    Ppair=np.eye(m)-B.T@np.linalg.inv(B@B.T)@B
    P1=np.ones((m,m))/m
    Pstd=np.eye(m)-Ppair-P1
    assert np.linalg.norm(Ppair@Ppair-Ppair)<1e-10
    assert np.linalg.norm(B@Ppair)<1e-10
    assert np.linalg.norm(Ppair@Pstd)<1e-10
    assert abs(np.trace(Ppair)-n*(n-3)/2)<1e-9
    for P in (P1,Pstd,Ppair):
        assert np.linalg.norm(P@P-P)<1e-10
    a,b,c=.73,-.19,.11
    K=a*np.eye(m)+b*share+c*disjoint
    eigenvalues=[
        a+2*(n-2)*b+(n-2)*(n-3)*c/2,
        a+(n-4)*b-(n-3)*c,
        a-2*b+c]
    for P,lam in zip((P1,Pstd,Ppair),eigenvalues):
        assert np.linalg.norm(K@P-lam*P)<1e-9
    assert abs(np.trace(K)-(eigenvalues[0]+(n-1)*eigenvalues[1]+n*(n-3)/2*eigenvalues[2]))<1e-9
    # CP-12 pair character equals pair-fixed count minus the N-dimensional
    # vertex permutation character. Test with sampled event permutations.
    if n<=7:
        for _ in range(20):
            perm=rng.permutation(n)
            T=np.zeros((m,m))
            for ix,(i,j) in enumerate(pairs):
                e=tuple(sorted((int(perm[i]),int(perm[j]))))
                T[lookup[e],ix]=1.
            fixed=sum(perm==np.arange(n))
            twocycles=sum(perm[int(perm[i])]==i and perm[i]!=i for i in range(n))//2
            character=math.comb(int(fixed),2)+twocycles-int(fixed)
            assert abs(np.trace(Ppair@T)-character)<1e-9
            assert np.linalg.norm(T@K-K@T)<1e-10
    # Comparator: mean adjacency of UNIFORM Q=N-1 simple graphs:
    # E[A_ij]=Q/binom(N,2). This is not the full quantum CP-12 state.
    q=n-1
    gap_uniform=n*q/m
    gap_ring=4*math.sin(math.pi/n)**2
    assert abs(gap_uniform-2.)<1e-12
    return dict(N=n,potential_edges=m,
       singlet_dimension=1,standard_dimension=n-1,
       pair_irrep_dimension=n*(n-3)//2,
       S_N_invariant_pair_kernel_eigenvalue=eigenvalues[-1],
       annealed_uniform_Q_N_minus_1_laplacian_gap=gap_uniform,
       prepared_cycle_laplacian_gap=gap_ring,
       ratio=gap_uniform/gap_ring)

def tt_3d_rank():
    # Conventional TT projection for a supplied physical wavevector along z.
    # Contrast: this physical k and spatial metric are NOT present in S_N pair labels.
    P=np.diag([1.,1.,0.])
    TT=np.zeros((9,9))
    for i,j,k,l in itertools.product(range(3),repeat=4):
        TT[3*i+j,3*k+l]=(
            .5*(P[i,k]*P[j,l]+P[i,l]*P[j,k])
            -.5*P[i,j]*P[k,l])
    assert np.linalg.norm(TT@TT-TT)<1e-12
    assert abs(np.trace(TT)-2)<1e-12
    return float(np.trace(TT))

def main():
    rng=np.random.default_rng(42012)
    samples=[analyze(n,rng) for n in (4,5,6,7,10,16)]
    assert [x["pair_irrep_dimension"] for x in samples[:4]]==[2,5,9,14]
    answer=dict(
        pinned_frozen_cp12_gravity_report="11fbaec7b28a1bbe88e6ec54808829ebbc804c48",
        theorem_one="Any S_N-invariant vertex-indexed two-point kernel has only two orbit coefficients; any S_N-invariant unordered pair-indexed kernel has only three: same pair, one shared endpoint, disjoint pairs.",
        theorem_two="On a single copy of [N-2,2], any S_N-equivariant pair operator acts as a scalar, a-2b+c; its irrep multiplicity N(N-3)/2 cannot be interpreted as helicity or wavevector.",
        theorem_three="Uniform fixed-Q annealed expected graph Laplacian = [2Q/(N-1)]*(I-J/N), has one zero mode and N-1 identical positive eigenvalues, not a spatial dispersive band.",
        qualification="Schur scalar applies to the single pair-label representation, not an entire CP12 many-body isotypic sector with potentially many copies; relational, graph-conditioned probes may evade label-only orbit constraints.",
        supplied_3d_TT_projector_rank=tt_3d_rank(),checks=samples,
        not_established=["No CP12 two-helicity tensor pole","No k-dependent physical spin 2","No Lorentzian light cone or universal matter coupling","No full thermal/ground phase inferred from annealed graph average"])
    Path(__file__).with_name("results.json").write_text(json.dumps(answer,indent=2)+"\n")
    print(json.dumps(samples,indent=2))
    print("CP12_PAIR_TENSOR_SYMMETRY_INDEPENDENT_PASS")

if __name__=="__main__":
    main()
