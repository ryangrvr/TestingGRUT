# R1: separating static source calibration from response shape

This diagnostic keeps the D1 source dictionary, Hamiltonian and boundary convention fixed. It can distinguish specified static-calibration explanations from additional response-shape differences. It cannot uniquely factor an arbitrary observable map from an arbitrary Hamiltonian: source vertices, spectral measures, ground preparation and boundaries jointly determine correlation functions. All tests and interpretation rules must be fixed before new responses are inspected.

## Static matrices and exact source sum rule

For registered source probes G_i, define A_i=G_i-<G_i> and |v_i>=A_i|Omega>, with H|Omega>=E0|Omega>. Let V have the |v_i> as columns. Then

    C = V^dagger V,
    K = V^dagger (H-E0) V,
    C(t) = V^dagger exp[-i(H-E0)t] V.

C is the connected equal-time covariance. K is the first energy-weighted moment, defined as an energy quadratic form. For the real ground state and real commuting probes in this study,

    S_ij = <[G_i,[H,G_j]]> = 2 K_ij.

For G_i=sum_l v_i,l y_l+c_i and H_hop,l=(T_l+T_l^dagger)/(2a), the exact lattice algebra [G_i,T_l]=v_i,l T_l gives

    S_ij(psi) = -sum_l v_i,l v_j,l <H_hop,l>_psi.

The source coefficients v_i,l should not be confused with the state vectors |v_i>. The lattice ground-state identity K=S/2 can be checked in two independent ways: apply H-E0 to the source-state vectors, and sum the bond hopping expectations.

For the matched continuum step sources f_i,a(z)=v_i,l/a on each interior cell and zero outside, with the declared constant calibration kept separately,

    C_ij(t) = sum_n f_i,n f_j,n exp(-i omega_n t)/(2 pi omega_n),
    K_ij = (1/(2 pi)) integral f_i,a(z) f_j,a(z) dz
         = (1/(2 pi a)) sum_l v_i,l v_j,l,
    S_ij = 2K_ij.

Use the exact Parseval expression for K. Truncating its mode sum would leave an O(1/M) tail, whereas C(t) has the covariance tail already controlled by the step-source variation bound. Source offsets cancel from connected covariances; they remain essential physical phases in the historical detector instrument.

The infinitesimal response to a unitary kick exp(i kappa G_j) is

    chi_ij(t) = i<[G_i(t),G_j]> = -2 Im C_ij(t),
    chi'_ij(0) = S_ij.

This fixes the response sign and the factor of two without a numerical derivative.

## Exact static-amplitude versus normalized-shape decomposition

For each scalar probe, set V_X=C_X,ii(0) and R_X(t)=C_X,ii(t)/V_X, where X is lattice or continuum. An affine scalar calibration G_i -> alpha_i G_i + beta_i leaves R_i unchanged. Therefore a significant normalized-response discrepancy cannot be repaired by that probe's gain or offset.

There are two exact decompositions of the complex correlation discrepancy:

    Delta C(t) = (V_L-V_C) R_C(t) + V_L [R_L(t)-R_C(t)],
    Delta C(t) = (V_L-V_C) R_L(t) + V_C [R_L(t)-R_C(t)].

Preserve both signed decompositions. Their terms can cancel, and they are not unique causal percentages. They separate amplitude factors from normalized spectral shape algebraically, not every possible source-operator correction from dynamics.

## Invariants under static mixing of the registered source span

Use the three sine probes as the primary span and the asymmetric beta probe as a separate witness. Under an invertible constant source mixing V' = V R,

    C' = R^dagger C R,
    K' = R^dagger K R,
    C'(t) = R^dagger C(t) R.

Since det(K'-nu C')=|det R|^2 det(K-nu C), the generalized eigenvalues of K x=nu C x are unchanged. They are energy-weighted centroids/Rayleigh-Ritz values of the fixed source span, not automatically physical excitation eigenvalues.

Likewise the generalized eigenvalues of

    B(t) x = lambda C x,
    B(t) = [C(t)+C(t)^dagger]/2,

are invariant. B(t) is the compressed cosine response. For the real probes used here, it equals the real symmetric part of C(t). Compare its sorted real eigenvalue sets, rather than attaching mode labels that can exchange order as cosine phases evolve.

A mismatch in either invariant excludes a constant invertible calibration confined to this span as a complete explanation. It does not exclude outside-span derivative/source corrections, a different endpoint prescription, or changed state matching; therefore it is not a Hamiltonian-only attribution. Agreement is compatibility, not proof that a single static map explains all response matrices jointly.

Whitening requires a positive, adequately conditioned covariance. If the frozen conditioning criterion fails, classify the invariant diagnostic as unresolved rather than regularizing or dropping a source after seeing the results.

## Prospectively fixed Galerkin response predictor

Let U=V C^(-1/2), so U^dagger U=I, and define

    A=C^(-1/2) K C^(-1/2).

The Hamiltonian-derived compressed predictor is

    C_Gal(t)=C^(1/2) exp(-i A t) C^(1/2),
    M_Gal(t)=exp(-i A t),
    M_true(t)=C^(-1/2) C(t) C^(-1/2).

No finite-time trajectory is fitted. C and K are computed from the registered Hamiltonian, sources and ground state; new times are held out from construction. Test lattice and continuum versions against their own full responses. A registered pass supports this finite predictive linear-response compression on those sources and times. It does not resolve the nonlinear detector-history task or supply a new GRUT law.

Under a static source mixing, the whitened matrices transform by unitary similarity: with Q=C^(1/2) R (C')^(-1/2), Q^dagger Q=I and A'=Q^dagger A Q. This explains their calibration invariance.

Step smearings have finite energy form, but the source-state vectors need not lie in the operator domain of H. Higher spectral moments and strong residual norms such as ||(I-UU^dagger)HU|| can diverge in the continuum. Do not turn the finite K calculation into an H^2 or operator-norm error certificate. The finite-time comparison is the registered validation.

## Reached-state and coincident-kick source algebra

The continuum S_ij=(1/pi)integral f_i,a f_j,a is state independent as a suitable quadratic-form identity. On the lattice, the bond expectation in S_ij(psi) can depend on the reached state. The diagonal ratio

    S_ii(psi)/S_ii(Omega)-1

is invariant under scalar gain. If full matrices are registered, the norm of

    S_vac^(-1/2) [S_psi-S_vac] S_vac^(-1/2)

also gives an invariant test of a common invertible mixing within the primary span. A finite-a state dependence cannot be removed by one state-independent normalization across those states. This locates a canonical source/current-algebra difference; it is not by itself a pure Hamiltonian attribution.

There is an especially clean no-evolution control. Let F be the unchanged historical detector and W=exp(i kappa F). For a real, zero-current ground state,

    <H_hop,l>_(W^m Omega) = cos(m kappa w_F,l) <H_hop,l>_Omega,
    S_ij(W^m Omega) = -sum_l v_i,l v_j,l cos(m kappa w_F,l) <H_hop,l>_Omega.

For m=1 and 2 these describe one and two coincident unitary kicks. They require no free propagation. A discrepancy there already exists without accumulated lattice evolution. The bond-cosine relation can independently check the direct state expectation. Delayed states test additional propagation/backaction effects, while retaining the same sources and calibration.

The prospective unconditional-kick scan uses kappa in {0.5,1,2,4}. On the primary three-source span, report two distinct full-matrix quantities:

    disturbance(kappa) = ||S_vac^(-1/2) [S_kappa-S_vac] S_vac^(-1/2)||_2,
    absolute(kappa) = ||S_C^(-1/2) [S_kappa-S_C] S_C^(-1/2)||_2.

The registered 0.01 thresholds define an unconditional source-current window only at the tested amplitudes. A small disturbance does not establish absolute continuum agreement when the vacuum baseline already differs. Neither criterion bounds selective detector histories, all intermediate amplitudes, or arbitrary states.

For the unchanged selective instrument K_s=(I+sW)/2, s=+1 or -1, two coincident recorded outcomes prepare

    K_t K_s|Omega> = [I+(s+t)W+st W^2]|Omega>/4.

Equivalently, use the ordered path basis B=(|Omega>,W|Omega>,W|Omega>,W^2|Omega>) and x_st=(1,s,t,st)/4. The norm x_st^dagger B^dagger B x_st is the branch probability; insert G, G^2 or the double-commutator operator between the same vectors and divide by that probability. The duplicate W columns make this four-path Gram matrix singular, which is harmless because the calculation requires no Gram inverse. Every registered branch is retained subject to the prospective numerical probability floor.

Compare the coincident conditional means and second moments to the unchanged D1 Gaussian reference evaluated with both backaction times and the probe time equal to zero. This uses D1's moment tolerance and reserve; it does not rewrite D1's delayed result. A zero-delay failure cannot be assigned to accumulated free propagation. A zero-delay pass followed by delayed failure identifies an additional response-shape/backaction issue, but still does not uniquely distinguish the Hamiltonian from frequency-dependent or outside-span source corrections. Coincident and delayed reached-state sum rules are reported separately from these conditional moment tests.

## Interpretation of prospective numerical classes

The proposed 1% static C/K and reached-state sum-rule tolerances, 0.5% centroid and normalized-response tolerances, 0.005 absolute cosine-spectrum tolerance, and 0.01 whitened Galerkin error are operational diagnostic choices. Each includes the prospective additive 1e-4 numerical guard before comparison with its allowance, on both standard and refined N22. The coincident moment/insertion allowance is 0.002 with its additive 1e-4 reserve. Normalized scalar responses retain all four probes; matrix whitening and invariant spectra use the three-sine span. The reached-state and kick-scan matrix errors use inverse-square-root metric whitening, not a ratio of unwhitened operator norms. These are not derived universal physical error bounds. Keep their conclusions separate rather than requiring every class to pass before any inference is allowed:

- Static covariance mismatch with normalized response and invariants passing is compatible with leading amplitude/static mismatch.
- Normalized scalar-response mismatch makes scalar gain/offset correction insufficient.
- Invariant mismatch makes constant mixing within the primary span insufficient.
- Coincident-kick sum-rule variation identifies state-dependent source/current algebra without free evolution.
- A Galerkin prediction pass establishes only the registered finite response-compression task.
- Ill-conditioned or numerically uncontrolled comparisons remain unresolved.

None of these classifications retrospectively changes B2, N22 or D1, authorizes a massive transfer, or proves a distinctive empirical GRUT prediction.
