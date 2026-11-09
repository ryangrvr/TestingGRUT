# R1-01: bounded source, response and backaction separation

This document and executable protocol.json are frozen before R1 physical evaluation. Parent is D1 result 7ba04fd35e47870cbdb41268b2d7ffd33f567f51. D1's existing results informed the question, not an unseen-outcome claim about static covariances. N22, B2 and D1 retain their historical verdicts. Massive transfer stays withheld.

## What this experiment can separate

For the same archived ground states and four source-native probes, form centered source vectors Xi_i=(G_i-<G_i>)|0>. Measure C=Xi^dagger Xi, K=Xi^dagger(H-E0)Xi and M(t)=Xi^dagger exp[-i(H-E0)t]Xi. The four scalar normalized responses M_ii/C_ii remove each source's gain and offset. The generalized spectra (K,C) and (Hermitian M(t),C), using only the first three sine probes, also remove invertible mixing within that fixed source span. Beta_left remains a scalar witness. No source eigenmode may be clipped or dropped.

An arbitrary new detector can change its excited spectral subspace. Consequently a normalized-response discrepancy rules out a simple gain/offset explanation, and an invariant-spectrum discrepancy rules out an invertible within-span mixing explanation, but neither uniquely identifies an error in H. The synthetic suite explicitly demonstrates this limit with fixed H, variance and centroid but changed response.

Archive both exact signed decompositions for each auto response:

    M_l-M_c = (C_l-C_c) R_c + C_l(R_l-R_c)
            = (C_l-C_c) R_l + C_c(R_l-R_c).

The terms are complex, can cancel, and depend on ordering; do not turn them into causal percentages. These identities and source invariants provide a bounded diagnostic, not an inter-Hilbert-space identification theorem.

## Fixed probes, states and no-evolution controls

All geometry, normalization, offsets, boundaries, bare mass, detector coupling kappa=2 and instrument phase=0 are unchanged from D1. Use N=14,18,20,22 and the archived refined N22 ground state. The refined run halves each propagation step. INPUTS.json identifies all five checkpoint hashes. The source construction is operationally calibrated; outer-cell contributions remain identity calibration, not assumed constant field cell averages.

Use new response times 0,.29,.83,1.57. Evaluate the source-current f-sum S_gh=<[G_g,[H,G_h]]>=-sum_j v_gj v_hj <Hhop,j>. Independently check S=2K in the real ground state. The continuum matched-source f-sum is S_c=(v v^T)/(pi a), independent of the state. Its finite-lattice state dependence is directly observable without evolving a branch.

Evaluate S on vacuum, W|0>, both K branches, W^2|0>, all four coincident K_s K_t branches, and the four D1 branches with operations at 0,.43 observed at1.31. Coincident branches use exactly the same instrument with zero intervening free evolution. Their dependent unitary basis is intentionally rank deficient; probabilities/moments use quadratic forms without any Gram inverse.

Evaluate all first/second probe moments and signed insertion matrices for the zero-delay one/two-instrument controls. A failure there cannot be caused by accumulated free evolution. It can still reflect the initial lattice state or source dictionary as well as finite-spacing algebra.

The prospectively fixed unconditional kick strengths .5,1,2,4 are diagnostic candidates. A current window requires BOTH ||S0^-1/2(S_k-S0)S0^-1/2||+1e-4<=.01 AND ||Sc^-1/2(S_k-Sc)Sc^-1/2||+1e-4<=.01 on the primary3x3 block. This says nothing by itself about conditional branches, delayed histories or massive transfer.

## Static-only predictor and gates

The three-source Galerkin predictor is C^(1/2) exp[-it C^(-1/2)K C^(-1/2)] C^(1/2). Build it only from static quantities, independently on lattice and continuum, then compare at the new times. It is a proposed finite predictive reduction; the source span need not be invariant under H. Step sources have finite first spectral moment but no assumed finite second moment, so no H^2, energy-variance or rigorous model-space error estimate is used.

Numerical eligibility requires source/synthetic identities, normalization and residuals <=1e-9; same-size refinement<=1e-8; reference mode doubling<=1e-7; conservative reference tail estimates<=1e-4; primary Gram eigenvalue>1e-10 and condition<1e6; single branch probability>=.01 and double branch>=.001. Missing or ineligible branches are unresolved, never omitted. Use 8192 modes with16384-mode audit. Bounds are evaluated in floating point, not interval certified.

Separate prospective hypothesis classifications use maxima over standard AND refined N22: relative primary static C and K mismatch<=.01; all-four normalized scalar response mismatch<=.005; primary generalized centroid relative mismatch<=.005; generalized cosine spectrum mismatch<=.005; each reached-state primary current mismatch<=.01; Galerkin own-response whitened operator norm +1e-4 reserve<=.01 on each side; zero-delay first/second moments and signed insertion errors +1e-4 reserve<=.002. Each static, scalar response, centroid, cosine-spectrum and current metric includes an additive1e-4 numerical guard before comparison to its stated limit. Numerical eligibility is also reported separately; these reserves are not universal physical error bounds. Archive all sizes for descriptive trends; do not fit a convergence rate or extrapolate from them.

`pass` in assessment.json means numerical eligibility only. Every physical hypothesis has its own label; there is no aggregate GRUT pass. The retained numerical reserve for Galerkin and instantaneous moments exceeds the eligible tail estimates but is not a rigorous continuum error bound.

## Resource and execution contract

One frozen suite: five points, 300 seconds per point,1200 seconds total,4GiB Windows Job Object committed-memory ceiling per point, one BLAS thread. Historical resource-limiter negative checks remain in N22. The runner requires clean committed source, records source and checkpoint hashes, refuses output overwrite and preserves partial failure. A resource failure terminates the suite rather than changing the protocol.

    python run.py --inputs PATH_TO_HASHED_GROUNDS --output NEW_EXTERNAL_DIRECTORY

No new physical response or source-current result is evaluated before the freeze. Synthetic oscillator, finite N6 operator algebra and experiment ideal gate checks are implementation checks, not R1 physical production. EX1-01 is a separate prospective hardware protocol. It may establish a calibrated finite operational realization under standard quantum mechanics after real data; it cannot promote this numerical suite to a new empirical law.
