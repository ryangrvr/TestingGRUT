# D1-01: physical source dictionary and independent observable checks

Completed 9 November 2026. **The detector has an explicit calibrated source realization; its finite matching on the expanded backaction-state tests remains unresolved.** The new source dictionary and protocol were frozen at `b6cdfda3196f75eea0fb2ab951b42a8093d7ac12` before production. All registered points and the N22 refinement completed in 105.265 seconds. The post-run audit passed. 14 of 16 gates pass; the misses are `moments_and_insertions` and `finest_step`.

N22-01 remains its completed finite massless numerical pass. Its lattice headroom after guard remains 5.9324e-7, distinct from its combined-budget headroom 2.5910e-4. D1 measures different quantities with separately frozen tolerances; these errors are not history TVs. B2's N20 component miss, the withheld massive transfer and P1's interacting-parent obligations are unchanged.

## What the source derivation establishes

In the even neutral sector, write epsilon_i=(-1)^i and y_j=L_j-epsilon_j/4. At m_lat=m_R-g^2a/8, direct finite-chain algebra gives

    H_a = H_hop + (g^2a/2) sum_j y_j^2
          + (m_R/2) sum_i(1+epsilon_i Z_i) - g^2a(N-1)/32.

Thus y is exactly the electric-energy variable of the shifted Hamiltonian. A calibrated flux source y_j+q_j gives (g^2a)^-1 dH_y/dq_j=y_j at zero source. For the explicitly prescribed continuum cell source q(z)=sum_j q_j 1_Ij(z), the same normalized derivative is (1/a)integral_Ij Phi. This identifies matched generalized forces and fixes their normalization.

The conventional L_j+q_j source also contains a known staggered source-only offset. In a controlled detector it is a relative qubit phase and must be calibrated, along with the quadratic source self-energy and B2's boundary offset. DERIVATION.md specifies a qubit-projector coupling and apparatus Hamiltonian that leave exactly the interaction -u(t)P F_B2. Its ideal impulsive limit produces the historical K+/− instrument without probability fitting.

This gives physical meaning to the operational dictionary. Parallel source derivatives do not establish exact finite-a equality across different Hilbert spaces. Spatial interpolation, boundary placement, higher-derivative source corrections and finite-duration switching remain distinct matching questions. Fixed outer-cell terms are declared identity-valued apparatus calibrations, not continuum cell averages implied by endpoint Dirichlet data.

The derivation also displays the classical point-versus-cell term a^2 Phi''/24 and the freedom to add finite-a source improvements. These are not a derivation of all quantum matching coefficients. The full algebra, signs, normalization and limitations are in DERIVATION.md.

## Independent checks that were frozen

Four source-native probe functions were used: three sine modes and an asymmetric left beta profile, all with L1 norm1. Their cell weights came from direct integrals, independently of the quadratic reconstruction. The probes retain their signed integrals, including zero for sine2 and1/3 for sine3.

At N14,18,20,22, the study measured vacuum moments/cross-covariances; first and second probe moments and complex insertion matrices after one detector measurement at new times0,0.43,1.31; and all four outcome branches of two measurements at0,0.43 observed at1.31. The full unitary-path matrix elements were retained before conditioning. Mean kick and conditional-branch energies were checked against an independent lattice hopping identity and a continuum Parseval expression.

Historical ground-state checkpoints were reused and hashed. The new continuum reference was implemented separately from the historical history-probability code. Branch probabilities were used for normalization, not as a new history-TV gate. The analysis includes conservative mode-tail bounds for conditional field moments and conditional energy, mode doubling, and same-size N22 propagation/ground-state refinement. Every choice and threshold preceded production.

## What passes and what still fails at N22

| Quantity | Measured discrepancy or guarded value | Frozen limit | Decision |
|---|---:|---:|---|
| Maximum first/second moment or insertion error, plus reserve | 0.00407959128033 | 0.002 | fail |
| Maximum vacuum covariance error, plus reserve | 0.000967022411203 | 0.001 | pass |
| Largest guarded relative backaction-energy error | 0.00720463483806 | 0.02 | pass |
| Maximum N20-to-N22 observable step | 0.00101064958986 | 0.001 | fail |
| Maximum N22 same-size refinement | 8.82627304577e-15 | 1e-8 | pass |

The worst error is the sine1 second moment in the two-measurement `--` branch:

    lattice:   0.35162179953004
    continuum: 0.35560139081038
    absolute:  0.0039795912803346

The branch probabilities are 0.0156471245385 and 0.0156032297649, above the frozen0.001 floor. This error is not attributed solely to rare-branch normalization: the unconditioned twice-kicked sine1 second-moment insertion also differs by 0.00356777404344. The largest conditional first-moment discrepancy is 0.00361976917304, in the asymmetric beta probe's `+-` branch.

The finest-step miss occurs in the beta probe's `--` second moment, at approximately0.00101065. It is distinct from the largest absolute source-matching discrepancy. The refined N22 solve preserves both misses. None of the tolerances was changed after the result.

The numerical controls pass: covariance/energy tails and doubling, normalization, source identities and same-size refinement. Maximum conservative conditional-moment tail estimate is 5.2283795456e-06; conditional-energy tail estimate is 1.65943679084e-07. These are far below the failed physical tolerances, but are floating-point evaluations of analytic estimates, not interval certificates.

## Resolution and state dependence

| N | Largest moment/insertion discrepancy | Discrepancy / a^2 | Largest vacuum covariance discrepancy | Largest relative energy discrepancy |
|---|---:|---:|---:|---:|
| 14 | 0.008691813651 | 0.0782263229 | 0.002385726288 | 0.01660010617 |
| 18 | 0.005706556599 | 0.0824026773 | 0.001368520066 | 0.01047780737 |
| 20 | 0.004710981466 | 0.0831017131 | 0.001078060698 | 0.008612586547 |
| 22 | 0.00397959128 | 0.0842081515 | 0.0008670224112 | 0.007203348346 |

The errors decline with refinement; the scaled coefficient also drifts. This observed trend does not prove an asymptotic law and was not used to subtract a correction or to authorize another lattice size.

At N22, the largest immediate one-measurement probe discrepancy is about0.000816. At time0.43 it is0.001241; at1.31 it is0.001925 before reserve. After two measurements the maximum is0.003980. That pattern makes accumulated dynamics and backaction worth separating from static source matching. It does not, by itself, identify which discretization or boundary term causes the discrepancy.

## What is now justified, and the next discriminating work

The generalized-force/ancilla derivation is exact for the specified finite lattice and calibrated source. The vacuum covariance and mean energy injection support its normalization on the registered probes. The stronger accuracy claim for all registered backacted moments fails at this resolution and tolerance. This is not a contradiction of N22, a violation of its combined probability budget, or a falsification of GRUT.

A useful next study would separate equal-time source matching from evolution: derive and test equal-time covariance and mixed-source quadratic-form sum rules on the same reached states, then compare low-frequency or mode-resolved response at fixed source. This avoids assuming finite higher energy moments of the step-smearing kick. The resulting Hamiltonian, boundary and source improvement terms should be derived and fixed before another comparison. Do not tune detector weights or declare a universal correction from the failed conditional moments. No larger lattice or additional physical trial was run here.

The source interpolation and boundary convention are part of the stated physical experiment. A different derivation would define a new experiment and preserve N22 and D1 as historical records. Finite-duration switching is also open: the source realization above is impulsive. Step smearings have finite mean kick energy but generally divergent energy variance, so D1 deliberately uses quadratic-form mean energies and does not claim an H^2 convergence result.

## Verification, archive and reproduction

Independent N6 full-operator/source/ancilla checks passed, with maximum residual1.94e-14. An independent dense oscillator verified one/two-measurement Gaussian insertions and conditional moments to3.55e-16, with negative controls for a missing source calibration and an incorrect coherent-shift sign.

The post-run audit independently rebuilt all numerical decisions with zero discrepancy. It verified13 source hashes,5 checkpoint hashes, unchanged historical detector weights and960 state records, including matrix contractions against direct state moments. Tampered moments, a forged pass and a wrong input hash were rejected. It did not rerun the full lattice evolution or implement an independent full-size physical solver.

Source, derivation, frozen protocol, raw observables, resource records and audit are on codex/grut-source-dictionary-20261109. The delivered ZIP also contains the five historical binary inputs identified by INPUTS.json. Its source files preserve the exact Windows bytes recorded in provenance. SHA256.json covers all payload files except itself.

For production reproduction, use a clean Git checkout of the source freeze, Python3.12, numpy2.3.5 and scipy1.17.0 on Windows, and the hashed input folder. Follow PROTOCOL.md. For an offline audit, preserve the original archive and use a copy of raw with its audit.json omitted, then run:

    python audit_d1.py --root source --inputs inputs --output RAW_COPY

The audit refuses to overwrite its result. No massive transfer, N24 extension, merge, release or Zenodo deposit occurred.
