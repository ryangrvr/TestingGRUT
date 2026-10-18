# GRUT Candidate 1: explicit theory proposal and demonstrated regime

**Research report · 18 October 2026**

## Decision

**Outcome: a narrower working model rather than the intended universal theory.**

Candidate 1 specifies one leading covariant effective field theory and executes a populated classical/test-field/weak-gravity regime. It adopts Einstein gravity, the Standard Model with Dirac neutrinos, quantum probabilities, a four-dimensional metric and leading local propagation. It adds a complex charge field X, a real organization field phi and a real reporter D with explicit couplings.

It does not claim validated physics, a nonperturbative quantum-gravity completion, a quantum soliton lifetime, global fission stability, autonomous production or an implemented microscopic detector. Internal automated challenge is disclosed and is not independent review.

## Candidate law

With signature `(-+++)` and `hbar=c=1`,

\[
S=\int d^4x\sqrt{-g}\left[\frac{M_P^2R}{2}-\rho_\Lambda
+\mathcal L_{\rm SM+\nu_R}-\frac12(\nabla\phi)^2
-|\nabla X|^2-\frac12(\nabla D)^2-U\right]+S_{\rm boundary},
\]

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+g_X^2\phi^2|X|^2+\lambda_X|X|^4
+\frac12(m_{D0}^2+h_D^2\phi^2)D^2+\frac{\lambda_D}{4}D^4.
\]

The postulated global `U(1)_X` conserves charge. High X occupation can deplete phi, reducing the X mass and binding the populated configuration; the altered phi profile changes D propagation. The structural profile is occupation-correlated and reversible. A historical memory surviving removal of X is not demonstrated.

The displayed action is a leading truncation. All symmetry-allowed counterterms, portals, nonminimal curvature terms and higher operators are admitted with an explicit matching prescription. Zero matched values are choices, not radiative protection. See [the full candidate](GRUT_CANDIDATE_1.md) for the state-space/observable prescription, constraints, time convention, field equations, preparation, measurement rules and calibrated-input definitions.

## Executed operating regime

Use `g_X=.1`, `lambda_phi=.01`, `lambda_X=.001`, `h_D=.01`, `m_D0=g_X v`, `lambda_D=.001`, and `v/M_P=10^-5`. The asymptotically flat benchmark sets renormalized `rho_Lambda=0`; this is not a cosmological prediction.

Prepare `X=e^(-i omega t)f(r)/sqrt(2)`, `phi=varphi(r)` and `D=0`. With `rho=g_X v r`, `F=f/v`, `S=varphi/v`, `Omega=omega/(g_X v)`, `a=1`, `b=.1`, solve

\[
F''+2F'/\rho=(S^2+bF^2-\Omega^2)F,\qquad
S''+2S'/\rho=[a(S^2-1)+F^2]S.
\]

The supplied charged preparation is regular, nodeless and localized. Six failed preliminary solver starts were retained. A successful bag-seed preparation was then frozen before four verification refinements; physical couplings never changed.

| Calculated result | Value | Claim boundary |
|---|---:|---|
| Classical/mean charge | 39,943.8462 | Not an exact integer-charge coherent state |
| Energy `E/v` | 3,639.2300 | Supplied initial matter resources |
| Energy per charge relative to free X mass | 0.91108652 | Blocks complete free-carrier dissolution at leading order; does not prove every fission comparison |
| Core `phi/v` | 0.049079979 | Organization response in this preparation |
| Energy radius `R90 v` | 51.03432 | 90% energy radius |
| Compactness `2GE/R90` | 5.67463×10^−10 | Weak gravity only |
| Central clock redshift | 6.68813×10^−10 | Full-stress first-order Einstein response |
| Reporter phase shift at momentum `5 g_X v` | 0.0040665116 radians | Elastic S phase is twice this number |
| Ideal monochromatic plus-output probability | 0.49593353 | Supplied balanced instrument with reference phase pi/2 |

The fixed-charge amplitude Hessian was checked in angular sectors 0–3 and the phase Hessian in sectors 0–2 across three grids. Identified translation and phase zero-mode errors decrease as mesh spacing squared. The positive radial finite-box level is not a universal gap. These results support local classical energetic stability; nonlinear persistence, every charged partition and quantum decay remain uncomputed.

Gravity uses the same profile's energy and anisotropic stresses. In dimensionless form `laplacian_rho Phi=4 pi Gv²(e+p_r+2p_t)` and `laplacian_rho Psi=4 pi Gv²e`. Stress conservation and the virial identity pass. A flat scalar profile plus its first-order metric is a weak-gravity expansion, not a solved nonlinear Einstein–scalar configuration.

The reporter phase was calculated using independent variable-phase and radial-wave equations. A finite packet was also evolved on three grids. At the finest grid the selected outgoing probability is about .992, joint plus probability about .492, and conditional plus probability about .49591. Their modest precision and separate normalization are retained; packet delivery is not memory retention. A microscopic source, splitter, recombiner and detector remain an explicit operational gap.

![Calculated scalar, stress, metric and reporter observables](REGIME_FIGURE.png)

## Compatibility with archived GRUT work

This is a **separately named successor**, selected over continuing the moving-word architecture. It replaces that mechanism rather than claiming to resolve its moving-channel bottleneck. No D1 finite-density, B1 positive-electric-energy, B2 fixed-path or MC-E1 moving delivery/retention theorem is inherited. X is not electric charge, and this is not a current-loop calculation.

Dimension, locality, time orientation, Born probabilities and Einstein dynamics are supplied openly. Their adoption is not relabeled emergence. The underlying charge-supported scalar mechanism is known prior work; its quartic extension and this parameter point have their own calculations and obligations. No originality claim is made for general relativity, quantum theory, Standard Model content or scalar solitons.

## Reproduction and verification

From this directory:

```bash
python code/run_all.py
```

The package includes source, pinned dependencies, prescribed preparations, failed attempts, raw results, convergence grids, execution notes and challenge records. All **25 internal numerical acceptance checks passed**, including an independent sampled check of the radial equations. This verifies the implemented leading regime, not the theory's realization in nature. See [code/README.md](code/README.md) and [verification evidence](code/verification.json).

No extra package installation or paid computation was performed. Runs are finite and require no account or background service. [EXECUTION_RECORD.md](EXECUTION_RECORD.md) records the original construction phase; [PUBLICATION_NOTE.md](PUBLICATION_NOTE.md) distinguishes this subsequent authorized repository posting from that earlier snapshot.

## Complete report packet

| Document | Content |
|---|---|
| [GRUT_CANDIDATE_1.md](GRUT_CANDIDATE_1.md) | Complete leading-law proposal, postulates, states, constraints, action, measurements and parameters |
| [WORKED_REGIME.md](WORKED_REGIME.md) | Preparation through field solution, reporter propagation, ideal records and gravitational observables |
| [COMPATIBILITY.md](COMPATIBILITY.md) | Integrated defined/derived/numerical/conjectured/missing table and challenge verdict |
| [PREDICTIONS_AND_TESTS.md](PREDICTIONS_AND_TESTS.md) | Calibration, baselines, held-out test design, rejection conditions and errors |
| [OPEN_CLAIMS.md](OPEN_CLAIMS.md) | Quantum/global stability, formation, apparatus, UV and empirical gaps |
| [audit/](audit/README.md) | Constructor and Challenger notes, including corrected and surviving objections |

The next decisive research test is a finite-charge branch/partition study and controlled nonlinear perturbation evolution for this same potential, followed by renormalized quantum stability. A failure may reject the operating regime. Better phase precision alone will not close the physical-apparatus gap.
