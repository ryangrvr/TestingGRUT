# Playing out the four next steps for GRUT

This records the original four-step probe and its later DMRG rerun. **Scope corrections by Codex, 7 November 2026:** original registrations, solver physics, outputs and missed windows are preserved. The initial ED pre-run chronology was reported by the builder; its registration and outputs were first committed together. The DMRG recipe and analysis have separate earlier commits. See `RESULTS.md` and the [review](../grut_qed_review_20261107/REVIEW.md) for the independently checked scope.

## Short answer

- **Interacting charged matter fits GRUT's detector-history framework.** In the lattice Schwinger model (1+1D QED), GRUT's two-time histories behave as its effective-description theorem predicts. All three preregistered history checks pass.
- **The initial ED spectrum check missed; the registered DMRG rerun passes both windows.** Historically, The zero-mass extrapolation lands 3.3% off the exact value against a 3% window. The massive holdout lands 9.9% off. Exact diagonalization at up to 24 sites is too small to reach the continuum to 3%. This is a limit of the method, not evidence against GRUT, but it still counts as a miss.
- **Finite-spacing dispersion is model dependent.** The massless Wilson and plaquette dispersions differ away from their common low-momentum limit. A conditional soft-Cherenkov estimate gives a much smaller illustrative spacing than the photon-only estimate; neither is a universal GRUT bound. See the assumptions below.
- **The chiral and spin-2 steps are stated as obligations, not solved.** No derivation was attempted. Below I write them in the form the GRUT program could adopt.

## Step 1: lattice QED in GRUT's framework

**Model.** The model uses Kogut–Susskind staggered fermions with U(1) links and open boundaries. Gauss's law is solved exactly, and θ=0. The Hamiltonian is \(W = x\sum(\sigma^+\sigma^-+{\rm h.c.}) + \tfrac{\mu}{2}\sum(1+(-1)^n\sigma^z) + \sum L_n^2\), with \(x=1/(ag)^2\) and \(\mu=2\sqrt{x}\,m/g\). The reference values come from the [exact massless result and Byrnes et al. (DMRG)](https://ar5iv.arxiv.org/html/hep-lat/0202014): vector binding energy 0.5641896 at \(m=0\) and 0.51918 at \(m/g=0.25\).

### Detector histories (P-C, P-D, P-E)

The instrument is GRUT's two-outcome Weyl instrument, \(K_\pm=(1\pm e^{i\varphi}W)/2\) with \(W=e^{icF}\). Here c=2, and F is a Gaussian-smeared electric field (σg=0.8 in a box with Lg≈5). Histories are computed at two times, out to T=6/g. The statistic is the maximum difference between the exact history probabilities and a quasi-free prediction built only from two-point functions.

| Lattice | \(m=0\) deviation | \(m/g=0.5\) deviation |
|---|---:|---:|
| N=10, x=4 | 0.00763 | 0.00788 |
| N=14, x=9 | 0.00551 | 0.00918 |
| N=20, x=16 | 0.00412 | 0.00895 |

- **P-C passes its finite-ladder trend diagnostic.** The Gaussian residual shrinks along the tested ladder. This is consistent with the massless bosonized limit, but does not establish continuum detector matching: the reference uses the lattice's own mean and covariance, and the box length is not exactly fixed across all three points.
- **P-D passes its declared finite comparison.** The selected Gaussian summary has a larger residual for the massive case on these three lattices, only slightly at x=4. Two exploratory matched-variance checks in the Codex review retain a larger massive residual. This supports information beyond this Gaussian summary for the tested histories; it excludes neither every reduced description nor an untested continuum behavior.
- **P-E passes.** The first instrument disturbs the second: by about 0.04 at \(m=0\) and 0.012–0.015 at \(m/g=0.5\). Probabilities stay normalized to machine precision.
- **Caveat.** The two masses have different field fluctuation sizes (Var(cF) ≈ 0.24 vs 0.12). The comparison is therefore not perfectly matched, though the smaller fluctuations at \(m>0\) should make deviations harder to see, not easier.

### Spectrum (P-A, P-B)

**Recipe.** I fixed the recipe in `ed/extrapolate.py` using only the \(m=0\) data, before the holdout run finished.

1. At each x, fit \(\omega^2=M^2+b/N^2+c/N^4\) on N=14–24.
2. Extrapolate in \(1/\sqrt{x}\) (quadratic) over x=2–32.
3. Report the min–max spread over 24 fit variants as fit-choice sensitivity, not a complete uncertainty band.

| Prediction | Target | Result (variant spread) | Error | Window | Verdict |
|---|---:|---:|---:|---:|---|
| P-A, \(m=0\) | 0.5642 | 0.583 (0.561–0.641) | +3.3% | ±3% | Narrow miss |
| P-B holdout, \(m/g=0.25\) | 0.5192 | 0.571 (0.523–0.609) | +9.9% | ±3% | Miss |

The infinite-length values decrease steadily with x but are still far from the continuum at x=32. Examples: \(m=0\) runs 0.785 → 0.613; \(m/g=0.25\) runs 0.737 → 0.582. At large x, the bound state spans more lattice sites than a 24-site chain can hold. The qualitative trend is right, since binding drops when mass is added, but the size is not. The published values were reached with DMRG/MPS on chains hundreds of sites long. That was the limitation of the initial ED implementation; the subsequent MPS implementation and rerun are now in this packet.

### Update: sealed DMRG rerun passes both

A DMRG/MPS solver was validated against ED (rel ≤ 2e-13), and its recipe was sealed before any production run (`mps/PREREGISTRATION_AMENDMENT_01.json`, commit `dd1e107`). It then reached chains of up to N=192 (physical length 24/g). P-A gives 0.5641906 against 0.5641896, and P-B (the previously examined massive benchmark) gives 0.519232 against 0.51918. Both are inside ±3%, with every fit variant within 0.25%. The ED misses above stay on record as the first test. See `RESULTS.md`.

## Step 2: the chiral obligation

The probe tests a finite-band Hermitian lattice Hamiltonian with smooth periodic momentum dependence and isolated Weyl nodes. In this setting the total Weyl charge on the Brillouin torus vanishes. The broader Nielsen–Ninomiya obstruction also depends on its locality, regularity and conventional chiral-symmetry assumptions; it is not a ban on every lattice realization of chiral physics. QED is vector-like and does not solve the Standard Model chiral-gauge problem.

[Ginsparg–Wilson symmetry](https://arxiv.org/abs/hep-lat/9802011) changes the chiral realization and therefore evades the naive hypotheses. Overlap, domain-wall and suitable interacting constructions require their own anomaly, locality and continuum checks. A complete chiral SM construction has not been supplied here. A declared effective theory may adopt a consistent regulator; deriving a unique choice from access or closure is a separate, stronger claim.

## Step 3: universal spin-2 coupling as the target

The calculation compares specified **linear, unscreened** exchange theories after Newtonian calibration. The massless universal spin-2 comparator has PPN gamma=1; the linear Fierz–Pauli small-mass limit has gamma=1/2 and 3/4 of GR's bending. The pure trace-coupled scalar comparator has no leading coupling to a traceless photon source. These alternatives are inconsistent with the quoted Cassini PPN result in their stated regime.

This does not exclude every massive or scalar–tensor gravity model. Nonlinear [Vainshtein screening](https://arxiv.org/abs/1304.7240) can invalidate application of the linear massive comparator to the Solar System. Cassini measured a frequency shift constraining PPN gamma, not the graviton's identity. Universal coupling, nonlinear constraints and dynamical geometry remain adopted or unconstructed here.

The stress-tensor part of [Weinberg–Witten](https://doi.org/10.1016/0370-2693(80)90212-9) concerns massless spin greater than one carrying nonzero energy-momentum in a Lorentz-invariant theory with an appropriate Lorentz-covariant conserved stress tensor and well-defined matrix elements. Its hypotheses must be checked in a proposed microscopic model; abandoning one hypothesis does not itself construct gravity. GR's local gravitational energy does not satisfy the theorem's covariant-stress premise in the required way.

## Step 4: spacing ceiling and shared light cone

`lightcone.py` expands both dispersion relations to \(O(k^2a^2)\), with \(v=1+C(ka)^2\):

- **Wilson fermion (r=1):** \(C_e=3/8-S/2\)
- **Plaquette photon:** \(C_\gamma=-S/8\)

Here \(S=\sum n_i^4\). Both formulas check numerically to about \(10^{-5}\).

- **Different finite-spacing group velocities, common limiting speed.** The difference is $C_e-C_\gamma=3(1-S)/8$, zero on axes and positive elsewhere. In the sampled directions, 77.6% have $C_e>0$ relative to the limiting speed. This is not a calculation of an exact finite-lattice front velocity.
- **Conditional estimate only.** The archived $a\lesssim3\times10^{-32}$ m body-diagonal value assumes a physical preferred lattice, its orientation, the stated Wilson and plaquette dispersions, a soft vacuum-Cherenkov emission model, the inverse-Compton interpretation of Crab gamma rays, and removal of the leading mass-dependent speed shift. LHAASO observes photons; the electron energy is inferred. This is not a bound on a regulator removed in the continuum limit.
- **The mass term matters.** For positive untuned Wilson bare mass, $v-1=-m^2/(2p^2)+ma/2+C_e(ap)^2+\cdots$. At the archived illustrative bound, $ma/2$ is 1.55 times the retained mass deficit. The old threshold cannot be used without stating the corresponding matching/tuning assumption. Hard-photon kinematics and radiative stability require separate treatment. See the [explicit derivation](../grut_qed_review_20261107/REVIEW.md).

## What this adds up to

| Step | Current scoped outcome |
|---|---|
| 1. Interacting charged matter | Three finite history diagnostics pass. The original ED spectrum misses remain; the registered DMRG rerun passes both spectrum windows. Interacting continuum histories remain unmatched. |
| 2. Chirality | Known conditional obstruction and regulator options identified; full chiral SM not constructed. |
| 3. Gravity | Specified linear comparators checked; universal coupling adopted and nonlinear recovery open. |
| 4. Propagation | Finite-spacing dispersion differences and conditional estimates; no universal GRUT spacing bound. |

The next connected calculation is independent detector-history matching at fixed physical parameters, states, boundaries and instruments, beginning with the massless interval reference. A further spectrum fit alone cannot establish that bridge.

## Scope and review record

These scripts and results are committed in draft PR #14. The Codex review reran four solver checks, independent finite-Hamiltonian checks, the spectrum analysis and one production point; it did not rerun the complete production grid or independently implement DMRG. The matched-variance follow-ups were exploratory. Broader external mathematical and physical review remains open.

The original ED misses and DMRG registration chronology are retained. The massive DMRG benchmark had already been examined and is not an untouched holdout. Twelve-fit spreads measure fit-choice sensitivity. See `RESULTS.md` for the repaired clean-checkout replay command and [the research report](../grut_solid_theory_research_20261107/RESEARCH_PASS.md) for the matching obligations. Zenodo remains paused.
