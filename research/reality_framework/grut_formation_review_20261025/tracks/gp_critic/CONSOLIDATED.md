# Consolidated review: Codex b5174e5 at the 7737487 gauge point (track gp_critic)

25 October 2026 · Claude (issue #13) · consolidates tracks gp_matching, gp_phase, gp_stability, gp_walls and gp_leftovers.
Codex packets were read only through `git archive` of b5174e5 (`grut_gauge_point_test_20261025/`) and 7737487
(`grut_z6_gauge_bridge_20261025/`). Own code is in `code/`, outputs in `results/`. No Codex code was run and nothing was
committed.

**This is one supplied tree-level working model. Nothing here is evidence that X, Y, Sigma, phi or the inflationary
inputs exist. It is not a validation, an exclusion or a GRUT prediction.** Comparisons with CMB, LVK and PTA data are
retrospective. The five tracks and this one are same-session AI work, not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Stated point [Q]: M_Y = 1e12, h = 1e11, mu = 5e11 GeV, kappa (Codex "y") = 1e-3, e = 1e-7, lambda_Y = lambda_Sigma = 1;
light target g_X = 0.1, lambda_X = 1e-3, lambda_phi = 0.03, v = 0.924089637 GeV, Q = 200000. The review inputs are
H_inf = 4.7e13 GeV, c_H = 1, the Hubble mass switched off at release, and instant reheating.

## 1. Verdicts on Codex's five handoff items

| # | Codex claim (b5174e5) | Consolidated finding | Verdict |
|---|---|---|---|
| 1 | Early matching: M_Y/H = 0.0213; kappa abs(X)^2/M_Y^2 about 1100; full minimum 1.1197e15 / 8.8895e13 / 1.2548e13 GeV; Higgs 177.5 times its vacuum value; phase masses 0.378H, 2.751H | Reproduced to 6 digits from 48 random complex starts with a full 6x6 Hessian (gp_matching) and to 9 digits by gp_phase [C] | CONFIRMED |
| 1 | Full potential bounded; truncated runaway is not a theorem | Global search plus the exact Y/Sigma-eliminated potential; lambda_eff = 8.81e-4 [C, D]. Stability survives one-loop running, but lambda_Y reaches 4 pi at 1.43e15 GeV, close to the inflationary abs(X) = 1.12e15 GeV [E] | CONFIRMED; perturbative control is marginal |
| 1 | "The phase masses do not supply a CMB isocurvature prediction" | The light mode is the source phase psi2 (2.654 rad per canonical H). A linear estimate gives S of about 0.036 cot(delta0) [E]. The chaotic release (item 2) overrides any linear map | CORRECTED |
| 2 | Real couplings plus own-potential alignment give exactly zero charge; the ±pi/24 releases and Gauss law | Proof [D]. Comoving charges (-16.261387466, 6.395904449, -0.487720980) reproduced to 1.2e-11 with an independent integrator [C] | CONFIRMED |
| 2 | The physical reading "own alignment supplies zero charge" | Zero mean only. The real subspace is unstable: growth of about e^29-e^30 by a of about 130-150 in both tracks. Local abs(eps_X): rms 0.14-0.23, median 0.11, maximum 0.81, sign chaotic. Inflation does not hold psi2 at its minimum (<cos psi2> = 0.26) [C] | CORRECTED |
| 2 | (our hypothesis) A moving minimum from higher operators supplies the charge | 67 CP-violating configurations, 8,832 runs: abs(<eps>) ≤ 0.034, chi^2/dof 0.94. DRT-type supplied mismatches also average to zero [C] | REFUTED at this point and release |
| 3 | Screened bound 6.56137e-20 GeV = 4.63e-23 of the 1417.26 GeV margin; m_A = 6eh; heavy masses exceed the object's energy; Q = 2 mod 6 | All reproduced with an independent profile and a real-space Yukawa route [C]. New unscreened bound: 3.40e-7 GeV = 2.4e-10 of the margin [C] | CONFIRMED |
| 3 | "A compact population needs compensating sectors; three cores" | The invariant is (Q_X + 3Q_Y) mod 6, on the whole closed slice only. Compensation can sit anywhere; the rule selects nothing [D, A] | CONFIRMED with clarification |
| 4 | Tilt Delta/4 [2 - 3x + x^3]: stationary ±v, energies 0 and Delta, positive potential, extra core energy 1.74e-17 GeV | Exact (sympy); the tilt equals Saikawa eq. 2.23 plus Delta/2; the first-order tension shift is exactly 0; 1.7382e-17 GeV [D, C] | CONFIRMED |
| 4 | Minimum annihilation temperature 4.26 MeV, above 1 MeV across controls | 4.2607 MeV reproduced, but it uses g* = 61.75 near 4 MeV. With g* = 10.75 the minimum is 6.596 MeV [C]. The central range is 13.9-24.3 MeV with a self-consistent g* (gp_walls), or 13.9-24.9 MeV in my fixed-g* recompute [C]. This holds only for t_ann ∝ 1/Delta. Readings of Babichev et al. span 5 keV-18 GeV, and wall–object capture could delay to about 50 keV [E] | CORRECTED (Codex's figure is the conservative one); the pass is conditional |
| 5 | Pair budget Y_pair ≤ 2.381e-11 (0.0049296 per net charge); eps ≥ 0.99175; quench debris 0.34-4.2 | Arithmetic and the six debris values reproduced exactly [C]. The 0.9(eps^-4/3 - 1) law overstates the frozen remainder by 1.2-4.2 times; thresholds are 0.990 (non-relativistic count) and 0.988 (adiabatic) [C] | CONFIRMED (arithmetic); law CORRECTED |
| 5 | "An existing ledger may already include the quench energy" | The ledger has no V0. Debris, binding pairs and remainder share one onset energy of 0.867 m_X/eps per net charge; adding them double counts for eps > 0.50 [C] | CORRECTED |

## 2. The six requirements at this gauge point

| Requirement | Codex | Consolidated | Basis |
|---|---|---|---|
| Potential stability | PASS (tree) | **PASS** at tree level and in a one-loop estimate. Marginal perturbativity near the inflationary amplitude; no matching scale is stated | gp_matching [C, E] |
| Fixed initial phase | NOT SUPPLIED | **FAIL at this release.** Supplying a phase does not help: the release is chaotic at a sensitivity of 1e-9 to 1e-15 rad, far below inflationary phase spreads. Super-horizon patches at a ≈ 100-300 therefore get independent signs. The domains are 30-90 m today and 150-450 GeV^-1 at an 86 GeV onset, against an object spacing of 116 GeV^-1. The ± charge is mixed on the scale of the objects themselves | gp_matching, gp_phase [C]; separate-universe argument here [A, E] |
| Bias to remove walls | ABSENT (variant passes) | **OPEN.** The algebra passes [D] and the timing is 6.6-39 MeV under the standard law [C]. The objects attract walls by about 1.40 TeV (99% of the binding margin), and wall-bound X has a gap of 0.828 m_X, below omega = 0.861 m_X. Capture dynamics are not computed | gp_walls |
| Symmetric leftover matter | OPEN | **FAIL on present estimates at H_inf = 4.7e13 GeV.** Net eps ≤ 0.03; the frozen remainder needs ≥ 0.990; the best equilibrium case (homogeneous) needs ≥ 0.19. Erasure is possible only through lumps, conversion and uncomputed 3-body capture, and only for eps ≥ 0.01. Two extra sources are each sufficient on their own [E]: the longitudinal vector (1.06e3 pairs per net charge) and the Y condensate (below) | gp_leftovers [C, E]; this track [E] |
| Gauged object stability | PARTIAL | **PARTIAL** (confirmed). The classical margin is huge whether screened or not. The fluctuation spectrum, quantum metastability and the wall interaction remain open | gp_stability, gp_walls |
| Corrected survival rate | PASS (included channels) | **PASS** for included channels. c6 = 2.1068e5, rate 4.618e-36 /s, rate times age 2.011e-18 [D, C]. The late c6 is robust to the dimension 5-7 operators that gp_phase lists (shifts ≤ 8.2e-4) [D] | this track |

Overall: two pass, two fail, two remain open or partial. **On homogeneous calculations plus estimates, the stated
point is not a viable genesis point under the supplied release.** A lattice calculation could change the details. It is not expected to
change the absence of a resolved net charge: abs(<eps>) ≤ 0.03 even with CP-violating operators [C], and super-horizon
signs are independent by the separate-universe argument [A].

## 3. Cross-track consistency and overreach

**Agreements.**
- Two independent integrators find the same release instability: gp_matching's rate of 0.263 per conformal time and
  gp_phase's envelope of 0.20 come from different fit windows, and both give about e^29-e^30 of total growth. The
  saturated spreads also match (median 0.11 against rms 0.177).
- Three object-profile solvers agree: Omega = 0.8611951 and margin 1417.255 GeV, with gp_stability and gp_walls
  each matching Codex to 3e-10 and 4e-11.
- The longitudinal-vector abundance matches: gp_matching 0.3 GeV and gp_leftovers 0.31 GeV, against my GMR
  re-evaluation of 0.307 GeV and 1061 pairs per net charge [C].

**Flags.**
1. *gp_matching, isocurvature "about 4,000 times S_max".* The linear formula assumes a smooth map from phase to
   charge, and the same track shows that map is chaotic. Read the figure as the linear-response size only. The
   operative problem is the stronger one in section 2: no correlated net charge forms.
2. *gp_phase, leftover figures.* They use the overstated 4/3 law. At abs(eps) = 0.17 and 0.03, the frozen excess is
   about 5.0 and 33 (non-relativistic count), not 8.6 and 95. The conclusion is unchanged.
3. *gp_phase, one-field "upper envelope 0.22".* The family omits the compulsory radial companion. For this
   construction, |beta6|/(2|alpha6|) = M_Y^2/(sqrt2 mu h) = 14.1, and it exceeds 1 whenever the Y mass matrix is
   positive [D]. The envelope applies only to the angular-only family. Our affleck_dine window (c6 from 137 to 1,610)
   omits the companion as well, so neither can be imported here (section 4).
4. *gp_stability, strings "plausible for high-scale fast reheating".* Visible reheating cannot restore the hidden
   U(1): the portal coupling of about 9e-10 keeps Gamma/H far below 1. During inflation the U(1) is broken by
   abs(X) of about 1.1e15 GeV and by abs(Sigma) at 177 times its vacuum value. Restoration would need the hidden
   sector's own energy to thermalize (T_h ≤ 2.2e14 GeV at release if it did [E]), or a non-thermal restoration
   during the release resonance. Gmu = 3.12e-15 [C] lies below the LVK O3 headline Gmu ≲ 4e-15 (abstract checked
   [Q]). The model C-1 range is quoted from gp_stability and not re-verified. String formation is open in both
   directions.
5. *gp_leftovers.* Its maximal-equilibrium thresholds assume full thermalization, which is optimistic. It does not
   include the Y-condensate source that gp_matching flagged.
6. *gp_walls.* The 50 keV snowplow is labelled crude and should stay so. Otherwise it is scoped correctly.

## 4. New consolidation results (this track)

**a. Survival recomputed [D, C].**
- alpha6 = 3.5533e-32 GeV^-2 and beta6 = -1.0050e-30 GeV^-2, so c6 = alpha6 M_red^2 = 210683.18, matching Codex.
- With Codex's unit-coefficient rates (coherent 3.92e-50 /s, pair 1.0400e-46 /s [Q]; our 2e9f296 total is
  1.04e-46 [Q]), the rate is 4.6179e-36 /s and rate times age is 2.0111e-18. The 1% ceiling on c6 is 1.486e13.
- Planck-suppressed operators shift the late alpha6 by 8.2e-4 (Sigma* Y X^3/M_red) and 1.4e-13 (Sigma* X^6/M_red^3).

**b. Y-condensate source [E, built on gp_matching C].**
- The Y condensate is excited even on the aligned, zero-charge trajectory. Its comoving energy grows linearly
  (0.733 a, in H_inf^4 units) until Y decays at a ≈ 3.5e5 (Gamma ≈ 378 GeV).
- At decay, E_Y/E_X ≈ 371, which is 3.2e-5 of the visible energy.
- If it never thermalizes: 0.50 pairs per net charge, about 100 times the 1% budget (f ≈ 1.0).
- If it thermalizes: about 180 units of hidden entropy per net charge, in gp_leftovers' failing regime (eps_eff ≲ 0.003).
- Both numbers assume the best-case uniform-sign charge, which alone would be 248 times Y_Q.

**c. Longitudinal vector [D, E].**
- m_A/H = 3.47e-6 at the inflationary point. The GMR estimate gives rho/s = 0.307 GeV and 1061 pairs per net charge.
- The GMR premise of a light vector during inflation fails once m_A(inflation) ≥ H, which needs e ≳ 0.029 at this
  point [D]. Other consequences of a larger e are not checked.

**d. Does the companion help in the regular regime? [D, C]**
- In the truncated EFT, the phase-independent companion caps the Hubble-held minimum at
  eta = |alpha6| X0^2/lambda_X ≤ r/12 = 5.9e-3 (at phi0 = pi/2), i.e. H_inf ≤ 5.8e11/sqrt(c_H) GeV [D].
- My own full three-field minimization shows this "edge" is an artefact of truncation [C]. The full minimum
  continues smoothly to 4.7e13 GeV with abs(X)/X0 = 1.04-1.07. By the EFT edge, lambda_Y abs(Y)^2/M_Y^2 ≈ 1.1, so the
  quadratic elimination has already failed there.
- Where the EFT is reliable (H ≤ 0.5 of the edge, lambda_Y abs(Y)^2/M_Y^2 ≈ 0.1), my release from that minimum with
  the best supplied phase gives eps = 0.010, 0.013, 0.024 and Y = 0.013, 0.047, 0.42 Y_Q for c_H = 1, 10, 100 [C].
- Controls: no-companion A_6 = 12.85 (affleck_dine 12.97); tolerance 1.5e-9; a = 1e4 at 3e-5; exact CP mirror.
- Closer to the edge the EFT gives eps ≤ 0.10-0.12 and Y up to 12 Y_Q, but those values sit outside its validity.
- At the GMR ceiling of 2e10 GeV, the best yield is 0.12 Y_Q (c_H = 100).
- **So no regime computed so far gives both Y ≥ Y_Q and eps ≳ 0.19.** The range from 0.5 of the edge up to
  4.7e13 GeV has been computed only at its top, where the release is chaotic.

## 5. Corrections for Codex to apply

1. Replace "own alignment generates zero charge" with: zero *mean*; the real subspace is unstable; local
   abs(eps) ≈ 0.1-0.2 with random sign; inflation does not prepare the aligned state.
2. Change "fixed initial phase: NOT SUPPLIED" to **FAILS at this release**. Neither a supplied DRT phase nor a
   moving minimum from dimension-5 operators gives a resolved net charge.
3. Drop "masses do not supply an isocurvature prediction". State the dictionary: the light mode is the source phase psi2.
4. State the phase count: 0 physical phases at dimension ≤ 4 with Z2_phi; 8 / 16 / 53 at dimension ≤ 5 / 6 / 7.
   The proposed wall tilt breaks Z2_phi, which makes phi Sigma* Y^2 renormalizable with one more phase. Re-audit
   kinetic mixing and CP accordingly.
5. Mod-6: use (Q_X + 3Q_Y) mod 6 for the whole closed slice, and drop the three-core grouping.
6. Add the Z6 string tension (Gmu ≈ 3.1e-15, conditional on formation) and its formation condition.
7. Walls: the minimum T_ann is 6.6 MeV with a consistent g*; the pass holds for the standard law only; credit
   Saikawa eq. 2.23; list wall–object capture and wall-bound X as open.
8. Leftovers: the ledger contains no V0, and the energies are alternatives, not additive; use thresholds of
   0.990/0.988; add the longitudinal-vector and Y-condensate sources.
9. Matching: state a matching scale (lambda_Y runs to 4 pi near 1.4e15 GeV). Note that the compulsory companion
   (|beta6| ≥ 2|alpha6|) forbids importing any angular-only AD window.

## 6. Status at the stated point

- **Pass:** potential stability (tree level and a one-loop estimate); cold classical binding with gauge screening;
  included-channel survival (rate times age 2.0e-18); the mathematics of the bias; the arithmetic of the mod-6 rule.
- **Fail (supplied release, homogeneous plus estimates):** net charge generation (chaotic, mean zero); leftover matter
  at H_inf = 4.7e13 GeV.
- **Open:**
  - wall annihilation with objects present;
  - string formation and its LVK-edge tension;
  - the coupled fluctuation spectrum and quantum metastability;
  - hidden-sector thermalization and erasure through lumps;
  - longitudinal-vector production in the real Higgsing history;
  - Y-decay CP and Y-decay products;
  - the H_inf range between the regular and chaotic regimes.

## 7. The single most decisive next calculation

**A full three-field (X, Y, Sigma) homogeneous release scanned in H_inf, from the EFT-reliable range
(about 3e11/sqrt(c_H) GeV) up to 4.7e13 GeV, at c_H = 1 and 10.** For each H_inf it should record:
- the net eps, and its change under a 1e-9 rad shift, as the chaos indicator;
- the yield against Y_Q;
- the Y-condensate energy at a = 1e3;
- the eaten-mode amplitude, as input to a longitudinal-vector Bogoliubov estimate.

The existing two-track code can run it in hours. The calculation decides whether this construction has any release
that is regular, productive and symmetric-safe:
- If no H_inf gives regular release with Y ≥ Y_Q and eps ≳ 0.19 (or ≳ 0.01 together with a demonstrated
  erasure through lumps), the one-mediator "creates and protects" link fails at these couplings, whatever the walls,
  strings or spectrum turn out to be.
- If a regular window appears, the longitudinal-vector population at that H_inf becomes the next gate.

## 8. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/gp_critic/code
python3 critic_checks.py ../../gp_matching/results/release.json ../results/critic_checks.json   # 1 s
python3 companion_ad.py ../results/companion_ad.json        # ~8 min (282-run phase/eta scan, controls, yields)
python3 regular_yields.py ../results/regular_yields.json    # ~3 min
python3 full_minimum_vs_H.py ../results/full_minimum_vs_H.json   # ~5 min
```

Credit: Affleck–Dine; Dine–Randall–Thomas; Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki (Q-ball
formation); Graham–Mardon–Rajendran (vector from inflation) [R]; Saikawa 1703.02576 and Babichev et al. 2504.07902
(walls) [Q, via gp_walls]; LVK O3, arXiv:2101.12248 (abstract read) [Q]; Krauss–Wilczek (discrete gauge charge) [R].
