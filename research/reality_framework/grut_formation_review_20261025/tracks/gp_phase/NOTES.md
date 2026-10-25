# Track gp_phase: handoff item 2 (phase) at the 7737487 gauge point

25 October 2026 · Claude (independent review track) · reads Codex b5174e5
(`origin/codex/grut-gauge-point-test-20261025`) and its parent 7737487
(`origin/codex/grut-z6-gauge-bridge-20261025`) through `git archive` only. Own code in `code/`, outputs in
`results/`, large scratch output in the session scratchpad. No commits. Codex's code was not run.

**This is a conditional calculation inside one model with supplied fields and a supplied release history. It is
not a GRUT prediction, and nothing here is evidence that X, Y, Sigma or any inflationary input exists.**

Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Questions and inputs

Handoff item 2 (b5174e5 CLAUDE_HANDOFF.md): "Constant real couplings plus own-potential alignment remain in a
real-field subspace and generate zero charge. A short full complex release verifies this, its conjugate alternatives
and compact Gauss conservation. The inflaton-dependent mismatch assumed in the review is not present in 7737487." [Q]

Tasks: (a) prove the statement and check it with our own complex-field release; (b) count the physical CP phases
of the 7737487 EFT and of its next order; (c) test whether a generic EFT supplies the charge without an assumed
initial mismatch, through a potential minimum whose phase moves as the radial amplitudes evolve.

Stated point [Q]: M_Y = 1e12 GeV, h = 1e11 GeV, mu = 5e11 GeV, kappa (Codex "y") = 1e-3, e = 1e-7,
lambda_Y = lambda_Sigma = 1, lambda_X = 1e-3. Supplied preparation [Q, as in b5174e5 point.py and our
affleck_dine W1]: H_inf = 4.7e13 GeV, Hubble mass -c_H H^2 |X|^2 with c_H = 1 during inflation, switched off at
release into radiation domination (instant reheating). Phi is held at 0 (|X| >> X_c). Units below: H_inf = 1.
Inflationary minimum (own code, matches b5174e5 POINT_RESULTS to 9 digits): |X|, |Y|, |Sigma| = 23.82417, 1.891386,
0.2669760 H; V = -283.79980 H^4; five gauge-fixed Hessian eigenvalues positive.

## 2. (a) The real-subspace statement

**Proof [D].** With kappa, mu > 0 after rephasing (possible, section 3), V(z, z*) has real coefficients, so
conj(dV/dz*_i (z, z*)) = dV/dz*_i (z*, z). If z(t) solves z'' + 3H z' + dV/dz* = 0 with any real backgrounds
(a(t), H(t), c_H(t), xi R, a real tilt), z*(t) solves the same equation. Real initial data (z, z') therefore give
z(t) = z*(t) for all t by uniqueness (V is a polynomial, locally Lipschitz). Then n_i = -2 Im(z_i* z_i') = 0 for X, Y,
Sigma; the gauge current vanishes, so A_mu = 0 stays consistent with Gauss. The phase potential
2 kappa|Y||X|^3 cos A + 2 mu |Sigma||Y|^2 cos B (A = 3 theta_X - theta_Y, B = 2 theta_Y - theta_Sigma, the two
gauge-invariant angles) is minimized at A = B = pi for every set of amplitudes, which is a gauge rotation of the real
point (X > 0, Y < 0, Sigma < 0). Charges are gauge invariant, so own-potential alignment with zero velocities gives
exactly zero charge for any real radial history. Quantum mechanically the same C symmetry gives <Q> = 0 exactly
for a C-invariant initial state [D].

**Numerical check [C]** (`code/release.c`: conformal variables chi = a z, a = tau, own Dormand-Prince 5(4); not
Codex's DOP853 physical-time code):
- Exactly real data: all three comoving charges are identically 0.0 to tau = 1000, with no extra operator and with
  real-coefficient extra operators (Y*X^3|X|^2 at delta = 0, pi; Sigma*YX^3; Sigma*X^6).
- Codex's misaligned alternatives (X phase -/+ pi/24) at a = 14.177 (Codex t = 100): our a^3 n = (-16.261387466,
  6.395904449, -0.487720980) H^3, against Codex's (-16.261387466, 6.395904449, -0.487720980): agreement 1.2e-11
  (relative). Conjugate runs are exact mirrors. Gauss residual < 3.3e-11 relative, conformal energy ledger < 6.1e-13.

**New finding: the real subspace is invariant but unstable [C].** A CP-odd seed e0 in the X phase grows. In the
linear regime a^3 n_X / e0 is identical for e0 = 1e-15, 1e-12, 1e-9 and grows by about e^29 (4e12) between
tau = 5 and 130 (fitted envelope rate 0.20 per unit conformal time, faster after tau ~ 70). A 1e-15 rad seed gives
|eps_X| ~ 0.1 by tau = 300 (`results/real_subspace.json`). Inflationary phase fluctuations are of order
H/(2 pi sqrt(2)|X|) ~ 5e-3 rad or larger, so the zero-charge trajectory is not the generic outcome: ensembles of
256 releases with 1e-3 rad random phase seeds give eps_X with mean -0.002 +- 0.011, rms 0.177, mean |eps| 0.141,
max 0.56 (section 4). The mean stays zero (C symmetry); the local charge does not. Physically this is random-sign
charge in regions about the release Hubble size, i.e. X and anti-X matter side by side, not "no charge". The
homogeneous k = 0 instability implies that small-k modes grow too; a lattice calculation is needed beyond
tau ~ 50-100 [A].

## 3. (b) Operators and physical CP phases

A monomial X^a X*^a' Y^b Y*^b' Sigma^c Sigma*^c' has rephasing vector n = (a - a', b - b', c - c'); gauge
invariance is n_X + 3 n_Y + 6 n_Sigma = 0. The solutions form a rank-2 lattice with the **unimodular basis**
u = (3, -1, 0) (Y* X^3, kappa) and w = (0, 2, -1) (Sigma* Y^2, mu) [D]. Every operator is u^A w^B times neutral
singlets, and its coefficient c has the invariant phase arg c - A arg kappa - B arg mu. With N complex couplings,
physical phases = N - 2 (the gauge direction (1, 3, 6) acts trivially) [D].

Enumeration (`code/operators.py`, non-derivative operators; singlets phi^2, D^2, H^dag H, R under Candidate 1's
Z2_phi x Z2_D) [C]:

| Dimension <= | Complex couplings | Physical phases | New operators at this dimension |
|---:|---:|---:|---|
| 3 | 1 | 0 | Sigma* Y^2 |
| 4 | 2 | **0** | Y* X^3 |
| 5 | 10 | 8 | Sigma* Y X^3; Sigma* Y^2 x {abs(X)^2, abs(Y)^2, abs(Sigma)^2, phi^2, D^2, H^dag H, R} |
| 6 | 18 | 16 | Y* X^3 x {abs(X)^2, abs(Y)^2, abs(Sigma)^2, phi^2, D^2, H^dag H, R}; Sigma*^2 Y^4 |
| 7 | 55 | 53 | includes the direct Sigma* X^6 (class kappa^2 mu) and Sigma* Y^3 X*^3 (class kappa^-1 mu) |

If Z2_phi is broken, as the proposed wall bias does, phi Sigma* Y^2 is renormalizable and the counts are 1, 10,
26, 71 phases at dimension <= 4, 5, 6, 7. Any inflaton coupling chi Sigma* Y^2 is also renormalizable and complex.
R Sigma* Y^2 / Lambda and R Y* X^3 / Lambda^2 are covariant Hubble-induced A-terms of Dine-Randall-Thomas type,
the analogues of Codex's proposed xi_X R |X|^2.

Strength at the inflationary point for c = 1, Lambda = M_red (`phase_scan.json`, natural_ratio) [C]: Sigma* Y X^3 is
**13.0 times** the mu term; Sigma* Y^2 |X|^2 is 1.03 times the mu term; Y* X^3 |X|^2 is 2.1e-4 of the kappa term;
Sigma* X^6 is 3.5e-5 of the mu term. The two Planck-suppressed dimension-5 operators therefore set the inflationary
Higgs-phase angle B, which shifts by 1.49 rad (Sigma*YX^3) and 0.80 rad (Sigma*Y^2|X|^2) at delta = pi/2. A direct
Sigma* X^6 / Lambda^3 matters only for Lambda ~ 1e14 GeV, below the inflationary |X| = 1.1e15 GeV [C, A].

## 4. (c) Does a moving minimum supply the charge?

**Method [C].** Full homogeneous X, Y, Sigma model plus one extra operator C O (ops: Sigma*X^6, Sigma*YX^3,
Sigma*Y^2|X|^2, Y*X^3|X|^2), C = s|C|e^{i delta}, delta = the irremovable phase (delta = 0 aligned). The field
starts at rest at its own inflationary minimum including C O, so there is no assumed mismatch; the minimum
genuinely moves (A or B shifts listed above; it relaxes to the low-field alignment as |X| falls). Strength r is set
relative to the kappa term (Y*X^3|X|^2) or the mu term (Sigma-linear ops). Release as in section 1, to tau = a = 1000
(|X| ~ 1.1e12 GeV). Because the dynamics are chaotic (section 2), each configuration is an ensemble of 128 releases
with 1e-3 rad random phase seeds plus the unperturbed run. eps_X = a^3 n_X / (2 sqrt(2 lambda_X)
(a^4 rho_X / 3 lambda_X)^(3/4)), the tracks/affleck_dine definition; n_X is the late conserved quantity (n_Y
oscillates about zero by Y-Ybar mixing through mu Sigma* Y^2).

**Full-model results** (`results/phase_scan.json`; 67 CP-violating configurations, r from 3.5e-5 to 100,
delta = pi/4, pi/2, 3pi/4, -pi/2; 8,832 runs) [C]:
- Ensemble means are all consistent with zero: chi^2/dof of mean/sem = 0.94 over 67 configurations, largest
  |mean|/sem = 2.1 (one beyond 2 sigma; about 3 expected). |<eps_X>| <= 0.034, with sem 0.013-0.020 per
  configuration. No CP-driven net charge is resolved at the +-0.03 (2 sigma) level, including the natural
  dimension-5 cases that shift B by 1.5 rad.
- Local spread: rms 0.14-0.23, mean |eps| 0.12-0.19; largest |eps_X| in any run 0.81.
- The part odd in delta, (mean(pi/2) - mean(-pi/2))/2, is below 0.016 in all 17 (operator, r) pairs, with sem
  0.010-0.014: no CP-driven bias is resolved.
- **No run reaches eps >= 0.99175** (0 of 8,832).
- Unperturbed runs give |eps| up to 0.44, but are not tolerance-converged (chaos); only ensemble statistics are
  reported as results.
- Controls (`results/controls.json`): (i) tau_end = 1e4 versus 1e3 (N = 32): baseline mean 0.009 -> -0.019, rms
  0.141 -> 0.145; natural Sigma*YX^3 mean -0.027 -> 0.008, rms 0.165 -> 0.178. (ii) rtol 1e-12 versus 1e-10
  (N = 64, same seeds): mean -0.046 +- 0.019 versus -0.006 +- 0.022, rms 0.148 versus 0.175; individual runs
  differ (chaos), the distribution does not. Pooling the four ensembles of that configuration gives
  -0.013 +- 0.010. (iii) Seed amplitude 1e-6, 1e-3, 1e-2, 0.1 rad: rms 0.138, 0.175, 0.144, 0.148, means within
  +-0.027 (sem 0.017-0.022). (iv) Ledgers: Gauss residual <= 1.4e-4 relative (<= 4e-5 in controls), conformal
  energy ledger <= 5.5e-7 to tau = 1e3 and <= 5e-4 to 1e4.

**DRT comparator [C]** (an A-term present only during inflation, e.g. xi R Y*X^3/Lambda^2, absent after release):
at the curvature-coupled natural size (r = 4.5e-6 of the kappa term, xi = 1, Lambda = M_red) and at
r = 0.1 and 0.3 (initial A shifted by 0.10 and 0.30 rad), delta = +-pi/2, N = 64: means -0.045 +- 0.021,
0.021 +- 0.019, 0.016 +- 0.020, -0.006 +- 0.022, 0.002 +- 0.022, 0.015 +- 0.019; rms 0.15-0.18. Unperturbed
runs give +-0.364 (r = 0.1) and +-0.414 (r = 0.3), exact CP mirrors, but tiny seeds randomize them. A supplied
DRT-type mismatch, the input of our affleck_dine window, therefore also leaves no resolved net charge at this
point and release: the multi-field dynamics forget which side of the real subspace they started on.

**Reduced one-field envelope [C]** (`code/single_field.py`; one complex field, global U(1), V = lambda|X|^4 +
2 Re(a X^6) + 2 Re(b |X|^2 X^6), no Y, Sigma or Gauss constraint, no chaos): starting at the own minimum,
the moving minimum does generate a deterministic, tolerance-converged charge (rtol 1e-12 and a = 1e4 agree to 1e-6; exact CP mirror; zero without the second term, 8e-8). It is bounded by the
existence of the Hubble-held minimum, eta < 1/12 [D: R^2 - 1 - 3 eta R^4 = 0 needs 1 - 12 eta >= 0]. Max |eps| over
480 coarse points (353 with a minimum) and a 729-point refinement: 0.014, 0.031, 0.049, 0.124, 0.138, 0.172 at
eta = 0.005, 0.01, 0.02, 0.04, 0.06, 0.08; best 0.216 (eta = 0.0786, r0 = 0.42, delta = 2.01). Even the regular,
ungauged form of the mechanism stays below 0.22, a factor 4.6 short of 0.99175. (A first pass of this scan accepted
a bound-constrained non-stationary start and reported 0.38; it was discarded and the minimizer replaced by a grid
search for interior minima.)

**Comparison with the leftover budget [D].** Codex's requirement eps >= 0.99175 (1% allowance, our remainder law
f_extra ~ 0.9(eps^-4/3 - 1)). At the full-model local rms 0.17 the law gives f_extra ~ 8.6; at a net |eps| <= 0.03,
f_extra >= 95, unless hidden number-changing processes erase the symmetric part (not computed).

## 5. Verdicts

| Claim | Codex | Mine | Verdict |
|---|---|---|---|
| Real couplings + own alignment stay in the real subspace and give zero charge (b5174e5) | proof sketch + 3 releases | proof [D]; exact 0.0 [C] | CONFIRMED |
| Misaligned +-pi/24 releases, conjugate sign flip, Gauss (RELEASE_RESULTS.json) | -16.261387466, 6.395904449, -0.487720980 | same to 1.2e-11 with own integrator | CONFIRMED (independent recomputation) |
| Physical reading "own alignment supplies zero charge" | zero | mean zero, but the real subspace is unstable: seeds >= 1e-15 rad give local eps_X rms 0.18, random sign | CORRECTED |
| 7737487 interaction phases removable by rephasing; inflaton mismatch absent | stated | 0 physical phases at dim <= 4 (Z2_phi kept) | CONFIRMED |
| Next order | not counted | 8 phases at dim <= 5, 16 at <= 6, 53 at <= 7; 1 already at dim 4 if Z2_phi is broken | NEW (counted) |
| Hypothesis: a generic EFT supplies the charge without an assumed mismatch | review proposal | the minimum moves by up to 1.5 rad for natural dim-5 terms, but net <eps_X> is consistent with 0 (<= 0.03) and no run reaches 0.99175 | REFUTED at this point and release (the regular one-field envelope also stays <= 0.22) |

## 6. Limits

Classical, homogeneous (k = 0) fields; no lattice, so fragmentation and resonance at k != 0 are absent although
the k = 0 instability implies them. One supplied release (Hubble mass off, instant radiation era); W3-type
persisting Hubble mass, reheating and thermal masses not computed. Phi held at 0; the phi transition, walls,
heavy-particle decays and Y-decay CP asymmetries are not included. Ensemble seeds (1e-3 rad, Gaussian) stand in for
inflationary fluctuations; a light-angle isocurvature calculation is not done. eps_X counts the X condensate only;
symmetric X from Y decays would lower it further. tau_end = 1000 (charge changes of order 0.01 afterwards).
Same-session AI agent, not outside replication.

## 7. Credit and read levels

Affleck & Dine, NPB 249 (1985) 361 [R]; Dine, Randall & Thomas, PRL 75 (1995) 398 and NPB 458 (1996) 291
(Hubble-induced A-terms whose phases differ from the low-energy ones) [R]. Kasuya, Kawasaki & Takahashi, JCAP
0810:017 (2008), AD isocurvature [R]. Kusenko & Shaposhnikov, PLB 418 (1998) 46; Enqvist & McDonald, PLB 425
(1998) 309; Kasuya & Kawasaki, PRD 61 (2000) 041301 (condensate fragmentation) [R]. Kofman, Linde & Starobinsky,
PRL 73 (1994) 3195 (parametric resonance after inflation) [R]. Krauss & Wilczek, PRL 62 (1989) 1221 (discrete gauge
symmetry) [R]. Jarlskog-style rephasing-invariant counting is the standard method [R]. Hairer, Norsett & Wanner,
Solving ODEs I (Dormand-Prince pairs) [R].

## 8. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/gp_phase/code
python3 real_subspace.py    # builds release.c into the scratchpad (gcc), -> results/real_subspace.json (~5 s)
python3 operators.py        # -> results/operators.json
python3 phase_scan.py       # ~20 min on 4 cores -> results/phase_scan.json
python3 controls.py         # ~15 min -> results/controls.json
python3 single_field.py     # ~25 min -> results/single_field.json
python3 summary.py          # -> results/summary.json
```
Set GP_SCRATCH to choose where the release binary is built.
