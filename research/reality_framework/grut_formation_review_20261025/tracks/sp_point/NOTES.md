# Track sp_point: the supplied-phase point of Codex 373d37e (H_i = 1.1211e13 GeV)

25 October 2026 · Claude (issue #13) · independent review of Codex 373d37e,
`research/reality_framework/grut_inflation_scan_20261025/` (read through `git archive` only; no Codex code was run;
nothing committed). Own code in `code/`, outputs in `results/`. Same-session AI work, not outside replication.

**One supplied tree-level working model.** Nothing here is evidence that X, Y, Sigma or the inflationary inputs exist.
The point is a supplied, CP-odd initial condition; no phase-selection interaction exists in the action. Not a
validation, exclusion or prediction. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued.

Couplings unchanged [Q]: M_Y = 1e12, h = 1e11, mu = 5e11 GeV, kappa = 1e-3, lambda_X = 1e-3, lambda_Y = lambda_Sigma = 1,
e = 1e-7; c_H = 1; g* = g*s = 106.75; Y_Q = 4.8294e-9. Fields in H_i units.

## Method (independent)

- `code/sprel.c`: own C integrator of the conformal release equations chi'' = -dU/dchi* (own derivation from
  z'' + 3Hz' + V' = 0 with a = H_i eta). Two integrators: an 8th-order Yoshida composition of the drift-kick-drift map
  in extended phase space (step c/Omega, Omega a Gershgorin bound evaluated at each species' quartic turning amplitude;
  Gauss law preserved to rounding), with the **exact discrete tangent map**; and an adaptive Dormand-Prince 5(4).
  Codex used DOP853. Convergence: Yoshida error ratio 150-209 per step halving (8th order: 256); the two
  integrators agree to 1.4e-10 in the state at a = 100.
- `code/spmodel.py`: own global radial minimum (elimination + full root scan; multi-start 6D complex minimisation
  agrees to 1e-15 at 0.8-1.2 H*, 3e11 and 4.7e13 GeV). Amplitudes at H*: (23.826064, 1.893762, 0.430868) H [C].
  At 4.7e13 GeV this reproduces 1.1197e15 / 8.8895e13 / 1.2548e13 GeV (our earlier gp_matching values).
  *Internal note:* a first fsolve version landed on a non-minimal stationary point (|Sigma| about 40 H) at 0.95 H*; it
  was replaced before any reported scan, and the H* amplitudes were the same in both.
- Window averages are trapezoidal in a over 600-1000 (Codex's window); eps = mean of N_X/N_circ(E_X).

## 1. Reproduction of the handoff targets [C]

| Quantity (a = 600-1000) | Codex | Ours (Yoshida c = 0.00625 / DP5 1e-14) | Verdict |
|---|---|---|---|
| eps_X | -0.24066124 | -0.2406612341 / -0.2406612342 | CONFIRMED |
| N_X | -238.64868 | -238.6486708 | CONFIRMED |
| E_X | 741.64988 | 741.6498837 | CONFIRMED |
| J_Y | 13.26903 | 13.2690292 | CONFIRMED |
| Y_X / Y_Q | -65.9823 | -65.982276 | CONFIRMED |
| d eps for delta + 1e-9 | about 2.6e-7 | 2.5722e-7 / 2.5725e-7; exact tangent 257.21 per rad | CONFIRMED |
| tolerance effect | 6.7e-9 (their DOP853 2e-12 vs 2e-14) | spread over six integrator settings 1.9e-9; DP5 1e-12 vs 1e-14: 1.8e-9 | CONFIRMED (in kind; integrator-specific number) |
| mean shift a = 2000-3000 vs 600-1000 | about 0.00046 | 0.000460 (both integrators) | CONFIRMED |
| sigma = 0.03, N = 128, seed 2026102600 | -0.24909 ± 0.00333 | -0.249090 ± 0.003333 (identical draws, to 1e-16; member-by-member max abs diff 6.9e-5 against Codex's ENSEMBLE_MEMBERS.csv) | CONFIRMED |

## 2. Window extent in H_inf and in the supplied phase [C]

`results/window.json`. Each point is a release from Codex's aligned-radii preparation to a = 1000, window 600-1000.
Points on the lines also carry the exact tangent S = |d eps/d delta|.

*Line in H at delta = pi/2* (301 points, 0.8-1.25 H*, spacing 0.15%):
- |eps| ≥ 0.19 at 6% of points. Signs are mixed and ln S is 6-27 at most points.
- A negative-sign band runs from 0.985 to 1.0075 H* (2.2% wide): mean eps -0.192, |eps| ≥ 0.19 in half its points,
  ln S from 2.6 to 17.
- The contiguous run with eps ≤ -0.19 around H* is only 0.997-1.000 H* (3 points). The next point up, 1.0015 H*,
  already has -0.178.
- **H* is the only point within the ±0.15% grid spacing with ln S < 10.**

*Line in delta at H = H*:*
- 315 points over (0, pi): |eps| ≥ 0.19 at 8%. A coherent band runs over delta ≈ 1.46-1.64 rad (0.18 rad wide):
  mean -0.240, 89% of points ≥ 0.19, ln S 1.4-16. Outside it there are only isolated points of alternating sign.
- Fine scan at pi/2 ± 0.02 (2e-4 rad steps): mean -0.255, std 0.029, 98.5% ≥ 0.19. The autocorrelation falls to 0.4
  at 1e-3 rad and to 0 at 4e-3 rad. The map is jagged on a 1e-3 rad scale, which is consistent with S = 257.

*Codex's own robustness metric* (sigma = 0.03 rad, N = 32 per H, delta = pi/2): eps_net(H/H*) =

| H/H* | 0.98 | 0.99 | 1.00 | 1.01 | 1.02 | 1.03 | 1.05 |
|---|---|---|---|---|---|---|---|
| eps_net | +0.032 ± 0.015 | -0.122 ± 0.007 | -0.245 ± 0.007 | -0.123 ± 0.022 | -0.110 ± 0.018 | -0.059 ± 0.016 | -0.034 ± 0.018 |

All other H in 0.90-1.10 H* give |eps_net| ≤ 0.061. **Only H/H* = 1.00 clears 0.19.** The robust window has a full
width at half maximum of about 1.5-2% in H, and the passing core is narrower than ±1%.

At H*, the same metric across delta0 = k pi/12 (k = 1..11) gives eps_net = -0.252 ± 0.007 at k = 6 (pi/2) and
|eps_net| ≤ 0.031 at every other k, including 5pi/12 and 7pi/12.

*Focused 2D map* (with tangent): H/H* = 0.96-1.04 in 33 steps × delta = 1.30-1.85 in 23 steps, 759 releases.
- eps ≤ -0.19 in 11.2% of cells. The connected eps ≤ -0.19 region containing H* spans H/H* 0.9925-1.030 and
  delta 1.40-1.65 rad. Its area is 0.0034 (fractional H × rad), 7% of the box, and its mean eps is -0.237.
- The region is a curved tongue: a core of about 1.5% × 0.2 rad, plus a one-cell-wide streak that reaches +3% in H at
  delta ≈ 1.6. A separate branch sits near delta ≈ 1.30-1.35 at H/H* = 1.008-1.018.
- Median ln S is 8.2 inside the tongue and 16 outside it. The least sensitive cell is (H*, delta = 1.55), with
  ln S = -0.3, i.e. S ≈ 0.8.

**Answer to (1).** The productive region is a narrow, CP-signed tongue of lower sensitivity:
- about 1.5% in H by 0.2 rad in the supplied phase, with thin streaks;
- connected area 0.0034 in (fractional H) × rad;
- inside a sea where the median S is about 1e7 and |eps| ≥ 0.19 occurs only at scattered points of either sign.

In Codex's own sigma = 0.03 metric, only H/H* = 1.00 (within ±1%) and delta0 = pi/2 (within ±pi/12) pass. It is a
narrow resonance or island, not a window in H_inf.

## 3. Mechanism: why this point is regular when 4.7e13 GeV is chaotic [C, A]

Exact tangent along the supplied phase, S = |d eps/d delta| (window 600-1000), 41 H values 3e11-4.7e13 GeV
(`results/lyap.json`):

- S ≲ 0.06 below about 2.3e12 GeV: the map is smooth and linear there, but |eps| ≤ 0.12 at delta = pi/2 and pi/4.
- Above about 3e12 GeV, ln S rises roughly linearly: **ln S ≈ 2.14 (H/1e12 GeV) - 6.3** (delta = pi/2, residual rms 7.4
  e-folds; delta = pi/4: 2.61 H/1e12 - 11.5). S = 1 at about 3e12 GeV. Trend value at H*: 17.6 (S about 4e7); at
  4.7e13 GeV: 94 (computed 113).
- **H* is an island:** ln S = 5.55 (S = 257), 12 e-folds below the trend. Neighbours: 1.032e13 GeV ln S = 22.8,
  1.171e13 GeV 16.8; 1.05 H* (delta = pi/2) S = 6.8e6; 0.95 H* (delta = 1.0) S = 2.2e5.
- **When the charge forms** (`results/formation_epoch.json`). At H* the X charge is created at a ≈ 25-40 (N_X: -12 at
  a = 20, +45 at 25, -164 at 30, -252 at 40). Gauss's law makes Sigma the partner: 6 N_Sigma ≈ +216 at the end.
  - Across 3.3e12-4.7e13 GeV, the epoch a75 at which |N_X|/N_circ reaches 75% of its maximum tracks the kappa
    resonance a_res = 3 omega_X H/M_Y. That is the epoch when the Y conformal mass a M_Y/H reaches the third harmonic
    of the quartic X oscillation, which Y* X^3 drives.
  - a75/a_res = 1.0-1.4 from 3.3e12 to 1.1e13 GeV (1.03 at H*: 38.5 vs 37.4) and 0.70-0.87 from 1.5e13 to 4.7e13 GeV.
  - At 2e12 GeV the release is quick (a75 = 10, a_res = 6.9) [C].
- **How the sensitivity grows.** At H*, S rises from 1e-2 (a = 20) through 34 (a = 90) to 250 (a ≥ 500).
  - The tangent norm completes half its growth by a50 = 82, close to the epoch when the Sigma vacuum mass passes
    omega_X, a_Sigma = omega_X/(sqrt(2 lambda_Sigma) h/H) = 88 [E].
  - a50 scales roughly ∝ H: 34 at 3.3e12 and 282 at 4.7e13 GeV, against a_Sigma of 27 and 338.
- **Mechanism [A].** At fixed couplings, the kappa resonance epoch, the Sigma crossing and the duration of the
  nonlinear X-Y-Sigma stage all scale ∝ H/M_Y in conformal time, while the X oscillation frequency is fixed
  (omega_X ≈ 1.1). The number of strongly coupled X oscillations, and with it the accumulated log-sensitivity, therefore
  grows ∝ H, as the fit above shows. Charge is created by passage through the kappa resonance (a CP-odd supplied phase
  biases it). The outcome depends on the resonance phase accumulated since release.
- Diagnostic rescalings at fixed H (not the stated model): at H*, mu × 2 gives ln S = 56.7 (chaotic), mu × 0.5 gives
  3.6; M_Y × 0.71 gives 38, M_Y × 0.5 overflows; M_Y × 1.41 gives 1.4. At 4.7e13 GeV, M_Y × 2 lowers ln S from
  113 to 21.5. Chaos is driven by the mediator couplings: at fixed couplings every mediator-induced X self-coupling
  falls as (H/(a M))^2 in conformal units, so the nonlinear stage lasts a number of X oscillations ∝ H [A].
- **Conclusion:** there is no regime change between 1.1e13 and 4.7e13 GeV. The log-sensitivity grows ∝ H with ±7
  e-fold scatter, and the leading point is a low-sensitivity island inside a region whose typical S at that H is
  about 1e7. It is regular because it is a rare fluctuation of that trend, not because the dynamics is regular there.

## 4. Which preparation is physical? [C, E, A]

Inputs: `results/prep.json`, `results/ensembles.json` (parts B, C, D).

*The two preparations.*
- Supplied (Codex): aligned radii, Sigma rotated by delta = pi/2 at fixed modulus, so Sigma = -i |Sigma|. It costs
  Delta V = 2 mu |Sigma||Y|^2 = 0.1378 H^4 above the aligned minimum, which is 3.6 × (3H^4/8pi^2) [C].
- Conditional (Codex's control): radii re-minimised at fixed phase. At delta = pi/2 the mu term vanishes, so Sigma
  falls exactly to its vacuum modulus h/sqrt2 = 0.0063072 H, 68 times smaller [D].
  - Release from the exact conditional state: eps = -0.0584 [C].
  - With Codex's own conditional radii (Sigma = 0.0063103, off by 3.1e-6 because of the L-BFGS-B tolerance): eps = -0.0484.
  - Codex reports |eps| = 0.0516. The conditional release is therefore sensitive at the 1e-6 level and is not a
    sharp number: **|eps_cond| ≈ 0.05-0.06 (CORRECTED; qualitative conclusion unchanged)**.
  - It also carries a larger Y excitation (J_Y = 29.6): 0.79 pairs per net charge, 160 times the allowance.
- Conditional preparations across delta (61 phases in (0, pi)): rms eps 0.066 with 16 sign changes. One isolated
  phase (delta = 1.115) reaches -0.27.

*Neither preparation is an inflationary state of this action.*
- The Sigma phase is not a Goldstone direction. |Sigma| = 0.43 H ≫ h/sqrt2 is held by the linear source mu Y^2, and
  the masses at the aligned minimum are 0.61 H (transverse, the light "source phase") and 1.05 H (radial) [C].
  "Phase held, radius relaxed" therefore has no adiabatic justification: the phase relaxes as fast as the radii do.
- *Inflationary relaxation test.* Start from the supplied state inside de Sitter (Hubble mass on, full potential),
  evolve N e-folds, then release with the acquired velocities:
  - N = 0: eps = -0.241; N = 0.02: -0.202; N = 0.05: -0.167; N = 0.1: +0.012; N = 1: -0.127; N = 8: +0.090.
  - At N = 0.1, Sigma has moved by only 1e-3 H but has acquired a velocity of 0.02 H^2.
  - **The supplied state must be imposed within about 0.05 e-fold of the end of inflation**, in effect by an
    interaction acting at release.
- Partial radial relaxation, r = r_aligned + f (r_cond - r_aligned): eps = -0.241, -0.211, -0.093, -0.005 at
  f = 0, 0.05, 0.10, 0.15. The window survives only within about 5% of the aligned radii.

*Fluctuations.*
- The light source mode receives a de Sitter kick of H/(2 pi sqrt2 |Sigma|) = 0.26 rad per e-fold [D].
- In the Gaussian Starobinsky–Yokoyama equilibrium around the aligned minimum [E], the rms values are
  Im Sigma 0.23 H and Re Sigma 0.13 H.
- The supplied state lies 3.8 sigma (Mahalanobis) from the aligned minimum, with a density ratio of 7e-4.
- Wider supplied phase spreads, net eps (N = 128 each, SEM sampling only):

| sigma_delta [rad] | 0.03 | 0.1 | 0.26 | 0.5 |
|---|---|---|---|---|
| eps_net | -0.249 ± 0.003 | -0.179 ± 0.010 | -0.092 ± 0.011 | -0.064 ± 0.010 |
| pairs per net charge | 0.098 | 0.19 | 0.42 | 0.59 |

- Gaussian de Sitter draws of all six real field components [E]. Chaotic members depend on the integrator; ensemble
  means do not: a rerun of B and C with the final step control changed eps_net by ≤ 0.003.
  - around the aligned minimum: eps_net = -0.001 ± 0.005 (N = 256), local rms 0.081, |eps| ≥ 0.19 in 1.2% of draws.
    This is CP-symmetric, as required.
  - centred on the supplied Sigma: -0.003 ± 0.006.
  - centred on the supplied Sigma, Sigma fluctuating only: -0.010 ± 0.009.
  - centred on the supplied Sigma, X and Y fluctuating only: -0.064 ± 0.009.

*Answer.*
- Neither preparation follows from the baseline action. The inflationary state is a CP-symmetric fluctuating state
  around the aligned minimum, and it gives zero net charge.
- The supplied aligned-radii rotation is the most favourable choice and the least physical one. The conditional
  minimum is not an adiabatic attractor either, since there is no mass hierarchy between phase and radius.
- A preparation that keeps |eps| ≥ 0.19 must:
  1. act within about 0.05 e-fold of release;
  2. hold the source-phase spread to ≲ 0.08 rad, which means making the transverse Sigma mode heavy (≳ 1.5 H, against
     0.61 H now) so that its de Sitter fluctuations are suppressed;
  3. also suppress the X/Y fluctuations (X/Y-only draws alone give -0.064; the X-radial mode, m^2 = 1.99 H^2, is at the
     edge of the light-field regime where the Gaussian estimate holds);
  4. leave the Y/Sigma radii within about 5% of the aligned values.

  This is a demanding target for the "phase-preparation interaction" that Codex lists as next [A].

## 5. Longer runs [C]

`results/longrun_y8.json` (Yoshida c = 0.025, a = 1e4), `results/longrun.json` (DP5 1e-12, a = 3e4):

| window | 600-1000 | 1000-2000 | 2000-3000 | 3000-5000 | 5000-10000 |
|---|---|---|---|---|---|
| eps_X | -0.24066123 | -0.24027393 | -0.24020157 | -0.24020170 | -0.24020169 |
| N_X | -238.64867 | -238.26981 | -238.19933 | -238.19937 | -238.19937 |
| N_eff = N_X + 3 N_Y | -232.52 | -237.19 | -238.00 | -238.15 | -238.19 |
| J_Y (Codex) | 13.26903 | 13.25950 | 13.25866 | 13.25848 | 13.25845 |
| J_Y,eig | 13.691547 | 13.691544 | 13.691544 | 13.691544 | 13.691544 |

- **The final homogeneous charge fraction is eps_X = -0.240202**: it is frozen to 2e-7 from a = 3000 to 1e4, and its
  magnitude is 0.19% below the 600-1000 value.
- The Y charge mean decays, so N_eff → N_X.
- The eigenmode-resolved Y number (each Majorana-split real mode divided by its own frequency) is the adiabatic
  invariant: it is constant to 1e-7. Codex's single-frequency J_Y drifts by -0.08% and is 3.1% low.
- DP5 1e-12 run to a = 3e4 (2.8e9 steps): eps_X = -0.24020169 in every window from 3000-5000 to 20000-30000, and
  N_X = -238.199372. It agrees with the Yoshida run to 1e-8. In this run J_Y,eig drifts by -1e-5 over 2e9 steps; the
  Yoshida run shows no such drift.
- Scope: this is the k = 0 mode only. Fragmentation (Q-ball formation) is expected long before a = 1e4 and is not
  modelled.

## 6. Y-excitation ledger [C, D, E]

`results/ledger.json`, built on our own integration:

- Gamma(Y → 3X) = 3 kappa^2 M_Y/(256 pi^3) = 377.95 GeV [D, matches Codex].
- Decay epoch: a_dec = sqrt(H_i/Gamma) = 1.7223e5. With t = 1/Gamma instead, a_dec = 2.44e5 [D].
- Y/X energy at decay: 274.85 (CONFIRMED). It is 283.6 with J_Y,eig and 388.7 with the t = 1/Gamma convention [C].
- Pairs per net charge: 1.5 J_Y/|N_X| = 0.083401, which is 16.92 times the quoted 0.0049296 allowance (CONFIRMED).
  - With J_Y,eig: 0.08606 (17.5 times).
  - With J_Y,eig and the final charge N_eff: 0.08832 (17.9 times).
  - Taking the ensemble ratio (sigma = 0.03) instead: 0.0977 (19.8 times).
- Energy-equivalent quartic fraction: 0.003555 (CONFIRMED arithmetic; it is a diagnostic only).
- Checks of the conserved-number premise that Codex did not state [E]:
  - The condensate-assisted two-body channel Y → X X through 3 kappa X_bg has Gamma_2/H ≈ 1.7e-3, constant in
    radiation. Integrated from a = 10 to decay, it is 1.7%.
  - The parametric resonance parameter of the Y condensate on X fluctuations is q ≈ 2e-5 at a = 1000, so there is no
    broad resonance.
  - So in the homogeneous model the Y number survives to perturbative decay.
  - As dark radiation the Y energy is only 1.4e-6 of the visible energy at decay (rho_X/rho_vis = 5.2e-9), so ΔN_eff
    is negligible.
- Thermalised branch: if the Y daughters equilibrate with X, the hidden energy per net charge rises 276-fold, and the
  energy-equivalent fraction 0.0036 lies far inside the failing regime of the gp_leftovers equilibrium benchmark
  (needs ≳ 0.19) [A]. This is not a solved mixed-plasma relic.

## 7. Verdicts

| # | Codex claim (373d37e) | Ours | Verdict |
|---|---|---|---|
| 1 | Point at H* = 1.1211167e13 GeV, delta = pi/2: eps -0.24066124, N_X -238.64868, E_X 741.64988, J_Y 13.26903, Y_X -65.9823 Y_Q | Same to all quoted digits with two independent integrators [C] | CONFIRMED |
| 2 | 1e-9 phase shift changes eps by 2.6e-7; tolerance 6.7e-9 | 2.572e-7; exact tangent 257.2 per rad; integrator spread 1.9e-9 [C] | CONFIRMED |
| 3 | a = 2000-3000 mean moves 0.00046; "late fraction remains about 0.241" | 0.000460. The final value is -0.240202, frozen to 2e-7 from a = 3000 to 3e4 (two integrators) [C] | CONFIRMED (final value -0.24020) |
| 4 | sigma = 0.03, N = 128, seed 2026102600: -0.24909 ± 0.00333 | -0.249090 ± 0.003333, identical draws [C] | CONFIRMED |
| 5 | "A small-noise ensemble preserves it" | True only for sigma ≲ 0.08 rad: sigma = 0.1 gives -0.179 ± 0.010, 0.26 gives -0.092, 0.5 gives -0.064. Gaussian de Sitter draws centred on the supplied Sigma give -0.003 ± 0.006 [C, E] | CORRECTED (scope) |
| 6 | Local map: structured and narrow (±2% in H, ±0.06 rad) | A tongue of about 1.5% in H by 0.2 rad (connected area 0.0034); with sigma = 0.03 ensembles only H/H* = 1.00 passes (0.99: -0.122, 1.01: -0.123), and only delta0 = pi/2 among k pi/12 [C] | CONFIRMED, narrower than "±2%" suggests |
| 7 | The point is "regular" | Locally (S = 257), but as an island 12 e-folds below the trend ln S ≈ 2.1 H/1e12 - 6.3; S = 2e5-7e6 at ±5% [C] | CORRECTED (island, not regime) |
| 8 | Conditional radial re-minimisation gives |eps| ≈ 0.0516 | -0.0584 from the exact conditional state (Sigma = h/sqrt2); -0.0484 from Codex's radii (Sigma off by 3e-6). Sensitive, ≈ 0.05-0.06 [C] | CORRECTED (value not sharp; conclusion holds) |
| 9 | J_Y subtraction: 0.0834 pairs per net charge, 16.9 times the 0.00493 allowance | Same. Eigenmode-resolved J_Y = 13.6915 (+3.2%, and the true adiabatic invariant): 0.0861 (17.5 times); with N_eff, 0.0883 (17.9 times); sigma = 0.03 ensemble, 0.098 (19.8 times) [C] | CONFIRMED, small upward correction |
| 10 | Gamma_Y = 377.95 GeV, a_dec = 1.72e5, E_Y about 275 E_X, energy-equivalent fraction 0.00356 | Same [D, C]. The t = 1/Gamma convention gives 389 E_X. Two-body and parametric channels do not change the number (1.7%, q ~ 2e-5) [E] | CONFIRMED |
| 11 | "Counterexample to a universal X-only no-window claim for supplied mismatches" | Holds as stated: a supplied, instantaneous, CP-odd, aligned-radii preparation [C] | CONFIRMED (conditional) |
| 12 | The physical preparation is unresolved | Neither preparation is an inflationary state. The preparation must act within about 0.05 e-fold of release, hold the phase spread ≲ 0.08 rad, and keep the radii within 5% of the aligned values. Natural de Sitter fluctuations erase the net charge [C, E, A] | CONFIRMED and sharpened |
| - | Phase-preparation interaction; thermal/mixed-plasma relic; fragmentation; vector population | not computed | NOT CHECKED |

## 8. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/sp_point/code
python3 repro.py ../results/repro.json                 # ~1 min on 4 cores
python3 ensembles.py ../results/ensembles.json A B C D # ~15 min
python3 prep.py ../results/prep.json                   # ~3 min
python3 lyap.py ../results/lyap.json                   # ~10 min
python3 window.py ../results/window.json lines ens map # ~60 min
python3 longrun.py ../results/longrun_y8.json Y8_c0.025_a1e4 8 0.025 10000   # ~3 min
python3 longrun.py ../results/longrun.json DP5_1e-12_a3e4 5 1e-12 30000      # ~40 min
python3 formation_epoch.py ../results/formation_epoch.json  # ~2 min
python3 ledger.py ../results/repro.json ../results/ensembles.json ../results/ledger.json
```
The C file compiles automatically (gcc -O3, -lm) into the scratchpad `bin/`. Python 3, NumPy, SciPy.

Credit: Affleck–Dine; Dine–Randall–Thomas; Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki (Q-ball formation and
fragmentation); Starobinsky–Yokoyama (de Sitter equilibrium); Yoshida (1990) composition; Dormand–Prince; Benettin et
al. (tangent-map Lyapunov method) [R].
