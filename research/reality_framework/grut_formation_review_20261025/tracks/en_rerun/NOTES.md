# en_rerun: reproduction and code audit of Codex 7d9c93b (matched early nonlinear X/φ pilot)

27 October 2026. Track `en_rerun` of Claude's independent review (issue #13). Packet:
`research/reality_framework/grut_early_nonlinear_20261027/` at 7d9c93b, run from the git-archive extraction in
`$SP/e7d` (SP = /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad).

**Kinds of evidence.** Three kinds are kept apart throughout:

- **Reproduction (R-run):** Codex's own scripts and C++ engine, re-executed.
- **Independent (own code):** checks under `code/` that never import Codex modules. They read only Codex's
  outputs, initial states and snapshots.
- **Reading:** a line-by-line audit of the code.

Nothing here validates the model or excludes anything. All of it concerns the reduced scalar-limit X/φ pilot.

## Bottom line

- **All seven histories reproduce.** Every Codex script was rerun in REPRODUCE.md order. All seven nonlinear
  histories were also rerun (2 threads; five of them in τ segments with restarts). Every published number is
  reproduced at round-off: max relative difference ≤ 3.5×10⁻¹² for any field-history quantity. The preparation and
  matching files are bitwise identical or within 2×10⁻¹⁶.
- **No code errors found.** The code audit found no error in the potential, forces, stencil, integrator, clock,
  sampling normalisation or export conventions.
- **New qualifications.** Three things change how the headline numbers should be read.
  - The coherent mean follows a very eccentric orbit, with |χ̄|² going from 314 to 20 344. Several headline numbers
    are single-instant samples of orbit-modulated quantities. All seven controls share the same pericentre clock.
  - The **massless** antiparticle reference oscillates between 1.25 and 2.64 over τ 600–800 (2.238 at τ = 800
    exactly). The isotropic Hartree reference stays between 1.08 and 1.22.
  - The factor-1.9 gap between them comes entirely from k < 1.8. In that range 26 box-scale modes carry 1.07 of
    the 2.24.
  - The seeded zero-point half-quanta contribute 5.3×10⁻⁴ per net charge. That is negligible.

## Reproduction command

Environment file: `$SP/review/en_rerun/env.sh`. It exports the GRUT_* work directories and
OMP/OPENBLAS = 2. Run from `$SP/e7d`:

```bash
source $SP/review/en_rerun/env.sh; cd $SP/e7d
for s in dense_response correct_clocks match verify_engine matching_check; do python3 -B $PKT/code/$s.py; done
python3 -B tracks/en_rerun/code/run_one.py $PKT/code N48            # Codex run(), threads=2, unsegmented
for t in N64 halfstep48 seed2_48 box16_64 N32 cut10_48; do          # segmented restarts, same engine
  for e in 280 540 800; do python3 -B tracks/en_rerun/code/chunked_run.py $PKT/code $t $e; done
  python3 -B tracks/en_rerun/code/chunked_run.py $PKT/code $t finalize; done
for s in linear_phi_control analyze export_covariance summarize; do python3 -B $PKT/code/$s.py; done
python3 -I tracks/en_rerun/code/reproduction_compare.py $SP/review/en_rerun/published $SP/review/en_rerun/regenerated out.json
python3 -I tracks/en_rerun/code/audit_checks.py $PKT $GRUT_EARLY_WORKDIR out.json
python3 -I tracks/en_rerun/code/ir_tangent_scan.py $PKT out.json
python3 -I tracks/en_rerun/code/orbit_and_projection.py $GRUT_EARLY_WORKDIR out.json N48 N64 halfstep48 N32 seed2_48 cut10_48 box16_64
python3 -I tracks/en_rerun/code/projection_history.py $GRUT_EARLY_WORKDIR/N48 out.json
```

**Housekeeping.**
- Published outputs were copied to `$SP/review/en_rerun/published/` (MD5SUMS) before anything ran. Regenerated
  outputs are in `$SP/review/en_rerun/regenerated/`.
- Afterwards the packet directory in `$SP/e7d` was restored from `git archive 7d9c93b`, and `diff -rq` shows no
  differences.
- Run times: N48 177 s on 2 threads; N64 460 s.
- The machine has 2 cores, not the 4 the brief assumed.

**Artefacts of the segmented runs, not physics:**
- `conformal_light_crossings` in FIELD_RESULTS is computed from segment-relative s.
- The decimated halfstep curve has one extra sample, 15 602 instead of 15 601.

## Reproduction table (R-run vs published)

| Quantity | Codex | Rerun | Δ |
|---|---|---|---|
| DENSE (3,.5) 10% warning | 56.45–56.475 | 56.45–56.475 | bitwise |
| max amplitude through 80 | 20.174% | 20.17424% | bitwise |
| CORRECTED_CLOCKS Y_L number / injection | 5.68e-6 / 1.01e-8 | 5.6808e-6 / 1.0148e-8 | bitwise |
| MATCH a_m, C | 8.8692, 1.98580 | same | ≤2e-16 |
| scalar symplectic / φ det | 7.26e-11 / 5.3e-14 | same | bitwise |
| C++ vs FFT-DOP853 | 4.958e-13 | 4.958e-13 | bitwise |
| 4th-order refinement | 1.30e-5, 8.12e-7, 5.07e-8 (ratio 16.0, 16.0) | same | bitwise |
| seed across grids | 1.137e-13 | same | bitwise |
| background / phase tangent / neutral / force | 0.547% / 0.977–1.031% / 1.482–2.717% / 0.206% | same | bitwise |
| N48 contrast = 1, X variance ½ | 60.60, 69.25 | 60.60, 69.25 | exact |
| N48 outside mean / φ+mix / Q0/Q | 62.722% / 0.3098% / 1.05462 | same | ≤9e-15 |
| N48 energy / current drift | 1.254e-8 / 1.08e-14 | 1.254e-8 / 2.18e-14 | current: reduction order (2 vs 4 threads) |
| halfstep energy drift | 7.818e-10 | 7.818e-10 | 1e-5 rel |
| all-7 ranges: contrast, variance | 60.50–60.90, 69.20–69.25 | identical | exact |
| all-7 outside / φ+mix / Q0/Q | 60.94–64.93% / 0.277–0.365% / 1.0543–1.1585 | identical | ≤2e-14 |
| LINEAR_PHI_CONTROL at 800 | 2.0917 | 2.091713 | ≤3e-15 |
| reference n₋/Q (N48) | 2.2375 / 1.1800 | 2.237515 / 1.179962 | ≤2e-15 |
| EXPORT_CHECKS charge / energy | 4.4e-16 / 2.2e-16 | ≤4.4e-16 / ≤2.2e-16 | round-off |
| VERIFICATION groups | 7 pass | 7 pass (summarize.py asserts) | — |
| frozen-Σ estimate max | 0.447% | 0.4469% | exact |
| work ledger max | 5.07e-7 | 5.07e-7 | 9e-9 rel |

## Verdicts on Codex's claims (BRIEF_EARLY numbering)

| # | Claim (Codex) | Verdict | Ours | Method / label |
|---|---|---|---|---|
| 1a | Charge-contrast RMS reaches the mean at τ ≈ 60.5–60.9 | CONFIRMED | 60.50–60.90. On orbit average, the cycle τ 56.65–69.6 is the first with mean contrast ≥ 1 (1.11–1.16) in all seven runs | R-run [C]; own orbit average [C] |
| 1b | X fluctuation variance reaches half at τ ≈ 69.2 | CONFIRMED as computed; see N1 | 69.20–69.25. This is an instantaneous crossing 0.4 before pericentre (69.6), where \|χ̄\|² ≈ 344. Orbit-averaged fraction in that cycle is 2.6–3.0%. The cycle average never reaches ½ by τ 800 (max 0.466–0.496) | R-run [C]; own [C] |
| 2 | At τ 800, 61–65% of energy is outside the spatial mean; mean charge 1.05–1.16× net | CONFIRMED (instantaneous) | 60.94–64.93%; 1.054–1.159. In the last partial orbit these swing over 47–75% and 1.05–1.21. τ 700–800 means: 60.1–63.0%, 1.07–1.18 | R-run [C]; own [C] |
| 3 | Neutral + interaction 0.28–0.36% at 800; fixed-orbit control gives 2.09× | CONFIRMED | 0.277–0.365% (700–800 means 0.279–0.367%); control 2.0917 (max 2.109; first ≥ 1 at τ 702.5). Same convention at start: 8.032e-4 (control) vs 8.038e-4 (nonlinear). The curves agree within 4% to τ 100 and separate after about τ 150. Ratio 675 at τ 800 | R-run [C]; own reading of conventions [A] |
| 4 | Reference projections 2.238 (massless), 1.180 (isotropic Hartree); not physical counts | CONFIRMED at τ = 800; see N2 | Own full-FFT code gives 2.237515 / 1.179962. Over τ 600–800 the massless value spans 1.25–2.64 (mean 1.98); Hartree spans 1.08–1.22 (mean 1.16) | R-run + own [C] |
| 5 | Matching 0.55% background, ~1% charged, 1.5–2.7% neutral, ~0.45% frozen estimate | CONFIRMED (scope: see A4) | 0.547%, 0.977/1.031%, 1.482/2.717%, 0.206%, 0.447% | R-run [C] |
| 6 | Seven groups pass; energy 1.3e-8, current 1.9e-14, export at round-off | CONFIRMED | 1.293e-8; ≤ 4.3e-14 (thread-dependent round-off); 4.4e-16 / 2.2e-16. Halfstep ratio 16.04 gives 4th order | R-run [C] |

## Code audit (reading, plus own numerical checks)

### (a) `early3d.cpp`

All CONFIRMED.

- **Potential and forces.** X force −(g²χφ² + λ_X r²)χ_a. φ force −(λ_φ(χφ² − ν²) + g² r²)χφ. ν = 0 is written by
  `write_initial`; the file header reads [0, .01, .001, .03]. These are exactly ∂/∂χ of MATHEMATICS.md's 𝓗. [D]
- **FD4 stencil.** (4/3 Σ nn − 1/12 Σ nnn − 7.5 f)/dx² is the standard (−1, 16, −30, 16, −1)/12 per axis. [D]
- **Matched gradient energy.** The forward-difference energy ½(4/3|Δ₁f|² − 1/12|Δ₂f|²)/dx² has symbol
  4 sin²(k/2)(1 + sin²(k/2)/3)/dx², which is the same as −Laplacian and as `en_common.symbol`. [D]
- **Energy check.** My own real-space (np.roll) FD4 energy matches the C++ trace to 3.9×10⁻¹⁵ at τ 20 and
  3.8×10⁻¹⁵ at τ 800. [C]
- **Yoshida.** w₁ = 1/(2 − 2^{1/3}) = 1.35121, w₀ = 1 − 2w₁. The composition is S(w₁h)S(w₀h)S(w₁h) with the
  correct half-kick merging. [D] Fourth order is confirmed: energy drift ratio 16.04 between h and h/2. [C]
- **Role of the clock.** a(s) never enters the equations of motion: the conformal massless quartic system is
  autonomous for R = 0. a(s) only
  - maps s to τ for output labels and thresholds;
  - sets the output schedule, and therefore the late steps: Δs = 0.05/a < h_max = 0.004 for τ > 39.6, down to
    8.9×10⁻⁴ at τ 800, so h_max is inactive after about τ 40;
  - sets the a⁻² of the frozen-heavy estimate.

  τ(s) = 20 + a_m s + C s²/2 satisfies dτ/ds = a(s) = a_m + C s. [D] Against the full quiet (3, .5) background
  (dense_b3.npz), a_clock(τ) = √(a_m² + 2C(τ − 20)) deviates by at most 2.1×10⁻⁶ over τ 20–240, with R = 1.9×10⁻⁷
  at matching. [C]
- **Step logic.** The inner loop steps h = min(h_max, target − s) up to each output time and then sets s = target.
  It is symplectic for each step. The rounding residual dropped is ≤ 10⁻¹³. One cosmetic effect: `next_tau +=
  sample` accumulates rounding, so the last record is a near-duplicate (Δs = 3×10⁻¹²). Harmless. [C]
- **Current conservation.** Δ(χ₁π₂ − χ₂π₁) under a kick is h(χ₁∇²χ₂ − χ₂∇²χ₁). The radial potential force drops
  out, and this sums to zero because the stencil is symmetric. Under a drift it is h(π₁π₂ − π₂π₁) = 0. The current
  is therefore exact for each sub-step; the observed ≤ 4×10⁻¹⁴ is accumulated round-off. [D] The work identity
  dE_X/ds = ∫g²χφπφ|χ_X|² is correct. The ledger residual is 4.8×10⁻⁷ at Δτ 0.05 and 1.2×10⁻⁷ at 0.025, i.e. ∝ Δ²
  (trapezoid). [D/C]

### (b) `match.py` / `initial.py`

All CONFIRMED.

- **Initial covariance and transfer.**
  - S = [U/√(2w), 0; 0, U√(w/2)] gives C_qq = U diag(1/2w)Uᵀ, C_pp = U diag(w/2)Uᵀ, C_qp = 0, with
    w² = k² + eig(M₀). The minimum w² is 0.41 at k = 1.8, using λ_min = −2.832.
  - Then evolved = T(0→20)·S for the full 12×12 scalar system, and likewise for the 2×2 φ system. [D]
  - Own check: the full 12×12 τ-20 covariance stays pure, with symplectic eigenvalues ½ ± 3.6×10⁻¹¹. [C]
  - The light X marginal has symplectic eigenvalues 0.500002–0.520 (largest at k = 1.81). Heavy-field
    entanglement therefore adds at most 0.02 quanta per mode (mean 4×10⁻⁴). [C]
- **Sampling.**
  - Each pair ±n gets 14 complex standard normals z with E|z|² = 1. The coefficients are c = A z/√V with
    c₋ₙ = c̄ₙ.
  - The rows of A are X q (0, 1), X p (6, 7) and φ (q, p). The zero mode is set to the matched mean.
  - The seed key is (seed, int(100L), i+128, j+128, k+128), independent of N. [D]
- **Variance normalisation and the factor 2.** For one pair, E[c cᴴ] = AAᵀ/V and E[c cᵀ] = 0, so
  ⟨χ²⟩ = (1/V)Σ_n C(k_n), which is the continuum normalisation. [D] Own numerical test on the N48 τ-20 field:
  - Whitened coefficients over 7 384 independent pairs give covariance diagonal 0.983–1.009 (σ = 0.012).
  - Max off-diagonal deviation is 0.024; max pseudo-covariance is 0.022, consistent with zero.
  - Mean whitened norm per component is 1.0007.
  - Real-space X fluctuation variance is 11.88 against 11.44 specified (+3.9%, about 0.75σ for N_eff ≈ 370).
  - No coefficient exists outside the band: max 1.9×10⁻¹⁴.
  - **No factor-2 error.** [C]
- **Continuum k for preparation, FD4 symbol for dynamics.** This is stated and consistent. At the band top
  (k = 8), k_eff²/k² is 0.589 (N32), 0.875 (N48) and 0.954 (N64). So the prepared state is squeezed relative to the
  lattice dispersion near the top of the band. This is the main grid effect in φ+mix: 0.277 / 0.310 / 0.312%. [C]

### (c) Export and reference projections

All CONFIRMED; Hartree item qualified.

- **n± formula.** The formula is the standard real-component form summed over ±k. Both references reproduce the
  nonzero-mode charge exactly: (n₊ − n₋)/Q = −0.0546, which equals 1 − Q₀/Q. [D/C]
- **k0 removal and weights.** The rfft half-space weights are right, including the Nyquist plane. Bins are 50
  equal bins in physical lattice |k|. The k²-weighted covariance uses the FD4 symbol. [D]
- **Hartree mass.** ⟨|χ_X|²⟩ **includes the zero mode**: 52.5% of ⟨|χ_X|²⟩ (5 051 of 9 622) at τ 800, N48.
  - The isotropic factor is correct: 2λ⟨r²⟩ is the mean of the radial 3λ and tangential λ curvatures, and the trace
    average of the O(2) Hartree matrix including the condensate. [D]
  - The count is insensitive to these choices (N48): 1.180 with the zero mode, 1.183 without, 1.180 tangential,
    1.231 radial. All controls and variants span 1.174–1.315. [C]
- **Zero-point.** The count includes the seeded half-quanta with no subtraction. Their share is 14 768 modes × ½ /
  (n_Q V = 1.402×10⁷) = **5.27×10⁻⁴ per net charge**. That is 0.024% of 2.238 and 0.045% of 1.180 (1.03×10⁻³ for
  cut10). [C]
  - At τ 20 the full projections give 1.85×10⁻³ (massless) and 2.41×10⁻³ (Hartree). These include squeezing from
    the 0→20 transfer.
  - If all 48³ lattice modes had half-quanta it would be 3.9×10⁻³.
- **Occupation per mode at τ 800 (N48).**
  - About 10⁵ per mode for k < 1.8, 15 at 6–8 and 2 at 8–10.
  - 0.12 for k > 10. That is below the half-quantum, where the classical-statistical description fails
    (Epelbaum–Gelis–Wu, arXiv:1402.0115).
  - Those k > 10 modes carry ≤ 0.06% of n₋. [C]

## New findings

- **N1 (orbit phase). [C]** The matched coherent mean follows a highly eccentric reduced orbit: |χ̄|² runs 314 →
  20 344 → 314. All seven controls keep the same pericentres (τ 26.0, 34.9, 45.1, 56.65, 69.6, 84.0, …) through the
  breakup window.
  - Instantaneous ratios such as (X2 − |χ̄|²)/X2, "energy outside the mean" and Q₀/Q are therefore strongly
    orbit-modulated. Over the last partial orbit, outside-mean goes 47→73%; Q₀/Q goes 1.05→1.08 (N48) and
    1.15→1.21 (seed2).
  - The cross-control spreads at the single instant τ = 800 are smaller than these swings. Orbit-phase-averaged
    numbers (τ 700–800 means, listed above) should accompany them.
  - The τ 60.5–60.9 / 69.2 thresholds are pericentre-timed operational crossings. Codex already says this.
  - The robust statement is: *orbit-mean charge contrast first exceeds 1 in the cycle 56.65–69.6.*
- **N2 (the reference gap is IR- and phase-dominated). [C]**
  - At τ 800, k < 1.8 holds 1.874 of the massless 2.238 and 0.836 of the Hartree 1.180. For k ≥ 1.8 the two agree
    (0.364 vs 0.344).
  - The 26 lowest box modes (k < 1) alone carry 1.066 massless; 1.403 with 80 modes in the L = 16 box.
  - The massless value oscillates with the orbit: 1.25–2.64 over the ten snapshots τ 600–800, mean 1.98.
    Hartree: 1.08–1.22, mean 1.16. The "1.9-fold" gap at τ 800 is an instantaneous value; its instantaneous ratio
    ranges 1.04–2.45.
- **N3 (equal-time projection floor). [C]** Over all quadratic projections with a common ω(k) per FD4 shell, the
  minimum n₋/Q is 1.060 at τ 800 (N48), 0.93–1.09 over τ 600–800, and 1.056–1.142 across controls. Per mode it
  is 0.960. The Hartree reference lies within about 11% of this floor.
  - In the IR (k = 0.52, 1.05), the equal-time ratio ⟨|π|²⟩/⟨|χ|²⟩ − k² ≈ 23–25, close to the Hartree m² of 19.3.
    This argues against ω = k there; it is a single-snapshot diagnostic, not a pole.
  - None of these is a physical count.
- **N4 (lower seed edge untested). [C/E]** The band starts at k = 1.8. The τ = 0 tachyonic threshold is
  k = 1.683, so the L = 12 shell at k = 1.737 (24 modes, ω² = 0.18 > 0) is excluded by choice.
  - On the reduced homogeneous orbit, own tangent maps from τ 20 to 55 grow almost uniformly across k = 0.2–2.3:
    14–22×. k = 1.737 grows 22.2× against a maximum of 22.4× at k = 1.9.
  - The X variance at τ 20 is dominated by the band bottom: k ≈ 1.96 holds 11%.
  - So a lower-cutoff (IR seed) control is missing, alongside the upper cutoff-10 control. Its effect on the onset
    clock is expected to be logarithmic in seed amplitude [E]; it was not computed.
- **N5 (minor).** The early3d dynamics are pointwise in OpenMP, so the trajectory does not depend on the thread
  count. Only the diagnostic reductions do: the recorded current drift is 1.1–4.3×10⁻¹⁴ depending on threads.
  Restarting from 64-bit snapshots at τ 280/540 (or 300/560) reproduces every published observable at ≤ 10⁻¹²: no
  chaotic amplification of round-off is visible to τ 800 at this precision.

## Not checked or out of scope

- **Engine independence.** The independent lattice engine is another track (`en_indep`). Here only the C++ code
  was read, plus Codex's short FFT-DOP853 test (reproduced) and my FD4 energy re-evaluation.
- **Matching scope.** The full-versus-reduced matching covers only two momenta (2.13 and 7.1), one Σ
  preparation-phase tangent and e = 0 on both sides, through τ 50, before the strong-inhomogeneity window. It does
  not test gauge, heavy-field or nonlinear errors after τ 50. Codex says this.
- **Neutral catch-up comparison.** The comparison with our critic's neutral catch-up (τ 665–697) used a different
  band and a zero-point-subtracted occupation convention. Here the control crosses 1 at τ 702.5 with no
  subtraction. Not recomputed.
- **Physical results.** Physical n̄, f_sym, kinetic threshold and particle matching remain open, as Codex states.

Credits: classical-statistical lattice method (Khlebnikov–Tkachev; Micha–Tkachev); UV limits of the
classical-statistical approximation (Epelbaum–Gelis–Wu, arXiv:1402.0115); Yoshida's fourth-order composition (1990).

## Files

- `code/`
  - Independent: audit_checks.py, ir_tangent_scan.py, orbit_and_projection.py, projection_history.py,
    reproduction_compare.py.
  - Reproduction drivers: run_one.py, chunked_run.py, jdiff.py, ndiff.py.
- `results/`: audit_checks.json, ir_tangent_scan.json, orbit_and_projection.json, projection_history_N48.json,
  occupancy_N48_tau800.json, reproduction_compare.json, en_rerun_verdicts.json.
