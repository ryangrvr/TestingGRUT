# pt_response SKEPTIC: adversarial check of track pt_response (Codex 5f5513b claim 4, issue #13)

1 November 2026. Claude, key `pt_response_skeptic`. Briefs: [BRIEF_TRANSPORT.md](../../BRIEF_TRANSPORT.md) and
[BRIEF.md](../../BRIEF.md). This checks [NOTES.md](NOTES.md) of track `pt_response`.

This is same-session AI work on one supplied, reduced, classical X/φ scalar model. It is not outside replication.
Nothing here is a validation, an exclusion, a GRUT prediction or new theory.

Claim labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## Provenance of the numbers

- **Codex engine (reproduction-type use).** Only for background states and tangents.
  - `skeptic/code/make_states.py` advances the regenerated base48 field to s = 20.12, 20.56 and 20.92.
  - `skeptic/code/probes_phase.py` runs Codex's exact tangent of the discrete Yoshida map, from our copy of the tree.
  - All six background states match the published `FIELD_TRACES.npz` mean field to ≤ 3×10⁻¹⁰ relative.
- **Independent (our code and our method).** Everything else:
  - the probe designs: family A, M = 8 random-phase probes on shell 1 only; family B, M = 8 on shells 2–4;
  - every estimator and channel algebra;
  - the closure models and the Goldstone-residual analysis;
  - the canonical-sum and free-map checks;
  - the J_r identity and the toy model.
- **An earlier run of this same key.** It was interrupted on 27 October before it wrote SKEPTIC.md. It left:
  - the code `analyse_phase_probes.py`, `circular_bare.py`, `counting_checks.py`, `toy_exchange.py`,
    `probes_phase.py` and `make_states.py`;
  - the projections at s0 = 20.12, 20.24 and 20.68.

  I reviewed that code and re-ran all of its analyses today. `counting_checks.json` regenerates identically. The s0 =
  20.92 probes are new today.
- **This run (resumed on 1 November after a second interruption).** Every number in this file was re-checked against
  `skeptic/results/*.json`. Two statements had been written without a stored computation; both were recomputed with
  new independent code and are now backed by `transform_nyquist.json` (canonical sum on the Nyquist grid, with
  Gaussian- and Lorentzian-damped test lines) and `initial_slope.json` (the O(Δ³) slope identity). Four wordings were
  tightened (R4, R6, N2 and the closure medians); no verdict changed.

Results are in `skeptic/results/*.json`. The large projection files (17–34 MB) are kept in
`$SP/review/pt_response_skeptic/probes/`.

## Verdicts on pt_response

| ID | pt_response | Skeptic recomputation | Verdict |
|---|---|---|---|
| R1 IR ratio | CORRECTED to 0.46–1.31 (two phases) | Six phases over one full radial period, shell 1: 1.309 (exact) / 0.522 / 0.283 / 0.501 / 0.923 / 1.258 at s0 = 20 / 20.12 / 20.24 / 20.44 / 20.68 / 20.92. Shells 1–3 span **0.28–1.31**. At 20.24 the bootstrap 16–84% is 0.26–0.31. The ratio is orbit-locked: 20.92 ≈ 20 + T_r gives 1.26 against 1.31. Our own reconstruction of the exact shell-1 data reproduces 1.309, the white-noise scatter 1.08–1.60 and the scattered/forward amplitudes 0.58/1.07/1.32/1.36 | **REFUTED** (the range). The per-phase values are confirmed |
| R2 canonical sum | identity; tests only the transform | Confirmed as an identity. The 0.1–0.4% residual is spectral weight beyond \|ω\| = 20, not transform error: on the Nyquist grid (±78.5) all 45 published shell sums give ≤ 3.4×10⁻⁵ and our low-noise traces ≤ 5.1×10⁻⁵. On Codex's grid, undamped test lines give ≤ 1.9×10⁻⁵ and Gaussian-damped lines ≤ 2×10⁻⁶, while Lorentzian-damped lines give 1.0–3.6%; the published residuals correspond to effective widths of 0.01–0.12 | CONFIRMED, refined |
| R3 free control | 2.161928e-10 | Exact discrete Yoshida map in 50-digit arithmetic: 2.1619489e-10, ω = 1.2541202921170822 | CONFIRMED |
| R4 definitions | stated; slopes 0.986–0.997 "consistent with ω² = k² + Hartree m²" | Our algebra matches in every script. The 20 low-noise initial slopes equal 1 − Δ²(κ̄² + ½tr H̄_X(s0))/6 − (½tr H̄_X)′Δ³/12 to ≤ 1.4×10⁻⁴ (≤ 3.3×10⁻⁴ without the drift term); Codex's single probes deviate by up to 1.3×10⁻³ in shell 1. Here ½tr H̄_X is exactly Codex's Hartree m², and it swings 15.7–25.3 with orbit phase. The slope is an O(Δ³) identity of the instantaneous mean Hessian, not a dispersion test | CONFIRMED; reading corrected |
| R5 matrix response | UV Floquet–HFB; IR fails; bare BdG unstable | Main picture holds; three specifics corrected (below) | CONFIRMED with corrections |
| R6 circular Bogoliubov | ratio 0.71–0.80 … 0.055–0.08; 97–99.8% antiparticle | A lab-frame DOP853 integration reproduces 0.705/0.804/0.800 and 0.054–0.080, with \|v/u\| 0.0188–0.0499. The radial-branch antiparticle fraction is **99.35–99.93%** (\|ω\|-weighted action) and 99.75–99.96% (amplitude²), not 97%; with \|v/u\| ≤ 0.05 it cannot be below 99.75% | CONFIRMED; minor correction |
| R7 J_r vs lab | +0.8–2.7% (0.6237/0.4798/1.6815) | Starting orbit, three routes: J_r = 1.7637252794 and lab = 1.7377883531 (Codex ✓). In the medium pt_response's weighting is sin², which is not like-for-like (below) | CONFIRMED in direction; numbers CORRECTED |
| N1 quasiparticle counting | 1.770–1.775; Floquet basis undefined | With like-for-like J_r: 1.7625–1.7668 (+0.6–0.9%). The toy shows the X↔X̄ conversion is an exchange with the condensate, so the total is conserved | CORRECTED (number); CONFIRMED (sector-wise) |
| N2 X̄ = (S − Q)/2 | S/Q 4.483–4.504; identity to 2e-4 | S/Q 4.4828–4.5035 ✓. The identity A₊ − A₋ = Q holds to **≤ 1.02×10⁻³** over the 18 windows: 1.2×10⁻⁴ for base48 s 14–24 (pt_response's case), up to 5.0×10⁻⁴ in the other 14–24 windows and 5.4×10⁻⁴–1.02×10⁻³ in every 24–34 window | CONFIRMED; tolerance CORRECTED |

## 1. The IR ratio is a phase- and window-dependent number [C]

The table gives the conjugate/normal RMS ratio over delays 0–6 ([meas_traces.json](skeptic/results/meas_traces.json),
[phase_probes_all.json](skeptic/results/phase_probes_all.json), [phase_probes_2092.json](skeptic/results/phase_probes_2092.json)).

| s0 | r (orbit) | Shell 1 | Shell 2 | Shell 3 | Shell 4 |
|---|---|---|---|---|---|
| 20 | 96.6 | 1.309 (exact) | 0.884 | 0.747 | 0.619 |
| 20.12 | 100.7 (apocentre) | 0.522 ± 0.03 | 0.566 | 0.464 | 0.407 |
| 20.24 | 87.4 | 0.283 ± 0.02 | 0.396 | 0.454 | 0.382 |
| 20.44 | 41.4 | 0.501 | 0.456 | 0.588 | 0.484 |
| 20.68 | 51.6 | 0.923 | 0.707 | 0.731 | 0.632 |
| 20.92 | 98.3 | 1.258 ± 0.06 | — | — | — |

**The ratio depends on the delay window.**
- Over delays 0–2 it is 0.16–0.58 across all phases and shells 1–4; over delays 2–6 it is 0.30–1.72.
- In the exact shell-1 data at s0 = 20 it rises from 0.25 (t < 1) through 0.69, 1.46 and 2.02 (t = 4–6).
- The conjugate channel grows relative to the normal one with delay. This is what a static X↔X̄ mixing would do. A
  splitting δω ≈ 0.4 gives an envelope ratio of tan(δω t/2), about 1.6 at t = 5.
- So both Codex's 0.44–1.70 and pt_response's 0.46–1.31 are properties of the chosen 6-unit window and the kick phase.
  They are not constants of the medium.

## 2. The closure conclusions: what survives [C/D]

Results are in [closures_phase.json](skeptic/results/closures_phase.json), [closures_uv.json](skeptic/results/closures_uv.json)
and [circular_bare.json](skeptic/results/circular_bare.json).

**Goldstone residuals on the recorded orbit.** We apply each closure's Hessian to u = Jq̄, the global-phase direction
of the coherent mode, and project ρ = u″ + H u on the tangent.
- **bare: −6.00** (−3.3 at apocentre, −9.7 at pericentre). This is tachyonic: the measured orbit is held by about 6
  units of extra central force from the fluctuations, which the bare Hessian omits.
- **meanH: +5.03** (+3.96 at apocentre, +0.09 at pericentre), with a radial residual of 2.35 RMS.
- So the Gaussian/HFB mean Hessian is **gapped**. This is the known conserving-versus-gapless problem of HFB
  (Hohenberg–Martin 1965; Griffin, PRB 53, 9341 (1996)) [R]. The anomalous pair density shifts the phase mode.

**(a) The bare-BdG instability is mostly this artefact.**
- A circular estimate with tangential mass −Σ = −6.0 gives a rate of 0.78 at k = 0, unstable for k² < 6.0. This
  matches pt_response's "k² ≤ 7, rates 0.19–0.82". A direct Lyapunov integration on the measured orbit gives
  0.86/0.76/0.54/0.12/0.19.
- A self-consistent eccentric orbit with the same r_max/r_min has a genuine parametric band at k² ≲ 4.0 with maximum
  rate 0.38. With a Hartree-dressed orbit matching the measured period of 0.90, the maximum rate is 0.30.
- So "condensate-only BdG is parametrically unstable for k² ≤ 7 at 0.19–0.82" overstates the band and the rate. The
  real parametric resonance is narrower and about 2× slower.

**(b) The IR pole offset is the HFB gap.**
- Over 5 phases × shells 1–4, the meanH ρ-peaks are a median 11.4% too high (maximum 24%, mean signed +12%). This is
  pt_response's "measured poles 8–20% low".
- A Goldstone-repaired closure (meanH minus the tangential defect) brings them to a median 3.5% (maximum 12%). The
  rank-2 repair gives a median 5.6% with mean signed +0.7%.
- The normal-channel L2 error over 0–6 drops from 0.93–1.41 to 0.52–1.22 over shells 1–4 (median 1.22 → 0.79; over
  shells 1–3 the median goes 1.24 → 0.84).
- **But no tested homogeneous closure fits everything.**
  - The repaired closures overpredict the UV ratio by ×1.9–2.7, and their conjugate-channel L2 errors stay at 0.9–2.7.
  - The measured momentum scattering (≥ 1.07× forward after t > 1) is model-independent evidence of inhomogeneity.
  - So pt_response's IR conclusion stands, but its "poles 8–20% low" is evidence of an HFB gap, not of non-mean-field
    physics.
- **meanH captures the orbital-phase pattern of the IR ratio.** Its correlation with the measurement is 0.87 over 15
  shell-phase points (model 0.30–0.70, measured 0.28–1.31). It matches at s0 = 20.12 and 20.24 and is about 2× low
  near apocentre. pt_response did not report this.

**(c) The UV poles agree less well than stated.**
- pt_response's own 2-pole fits ([pole_fits.json](results/pole_fits.json)) give a meanH error of **5.5%** (shell 5,
  s0 = 20.44) and 4.4% (shell 5, s0 = 20, seed 7).
- So the bound is ≤ 2.4% only for k ≥ 3.1, and ≤ 0.6% for k ≥ 4.7 ✓. Our ρ-peak method gives ≤ 5.8% and ≤ 1.25%.
- The UV ratio errors, −27% to +21% with a median |error| of 9.4%, are confirmed.

## 3. Counting [D/C]

**Identity [D].** For a quasi-periodic orbit with ω_n = ω_θ + nΩ_r:

  J_r = Σ n ω_n |c_n|² = (⟨|ψ′|²⟩ − ω_θ q)/Ω_r.

So, when ω₀ > 0 > ω₋₁,

  J_r − A₋ = Σ_{n≥1} n ω_n |c_n|² + Σ_{n≤−2} |n+1| |ω_n| |c_n|² ≥ 0.

The radial action therefore always exceeds the lab statistic, by the weight of the higher harmonics. For the starting
orbit this is 0.02594 = 0.01653 + 0.00940 per charge, which reproduces Codex's 1.763725 − 1.737788
([jr_identity.json](skeptic/results/jr_identity.json)).

**Weighting in the medium.**
- The lab statistic is an h²-weighted (Hann²) time average of line power: by Parseval, Σ_ω P(ω) = Σ_j h_j²|ψ_j|²/Σh_j²,
  so a slowly drifting line amplitude enters with weight h² = sin⁴. pt_response weighted the per-period J_r with sin²
  (Hann), not sin⁴. (In `counting_checks.json` the key `Jr_hann2_per_Q` is, despite its name, the sin²-weighted
  value, 0.6243/0.4799/1.6745, and reproduces pt_response's choice; the like-for-like sin⁴ values are the `hann2`
  entries of `jr_identity.json`.)

| Window (s) | Lab | J_r, Hann² (like-for-like) | J_r, Hann (pt_response) | J_r, uniform | Per-period J_r drift in window |
|---|---|---|---|---|---|
| 0–6 | 1.6657 | 1.6911 (+1.5%) | 1.6745 (+0.5%; pt_response 1.6815) | 1.5502 (−6.9%) | 1.761 → 1.043 (−41%) |
| 14–24 | 0.6075 | 0.6160 (+1.4%) | 0.6243 (+2.8%) | 0.6424 (+5.8%) | 0.782 → 0.588 (−25%) |
| 24–34 | 0.4760 | 0.4783 (+0.5%) | 0.4799 (+0.8%) | 0.4872 (+2.4%) | 0.566 → 0.407 (−28%) |

- J_r is not conserved in the medium: it drains by 25–41% within each window. The 0.5–1.5% J_r/lab distinction is
  smaller than the weighting systematic.
- The radial-action/lab distinction is exact for the isolated orbit, and only indicative in the medium.
- With the like-for-like coherent J_r, the circular/adiabatic total is **1.7625–1.7668** (+0.6–0.9% over 1.7517),
  not 1.770–1.775.

**The X↔X̄ conversion is an exchange [C].** This uses a closed 3-mode Galerkin truncation of the quartic field with an
eccentric condensate (r_max/r_min ≈ 3.9) and waves started as pure X̄ ([toy_exchange.json](skeptic/results/toy_exchange.json)).
- At k² = 20, the waves' A₋ swings 20.7–38.6 and the wave charge runs from −31.4 to +4.3.
- The waves' A₋ plus the condensate's change in A₋ stays at 48.0–50.1 (k² = 20) and 65.12–65.32 (k² = 45).
- The process is X̄_w + X_c ↔ X_w + X̄_c. So "X̄ is not a quasiparticle quantum number" holds for the wave sector alone.
  The total coherent-plus-waves budget (S − Q)/2 is not lost, which is consistent with pt_response's "only Q and S carry
  over".

**The spectral charge identity.** Its tolerance is 1.02×10⁻³ (IR12_s1_w3, s 24–34), not 2×10⁻⁴
([counting_checks.json](skeptic/results/counting_checks.json)). The 2×10⁻⁴ holds only for the base48 and base48_half
14–24 windows (1.2×10⁻⁴); eleven of the 18 windows exceed it. The residual is the Hann-window leakage of the k ≤ 8
band and grows with the wave share, so it is a property of the finite transform, as for the canonical sum.

## 3b. Two things that are identities, not physics [D/C]

Results are in [transform_nyquist.json](skeptic/results/transform_nyquist.json) and [initial_slope.json](skeptic/results/initial_slope.json).

**The canonical sum.** On the full Nyquist grid of the Δ = 0.04 sampling (±78.5), Σ_ω ω sin(ω t_j) is the exact
discrete kernel −(−1)^j/(jΔ²), so Z₊ + Z₋ is an identity of the sampled data. All 45 published shell sums then give
|Z − 1| ≤ 3.4×10⁻⁵ (our low-noise traces ≤ 5.1×10⁻⁵), against 0.1–0.4% on Codex's ±20 grid. The ±20 residual is the
weight of the response beyond |ω| = 20: on that grid undamped test lines give ≤ 1.9×10⁻⁵, Gaussian-damped lines
≤ 2×10⁻⁶, and Lorentzian-damped lines e^{−Γt/2} give −(2/π)Γ/20, i.e. 1.0–3.6% for Γ = 0.3–1. Read that way, the
published residuals correspond to effective widths of 0.01–0.12 per shell. None of this says anything about
quasiparticles, as pt_response concluded.

**The initial slope.** For q″ = −Ω²(s)q with q(0) = 0 and q′(0) = 1, q(Δ)/Δ = 1 − Ω₀²Δ²/6 − (Ω²)′₀Δ³/12 + O(Δ⁴), and
for the shell-trace estimator with a white source the expected Ω₀² is κ̄²(shell) + ½tr H̄_X(s0), with H̄_X the
spatially averaged X Hessian. Our 20 low-noise slopes (five phases, shells 1–4) match this to ≤ 1.4×10⁻⁴ (≤ 3.3×10⁻⁴
without the drift term). Codex's 45 single-probe slopes match to ≤ 1.3×10⁻³, with the largest misses in shell 1 (18
modes), where the spatial inhomogeneity of H enters the trace as noise. Our shell mode counts equal Codex's (18, 62,
98, 210, 350, …) exactly. Two consequences:
- ½tr H̄_X = 2λ⟨r²⟩ + g²⟨φ²⟩ is *identically* Codex's Hartree m², so "slopes consistent with ω² = k² + Hartree m²" is
  an O(Δ³) tautology, not evidence for a dispersion relation.
- It swings with the orbit phase: 25.0 (s0 = 20), 25.3 (20.12), 22.0 (20.24), 15.7 (20.44), 17.7 (20.68). A
  mass read from the first sample is a phase sample, not a medium constant.

## New issues

1. The IR ratio range is 0.28–1.31 over one radial period. It is orbit-locked and depends on the delay window (0–2:
   0.16–0.58; 2–6: 0.30–1.72).
2. The meanH (HFB) closure is gapped, with a Goldstone defect of +5.0. This explains the IR poles being 12% high. The
   bare closure is tachyonic (−6.0), which explains most of the "parametric" instability; the genuine band is k² ≲ 4
   at a rate ≤ 0.30–0.38.
3. Goldstone-repaired closures fix the IR poles (median 3.5%) but break the UV ratio (×1.9–2.7). No homogeneous closure
   tested fits both.
4. meanH reproduces the phase pattern of the IR ratio (correlation 0.87) but not its apocentre amplitude.
5. The UV pole agreement is ≤ 5.5% at k = 2.6 by pt_response's own fits.
6. The J_r comparison used a non-like-for-like weight. The like-for-like excess is +0.5–1.5%, J_r drifts 25–41% per
   window, and the quasiparticle-basis total is 1.7625–1.7668.
7. The identity J_r − A₋ ≥ 0 is derived; the excess is the harmonic weight.
8. The canonical-sum residual is weight beyond |ω| = 20; on the Nyquist grid it is ≤ 3.4×10⁻⁵ (an identity of the
   sampled data). Lorentzian widths of 0.3–1 would give 1–3.6%, so the published 0.1–0.4% sizes the effective widths
   at 0.01–0.12.
9. The initial slopes are the O(Δ³) identity 1 − Δ²(κ̄² + ½tr H̄_X(s0))/6 to 1.4×10⁻⁴ (low-noise) and 1.3×10⁻³ (Codex's
   single probes). ½tr H̄_X is exactly the Hartree m² and swings 15.7–25.3 with phase, so "consistent with
   ω² = k² + Hartree m²" is not a dispersion test.
10. The charge identity holds to 1.0×10⁻³ (2×10⁻⁴ only in two of 18 windows). The radial antiparticle fraction is
    99.35–99.93% (action-weighted), 99.75–99.96% (amplitude²).

## Limits

- **Coverage.**
  - One box (L12, N48) and one radial period of phases, using shells 1–4 for the low-noise traces.
  - The closures are homogeneous, X block only (φ mixing neglected), and use the recorded spatial-mean Hessian up to
    s = 26.8.
  - The Goldstone repairs are diagnostic constructions, not a derived 2PI closure.
- **Toy.** The 3-mode toy is a bare-condensate truncation; its k² = 3 case is inside the parametric band and is not
  informative.
- **Physics not covered.** The quantum half, physical n̄, f_sym, binding and Candidate 1 are untouched.

## Reproduction

```bash
source $SP/review/pt_response_skeptic/env.sh
T=research/reality_framework/grut_formation_review_20261025/tracks/pt_response/skeptic
P=$SP/f55/research/reality_framework/grut_particle_transport_20261027; A=$SP/review/pt_response; B=$SP/review/pt_response_skeptic
(cd $PKC/code && python3 -B $OLDPWD/$T/code/make_states.py && python3 -B $OLDPWD/$T/code/make_states.py late)
# probes, chunked (about 4 s per sample per 16 tangents); repeat each until FINISHED
for s in 20.12 20.24 20.68; do (cd $PKC/code && python3 -B $OLDPWD/$T/code/probes_phase.py $s 8 $B/probes 75); done
(cd $PKC/code && python3 -B $OLDPWD/$T/code/probes_phase.py 20.92 8 $B/probes 75 424242 A)
cd $B/an
python3 -I $T/code/analyse_phase_probes.py phase_probes_all.json $B/probes/PROJ_skp_s20.{12,24,68}_M8_seed424242.npz $A/probes/PROJ_hutch_s20{,.44}_M8.npz
python3 -I $T/code/analyse_phase_probes.py phase_probes_2092.json $B/probes/PROJ_skp_s20.92_M8_seed424242_A.npz
python3 -I $T/code/exact_shell1.py $A/probes/PROJ_exact1_s20.npz exact_shell1.json
python3 -I $T/code/emulate_whitenoise.py $A/probes/PROJ_exact1_s20.npz emulate_whitenoise.json
python3 -I $T/code/measured_traces.py $A/probes $B/probes meas_traces.npz meas_traces.json
python3 -I $T/code/closures_phase.py $A/hessian_s20_26.8.npz meas_traces.npz closures_phase.json
python3 -I $T/code/closures_uv.py $A/hessian_s20_26.8.npz $P closures_uv.json
python3 -I $T/code/circular_bare.py $A/hessian_s20_26.8.npz $P circular_bare.json      # about 100 s
python3 -I $T/code/transform_free_checks.py $P transform_free.json
python3 -I $T/code/transform_nyquist.py $P meas_traces.npz transform_nyquist.json
python3 -I $T/code/initial_slope_check.py $A/hessian_s20_26.8.npz $P meas_traces.npz initial_slope.json
python3 -I $T/code/jr_identity.py $P jr_identity.json
python3 -I $T/code/counting_checks.py $P counting_checks.json
python3 -I $T/code/toy_exchange.py toy_exchange.json
```
