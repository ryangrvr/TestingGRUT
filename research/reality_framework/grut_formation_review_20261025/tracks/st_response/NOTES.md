# st_response: independent linear phase-seed response on the quiet reheating histories

Track of Claude's review of Codex b73de64 (`grut_source_timing_20261026/`), issue #13. 27 October 2026.
Scope: claim 1 of BRIEF_TIMING.md (linear response crossing times) and the response-side items in Codex's handoff
(direct-vector equations, phase-current sign convention). The decay clocks, the quartic Floquet spectrum and the
cube charges are not covered here.

**Bottom line.** Our own code reproduces Codex's archived linear responses to 8–9 significant digits. It confirms every
unit-level crossing bracket that it recomputes, and finely sampled runs place each first continuous unit crossing
inside its sampled bracket. One secondary number is wrong. The quiet (3, 0.5) 10% warning bracket for the 0.03 rad
seed, sampled τ = [82.5, 83.0], comes from the sparse archived grid below τ = 80. The response first reaches 10% at
τ ≈ 56.46, between archived samples 55.29 and 57.61. All of this is a linear extrapolation on a supplied model. It is
not a nonlinear history and not evidence that these fields exist.

## 1. Method (independent code; nothing from Codex's modes.py, response.c or timing_common.py is imported)

- **Background [C].** These are our rp_reheating engine's equations (`tracks/rp_reheating/code/rpreheat.c`, mode 0:
  an oscillating inflaton with Γ friction and a nonminimal ξ_X = −5/6). We transcribed them to Python in
  `code/st_tangent.py` and integrate them together with the tangent. H is rebuilt from the Friedmann constraint at
  every evaluation, while Codex evolves H. Raychaudhuri is integrated into H_aux as a check. Parameters and the
  preparation come from `rp.py` (`rp.P`, `rp.XI`, `rp.ETA2`, `rp.prepare`, `rp.potential`).
  Transcription check against the rp C engine over τ = 80–800 at 0.5 spacing: relative differences are 2.5e-13 in a,
  1.0e-14 in H, 1.4e-11 in N_X and 2.7e-13 in |X|. The rp `extra` yield yr_b3_g0.5 reruns as 1.04542056 (committed
  value 1.0454206).
- **Mode equations [D].** For χ = a r, r = √2(Re, Im)(X, Y, Σ), conformal canonical q and p, g = Tχ, g' = Tχ' and
  T = diag(J, 3J, 6J):
  H₂ = p²/2 + qᵀ(k²+M)q/2 + (e²/2k²)(g·p − g'·q)², with M = a²(V_rr + ξ_X R P_X) − (a²R/6)·1.
  The equations are integrated in τ with a factor 1/a.
  - We rederived the A0 elimination from D₀χ = χ' − (a'/a)χ − eA₀Tχ. The (a'/a) cross terms cancel, so the
    constraint current is g·q' − g'·q, matching MATHEMATICS.md.
  - V_rr is our own Wirtinger Hessian of `rp.potential`, with A = ∂²V/∂z*∂z, B = ∂²V/∂z*∂z* and block
    [[Re(A+B), −Im(A−B)], [Im(A+B), Re(A−B)]]. It matches a Richardson finite difference of rp.potential's gradient
    to 2.3e-12. It is used directly, with no Ward-identity projection, which is a different numerical route from
    Codex's C engine.
  - e = 1e-7 is quoted [Q] from `tracks/hw_gates/code/spotmodel.py`.
- **Seed [C, D].** q₀ = dχ/dδ at a = 1, from a central difference of `rp.prepare(b, delta=±1e-5)`. Here Σ has phase
  −π/2 + δ, X stays real and positive, and Y and the radii are re-solved. We set p₀ = q₀.
  - Why p₀ = q₀ [D]: every prepared state has zero velocity, so δṙ = 0. There is no metric perturbation, so
    δχ' = a²(δṙ + Hδr) = δr = q₀ at a = H = 1. At τ = 0, g' = g, so the seed current g·p₀ − g'·q₀ vanishes and
    q' = p₀.
  - Our q₀ matches Codex's archived seed to 2.5e-10 absolute (6e-12 of its norm). Codex's archive also has p₀ = q₀.
  - The fractional amplitude starts at A(0) = −0.0430432.
- **Observables [C].**
  - Amplitude per radian: A = χ_X·q_X/|χ_X|².
  - Charge: N_X = χ₀χ₁' − χ₁χ₀', which equals q₀p₁ − q₁p₀ for the background. Its perturbation is
    δN_X = δq₀χ₁' + χ₀δp₁ − δq₁χ₀' − χ₁δp₀, using the canonical, gauge-covariant p.
  - Charge fraction = δN_X / |⟨N_X⟩₆₄–₈₀|, a trapezoid over the archived sample grid, as Codex does.
- **Integrator.** scipy DOP853 with rtol 1e-10, atol 1e-13 and max step 0.05, on Codex's sample grid to τ = 800 so that
  the samples coincide. Own controls:
  - rtol 1e-12: change ≤ 2.6e-13 of the running envelope through τ = 240.
  - RK45 at rtol 1e-11: change 3.0e-13 through τ = 160.
  - k = 2.13 alone versus inside a three-k run: 7.4e-13 through τ = 800.
  - Gauss charge ≤ 7.5e-13 and |H_aux/H − 1| ≤ 6.3e-13 in every run.

## 2. Comparison with Codex's archive (quiet3.npz, quiet5.npz, gauge3.npz; read-only) [C]

"Envelope" means max over τ of |Δ(τ)| / max_{τ'≤τ}|A_Codex(τ')|. "Pointwise" excludes near-zeros, where
|A| < 1e-3 × envelope.

| History, k | Amplitude: envelope / pointwise | Charge fraction: envelope | max over τ of |A|, ours vs Codex |
|---|---|---|---|
| (3, .5), 2.13 | 5.5e-9 (at τ 772) / 9.7e-8 | 3.2e-10 | 1583890.118 vs 1583890.119 |
| (3, .5), 0.02 | 9.5e-9 / 1.3e-6 | 2.8e-11 | 12.4973349 vs 12.4973349 |
| (3, .5), 5.8 | 8.2e-9 / 3.1e-7 | 8.4e-10 | 7.9308932 vs 7.9308932 |
| (5, 1), 2.0 | 2.1e-8 / 3.6e-7 | 1.1e-9 | 1441139.089 vs 1441139.090 |
| (5, 1), 2.13 | 2.0e-8 / 1.1e-7 | 1.2e-9 | 1085865.639 vs 1085865.638 |
| (5, 1), 0.02 | 3.6e-8 / 4.4e-7 | 4e-11 | 14.017174 vs 14.017174 |
| (3, .5) at e = 1e-6, k = 0.02 / 2.13 / 5.8 | 8.3e-8 / 4.8e-8 / 7.2e-8 | ≤ 7.4e-9 | (gauge3.npz) |

Our controls put our own error near 1e-12, so the residual of 5e-9 to 8e-8 belongs to Codex's integration. Codex
reports its own tolerance spread as 1.8e-8.

Background, ours vs Codex:
- a agrees to 3.5e-11 (3, .5) and 1.2e-10 (5, 1); H to 1.4e-10.
- ⟨N_X⟩₆₄–₈₀: 8156.643022496 vs 8156.643022495, and 5236.090504367 vs 5236.090504368.
- Late C = a²H at τ = 800: 1.98580528 vs 1.98580528, and 1.47586779 vs 1.47586779.
- Final a: 56.360618 and 48.599254, both matching Codex.
- Change of the unperturbed N_X from its 64–80 mean to τ = 800: −0.497% and −0.504%.
- The late envelope growth rate over τ = 200–700 at k = 2.13 on (3, .5) is 0.01545368 /τ in both codes.

### Crossing brackets (first sampled seed·|value| ≥ level; our brackets equal Codex's per-column brackets in every case)

"Fine" is the first continuous crossing found at 0.002 spacing (`code/fine.py`).

| History, k | Seed | Amplitude 0.1 | Amplitude 1 (a) | Local charge 1 | Fine first crossing (amplitude 1 / charge 1) |
|---|---|---|---|---|---|
| (3, .5), 2.13 | 0.03 | [82.5, 83.0]* | [133.5, 134.0] (23.0527) | [119.5, 120.0] | 133.602 / 119.580 |
| (3, .5), 2.13 | 1e-6 | [563.0, 563.5] (fine 563.468) | [771.5, 772.0] (55.3653) | [737.0, 737.5] | 771.942 / 737.422 |
| (5, 1), 2.0 | 0.03 | [80.0, 80.5]* | [147.0, 147.5] (20.8777) | [134.0, 134.5] | 147.322 / 134.414 |
| (5, 1), 2.0 | 1e-6 | [576.5, 577.0] (fine 576.812) | [781.5, 782.0] (48.0495) | [751.0, 751.5] | 781.812 / 751.072 |
| (5, 1), 2.13 | 0.03 | [80.0, 80.5]* | [147.0, 147.5] | [134.0, 134.5] | 147.342 / 134.456 |
| (5, 1), 2.13 | 1e-6 | [614.5, 615.0] (fine 614.840) | [782.0, 782.5] | [759.0, 759.5] | 782.440 / 759.418 |
| all of the above | 1e-9 | none by τ 800 | none | none | the largest seed·|A| is 1.6e-3 |

Sampled values at the unit crossings:
- 0.03 seed: 1.48497897421 vs Codex 1.48497897416 on (3, .5), and 1.2154560363 vs 1.2154560366 on (5, 1) at k = 2.0.
- 1e-6 seed: 1.0320225039 vs 1.0320224982, and 1.1035231758 vs 1.1035231947.

\* **The sparse grid below τ = 80 misses the 10% level.** Codex's sample grid there is the old MODE_PROFILES grid,
spaced about 2.3–3 in τ. On a 0.002 grid, the 0.03 rad seed first reaches 0.1:
- on (3, .5), at τ = 56.458 (k = 2.13) and 56.544 (k = 2.4), between archived samples 55.293 and 57.609. It reaches
  0.20 by τ ≈ 69.6. The archive's max-over-k bracket is [82.5, 83.0] (k = 2.4).
- on (5, 1), at τ = 67.11 (k = 2.0 and 2.13). The archive's max-over-k bracket, k = 0.02 at [53.07, 55.29], still
  holds: the fine crossing for that mode is τ = 55.168.

The unit crossings are not affected. Each fine first crossing lies 0.04–0.43 before its sampled upper end. Before
τ = 80 the largest value of 0.03|A| is 0.20, so no earlier excursion reaches unity.

### Sign convention [C]

- The archive's charge fraction agrees with our canonical q₀p₁ − q₁p₀ at the **same** sign: the dot product is +1, and
  flipping the sign would give a relative error of 2.0. Normalisations are positive (+8156.64 and +5236.09).
- The old 2cbc7f1 PHASE_PROFILES archive at k = 2.122989 through τ = 80 agrees only with the **opposite** sign. The
  least-squares ratio is −0.9999999986, with residuals of 5e-8 at float32 precision. The amplitude row agrees to 4.8e-8.
- Codex's stated conversion is therefore correct.
- For the record, three conventions are in play:
  - 2cbc7f1 MATHEMATICS.md writes N_X = −(Jχ)·χ', which is −(q₀p₁ − q₁p₀).
  - Our rp engine's NX = −2a³Im(X*Ẋ) has that same sign.
  - The new archive uses the canonical sign.
- The signed charge fraction flips with the convention. The |·| crossing brackets do not.

### Grid and gauge sensitivities [C]

- **Window mean.** On the archived grid the "64–80" trapezoid uses 6 samples and actually spans τ = 65.16–80. On a
  0.002 grid over 64–80 it changes by +1.1e-5 on (3, .5) and −1.2e-6 on (5, 1). No bracket moves: the smallest margin among the sampled unit charge values is
  1.0015 − 1 = 1.5e-3.
- **e = 1e-7 against e = 0.** The response changes by 4.2e-9 of the envelope (k = 2.13, τ ≤ 240) and 8.8e-9
  (k = 0.02, τ ≤ 800).
- **e = 1e-6 against e = 1e-7.**
  - Size of the effect, ours vs Codex: 8.7e-7 vs 8.9e-7 at k = 0.02, 4.14e-7 vs 4.14e-7 at k = 2.13, and 1.45e-7 at
    k = 5.8 in both.
  - Our effect curve agrees with Codex's to 8% at k = 0.02 and 15% at k = 2.13. At k = 5.8 the curves are at the
    numerical noise level and are not resolved.
  - So at these couplings the gauge/longitudinal term is barely exercised. The agreement tests it only to roughly 10%
    of an effect of size 1e-6.

## 3. Verdicts

| # | Codex claim | Codex | Ours | Verdict |
|---|---|---|---|---|
| 1 | 0.03 rad seed, unit fractional modulation of \|X\|, (3, .5) | sampled τ 134, a 23.053 | [133.5, 134.0], a 23.0527; fine 133.602 (k 2.13) | CONFIRMED [C] |
| 2 | Same, (5, 1) | 147.5, a 20.878 | [147.0, 147.5], a 20.8777; fine 147.322 (k 2.0) | CONFIRMED [C] (only k = 2.0 and 2.13 rerun) |
| 3 | 1e-6 rad seed, unit modulation | 772 / 782 | [771.5, 772.0] / [781.5, 782.0]; fine 771.942 / 781.812 | CONFIRMED [C] |
| 4 | Local charge modulation reaches 1 earlier (0.03 rad) | 120 / 134.5 | [119.5, 120.0] / [134.0, 134.5]; fine 119.580 / 134.414 | CONFIRMED [C] |
| 5 | 1e-9 rad seed crosses neither threshold before τ 800 | none | none; largest seed·\|A\| 1.6e-3 | CONFIRMED [C] (for the k values computed) |
| 6 | 10% warning bracket, 0.03 rad, (3, .5) (RESPONSE_RESULTS.json) | [82.5, 83.0] (k 2.4) | first continuous crossing τ 56.46 (k 2.13), 56.54 (k 2.4) | CORRECTED [C]: a sampling artefact of the sparse pre-80 grid |
| 6b | 10% bracket, 0.03 rad, (5, 1); 1e-6 rad 10% brackets | [53.07, 55.29]; 563 / 576.5 | fine 55.168; 563.468 / 576.812 | CONFIRMED [C] |
| 7 | Phase-current sign: canonical now, negative in the old archive | stated | ratio to the old archive −0.9999999986 | CONFIRMED [C] |
| 8 | The direct-vector equations reproduce the response (profiles) | 7 internal checks | 5.5e-9 to 8e-8 of the envelope, independent code | CONFIRMED [C] |
| 9 | Late C = a²H/H_i | 1.986 / 1.476 | 1.98580528 / 1.47586779 (own background) | CONFIRMED [C] |
| 10 | Unperturbed N_X changes about −0.5% from the 64–80 mean to τ 800 | ≈ −0.5% | −0.497% / −0.504% | CONFIRMED [C] |
| 11 | e = 1e-6 has little effect at the peak | "little" | 4.1e-7 of the envelope at k 2.13 | CONFIRMED [C], with the exercise caveat above |
| 12 | A tighter tolerance leaves the brackets unchanged | stated | own rtol 1e-12 and RK45 controls ≤ 3e-13 | CONFIRMED [C] |
| – | Decay medians, injected fractions, quartic bands, cube charges | – | – | NOT CHECKED (other tracks) |
| – | The earliest crossing over all 18 k values | k 2.13 / 2.0 | only k = 0.02, 2.13, 5.8 (3, .5) and 0.02, 2.0, 2.13, 2.4 (fine) | NOT CHECKED beyond these k |

**Interpretation [A].** These crossings are breakdown times of a linear extrapolation for supplied seeds. We agree with
Codex's wording that they are not nonlinear breakup times. At the unit level, the first continuous crossing lies inside
each sampled bracket. Comparing these times with the decay clocks is outside this track. Nothing here selects a seed,
a covariance, a fragment population or an abundance.

## 4. Limits

- Same model and inputs as Codex: supplied Higgs-portal fields, couplings and the preparation family. The independence
  is in code and formulation only (constraint-H background, direct Hessian, Python/scipy integrator).
- The A0 elimination was rederived only in outline, as above. The previous review bca4927 checked it in full.
- At the physical e, the longitudinal term contributes at about 1e-8 to 1e-6. Agreement with Codex therefore does not
  test it strongly.
- We reran 7 of Codex's 42 histories: (3, .5) at k = 0.02, 2.13, 5.8 with e = 1e-7 and 1e-6; (5, 1) at k = 0.02, 2.0,
  2.13; and the old-archive k. The k = 2.4 mode was run only on the fine early grid.

## 5. Reproduction

```bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 RP_SCRATCH=$SP/review/st_response/rp
cd tracks/st_response/code
python3 -B st_tangent.py                       # potential/force/Hessian self-test
python3 -B runs.py q3_k213 q3_k3 q5_k2 q3_k213_tight q3_k213_rk45 q3_k213_e0 q3_gauge6 q3_k002_e0 q5_k002 q3_old
python3 -B compare.py                          # -> results/st_response.json (also reruns the rp C engine checks)
python3 -B fine.py q3 q5 q3e q5e               # -> $SP/review/st_response/fine/*.json (collected in results/st_fine_crossings.json)
```

Each 800-τ run takes about 75–100 s on one core. Raw samples are saved in `$SP/review/st_response/runs/*.npz`.
