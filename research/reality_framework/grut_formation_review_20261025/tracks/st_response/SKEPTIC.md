# st_response_skeptic: adversarial check of track st_response (Codex b73de64, linear phase-seed response)

Claude's review of Codex b73de64 (`grut_source_timing_20261026/`), issue #13. 27 October 2026.
The target is st_response's NOTES.md and results. Code and results are in `tracks/st_response/skeptic/` and scratch
output in `$SP/review/st_response_skeptic/`.

**Bottom line.** We found no concrete basis to refute any st_response number or verdict. Our own third formulation
reproduces every bracket and first continuous crossing that st_response reports. It also closes st_response's main
NOT CHECKED item: over all 18 of Codex's momenta, the earliest continuous unit crossing of each history lies inside
Codex's sampled bracket.

The 10% correction is confirmed over all 18 k. That artefact also propagates into two places st_response did not
flag: SOURCE_CLOCKS.json and TIMING.png. Two of st_response's statements are worded too strongly or too weakly:

- **Too strong.** The attribution "the residual is Codex's integration noise" is not established. A ~5e-9 to 1e-8
  floor is shared by all three implementations.
- **Too weak.** The gauge-term caveat is too pessimistic. In the charge channel, and at k = 0.05–0.8, the e = 1e-7
  term is tested to about 1e-4 of its size.

Everything here is a linear extrapolation on a supplied model. It is not a nonlinear history and not evidence that
the fields exist.

## 1. Methods

### Independent code (sk_tangent.py)

A third formulation, deliberately different from both Codex's and st_response's:

- **Variables.** Cosmic time τ = H_i t and physical complex fields z = (X, Y, Σ), in second-order form:
  δz̈ = −3Hδż − (k²/a²)δz − dF[z](δz) − ξR(δX, 0, 0).
- **No hand-written Hessian.** dF[z](δz) is a forward-mode dual-number derivative of our own force F = ∂V/∂z*. F
  matches a numerical gradient of V to 1.2e-12, and the dual Jacobian matches finite differences to 6e-11.
- **Simplified geometry.** H = √(ρ/3) and R = ρ − 3p. We drop f = 1 − 2ξη²|X|² and the ξη² terms in R, which are
  O(1e-11). The background still matches Codex's archive: a to 4.3e-11 / 1.1e-10 and H to 1.4e-10 / 4.7e-10.
- **e = 0.** The gauge term is omitted. Its size is measured in §3.
- **Seed.** A 4-point difference of rp.prepare. It matches Codex's archived seed to 4.8e-10.
- **Observables.** A = Re(X*δX)/|X|², N_X = 2a³Im(X*Ẋ) (canonical sign), and
  δN_X = 2a³Im(δX*Ẋ + X*δẊ).
- **Model inputs.** Only the model constants and z₀ come from rp.py [Q].
- **Runs.** All 18 of Codex's k (plus the old-archive k = 2.122989) on both histories, τ = 0–800. Samples are on a
  0.002 grid plus Codex's sample times. DOP853 at rtol 1e-11.

### Controls (own code)

| Control | Result |
|---|---|
| rtol 1e-12 (k = 2.0, 2.13, 2.4, to τ 160) | changes ≤ 9e-13 of the envelope |
| Nonlinear two-point antisymmetric finite difference (eps 1e-6; F(z+d) − F(z) evaluated directly, no linearisation at all; k = 2.13, 2.4 to τ 200) | agrees to 6.4e-9 / 6.8e-9; identical first crossings at 0.002 resolution |
| Swapping in st_response's or Codex's seed | changes our differences by < 4e-11 |

### Reproduction (Codex's response.c, run in our own work dir)

k = 2.13 on (3, .5), at rtol 2e-10 and 2e-12.

### Tool use of st_response's code (labelled)

st_tangent.py at e = 0 and e = 1e-7 for k = 0.1, 0.2, 0.8, used only to attribute the low-k differences in §3.

## 2. Checks of st_response's claims

Agreement is measured as max |Δ| / the running envelope of the archive, with all 18 k.

### Codex's archive against our e = 0 code

**Amplitude:**
- k ≥ 1.2: ≤ 3e-8 on (3, .5) and ≤ 8.2e-8 on (5, 1).
- k = 0.02: 1.0e-8 and 3.6e-8.
- k = 0.05–0.8: 5.6e-8 to 1.1e-6. This is the gauge term (§3).

**Charge:**
- k ≥ 1.2: ≤ 2.4e-8 and ≤ 7.8e-8.
- k ≤ 0.8: 2.3e-7 to 4.2e-6 (gauge term).

**Background:**
- Late C: 1.98580528 and 1.47586779.
- Final a: 56.360618 and 48.599254.
- ⟨N_X⟩ over 64–80 on Codex's grid: 8156.643021831 and 5236.090504012.
- Change of N_X to τ 800: −0.4974% and −0.5037%.

### Fine first continuous crossings (max over all 18 k)

Each is compared with Codex's bracket, and our sampled brackets are identical to Codex's.

| History | Seed | Amplitude 0.1 | Amplitude 1 | Local charge 1 |
|---|---|---|---|---|
| (3, .5) | 0.03 | **56.458 (k 2.13); Codex [82.5, 83.0]** | 133.602 (k 2.13) ∈ [133.5, 134.0] | 119.580 ∈ [119.5, 120.0] |
| (3, .5) | 1e-6 | 563.468 ∈ [563.0, 563.5] | 771.942 ∈ [771.5, 772.0] | 737.422 ∈ [737.0, 737.5] |
| (5, 1) | 0.03 | 55.168 (k 0.02) ∈ [53.07, 55.29] | 147.322 (k 2.0) ∈ [147.0, 147.5] | 134.414 ∈ [134.0, 134.5] |
| (5, 1) | 1e-6 | 576.812 ∈ [576.5, 577.0] | 781.812 ∈ [781.5, 782.0] | 751.072 ∈ [751.0, 751.5] |
| both | 1e-9 | none | none (max 1e-9·\|A\| = 1.59e-3 / 1.44e-3) | none (max 2.6e-3 / 1.9e-3) |

- **Runners-up on (3, .5), 0.03 seed.** Amplitude 1: k 2.0 at 133.678, k 2.4 at 133.708. 10% level: k 2.0 at
  56.496, k 2.4 at 56.544.
- **Pre-80 maxima (fine grid).**
  - (3, .5): amplitude 0.03·max|A| = 0.2018, against 0.0887 on Codex's grid (under-read 2.3×); charge 0.343.
  - (5, 1): 0.118 for the amplitude and 0.188 for the charge.
- **Grid.** Codex's grid spacing exceeds 0.5 from τ ≈ 12.1 to 80, growing from 2.1 to 3.2 between τ 50 and 80.

### Verdicts on st_response

| Claim | st_response | Skeptic recomputation (own code unless stated) | Refuted? |
|---|---|---|---|
| R1 | [133.5, 134], fine 133.602 | same; earliest over all 18 k is k 2.13 | no |
| R2 | [147, 147.5], fine 147.322 | same; earliest over 18 k is k 2.0 | no |
| R3 | 771.942 / 781.812 | same | no |
| R4 | 119.58 / 134.414; 1e-6: 737, 751 | same | no |
| R5 | none; 1.6e-3 (3 k) | none over all 18 k; 1.59e-3 / 1.44e-3 (charge 2.6e-3 / 1.9e-3) | no |
| R6 (CORRECTED) | 56.458 (k 2.13) | 56.458 is the earliest of all 18 k; the nonlinear FD gives the same | no (correction stands; see the framing note below) |
| R6b | 55.168; 563.468; 576.812 | same, over all 18 k | no |
| R7 sign | ratio −0.9999999986 | −1.0000000017; same-sign error 1.84, opposite 3.2e-8 | no |
| R8 | 5.5e-9 to 3.6e-8; "residual is Codex's noise" | numbers reproduced; the attribution is overstated (below) | no (wording) |
| R9 | 1.98580528 / 1.47586779 | same (simplified background) | no |
| R10 | −0.497% / −0.504% | −0.4974% / −0.5037% | no |
| R11 | gauge test "~10% of a 1e-6 effect" | the caveat is too pessimistic (§3) | no (refined) |
| R12 | own controls only | Codex's own code at rtol 2e-12: 1.57e-8 change, identical k 2.13 brackets [Reproduction] | no |
| N1 | seed 2.5e-10 | our 4-pt seed vs Codex 4.8e-10; the seed choice is immaterial | no |
| N2 | 65.16–80, 6 samples; +1.1e-5 / −1.2e-6 | 65.158–80, 6 samples; +1.12e-5 / −1.17e-6 | no |
| X1 | all 18 k NOT CHECKED | now checked: no earlier max-over-k unit crossing outside the brackets | – |

Notes on the verdicts:

- **R6 framing [A].** Codex's PROTOCOL and MATHEMATICS define the bracket as the first *sampled* crossing. They warn
  explicitly that narrow excursions between samples can precede it, so [82.5, 83.0] is correct by Codex's own
  definition. It is nevertheless not a bracket of the first crossing, which is 26.5 τ earlier.
  - st_response also "corrects" the per-column (5, 1) k = 2.0 and 2.13 10% brackets ([80.0, 80.5] against 67.11).
    Our numbers agree: 67.114 and 67.092. But Codex never reports per-column brackets, so that is a correction of an
    unreported quantity.
- **R8 attribution.** We reran Codex's own response.c for k = 2.13.
  - **Codex against itself.** Codex's code at rtol 2e-12 differs from its archive by 5.7e-9, and from its rtol 2e-10
    one-k run by 1.6e-8.
  - **Against st_response.** Codex's code at rtol 2e-12 differs from st_response by 1.1e-8, which is *more* than the
    archive does (5.5e-9).
  - **Our code (e = 0) against the others.** It differs from st_response's own e = 0 run by 4.9e-9 (τ ≤ 240) and
    1.2e-8 (τ ≤ 800). Both codes have integration errors of about 1e-12.
  - **Conclusion.** A floor of ~5e-9 to 1e-8 is common to all three implementations. Its source is not seeds (swap
    test) and not isolated. Our simplified geometry may contribute to our share. The remark "our own error is ~1e-12,
    so the residual belongs to Codex's integration" therefore rests only on st_response's integration controls, not on
    formulation differences. This is immaterial: the smallest bracket margin is 1.5e-3.
- **Growth-rate remark.** The "0.01545368 /τ in both codes" is a property of Codex's sample grid. On a 0.002 grid
  the running-envelope rate over 200–700 is 0.015469 (3, .5) and 0.016129 (5, 1). On Codex's grid we get 0.0154537 and
  0.0162003, matching the archive.

## 3. Gauge/longitudinal term (refinement, in st_response's favour)

The low-k differences are the gauge term. Running st_response's code at e = 0 against e = 1e-7 for k = 0.1, 0.2, 0.8,
the effect of e = 1e-7 on the amplitude is 2.4e-7, 7.3e-7 and 1.1e-6 of the envelope. Our e = 0 code matches
st_response's e = 0 to 1.9e-8, 1.1e-8 and 9.7e-9.

The charge channel feels it more. The e = 1e-7 effect is 2.4e-7 (k 0.02), 4.2e-6 (0.1), 3.5e-6 (0.2) and 2.0e-6
(0.8). st_response's e = 1e-7 charge matches Codex's archive to 2.8e-11, 2.8e-11, 3.7e-11 and 1.4e-10.

So the gauge term in Codex's direct-vector code is independently confirmed to about 1e-4 of its own size. That is far
better than st_response's "roughly 10% of a 1e-6 effect", which is limited by the amplitude channel at k = 0.02 and
2.13, where the effect sits at the ~1e-8 inter-code floor. This bears on Codex's request to check the direct-vector
adaptation, but it remains a test at a tiny coupling: e²/k² ≤ 2.5e-11.

## 4. New issues with the Codex packet (not raised by st_response)

1. **SOURCE_CLOCKS.json inherits the 10% artefact.**
   - The quiet3 `at_response_crossings` entry for the 0.03 rad amplitude_10percent is evaluated at τ = 83.0
     (a = 18.13). There the Y_L rest-clock number fraction is 8.35e-6 and the injected conformal-energy fraction
     1.81e-8.
   - At the first continuous crossing, τ = 56.458, these are 5.68e-6 and 1.01e-8 [C, own a(τ); γ = 1.00590e-7 quoted
     from SOURCE_CLOCKS].
   - The error goes in the direction that *strengthens* Codex's ordering. The REPORT does not quote these numbers.
2. **The upper-right panel of TIMING.png is aliased before τ = 80.** plot.py draws the max-over-k response on the
   archived samples. On (3, .5) that curve peaks at 0.089 before τ 80 against a true 0.202, so the dotted 10% line
   appears crossed near τ 83 instead of 56.5. The figure is visual only, and the unit crossings are unaffected.
3. **The source-clock values at the unit samples are reproduced** [C, own a(τ) and closed-form tail; widths quoted].
   - Values: 1.3479e-5 / 3.7161e-8 (τ 134) and 1.4837e-5 / 4.3067e-8 (τ 147.5). Codex's PCHIP-on-sparse-grid values
     agree to 2.3e-7 relative. Mean a per decay: 5568.6529 and 4800.7110.
   - At the continuous crossings the values are 1.3439e-5 / 3.700e-8 and 1.4819e-5 / 4.299e-8. This changes nothing,
     but the REPORT's "1.35–1.48e-5, 3.72–4.31e-8" are sampled-time values.

## 5. Limits

- **Same model and inputs.** The independence is in formulation and code only.
- **Our approximations.** e = 0 and the simplified geometry: their effects are measured above.
- **Not checked.** The decay medians (except the Y_L fractions above), the vector clocks, the quartic bands and the cube
  charges, as well as the gauge3 and gauge5 histories beyond the archive differences.
- **Third-party code.** Codex's code was run only for k = 2.13 on (3, .5).

## 6. Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RP_SCRATCH=$SP/review/st_response_skeptic/rp
cd tracks/st_response/skeptic/code
python3 -B sk_tangent.py                      # force / dual-Jacobian self-test
python3 -B sk_runs.py q3 & python3 -B sk_runs.py q5 & wait     # ~215 s each; rerun until DONE
python3 -B sk_runs.py q3fd; python3 -B sk_runs.py q3tight
python3 -B sk_analyze.py                      # -> results/sk_response.json
python3 -B sk_seedtest.py; python3 -B sk_clock.py
python3 -B sk_lowk_probe.py 0; python3 -B sk_lowk_probe.py 1e-7   # uses st_response's st_tangent.py as a tool
# Reproduction of Codex response.c (own work dir, GRUT_TIMING_WORKDIR=$SP/review/st_response_skeptic/twork): python3 -B -I sk_codex_tight_reproduction.py
```

Results: `skeptic/results/sk_response.json`, `sk_seedtest.json`, `sk_codex_tight_k213.json`, `sk_clock.json`,
`sk_lowk_probe.json`, `sk_gauge_charge.json`, `sk_charge_agreement.json`, `sk_pre80.json`, `sk_verdicts.json`.
