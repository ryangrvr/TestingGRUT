# st_rerun: rerun of Codex b73de64 (source timing) with Codex's own code, plus a code audit

Track key `st_rerun`, 27 October 2026. This review covers Codex packet
`research/reality_framework/grut_source_timing_20261026/` at b73de64, using the extracted tree in
`$SP/b73/research/reality_framework/` (SP = /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad).

**What kind of check this is.** Almost every number here is a **reproduction**: Codex's scripts and C integrator,
run on our machine (Python 3.13, NumPy 2.5.3, SciPy 1.18.1, gcc 13.3 -O3). Codex used NumPy 2.3.5. Two things are
not plain reproduction:

- **Sensitivity reruns.** These still use Codex's own integrator, but with a dense time grid, a dense k grid, a tighter
  tolerance and a rescaled seed.
- **Own arithmetic.** Our own code computes the decay clocks and growth slopes, but it works on Codex's background
  arrays. It is not an independent solver.

No independent physics engine was written for this track. Claim labels follow BRIEF.md: [D] derived, [C] computed,
[Q] quoted, [E] estimate, [A] argued, [R] recalled. Nothing here validates or excludes the model. All of it is
same-model numerical work on supplied fields.

## 1. Housekeeping

- Before any run, we copied every published JSON/NPZ the scripts can overwrite to `$SP/review/st_rerun/published/`.
  The copy also took RESULT_TABLE.md, MANIFEST.json and TIMING.png, with SHA256SUMS. The published `quiet3.npz`
  matches `git show b73de64:` byte for byte.
- After the comparisons, we copied the published files back into `$SP/b73`, and all 16 checksums verify. Other tracks
  therefore see Codex's published outputs. We did not run plot.py. package.py needs a git checkout and stopped before
  writing anything.
- We wrote no git objects and made no commits. Raw binaries live in `$SP/review/st_rerun/{twork,mwork,dense}`.

## 2. Reproduction of verify_response.py, run_response.py, quartic_scan.py, source_clocks.py, check_late.py, summarize.py

Run times on our machine were 8.7 s, 5.2 s, 5.2 s, 0.7 s, 1.5 s and 0.7 s. The 42 tangent histories to τ = 800
take about 5 s of wall time with two threads, so the full run was used and no subset or reused npz was needed.

| Output | Leaves | Differences from the published copy |
|---|---:|---|
| quiet3/quiet5/gauge3/gauge5.npz (every array) | – | **bitwise identical**; `raw_sha256` of all four raw binaries also identical |
| RESPONSE_RESULTS.json | 313 | only the 4 wall-clock `seconds` |
| SOURCE_CLOCKS.json | 434 | none (bitwise) |
| RESULTS.json, RESULT_TABLE.md | 206 | 4 quartic μ values at 1e-15 relative; table identical |
| RESPONSE_CHECKS.json | 10 | short Radau: 3.397e-10 against 3.488e-10, nfev 20559 against 20407 (SciPy version) |
| QUARTIC_SCAN/CHECKS.json | 515/25 | stable-mode μ ~3e-13 (round-off noise); IR μ/k to 1e-11 relative |
| VERIFICATION_RESULTS.json | 65 | the same Radau and quartic noise; all 7 groups pass |

verify_response.py has three groups:

- **Direct vector against full transfer.** Scaled error 2.6613e-9 for the vector and 2.0111e-11 for the background.
  These are identical to Codex's.
- **Short Radau run.** 3.40e-10 here against Codex's 3.49e-10.
- **Published phase archive.** Amplitude 1.994e-8 and charge 3.377e-8, both identical to Codex's.

check_late.py reproduces three Codex results:

- The tight-tolerance response error is 1.789e-8 for β = 3 and 3.24e-9 for β = 5, with identical brackets.
- The clock checks give 1.1e-16, 2.2e-16 and 3.3e-16.
- The largest Friedmann residual is 3.73e-9 (in gauge3) and the largest Gauss residual is 6.94e-10. REPORT.md says
  "3.8e-9", which rounds 3.73e-9 up. That does not matter.

### Requested numbers (published = reproduced, all digits)

| History | Seed | 0.1 bracket (τ) | Unity bracket (τ, a) | Local-charge unity (τ) |
|---|---|---|---|---|
| quiet3 | 0.03 | [82.5, 83.0] k = 2.4 | [133.5, 134.0], a = 23.0527 | [119.5, 120.0] |
| quiet3 | 1e-6 | [563.0, 563.5] | [771.5, 772.0], a = 55.365 | [737.0, 737.5] |
| quiet5 | 0.03 | [53.069, 55.293] k = 0.02 | [147.0, 147.5], a = 20.8777 | [134.0, 134.5] |
| quiet5 | 1e-6 | [576.5, 577.0] | [781.5, 782.0], a = 48.050 | [751.0, 751.5] |
| both | 1e-9 | none | none | none |

The gauge3/gauge5 rows (e = 1e-6) have the same brackets as quiet3. For quiet5, k = 2.13 gives a 1e-6 unity bracket
of [782.0, 782.5], while k = 2.0 crosses first at 782.0.

Other reproduced values:

- **Late C and matched offsets.** late_C = 1.9858052758 (quiet3) and 1.4758677865 (quiet5). The matched offsets
  τ_* are 0.19367410 and −0.16906188.
- **Relative range of C over τ 720–800.** 4.0e-11 and 1.3e-10.
- **Change in the homogeneous light charge.** N_X changes by −0.497% (quiet3) and −0.504% (quiet5) between its
  64–80 mean and τ = 800.
- **Medians.** Y_L τ = 6.8908e6, with a = 5231.4 and 4510.0. Vector (e = 1e-6) τ = 4.3552e17, with a = 1.3152e9 and
  1.1338e9. Vector (e = 1e-7) τ = 4.3552e20, with a = 4.159e10 and 3.585e10.
- **Y_L fractions at the 0.03-rad unity samples.**
  - Decayed number: 1.3479e-5 (quiet3, τ = 134) and 1.4837e-5 (quiet5, τ = 147.5).
  - Conformal energy: 3.7161e-8 and 4.3067e-8.
- **Mean a per decay.** 5568.65 and 4800.71, against the asymptotic form √(πC/2γ). The old C = 1 value was 3951.68.
- **Quartic orbit.**
  - μ/k = 0.43980, 0.43979 and 0.43972 at k = 0.01, 0.02 and 0.05.
  - Bands: [0.01*, 3.135–3.140] with a peak μ of 0.63178 at k = 2.13 (*the lower edge is not bracketed, because μ → 0 as k → 0), and [5.57–5.58, 5.97–5.98] with a peak of 0.1602
    at k = 5.77.
- **Wavelength cube.** Q_cube at k = 2.13 is 209368.9 and 134402.7.

## 3. Code audit

**(a) response.c against modes.c [C, A]**

The diff contains three kinds of change: the stride change 12 → `NCOL` = 1, the removal of the `NEUTRAL_SECTOR`
branch, and the seed read in place of the identity initialization. The perturbation matrix, the Ward vector, the
Coulomb-gauge current terms (c = e²/k²), the 1/a factors and the DOP853 driver are character-identical. For column
`col`, the column operations read only that column, so the direct 12-vector right-hand side is exactly the
full-transfer generator applied to the seed column.

The Ward projection is kept: `mq = common*ward + M(T − common*g)`. Two numerical checks back this up:

- **Ward identity.** It holds to 3.9e-16 relative on all 227 sampled backgrounds per history over τ = 0–800, so the
  projection is an exact rewrite and not a change of equations.
- **Projected right-hand side.** It equals Codex's Python `generator(v,β,k)@x` to 6.0e-16 for random x at
  k = 0.02, 2.13 and 5.8.

The seed layout `np.tile(seed,(nk,1))` matches `T[i*NCOL]`.

DOP853 step control does include the tangent components. The RMS norm runs over all 19 + 12·n_k components, with
scale_j = atol + rtol·max(|y_j|, |y_new_j|) and atol = 1e-3·rtol. This matches SciPy's `_estimate_error_norm`, and the
step factor is clipped to [0.2, 5]. The k = 2.13 tangent grows by e^2.84, e^4.25, e^8.86 and e^14.58 by τ = 80, 134,
400 and 772. Once it has grown, rtol·|y| dominates atol by many orders, so control is relative and invariant under the
linear rescaling. That growth therefore does not loosen the tangent's relative accuracy. Background components share
the RMS with 216 tangent components, which dilutes their weight by up to √(235/19) ≈ 3.5.

Tests, using Codex's integrator:

1. **Seed scaling.** Seed ×1e-6 and ×1e6 reproduce response/s to 2.1e-9 and 9.7e-8 of the peak.
2. **Tolerance, dense grid.** rtol 2e-12 against 2e-10 on a dense grid to τ = 800 for all 18 k gives 6.9e-8 (amplitude)
   and 6.5e-9 (charge) of the per-k peak. All crossings agree to 0.01 τ. The 2e-12 run is at the round-off floor:
   12.9M steps and 2.3M rejections, against 0.19M steps at 2e-10.

**(b) responses() [D, C]**

The amplitude row χ_X·δq_X/|χ_X|² is the first-order δ|X|/|X|.

The local charge δN = δq0·p1 − δq1·p0 + q0·δp1 − q1·δp0 equals g_X·δp_X − g'_X·δq_X, with g = Jχ and
J = [[0,−1],[1,0]]. This agrees to 4.1e-15 on the archived vectors, and it is exactly the X part of the Gauss-law
current in response.c.

In H₂, q' = p + (e²/k²)g(g·p − g'·q), so the canonical p = q' − eδA₀g is the covariant D₀χ at linear order. The
canonical momentum therefore already contains the eliminated A₀, and δN is the gauge-covariant X-charge density. The
background uses p = χ', which is valid because the background's total gauge charge g·χ' ≤ 6e-10 while X alone is 8.1e3,
so A₀ = 0 there. Numerically, p and q' differ by ≤ 9e-7 relative at e = 1e-7.

Sign convention:

- The old 2cbc7f1 row is the exact negative of the new one: charge[:2] = [−p1, p0] and charge[6:8] = [q1, −q0].
- The old background charge −2a³Im(z*ż) equals −(q0p1 − q1p0) = −(Jχ)·χ'.
- Both conventions divide by |mean|, so the new fraction is −(old). verify_response adds the two and gets 3.38e-8.
- The crossings use |·|, so the sign does not affect them.

At τ = 0 the seed current is exactly 0, and dv[0]/dδ = 0, so the inflaton is not perturbed by the preparation phase.

**(c) crossings() [C]**

The function takes the first sample with max_k ≥ level and returns the bracket [t_{i−1}, t_i] and k at t_i. It has no
off-by-one error. Sample spacing is 0.5 τ for τ ≥ 80, but 1.89–3.22 τ for τ 40–80, where the inherited geomspace grid
is used.

The responses are bursty: δ|X|/|X| spikes at the background radial minima. For example, the 0.03-rad excursion above
unity lasts 1.36 τ (quiet3) and 1.08 τ (quiet5), and a 10% excursion near τ = 56 lasts only 0.52 τ.

We reran Codex's integrator on a dense time grid with Δτ = 0.01, at both rtol values, with identical results:

- **First crossings.**
  - quiet3: 0.03 rad gives unity 133.61 and charge 119.58; 1e-6 gives 0.1 at 563.47, unity 771.95 and charge 737.43.
  - quiet5: 0.03 rad gives 0.1 at 55.17, unity 147.33 and charge 134.42; 1e-6 gives 0.1 at 576.82, unity 781.82 and
    charge 751.08.

  All of these lie inside the published brackets. The quiet5 spike at τ ≈ 129.2 peaks at 0.842, below 1.
- **One miss.** The quiet3 0.03-rad 10% warning is first reached at **τ = 56.46** (k = 2.13), not [82.5, 83.0]. A
  0.52-τ excursion fell between the geomspace samples 55.29 and 57.61.

We also reran on a dense k grid (1.5–3.0 in steps of 0.02, dense in time):

- Amplitude-unity crossings move by ≤ 0.02 τ: 133.60 at k = 2.20 and 147.32 at k = 2.04.
- Local charge moves to 119.40 (k = 2.18), 134.15, 737.28 and 750.80, up to 0.28 τ earlier. Both 0.03-rad
  local-charge "by 120 / by 134.5" statements still hold.

**(d) source_clocks.py [C, D]**

- **Tail formula.** ∫_t^∞ γ√(2C(s−τ_*)) e^{−γs} ds = √(2C/γ) e^{−γτ_*} Γ(3/2, γ(t−τ_*)), coded as
  `gamma(1.5)*gammaincc(1.5,u)`. Re-derived.
- **Matching.** τ_* = τ_m − a_m²/(2C) follows from a² = 2C(τ−τ_*) and C = a²H/H_i. Re-derived, and it is continuous
  at τ_m.
- **Clock class.** It applies PCHIP with 8-point Gauss per interval below τ_m and the closed-form tail above it, and
  finds quantiles with brentq on u ∈ [1e-8, 30].
- **Comparison with our own arithmetic.** Our code uses trapezoid integration on the dense Δτ = 0.01 a(τ) and
  numerical quadrature for the tail, with no PCHIP, Gauss or closed form. It matches the medians and median a exactly,
  the mean a per decay to ≤ 8e-15, the number fractions exactly and the energy fractions to ≤ 4.6e-7 relative.
- **Crossing time used.** The fractions are evaluated at the later bracket sample, which is the conservative choice.
  At the corrected τ = 56.46, the Y_L number fraction would be 5.68e-6, against 8.35e-6 at τ = 83.

## 4. Verdicts on the b73de64 claims

| Claim | Codex | st_rerun | Verdict |
|---|---|---|---|
| 0.03 rad unity by sampled τ | 134 (a 23.053) / 147.5 (a 20.878) | identical; dense 133.61 / 147.33 | CONFIRMED [C] |
| 1e-6 unity by sampled τ | 772 / 782 | identical; dense 771.95 / 781.82 | CONFIRMED [C] |
| 1e-9: no 0.1 or 1 crossing before 800 | none | none; max 1.6e-3 (amp), 2.6e-3 (charge) | CONFIRMED [C] |
| Local charge unity, 0.03 rad | 120 / 134.5 | identical; dense 119.58 / 134.42; dense-k 119.40 / 134.15 | CONFIRMED [C] |
| 10% warning, quiet3 0.03 rad | [82.5, 83.0] | **[56.45, 56.46]** | CORRECTED [C] |
| 10% warnings, other cases | 53.1–55.3, 563.5, 577 | 55.17, 563.47, 576.82 | CONFIRMED [C] |
| late C; matched τ_* | 1.986 / 1.476 | 1.98581 / 1.47587; 0.1937 / −0.1691 | CONFIRMED [C, D] |
| Medians Y_L / A(1e-6) / A(1e-7) | 6.89e6 / 4.36e17 / 4.36e20 | 6.8908e6 / 4.3552e17 / 4.3552e20 | CONFIRMED [C] |
| Y_L number and energy fractions at the 0.03 unity samples | 1.35–1.48e-5; 3.72–4.31e-8 | 1.3479–1.4837e-5; 3.7161–4.3067e-8 (own ≤ 4.6e-7 relative) | CONFIRMED [C] |
| Quartic μ/k, edges, second band | 0.440; 3.135–3.140; 5.57–5.98 | identical (1e-11) | CONFIRMED [C] (reproduction) |
| Q_cube at k = 2.13 | 1.34–2.09e5 | 1.3440e5–2.0937e5 | CONFIRMED [D] |
| Seven verification groups pass | pass | pass; numbers as in §2 | CONFIRMED [C] |
| Direct vector = full-transfer column; Ward kept | stated | diff + 6e-16 / 3.9e-16 checks | CONFIRMED [C, A] |
| Phase-current sign conversion | stated | exact negative, 4.1e-15 identity | CONFIRMED [D, C] |
| Conditional ordering: linear control lost before prescribed heating | stated | 772–782 ≪ 6.89e6 (number fraction ≤ 7.9e-5) | CONFIRMED as stated (conditional, linear) [A] |

## 5. New items (beyond the claims)

1. **Inherited 2cbc7f1 table.** [C, reproduction with Codex's integrator] Its "amplitude through τ80" for quiet (3, 0.5)
   is undersampled in time.
   - On the dense grid the 0.03-rad maximum is **20.2%** (τ = 69.61, k = 2.13), against the published 8.88%.
   - On the 18-k grid, the τ40 value is 7.42% against 6.97% (quiet3), and 8.24% against 5.44% (quiet5, k = 0.02).
     These τ40 comparisons are only approximate, because the 2cbc7f1 k grid differs.
   - quiet5 through τ80 (11.85% against 11.9%) and both local-charge values through τ80 (34.3% / 18.8% against
     34.5% / 18.8%) are unchanged.
   - The new packet's archive check compares only archived samples, so it cannot detect this.
2. **Growth rate.** [C] The response envelope (k ≈ 2) grows at d ln/dη ≈ 0.45–0.58 in conformal time over τ 150–800.
   This is of the order of, but below, the auxiliary quartic μ_peak of 0.632.
3. **The 1e-9 seed.** [E, extrapolation only, not computed] If the τ 400–800 envelope slope persisted, the seed would
   reach unity near τ ≈ 1.4–1.7×10³, still well before the Y_L median. Neither the persistence of the quartic regime
   nor the nonlinear evolution is established.

## 6. Limits

- These are reruns of Codex's code; we built no independent tangent solver.
- The dense-k scan covers k = 1.5–3.0 only.
- The dense sampling resolves crossings to 0.01 τ. It does not resolve the first continuous crossing in k-continuous
  space beyond a 0.02 k grid.
- The questions of physics scope remain as Codex states them and are not tested here: no metric, inflaton or radiation
  perturbations, no φ backreaction, no nonlinear release, no thermalization, and no production normalization.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
export GRUT_KINETIC_INPUT_ROOT=$SP/pinned/research/reality_framework/grut_formation_review_20261025
export GRUT_KINETIC_WORKDIR=$SP/kwork GRUT_ONSET_WORKDIR=$SP/owork
export GRUT_TIMING_WORKDIR=$SP/review/st_rerun/twork GRUT_MODE_WORKDIR=$SP/review/st_rerun/mwork
C=$SP/b73/research/reality_framework/grut_source_timing_20261026/code; T=tracks/st_rerun/code
cd $SP/b73 && for s in verify_response run_response quartic_scan source_clocks check_late summarize; do python3 -B $C/$s.py; done
python3 -I $T/compare_outputs.py $SP/review/st_rerun/published <regenerated dir> out.json
for j in time_b3_r10 time_b5_r10 time_b3_r12 time_b5_r12 kgrid_b3 kgrid_b5; do python3 -B $T/dense_crossings.py $C $SP/review/st_rerun/published/RESPONSE_RESULTS.json $SP/review/st_rerun/dense $j; done
python3 -B $T/ward_and_generator.py $C $SP/review/st_rerun/published ward.json
python3 -B $T/tolerance_scaling.py $C $SP/review/st_rerun/published $SP/review/st_rerun/dense tolerance.json
python3 -I $T/clocks_check.py $SP/review/st_rerun/published $SP/review/st_rerun/dense clocks.json
```

Results are in `results/`: reproduction_diff_vs_published.json, crossings_dense_vs_published.json,
ward_generator_charge_audit.json, tolerance_audit.json, clocks_own_arithmetic.json and st_rerun_summary.json.
