# st_rerun_skeptic: adversarial check of the st_rerun findings on Codex b73de64 (source timing)

Track key `st_rerun_skeptic`, 27 October 2026. This note tries to refute the st_rerun verdicts listed in its NOTES.md
and `results/st_rerun_summary.json`. The subject is Codex's packet `research/reality_framework/grut_source_timing_20261026/`
at b73de64, read through the extracted tree `$SP/b73/research/reality_framework/` and
SP = /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad.

**Bottom line.** None of st_rerun's verdicts is refuted. Its crossing times, the corrected 10% warning, the clock
numbers and the quartic spectrum all come out the same with our own code, which shares no code with Codex's tangent
engine or with st_rerun's reruns of it. Several of st_rerun's verdicts rested only on reproduction, so this check makes
them stronger. We found four new items:

- The "unit modulation" times depend on how the threshold is defined.
- The gauge3 row has the same 10% sampling miss as quiet3.
- Codex's charge normalization is taken on a sparse time grid. The effect is 1e-5, which is negligible.
- The "upper bound" wording for the rest clock holds only with respect to time dilation.

None of these changes the conditional time ordering. All of this is same-model numerical work on supplied fields. It
validates nothing and excludes nothing.

## 1. What we built, and what is independent of what

Labels: [C] computed, [D] derived, [A] argued, [E] estimate, [Q] quoted.

- **Background [C], independent of Codex.** We used Claude's own reheating C engine,
  `tracks/rp_reheating/code/rpreheat.c` (mode 0). It reconstructs H from the Friedmann constraint, where Codex evolves
  H. Its preparation is our own `rp.prepare`. We ran it to τ = 800 at Δτ = 0.01
  (`skeptic/code/bg_own.py` → `skeptic/results/bg_own.json`). Constraint checks: Gauss ≤ 7e-14 and |H_aux/H − 1| ≤ 3.6e-12.
- **Tangent engine [C], independent formulation** (`skeptic/code/skeptic_tangent.c`, driver `run_tangent.py`):
  - It integrates physical cosmic-time perturbations of the six canonical reals:
    δφ̈ + 3Hδφ̇ + (k²/a² + V_rr + ξ_X R P_X)δφ = 0.
  - Codex instead uses conformal canonical (q, p), a Ward-projected Hessian and an a''/a term.
  - R = 6(Ḣ + 2H²) comes from the engine's Raychaudhuri Ḣ.
  - V_rr is a sympy Hessian of our potential, emitted as C. It matches a finite difference of `rp.potential`'s force
    to 4e-11.
  - The gauge term is dropped (e = 0). Its size, c|g|² ≈ e²|g|²/k², is ≤ 4e-7 at k ≥ 2.
  - The seed is a central difference of `rp.prepare(β, δ = ±1e-5)` with zero velocities.
    It gives δ|X|/|X|(0) = −0.0430432 per radian.
  - The charge is δN = a³(δφ₁φ̇₂ + φ₁δφ̇₂ − δφ₂φ̇₁ − φ₂δφ̇₁), normalized by our dense ⟨N_X⟩ over τ 64–80.
  - Integration is Dormand–Prince 5(4) with a **max-norm** error over all components, where Codex uses DOP853 with an RMS
    norm. The rtol is 1e-10, output is every 0.005 τ, and crossings are linearly interpolated between outputs.
  - The k sets are Codex's 18 values plus 0.005, 0.01, 2.04, 2.16, 2.18 and 2.2, and a band from 1.80 to 2.60 in steps
    of 0.01.
  - Each τ = 800 history takes 5–14 s.
- **Prototype cross-check.** A Python/LSODA version (`tangent_own.py`) uses Codex's Python background right-hand side.
  It agrees with the C engine to 5e-13 in the maximum 0.03-rad response through τ = 20. That cross-checks both
  backgrounds in the tangent context.
- **Floquet [C], our own method** (`skeptic/code/floquet_own.py`). The orbit is built from the invariants (E, L) of
  Codex's matched χ. The period comes from quadrature, and the monodromy from Radau starting at apoapsis. Edges are
  found by bisection.
- **Clocks [D, C]** (`skeptic/code/analyse.py`):
  - Widths come from our own formulas and couplings: m_L = √(M² − √2μh) = 22.2486, γ_L = 3κ²m_L/(256π³), and
    γ_A = e²(6eh)/(48π).
  - a(τ) comes from our engine, with a closed-form Γ(3/2, u) tail.
  - Our γ_L = 1.0059033e-7 equals the pinned value.

**Direct agreement with Codex's published arrays** (quiet3.npz and quiet5.npz, at Codex's sample times):

- The amplitude per radian agrees to 4.8e-9–5.1e-7 of the running envelope for k = 0.02–8.
- The charge fraction differs by a constant 1.123e-5 (quiet3) and 1.17e-6 (quiet5). This is entirely normalization: see
  new item 2.
- Over τ = 0–800, our rtol 1e-12 run differs from our 1e-10 run by 3e-11 at k = 2.13. Crossings move by less than 1e-11.

## 2. Claim-by-claim result (st_rerun value → our independent value)

| id | st_rerun verdict | Our recomputation | Refuted? |
|---|---|---|---|
| R1 | CONFIRMED (repro) | Not rerun. We did check the extracted packet: all 22 output files and 10 code files are hash-identical to `git show b73de64`, so the restoration holds. Our independent engine agrees with the published arrays (§1). | no |
| R2 | npz bitwise | The same hash check. We made no regeneration of our own. | no |
| R3 | dense 133.61 / 147.33 | 133.601 (k 2.13, a 23.018) / 147.320 (k 2.0, a 20.865). On the band: 133.592 (k 2.19) / 147.312 (k 2.05). | no |
| R4 | 771.95 / 781.82 | 771.942 / 781.810. On the band: 771.933 / 781.792 (k 1.97). | no |
| R5 | max 1.6e-3 / 2.6e-3 | quiet3 1.585e-3 / 2.588e-3; quiet5 1.441e-3 / 1.944e-3 | no |
| R6 | 119.58 / 134.42; 737.43 / 751.08; dense-k ≤ 0.28 τ earlier | 119.579 / 134.414; 737.422 / 751.071. On the 0.01 band: 119.393 / 134.148 / 737.259 / **750.758**, which is up to **0.31 τ** earlier (st_rerun: ≤ 0.28). | no (minor) |
| R7 | CORRECTED → 56.46 (k 2.13) | **56.457** (k 2.13); band 56.456 (k 2.15). The excursion runs 56.46–57.00 (0.535 τ). The Y_L number fraction there is 5.679e-6. | no |
| R8 | 55.17; 563.47; 576.82 | 55.167 (k 0.02; k = 0.005 gives 55.166); 563.468; 576.811 | no |
| R9 | C 1.98581 / 1.47587; τ* 0.19367 / −0.16906 | Our engine: 1.985805275 / 1.475867788; τ* 0.193674 / −0.169062. a and H agree with Codex to ≤ 1.2e-10 and 4.7e-10 over τ 80–800. | no |
| R10 | medians as Codex | Y_L 6.890794e6 (a 5231.40 / 4509.97); A 4.355172e17 (1.31518e9 / 1.13381e9); 4.355172e20 (4.15897e10 / 3.58543e10) | no |
| R11 | 1.3479–1.4837e-5; 3.7161–4.3067e-8 | The same to 5 digits. At our first crossings (133.60 / 147.32): 1.3439e-5 / 1.4819e-5. Mean a per decay: 5568.653 / 4800.711. | no |
| R12 | CONFIRMED (reproduction only) | **Independent:** μ/k = 0.439798 / 0.439789 / 0.439722 at k = 0.01 / 0.02 / 0.05. Peak 0.631777 at k = 2.13. Edge at 3.13949. Second band 5.57118–5.97645, peak 0.16027 at 5.775. Period 0.81403767, matching Codex. | no |
| R13 | 1.3440e5–2.0937e5 | 2.09371e5 / 1.34403e5 with our dense ⟨N_X⟩ (Codex 209368.9; the 1.1e-5 difference is normalization) | no |
| R14 | 7 groups pass (repro) | Not rerun. These are same-author checks, and our engine is the external check. | no |
| A1 | RHS = full-transfer column | Our different formulation (no Ward projection, cosmic time) reproduces the published responses to ≤ 5e-7 of the envelope, so the equations are equivalent. | no |
| A2 | tolerance adequate | Our 1e-12 against 1e-10 differs by 3e-11. Max-norm control gives the same crossings. | no |
| A3 | δN is the covariant X charge; sign | With e = 0, the φ̇-based δN equals Codex's canonical-p δN up to the normalization constant: 1.1232e-5 = 8156.7346/8156.6430 − 1. | no |
| A4 | crossings() ok; spacing 1.89–3.22 | Logic ok. Actual spacing over τ 40–80 is **1.49**–3.22: the inserted τ = 40 knot gives 1.49. Unity excursions last 1.37 / 1.10 τ (st_rerun: 1.36 / 1.08). | no (trivial) |
| A5 | clocks re-derived | We re-derived the tail and τ*. Our own background and widths give the same numbers. | no |
| N1 | 2cbc7f1 8.88% → 20.2% | **20.18%** at τ 69.61, k 2.13. The band gives 20.36% (k 2.19). quiet5 11.85%. Charge 34.26% / 18.78% on 18 k, and 34.57% / 18.87% on the band. | no |
| N2 | 0.45–0.58 per η (said k≈2) | Max over k: 0.554 / 0.478 / 0.546 (quiet3) and 0.535 / 0.572 / 0.448 (quiet5) over τ 150–400 / 400–800 / 150–800. k = 2.13 alone gives 0.433 for quiet5 over 150–800. st_rerun's code takes the max over all k, although the note says k ≈ 2. | no (label) |
| N3 | [E] ~1429 / 1704 | [E] 1429 / 1702 with the same assumption | no |

## 3. New items with the Codex packet that st_rerun missed

1. **The unit-threshold times are set by periapsis spikes; the result is definition-sensitive [C, A].**
   - δ|X|/|X| peaks where the X orbit passes its radial minimum. There |X| is only 0.11–0.12 of the local orbit
     maximum at τ 100–150, and 0.22–0.24 at τ 700–800.
   - We also measured δ|X| against the local orbit maximum, max |X| within ±8 τ (`orbit_scale.json`):
     - The 0.03-rad seed reaches unity at τ = **196.3 / 234.0**, not 133.6 / 147.3. It reaches 0.3 at 153.3 / 188.0.
     - The 1e-6 seed reaches only **0.41 / 0.34** by τ = 800, and 0.3 at 772.3 / 782.4.
   - So "a 1e-6 seed reaches order-one modulation by 772–782" is true for Codex's fractional-|X| diagnostic. At that time
     the field-scale perturbation is about 0.3 of the orbit.
   - MATHEMATICS.md says the solution has "certainly lost small-perturbation control" at fractional unity. That is a
     statement about |X| at periapsis, not about the Cartesian field dynamics, whose nonlinearity scale is the orbit
     amplitude.
   - The ordering against Y_L (median 6.89e6) is unaffected. [E] At the measured growth, orbit-scale unity for 1e-6 would
     come at τ ≈ 900 (extrapolation; not computed).
2. **The gauge3 row has the same 10% sampling miss [C].** The gauge3 0.03-rad 10% bracket [82.5, 83.0] (k 2.13) is
   time-undersampled in the same way as quiet3. Our engine on {0.02, 2.13, 5.8} gives 56.46. Changing e from 1e-7 to
   1e-6 moves the published responses by ≤ 1.1e-6, and our e = 0 run agrees with e = 1e-7 to ≤ 5e-7. The R7 correction
   therefore applies to gauge3 as well. The gauge5 brackets (614.84, 782.44 and 759.42 by our dense run) are fine.
3. **Codex's charge normalization uses a sparse grid [C].** `late_NX_normalization_through80` is a trapezoid over the 6
   samples in [64, 80]. A dense time mean differs by 1.12e-5 (quiet3) and 1.1e-6 (quiet5). This shifts every
   charge fraction, and Q_cube (209368.9 against 209371.3), by that factor, which is negligible.
4. **"Upper bound" holds only against time dilation [A].** REPORT.md calls the early rest-clock fraction an upper bound.
   That holds only with respect to time dilation and with the vacuum width. In-medium widths are listed as uncomputed,
   and with the very large X occupation, stimulated Y_L → XXX is not excluded. To bring Y_L decays to order unity by the
   linear crossings, the width would need an enhancement of about 7e4 (0.03 rad, τ 134) or 1.3e4 (1e-6, τ 772).
   Nothing here computes such an effect.
5. **Gauge comparison.** At these e the gauge comparison is close to a null test: c|g|² ≲ 1e-6 at k ≥ 2. "Changing e to
   1e-6 has little effect" is true, but it tests little of the gauge physics. Codex already says so in part.

## 4. Limits

- The background right-hand side and the tangent equations are ours. They still encode the same supplied model:
  parameters [Q] from rp.P, the same as Codex's.
- Neither our engine nor Codex's includes metric, inflaton or radiation perturbations, φ backreaction, the gauge
  dynamics beyond e = 0, nonlinear release or thermalization.
- The fine k band covers 1.80–2.60 only; elsewhere only the 18 published k plus 0.005, 0.01 and 2.04–2.2 were run.
- We did not rerun Codex's verify_response.py or check_late.py. R1, R2 and R14 rest on st_rerun's reproduction, plus
  our hash check and the independent agreement in §1.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad; L=$SP/review/st_rerun_skeptic
K=research/reality_framework/grut_formation_review_20261025/tracks/st_rerun/skeptic
CX=$SP/b73/research/reality_framework/grut_source_timing_20261026
python3 -B $K/code/bg_own.py $CX $L $K/results/bg_own.json
for h in quiet3 quiet5; do for s in main band; do python3 -B $K/code/run_tangent.py $L $h 800 $s 1e-10 0.005; done
  python3 -B $K/code/run_tangent.py $L $h 800 k213 1e-12 0.005; done
python3 -I $K/code/floquet_own.py $CX/QUARTIC_SCAN.json $K/results/floquet_own.json
python3 -B $K/code/analyse.py $L $CX $K/results/analysis_own.json
python3 -B $K/code/orbit_scale.py $L $K/results/orbit_scale.json
```

The results are in `skeptic/results/`, with the summary in `skeptic/results/skeptic_summary.json`.
