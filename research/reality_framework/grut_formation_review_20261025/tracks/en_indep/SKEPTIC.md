# Skeptic review of track en_indep (Codex 7d9c93b, early nonlinear X/φ pilot, issue #13)

Track key en_indep_skeptic. 27 October 2026. I tried to refute en_indep's findings by recomputing each key number with a
different code or method, or by rechecking the published data. Codex's packet was read only from a fresh git archive of
7d9c93b (`$SP/review/en_indep_skeptic/pub`). Nothing in Codex's tree or in other tracks was changed, and nothing was
committed. My code and results are in `skeptic/` (`code/`, `results/`).

**What each number comes from**

- **Independent code [C].** `code/linear_phi_moments.py` recomputes Codex's fixed-orbit linear φ control by its own
  method:
  - second-moment equations for each shell;
  - its own fixed-step RK4;
  - its own full-FFT binning.
- **Reproduction-type checks [C].**
  - `code/refcycle_codex.py` runs Codex's own `analyze.spectral()`, executed from the archived text, on Codex-engine
    fields.
  - The twin runs use Codex's `early3d.cpp`, copied verbatim and compiled with Codex's flags, on 2 threads. They restart
    from the read-only en_rerun N48 snapshots at τ = 100 and τ = 740. The twin design and analysis are my own
    (`code/snap2init.py`, `code/twin_analyze.py`).
  - None of these numbers uses en_indep's code or en_indep's fields.
- **Own analysis of Codex's published files [C].** `code/control_cycles.py` and `code/published_checks.py`, applied to
  FIELD_CURVES.npz, LINEAR_PHI_CONTROL.json, FIELD_RESULTS.json and RUNS.json.

## Verdicts on en_indep's findings

| # | en_indep finding | Skeptic verdict | Skeptic number vs en_indep / Codex |
|---|---|---|---|
| 1a | Charge contrast crosses at 60.60 (N48); controls 60.5–60.9 NOT CHECKED because "no control curves are published" | **N48 CONFIRMED. The reason for NOT CHECKED is wrong** | FIELD_CURVES.npz holds all seven histories: controls at Δτ = 0.5, halfstep48 at 0.05. Interpolated control crossings are 60.45–60.86; Codex's dense thresholds are 60.5–60.9. N48 is 60.5503 / 60.5503. |
| 1b | X variance reaches half at 69.25 | **CONFIRMED.** The controls are checkable. | Controls 69.21–69.28 from the Δτ = 0.5 grid; reported 69.2–69.25. Linear interpolation on that grid runs late by up to 0.04 on this convex rise. |
| 2a | 62.72176% outside the mean at τ 800 | **CONFIRMED** | Codex engine restarted from τ 740: 0.6272176034. Range over τ 600–800: 46.71–72.70%. Last-cycle mean: 60.74% (en_indep 60.7%). |
| 2b | Q₀/Q = 1.0546248 | **CONFIRMED** | Range 1.0467–1.1421. Last-cycle mean 1.0659 (en_indep 1.066). |
| 3a | φ + interaction 0.3098%; E_φ 0.2256% | **CONFIRMED** | Range over τ 600–800: 0.2694–0.3440%. The maximum after τ 200 is 0.34397%. |
| 3b | 2.09× linear φ control: NOT CHECKED | **Skeptic CONFIRMS Codex** (independent code) | Own result 2.091701 at h = 0.002 (Richardson about 2.09171) vs Codex 2.0917129. First unit crossing 702.5 / 702.5. Maximum 2.10943 / 2.10944. |
| 4 | Projections at τ 800 confirmed; "1.9-fold" CORRECTED to 1.04–2.51 | **The numbers are CONFIRMED. The CORRECTED label is overstated.** | Codex's spectral() on Codex-engine fields over τ 752–800 gives massless 1.2626–2.7475 (mean 1.978), Hartree 1.0923–1.2357 (mean 1.174) and a ratio of 1.0375–2.5136. The largest difference from en_indep is 1.4×10⁻⁷. Codex's τ = 800 ratio is exactly 1.896 as stated, and Codex never claimed it was phase-independent. The right verdict is CONFIRMED at τ 800 plus a NEW qualification. The one-period trapezoid over 740–788 gives massless 1.91, Hartree 1.175 and a mean ratio of 1.65. |
| 5 | NOT CHECKED; only the clock checked | **Agree** | C = a²H agrees bitwise because the arithmetic is the same. That is not an independent check of the matching. |
| 6 | Conservation CONFIRMED; work ledger 4.9×10⁻⁷ on τ ≤ 150 only | **CONFIRMED** | On the Codex-engine N48 trace the ledger residual peaks at 4.82×10⁻⁷ at τ 22.65 and is only 1.36×10⁻⁷ for τ ≥ 150. At half resolution it is 1.95×10⁻⁶ (4×, as expected for the second-order trapezoid rule). It is therefore limited by quadrature, and en_indep's τ ≤ 150 window contains the maximum. |
| N1 | No divergence between the codes; "no exponential growth is resolved" | **Headline CONFIRMED. The growth statement is CORRECTED to "slow growth is resolved".** | See N1 below. |
| N2 | The τ = 800 values sit on a large orbit oscillation, so control effects "cannot be separated from orbit phase" | **Oscillation CONFIRMED. The inference and its premise are REFUTED.** | See N2 below. |

### N1: slow growth is present but hidden by a test dominated by short wavelengths [C]

- **Codex's own published N48 and halfstep48 curves.** Their differences are flat from τ 40 to 800. For E_φ the
  difference is 6.4–6.5×10⁻⁶ in every segment. This reproduces en_indep's own-versus-Codex pattern without their code.
- **Why the differences are flat.** After τ ≈ 40 Codex's step is set by the sampling (Δs = 0.05/a, falling to
  8.9×10⁻⁴ at τ 800), not by hmax = 0.004. Truncation error therefore accumulates almost only before τ ≈ 40–60.
- **Twin runs on the Codex engine** from the same τ = 100 state, over Δη = 18.4. Growth is measured relative to the
  base run's fluctuation norm.

  | Perturbation | X | φ |
  |---|---|---|
  | IR band, \|k\| ≤ 2.5, size 10⁻⁸ | 1.0 → 5.3×10⁻⁸ (×5.3) | ×20.5 |
  | Mid band, size 10⁻⁸ | ×1.0 | ×4.2 |
  | White noise, size 10⁻⁶ (en_indep's test) | total ×0.71 (en_indep ×0.72); its IR part grows ×4.6, from 6.5×10⁻⁸ to 3.0×10⁻⁷ | ×2.2 (en_indep 2.2) |

- **Why en_indep's test missed it.** In the white-noise twin, 97% of the initial norm sits in short-wavelength modes
  that do not grow. That norm hides the long-wavelength growth.
- **Size of the growth.** The effective rates are at most 0.10 (X) and 0.15 (φ) per unit η. The orbit frequency is 7.0
  per unit η. The X growth fits a linear (secular) law slightly better than an exponential one.
- **So:** there is no rapid chaotic divergence, and the pointwise agreement between the codes stands. But "no
  exponential growth is resolved" should read "growth of up to ×5 (X) and ×20 (φ) over τ 100–800 is resolved".

### N2: control spreads at τ 800 are phase-matched [C]

- **Premise.** FIELD_CURVES.npz does contain all the control curves.
- **All seven histories reach τ 800 at the same orbit phase.** Each is 0.248–0.267 of a cycle past its last maximum of
  energy outside the mean, at τ 786.5–787.5. The period is 0.889–0.893 in η.
- **The single-instant spread and the cycle-mean spread agree:**

  | Diagnostic | τ = 800 spread | Last-cycle mean spread |
  |---|---|---|
  | Energy outside the mean | 60.94–64.93% | 60.28–63.09% |
  | Q₀/Q | 1.0543–1.1585 | 1.0655–1.1759 |
  | φ + interaction | 0.277–0.365% | 0.280–0.369% |
  | E_φ | 0.188–0.267% | 0.189–0.268% |

- **Ordering.** The extremes keep their places in both columns: box16 highest and seed2 lowest outside the mean. The φ
  ordering is the same throughout.
- **Real control effects.** The seed2 (+0.110) and box16 (+0.075) offsets in Q₀/Q persist in the cycle means, so they
  are not phase artefacts.
- **So:** Codex's τ = 800 control ranges are valid same-phase comparisons. They share a common bias against the cycle
  means: +2.0 points for energy outside the mean, −0.011 for Q₀/Q.

### Other

- **Stale environment warning.** en_indep's warning about the overwritten copy no longer applies:
  `$SP/e7d/.../MATCH_L12.npz` now has SHA-256 4462a2a3…, the published hash.

## New issues with the Codex packet that en_indep missed

1. **N32 is under-resolved and sets the low end of both neutral ranges.**
   - Its boundary-shell energy is 1.21%, against 0.18% for N48.
   - φ + interaction is 0.277% for N32 against 0.310–0.365% for the other six histories.
   - E_φ is 0.188% for N32 against 0.226–0.267%.
   - Going from N48 to N64 still changes E_φ by +1.8%.
2. **The orbit oscillation and the phase bias are not reported.** All τ = 800 values are taken at the same orbit phase.
   The packet does not state the ±13-point oscillation in energy outside the mean, or the cycle means. Its own
   FIELD_CURVES.npz contains everything needed to compute them.
3. **When the neutral forecast actually changes.**
   - Linear and nonlinear φ + interaction agree within 0.6% up to τ 40 and within 3% up to τ 60.
   - They separate from τ ≈ 70–90: 18% apart at τ 91, a mean ratio of 1.23 over τ 100–120 and of 2.38 over τ 150–200.
   - The forecast change is therefore established by τ ≈ 200.
   - The 2.09× value at τ 800 does not conserve energy: it is an extrapolation that exceeds the total active energy from
     τ 702.5. It is also not the review's 665–697 calculation.
4. **The late time step is set by the sampling, not by hmax** (see N1 above). Codex's packet describes the step as
   hmax = 0.004 without stating this.
5. **The reference projections also vary across controls.**
   - Hartree at τ 800: 1.174–1.266 across controls (seed2 1.212, box16 1.266). This is comparable to its own variation
     over one cycle, 1.08–1.24.
   - Massless at τ 800: 2.06–2.41 across controls, against 1.26–2.75 over one cycle.

## Limits

- I did not check claim 5 (the matching errors), the spectral-export binning, the N64 or box16 dynamics beyond their
  published curves, or the Y/Σ/gauge and quantum scope.
- The twin runs cover N48 only, from τ = 100.
- The growth rates are effective rates over a finite window, not Lyapunov exponents. I did not orthonormalise
  perturbation directions.
- The reference projections remain non-physical. Agreement between codes on the reduced Hamiltonian is not
  validation. No validation, exclusion or new-theory claim is made.
- Credit: classical-statistical lattice methods follow Khlebnikov–Tkachev, Micha–Tkachev and Epelbaum–Gelis–Wu
  (arXiv:1402.0115). The integrators are Yoshida (1990), used in Codex's engine, and Blanes & Moan (2002), used by
  en_indep.

## Reproduction

```bash
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
K=research/reality_framework/grut_formation_review_20261025/tracks/en_indep/skeptic
PUB=$SP/review/en_indep_skeptic/pub/research/reality_framework/grut_early_nonlinear_20261027   # git archive 7d9c93b
R=$SP/review/en_indep_skeptic/runs
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 -I $K/code/control_cycles.py $PUB/FIELD_CURVES.npz $K/results/control_cycles.json
python3 -I $K/code/published_checks.py $PUB $SP/review/en_rerun/ework/N48/state.tsv $K/results/published_checks.json
python3 -I $K/code/linear_phi_moments.py $SP/review/en_indep/ework/N48_initial.bin $K/results/linear_phi_control_own.json  # sha abe46838... = RAW_FILES
# Codex engine restarts: snap2init.py SNAP OUT KIND EPS SEED; run_engine.sh INIT PREFIX TAU0 A0 SAMPLE SNAP
#   cyc740: state_s0000036 none, sample .05, snap 4;  t100_{base,ir8,mid8,white6}: state_s0000004 {none,ir 1e-8,mid 1e-8,white 1e-6}, snap 50
python3 -I $K/code/refcycle_codex.py $PUB $R/cyc740 $K/results/refcycle_codex_engine.json $SP/review/en_rerun/ework/N48
python3 -I $K/code/twin_analyze.py $R/t100_base $R/t100_ir8 $R/t100_mid8 $R/t100_white6 $K/results/twin_growth_codex_engine.json
```
