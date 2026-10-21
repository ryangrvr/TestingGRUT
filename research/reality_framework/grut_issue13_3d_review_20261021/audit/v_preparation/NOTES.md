# Verifier of audit track "preparation" (v_preparation)

21 October 2026. Adversarial check of `audit/preparation/` (NOTES.md, code, results) in Claude's audit of Codex `a69bd86`,
`research/reality_framework/grut_expanding_charge_20261021`, issue #13.

Labels:
- [C] computed here;
- [E] estimate or model;
- [P] preparation input;
- [Q] quoted.

Nothing here validates, excludes or identifies Candidate 1.

## 0. Verdicts

| # | Claim of the track | Verdict | Why (independent value) |
|---|---|---|---|
| 1 | s101 → almost a single positive vacuum. Reproduced. Robust to dt/2, 256^3, 384^3 | CONFIRMED | The audit's repro csv vs Codex: max \|d mean_S\| = 5.15e-14, fractions bit-identical. Codex refinements: end f+ 0.9777-0.9789, f- <= 2e-6; mean_S(1.54) within 0.08% of main. These test discretization of ONE polynomial, not the preparation. |
| 2 | large_s101 is the same preparation, not an independent larger volume. The one-sign share does not fall with L. | CONFIRMED (two supporting numbers corrected) | Coefficients identical (max diff 0.0). Physical k re-rounded by up to **0.077** (median 0.047), not <= 0.055; \|d\|k\|\| up to 0.065. mean_S tracks main_s101 within **0.1-0.2% for 1.39 <= a <= 1.60 (0.5% at 1.30)**, not 0.06%. Family scan (own generator transcription, bit-identical tables): 12.7% / 15.3%, sd C_uS 0.076 / 0.076, class agreement 94.7%. But only 68% of seeds reuse identical coefficients across L, and 12.7 vs 15.3 is within binomial noise (sd about 2.8 points). The L-independence is structural: always 64 unit-power modes. |
| 3 | The parameter-free quasi-linear model reproduces the nonlinear mean_S, sign changes included. 64% of s101's S0 is from w. The kappa4 sign flip is predicted. | CONFIRMED | Re-derived from expanding.cpp with the solver's own KDK leapfrog (dt = 0.06, lattice khat, leapfrog background orbit). This is a different integrator from the audit's DOP853 code. Max \|rec/pred - 1\|:<br>- s101 0.22% (1.12-1.60);<br>- large_s101 0.18%;<br>- s303 1.4%;<br>- kappa4 1.3%;<br>- s202 0.6% to 1.54, 3% at 1.57.<br>J0 = -2.6247e-6. s101 S0 (rms 0.5) = 0.349 = u 0.126 + w 0.223 (64% w). kappa4: u +0.023, w -0.212. kappa4 a = 1.75: rec -3.554e-4 vs -3.550e-4. **New pre-registered dose test** (X scaled by 0.1 and 0.3, S unchanged): mean_S scales exactly (a = 1.54: 1.2627e-3 / 3.788e-3 rec vs 1.2607e-3 / 3.782e-3 pred). |
| 4 | Selection needs both the S-X cross-correlation and the infrared gap (flipX / noX / gapfill) | CONFIRMED, with a caveat on what gapfill isolates | flipX was rerun from an independently built table and is **bit-identical** to the audit's record: f- 0.978047, genuine - 0.964622, + 1e-6, winding 0/7. noX record: f+ - f- about 1e-5. There is no exact lattice antisymmetry: no translation d with n.d = N/2 mod N exists for s101/s202/s303, so the balance is dynamical, not imposed. Pre-registered dose tests on the same S field: Xscale0.1 (ratio 0.096) keeps a + biased network, winding 7/7, m_end +0.196. Xscale0.3 (ratio 0.285) removes the minority by a = 2.2. Both outcomes were as predicted. Caveat: gapfill adds a first-order mode of average power (1/64). It shows that a populated near-k = 0 mode swamps a second-order zero mode. It does not measure how much gap is needed. |
| 5 | Seed ordering is predicted by initial-condition diagnostics. Sign right in 22/22 runs with \|m_lin\| > 0.005. | CORRECTED | S0/sigma_f and m_lin are reproduced: 0.912/0.524, -0.045/-0.021, 0.131/0.073, -0.387/-0.211, 1.332/0.804. (a) The quoted head starts (3.15 / 1.61 / 0.82) are WKB-from-crossing. The actual leapfrog amplitude ratio, k = 0 vs the best of the 6 lowest modes, is 3.00 / 2.69 / 1.06 (s101 / s303 / s202). So s101 vs s303 differ by 0.3 e-fold, not 1.5. That contrast is mainly the 3x larger \|C\|, not the gap. (b) By the track's own stated observable (first sample with f+ + f- >= 0.5), gaussian_s101 has m_f = **+0.0046** (a = 1.57, the track's outcomes.json) against m_lin = -0.017, a sign miss. The NOTES table uses a = 1.60 (-0.018) instead. A like-for-like linear m05 at a = 1.57 is +0.0058 and does match. (c) The dt/grid/repro/large copies are not independent sign tests. Dense signs sit at 0.2-0.8 ensemble sd, so they are weak tests. |
| 6 | Dense Gaussian imbalance is finite-volume statistics, not a sourced bias | CONFIRMED | Noise regenerated (covariance deviation 0.0 for all three). Linear-field m05 at a = 1.60 is -0.0166 / -0.0079 / -0.0126 vs recorded -0.0180 / -0.0076 / -0.0115. At a = 1.57: +0.0058 / -0.0121 / -0.0119 vs +0.0046 / -0.0115 / -0.0096. Arcsine-law ensemble sd: 0.025 / 0.025 / 0.014. C_uS (realized covariance) -0.0010 / -0.0059 / +0.0059. Analytic ensemble sd (2 pi)^(3/4)/L^(3/2) = 0.0059 / 0.0032. WKB k = 0 vs lowest shell: 0.090 / 0.040 e-folds. |
| 7 | Bias law in window: paired gamma_eff = 1.43 at L = 115.2 "matches gamma = 1.41", 2.05 at 76.8, 2.97 at 51.2. The excess at small L is finite-box. | CORRECTED | Reproduced at the track's start a0 = 1.60 (1.425 / 2.05 / 2.97). But gamma_eff depends on a0. At L = 115.2 it is **0.95 / 1.43 / 1.67 / 1.81** for a0 = 1.57 / 1.60 / 1.72 / 1.84. The biased m dips just after formation (0.33 -> 0.26) before growing. So "matches 1.41" and "consistent with gamma < d/2 in the largest box" are window choices. For a0 >= 1.72, gamma_eff exceeds d/2 = 1.5 even at L = 115.2. Robust parts: the ordering 51.2 > 76.8 > 115.2 at every a0, and a response ratio of about 2.0-2.2 at a0 = 1.60-1.72. The asymptotic gamma = 1.41 is NOT CHECKED (as the track says). |
| 8 | Cosmological reading (A1-A3) | NOT CHECKED | The track itself labels it argued, not shown. Nothing here tests it. |

## 1. What was done

- Codex's packet was extracted read-only (`git archive a69bd86`) to scratch. expanding.cpp sha256 e62701e5... matches
  SOURCE_FREEZE. Compiled as run.py does.
- **Mode tables** (code/mode_tables.py -> results/mode_tables.json) [C].
  - Coefficient identity across L.
  - k re-rounding.
  - An exhaustive search over all 192^3 translations for an exact S antisymmetry, n.d = 96 mod 192 for all 64 modes:
    none for s101/s202/s303.
- **Leapfrog quasi-linear model** (code/leapfrog_ql.py -> results/leapfrog_ql.json) [E, validated against Codex's
  records].
  - Written from expanding.cpp `initial()` and `kick()`: du = 1e-3 F f_u, dw = 1e-3 F f_w, S = 1e-2 f_S; p_u = omega q_w,
    p_w = -omega q_u, p_S = q_S/40.
  - Everything (orbit, X perturbations, S modes, zero mode with kick (a^2 - \|X_bg\|^2) chi0 - <dF^2 chi_S>) is advanced
    with the solver's KDK map at dt = 0.06.
  - The u/w split is by linearity.
  - Its `head_start_efolds_k0_vs_kmin` field uses the signed k = 0 amplitude (NaN for kappa4). Use
    results/head_start.json instead.
- **Head start** (code/head_start.py -> results/head_start.json). WKB-from-crossing (the track's measure, reproduced to
  0.01) vs the actual linear amplitude ratio under the preparation's initial velocity.
  - s303's lowest mode (khat 0.337) crosses with an unfavourable phase. Its third-lowest mode is the effective
    infrared competitor, at -2.69 e-folds vs k = 0.
  - kappa4: WKB 2.75, amplitude 1.95.
- **Family scan** (code/family_scan.py -> results/family_scan.json). Own transcription of run.py `modes()`, bit-identical
  to Codex's four tables. 150 seeds x 2 boxes, leapfrog S0/sigma_f at rms 0.5.
- **Dense check** (code/dense_check.py -> results/dense_check.json).
  - Codex's noise is regenerated.
  - S is evolved linearly per lattice khat with the leapfrog.
  - The sign statistics are taken at the recorded samples.
  - The ensemble sd comes from the arcsine law with the ensemble spectrum exp(-k^2) g(khat)^2.
- **Paired exponent** (code/paired_check.py -> results/paired_check.json). Recomputed from the track's fb_* records and
  Codex's gaussian records. Start sample varied.
- **Codex-binary runs** (code/run_variants.py, 192^3, Codex main settings, 2 threads, about 190-230 s each). Predictions
  were written before the runs (code/predict_dose.py -> results/dose_predictions_preregistered.json, 16:51 UTC; the runs
  were launched after).

  | Run | Pre-registered | Observed | Result |
  |---|---|---|---|
  | flipX (rerun) | bit-identical to the track's int_flipX | all csv columns diff 0, topology identical | pass |
  | Xscale0.1 (X x 0.1) | S0/sigma_f = 0.096, m_lin 0.047. Network retained, + bias, m_f in [0.066, 0.117] | mean_S(1.54) = 1.2627e-3 (pred 1.2607e-3). First f+ + f- >= 0.5 at a = 1.66 with m = +0.0725. End m = +0.196, f- = 0.351, both phases winding 7 | pass |
  | Xscale0.3 (X x 0.3) | S0/sigma_f = 0.285 (the track's untested band 0.13-0.34). Guess: minority removed by 2.2, LOW confidence | mean_S(1.54) = 3.788e-3 (pred 3.782e-3, +0.16%). First saturation a = 1.66 with m = +0.229 (in the band [0.198, 0.353]). Minority 0.26-0.36 until a = 1.87, then removed: f- = 1.0e-4 at 2.2, minus winding 0 | pass (guess right) |

**Threshold refined [C].** The kappa = 2 sparse one-sign threshold in S0/sigma_f(rms 0.5) lies in (0.13, 0.285]:
- retained: Xscale0.1 (0.096), s303 (0.131);
- removed: Xscale0.3 (0.285), noW (0.34).

The track had (0.13, 0.34]. With 0.285 the family fraction above threshold is 13.3% (L = 76.8) and 18.0% (L = 115.2).
With 0.13 it is 28.7% / 32.0%. The bracket is set by one preparation's scalings plus three other runs. The new runs
remove the minority at m_f = 0.23. The track's dense ladder retains it at m_f = 0.49 (L = 51.2), so the sparse family's
bias growth is far stronger than Claude's law with the sparse K_ref of about 1.0 [C].

## 2. Errors and overclaims found in the track

1. Physical-k re-rounding bound: 0.077, not <= 0.055.
2. The large_s101 mean_S agreement is 0.1-0.2% (0.5% at a = 1.30), not 0.06%.
3. **Head-start attribution.** The WKB-from-crossing e-folds overstate the s101/s303 infrared contrast by 5x (1.54 vs
   0.31 e-fold against the best low mode). The seed ordering follows S0/sigma_f, which is reproduced. The headline
   attribution to k_min should be weakened: \|C_wS\| is 0.248 vs 0.085, and S0 at rms 0.5 is 0.349 vs 0.066.
4. **Observable switching.** NOTES §4 defines the observed m_f at the first sample with f+ + f- >= 0.5. For the Gaussian
   rows it quotes a = 1.60 values. With the stated definition gaussian_s101 is a sign miss (+0.0046 vs m_lin -0.017).
   gapfill's first-saturation m is +0.046 (a = 1.54), not +0.013. The count "22 of 22" cannot be reproduced as stated:
   the listed qualifying runs number 23-24.
5. **Window choice in the paired exponent** (see verdict 7). The statement "consistent with gamma < d/2 in the largest
   box" (§7 A1) does not survive a later start.
6. The dense-equivalent mu of s101 (0.065, "6.6 sigma") depends on the exact a used for the dense ratio: 14.7 per unit
   mu at a = 1.558 and 17.3 at a = 1.5625, against the track's 14.1. That is about ±15%; the order of magnitude stands.

## 3. Limits of this verification

- main_s101, gapfill, noU, noW, noX and the forced-bias ladder were not rerun. Their records were read, and flipX's
  bit-identity shows that the track's pipeline reproduces.
- One preparation (s101) carries every intervention, and the dose test adds two scalings of the same preparation. No
  seed other than s101 was intervened on.
- The quasi-linear model (both versions) is second order. It is valid until total rms is about 0.3 (a <= 1.60). The
  formation magnitude m_f is under-predicted by 1.4-2.5x in the sparse runs (Xscale0.1: 0.047 -> 0.073).
- The forced-bias ladder runs at L = 51.2 and the paired boxes have one seed each. gamma_eff spans more than 0.8 over
  plausible start samples, so these windows cannot discriminate gamma = 1.41 from d/2 or above.
- No thermal bath, no D, prescribed expansion, H/m_X about 10^7 too large [Q Codex]. No cosmological persistence or
  survival statement follows.

## 4. Reproduce

    git archive a69bd86 research/reality_framework/grut_expanding_charge_20261021 | tar -x -C <scr>/codex
    tar -xzf <scr>/codex/.../raw_diagnostics.tar.gz -C <scr>/raw       # -> <raw> = <scr>/raw/raw
    g++ -std=c++17 -O3 -fopenmp -ffp-contract=off expanding.cpp -o <scr>/build/expanding
    python3 -I code/mode_tables.py <raw> results/mode_tables.json
    python3 -I code/leapfrog_ql.py <raw> results/leapfrog_ql.json            # about 10 s
    python3 -I code/head_start.py <raw> results/head_start.json
    python3 -I code/family_scan.py <raw> results/family_scan.json 150        # about 4 min
    python3 -I code/dense_check.py <raw> results/dense_check.json            # about 3 min
    python3 -I code/paired_check.py <raw> ../preparation/results/rerun_series results/paired_check.json
    python3 -I code/predict_dose.py <raw> results/dose_predictions_preregistered.json
    OMP_NUM_THREADS=2 python3 -I code/run_variants.py <raw> <scr>/runs <scr>/build/expanding Xscale0.1 flipX Xscale0.3
    python3 -I code/summarize_variants.py <scr>/runs ../preparation/results/rerun_series results/variant_runs.json
