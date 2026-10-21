# Audit track "preparation": why sparse seed s101 selects one vacuum

21 October 2026. Claude audit of Codex `a69bd86`, `research/reality_framework/grut_expanding_charge_20261021`, issue #13.
Claude's packet: `c7b8480`, `modules/bias_threshold`.

Labels:
- [C] computed here;
- [E] estimate or model;
- [P] preparation input;
- [Q] quoted from Codex or Claude.

Nothing here validates, excludes or identifies Candidate 1.

## 0. Result

1. **Mechanism [C].** The sparse 64-mode preparations select a vacuum through a second-order S zero mode. It is
   sourced by the realization-level S-X cross-correlation, `chi0'' = (a^2 - F^2) chi0 - <dF^2 chi_S>_box`. It then
   outgrows every populated mode because the preparations have an infrared gap: k_min = 0.25-0.51, while the k = 0
   mode is tachyonic first, at a = 2^(1/3) = 1.26.
   - The quasi-linear model has no free parameters. Its inputs are the mode table, linear X response about the
     rotating orbit, and the lattice k-hat.
   - It reproduces Codex's recorded nonlinear `mean_S`, sign changes included, in all five sparse charged runs
     (s101, s202, s303, kappa4 s101, large s101) and in four interventions:
     - within 1.6% for 1.12 <= a <= 1.54 (3.8% at one point within 1e-6 of a sign change);
     - within 3% at a = 1.57;
     - s101 and large s101 within 0.4% through a = 1.60.
   - Gapfill matches to 4% at a = 1.54.
2. **Causality [C].** Interventions with Codex's unchanged binary, at 192^3 and Codex's main settings:
   - **flipX:** the X perturbation is negated and S is unchanged. The run selects the NEGATIVE vacuum: genuine minus
     0.9646, plus 0 at a = 2.2. That mirrors main_s101.
   - **noX:** X is unperturbed. Mean S stays <= 1e-9 before formation and <= 6e-6 through a = 2.2. The network is
     balanced to 3e-5 (m_end = 0.0000, both phases winding 7).
   - **gapfill:** one populated box-scale S mode is added, n = (1,0,0), power 1/64, with no X partner. Selection is
     removed: m_end = +0.014, two slabs, winding 6/6. The zero-mode source is unchanged.
   - Both ingredients are therefore necessary for s101's selection: the cross-correlation and the infrared gap.
3. **Seed ordering [C].**
   - s101: k_min = 0.477, WKB head start of k = 0 = 3.15 e-folds at rms 0.5, S0/sigma_f = 0.91 -> one-sign.
   - s303: k_min = 0.337, 1.61 e-folds, ratio 0.13 -> network biased +.
   - s202: k_min = 0.245, 0.82 e-folds, ratio -0.045 -> near-unbiased slabs.
   - kappa4: ratio -0.385, late formation at a = 1.96 -> network biased -. The sign flip relative to kappa = 2 is
     predicted.
   - Codex's initial cross term (-2.6e-6 for s101) has the WRONG sign as a predictor. In s101, 64% of S0 comes from
     the phase (w) component, which is invisible in the initial dF^2. For kappa4 the w component alone fixes the
     negative sign.
4. **Dense Gaussian [C].**
   - The cross-correlation is 10-80x smaller (C_uS = -0.001, -0.006, +0.006; ensemble sd 0.0059 at L = 76.8 and
     0.0032 at L = 115.2). The k = 0 head start over the lowest shell is 0.09 / 0.04 e-folds.
   - The formation imbalance equals the linear field's statistical imbalance:

     | Run | Predicted (linear field) | Observed (a = 1.60) | Ensemble sd |
     |---|---|---|---|
     | gaussian_s101 | -0.017 | -0.018 | 0.025 |
     | gaussian_s202 | -0.0076 | -0.0076 | 0.025 |
     | gaussian_large_s101 | -0.0099 | -0.0115 | 0.014 |

5. **The sparse family's bias does not shrink with the box [C].** Over 150 Codex-generator seeds:
   - |S0/sigma_f| >= 0.34 (the empirical one-sign class here) in 12.7% at L = 76.8 and 15.3% at L = 115.2.
   - sd(C_uS) = 0.076 at both box sizes.
   - large_s101 reuses s101's 64 coefficients exactly; its physical k vectors are re-rounded by <= 0.055. Its
     mean-S history equals main_s101's to 0.1% through a = 1.60.
   - So "s101 survives a larger box" is a property of the fixed-mode family, not finite-volume robustness.
6. **Forced box bias in a dense spectrum [C, reduced box].** 128^3 with dx = 0.4 (L = 51.2, Codex's generator,
   seed 101) and S noise + mu:
   - Linear m_f predictions are 0.22 / 0.43 / 0.66 / 0.80 / 0.94. The observed values at a = 1.60 are 0.27 / 0.49 /
     0.72 / 0.85 / 0.96.
   - The minority is removed by a = 2.2 for m_f >= 0.85 (mu >= 0.10). At m_f = 0.72 it ends at 1.3% volume with no
     winding. At m_f <= 0.49 the network is retained.
   - The gamma = 1.41 law with the unbiased reference K = 1.27 predicts a threshold m_f = 0.88. That has the right
     order but under-predicts growth in this small box: paired gamma_eff = 2.9-3.4. At mu = 0.02 the paired
     gamma_eff falls to 2.05 (L = 76.8) and 1.43 (L = 115.2) (§6b).
   - s101's sourced zero mode is equivalent to a forced mu of about 0.065 (ratio 0.91, against a dense ratio of
     14.1 per unit mu at rms 0.5). That is 6.6 sigma of the natural k = 0 fluctuation of a dense box at L = 76.8
     (sd 0.0099), and 17 sigma at L = 115.2 (sd 0.0054, large_s101 ratio 1.33) [E].
7. **Bias law (Claude) against these runs [C/E].**
   - Fed with the observed first-saturation bias and the dense reference K_ref = 1.48 (Codex's gaussian_s101), the
     takeover criterion erfinv(m_f) K^1.41 >= 1.54 (threshold m_f >= 0.78) sorts all 12 Codex runs correctly. The
     cold control is outside the law: it has no exit and is wall-free.
   - It also sorts the 9 forced-bias runs and flipX, noX and gapfill. It FAILS for the sparse interventions noU
     (m_f = 0.54) and noW (m_f = 0.28), which both lost the minority by a = 2.2.
   - With the sparse family's own unbiased reference (K_ref = 1.00-1.01: quasi-planar walls that do not coarsen), it
     would also mis-sort main_s101, dt_half, grid256, grid384 and flipX.
   - Paired responses (biased minus unbiased partner with identical fluctuations) give gamma_eff = 2.9-3.4 at
     L = 51.2, 2.05 at L = 76.8 and **1.43 at L = 115.2**. The last matches gamma = 1.41 over K = 1.62, but it is
     one seed.
   - The response ratio is nearly constant (2.0-2.2) while K_ref grows with L. So the excess at small L is a
     finite-box effect on the reference coarsening.
   - The sparse survivors give 2.6 (s303) and 3.9 (kappa4) in L = 76.8 boxes with box-coherent correlations.
   - This window (K <= 1.6) cannot test the asymptotic exponent or the bound gamma <= d/2.
8. **Cosmology.** Shown versus argued is in section 7. In short:
   - A periodic box's k = 0 mode is dynamically a systematic bias for that box. That is shown by the forced ladder
     and flipX.
   - Whether it is systematic physically depends on how its amplitude scales with the averaging volume. In the
     sparse family it does not scale. For a dense short-range spectrum it falls as (xi/L)^(3/2). That is computed
     for the ensemble, and the realized values at L = 76.8 and 115.2 (one seed each) are consistent with it.
   - For a thermal universe with horizon >> xi, a region-scale imbalance is a statistical fluctuation, so it cannot
     remove the network IF gamma <= d/2 holds asymptotically. That last step is argued (Claude's §5); it is not shown
     by any run here.

## 1. What was read and reproduced

- Codex files were extracted read only (`git archive a69bd86 ... | tar -x`) to scratch. They include run.py
  `modes()`, gaussian_run.py `noise()`, the expanding.cpp `initial()` and the kick (`(s*s-a*a+f)*s`), gaussian.cpp,
  RESULTS.json, RUN_TABLE.md, the protocols, and raw_diagnostics.tar.gz (all csv/topology series and the mode
  tables).
- Source hashes match SOURCE_FREEZE: expanding.cpp e62701e5..., gaussian.cpp 618668f6...
- Compiled exactly as run.py does (`g++ -std=c++17 -O3 -fopenmp -ffp-contract=off`).
- **Reproduction [C].** main_s101 at 192^3 with 2 threads in 216 s matches Codex's csv: |d mean_S| <= 5.2e-14, raw
  fractions identical, |dQ| = 2e-8 out of 2.9e5.
- **Mode generator replicated [C].** `codex_modes()` reproduces modes_L76.8_s{101,202,303} and modes_L115.2_s101
  bit for bit (max |diff| = 0).
- **Gaussian noise regenerated [C].** The covariances match Codex's records exactly (deviation 0.0; numpy 2.5.3 here
  vs Codex's 2.3.5).

## 2. Quasi-linear model (code/quasilinear.py) [E, validated against the records]

Conformal variables of expanding.cpp. Radiation a = 1 + tau/40, kappa = 2 (and 4), F^2 = kappa^(2/3).

- **S modes.** chi_k'' = -(khat^2 + F^2 - a^2) chi_k, with chi_k'(0) = chi_k(0)/40 (zero physical S velocity, as in
  `initial()`).
- **X modes about the exact orbit** U = F cos(w t), W = -F sin(w t), w = sqrt(0.1) F:
  - du'' = -khat^2 du - 0.1[F^2 du + 2U(U du + W dw)], and the same for w.
  - Initial momenta are the preparation's rigid rotation: du'(0) = w dw(0), dw'(0) = -w du(0).
  - dF^2 = 2(U du + W dw).
- **Zero mode.** Initially absent: <S> = 0 exactly on the grid for every preparation, since modes with
  0 < |n_i| < N sum to zero.
  - Equation: chi0'' = (a^2 - F^2) chi0 + J, with J = -<dF^2 chi_S>_box.
  - In the sparse case J = -eps_X F eps_S sum_m g_m (A_u,m C_uS,m + A_w,m C_wS,m), with C_fS,m = (c_f c_S + s_f s_S)/2
    per mode.
  - At t = 0, J reduces to Codex's `leading_mean_chiS_acceleration_cross_term` (-2.6247e-6 for s101, reproduced).
- **Cubic terms.** The bare -<chi_f^3> was computed from triads for the sparse case: s101 at L = 76.8 has no triads,
  and the ratio cubic/cross is < 0.005 for s202 and s303. For the Gaussian runs it was computed from the field.
  - The bare cubic over-predicts the Gaussian `mean_S` by about 6x at a = 1.54-1.57, because the S^2-driven X
    response is not modelled.
  - So for the dense runs only the sign/fluctuation statistics are used, not the cubic zero mode.
- **Validity.** Until total S rms is about 0.15-0.3 (a <= 1.57-1.60 at kappa = 2). After that the record saturates
  (e.g. s101: ratio 0.94 at a = 1.63, 0.54 at 1.66).

Comparison with Codex's nonlinear `mean_S` (results/sparse_diagnostics.json and intervention_quasilinear_*.json):

| Run | a range | Max |rec/pred - 1| | Example |
|---|---|---|---|
| main_s101 | 1.12-1.60 | 0.4% | a = 1.54: rec 1.263e-2, pred 1.262e-2 |
| main_s202 | 1.12-1.57 | 0.4% to a = 1.54; 3% at 1.57 | a = 1.54: -6.218e-3 / -6.200e-3 |
| main_s303 | 1.12-1.60 | 1.6% | a = 1.54: 4.081e-3 / 4.071e-3 |
| kappa4_s101 | 1.12-1.75 | <= 1.4%; 3.8% at a = 1.48, where |mean_S| < 1e-6 near a sign change | a = 1.75: -3.554e-4 / -3.552e-4 |
| large_s101 | 1.12-1.60 | 0.3% | a = 1.60: 1.083e-1 / 1.085e-1 |
| flipX / noU / noW | 1.12-1.57 | 0.3% / 0.6% / 1.0% | |
| noX | all | predicted 0; recorded <= 1e-9 before formation and <= 6e-6 after (roundoff/high order) | |

## 3. Initial-condition diagnostics per seed [C] (results/sparse_diagnostics.json, gaussian_diagnostics.json)

The predictions are for the linear field at total S rms = 0.5. m_lin is <sign S> including the sourced zero mode.
The lowest-mode share is that mode's fraction of the growth-weighted power.

| Run | k_min | S power k < 0.5 | C_uS | C_wS | Head start of k = 0 (e-folds) | a(rms 0.5) | S0 | sigma_f | S0/sigma_f | m_lin | Lowest-mode share |
|---|---|---|---|---|---|---|---|---|---|---|---|
| main_s101 (= dt_half, grid256, grid384) | 0.477 | 0.004 | +0.083 | -0.248 | 3.15 | 1.630 | +0.350 (u 0.127, w 0.223) | 0.383 | +0.91 | +0.52 | 0.89 |
| main_s202 | 0.245 | 0.023 | -0.063 | +0.082 | 0.82 | 1.578 | -0.023 | 0.513 | -0.045 | -0.021 | 0.94 |
| main_s303 | 0.337 | 0.074 | +0.023 | -0.085 | 1.61 | 1.616 | +0.066 | 0.502 | +0.13 | +0.073 | 0.19 |
| kappa4_s101 | 0.477 | 0.004 | +0.083 | -0.248 | 2.74 | 1.928 | -0.182 (u +0.022, w -0.204) | 0.474 | -0.385 | -0.21 | 0.82 |
| large_s101 (L = 115.2) | 0.506 | 0 | +0.083 | -0.248 | 3.53 | 1.634 | +0.406 | 0.306 | +1.33 | +0.80 | 0.81 |
| gaussian_s101 | 0.082 (dense) | 0.078 | -0.0010 | -0.0023 | 0.09 | 1.558 | ~1e-4 rec. | 0.514 | <= 0.01 | -0.017 (fluctuation) | 0.047 (lowest shell) |
| gaussian_s202 | 0.082 | 0.081 | -0.0059 | -0.0068 | 0.09 | 1.560 | ~1e-4 | 0.536 | <= 0.01 | -0.0076 | 0.064 |
| gaussian_large_s101 | 0.055 | 0.080 | +0.0059 | -0.0018 | 0.04 | 1.558 | ~1e-4 | 0.517 | <= 0.01 | -0.0099 | 0.027 |

- **Fluctuation-only sign mean.** For the sparse tables, <sign S_f> over the box is <= 2e-4. Sums of
  non-resonant plane waves have almost exactly symmetric sign fractions. The sparse family's statistical
  imbalance is therefore anomalously SMALL, while its sourced zero mode is anomalously LARGE.
- **Dense spectra.** The growth-weighted spectrum at formation peaks at k_eff = 0.21. 85% of its power is at
  k < 0.25, but only 3-6% is in the lowest box shell.
- **Mean S in the dense runs.** It stays <= 1e-4 until a = 1.48, and Codex records |mean_S| / sigma of about
  0.005-0.015 at a = 1.57-1.60.

## 4. Linear prediction of formation sign fractions, and outcome tests

The observed m is the raw-vacuum m at the first sample with f+ + f- >= 0.5, or at the formation sample.

| Run | m_lin (rms 0.5) | Observed m_f | Observed end (a = 2.2) | Outcome class |
|---|---|---|---|---|
| main_s101, dt_half, grid256, grid384 | +0.52 | +0.86 (a = 1.66) | f- <= 2e-6 | one-sign, no 10/10 formation |
| large_s101 | +0.80 | +0.997 | f- = 0 | one-sign |
| main_s202 | -0.021 | -0.030 | m = -0.034, K_basin = 1.00 (parallel slabs) | unbiased network |
| main_s303 | +0.073 | +0.12 to +0.16 | m = +0.50 | biased network |
| kappa4_s101 | -0.21 | -0.31 to -0.53 | m = -0.73, + phase not winding | biased network |
| gaussian_s101 | -0.017 | -0.018 | -0.094 | network |
| gaussian_s202 | -0.0076 | -0.0076 | +0.0008 (sign changed) | network |
| gaussian_large_s101 | -0.0099 | -0.0115 | -0.030 | network |
| cold_s101 | no tachyonic exit | | f+ = 1 | wall-free |
| int_flipX | -0.52 | -0.86 | f+ = 0 | one-sign, negative |
| int_noX | 0.000 | 0.000 | 0.0000 | exactly balanced network |
| int_gapfill | +0.012 | +0.013 | +0.014 | network (two slabs) |
| int_noU (w only) | +0.31 | +0.54 | f- = 0.0008 | one-sign by 2.2 |
| int_noW (u only) | +0.17 | +0.28 | f- = 0.0009 | one-sign by 2.2 |

- **Sign.** The predicted sign is correct in every run with |m_lin| > 0.005: 22 of 22, counting the dt/grid copies
  and the forced-bias runs.
- **Magnitude.** The linear magnitude under-predicts |m_f| by 1.0-1.2x (dense) to 1.4-2.5x (sparse). The extra
  nonlinear amplification is consistent with the k = 0 head start persisting past rms 0.5.
- **Pre-registered class rule.** The rule was |m_lin| >= 0.45 -> one-sign
  (results/intervention_predictions_preregistered.json, written before the runs).
  - It is right for flipX, noX and gapfill.
  - It is WRONG for noU and noW: predicted "biased two-sign", observed one-sign.
  - The empirical sparse threshold at kappa = 2 lies between S0/sigma_f = 0.13 (s303, network) and 0.34 (noW,
    removed).

## 5. Bias-law test (Claude's m = erf(y_f K^1.41), x10 = 1.54) (code/outcomes.py, results/outcomes.json)

K is the comoving V/A growth from Codex's basin interface estimator.
- **Unbiased references** (a = 1.6 -> 2.2):
  - Gaussian: K_ref = 1.48 (L = 76.8, s101), 1.40 (s202), 1.62 (L = 115.2), 1.27 (L = 51.2, this audit).
  - Sparse unbiased (s202, noX, gapfill): K_ref = 1.00-1.01. These are quasi-planar walls, so there is no
    coarsening in the window.
- **Takeover threshold within the window.** erfinv(m_f) K^1.41 >= 1.54 gives m_f >= 0.78 (K = 1.48), 0.72
  (K = 1.62) and 0.88 (K = 1.27).
- **Codex's 12 runs.** All are sorted correctly when fed the observed m_f with the dense K_ref = 1.48:
  - s101 x5: m_f 0.86-0.997 -> takeover, observed.
  - s303, kappa4, s202 and 3 Gaussian runs: no takeover, observed.
  - cold: no exit.
  - These are post-hoc inputs, not blind predictions. With the sparse K_ref of about 1.0 the threshold becomes
    0.97, and 4 of the 5 s101 runs would be mis-sorted.
- **Failures.**
  - noU (m_f = 0.54) and noW (m_f = 0.28) lost their minority by a = 2.2. With either reference K (1.48 or about
    1.0) the law predicts survival.
  - In biased survivors the law under-predicts m_end. For the dense runs gamma_eff is the paired response against
    the mu = 0 partner, with K_ref = 1.27:

    | Run | m_f -> observed m_end | Law prediction | gamma_eff |
    |---|---|---|---|
    | s303 | 0.17 -> 0.50 | 0.31 | 2.55 |
    | kappa4 | -0.31 -> -0.73 | -0.43 | 3.9 |
    | dense mu = 0.02, L = 51.2 | 0.27 -> 0.55 | 0.37 | 3.0 |
    | dense mu = 0.04, L = 51.2 | 0.49 -> 0.85 | 0.65 | 3.1 |
    | dense mu = 0.07, L = 51.2 | 0.72 -> 0.974 | 0.87 | 2.9 |

  - The paired larger-box results are in §6b.
- **Interpretation [E].**
  - These are early (K <= 1.6) and finite-box (L/xi_end = 3-6) growth rates. They exceed d/2 = 1.5, which the
    Fisher-Huse bound forbids only asymptotically under scaling with short-range initial correlations.
  - In the sparse family the initial correlations are box-coherent, so the bound does not apply at all.
  - The runs neither test nor refute gamma = 1.41 at the 10^11-10^17 extrapolation. They show three things:
    - the law under-estimates bias growth in boxes with L/xi of about 3-5;
    - in the largest dense box (L = 115.2) the paired response follows K^1.43 over this window (§6b);
    - the law does not describe the sparse family at all.

## 6. Forced-bias dense ladder (code/forced_bias.py; predictions pre-registered before the runs)

128^3 with dx = 0.4, L = 51.2, Codex's Gaussian generator (seed 101 at N = 128), Codex's unchanged gaussian.cpp,
dt = 0.06 and T = 48. The S noise is shifted by mu, giving a uniform S offset of 0.01 mu [P].

| mu | Pred. m_lin (rms 1, no bare cubic) | Obs. m(1.60) | Obs. m_end | f-_end | Basin area end / unbiased end | - winding at end |
|---|---|---|---|---|---|---|
| 0 | -0.003 | +0.007 | +0.061 | 0.441 | 1 | 7 |
| 0.02 | 0.215 | 0.265 | 0.549 | 0.214 | 0.72 | 6 |
| 0.04 | 0.426 | 0.490 | 0.851 | 0.072 | 0.40 | 4 |
| 0.07 | 0.655 | 0.717 | 0.974 | 0.0128 | 0.15 | 0 |
| 0.10 | 0.795 | 0.845 | 0.999 | 0.0005 | 0.03 | 0 |
| 0.15 | 0.941 | 0.964 | 0.9999 | 3e-5 | 0.002 | 0 |
| 0.25 | 1.000 | 1.000 | 0.9999 | 4e-5 | 0.002 | 0 |

- **Takeover threshold.** Within a = 2.2 in this box the threshold is m_f of about 0.72-0.85; the law predicts 0.88.
- **Natural k = 0 fluctuation.** A dense box of this spectrum has a natural k = 0 sd of 6.67/L^1.5: 0.018 at
  L = 51.2, 0.0099 at 76.8, 0.0054 at 115.2. Codex's generator removes it.
- **Its effect.** At L = 51.2 it corresponds to m_f of about 0.2, which is not enough for takeover by 2.2.

### 6b. Paired larger boxes (mu = 0.02; the unbiased partners are Codex's recorded gaussian_s101 / gaussian_large_s101)

This is the same noise as Codex's runs plus mu = 0.02, so the paired difference y_b - y_0 = erfinv(m_b) - erfinv(m_0)
isolates the forced bias from the realization's statistical imbalance (results/paired_boxes.json). K_ref is the
unbiased partner's basin-area growth from a = 1.6 to 2.2.

| L | K_ref | m_f biased (partner) | m_end biased (partner) | f-_end | Response ratio y(2.2)/y(1.6) | Paired gamma_eff | Law m_end (gamma 1.41) |
|---|---|---|---|---|---|---|---|
| 51.2 (N = 128, this audit's partner) | 1.27 | 0.265 (0.007) | 0.549 (0.061) | 0.214 | 2.05 | 2.97 | 0.37 |
| 76.8 (N = 192, Codex gaussian_s101) | 1.48 | 0.261 (-0.018) | 0.498 (-0.094) | 0.238 | 2.22 | 2.05 | 0.44 |
| 115.2 (N = 288, Codex gaussian_large_s101) | 1.62 | 0.261 (-0.011) | 0.488 (-0.030) | 0.243 | 1.99 | **1.43** | **0.49** |

- **Box dependence.** The response ratio is nearly box-independent (2.0-2.2), while the unbiased coarsening K_ref
  grows with L. So the excess gamma_eff in smaller boxes is a finite-box effect: the small box slows the reference
  coarsening, not the bias growth.
- **Largest box.** In the L = 115.2 box the paired response grows as K^1.43 over K = 1.62, matching Claude's
  gamma = 1.41 within this window. The law's m_end (0.49) equals the observed value (0.488). This is one seed and a
  short K range [C, E]. It is consistent with the law; it does not validate the 10^11-10^17 extrapolation.
- **Linear prediction.** The pre-registered m_lin (rms 1, no bare cubic) was 0.216 (L = 76.8) and 0.220
  (L = 115.2), against an observed m_f of 0.261 in both.
- **Takeover.** None of the three biased boxes takes over by a = 2.2 (minority 21-24%), as the law predicts for
  m_f = 0.26 (threshold 0.73-0.88).

## 7. Cosmological reading: what is shown and what is argued

**Shown here (numerically, in these boxes):**
- (S1) A box k = 0 S offset acts as a systematic bias for that periodic box. It is amplified by the k = 0 tachyonic
  rate, and above m_f of about 0.8 it removes the minority by a = 2.2 (dense ladder). A sourced zero mode acts the
  same way: flipX reverses the outcome, noX removes it.
- (S2) In Codex's sparse family the box-scale bias is sourced by a Z2-odd S-X cross-correlation of 64 global modes.
  - Its ensemble mean is zero: the signed mean of S0/sigma_f over 150 seeds is +0.008, and outcomes have both signs
    (s101 +, kappa4 -, flipX -).
  - Its realization amplitude does NOT fall with box volume: sd C_uS = 0.076 at L = 76.8 and 115.2, and the one-sign
    fraction is 12.7% -> 15.3%.
  - The infrared gap (k_min of about 0.25-0.5) gives it a head start of 1-3.5 e-folds.
  - So for this family it acts as a systematic bias at every box size, which is why the larger box did not restore
    two signs.
- (S3) For a dense Gaussian spectrum the same quantities scale as finite-volume fluctuations:
  - C_uS ensemble sd 0.0059 -> 0.0032 (proportional to L^-1.5);
  - the formation imbalance sd 0.025 -> 0.014, with realized values within 1 sigma;
  - the k = 0 head start over the lowest shell is 0.09 -> 0.04 e-folds.
  - The s101 bias equals a dense-box k = 0 fluctuation of 6.8 sigma (L = 76.8) or 17 sigma (L = 115.2).

**Argued, not shown:**
- (A1) In a universe whose horizon is >> xi, with a dense short-range (thermal, screened) spectrum, no coherent zero
  mode exists.
  - The analogue of a box mean over a region R is a random variable with zero ensemble mean and sd proportional to
    (xi/R)^(3/2). That covers both the S field itself (white-noise infrared tail) and the Z2-odd S-X correlation,
    which is local and averages the same way.
  - A box-scale infrared mode is therefore systematic for a box but statistical for the universe. It is a bias that
    varies from region to region, with no common sign.
  - The step from "statistical" to "cannot remove the network" is Claude's §5 argument. A region-R imbalance grows
    at most like (xi/R)^(3/2) K^gamma until domains reach R, i.e. to O((R/xi)^(gamma - 3/2)) <= O(1) if
    gamma <= d/2. It never exceeds the scaling network's own domain-scale imbalance.
  - That needs dynamical scaling and short-range initial correlations asymptotically. These runs do not show it.
    Their K <= 1.6 windows show box-limited excess growth at L/xi of about 3-5 (paired gamma_eff 2.0-3.4) and
    gamma_eff = 1.43 in the L = 115.2 box. That is consistent with gamma < d/2 in the largest box, but it neither
    confirms nor refutes the asymptotic bound.
- (A2) A cosmological analogue of the s101 selection would need one of two things.
  - A genuinely systematic Z2-odd ensemble asymmetry. That is Claude's m_f, whose size in any history is not computed.
  - Infrared-correlated (super-horizon) fluctuations with more than white-noise power at k -> 0, e.g. inflationary
    spectator modes. That is the spatially varying random-field bias Claude's module lists as not done.
  - A thermal dense spectrum supplies neither, but this argument is not a calculation of the hidden bath.
- (A3) Physical quenches are about 10^6.5-10^7.5 times slower: H/m_X of about 0.016 here vs 5e-10 to 5e-9 [Q Codex].
  - That widens the k = 0 head start for any FIXED k, but in a dense spectrum only modes with k < k_KZ compete, so a
    region's zero mode matters only if the region is smaller than the Kibble-Zurek length.
  - A horizon holds 10^13.7-10^15.7 KZ domains [Q Claude's bias_threshold], so this does not apply.
  - This is a scaling argument [E] and was not simulated.

## 8. Verdicts on claims in this track

- Codex: "s101 evolves toward almost a single positive vacuum, robust to dt halving, 256^3/384^3 and a larger box."
  **CONFIRMED** as an execution: reproduced main_s101, and the records agree. **CORRECTED** in interpretation:
  large_s101 uses the identical 64 coefficients (physical k re-rounded by <= 0.055), so it is not an independent
  larger-volume test. The selection probability of the family does not decrease with L (12.7% -> 15.3%).
- Codex: "finite mode sets can have appreciable S-X cross-correlations which may source a coherent S zero mode; its
  sign need not fix the eventual vacuum." **CONFIRMED and quantified.** The sourced zero mode is the measured
  mean S (<= 1% for a <= 1.57). The interventions establish causality. The initial cross term has the opposite sign
  to the s101 outcome because of the phase (w) component and the X response.
- Codex: "sparse finite global-mode ensemble is not a thermal ensemble; not a refutation of Claude's hot history."
  **CONFIRMED**, with mechanism: infrared gap plus a box-coherent Z2-odd cross-correlation.
- Codex: "Gaussian preparations form both vacua; box-winding through a = 2.2." **CONFIRMED** as consistent. Their
  formation imbalances are the statistical imbalances of the linear field (-0.017 / -0.0076 / -0.0099 predicted vs
  -0.018 / -0.0076 / -0.0115 observed).
- Claude: "systematic bias grows as m0 K^gamma, gamma = 1.41 in 3D." **NOT CHECKED** asymptotically.
  - Within the window (K <= 1.6, charged walls, dense forced bias) the paired gamma_eff is 1.43 in the L = 115.2 box
    (consistent, one seed) and 2.05 / 2.97 in the 76.8 / 51.2 boxes (finite-box excess).
  - It sorts Codex's 12 outcomes post hoc. It does not describe the sparse family: it fails for noU and noW, and with
    the sparse reference K it fails for the s101 runs too.
- Claude: "finite-volume statistical fluctuations of a short-range Gaussian field cannot remove a network because
  lambda_bar >= d/2." **NOT CHECKED** (argued). The runs show the prerequisite scaling of the statistical amplitude
  with L ((xi/L)^(3/2), dense) and no takeover in three dense boxes. They do not test the asymptotic bound.

## 9. Limits, failures, not done

- The quasi-linear model ignores the S^2-driven X response. It is validated only until total rms is about 0.15-0.3.
  The formation magnitudes beyond that are under-predicted by 1.2-2.5x.
- The pre-registered class rule failed for noU and noW. The gamma = 1.41 law failed for those two. It under-predicts
  the biased survivors in boxes with L <= 76.8 and matches only in the single L = 115.2 paired run.
- The reduced box L = 51.2 (N = 128) has only about 3 domains across at the end, so finite-size effects are strong.
  The paired larger boxes have one seed each (§6b). There is a single seed per mu.
- No large_s202 or large_s303 runs (Codex's unexecuted rows were not run here). No convergence study of the
  interventions: they use Codex's main 192^3, dt = 0.06, whose s101 outcome is unchanged under Codex's dt and grid
  refinements.
- The genuine-vacuum mask was not varied: F^2 < 0.5 as in Codex. The basin-area K includes the transient right
  after formation.
- No thermal bath, no D, prescribed expansion, and H/m_X about 10^7 too large. No cosmological persistence or
  survival probability follows.
- The 150-seed family scan uses the quasi-linear class threshold (0.34), which is empirical from three runs at
  kappa = 2. It is a fraction estimate [E], not run outcomes.

## 10. Files

- `code/quasilinear.py`: linear S growth, linear X response about the rotating orbit, the sourced zero mode, and the
  WKB exponent.
- `code/sparse_diagnostics.py`: per-seed diagnostics, comparison with Codex's mean_S, and sign fractions.
  Output: `results/sparse_diagnostics.json`.
- `code/gaussian_diagnostics.py`: regenerates Codex's dense noise and computes spectra, cross-correlations,
  predictions and the ensemble sd of m. Output: `results/gaussian_diagnostics.json`.
- `code/predict_interventions.py` -> `results/intervention_predictions_preregistered.json` (written before the runs).
  `code/interventions.py` runs Codex's binary. `results/intervention_quasilinear_*.json` holds the model against each
  record.
- `code/forced_bias.py` -> `results/forced_bias_predictions_preregistered.json` and
  `forced_bias_predictions_box{192,288}_preregistered.json`, plus the runs.
- `code/sparse_family_scan.py` -> `results/sparse_family_scan.json` (150 seeds x 2 boxes).
- `code/outcomes.py` -> `results/outcomes.json` (classes, K, gamma_eff for all runs).
- `results/rerun_series/`: csv and topology series of this audit's reruns (the inputs to outcomes.json).

**Reproduce.** Extract Codex's packet and the raw archive into scratch, compile expanding.cpp and gaussian.cpp as in
run.py, then run (OMP_NUM_THREADS = 1-2):

    python3 -I code/sparse_diagnostics.py <raw> results/sparse_diagnostics.json
    python3 -I code/gaussian_diagnostics.py <raw> results/gaussian_diagnostics.json
    python3 -I code/interventions.py <raw> <scratch>/runs <build>/expanding
    python3 -I code/forced_bias.py predict|run <scratch>/runs/fb ...
    python3 -I code/sparse_family_scan.py <raw> results/sparse_family_scan.json 150
    python3 -I code/outcomes.py <raw> <scratch>/runs results

Wall times with 2 threads: about 230 s per 192^3 run and about 75-110 s per 128^3 run.
