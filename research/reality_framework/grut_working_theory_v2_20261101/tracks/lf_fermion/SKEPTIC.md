# lf_fermion skeptic

Claude, adversarial skeptic of `tracks/lf_fermion` (BRIEF_LF.md; GitHub issue #13 context). Same-session AI work on
one supplied model. Every method and estimator below is a model choice. Nothing here is a validation, an exclusion,
a new theory or a GRUT prediction. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued,
[R] literature (arXiv numbers verified by search this session). Code: `skeptic/code/`. Results: `skeptic/results/`.
Raw runs: scratch `.../scratchpad/wt2/lf_fermion_skeptic/`. lf_fermion's code and recorded lines were only read.

## Verdict in one paragraph

The method stands where it was regressed: an independent exact mode-function solver reproduces the homogeneous
control to ≤ 0.1% through the production burst and to ≤ 1.3% (χ) / ≤ 2.3% (N1) later. The exported history is
energy-consistent. Three things need correcting. (1) The reported line energy Q carries a grid-dependent bias, and the
quoted "fermion grid ×2: +2–5%" check measures that artefact rather than convergence. With a clean estimator the χ
species are unbiased after fragmentation (≤ 2%), but N1 is high by 14–23% and N2 by 9–19% at t = 0.4–0.5.
(2) λ_out mixes a common-mode closure deficit, which already exists without inhomogeneity, with the locality effect
the brief asks about. Split species by species, the locality part is λ_loc ≈ −1.5 to −1.7 (N2), −0.4 to −0.5 (N3),
≈ 0 (N4) and −0.7 to −1.4 (N1, bias-corrected). (3) The post-fragmentation enhancement is carried by hot spots one to
two lattice cells wide, so κ_frag depends on the lattice resolution until dx is varied. Nothing refutes Q1's headline:
the output-level closure lies beyond P0 for every species.

## S1. Independent homogeneous solver [C] (`sk_hom_exact.py`, `results/sk_hom_REC_P0_l*_f1_hermite.json`)

The independent solver works per k, in cosmic time, with a 4th-order Gauss–Magnus integrator. It takes the exact SU(2)
exponential of the Magnus generator, including the commutator term (Blanes, Casas, Oteo, Ros, arXiv:0810.5488 [R]). It
uses no FFT, no line grid and no Strang splitting, and it reads lines.bin with its own reader. It ran on the line-mean
histories of χ line 0 and σ line 3.

| t | 0.20 | 0.22 | 0.30 | 0.50 | 0.70 | 1.10 |
|---|---|---|---|---|---|---|
| N2 sk/lf | 1.0017 | 1.0007 | 1.0007 | 0.991 | 0.989 | 0.990 |
| N3 | 1.0010 | 1.0002 | 1.0001 | 0.997 | 0.997 | 0.997 |
| N4 | 1.0009 | 1.0001 | 1.0000 | 0.999 | 0.999 | 0.999 |
| N1 | 1.0045 (t = 0.2) | 1.0004 (0.25) | 0.9999 | 0.99 (0.4) | 0.98 (0.6) | 0.98 (1.0) |

The 1–2% differences after t = 0.34 may come from my k_max (3.06×10⁴, against their 6.1×10⁴ after the refinement); this
was not separated. **Q_hom and the ratios Q_hom/closure_hom = 1.54 / 2.11 / 2.48 / 1.67 (N2 / N3 / N4 / N1) stand.**

**Mechanism [C on the numbers, A on the reading].** I resolved Q_exact/closure in time.
- **N2:** the excess forms in the first swing, t = 0.19 → 0.21 (ratio 0.35 → 1.82). Q there exceeds the closure's own
  capacity F_sat(E)·E by ×1.2–1.5. This matches the constant-amplitude overshoot F_tot(2)/F_sat = 1.39 in R1.
- **N4:** the ratio grows steadily, 1.12 → 2.47, over t = 0.20–0.50. Over that interval the closure is already
  saturated (TE ≈ cap).

So NOTES §3, reading 1 ("the deficit is the changing amplitude, the band edge shrinking") is at most half the story.
The non-monotonic onset F_tot(n), an overshoot the exponential-approach closure cannot represent, carries N2 and
part of N3.

## S2. λ_PB: common mode against locality [C] (`sk_lambda_split.py`, `results/sk_lambda_split.json`)

I recomputed λ_out = (P0 − Q)/(P0 − GM) from `cmp_REC_P0.json`; it reproduces the NOTES values. Q/P0 is already
1.6–2.1 at t = 0.3, where P0 = GM. That common-mode part comes from the onset and is shared by both closures.
Dividing it out (r_hom = Q_hom/closure_hom applied to both P0 and GM) gives a locality-only parameter λ_loc:

| species | extra factor Q/(r_hom·P0), t = 0.7 / 1.1 | λ_out (lf) | λ_loc |
|---|---|---|---|
| N2 | 2.13 / 2.30 | −3.2 / −3.5 | −1.50 / −1.70 |
| N3 | 1.25 / 1.36 | −2.3 / −2.6 | −0.36 / −0.49 |
| N4 | 0.98 / 1.04 | −2.2 / −2.4 | +0.04 / −0.06 |
| N1 | 1.86 / 2.51 (bias-corrected ≈ 1.6 / 2.2) | −2.5 / −3.8 | −1.00 / −1.79 (corrected −0.72 / −1.40) |

**Correction.** "Fragmentation adds a further ×1.5–2.5 over P0" (NOTES §3) holds for N2 and N1 only. For N3 the
factor is ×1.25–1.36. For N4 it is ×1.0: once the common-mode onset deficit is fixed, **local blocking (P0) is right
for N4.** Its "beyond P0" result comes entirely from the homogeneous onset deficit.

Q1 should be reported per species, with the two parts kept separate. The r_hom used here comes from the line-mean
history, which dephases after fragmentation, so the split is approximate [A].

## S3. Energy-measure and grid-refinement bias [C] (`sk_line_drv.py`, `sk_measure_bias.py`, `results/sk_measure_bias.json`, `sk_n1_grid.json`, `sk_chi_grid.json`)

`closure_cmp.py` uses Q = e_proj_cut minus the mean over t < 0.15. I reran σ line 3 and χ line 1 with lf_fermion's own
`lfline.main`, read-only, with the output redirected. Three effects show up:

1. **Projection measure.** It projects at the line-mean phase. Once θ(x) is inhomogeneous, it carries an additive
   residual ∝ k_Nyq². The proj − chir excess on N1 line 3 at t = 0.3 is 6.5 / 28.5 / 117 (ρ/v² units) on grids
   4096 / 8192 / 16384, a factor ×4.4 and then ×4.1. The phase spread is |⟨e^{iΔθ}⟩| = 0.26 at t = 0.3.
2. **Discrete-map vacuum residual.** The chiral measure's residual scales as k_Nyq²/a²: baselines are 2.95 / 13.2 /
   56.1 on the three grids, and the early-time drift follows a⁻². A constant baseline taken at t < 0.15 therefore
   misses part of it.
3. **Refinement at t_ref = 0.34.** The new modes are filled with the homogeneous vacuum of the mean phase, while the
   line already has |⟨e^{iΔθ}⟩| = 0.07–0.33. This puts spurious particles at high k, and both measures count them.
   lf_fermion's `nf16k` test inherits the effect: on line 1 it gives +33% (chir) and +54% (proj) for N2 at t = 0.4.

**Clean estimator [model choice].** Run without refinement (Nc = Nf), use the chiral measure, and subtract a fitted
C/a² residual. It is grid- and seed-converged for N1: grid 8192 (seed 34), 8192 (seed 134) and 16384 give
183.7–185.7 at t = 0.5 and 198.1–199.2 at t = 0.7. Bias of the reported Q against this clean estimate:

| t | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 |
|---|---|---|---|---|---|---|
| N1 (σ line 3) | 1.04 | 1.23 | 1.18 | 1.16 | 1.14 | 1.15 |
| N2 (χ line 1) | 1.04 | 1.19 | 1.09 | 1.005 | 1.005 | 1.00 |
| N3 | 0.99 | 1.06 | 1.02 | 1.00 | 0.99 | 0.99 |
| N4 | 1.00 | 1.05 | 1.01 | 0.98 | 0.99 | 1.00 |

Seed-to-seed spread on χ line 1 is up to 12% (N2 at t = 0.5), consistent with lf_fermion's realisation scatter. Only
one line was tested per field. The size of the bias on the other lines depends on their phase spread at t_ref.

**Consequences.**
- χ results after fragmentation stand.
- N1: Q/P0 at t = 1.1 goes from 4.2 to ≈ 3.6, and κ_frag(N1) = 29 is biased high (its "lo" value of 7 is closer).
- Part of the pre-fragmentation self-consistency miss for N2 (Q/LF = 1.05–1.45 at t = 0.3–0.4) is this artefact.
- The convergence entry "fermion grid ×2: +2–5% (+12% N2 transient)" is not a UV-convergence test. With the clean
  estimator, N1 grid ×2 gives −1.0% to +0.1%. A clean χ check at 16384 was not run (budget).

## S4. Interpolation and hot spots [C] (`sk_hotspots.py`, `results/sk_interp_mu4.json`, `sk_kbin_fraction.json`)

- **Interpolation.** Band-limited interpolation of θ̇ from the 64 lattice sites to the fermion grid does not inflate
  the hot spots. The ratio of ⟨μ5⁴⟩ on the fine grid to that on the samples is 1.00 / 1.04 / 0.93 / 0.88 / 0.79 at
  t = 0.3 / 0.5 / 0.6 / 0.7 / 1.0, so the interpolation is not a source of excess Q. Linear interpolation gives
  0.54–0.74 instead, which means the hot spots are one to two lattice cells wide.
- **Energy above the homogeneous band edge.** At t = 0.7–1.0, the share of energy at k > 1.5×10⁴ (above the
  homogeneous a·μ5 ≈ 1.1×10⁴) is 81–83% (N2), 62–65% (N3), 44–46% (N4) and 69–74% (N1).
- **Concentration in the closure.** The top 5% of line sites carry 53–82% of the P0 tracer energy (5% at t = 0.3).
- **Streaming time.** Crossing one cell (dx = 0.0066 comoving) takes about one half-period, so blocking at these
  sites is relieved by streaming.

So the locality enhancement (S2) is a property of structure at the lattice cutoff. κ_frag and λ_loc are not
continuum-controlled until dx is varied with fermions (lf_fermion lists this as "not done"; it is upgraded here to the
main systematic, next to the ray reduction).

## S5. Exported history [C] (`sk_hist_check.py`, `results/sk_hist_check.json`)

For hist_HFa_LF.json, primary and four alternates:

| check | result |
|---|---|
| Σf − 1 | ≤ 8×10⁻⁶ |
| f_rad = Σf_N + f_SM | ≤ 5×10⁻⁸ |
| ln H + ∫ε_H H dt | ≤ 0.47% (0.40% for H84) |
| ln a − ∫H dt | ≤ 1.1×10⁻⁵ |
| window means of −Ḣ/H² against ε_H and 1.5(1+w) | ≤ 0.9% apart |
| lattice → fluid splice | continuous to 10⁻⁵ in f and 4×10⁻⁶ in H |

The largest jump of f_rad between samples (0.05, at t = 0.197) is the onset burst. The H84 alternate's
budget_absmax = 0.97%, above the ≤ 0.6% stated for the primary; this should be stated. **No continuity or energy
error.**

## S6. K estimate [E recomputed]

The χ-weighted rule K ∝ (comoving χ energy)^{2/3} gives:

| run | K (P0 or GM anchor; χ or χ+grad) | lf_fermion |
|---|---|---|
| KFmid | 5.06–5.37 | 5.4 |
| KFlo | 5.35–5.66 | 5.7 |
| KFhi | 4.54–4.83 | 4.9 |

The rule reproduces GM from P0 to 3.5–4.8%, not 2%. LF's χ fraction at t = 3 (0.164) lies outside the anchored range
(0.29–0.50), and its σ remnant (0.023, against 0.097 for P0) is ignored by the rule. So K ≈ 5.0–5.4 is an
extrapolation; only lf_chain's hand-over decides it.

## Corrections for the chain

1. Report Q1 per species and split into a common-mode part (r_hom = 1.5–2.5) and a locality part
   (λ_loc ≈ −1.6 / −0.45 / 0 / −1.0 for N2 / N3 / N4 / N1 at t ≈ 1). Do not quote λ_out ≈ −2 to −4 as a locality
   measure.
2. Lower N1's Q by ≈ 15% after t = 0.4, which lowers κ_frag(N1). Treat N2's Q at t = 0.4–0.5 as ≈ 10–19% high.
   Replace the "grid ×2: +2–5%" entry with "artefact; clean N1 grid ×2 ≤ 1%". Use the chiral measure without
   refinement, or fill new modes with the local vacuum, in any rerun.
3. Reword NOTES §3, reading 1: the onset overshoot (N2, N3) and band-edge shrinkage (N4) are both at work.
4. State that κ_frag is a lattice-cutoff (dx) quantity until a fermion run on a finer lattice exists.
5. Use K ≈ 5.0–5.7 [E], marked as an extrapolation; lf_chain must take K from the kinetic hand-over and the m_s²
   column of hist_HFa_LF.json.
6. Note the H84 budget of 0.97%.

## Not done

- An exact many-body small-box check of the stochastic estimator. It is unbiased by construction; seed checks gave
  0.2–1% (N1) and up to 12% (N2, one line).
- A clean 16384 run for χ.
- The ray reduction and 3D helicity mixing.
- The time step on inhomogeneous lines, beyond lf_fermion's own check.
- Variation of the lattice spacing.
- Q2, Q3 and H1 (lf_chain). Q4 was not attacked.

## Literature

- S. Blanes, F. Casas, J. A. Oteo, J. Ros, *The Magnus expansion and some of its applications*, Phys. Rep. 470 (2009)
  151, arXiv:0810.5488 (Gauss–Magnus integrator of S1).
- Sz. Borsanyi, M. Hindmarsh, *Low-cost fermions in classical field simulations*, PRD 79 (2009) 065010,
  arXiv:0809.4711 (the estimator under test).
