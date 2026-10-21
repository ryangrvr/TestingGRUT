# Adversarial verification of audit track `estimators`

21 October 2026. Verifies `audit/estimators/` (NOTES.md, results, code), which audits Codex `a69bd86`
(`research/reality_framework/grut_expanding_charge_20261021/`, DERIVATION.md §2, REVIEW.md, RESULTS.json,
raw_diagnostics.tar.gz) against Claude's packet `c7b8480` (`research/reality_framework/grut_issue13_wall_history_20261020/`).
Both packets were read from `git archive` copies in the scratch directory. Nothing outside this directory was edited.
None of the audited track's code was reused. Every check below uses a different method from the one the track used.

Labels: [D] derived, [C] computed here, [E] estimate, [P] preparation input.

## Files

| File | Method (differs from the audited track) |
|---|---|
| `code/verify_bounds_Q.py` → `results/bounds_Q.json` | Bounds by Monte-Carlo over 4×10⁵ random discrete energy-measure distributions (the track used LP). Capped sharp bound by brute-force two-point search (the track used the chord formula). Thin-wall stress and the Q identity derived symbolically with sympy by boosting the rest-frame surface stress. Kernel fractions by quadrature, including from the z-profile itself. |
| `code/verify_estimator.py` → `results/estimator.json` | (A) 3D planar kinks at arbitrary orientation, offset-averaged continuum sampling, isotropic orientation average (the track used 1D only). (B) Own leapfrog for the exact lattice reductions of 3D walls with normals (1,0,0), (1,1,0) and (1,1,1), i.e. 1D steps h = dx, dx/√2, dx/√3 at network_sim's dt = dx/3. (C) network_sim run series read directly. |
| `code/verify_tpar_proxies.py` → `results/tpar_proxies.json` | T_par rerun with Claude's `cosmo_observable.core.history`, imported read-only from the scratch copy. Codex proxies parsed directly from the raw CSVs with my own parser, and checked against RESULTS.json. |
| `code/realistic_R.py` → `results/realistic_R.json` | New: where smooth unimodal speed distributions sit inside the rigorous interval. |

Reproduce, with the scratch extraction per BRIEF.md (`<P>` = Claude packet, `<X>` = Codex packet, `<raw>` = extracted `raw/`):

```
python3 -I code/verify_bounds_Q.py results/bounds_Q.json
python3 -I code/verify_estimator.py <P>/modules/network_sim results/estimator.json      # ~10 min, 1 core
python3 -I -B code/verify_tpar_proxies.py <P> <X> <raw> results/tpar_proxies.json
python3 -I code/realistic_R.py results/realistic_R.json
```

## 1. Lorentz-weighting bounds

- Monte-Carlo, 4×10⁵ distributions with up to 6 speeds [C]: no violation of (1−s_E)^(−1/2) ≤ R ≤ (1−s_E)^(−1).
  min R/lower = 1 + 1.6×10⁻¹², approached by near one-speed draws. max R/upper = 0.9999995, approached by a rest + near-luminal draw. **CONFIRMED.**
- Brute-force two-point maximum with cap v_max: 1.14937, 1.16934, 1.18829, 1.22284, 1.24850 for v_max = 0.8, 0.9, 0.95,
  0.99, 0.999. The maximum always has support {0, v_max}. This equals the track's chord formula to machine precision [C]. **CONFIRMED.**
- Area-weighted mixture: R = 2.291790 (v_f = 0.99), 5.4517, 15.499, 147.83. s_E of the 0.99 mix = 0.64317 [C]. **CONFIRMED.**
- Interval at s = 0.207936: 1.123621–1.262524 [C]. **CONFIRMED.**
- Addition [D]: (1−s)^(−1/2) is a lower bound under **either** weighting (Jensen for concave √(1−x) in the energy measure,
  convex γ in the area measure). So Claude's "up to ≈ 12% at the one-speed γ" points the wrong way whichever weighting the
  estimator has. This strengthens the track's REFUTED verdict on REPORT.md:345–346.
- Addition, where realistic distributions sit [C, E]: at s_E = 0.207936, a half-normal in γv gives R = 1.1335.
  A Maxwell-like distribution in γv gives 1.1273, uniform in v gives 1.1311, β(2,2) gives 1.1294 and β(½,½) gives 1.1321. All lie within 3–7% of the interval width above
  the lower edge. The upper half of the interval needs about 21% of the energy in near-luminal walls. The rigorous statement "12–26%"
  stands, but a smooth population gives ≈ 13% (T_par factor ≈ 1.063). Claude's number happened to lie near the probable value.
  Its logic (one speed = maximum) is wrong.

## 2. What network_sim's 0.456 is

Recomputed from `results/runs/*.json` [C]:
- v_rms_win = √(window mean of v2_wall) over h ∈ [1/8, 1/4]. Seeds 0.455920, 0.462422, 0.449998. The mean is 0.456113 and the ddof=1 std is 0.006214.
  The mean of v_rms² is 0.208065. Window v2_wall values: 0.1941–0.2178. The 384³ seed gives 0.4680. **CONFIRMED.**
- 3D PRS max v2_wall 0.689–0.727 at τ ≈ 7 (formation), before the window τ ∈ [36, 62]. **CONFIRMED.**
- Physical-equation runs: max v2_wall 2.09–2.78 (2D, τ ≈ 3.6–3.75) and 2.11–2.22 (3D, τ = 3.875). 2D PRS α = 2: 1.003–1.495.
  The track's "2.1–2.8" and "1.0–1.5" are **CONFIRMED.** network_sim NOTES:215 already uses v2_wall > 1 as a flag.

## 3. Estimator weighting and its calibration at network_sim resolution (dx = 1 = δ/2)

**Weighting [D].** For φ = f(γ(n·x − vt)), φ̇ = −vγf′ and ∇φ = γf′n pointwise. The core-integrated |∇φ|² is ∝ γ dA ∝ energy.
Fine-grid results [C]: sampled 3D isotropic at dx = 0.1 gives est/s_E − 1 = +0.1% to +0.3% and est/s_A − 1 = +11% to +40% for mixtures.
Evolved (1,1,1) at dx = 0.5 (h = 0.289 ≈ δ/7) gives +1.3% to +6.0% vs s_E and +20% to +46% vs s_A. **CONFIRMED** (energy weighting for resolved
isolated walls).

**Discretisation at dx = δ/2. CORRECTED (refined).** Only the direction is one-sided: the estimator reads HIGH.
- Sampled continuum (pure estimator error), one speed 0.456 [C]: axis (1,0,0) gives +16.7%, (1,1,0) gives +8.2% and (1,1,1) gives +5.5%.
  The **isotropic average is +9.7%** (dx = 1) and +2.4% (dx = 0.5). It grows with speed: +8.0%, +9.7%, +15.1% and +41% at v = 0.2, 0.456, 0.7 and 0.9.
  So at dx = δ/2 the estimator is not a clean energy mean: fast walls are under-weighted in the denominator. A rest + 0.9
  mixture reads +30% (isotropic, dx = 1) above s_E.
- Lattice-evolved, exact 3D reductions [C]. Each entry is the estimator relative to the measured lattice speed / relative to the continuum speed that was set:
  (1,0,0) +20.6% / +10.5% (kinks slowed 0.456 → 0.436), (1,1,0) +8.7% / +5.3%, (1,1,1) +5.4% / +3.3%. At dx = 0.5:
  (1,0,0) +4.2% / +2.5%, (1,1,1) +1.3% / +0.8%. This confirms the track's claim that the axis-aligned 1D case is the worst case. The
  effective step is exactly dx/√2 and dx/√3 for those normals.
- 3D PRS same-field pair (dx = 1 vs 0.5, L = 128) [C]: v2_wall ratio −3.57% to +1.46% over τ = 22.3–64 (6 outputs).
  The mean is −0.8% ± 0.8%. **CONFIRMED** as stated. But the isolated-wall models predict that dx = 1 reads ≈ +4–5% (evolved, isotropic interpolation by Σn_i⁴) to +7%
  (sampled) above dx = 0.5, and the network shows no such trend. So the network-level discretisation effect is at most a few
  percent, or it is compensated by slower lattice walls. The track's "a few percent to ≈ 20%" upper end is the axis-aligned
  worst case relative to slowed lattice kinks, not an isotropic network estimate.
- **Missing case.** The track propagated only upward variations of s (seeds, 384³: rms 0.450–0.476). The calibration bias
  pushes the continuum ⟨v²⟩_E **down**, to ≈ 0.19–0.21 [E]. For s_E/1.10 = 0.189: R ∈ [1.1104, 1.2331], Q t/ρ_w = 0.2165, and
  T_par (A = 0.8, t) is 6.135–6.465 keV without backreaction and 6.709–7.070 keV with it. For s_E/1.20 = 0.173: Q = 0.240.
  The Q range quoted in the track's replacement text, "0.14–0.19", should read ≈ 0.14–0.22 [E].
  (Indicative only: Leite–Martins' 3D PRS γv = 0.46–0.48 would mean ⟨v²⟩ ≈ 0.175–0.187 for one speed. Their estimator
  and resolution differ, so this is consistent with a read-high bias but is not evidence for one.)

## 4. Thin-wall stress and Q identity

sympy, boosting the rest-frame surface stress diag(σ, −σ, −σ, 0) along the normal and dividing by γ for the lab-frame
delta function [D]: T⁰⁰ = γσ, T^zz = γσv² and T^xx = T^yy = −σ/γ. Hence w = v² − 2/3, and for a population w = ⟨v²⟩_E − 2/3.
With ρ̇ = −3H(1+w)ρ − Q and ρ = Aσ/t: **Q t/ρ_w = 1 − Ht(1 + 3s)**. In radiation this is (1 − 3s)/2; off-scaling it is (−μ − 3s)/2.
Numbers [C]: 0.188096 (s = 0.456²); −0.184 (s = 0.456 literal); sign change at s = 1/3; 0.1601–0.1963 (rms 0.450–0.476);
0.1361 and 0.1756 (μ = −0.896, −0.975). **CONFIRMED.** Claude's REPORT:338–339 has the right number (0.19) with the wrong label (⟨v²⟩ ≈
0.456). The track's CORRECTED verdict is **CONFIRMED**.

## 5. T_par (Claude's integrator, rerun)

[C] A = 0.8, convention t: area 5.82206 → 6.17139 (R_lo) – 6.54167 (R_hi) keV without backreaction; 6.36686 → 6.74885–7.15374
keV with backreaction. Full bracket 5.44609–7.34400 → 5.77286–8.25186 keV. T ratio 1.05999–1.12362 = R^(1/2) in every row,
with backreaction and both conventions included. **CONFIRMED.**

**HKS A is an area count. NOT CHECKED in the track; CHECKED here.** I read the arXiv:1309.5001 text (fetched PDF, data only).
Their A (Eq. 13) is A·t/(aV) with the comoving area from PRS link counting weighted by |∇φ|/(|φ,x|+|φ,y|+|φ,z|) (Eq. 14),
times the rest tension σ_wall (Eq. 12). They also compute a field-energy ρ_wall in wall cores (Eqs. 15–16), but only to check
the t⁻¹ scaling. A ≃ 0.8 ± 0.1 (Eq. 29) comes from the area count. So applying R is appropriate.
**Unstated assumption (missing case).** HKS evolve the physical equation. The R applied to their A uses s_E from 3D PRS α = 3. The wall speed
of a physical-equation network in the scaling regime was not measured (network_sim phys-3D v_rms ≈ 0.99 is a rolling-field
value). The transfer relies on PRS α = 3 reproducing wall momentum dynamics [E].

## 6. Codex charged-state proxies (raw CSVs)

My parser reproduces RESULTS.json final v2_core for all 12 runs exactly [C]. All 11 non-cold runs exceed 1 before formation.
The first exceedance is at a = 1.03–1.06, with |S|max = 0.019–0.032 at that sample. The pre-formation maximum is 3.24–10.31 at a = 1.51–1.90. The first post-formation
sample of the Gaussian runs (a = 1.57) is 2.149–2.269. Post-formation values swing between 0.040 and 0.880. The final core charge fraction is 0.300–0.467
in forming runs and 0.187–0.236 in the s101 family. main_s101 has 126 interface links at a = 2.2. The cold run's max v2_smooth is 52.33.
The kink controls give +1.38% (dx = 0.2) and +0.32% (dx = 0.1). All **CONFIRMED.**
**One detail CORRECTED:** the track's "box K_φ/G_φ oscillates between 0.5 and 3.3" holds only for main_s202. Across forming runs,
post-formation K_φ/G_φ spans 0.27–7.83: 0.27–3.27 in the sparse runs and 0.69–7.83 in the Gaussian runs. This strengthens the conclusion that the core mask holds
non-wall kinetic content.
The linear-mode argument [D] holds. A travelling massive mode gives exactly ω²/k² > 1. A standing mode gives ω²/k² as a ratio of time averages.
Tachyonic growth gives Γ²/k², and a homogeneous mode gives an unbounded ratio. **CONFIRMED.**

## 7. Kernel calibrations

Quadrature [C]: ∫¾(1−S²)dS = 1. The sharp |S| < ½ fraction is 0.6875, also obtained from the z-profile of tanh(z/√2). The (1−S²)-weighted fraction is 0.8.
**CONFIRMED.**

## 8. Replacement text in the track's §6

The text is accurate apart from these points. (i) README "T ≈ 6 keV" is at README.md:14, not :13. (ii) Item 1's "(0.14–0.19 …)" should
add "≈ 0.22 if the estimator reads 10% high at dx = δ/2". (iii) Item 2 could add "smooth unimodal distributions give ≈ 1.13; the upper edge needs
≈ 21% of the wall energy in near-luminal walls" and "R is taken from 3D PRS speeds" [E]. All other line citations checked (REPORT
338/345–346/61, network_sim NOTES 33/187/206/215/303, cosmo_observable core.py:189, charged_exit_sim analyze.py:160 and NOTES:200).

## Limits

- No 3D network was rerun. The calibration uses isolated planar walls (sampled and lattice-evolved). The network-level
  evidence is one same-field resolution pair with 6 outputs.
- Realistic-R families are illustrations. The 3D PRS speed distribution was not measured.
- HKS was read from the arXiv text only. No figure data were digitised.
- Nothing here is an exclusion, validation or identification of Candidate 1.
