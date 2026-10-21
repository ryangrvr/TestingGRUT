# Audit track: estimators (Lorentz weighting, wall speed, energy-loss identity, velocity proxies)

21 October 2026. Audits Codex `a69bd86` (`research/reality_framework/grut_expanding_charge_20261021/`, DERIVATION.md
§2, REVIEW.md, RESULTS.json, raw_diagnostics.tar.gz) against Claude's packet `c7b8480`
(`research/reality_framework/grut_issue13_wall_history_20261020/`; working tree identical). Codex files were read
from a `git archive` copy in the scratch directory. Neither packet was edited.

Labels: [D] derived here, [C] computed here, [E] estimate, [P] preparation input.

## Files

| File | Content |
|---|---|
| `code/lorentz_bounds.py` | Bounds, linear-programming extremal distributions, Codex counterexample, Q identity, sensitivity, T_par with Claude's `cosmo_observable.core.history`, sampled-profile test of network_sim's estimator (also through `wallsim.measure` itself) |
| `code/lattice_kinks.py` | Lattice-evolved 1D kinks with chosen speeds, network_sim discretisation; estimator vs measured area- and energy-weighted ⟨v²⟩ |
| `code/tpar_scan.py` | T_par for area A and energy-weighted A over A = 0.7/0.8/0.9, both conventions, with/without backreaction |
| `code/proxy_series.py` | Codex's v2_core_proxy / v2_smooth_proxy time series from their raw CSVs |
| `results/lorentz_bounds.json`, `results/lattice_kinks.json`, `results/tpar_scan.json`, `results/codex_proxy_series.json`, `results/verdicts.json` | Outputs |

Reproduce (scratch extraction of Codex raw data first, see BRIEF.md):

```
python3 -I -B code/lorentz_bounds.py <claude_packet> results/lorentz_bounds.json
python3 -I code/lattice_kinks.py results/lattice_kinks.json
python3 -I -B code/tpar_scan.py <claude_packet> results/tpar_scan.json
python3 -I code/proxy_series.py <scratch>/raw/raw results/codex_proxy_series.json
```

## 1. Lorentz-weighting bounds (DERIVATION.md §2)

Thin, unoccupied wall elements with common rest tension σ and normal speed v. Energy per area γσ. R = E/(σ A_area).
Write x = v² ∈ [0, 1).

**Energy-weighted s_E = ⟨x⟩_E** (weights γ σ dA). Then A_area = ∫dE/(γσ), so R = 1/⟨√(1−x)⟩_E [D].
- Lower: √ is concave, so ⟨√(1−x)⟩_E ≤ √(1−s_E), giving R ≥ (1−s_E)^(−1/2). Equality iff x is constant (one speed).
- Upper: √(1−x) ≥ 1−x on [0, 1], equality only at x ∈ {0, 1}. So ⟨√(1−x)⟩_E ≥ 1−s_E and R ≤ (1−s_E)^(−1).
  The supremum is approached by energy fraction 1−s_E at rest plus energy fraction s_E at v → 1 (whose area → 0).
  It is not attained for v < 1. With a speed cap v_max the sharp bound is
  R ≤ 1/[1 − (s_E/v_max²)(1 − √(1−v_max²))]: 1.149, 1.169, 1.188, 1.223, 1.248 for v_max = 0.8, 0.9, 0.95, 0.99, 0.999 [C].
- An LP over 20001 speeds with ⟨x⟩_E fixed returns R_min = 1.123621 (support at v = 0.456) and R_max = 1.262193
  (support {0, 0.9999995}, energy weights 0.79206/0.20794) [C].

**Area-weighted s_A.** R = ⟨(1−x)^(−1/2)⟩_A. The function is convex, so R ≥ (1−s_A)^(−1/2) [D]. It diverges as
x → 1, so no upper bound exists. Codex's mixture (rest + v = 0.99, area fraction s_A/0.99² = 0.212158) gives
R = 2.291790 [C]. With v_fast = 0.999, 0.9999, 0.999999: R = 5.45, 15.5, 147.8 [C]. That mixture has s_E = 0.643,
so it is not compatible with an energy-weighted measurement of 0.208.

**Numbers at s = 0.456² = 0.207936:** (1−s)^(−1/2) = 1.123621, (1−s)^(−1) = 1.262524 [C]. A_E/A_area is 12.4–26.3%
above one. Seed scatter 0.456 ± 0.006 gives R_lo 1.120–1.128, R_hi 1.254–1.271. The 384³ seed (0.468) gives
1.132–1.280. The window time variation of s (pooled variance 3.5×10⁻⁵) changes the time-averaged bounds by ≤ 10⁻⁴.

Codex's kernel calibrations were also checked analytically: ∫₋₁¹ ¾(1−S²) dS = 1; sharp |S| < ½ energy fraction
∫(1−f²)df over |f|<½ / ∫ over |f|<1 = (11/12)/(4/3) = 11/16; (1−S²)-weighted fraction (16/15)/(4/3) = 4/5 [D].

**Verdicts.** Energy-weighted bounds CONFIRMED. No upper bound from s_A CONFIRMED. 1.12362–1.26252 CONFIRMED.
"One-speed value is the LOWER edge, not an upper limit of 12 percent" CONFIRMED. Claude's "up to ≈ 12% at the
one-speed γ = 1.124" and "parity earlier by ≲ 6%" REFUTED.

## 2. What network_sim measured as v

`wallsim.measure` (lines 259–277): v2_wall = Σ_{|φ|<½} π² / Σ_{|φ|<½} |∇_c φ|², where π = dφ/dτ (conformal, the
leapfrog half-step momentum) and ∇_c is the central difference. It is one ratio of sums per output. `analyze.py:82–83`
averages v2_wall over the window h = τ/L ∈ [1/8, 1/4] and stores v_rms_win = √(mean). The quoted 0.456 ± 0.006 is the
mean ± seed std of v_rms_win over three 256³ PRS seeds (0.4559, 0.4624, 0.4500; window values of v2_wall 0.194–0.218).
So **0.456 is an rms speed; the measured mean square is 0.208.**

Weighting [D]: for a boosted kink φ = f(γ(x−vt)), φ̇ = −vγf′ and ∇φ = γf′ pointwise, so each wall element contributes
v² × ∫_{core}|∇φ|² = v² × γ × const × dA (contracted width 1/γ, gradient² γ²). The weight is ∝ γ dA ∝ energy. The
ratio is the **energy-weighted** ⟨v_n²⟩_E for thin, unoccupied, resolved walls.

Lattice tests [C]:
- Lattice-evolved 1D kinks (`lattice_kinks.py`, network_sim discretisation, measured kink speeds). At dx = δ/8 the
  estimator matches the measured s_E to +0.4% to +4.8% and, for the mixtures, differs from s_A by +11% to +44% (rest + 0.9 mixture:
  est 0.573, s_E 0.547, s_A 0.397; uniform 0–0.8: 0.269, 0.265, 0.227). Energy weighting CONFIRMED.
- At network_sim's resolution dx = δ/2 (axis-aligned), the estimator reads +20% to +27% above the lattice kinks'
  measured s_E, and lattice kinks are slowed (set 0.456 → measured 0.437; 0.9 → 0.63).
- The 3D PRS resolution pair with the same continuum field (dx = 1 vs 0.5, L = 128) agrees to −3.6% to +1.5% in v2_wall
  over τ = 22–64. Oblique walls see an effective spacing between dx/√3 and dx, so the 1D axis-aligned number is a
  worst case. **The discretisation uncertainty in ⟨v²⟩_E = 0.208 is not established; a few percent to ≈ 20% [E].**
- Sampled continuum profiles fed to `wallsim.measure` give the same weighting conclusion and the same O(dx²) bias.

## 3. Energy-loss identity

Stress of a thin wall element [D]: energy γσ, normal stress γσv², two tangential −σ/γ. Hence p = ρ(⟨v²⟩_E − 2/3) for
a population. With ρ̇_w = −3H(1+w)ρ_w − Q and scaling ρ_w = A_E σ/t:

    Q = [1 − Ht(1 + 3⟨v²⟩_E)] ρ_w/t.

Radiation (Ht = ½): Q t/ρ_w = (1 − 3⟨v²⟩_E)/2.
- ⟨v²⟩_E = 0.207936 → **0.188096** [C].
- 0.456 used literally as ⟨v²⟩ → **−0.184** (walls would gain energy) [C]. The sign changes at ⟨v²⟩_E = 1/3.
- Seed/resolution spread (rms 0.450–0.476): 0.160–0.196 [E].
- Off exact scaling: if ρ_w ∝ t^((μ−1)/2) with area-density exponent μ, Q t/ρ_w = (−μ − 3⟨v²⟩_E)/2. With the measured
  3D PRS μ = −0.896 to −0.975 this is 0.136–0.176 [E].

**Verdict.** Codex's statement CONFIRMED. Claude's REPORT computed Q ≈ 0.19 correctly but labelled the input
⟨v²⟩ ≈ 0.456 (CORRECTED: rms 0.456, ⟨v²⟩_E ≈ 0.208).

## 4. Effect on T_par (Claude's cosmo_observable integrator, area-count A)

T_par ∝ (Aσ)^(1/2) without backreaction. With R = 1.1236–1.2625 [C, tpar_scan.json]:
- A = 0.8, convention t: 5.822 → 6.171–6.542 keV (no backreaction); 6.367 → 6.749–7.154 keV (backreaction).
- Full bracket (A 0.7–0.9, both conventions, ±backreaction): 5.446–7.344 keV (reproduces REPORT) → 5.773–8.252 keV.
- T_par factor 1.060–1.124; t_par factor 0.890–0.792.

This assumes the A = 0.8 ± 0.1 of HKS is an area count ("same definition" per network_sim NOTES:205; NOT CHECKED here
against the paper). ΔN_eff ≥ 7.4 is set at parity by ρ_hidden ≥ ρ_rad and does not depend on this factor.

## 5. Codex's charged-state velocity proxies

Codex's v2_core_proxy (expanding.cpp:166, 178) is Σ_{|S/a|<½}(χ_S′ − (a′/a)χ_S)² / Σ|∇_c χ_S|² = S′²/|∇S|², the same
estimator as network_sim. v2_smooth_proxy uses weight 1 − (S/a)² on |S/a| < 1. Their boosted planar kink control
reproduces 0.456² to +1.4% (dx = 0.2) and +0.32% (dx = 0.1) [C, from raw CSVs].

Final values (RESULTS.json, confirmed from raw CSVs) [C]:

| Run | v2_core | v2_smooth | outcome |
|---|---|---|---|
| main_s101, dt_half, grid256, grid384 | 1.263, 1.270, 1.249, 1.263 | 1.73–1.79 | one vacuum (genuine− ≈ 10⁻⁶) |
| large_s101 | 3.188 | 2.842 | one vacuum |
| main_s202, main_s303 | 0.175, 0.160 | 0.263, 0.250 | two vacua |
| kappa4_s101 | 0.404 | 0.455 | two vacua |
| gaussian_s101, s202, large_s101 | 0.301, 0.280, 0.336 | 0.433, 0.394, 0.476 | two vacua |
| cold_s101 | 0 (empty core mask) | 8.52 (max 52.3) | wall-free |

Time series [C]: every charged run exceeds one in the restored phase before formation (first at a = 1.03–1.06;
maxima 3.2–10.3 at a = 1.51–1.90). There the mask is the whole box. S is a small, nearly homogeneous growing or
oscillating fluctuation (|S| ≲ 0.05), and for a linear mode S′²/|∇S|² = ω²/k² > 1 (massive) or Γ²/k² (tachyonic).
In the Gaussian runs the first post-formation sample (a = 1.57) is 2.15–2.27. Later samples in the forming runs swing between 0.04 and 0.88
from one output to the next, while the box K_φ/G_φ oscillates between 0.5 and 3.3. At the end the mask holds 30–47% of the
charge in the forming runs (X-rich restored cores).
In the one-vacuum s101 family the mask contains almost no interface (main_s101: 126 interface links) but 19–24% of
the charge. Its oscillating cores keep the ratio above one to a = 2.2, independent of dt and of refinement to 384³.

**Verdict.** CONFIRMED: charged-state core proxies exceed one. This diagnoses that the mask holds non-wall content
(oscillating restored cores, bulk S oscillation, linear modes), not superluminal motion; a homogeneous oscillation has
∇S = 0 and an unbounded ratio. Values below one in charged runs are not calibrated speeds either. Claude's own
network_sim shows the same failure: v2_wall reaches 2.1–2.8 at formation in physical-equation runs and 1.0–1.5 in 2D
PRS with alpha = 2. network_sim NOTES §5.3 already uses v2_wall > 1 to flag non-wall sign changes, and calls phys-3D v_rms ≈ 0.99
"rolling-field values". The 3D PRS series behind 0.456 never exceeds 0.73 and its window is post-formation.

Claude's packet does not apply the π²/|∇φ|² proxy to the charged network. charged_exit_sim reports a coarsening rate
dξ_c/dτ (0.2–0.8). analyze.py:160 calls it a "flat-space wall velocity proxy". The neutral Lorentz factor of §1 is not
calibrated for the charged network.

## 6. Places in Claude's packet and replacement text

Required:

1. REPORT.md:338–339. Replace "A scaling network with ⟨v²⟩ ≈ 0.456 loses energy at\n  Q ≈ 0.19 ρ_w/t, and
   cosmo_observable evolves no receiver for it." with:
   "For a scaling network ρ_w = A_E σ/t with wall stress w = ⟨v²⟩_E − 2/3, energy conservation gives
   Q = [1 − Ht(1 + 3⟨v²⟩_E)] ρ_w/t. network_sim's 3D PRS estimator gives an rms speed of 0.456 ± 0.006, i.e. an
   energy-weighted mean squared speed ⟨v²⟩_E ≈ 0.208 (calibrated for thin unoccupied walls only). In radiation
   (Ht = 1/2) the loss is Q ≈ 0.19 ρ_w/t [E] (0.14–0.19 over the measured area-density exponents −0.90 to −0.975);
   reading 0.456 itself as ⟨v²⟩ would reverse the sign. cosmo_observable evolves no receiver for it."
2. REPORT.md:344–346. Replace the bullet text after "*Area A versus energy A.*" with:
   "Agreed that they differ. network_sim measures wall area (as does, per network_sim NOTES, the HKS A = 0.8 used
   for T_par); a wall element moving at normal speed v carries γσ per area. For thin unoccupied walls network_sim's
   estimator Σπ²/Σ|∇φ|² over |φ| < 1/2 is the energy-weighted ⟨v²⟩_E ≈ 0.208 (rms 0.456). Then
   E/(σA) = 1/⟨(1 − v²)^(1/2)⟩_E lies between (1 − ⟨v²⟩_E)^(−1/2) = 1.124 (all walls at one speed) and
   (1 − ⟨v²⟩_E)^(−1) = 1.263 (walls at rest plus walls with v → 1). The one-speed γ = 1.124 is the lower edge, not an
   upper limit: the energy-weighted A is 12–26% larger than the area A [E]. Since T_par ∝ A^(1/2), T_par rises by a
   factor 1.06–1.12 (5.82 → 6.17–6.54 keV without backreaction, 6.37 → 6.75–7.15 keV with it; 5.45–7.34 →
   5.77–8.25 keV over the full bracket), so parity comes earlier. The direction strengthens the rejection; ΔN_eff ≥ 7.4
   is unchanged. The exact factor was not measured, and the estimator is not calibrated for the charged network."
3. network_sim/NOTES.md:33, row label and value. Replace "| Wall velocity v_rms (estimator sum pi^2 / sum grad^2 over
   abs(phi) < 1/2) | 3D PRS 0.456 +- 0.006;" with "| Wall rms speed v_rms = sqrt(<v^2>_E): square root of the window
   mean of sum pi^2 / sum (central grad phi)^2 over abs(phi) < 1/2, which for resolved thin unoccupied walls is the
   energy-weighted mean squared normal speed (a proxy elsewhere; values > 1 occur at formation) | 3D PRS 0.456 +- 0.006
   (<v^2>_E = 0.208); dx = 1 vs 0.5 same field within -3.6 % to +1.5 %; a 1D axis-aligned lattice kink at dx = delta/2
   reads up to 20 % high in v^2 [E];"

Recommended clarifications (formulas or labels are not wrong):

4. REPORT.md:61 "w ≈ −2/3 + ⟨v²⟩" → "w = ⟨v²⟩_E − 2/3 (⟨v²⟩_E energy-weighted)".
5. REPORT.md:57 append after "(5.45–7.34 keV across A, backreaction and convention)": " for an area-count A; energy
   weighting raises T_par by a factor 1.06–1.12 [E] (§6)". README.md:13 "T ≈ 6 keV" → "T ≈ 6–7 keV". README.md has
   no ⟨v²⟩, 12% or 6% wording.
6. network_sim/NOTES.md:187 "(A = 1.04, v_rms = 0.99 there are rolling-field values, not wall properties)" →
   "(A = 1.04 and the velocity proxy v_rms = 0.99, with single outputs up to v2_wall = 2.2, are rolling-field values
   outside the estimator's thin-wall calibration, not wall properties)". Same wording at :303.
7. network_sim/NOTES.md:206 "v_rms = 0.456 (gamma v ~ 0.51; estimators differ)" → "v_rms = 0.456, <v^2>_E = 0.208
   (gamma v ~ 0.51 only if every wall had the rms speed; estimators and weightings differ)".
8. network_sim/code/wallsim.py:259 comment → "# energies ... and wall-speed estimator: v2_wall is the energy-weighted
   <v_n^2> for resolved thin unoccupied walls; a proxy otherwise (> 1 flags non-wall content)". analyze.py:83 add
   comment "# rms speed; v_rms_win**2 (not v_rms_win) is <v^2>_E".
9. cosmo_observable/NOTES.md:176–180, after "n = 3(1 + w) = 1 + 3 v^2": "For a population v^2 is the energy-weighted
   <v^2>_E; network_sim's radiation-era value 0.208 (rms 0.456) would give n = 1.62 without energy loss."
   cosmo_observable/code/core.py:189 "n_post = 2 v_rms^2 = 1/3 network" → "n_post = 2 (<v^2>_E = v_rms^2 = 1/3)
   network" (current wording is consistent; clarification only).
10. Outside the four named locations: charged_exit_sim/code/analyze.py:160 "flat-space wall velocity proxy" →
    "coarsening rate d xi_comoving / d tau (domain-scale growth; not a wall speed)"; charged_exit_sim/NOTES.md:200
    "peculiar coarsening speed in units of c" → "dimensionless coarsening rate (not a wall speed)".

No occurrence in code/history.py, code/verify.py, MANIFEST.json or the result JSONs (number matches there are
coincidental). verification/v_network/VERIFY.md a4/b3 label v_rms consistently and need no change.

## Limits

- The lattice test is 1D, flat-space and axis-aligned. It uses no PRS friction and is not a 3D network rerun.
- The bounds assume thin, unoccupied walls of common rest tension. They do not apply to occupied or charged walls,
  restored cores or radiation, for which no calibrated E/(σA) exists.
- Whether HKS A = 0.8 is an area count was taken from network_sim NOTES, NOT CHECKED against the paper.
- The Q identity assumes exact scaling and the thin-wall stress. The scaling window is box-limited (μ ≠ −1).
- Codex's proxies were read from their raw CSVs; no Codex run was repeated here.
- Nothing here is an exclusion, validation or identification of Candidate 1.
