# Verification of track `bound` (independent, adversarial)

21 October 2026. Verifier track `v_bound` of the stage-2 construction pass (`../../BRIEF.md`, `../../BRIEF_STAGE2.md`).
Object of verification: `../bound/NOTES.md`, `../bound/code/`, `../bound/results/` as they stood on this date.

**Labels.** [D] derived here, [C] computed by `code/` (outputs in `results/`), [Q] stated by a source at the read
level given, [E] estimate, [A] argued, [R] recalled. Verdicts: CONFIRMED / CORRECTED (with the corrected value) /
REFUTED / NOT CHECKED.

**Independence.** No code from `../bound/code/` is imported or copied. The checks use different formulations:
- symbolic checks start from minisuperspace actions, not from bound's intermediate formulas;
- the FRW integrations use N = ln a and solve H² from the constraint, where bound evolves H in t;
- the Hu-Sawicki background is solved in scalar-tensor form (φ = f_R), where bound uses HS's reduced (y_H, y_R) system;
- the Eöt-Wash curve is re-digitized by a separate SVG pipeline from my own download of the arXiv source.

Downloads are in the scratchpad (`.../qtheory/v_bound_dl/`, `.../qtheory/v_bound_work/`), none in the repository.

---------------------------------------------------------------------------------------------------------------------

## 0. Summary

Every load-bearing conclusion of `bound` survives:
- the class-wide adiabatic relation and its coefficient;
- the numerical validation;
- the Eöt-Wash bound;
- the 3.5e-62 maximum;
- the regime (ii) exclusion;
- the f(R) map and the meaning of P_q;
- the vacuum structure of the three published models;
- the threshold criterion;
- the Hu-Sawicki solution.

**Corrections (6):**
1. The relic growth factor (A.7) is low by a factor 3.50. The overproduction becomes 7.8e5–2.3e9 (bound: 2e5–7e8).
2. The Appleby-Battye bound G_Mink/G ≳ 1e59–1e122 holds only under AB's own sufficient conditions. The exposed f_R0
   data require only ≳ 5e3 to ≳ 2e6.
3. The affine-gauge sentence in B.5 ("holds automatically … ~1e113") is wrong. The gauge needs the same tuning.
4. DESI diagonal distances: 4.1/4.1/5.2/5.7, not 4.1/4.2/5.3/5.9.
5. An un-noted stage-1 update: ledger row 8 gives 9e-39 ρ_Λ. That figure is the gravitating ρ_V. The effective dark
   energy is 8.1e-38 ρ_Λ.
6. A label: the verdict-table entry A5 needs [E] on "for any m_q ≲ H(z_eq)".

**Missing credit (1).** Wang, Hui, Khoury, PRL 109, 241301 (2012), arXiv:1208.4612 (abstract read). They prove in
general that chameleon-like scalars cannot self-accelerate, which is B.6's "ΛCDM-mimic" conclusion.

**Qualifications, no change of verdict:**
- the CPL small-field slope depends on the fit weighting;
- the error coefficient c can take either sign (a fourth family gives c = −42);
- the 14.5σ G_cos/G_lab comparison is not like-for-like.

No REFUTED item.

## 1. Verdict table

Full list with sources: `results/verdicts.json` (V01–V37).

| id | claim (bound section) | bound | v_bound | verdict |
|---|---|---|---|---|
| V01 | ρ_dyn = [(9/4)ρ_d² + 2ρ_dρ_r − 4Λρ_r]/(3Km²) (A.2) | [D] | identical; the 00 equation was re-derived from S = ∫N a³[f(R)/2 − ρ] and time derivatives were taken through the background equations (§2.1) | CONFIRMED |
| V02 | Ω_dyn(t0) = (9/4)Ω_m²(H0/m)², w = +1, no Λ dependence (A.2) | [D] | identical | CONFIRMED |
| V03 | Λ_cal = −q0δμ + [−1/(2eps'') + (4/3)q0²/(Km²)]δμ² (A.2) | [D] | identical with Taylor data up to eps'''' and K'''. It equals stage-1 equations §5.2 through c*(δμ) = q0δμ + δμ²/(2eps'') | CONFIRMED |
| V04 | effective DE = 9 × stage-1 ρ_V; formal m_q = 0.57 H0 (A.2) | 9; 0.57 | 9; 0.5709 (stage 1: 0.1903) | CONFIRMED |
| V05 | full-FRW ratio and c (A.3) | F1 0.9832/0.99848/0.99983, c = 15.2; F2 c = 53; F3 c = 27 | F1 0.98307/0.998481/0.999831, c = 15.05–15.18; F2 0.94344/0.994710/0.999409, c = 50.1–53.2; F3 0.97008/0.997267/0.999696, c = 26.5–27.4; new F4: c = −41.7 | CONFIRMED (\|c\| = 15–53, sign family dependent) |
| V06 | fails below m ≈ 10 H0 | F1 0.832, F2 0.595 | F1 0.832, F2 0.620 (depends on initial data) | CONFIRMED |
| V07 | Eöt-Wash α = 1/3: λ = 54.0 ± 0.8 μm, m_q ≥ 3.65 meV (A.4) | 53.68 raw → 54.0 μm; 3.65 (3.60–3.71) | 53.99 raw → 54.02 μm (α = 1 crossing at 38.575 against 38.6 stated); **3.653 (3.62–3.68) meV** | CONFIRMED |
| V08 | Cassini m_q ≥ 2.9e-15 eV; light-q excess 2.2e4σ | 2.9e-15 | 2.891e-15 (1.81e-15 at b = 1.6 R_sun); 21739 | CONFIRMED |
| V09 | max Ω_dyn, ρ_dyn, 1+w; values at the Cassini bound and at natural m_q (A.5) | 3.5e-62; 1.3e-108; 1.0e-61; 5.5e-38; 1.2e-121 | 3.459e-62; 1.2725e-108 GeV⁴; 1.010e-61; 5.523e-38; 1.167e-121 | CONFIRMED |
| V10 | K_lab = f'(0) to 1e-32 | 1e-32 | 1.8e-32 | CONFIRMED |
| V11 | gaps to H0; 14.5σ from G_cos/G_lab = 3/4 (A.6) | 2.5e30, 2.0e18; 14.5σ | 2.54e30, 2.01e18; 14.5σ (14.3σ with the unrounded ledger 1.01123 ± 0.01830) | CONFIRMED; the 14.5σ is qualified (§3.4) |
| V12 | φMDE runs (A.6 table) | see bound | m = 0.1: 0.001987, 1.4303, 3.06e-4, p = 0.500; m = 10: 0.03077, 28.117, 1.23e-5, 0.503; m = 100: 0.1076, 146.82, 1.57e-6, 0.704; m = 1 sits on another calibration branch: 0.0111, 5.5 (bound 0.013, 7.0) | CONFIRMED |
| V13 | ω_BD = 0 dust era a ∝ t^{1/2} (APT07) | [R] | exact solution: a ∝ t^{1/2}, K ∝ t^{1/2}, R = 0 | CONFIRMED |
| V14 | "any m_q ≲ H(z_eq) ≈ 1.6e5 H0" (A5 verdict) | without [E] in the table | H(z_eq)/H0 = 1.574e5; the runs stop at 1e3 H0 | CONFIRMED as [E]; the label needs adding |
| V15 | ρ_V = (3/4)m²K(δK/K)², δK/K = 2√Ω_Λ H0/m (A.6) | [D] | identical | CONFIRMED |
| V16 | relic toy ≈ 0.10 κ² ln²(H_i/m) (A.7) | 0.091–0.102 | 0.099–0.105 | CONFIRMED |
| V17 | κ_QCD = 1.7e-3, T_osc = 1.6 TeV, τ < 1 s needs m ≳ 1.3e5 GeV | 1.68e-3; 1612.7; 1.33e5 | 1.68e-3, with (ε−3p)/T⁴ = (49/6)α_s² derived for n_f = 6; 1612.3; 1.33e5 | CONFIRMED |
| V18 | relic/radiation growth T_osc → today; overproduction (A.7) | 2.3e15; Ω_relic 6e4–2e8; ×2e5–7e8 of Ω_c | **7.98e15; Ω_relic 2.0e5–6.0e8; ×7.8e5–2.3e9** | **CORRECTED** (§3.1) |
| V19 | f(R) map, f'' and the affine inverse for any C² f (B.1) | [D] | identical. The q-stationarity factorizes as K'(q)(R − λ(q−q0)), so the root is unique where K' ≠ 0 | CONFIRMED |
| V20 | P_q ⟺ f(0) = 0, f'(0) > 0, f''(0) > 0 (B.2) | [D] | identical | CONFIRMED |
| V21 | HS n=1: Minkowski ghost; R_*/R_1; R_1; stability (B.3) | 0.090/0.011/0.0011; 25.12/26.09/26.10; 40/4.0e3/4.0e5 | identical in every digit shown | CONFIRMED. Novelty: a narrow arXiv abstract search found nothing earlier, which does not prove novelty |
| V22 | Starobinsky n=2: f''(0) < 0, x1 = 1.637, 3.959, stable | as stated | identical; Starobinsky's own stability inequality and λ(x1) relation are reproduced | CONFIRMED |
| V23 | AB: G_Mink/G = 1 + e^{2b} ≳ 1e59–1e122 (B.3) | 1.3e122 / 4e59 | reproduced as readings of AB's conditions (0705.3199 §4 TeV: "F''(R0) ≪ 1", "F''(R0) ≪ F'(R0)x²"). **From the exposed \|f_R0\| bounds: ≳ 5.4e3 (DESI FS), 7.9e3 (clusters), 2.4e6 (galaxies)** | **CORRECTED** as a lower bound (§3.2) |
| V24 | threshold δf = 2c(q_R−q0)/q0, δR = 2c/(q0K'), δR/R_dS = 2.4e54 (B.5) | 4.0e-29; 1.7e-83; 2.4e54 | sympy to R³. Exact linear-K form f_new(R) = f(R+d) − d f'(0) + 2(μ'q0 − c), d = 2μ'/K1, holds to 1e-15. Values 4.01e-29, 1.700e-83, 2.36e54 | CONFIRMED |
| V25 | "affine gauge: holds automatically, K'q0/K ~ 1e113" (B.5) | 1e113 (not in code) | δK/K = 2c/(λq0K); survival needs λq0 ≳ 2Δ/K = 2.4e54 R_dS | **CORRECTED** (§3.3) |
| V26 | HS background, CPL table, Ω_m(true), phantom crossing, attractor (B.6) | (−0.776,−0.295) … (−0.999668,−0.000354); z_x = 1.59–1.60 | identical to ≤ 1e-6; z_x = 1.588, 1.592; R(200)/R_1 − 1 = 5.1e-8, 1.4e-8 | CONFIRMED |
| V27 | (1+w0, wa) ≈ (3.32, −3.53)\|f_R0\| | (3.32, −3.53) | (3.318, −3.536) with bound's weighting (uniform in ln a); uniform in a: (2.82, −2.03); uniform in z: (3.96, −5.17) | CONFIRMED; depends on the weighting |
| V28 | values at the published f_R0 bounds; Compton λ | e.g. 8.5e-5; 12.3 Mpc/h | 8.53e-5; 12.30 Mpc/h (18.2 Mpc) | CONFIRMED |
| V29 | linear growth table; background shift ≤ 1e-4 | e.g. +6.8% | identical; shift on the true background −6.1e-5 (\|f_R0\| = 1e-4, k = 0.1) | CONFIRMED |
| V30 | DESI DR2 diagonal distances of the ΛCDM point | 4.1/4.2/5.3/5.9 | **4.1/4.1/5.2/5.7**, using the wa error on the model's side | **CORRECTED** (minor) |
| V31 | screened QV2 viable only as a ΛCDM mimic (B.6) | [C][A] | agrees | CONFIRMED; **missing credit**: Wang-Hui-Khoury 2012 |
| V32 | stage-1 updates (D) | eq. §6.4, ledger row 6 | ledger row 8's "9e-39 ρ_Λ at the Cassini bound" is ρ_V; the effective DE is 8.1e-38 ρ_Λ | **CORRECTED** (un-noted; negligible) |
| V33 | quotes from KV 0807.3896 §5–6, Starobinsky, AB, HS, Lee+ | [Q] | all found at the places stated | CONFIRMED |
| V34 | AB vacuum structure | [Q, C] | middle vacuum unstable; R_1 stable; F'(0) = 1/(1+e^{2b}) | CONFIRMED |
| V35 | α = 1/3 excluded from 54 μm to the Cassini regime | [Q-figure] | below 1/3 from 54 μm to the 1 mm frame edge (own digitization); not re-read beyond | CONFIRMED to 1 mm; NOT CHECKED beyond |
| V36 | signed α (q-mode α = +1/3; the figure gives \|α\|) | Supplemental not obtained | not obtained | NOT CHECKED |
| V37 | retrospective/held-out labelling | correct | correct; the "check at use" caveat matters, since Euclid DR1 is due near this date | CONFIRMED |

## 2. Derivations checked [D][C] (`code/v1_symbolic.py` → `results/v1_symbolic.json`)

### 2.1 The f(R) 00 equation (V1a)

Vary S = ∫N a³[f(R)/2 − ρ(a)]dt with respect to the lapse. Take R(N) = 6[ä/(aN²) + ȧ²/(a²N²) − ȧṄ/(aN³)] and
f = c0 + c1R + c2R² + c3R³ + c4 sin(c5R) as a generic test function. At N = 1 the result equals a³ × [3f'H² − (Rf' − f)/2
+ 3Hḟ' − ρ], with ratio exactly 1. This is the equation bound uses in A.1 and in its HS cross-check.

### 2.2 q-theory minisuperspace

From S = ∫N a³[K(q)R/2 − U(q) − ρ]:
- the constraint 3H²K + 3HK̇ = ρ + U;
- the q-equation K'R = 2U';
- the trace equation K̈ + 3HK̇ = (ρ_d + 4U − KR)/3. The trace equation follows from the a-equation up to (1/3) × constraint,
  for dust and for radiation.

The signs agree with stage-1 equations §1–3.

### 2.3 Adiabatic relation (V1b)

The exact identity 3K_labH² − ρ = Λ + (RF' − F)/2 − 3F'H² − 3HF''Ṙ holds to all orders in F. With F = KR²/(6m²) on
the GR + Λ background, where ρ̇_d = −3Hρ_d, ρ̇_r = −4Hρ_r and R = (ρ_d + 4Λ)/K:

    ρ_dyn = [(9/4)ρ_d² + 2ρ_dρ_r − 4Λρ_r]/(3Km²)

The Λρ_d and Λ² terms cancel identically. This gives Ω_dyn = 9H0²Ω_m²/(4m²) and w = 1. The gravitating ρ_V from the q-variables
is K1²R²/(8e2) = K0R²/(12m²), so the ratio is exactly 9 (V1c).

### 2.4 Static root (V1d)

The δμ series of the dS/AdS root, with Taylor data up to E4 and k3, reproduces bound's coefficient. It reproduces
stage-1 §5.2 through c*(δμ) = q0δμ + δμ²/(2eps''). There is no conflict with stage 1.

### 2.5 Map, inverse and P_q (V1e)

Implicit differentiation gives f' = K and f'' = K'²/(2eps'' − RK''). In the affine inverse the residual of the q-equation
is 0 and f is recovered exactly.

### 2.6 Threshold (V1f, plus v5 B4)

Explicit Taylor data (e2, e3, K0, K1, K2) were used. The P_q re-selection has dμ0/dc = 1/q0, and the O(c) change of f
equals 2c(q_R − q0)/q0 through R³.

In the linear-K realization I derived the **exact** relation

    f_new(R) = f(R + d) − d f'(0) + 2(μ'q0 − c),    d = 2μ'/K1.

It was confirmed numerically to 1e-15 on an HS-type f with c1 < 1. Since μ'q0 − c = O(c²), bound's shift δR = 2c/(q0K1)
follows with an O(c) relative error (3e-4, 3e-3, 3e-2 at c = 1e-4, 1e-3, 1e-2).

### 2.7 Other checks

- φMDE (V1g): with U = 0 and dust, a ∝ t^{1/2} and K = (2ρ0/3)t^{1/2} solve the constraint and the trace equation exactly.
- V1h: ρ_V = (3/4)m²K(δK/K)². The QCD interaction measure is (ε − 3p)/T⁴ = (1/3)(1 + 5n_f/12)(11 − 2n_f/3)α_s², which
  is 49/6 α_s² for n_f = 6. This is the "2940/360" in bound_numbers.py.
- V1i: AB F' = [1 + tanh(aR − b)]/2 and 1/F'(0) = 1 + e^{2b}.

## 3. Corrections in detail

### 3.1 A.7 relic growth factor (V18) [C]

`bound_numbers.py` sets growth = (T_osc/T0)(3.91/106.75)^{1/3}.

Correct derivation. Entropy conservation gives g_s a³T³ = const, and ρ_rad = (π²/30)g*T⁴. Then

    [ρ_relic/ρ_rad]_0 / [ρ_relic/ρ_rad]_osc = (T_osc/T0) · (g*_osc/g_s,osc) · (g_s0/g*0) = 1.163 · T_osc/T0

**Values.**
- Growth factor: 7.98e15, against bound's 2.28e15 (factor 3.50).
- Ω_relic, using bound's own onset fractions (2.82e-7 and 8.32e-4): **2.0e5** (one e-fold of drive) and **6.0e8** (drive from
  1e15 GeV), against bound's 6e4 and 2e8.
- Overproduction relative to Ω_c = 0.264: **7.8e5–2.3e9**, against bound's 2e5–7e8.

The A6/A.7 verdict (tuned misalignment or an extra decay channel is needed) is unchanged, and stronger. This is an [E]
that A.5 does not use.

### 3.2 Appleby-Battye G_Mink (V23) [Q][C][E]

The TeX of 0705.3199 §4 was checked. Bound transcribes AB's two conditions correctly:
- "F''(R0) ≪ 1" with R0 = 12H0² ≪ 1 in M = 1 units, i.e. 8be^{−2b} ≪ 2R0, gives b ≈ 140.6;
- "F''(R0) ≪ F'(R0)x²" with x = 1 mm gives b ≈ 68.7.

Both are AB's sufficient conditions for their Newtonian limit, expanded about R0. AB themselves say the expansion
should be about the local density, where F'' is exponentially smaller. Neither condition is required by data. The
observational requirement is on |f_R0|.

For the AB form, |f_R0| = 1/(1 + e^{2α0}) with α0 = aR_today − b = b(1 + Ω_m/(2Ω_Λ)) = 1.23b. Then

    G_Mink/G = 1 + e^{2b} ≈ |f_R0|^{−0.813}

This gives ≥ 5.4e3 (|f_R0| ≤ 2.6e-5), 7.9e3 (1.6e-5) and 2.4e6 (1.4e-8) [E: HS-calibrated bounds applied to AB].

**Corrected statement.** AB's Minkowski vacuum is healthy and has G_Mink/G ≳ 1e4–1e6 on exposed data. It is ≳ 1e59–1e122
only under AB's own stronger conditions. Bound's conclusion stands: the vacuum P_q describes has a Newton constant very
different from ours and is never reached.

### 3.3 B.5 affine-gauge sentence (V25) [D]

Bound writes "In the affine gauge this holds automatically: there K'q0/K ~ λq0|f_R0|/R0 ~ 1e113." No code or note
derives the 1e113.

In the affine gauge (q − q0 = R/λ) the threshold changes f by (2c/(λq0))R (V1f). This is a renormalization
δK/K = 2c/(λq0K), and it moves R_dS ∝ Λ_eff/K by the same relative amount. Survival therefore needs
λq0 ≳ 2Δ/K = 4.0e-29 GeV² = 2.4e54 R_dS. That is, q may change by at most ~4e-55 q0 across the feature.

This is exactly the "feature of relative q-width ≲ ρ_Λ/Δ" that bound states in the next sentence. The affine gauge is
not an escape. **Corrected:** "In every realization the feature survives only if its q-width is ≲ ρ_Λ/Δ relative to q0.
In the affine gauge this reads λq0 ≳ 2.4e54 R_dS." The B4/B5 verdict (q-theory protects f(0) = 0, not the meV scale)
is unchanged.

### 3.4 Smaller items

- **DESI distances (V30).** Bound averages the asymmetric wa errors. Using the side toward wa = 0 (+0.22, +0.31, +0.23)
  gives 4.09/4.08/5.17/5.74. The statistic is labelled indicative in any case.
- **Ledger row 8 (V32).** The row's "ρ_V shift today 9e-39 ρ_Λ even at the Cassini bound" is the gravitating ρ_V. The
  w-relevant effective DE is 9 × that = 8.1e-38 ρ_Λ (5.5e-38 ρ_c). Section D of bound should list this beside the
  equations §6.4 update. Both are negligible.
- **G_cos/G_lab = 3/4 against Lamine+ (V11, qualified).** For a light q the background (Friedmann) G is G* and
  sub-horizon perturbations see G_eff = (4/3)G* = G_lab. Lamine+'s λ_G rescales G uniformly, so "14.5σ" is not a
  like-for-like test. Regime (ii) is excluded regardless, by Cassini (2.2e4σ) and by the φMDE runs.
- **CPL slope (V27).** At |f_R0| = 0.1 the fitted wa runs from −0.26 to −0.33 across fit weightings, while the local slope
  at z = 0 has the opposite sign (wa_loc = +0.011). "(3.3, −3.5)|f_R0|" is specific to uniform-in-ln a weighting over
  z ≤ 2. The conclusions do not depend on it.
- **Error coefficient (V05).** A fourth family, F4 (eps = μ0q + e[cosh(q−1) − 1], K = 1/(1 − 0.5(q−1))), gives
  c = −40.0 (m = 30) to −41.7 (m = 300), so the formula can over- or under-predict. "|c| ≈ 15–53 in four families" is the accurate wording.
- **f_R0 label.** In the |f_R0| = 0.1 row the actual f_R today is −0.0935, because the HS relation for c1/c2² is
  approximate. The row is labelled by its nominal value. This is immaterial.

## 4. Missing credit

**Wang, Hui, Khoury**, "No-Go Theorems for Generalized Chameleon Field Theories", PRL 109, 241301 (2012),
arXiv:1208.4612 (abstract read via the arXiv API). Under general conditions they prove two things:
- the present Compton wavelength is ≲ Mpc;
- there is no self-acceleration. Chameleon-type cosmic acceleration "rel[ies] on a form of dark energy".

B.6's conclusion, that the screened QV2 is viable only as a ΛCDM mimic equivalent to QV1-cal plus an invisible
parameter, is a special case and should cite it. The HS-specific numbers in B.6 remain bound's own computation.

Everything else that bound credits was checked at the stated places:
- KV 0807.3896 §5–6 (scratchpad text, lines ~746–775);
- Starobinsky 0706.2041 ("should be zero due to some other symmetry"; scalaron overproduction);
- AB 0705.3199 ("we must still fine tune"; conditions of §4);
- HS 0705.1158 (lim_{R→0} f = 0; freezing; phantom crossing; eq. fR0bound);
- Lee+ 2020 (38.6 μm, 95%/2σ).

## 5. Numerics: convergence and independence [C]

- **v2** (`v2_frw_adiabatic.py`).
  - H² comes from the constraint, so the constraint holds exactly.
  - The q-equation is not imposed and is violated by ≤ 7.4e-10 relative.
  - The slow-manifold fixed-point iteration converges to ≤ 1e-11 in q within 4–5 steps.
  - Averaging uses a Hann window over 2 periods, with residual oscillation ≤ 1e-3 of ⟨P⟩ at m ≥ 100.
  - The ratio agrees with bound's one-period top-hat average to ≤ 2e-6 at m ≥ 100 and to ≤ 4e-4 at m = 30, where
    initial data and window matter at O((H/m)⁴). At m = 10 the slow-manifold iteration
    diverges, as expected for an asymptotic series at H/m ≈ 0.1, so one step is used there.
- **v4.** At m = 0.1, 10 and 100 the φMDE runs reproduce bound's numbers to ≤ 1e-9 relative. This shows that the
  formulations are equivalent and the integrations converged.
- **v5.** Starts at z = 10 and z = 20 agree to ≤ 1e-6 in w0. The HS background agrees with bound to ≤ 1e-6.

## 6. Inputs and provenance

| input | value | provenance |
|---|---|---|
| H0, Ω_m, Ω_Λ, Ω_r | 67.4 km/s/Mpc, 0.315, 0.685, 9.1e-5 | as in bound (Planck 2018 via ledger row 9; Ω_r [R]) |
| M, ρ_Λ | 2.435e18 GeV, 2.52e-47 GeV⁴ | BRIEF.md |
| ħ, ħc, k_B, Mpc, R_sun, T0 | CODATA / IAU values | [R] |
| g*(T_osc), g_s0, g*0 | 106.75, 3.909, 3.363 | [R] (SM; massless ν with T_ν/T = (4/11)^{1/3}) |
| α_s(1.6 TeV) | 0.085 | bound's value [R] |
| Eöt-Wash curve | own digitization of fig5b1.pdf (arXiv:2002.11761 source) | [C from Q-figure]; read level: figure vector data plus the constraint paragraph (TeX) |
| Cassini, Lamine+, DESI DR2, Ω_ch² | ledger rows 7, 5, 8, 13 | [Q via ledger] |
| f_R0 bounds | as in bound B.6 | [Q via bound]; not re-read |
| m_H, v (Δ_EW) | 125.25, 246.22 GeV | [R] (equations §7) |
| AB conditions | 0705.3199 §4 | TeX read (own download) |
| (eps, K) families F1–F3 (re-implemented from bound's NOTES), F4 (new), HS-type f with c1 = 0.5 | illustrative | [P] |

**Calibrations.** These are the same as in bound and are named as such:
- δμ is calibrated to Ω_Λ = 0.685 (the static root);
- HS c1/c2 = 6Ω_Λ/Ω_m = 13.05.

No number is presented as a prediction.

## 7. Limits

- The Supplemental Material of Lee+ 2020, which has the signed α limits, was not obtained by bound or by me (V36).
- The long-range exclusion figures (hep-ph/0307284) were not re-read beyond the 1 mm edge of the Eöt-Wash figure (V35).
- The f_R0 bounds were not re-read. Neither were Lamine+ or the DESI numbers beyond the ledger.
- The novelty of the HS Minkowski ghost was searched only at abstract level, with two arXiv API queries.
- The A.7 relic estimate remains an order-of-magnitude [E], even after correction (SM trace history, single-channel
  decay width).
- The regime-(ii) threshold m_q ≈ H(z_eq) was not run beyond 1e3 H0, by bound or by me.
- My m = 1 light-q run used a different de Sitter calibration branch. Only qualitative agreement is claimed there.

## 8. Files and reproduction

| file | does | runtime |
|---|---|---|
| `code/v1_symbolic.py` → `results/v1_symbolic.json` | §2 | ~10 s |
| `code/v2_frw_adiabatic.py` → `results/v2_frw_adiabatic.json` | V05, V06 | ~30 s |
| `code/v3_eotwash_digitize.py` → `results/v3_eotwash_digitize.json` | V07, V35 | — |
| `code/v3_numbers.py` → `results/v3_numbers.json` | V07–V11, V14, V17, V18, V23, V24 numbers | < 1 s |
| `code/v4_light_and_relic.py` → `results/v4_light_and_relic.json` | V12, V16 | ~50 s |
| `code/v5_fR_models.py` → `results/v5_fR_models.json` | V21–V30, V34 | ~30 s |
| `results/verdicts.json` | all verdicts with sources | — |

`v3_eotwash_digitize.py` takes an argument: the SVG made from the arXiv figure. `v3_numbers.py` reads its JSON
output.

    cd research/reality_framework/grut_q_vacuum_construction_20261021/tracks/v_bound/code
    # figure (scratchpad only): curl -L https://arxiv.org/e-print/2002.11761 | tar xz -C $SCRATCH/v_bound_dl/src_2002.11761
    #                           pdftocairo -svg $SCRATCH/v_bound_dl/src_2002.11761/fig5b1.pdf $SCRATCH/fig5b1.svg
    python3 -I v3_eotwash_digitize.py $SCRATCH/fig5b1.svg > ../results/v3_eotwash_digitize.json
    python3 v1_symbolic.py && python3 v2_frw_adiabatic.py && python3 v3_numbers.py && python3 v4_light_and_relic.py && python3 v5_fR_models.py
