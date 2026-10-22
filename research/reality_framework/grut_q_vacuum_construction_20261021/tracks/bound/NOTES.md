# Track `bound`: what any q-sector remnant can contribute, and screening as disappearing-Λ f(R)

21 October 2026. Stage-2 track of the construction pass (`../../BRIEF.md`, `../../BRIEF_STAGE2.md`).

**Labels.** [D] derived here; [C] computed by `code/` (outputs in `results/`); [Q] stated by the cited source at the read
level of §F; [E] estimate; [A] argued; [R] recalled, not checked here.

**Scope of the claims.** No new theory, no empirical validation, no confirmatory test. Every observational comparison
is **retrospective** (the data are exposed). Where a future measurement is named, it is named as held out. The
screened version of part B is a **published model class** (Starobinsky 2007; Hu-Sawicki 2007; Appleby-Battye 2007).
Writing it in q-theory variables does not make it new, and Klinkhamer-Volovik themselves made the identification and
discussed these models (§B.1, §B.3).

**Sibling records used.** `../equations/NOTES.md` (cited "equations §n"), `../literature/NOTES.md` ("literature §n"),
`../ledger/NOTES.md` and `../ledger/results/ledger.json` ("ledger row n").

**Conventions.** As in BRIEF.md: signature (−,+,+,+), K = M² = 1/(8πG), M = 2.435e18 GeV. At fixed μ the theory is
f(R) gravity with S = ∫√−g f(R)/2 (equations §2.2). H0 = 67.4 km/s/Mpc = 1.438e-33 eV, Ω_m = 0.315 (Planck 2018,
ledger row 9) unless stated.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdicts in one screen

| # | Result | Label |
|---|---|---|
| A1 | For any eps(q), K(q) and any μ, in the adiabatic regime (m_q ≫ H) the q-sector contributes ρ_DE(t0) = Λ_cal(μ) + ρ_dyn. Here Λ_cal(μ) = −q0(μ−μ0) + [−1/(2eps'') + (4/3)q0²/(K m_q²)](μ−μ0)² + … is the μ-offset, i.e. the calibration. The dynamically generated part is ρ_dyn = [(9/4)ρ_d² + 2ρ_dρ_r − 4Λρ_r]/(3K m_q²), so **Ω_dyn(t0) = (9/4) Ω_m² (H0/m_q)²**, independent of Λ and of the shapes of eps and K at this order. It is a **stiff** component (w = +1, ∝ a⁻⁶), not dark energy. | [D][C] |
| A2 | Full FRW integrations of three (eps, K) families reproduce A1. The relative error is c(H/m)², with c = 15, 53 and 27 for the three families: 0.02–0.06% at m = 300 H0, 0.15–0.5% at m = 100 H0, and the formula fails below m ≈ 10 H0. The constraint is held to 1e-12. | [C] |
| A3 | An unscreened q-mode has α = 1/3. The Eöt-Wash 2020 95% exclusion curve crosses α = 1/3 at λ = 54.0 ± 0.8 μm (digitized from the source figure). Hence **m_q ≥ 3.65 meV** (3.60–3.71). α = 1/3 is also excluded at every longer range up to the Cassini regime. | [Q-figure][C] |
| A4 | **Maximal non-calibrated Ω_q(t0) of an unscreened q-mode: (9/4)Ω_m²(H0/m_q,min)² = 3.5e-62** (ρ = 1.3e-108 GeV⁴), with 1 + w ≤ 1.0e-61. It is 5.5e-38 at the Cassini-only bound and 1.2e-121 at natural m_q. | [C] |
| A5 | Regime (ii) (Hubble-frozen thawing q) needs m_q ≲ H0 = 1.4e-33 eV. That is 2.5e30 below the Eöt-Wash bound and 2.0e18 below the Cassini bound. Locally a light q has γ = 1/2 (excess 2.2e4 σ over Cassini) and G_cos/G_lab = 3/4 (14.5 σ from the Lamine+ 2025 measurement). Within the theory itself, the ω_BD = 0 coupling to the dust trace drives K by factors 9–500 between z = 1e4 and today (m_q ≤ 100 H0 in the runs), for any m_q ≲ H(z_eq) ≈ 1.6e5 H0, and turns the matter era into a ∝ t^{1/2} (APT07's φMDE). No viable regime-(ii) branch exists. The DE such a branch would carry is in any case an initial displacement δK/K ≈ 2√Ω_Λ H0/m_q, i.e. calibrated initial data. | [D][C][Q] |
| A6 | Regime (iii) relics are dust (w = 0), not DE. Their amount is state data. A toy run gives a relic fraction at onset ≈ 0.10 κ² ln²(H_i/m), where κ is the trace fraction of the plasma. For m_q = 3.65 meV this implies overproduction by ≈2e5–7e8 relative to Ω_c unless the misalignment is tuned (a moduli/scalaron-overproduction-type problem; Starobinsky 2007 §2, §4). | [C] toy, [E] SM numbers |
| B1 | Map: f(R) = K(q_R)R − 2eps(q_R) + 2μq_R, with f' = K, f'' = K'²/[2(eps'' − K''R/2)], ∂f/∂μ = 2q_R and f(0) = 2eps*(μ). Inverse: **every C² f is realizable**, e.g. K(q) = f'(λ(q−q0)), eps = [K R(q) − f(R(q))]/2 + μq. The realization is unique only up to reparametrizations of q and the shift (eps + d·q, μ + d). Credit: KV 0806.2805 §II, 0807.3896 §4. | [D][C][Q] |
| B2 | **P_q in f(R) language: f(0) = 0 with f'(0) > 0 and f''(0) > 0** for the full final f. This is "Minkowski is a stable vacuum solution" (checked: R f' − 2f = 0 at R = 0 ⟺ f(0) = 0). The final-state reading of P_q needs in addition that the cosmological attractor be R → 0. f(0) = 0 alone does not imply it. | [D] |
| B3 | All three published models have f(0) = 0, so they satisfy P_q **as a label** (its stability part, f'(0), f''(0) > 0, only in Appleby-Battye). All three end in a **stable de Sitter attractor** R_1 ≈ 4Λ_eff, not in Minkowski. Their Minkowski point: Starobinsky n=2 has f''(0) < 0, tachyonic, as the source states. Hu-Sawicki n=1 has f'(0) = 1 − c1 < 0, a **ghost**, for every \|f_R0\| ≲ 0.2 (c1 ≈ 0.2/\|f_R0\|), with an f' = 0 singular point between Minkowski and de Sitter. Hu-Sawicki n=2 is tachyonic. Appleby-Battye has a healthy Minkowski vacuum, but there G_Mink/G = 1 + e^{2b} ≳ 1e59–1e122. So **P_q's final-state reading fails in all three**. | [Q][D][C] |
| B4 | QV1 vs QV2: a QV2 built on these models must replace P_q's final-state reading by the label condition f(0) = 0. It must also put a calibrated meV⁴ feature into eps(q). Its observed Λ is that calibrated feature. Gate Q3: under P_q re-selection a threshold c shifts the curvature structure by δR = 2c/(q0K'). With natural q0K' ~ K this is δR/R_dS = 2.4e54 for Δ_EW, so the feature survives only if q0K'/K ≳ 2.4e54 across it. | [D][C] |
| B5 | Hu-Sawicki n=1 solved. The calibrated number is c1/c2 = 6Ω_Λ/Ω_m = 13.05 (the effective Λ, calibrated to Ω_m = 0.315). For \|f_R0\| ≤ 1e-3, (1 + w0, wa)_CPL ≈ (3.3, −3.5)·\|f_R0\|. At the published bounds this gives \|1+w0\| ≤ 8.5e-5 (DESI DR1 full shape) and ≤ 4.6e-8 (galaxy morphology), i.e. the ΛCDM point. It shares ΛCDM's retrospective DESI DR2 tension (2.8–4.2σ). Approaching the DESI direction needs \|f_R0\| ~ 0.1, which gives (w0, wa) = (−0.78, −0.30) and is excluded by 4 (clusters, full shape) to 7 (galaxies) orders of magnitude. | [C][Q] |
| B6 | Viability: on exposed data the screened QV2 survives only as an f_R0 ≲ 1e-5 (or ≲ 3e-8) ΛCDM-mimic. It is observationally equivalent to QV1-cal (one calibrated Λ) plus a parameter bounded to be invisible. Held out: a ≥ 5σ detection of w ≠ −1 by data not yet released would falsify it, since it predicts \|1+w0\|, \|wa\| ≲ 1e-4. Its growth signal (≤ 0.05% in fσ8 at k ≤ 0.2 h/Mpc for \|f_R0\| ≤ 1.4e-8) is below Stage-IV reach. Neither test probes the q-content (μ). | [C][A] |
| B7 | GRUT contribution in this track: none. Every element used is KV, the cited f(R) literature, or standard scalar-tensor physics. | [A] |

---------------------------------------------------------------------------------------------------------------------

## A. Class-wide relation for an unscreened q-sector remnant

### A.1 Setting and the quantity bounded [D]

Fix μ, so that the theory is f(R) gravity (equations §2.2; adiabatic_relation.py S4). Write
f(R) = K_lab R − 2Λ + F(R), where:
- K_lab = f'(R_lab) is the coupling a laboratory measures;
- Λ is the constant part;
- F collects the q-dynamics.

For an unscreened q-mode, R_lab ≪ m_q², so K_lab = f'(0) = K(q_{R=0}(μ)) to 1e-32 (ρ_lab/(K m_q²) at m_q = 3.65 meV).
The effective dark-energy density an observer infers with G_lab is

    ρ_DE,eff ≡ 3 K_lab H² − ρ_m = Λ + (R F' − F)/2 − 3 F' H² − 3 H F'' Ṙ          (exact in F for the 00 equation)

Evaluated on the GR + Λ background, this is first order in F, i.e. in 1/m_q² (sympy, S1). It is also the formula used in
§B.6 for Hu-Sawicki, where it agrees with the full integration to 4e-6 in y_H(0) at |f_R0| = 1e-3. A perturbative
treatment of Hu-Sawicki-type models in a deviation parameter is due to Basilakos, Nesseris, Perivolaropoulos (1302.6051,
abstract read) [Q]. The formula here was derived independently and is not claimed to be new.

### A.2 Regime (i): adiabatic, m_q ≫ H [D]

The q-mode follows the slow manifold. To first order, F = f2R²/2 with f2 = f''(0) = K'²/(2eps'') = K/(3m_q²) (S2: the
same from the Taylor data of eps and K). Substituting R = (ρ_d + 4Λ)/K, H² = (ρ_d + ρ_r + Λ)/(3K) and Ṙ = −3Hρ_d/K:

    ρ_dyn ≡ ρ_DE,eff − Λ_cal = [ (9/4) ρ_d² + 2 ρ_d ρ_r − 4 Λ ρ_r ] / (3 K m_q²)   + O(H⁴/m_q⁴)
    Ω_dyn(t0) = (9/4) Ω_m² (H0/m_q)²            (radiation term < 1e-4 of it today)
    w_dyn = +1 in the dust era (ρ_dyn ∝ a⁻⁶)

The Λ dependence cancels exactly at this order (S1: ∂ρ_dyn/∂Λ = 0 without radiation). The coefficient depends on
nothing but m_q, because at first order only f''(0) enters, and that is m_q. This holds for **any** eps and K with
K > 0, eps'' > 0, K' ≠ 0, and for any μ.

**The μ-offset part (calibration).** The static (dS/AdS) root as a series in δμ = μ − μ0 (S3) is

    ρ_V,s = −q0 δμ + [ −1/(2 eps'') + (4/3) q0²/(K m_q²) ] δμ² + O(δμ³)

Λ_cal ≡ 3K_lab H_s² differs from ρ_V,s only at O(δμ²) (normalization). This is the calibration of QV1-cal: δμ is chosen
from ρ_Λ,obs.

**Relation to stage 1.** Equations §6.4 quoted the gravitating vacuum energy ρ_V ≈ ρ_m²/(12 K m_q²). The effective
dark energy that enters the Friedmann equation with G_lab also contains the −3HK̇ and the 3(K − K_lab)H² terms, and is
(3/4)ρ_m²/(K m_q²), **9 times larger**. The two are different quantities, not an inconsistency. Equations §6.4 also says
"raising it to ρ_Λ needs m_q ≈ 0.2 H0". With the effective quantity the formal value is m_q = 0.57 H0. That point lies
outside the regime of the formula, and the component would have w = +1 in any case.

### A.3 Numerical verification [C] (`code/frw_regimes.py`, `results/frw_regimes.json` part I)

**Families [P].** Units: H0 ≈ 1, K(q0) = 1.
- F1: eps = (e2/2)(1+q²), K = 1 + 0.4(q−1). Here f is exactly quadratic.
- F2: the KV 0806.2805 form, eps = (E/2)(−q² + q⁴/3), K = q (literature §4.1).
- F3: eps = e2x²/2 + e3x³/6 + e4x⁴/24 + μ0q with e3 = 3e2, e4 = 10e2; K = e^{0.7x}, x = q − 1.

**Procedure.**
- e2 (or E) is set from m_q.
- δμ is calibrated so that the static root gives Ω_Λ = 0.685.
- Full equations (equations §3.2) with dust and radiation.
- Second-order slow-manifold initial data, restarted at a_e/1.5 for each comparison epoch a_e.
- One-period average of Δ = 3K_labH² − ρ_d − ρ_r − Λ_cal, compared with ρ_dyn of A.2.
- The maximum relative constraint violation is ≤ 4e-12.

Ratio ⟨Δ⟩/⟨ρ_dyn⟩ at a = 1:

| family | m/H0 = 10 | 30 | 100 | 300 | c = (1 − ratio)/(H/m)² |
|---|---|---|---|---|---|
| F1 | 0.832 | 0.9832 | 0.99848 | 0.99983 | 15.2 |
| F2 (KV) | 0.595 | 0.9430 | 0.99471 | 0.99941 | 53 |
| F3 | 0.741 | 0.9700 | 0.99727 | 0.99970 | 27 |

At m/H0 = 300: Ω_dyn measured 2.4798e-6 (F1), 2.4789e-6 (F2), 2.4795e-6 (F3), against the formula (9/4)Ω_m²(H/m)² =
2.4806e-6. The error scales as (H/m)², as it should for a first-order result. The F2/F3 coefficients are larger
because f''' ≠ 0. For physical m_q the relative error is c(H0/m_q)² < 1e-58. The formula is not usable below
m ≈ 10 H0 (errors > 15%), which is the overlap region with (ii) and (iii) (A.6, A.7).

### A.4 Fifth-force inputs (exposed, retrospective) [Q][C] (`code/bound_numbers.py`)

The unscreened coupling is fixed: G_eff(r) = G_*(1 + e^{−m r}/3), so α = 1/3 and γ(r) = (1 − e^{−mr}/3)/(1 + e^{−mr}/3)
(equations §4.3). There is no free coupling, because ω_BD = 0.

**Eöt-Wash at α = 1/3.**
- Source: Lee, Adelberger, Cook, Fleischer, Heckel, PRL 124, 101101 (2020), arXiv:2002.11761, Fig. 5 bottom, 95% CL.
- Method: the exclusion curve was digitized from the vector paths of the arXiv source figure (scratchpad, not in the
  repository).
- Calibration: the same extraction puts the α = 1 crossing at 38.37 μm (paper text: 38.6 μm [Q]). For the Eöt-Wash
  2007 curve it gives 55.2 μm (Kapner+ 2007 abstract: 56 μm [Q]). The digitization uncertainty is 1.4%.
- **α = 1/3 crossing: 53.7 μm digitized, 54.0 ± 0.8 μm after rescaling to the stated α = 1 point.**
- Other curves at α = 1/3: Eöt-Wash 2007 gives 72.9 μm and the HUST curve 57.2 μm. Eöt-Wash 2020 is the strongest.
- **m_q ≥ ħc/λ = 3.65 meV (3.60–3.71 meV).** The ledger's "~4 meV [E]" is replaced by this value.

**Longer ranges.** α = 1/3 is excluded at every range from 54 μm up to the Cassini regime:
- 54 μm to 1 cm: Eöt-Wash 2004/2007 and Irvine, |α| < 0.1.
- ~4 mm to ~3e14 m: laboratory, geophysical, LAGEOS, LLR and planetary, |α| < 0.1.
- Source for both: Adelberger, Heckel, Nelson, Ann. Rev. Nucl. Part. Sci. 53, 77 (2003), hep-ph/0307284, figures
  "midlam" and "longlam" rendered from the arXiv source; figure-level read.
- λ ≫ AU: Cassini (below).

**Cassini.** γ − 1 = (2.1 ± 2.3)e-5 (Bertotti+ 2003, ledger row 7 [Q]). At the 2σ lower edge the bound is
m b ≥ 10.19, so m_q ≥ 2.9e-15 eV for b = R_sun. A light q (m r ≪ 1) has γ − 1 = −1/2, a 2.2e4 σ excess.

### A.5 The maximal dynamically generated Ω_q(t0) [C]

    Ω_dyn(t0) = (9/4) Ω_m² (H0/m_q)²  ≤  (9/4) Ω_m² (H0 λ_{1/3}/ħc)²  =  3.5e-62       (ρ_dyn ≤ 1.3e-108 GeV⁴)
    1 + w_DE(t0) = 2 Ω_dyn/Ω_Λ ≤ 1.0e-61

Inputs: H0 = 1.438e-33 eV, Ω_m = 0.315, λ_{1/3} = 54.0 μm.

| m_q | Ω_dyn(t0) |
|---|---|
| Eöt-Wash α = 1/3 bound, 3.65 meV | 3.5e-62 |
| Cassini-only bound, 2.9e-15 eV | 5.5e-38 |
| natural √(2/3) M | 1.2e-121 |

**Every non-negligible present-day vacuum energy of an unscreened q-sector in regime (i) is the μ-offset, i.e. the
calibration** (QV1-cal), and it has w = −1 to 1e-61.

### A.6 Regime (ii): Hubble-frozen, displaced q [D][C][Q]

**Mass needed.** The q-mode must be frozen today, so m_q ≲ H0 = 1.4e-33 eV. That is a factor 2.5e30 below the
Eöt-Wash bound and 2.0e18 below the Cassini bound.

**Local exclusion.**
- γ = 1/2: excess 2.2e4 σ.
- G_cos/G_lab = 3/4 against Lamine+ 2025, G/G_N = 1.011 ± 0.018 (ledger row 5): 14.5 σ.
- Ġ/G = O(H0) against LLR (ledger row 1).

**Cosmological exclusion within the theory** (`frw_regimes.json` part II). Family F1, Λ_cal calibrated, q started at
its vacuum (de Sitter) value at a = 1e-4. For m ≲ H(z_eq) the ω_BD = 0 coupling makes the dust trace drive K at O(1)
per e-fold. This is the φMDE of Amendola, Polarski, Tsujikawa (PRL 98, 131302 (2007), astro-ph/0603703 [R]; cited by
Starobinsky 2007 as APT07 [Q]).

| m/H0 | p = dln a/dln t at z = 100 | K(a=1e-4)/K(a=1) | H(a=1) (ΛCDM: 1) | Ω_m(a=1) |
|---|---|---|---|---|
| 0.1 | 0.500 | 0.0020 | 1.43 | 3e-4 |
| 1 | 0.500 | 0.013 | 7.0 | 8e-5 |
| 10 | 0.503 | 0.031 | 28 | 1e-5 |
| 100 | 0.70 | 0.11 | 147 | 2e-6 |
| 1000 | (oscillation-dominated) | 1.00 | 56 | 1e-4 |

Every case with m ≲ 1e3 H0 is grossly non-ΛCDM. The threshold is m_q ≈ H(z_eq) ≈ 1.6e5 H0 = 2.3e-28 eV [E].

**The amount is initial data.** For eps quadratic near q0, ρ_V = (3/4)m_q²K(δK/K)² [D]. To carry Ω_Λ the displacement
must be δK/K ≈ 2√Ω_Λ H0/m_q, which is O(1) for m_q ~ H0. That displacement is state data, so matching it to Ω_Λ is a
calibration.

**Verdict.** Regime (ii) is not a viable branch for an unscreened q-mode, for local and cosmological reasons alike
[C][Q]. It becomes available only if the scalar is heavy in dense environments and light cosmologically, which is
screening (part B).

### A.7 Regime (iii): oscillation relics [D][C][E]

Once H < m_q the displaced q oscillates. The relic is dust, with ⟨w⟩ = 0, so it contributes **nothing** to dark energy.
It adds to Ω_c, bounded by Planck Ω_ch² = 0.120 (ledger row 13).

**Sources of the amount.**
- Transitions: suppressed by (x/sinh x)², i.e. zero at natural m_q (equations §5.3).
- Initial data: state, uncomputed.
- The trace of the plasma before onset: computed here.

**Toy computation** (`frw_regimes.json` part III). One fluid with trace fraction κ, w_f = (1−κ)/3; q starts at its
vacuum value with H_i/m = 1e2 and 1e3. The relic fraction at onset (H = m) is

    E_osc/ρ_tot |_onset ≈ 0.10 κ² ln²(H_i/m)        (coefficient 0.091–0.102 for κ = 1e-3 … 1e-2)

i.e. log-enhanced misalignment of order κ, as for moduli [R: the cosmological moduli problem]. The physical scaling is
quadratic in κ [C].

**SM estimate at m_q = 3.65 meV [E].**
- Onset at T_osc ≈ 1.6 TeV.
- κ from the QCD trace anomaly ≈ 1.7e-3 (leading order, α_s = 0.085 [R]). The electroweak phase is symmetric there,
  so mass terms do not contribute.
- Relic-to-radiation growth to today: 2.3e15.
- If the relic is stable, Ω_relic ≈ 6e4 (drive over one e-fold) to 2e8 (drive from 1e15 GeV), against Ω_c = 0.26.
- With one scalar channel, Γ ~ m³/(192πM²) gives decay before 1 s only for m_q ≳ 1.3e5 GeV [E].

So an unscreened q-mode in the window 3.65 meV ≲ m_q ≲ 1e5 GeV also needs tuned initial misalignment or an extra decay
channel. This is the q-theory version of the scalaron-overproduction problem Starobinsky raises (0706.2041 §2 end, §4 [Q]).
Status: an estimate. It is not used in A.5, and its next dependency is the SM trace history and the decay widths.

### A.8 The class-wide relation (summary)

For the action of BRIEF.md with any eps(q), any K(q) (K > 0, eps'' > 0, K' ≠ 0) and any μ, with an unscreened q-mode:

    ρ_DE(t0) = Λ_cal(μ)  +  (9/4) Ω_m² (H0/m_q)² ρ_c [w = +1]  +  ρ_relic [w = 0]
    Λ_cal(μ) = −q0 (μ − μ0) + O((μ−μ0)²)          ← calibration (QV1-cal)
    m_q ≥ 3.65 meV (Eöt-Wash, α = 1/3, 95%)  ⇒  dynamical part ≤ 3.5e-62 ρ_c ;  regime (ii) requires m_q ≲ H0 (excluded)

**The only way an unscreened q-sector carries today's dark energy is the calibrated μ-offset. Nothing dynamical
survives.** [D][C]

---------------------------------------------------------------------------------------------------------------------

## B. Screening as the escape, and its identity with disappearing-Λ f(R)

### B.1 The map (eps, K, μ) → f(R) and its inverse [D][C] (`code/adiabatic_relation.py` S4)

**Forward map.** At fixed μ, eliminate q with the q-equation eps'(q) − K'(q)R/2 = μ:

    f(R) = K(q_R) R − 2 eps(q_R) + 2 μ q_R
    f'(R) = K(q_R)   (envelope)        f''(R) = K'²/[2(eps'' − K''R/2)]        ∂f/∂μ = 2 q_R
    f(0) = −2 [eps(q0) − μ q0]_{eps'(q0) = μ} = 2 eps*(μ)  (Legendre transform)

**Inverse (canonical, affine gauge).** For any C² f, take R(q) = λ(q − q0) and

    K(q) = f'(λ(q − q0)),   eps(q) = [K(q) λ(q−q0) − f(λ(q−q0))]/2 + μ q

The q-equation then holds identically, and the reconstructed f equals f (sympy residual 0). Any reparametrization R(q)
works too, e.g. linear K with R(q) = (f')⁻¹(K(q)) where f'' ≠ 0.

**Realizable f.** Every C² f(R) is realizable. Where f'' = 0 one has K' = 0, so the q-mode is not dynamical there; this
is Starobinsky's "weak singularity". K > 0 ⟺ f' > 0. The realization is not unique: reparametrizations of q, and
(eps + d q, μ + d), give the same f at a given μ. They differ only when μ changes, which is what §B.5 exploits.

**Credit.**
- KV 0806.2805 (end of §II, TeX read): expressing F in terms of R "gives an equation ... identical to the one of f(R)
  cosmology". q-theory "may give a microscopic justification for the phenomenological f(R) models" [Q].
- KV 0807.3896 §4: near equilibrium f̃ = −R + 16πGχ̃R² − 2λ with λ = −8πG F0 δμ (Weinberg sign convention), i.e. the
  μ-offset is the f(R) constant [Q].
- ω_BD = 0 ⟺ f(R): O'Hanlon 1972, Teyssandier-Tourrenc 1983, Sotiriou-Faraoni 0805.1726 [R].

What is derived here beyond those: the exact nonlinear inverse, the non-uniqueness, and §B.5.

### B.2 What P_q means in f(R) language [D]

P_q requires: eps_∞'(q0) = μ, eps_∞(q0) − μq0 = 0, eps_∞'' > 0, K(q0) > 0. With B.1 this is equivalent to

    f(0) = 0,  f'(0) > 0,  f''(0) > 0      (f = the full, final, renormalized f of the selected μ)

**Check (expected f(0) = 0 with Minkowski a vacuum solution).** Constant-curvature vacua solve R f' − 2f = 0. At R = 0
this is f(0) = 0. So **P_q ⟺ Minkowski is a ghost-free, non-tachyonic vacuum solution of the selected f.** ✓

The two conditions f'(0) > 0 and f''(0) > 0 are KV's own stability conditions (0807.3896 eq (26) [Q]).

**Two readings of P_q.**
- **Label reading (P_q^L):** μ is the value for which f(0) = 0. In f(R) language this is exactly Starobinsky's
  "disappearing cosmological constant", imposed there by hand ("should be zero due to some other symmetry",
  0706.2041 §1 [Q]). In q-theory the condition is carried by the state.
- **Final-state reading (QV1's P_q; KV 1102.3152 §6 via literature §3.2):** in addition, the actual cosmological
  attractor is R → 0. When U_μ has a single stationary point (QV1's natural eps) the two readings coincide. When f has
  further vacua R_1 > 0 (R_1 f'(R_1) = 2f(R_1)), they do not.

### B.3 The published models as q-theories [Q][D][C] (`code/hu_sawicki.py` B2)

| model (source, read level) | f(0) | Minkowski (P_q's equilibrium) | final attractor | source statement |
|---|---|---|---|---|
| Starobinsky 2007, f = R + λR0[(1+R²/R0²)^{−n} − 1] (0706.2041, full TeX) | 0 | f''(0) = −2nλ/R0 < 0: **tachyonic**. For n = 2, λ = 1, 2: m² = −R0/12, −R0/24 [C] | de Sitter R_1 = x_1R0, x_1 = 1.637 (λ=1), 3.959 (λ=2); stable [C] | "f''(0) is negative -- flat space-time is unstable"; the trajectory runs "to R = R_1 in the infinite future"; "flat space-time becomes unstable with a characteristic time of the order of the present Universe age" [Q] |
| Hu-Sawicki 2007 n = 1, f = R − m²c1x/(c2x+1) (0705.1158, sections) | 0 | f'(0) = 1 − c1: **ghost** for c1 > 1, i.e. \|f_R0\| ≲ 0.2. c1 = 20, 2.0e3, 2.0e5 for \|f_R0\| = 1e-2, 1e-4, 1e-6. f' = 0 at R_*/R_1 = 0.090, 0.011, 0.0011 [D][C] | de Sitter R_1 = 25.12, 26.09, 26.10 m² (= 0.963, 0.9996, 0.999996 × 4Λ_eff). Stability f'/(Rf'') = 40, 4.0e3, 4.0e5. Reached: R(a=200)/R_1 − 1 = 4e-8 (\|f_R0\|=1e-2), 1.4e-8 (1e-4) [C] | "the curvature freezes into a fixed value and ceases to decline" [Q]. The models are built so that lim_{R→0} f = 0 [Q]. The ghost at R = 0 is derived here; HS note a 1 + f_R < 0 regime only for the mCDTT models |
| Hu-Sawicki n = 2 | 0 | f''(0) = −2c1/m² < 0: **tachyonic** [C] | de Sitter (R_1 = 0.9997 × 4Λ_eff, \|f_R0\| = 1e-4) [C] | — |
| Appleby-Battye 2007, F = R/2 + (1/2a)log[cosh aR − tanh b sinh aR] (0705.3199, full TeX) | 0 | F'(0) = 1/(1+e^{2b}) > 0, F''(0) > 0: **healthy but G_Mink/G = 1 + e^{2b}**. Their requirement b ≫ 1 gives G_Mink/G ≳ 1.3e122 as written (8be^{−2b} ≪ 2R0, M = reduced Planck mass [A]) or ≳ 4e59 for F''(R0) ≪ (1 mm)² [C] | stable de Sitter minimum V = V0 > 0; Minkowski is the other minimum (V = 0); a third, unstable vacuum lies between [Q, C] | "two minima, one at V = 0 (Minkowski) and one at V0 > 0 ... stable de Sitter"; "we must still fine tune ... Λ = b/2a" [Q] |

**Answer to the central question.** In all three models the label part of P_q holds, since f(0) = 0. The stability part
(f'(0) > 0, f''(0) > 0) holds only in Appleby-Battye. The **final state is de
Sitter** with Λ_eff = K R_1/4 ≈ the calibrated scale, so **P_q's final-state reading is violated**. The Minkowski point
that P_q describes is never reached, and:
- in Starobinsky n=2 and Hu-Sawicki n=2 it is tachyonic;
- in Hu-Sawicki n=1 it is a ghost, separated from the cosmological branch by an f' = 0 (K = 0) singular point;
- in Appleby-Battye it is healthy but carries a Newton constant ≥ 1e59 times ours.

**KV reached the Starobinsky part of this first.** In 0807.3896 §5 (full text read): "For the particular models
suggested in [Starobinsky 2007], the functions f(R) expanded around R = 0 have negative M² (making Minkowski spacetime
unstable) ... can, in principle, be obtained by using special choices of ε and K. But such choices require a careful
tuning of the parameters, which is unnatural from our point-of-view". §6: the model "in which the cosmological constant
appears just at the latest stage does not follow from the macroscopic quantum-vacuum approach" [Q]. KV 0907.4887 §2
cites HS, AB, Starobinsky and Brax+ 2008 and locates "the crucial difference" in the integration constant μ [Q]. The
Hu-Sawicki ghost, the Appleby-Battye G_Mink and the threshold criterion of §B.5 are additions here.

### B.4 QV1 versus a QV2 built on them [D][A]

**QV1** (P_q in its final-state reading; natural eps with a single Minkowski equilibrium) **is kept unchanged**:
- ρ_V,∞ = 0;
- today's ρ_V is either zero or the calibrated QV1-cal offset;
- the dynamical part is ≤ 3.5e-62 ρ_c (A.5).

**QV2-screened** (new labelled version) = QV1 with two changes:
1. P_q weakened to P_q^L (f(0) = 0 for the final full f);
2. eps(q), K(q) given a feature at the curvature scale R ~ 4ρ_Λ/K, i.e. a meV⁴ structure in eps with a Hubble-light
   scalaron at cosmological densities, so that f is one of the published forms.

Consequences:
- **Its observed Λ is the calibrated height of that feature**: Λ_eff = m²c1/(2c2) in HS, λR0/2 in Starobinsky, b/2a in
  AB. The model is still one calibrated number, now placed in the shape of eps rather than in δμ. f(0) = 0 removes the
  bare constant, as P_q did in QV1. It does not produce the observed one.
- The **final state is de Sitter**, so QV2 is not a closure of "the equilibrium of the whole is pressure-free". The
  pressure-free equilibrium exists only as an unreached, and in HS/Starobinsky unstable, vacuum.

### B.5 Gate Q3 for QV2: threshold stability of the curvature feature [D][C]

Under eps → eps + c (any heavy threshold, any loop), P_q re-selects μ0 → μ0 + c/q0 (equations §5.1). From ∂f/∂μ = 2q_R:

    δf(R) = 2c (q_R − q0)/q0 + O(c²)

- **Linear K** (q_R − q0 = (f'(R) − f'(0))/K'): f_new(R) = f(R + δR) − δR f'(0), i.e. **the curvature structure moves
  by δR = 2c/(q0K')**.
- **Affine gauge** (q − q0 = R/λ): f_new = f + (2c/λq0)R, a renormalization of G only.
- General local statement: **the feature at R moves by δR = 2c/(q0K'(q_R))**.

Exact numerical check (Appleby-Battye form, compressed units; `hu_sawicki.json` B4): the P_q re-selection is solved
exactly. The first-order predictions hold to a relative 2e-5 (linear K, c = 1e-3), 2e-4 (c = 1e-2), 1e-5 (affine,
c = 1e-7) and 1e-4 (c = 1e-6), i.e. errors of O(c) as expected [C].

**Physical numbers** [C]:
- c = Δ_EW = 1.19e8 GeV⁴ and natural q0K' = K give δR = 4.0e-29 GeV², against R_dS = 4ρ_Λ/K = 1.7e-83 GeV², i.e.
  **δR/R_dS = 2.4e54**.
- The feature survives a threshold only if q0K'/K ≳ 2.4e54 across it, and more for heavier thresholds. In the affine
  gauge this holds automatically: there K'q0/K ~ λ q0 |f_R0|/R0 ~ 1e113.

So QV2 must abandon QV1's natural-scale assumption (q0K' ~ M², equations §7). Its q-functions need a feature of
relative q-width ≲ ρ_Λ/Δ. Otherwise every threshold retunes the dark-energy scale. This is the cosmological-constant
problem reappearing as the position of the curvature feature, as AB say of their own model ("we must still fine tune")
[Q]. The q-theory protects f(0) = 0, not the meV scale.

### B.6 Representative model solved: Hu-Sawicki n = 1 [C] (`code/hu_sawicki.py` B3; `results/hu_sawicki.json`)

**Setup.**
- Equations: HS's reduced Friedmann system (y_H, y_R) (0705.1158 §II.B [Q]), solved with Radau and an analytic
  Jacobian, from quasi-static initial data at z = 10 and z = 20. The two starts agree to 1e-6 in w0.
- Cross-check: the first-order formula of A.1 agrees with the full solution to 4e-6 in y_H(0) at |f_R0| = 1e-3, and to
  1% (1e-3) and 0.1% (1e-4) in 1 + w0. It is used below 1e-4.
- Units: m² = κ²ρ_m0/3 = 1.
- **Calibrated (single number):** c1/c2 = 6Ω_Λ/Ω_m = 13.05, i.e. the effective Λ, calibrated to Ω_m = 0.315 (Planck
  2018).
- **Free, bounded, not calibrated:** f_R0.

**Background and CPL projection.** Least-squares fit of w0 + wa(1−a) over z ≤ 2:

| \|f_R0\| | w0 | wa | w(z=0) | w(z=1) | Ω_m (true) |
|---|---|---|---|---|---|
| 0.1 | −0.776 | −0.295 | −0.828 | −0.906 | 0.337 |
| 0.03 | −0.898 | −0.136 | −0.937 | −0.955 | 0.327 |
| 0.01 | −0.9655 | −0.0413 | −0.9834 | −0.9794 | 0.319 |
| 1e-3 | −0.99667 | −0.00360 | −0.99846 | −0.99750 | 0.3155 |
| 1e-4 | −0.999668 | −0.000354 | −0.999847 | −0.999745 | 0.31505 |

- Small-field scaling: (1 + w0, wa) ≈ (3.32, −3.53)·|f_R0|.
- w_eff crosses −1 at z ≈ 1.6 (1.59 for |f_R0| = 1e-2, 1.60 for 1e-4) and is phantom above it (at |f_R0| = 1e-2:
  w(2) = −1.010, w(3) = −1.017), matching HS's statement that the models are phantom at high redshift and cross the
  divide [Q][C].
- **Future attractor:** de Sitter, with R → R_1 to 1e-8 and 1 + f_R(R_1) = 0.987 (1e-2) and 0.99988 (1e-4) > 0 [C].

**At the published f_R0 bounds** (exposed, retrospective; first-order formula):

| bound on \|f_R0\| (source, CL, read level) | value | 1 + w0 | wa | Compton λ today [Mpc/h] |
|---|---|---|---|---|
| Schmidt, Vikhlinin, Hu 2009, clusters (0908.2457, 95% stat; TeX section) | 1.3e-4 | 4.3e-4 | −4.6e-4 | 27.6 |
| Cataneo+ 2015, clusters + CMB (1412.0133, 95.4%; TeX section; the arXiv abstract prints "< 4.79", the text −4.79) | 1.6e-5 | 5.4e-5 | −5.7e-5 | 9.8 |
| DESI DR1 full shape, Gonzalez+ 2026 (2604.26915, 95%; abstract + TeX) | 2.6e-5 | 8.5e-5 | −9.1e-5 | 12.3 (paper: λ < 17.81 Mpc) |
| BOSS DR12 full shape, Aviles 2024 (2409.15640, 95%; abstract) | 1.53e-5 | 5.1e-5 | −5.4e-5 | 9.5 |
| Milky Way static thin shell, Hu-Sawicki 2007 eq (fR0bound), v = 300 km/s, "suggestive but not definitive" [Q] | 2e-6 | 6.6e-6 | −7.1e-6 | 3.4 |
| Jain, Vikram, Sakstein 2013, distance indicators (1204.6044, 95%; abstract) | 5e-7 | 1.7e-6 | −1.8e-6 | 1.7 |
| Desmond, Ferreira 2020, galaxy morphology (2009.08743, 1σ; 2σ: 2.8e-8; TeX) | 1.4e-8 | 4.6e-8 | −5.0e-8 | 0.29 |

Secondary quotes, via the DESI paper only [Q-secondary]: SPT clusters + Planck |f_R0| < 4.8e-6; Landim 2024
|f_R0| < 1e-8 (95%).

**Linear growth** (quasi-static G_eff(k,a) = [1 + 4u]/[(1 + 3u)(1 + f_R)], u = (k/a)² f_RR/(1+f_R); Pogosian-Silvestri
0709.0296, abstract [R]). Ratio of the fσ8 proxy (fD) to ΛCDM at z = 0.5:

| \|f_R0\| | k = 0.02 h/Mpc | 0.1 | 0.2 |
|---|---|---|---|
| 1e-4 | +3.3% | +21% | +29% |
| 1e-5 | +0.38% | +6.8% | +15% |
| 1e-6 | +0.04% | +0.93% | +3.3% |
| 1e-7 | +4e-5 | +0.10% | +0.38% |
| 1.4e-8 | +5e-6 | +0.014% | +0.055% |

**Against DESI DR2 w0wa** (ledger row 8, retrospective). The diagonal-Gaussian distance in (w0, wa) ignores the
reported w0–wa correlation and is indicative only. It is 4.1/4.2/5.3/5.9 for the ΛCDM point against DESI+CMB,
+Pantheon+, +Union3 and +DESY5, and the model at |f_R0| ≤ 1e-4 sits at the same numbers to 3 digits. DESI's own
statement of the ΛCDM disfavour is 2.8–4.2σ (Δχ²_MAP) [Q]. Even at |f_R0| = 0.1 the model gives wa = −0.30 and is
2–3 "σ" (diagonal) from the DESI central values. That is the right quadrant, but it is excluded by every growth and
local bound by ≥ 4 orders.

**Verdict on this QV2 (retrospective).** A failed comparison would be reported as a failure of the version and not
repaired inside it.
1. Allowed region: by the exposed growth and local data the screened QV2 survives only for |f_R0| ≲ 1.6e-5 (clusters,
   full shape), and for ≲ 3e-8 if the galaxy-morphology bound is accepted.
2. The DESI tension: in the allowed region its expansion history is the ΛCDM point (|1+w0| ≤ 8.5e-5, or ≤ 5e-8). It
   inherits ΛCDM's retrospective DESI DR2 tension unchanged and cannot account for the DR2 (w0, wa) preference.
3. Internal status: it violates P_q's final-state reading (B.3) and needs q0K'/K ≳ 1e54 to be threshold-stable (B.5).
4. Equivalence: it is viable only as a ΛCDM-mimic, observationally equivalent to QV1-cal with one calibrated number
   (Λ_eff) plus one parameter (f_R0) bounded to be invisible.

**The single calibrated number:** Λ_eff (c1/c2 = 13.05 from Ω_m = 0.315; equivalently R0 in Starobinsky, b/2a in AB).

**Held-out tests.**
- **(a)** Any ≥ 5σ detection of w ≠ −1 by data not yet released when this was written would falsify QV2 in its allowed
  region, since it predicts |1+w0|, |wa| ≲ 1e-4: for example the final DESI BAO release, Euclid DR1 BAO/weak lensing, or
  Rubin/LSST supernovae. Whether a given release is already public has to be checked at use.
- **(b)** Scale-dependent growth at the 1% level at k ~ 0.1 h/Mpc tests |f_R0| ~ 1e-6. Below the galaxy bound
  (1.4e-8) the signal is ≤ 0.06% and beyond Stage-IV reach (Desmond-Ferreira reach the same conclusion [Q]).
- **(c)** Screened fifth-force searches in the laboratory or in void dwarfs.

None of these tests the q-content (μ, P_q). They test the f(R) shape.

---------------------------------------------------------------------------------------------------------------------

## C. Inputs and provenance

| input | value | provenance | role |
|---|---|---|---|
| action, (q), (g) equations | BRIEF.md | KV, credited; verified in equations §1, literature §2 | all |
| m_q² = (2/3)K eps''/K'², α = 1/3, γ(r) | — | equations §4 [D] | A |
| H0 | 67.4 km/s/Mpc = 1.438e-33 eV | Planck 2018 [Q via ledger]; ħ [R] | A.4–A.6 |
| Ω_m, Ω_Λ | 0.315, 0.685 | Planck 2018, ledger row 9 [Q] | A, B (HS calibration) |
| Ω_r | 9.1e-5 | [R] | A numerics |
| ρ_Λ, M | 2.52e-47 GeV⁴, 2.435e18 GeV | BRIEF.md | A.5, B.5 |
| λ(α = 1/3), Eöt-Wash 2020 | 54.0 ± 0.8 μm | digitized from arXiv:2002.11761 Fig. 5 [C from Q-figure] | A.3–A.5 |
| α-limits, 54 μm … 3e14 m | \|α\| < 0.1 | hep-ph/0307284 figures [Q-figure] | A.4 |
| Cassini γ − 1 | (2.1 ± 2.3)e-5 | ledger row 7 [Q] | A.4 |
| G/G_N at recombination | 1.011 ± 0.018 | ledger row 5 [Q] | A.6 |
| Ω_ch² | 0.120 | ledger row 13 [Q] | A.7 |
| α_s(1.6 TeV), R_sun, T0, g* | 0.085, 6.957e8 m, 2.7255 K, 106.75 | [R] | A.4, A.7 [E] |
| Δ_EW | 1.19e8 GeV⁴ | equations §7 [C] | B.5 |
| HS c1/c2 | 13.05 | **calibrated** to Ω_m = 0.315 | B.6 |
| f_R0 bounds | table B.6 | sources as listed [Q] | B.6 |
| DESI DR2 (w0, wa) | ledger row 8 | [Q] | B.6 |
| (eps, K) families F1–F3, toy κ fluid | — | supplied [P], illustrative | A.3, A.6, A.7 |
| δμ in the runs | e.g. −2.0547 | **calibrated** to Ω_Λ = 0.685 (static root) | A.3 |

No parameter is tuned to ρ_Λ,obs and then called a prediction. Λ_cal and c1/c2 are calibrations and are named as such.

## D. Updates to stage-1 records (not edits of them)

1. **ledger row 6 / ledger §4.2:** "for α = 1/3 about 4 meV [E]" becomes **3.65 meV (3.60–3.71), 95% CL**, from the
   digitized curve.
2. **equations §6.4:** the matter-induced gravitating ρ_V = ρ_m²/(12Km²) is correct as a ρ_V. The effective dark energy
   in the Friedmann equation with G_lab is 9 times larger and stiff (A.2). The "m_q ≈ 0.2 H0" figure becomes 0.57 H0
   (formal; outside the regime; w = +1).
3. **ledger §4.2** said "screening would need eps and K designed as in chameleon f(R) ... a design KV do not make". KV
   do discuss and reject it (0807.3896 §5–6, B.3). The rest of that statement stands.

## E. Limits

**Classical analysis only.** No loop corrections to f, no quantum stability of the scalaron, no membrane nucleation.

**A.2.** The relation is first order in H²/m_q². The second-order coefficient c is family dependent (15–53 in the
tested families). For physical m_q it is irrelevant.

**A.4.**
- The Eöt-Wash α = 1/3 value is a digitization (1.4% uncertainty), not a number printed in the paper. The Supplemental
  Material tables were not obtained.
- The long-range exclusion uses 2003-era figures, which are conservative.
- The Cassini bound uses e^{−mb} at the impact parameter (factor < 2 in m).

**A.6.** The φMDE runs start q at its vacuum value at z = 1e4 with v = 0. Other starts change the numbers but not the
O(1)-per-e-fold drive (an attractor of the ω_BD = 0 dust system). The threshold m ≈ H(z_eq) is an estimate.

**A.7.** The relic coefficient comes from a constant-κ toy. The SM numbers are leading-order estimates. Decay widths
are single-channel recollections [R/E].

**B.**
- The HS analysis neglects radiation (z ≤ 10).
- CPL projections are fits over z ≤ 2, not distance-matched.
- The DESI comparison is diagonal-Gaussian and indicative.
- Growth is linear and quasi-static on the ΛCDM background. The O(f_R0) background shift is neglected (≤ 1e-4).
- Nonlinear chameleon screening, which is what the astrophysical bounds rely on, is not computed here. Those bounds are
  quoted.
- The Starobinsky and AB backgrounds were not integrated. Only their vacua and stability were computed, plus their
  published statements.

**Read levels** are listed in §F. The Hu-Sawicki and Starobinsky "[Q]" statements come from the TeX sections named
there.

## F. Sources and read levels (this track; downloads in the scratchpad `.../qtheory/bound_dl/`, none in the repository)

| source | identifier | read level | used for |
|---|---|---|---|
| Starobinsky, "Disappearing cosmological constant in f(R) gravity" | JETP Lett. 86, 157 (2007), arXiv:0706.2041 | full (TeX main text) | B.2, B.3, A.7 |
| Hu, Sawicki | PRD 76, 064004 (2007), arXiv:0705.1158 | sections (TeX §II.A–E, §III.D–E) | B.3, B.6 |
| Appleby, Battye | PLB 654, 7 (2007), arXiv:0705.3199 | full (TeX main text) | B.3, B.5 |
| Klinkhamer, Volovik, "f(R) cosmology from q-theory" | JETP Lett. 88, 289 (2008), arXiv:0807.3896 | full (text, literature-track download) | B.1, B.3 |
| Klinkhamer, Volovik | 0806.2805 (end of §II); 0907.4887 §2 | section (TeX) | B.1, B.3 |
| Lee, Adelberger, Cook, Fleischer, Heckel | PRL 124, 101101 (2020), arXiv:2002.11761 | abstract + constraint paragraph + Fig. 5 (vector data) | A.3, A.4 |
| Kapner et al. | PRL 98, 021101 (2007), hep-ph/0611184 | abstract | A.4 calibration |
| Adelberger, Heckel, Nelson | Ann. Rev. Nucl. Part. Sci. 53, 77 (2003), hep-ph/0307284 | figures midlam/longlam and captions | A.4 |
| Bertotti, Iess, Tortora; Lamine+; Planck 2018; DESI DR2 | via ledger rows 5, 7, 9, 8, 13 | as in ledger | A, B.6 |
| Schmidt, Vikhlinin, Hu | PRD 80, 083505 (2009), arXiv:0908.2457 | section (TeX: the 95% statement) | B.6 |
| Cataneo et al. | PRD 92, 044009 (2015), arXiv:1412.0133 | section (TeX: result and table) | B.6 |
| Gonzalez et al. (DESI) | arXiv:2604.26915 (2026) | abstract + TeX discussion | B.6 |
| Aviles | arXiv:2409.15640 | abstract | B.6 |
| Jain, Vikram, Sakstein | ApJ 779, 39 (2013), arXiv:1204.6044 | abstract | B.6 |
| Desmond, Ferreira | PRD 102, 104060 (2020), arXiv:2009.08743 | sections (TeX) | B.6 |
| Lombriser review | Ann. Phys. 526, 259 (2014), arXiv:1403.4268 | abstract | context |
| Basilakos, Nesseris, Perivolaropoulos | PRD 87, 123529 (2013), arXiv:1302.6051 | abstract | A.1 credit |
| Pogosian, Silvestri | PRD 77, 023503 (2008), arXiv:0709.0296 | abstract; G_eff form [R] | B.6 |
| Frolov | PRL 101, 061103 (2008), arXiv:0803.2500 | abstract | context (f(R) singularity) |
| Amendola, Polarski, Tsujikawa | PRL 98, 131302 (2007), astro-ph/0603703 | recalled; cited by Starobinsky 2007 | A.6 |
| O'Hanlon; Teyssandier-Tourrenc; Sotiriou-Faraoni | 1972; 1983; RMP 82, 451 (2010) | recalled (via equations §2.2) | B.1 |

## G. Files and reproduction

- `code/adiabatic_relation.py` (sympy, < 1 min) → `results/adiabatic_relation.json`: A.1–A.2 (S1–S3) and B.1–B.2, B.5
  (S4–S5).
- `code/frw_regimes.py` (numpy/scipy, ~2 min) → `results/frw_regimes.json`: A.3 (part I), A.6 (part II), A.7 (part III).
- `code/hu_sawicki.py` (~1.5 min) → `results/hu_sawicki.json`: B.3 vacua, B.6 background, CPL, f_R0 bounds, DESI distance,
  future attractor, growth, B.5 threshold check (B4).
- `code/bound_numbers.py` (< 1 s) → `results/bound_numbers.json`: the Eöt-Wash digitization, Cassini, Ω_dyn, the
  regime (ii) gaps, the regime (iii) estimate, the B.5 physical numbers, and AB's G_Mink.

Reproduce:

    cd research/reality_framework/grut_q_vacuum_construction_20261021/tracks/bound/code
    python3 adiabatic_relation.py && python3 frw_regimes.py && python3 hu_sawicki.py && python3 bound_numbers.py
