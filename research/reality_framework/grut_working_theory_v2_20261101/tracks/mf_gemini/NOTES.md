# Track mf_gemini: audit of the external two-field waterfall proposal

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md; GitHub issue #13 context).

This is same-session AI work on one supplied model, not outside replication. It audits a collaborator's design (Gemini,
relayed by the owner) exactly as summarised in BRIEF_MULTIFIELD.md, then its minimal variants. Nothing here is a
validation, an exclusion of the real world, a new theory or a GRUT prediction. **Every potential shape, coupling, slice
convention, threshold and ledger used below is a model choice, not a consequence of closure.** Data comparisons (Planck
2018, BICEP/Keck 2021, ACT DR6) are retrospective. Git was read-only and nothing was committed. All runs were in the
foreground, one process; scratch is `$SP/wt2/mf_gemini/`.

Labels: [D] derived, [C] computed with this track's code (`code/`, output in `results/`), [E] estimate, [Q] quoted from a
reviewed track, [A] argued, [R] recalled or literature (arXiv numbers verified against the arXiv API on 2 Nov 2026).

Units: reduced M_P = 2.435×10¹⁸ GeV. **H0 ≡ (V0/3M_P²)^{1/2}** is the Hubble rate on the valley; it is the proposal's
H_i (its λ_χv⁴/4 = 3M_P²H_i² and m_χ = 3H_i). M1r's H_i is the Hubble rate at the end-of-inflation slice; on the
proposal's own background **H_end = 0.7126 H0** [C], so m_χ = 3H0 = 4.21 H_end. Representative H0 = 1.3×10¹¹ GeV; the
ends of the range are run where it matters.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

| Item | Verdict | Deciding number | Label |
|---|---|---|---|
| G1 amplitude Φ = √6 HM_P/m | **verified** | 1.22 / 0.82 / 0.61 M_P at m = 2 / 3 / 4 H. After the mild waterfall the first swing about v is 0.24 M_P with H there ≈ 0.32 H0, consistent with the identity at the instantaneous H | D/C |
| G1 blocking and loops for a quiet dominant field | **verified** | Γ = 0.5H at m = 3H needs y = 2.89 (Majorana); yΦ/(m/2) = 3.0×10⁷; y⁴Φ⁴/32π² = 1.2×10¹³ × the tree energy | C |
| G2 v_χ, λ_χ | **verified** | v = √(8/3) M_P = 1.633 M_P; λ = 4.45 / 4.81 / 6.84×10⁻¹⁵ at H0 = 1.25 / 1.3 / 1.55×10¹¹ GeV | D |
| G2 CW ratio at y = 2 | **verified** | y⁴/16π²λ = 2.1×10¹³; after the best polynomial counterterms the residual is 1.1×10¹³ V0 (Z2 form) or 1.7×10¹² V0 (bare mass tuned). The 1% cap is y ≤ 3.5–5.5×10⁻⁴ | C |
| G2 blocking ratio | **verified** | M_N = yv = 4.1×10⁷ × m_χ/2 (y = 2); 3.1×10⁷ (y = 1.5). Even unblocked, y = 1.5–2 give Γ/H = 0.13–0.24 (Majorana, one flavour) | C |
| G2 mild waterfall "about 17 e-folds" | **verified as a floor** | N_wf = 16.9–17.5 for a fast σ crossing (Δ ≥ 0.3), rising to 21 (Δ = 0.1), 28 (0.03) and 41 (0.01). Δ = fraction of σ_c crossed per e-fold | C |
| G2 curvature perturbation / PBH | **verified, sharpened** | P_ζ,peak = 0.17–0.70 (linear δN); a nonlinear exponential δN tail gives P(ζ > 0.67) = 0.20–0.27; PBH 2×10¹⁶–10³⁶ g, f_PBH ≈ 10⁶–3×10¹⁶ | C |
| G2 domain walls (real χ) | **verified** | v = 1.633 M_P is 0.99 of the topological-inflation threshold 0.33 m_Pl (gr-qc/9506068); wall tension ≈ 5.3 H0 M_P², which dominates on horizon entry | D/Q |
| G2 "holder on χ is off during inflation" | **refuted as stated** | the holder −(g e^{iα}(χ − v)Σ*Y² + h.c.) has u = μ_late − Gv during inflation (on, stationary), g = 4.0×10⁻⁷. But ±v domains get different late μ and opposite CP phase | D |
| Any choice meeting all the proposal's requirements | **none found** | the m_χ = 3H0 plus "χ carries V0" pair forces a Planckian quiet field. No radiatively acceptable channel, perturbative or not, gives Γ_eff ≥ 0.5 H; real χ fails the no-defect/PBH requirement for every g_hyb and V_SR | C/A |
| Claimed n_s 0.967–0.973 | **conditionally verified** | convex V_SR gives n_s > 1 (1.004–1.040 for η = 0.002–0.02). Hilltop: |η| = 0.0136–0.0166, α_s ≈ −4×10⁻⁹. Log (CW-like): needs N_σ = 30–37, i.e. N_wf ≤ 23, α_s = −(0.7–1.1)×10⁻³. r = 2.5–3.9×10⁻⁷ | C |
| Minimal change | **a third field and a heavy waterfall** | heavy χ (λ ≈ 0.013, m ≈ 3.9×10³ H, v ≈ 1.3×10⁻³ M_P, y ≤ 0.08, loop 2×10⁻⁵ λ: G3 verified) plus a separate quiet subdominant holder field (G5) plus a defect fix | C/A |

**Bottom line.** The proposal's two requirements m_χ = 3H_i and λ_χv⁴/4 = V0 put the holder on a quiet field that also
carries the vacuum energy, so the field has a Planckian vev and amplitude (G1). That is the single-field obstruction of
the inflaton pass, moved to the waterfall field.
- **Fast reheating.** None of the routes tried, perturbative (y_N, λ_χH) or non-perturbative (fermionic, bosonic or
  tachyonic preheating), gives Γ_eff ≥ 0.5 H_i without spoiling the potential or destabilising the Higgs.
- **Defects and PBHs.** With a real χ the waterfall is mild (≥ 17 e-folds, quantum-seeded). It produces Planck-scale walls
  and an O(1) non-Gaussian curvature perturbation on scales inside today's horizon, and PBHs overclose by ≥ 10⁶.
- **What survives:**
  - the low-scale V0-dominated structure, with r ≈ 3×10⁻⁷;
  - the n_s window, but only with a non-convex V_SR;
  - a holder on (χ − v), as an algebraic possibility;
  - a Higgs kept heavy in inflation, but only with the opposite portal sign to the one written, or with a σ portal.

---------------------------------------------------------------------------------------------------------------------

## 1. The proposal as tested

V(σ, χ) = λ/4 (χ² − v²)² + ½ g² σ² χ² + V_SR(σ), with:
- λv⁴/4 = 3M_P²H0² and m_χ = √(2λ) v = 3H0;
- χ held at 0 for σ > σ_c, with g²σ_c² = λv²;
- y_N χ N̄N and λ_χH χ²|H|² for reheating;
- V_SR chosen for n_s.

Code units in `code/mfg_common.py`: M_P = H0 = 1, so v² = 8/3 and λ = 27/16.
- The quantum noise H/2π enters through hub = H0/M_P = 5.34×10⁻⁸.
- A_s = 2.10×10⁻⁹ fixes the pivot slope √(2ε_*) = 1.85×10⁻⁴.

**V_SR near the waterfall is not specified by the proposal.** I parametrise it by the slope at σ_c, p = √(2ε_c), and
Δ = p/σ_c, the fraction of σ_c crossed per e-fold. Two local shapes, both model choices:
- `lin`: a constant slope through the band |σ| < σ_c. σ can leave the band on the far side.
- `quad`: convex near σ_c, with σ decelerating towards 0.

The default slope is p_c = 1.85×10⁻⁴ e^{0.015×36} = 3.18×10⁻⁴, a hilltop η = −0.015 rolled over 36 e-folds. The pivot
slope p_* is run as a variant. Since g = √λ v/σ_c, a scan in Δ at fixed p is a scan in g_hyb: g = 3.6×10⁻⁶ to 1.1×10⁻³.

## 2. G1 [D/C] (`code/g1g2_checks.py` → `results/g1g2.json`)

- **The identity [D].** ½m²Φ² = 3M_P²H² gives Φ = √6 HM_P/m: 1.225 / 0.816 / 0.612 M_P at m = 2 / 3 / 4 H.
  - It holds with the instantaneous H.
  - On the proposal's own background the first overshoot about v is 0.147 v = 0.24 M_P. The Hubble rate at that turning
    point is ≈ 0.32 H0, which is what the identity predicts [C, `waterfall_bg.py`, `vscan.py`].
- **Fast perturbative decay [C].**
  - Γ = 0.5 H at m = 3H needs y = 2.894 (Majorana, Γ = y²mβ³/16π, inf_reheat R4 [Q]) or 2.05 (Dirac).
  - yΦ/(m/2) = 3.0×10⁷, so the decay is kinematically blocked.
- **One loop [C].** y⁴Φ⁴/32π² is 1.2×10¹³ times the tree energy.
- **Conclusion [A].** A quiet field that carries V0 cannot decay fast through any perturbative, radiatively acceptable
  coupling. The holder must sit on a subdominant field, or the dominant field must not be quiet. §4 extends this to the
  non-perturbative routes.

## 3. G2: the proposal as written

### 3.1 Parameters, loops, blocking [D/C]

| | H0 = 1.25×10¹¹ | 1.30×10¹¹ | 1.55×10¹¹ GeV |
|---|---|---|---|
| v/M_P; λ | 1.633; 4.45×10⁻¹⁵ | 1.633; 4.81×10⁻¹⁵ | 1.633; 6.84×10⁻¹⁵ |
| V0^{1/4} (GeV) | 7.26×10¹⁴ | 7.40×10¹⁴ | 8.09×10¹⁴ |
| y⁴/16π²λ at y = 2 (1.5) | 2.3×10¹³ (7.2×10¹²) | 2.1×10¹³ (6.7×10¹²) | 1.5×10¹³ (4.7×10¹²) |
| CW residual / V0 after counterterms, y = 2: Z2 form m_N = yχ; bare mass tuned m_N = y(χ − v) | 1.2×10¹³; 1.9×10¹² | 1.1×10¹³; 1.7×10¹² | 7.6×10¹²; 1.2×10¹² |
| y cap (residual ≤ 1% V0): Z2; tuned | 3.4×10⁻⁴; 5.4×10⁻⁴ | 3.5×10⁻⁴; 5.5×10⁻⁴ | 3.8×10⁻⁴; 6.0×10⁻⁴ |
| Γ/H at the cap (unblocked, Majorana) | 7.0×10⁻⁹ | 7.3×10⁻⁹ | 8.7×10⁻⁹ |
| M_N = yv over m_χ/2, y = 2 | 4.2×10⁷ | 4.1×10⁷ | 3.4×10⁷ |

- **The counterterm fit.** It is a least-squares fit of {1, χ², χ⁴}, or of {1, χ, …, χ⁴} when the bare mass is tuned, over
  χ ∈ [0, 1.5v]. The χ⁴ ln χ² shape cannot be absorbed.
- **The proposal's y_N = 1.5–2 would also be too slow if unblocked.** It gives Γ/H = 0.13–0.24 (Majorana) or 0.27–0.48
  (Dirac) per flavour, so Γ_eff/H = 0.5–1 needs N_f ≈ 2–4 even before blocking.

### 3.2 The mild waterfall: duration and background (`code/waterfall_dn.py`, `waterfall_bg.py`, `vscan.py`) [C]

**Method.**
- A σ-only background (χ = 0 exactly). At linear order the χ modes then decouple from the metric.
- Bunch–Davies χ modes, WKB-started at k/aH = 100, or where χ is heavy enough to be adiabatic.
  - Regression: in exact de Sitter the massless and conformal (m² = 2H²) spectra are reproduced to 5×10⁻⁵ and 4×10⁻⁵.
- Coarse-grained χ_rms from the super-horizon modes.
- Matching at χ_rms = 10⁻² χ_trap, with χ_trap = (2pσ_c/3)^{1/2}, the value at which g²χ²σ balances the slope. There the
  back-reaction on σ is ≤ 10⁻⁴.
- After matching, the full two-field equations, then a WKB-averaged σ (adiabatic invariant) once gχ > 30H, until ε_H = 1.

The tachyonic rate is s = 1.098 per e-fold at most, i.e. m² = −4.5H² at σ = 0.

**One bug found and fixed during the run.** The WKB amplitude was first coded as 1/(a√2ω); the correct form is
1/(a^{3/2}√2ω).
- The normalisation-free quantities were unaffected: ΣP_ζ dlnk = 1/s² and the δN tail.
- N_wf and the spectral shape were affected, and all reported numbers come from the corrected runs.

| case | Δ | σ_c (M_P) | g | N_wf from σ_c (classical / mean) | P_ζ peak (linear δN) | peak exits (e-folds before end) | M_PBH (g) | P(ζ > 0.67), nonlinear δN | f_PBH (tail) |
|---|---|---|---|---|---|---|---|---|---|
| quad | 0.01 | 0.032 | 3.6×10⁻⁶ | 39.8 / 40.7 | 0.25 | 38.6 | 1.0×10³⁶ | 0.274 | 5.5×10⁶ |
| quad | 0.03 | 0.011 | 1.1×10⁻⁵ | 27.4 / 28.1 | 0.22 | 26.3 | 1.9×10²⁵ | 0.233 | 1.1×10¹² |
| quad | 0.1 | 3.2×10⁻³ | 3.6×10⁻⁵ | 20.6 / 21.1 | 0.25 | 19.5 | 2.5×10¹⁹ | 0.206 | 8.5×10¹⁴ |
| quad | 0.3 | 1.1×10⁻³ | 1.1×10⁻⁴ | 17.9 / 18.5 | 0.36 | 16.9 | 1.5×10¹⁷ | 0.201 | 1.1×10¹⁶ |
| quad | 1 | 3.2×10⁻⁴ | 3.6×10⁻⁴ | 16.9 / 17.5 | 0.56 | 15.9 | 2.0×10¹⁶ | 0.200 | 2.9×10¹⁶ |
| quad | 3 | 1.1×10⁻⁴ | 1.1×10⁻³ | 17.0 / 17.6 | 0.70 | 16.0 | 2.4×10¹⁶ | 0.200 | 2.7×10¹⁶ |
| lin | 0.03 | 0.011 | 1.1×10⁻⁵ | 25.8 / 26.5 | 0.17 | 24.8 | 1.1×10²⁴ | 0.216 | 4.3×10¹² |
| lin | 0.06 | 5.3×10⁻³ | 2.1×10⁻⁵ | 21.4 / 22.0 | 0.19 | 20.4 | 1.5×10²⁰ | 0.200 | 3.4×10¹⁴ |
| lin | 0.1 | 3.2×10⁻³ | 3.6×10⁻⁵ | 19.7 / 20.5 | 0.32 | 18.7 | 5.0×10¹⁸ | 0.214 | 2.1×10¹⁵ |
| lin | 0.2 | 1.6×10⁻³ | 7.1×10⁻⁵ | **no waterfall in this crossing**: σ leaves the band first | — | — | — | — | — |

Notes on the table:
- The p_* variants (Δ = 0.1, 1) and H0 = 1.25 / 1.55×10¹¹ (Δ = 0.3) agree with these rows to < 3% in N_wf. The PBH mass
  scales with H.
- `results/requirements.json` holds Monte Carlo Var(δN) from each N(c) curve. Its agreement with π²/(8s²) is 1–5% for
  Δ ≥ 0.03.

**Readings [C].**
- **Duration.**
  - The waterfall lasts at least ≈ 17 e-folds. That is set by ln(v/(H/2π))/s with s ≤ 1.1, which no choice of g_hyb or
    V_SR can shorten for a quantum-seeded real χ.
  - A slow crossing (Δ ≲ 0.1) lengthens it.
  - A constant-slope V_SR with Δ ≳ 0.15 does not complete the waterfall at all: χ grows by only ≈ e^{1.5/Δ} per band
    crossing.
- **End of inflation.**
  - It is universal: H_end = 0.7126 H0 and χ = 0.647 v.
  - σ carries ≤ 6×10⁻⁷ of the energy (m_σ = gv = 99–2964 H_end).
  - **G4's "σ, not χ, ends up holding V0" therefore does not apply to the proposal as written.** χ holds the energy, which
    is exactly the G1 obstruction.
- **After the end** (`waterfall_bg.json`).
  - χ oscillates about v with overshoot 0.147 v, then return 0.112 v, and decays. It never comes back near 0.
  - The oscillation-averaged w = −0.03 → −10⁻⁴, i.e. matter.
  - A sudden decay at t = 0.5 / 1 / 2 / 4 / 8 / 16 H_end⁻¹ would give K = a²H/H_end = 0.98 / 1.21 / 1.31 / 1.66 / 1.97 /
    2.48.
- **Slice convention.**
  - m_χ = 4.21 H_end on the slice M1r uses. inf_reheat found m_I = 4 already "marginal" (N_YL = 2.3) [Q].
  - m_χ = 3 H_end would need v ≈ 3.2 M_P and ≈ 55 quantum-seeded e-folds (`vscan.json`), which would put the CMB scales
    inside the waterfall.

### 3.3 Curvature perturbation and PBHs [C/E]

- **Linear δN** (Clesse-type, on the classical trajectory from χ_rms).
  - P_ζ(k) = (∂N/∂χ)² P_χ(k) peaks for modes exiting ≈ 1 e-fold after σ_c, at 0.17–0.70.
  - ΣP_ζ dlnk = 1/s² = 0.83 (Δ ≥ 0.3), independent of normalisation [D/C].
- **Nonlinear δN.**
  - χ_L is Gaussian with zero mean, and N(χ_L) ≈ N_0 − ln|χ_L/χ_rms|/s. So Var(δN) = π²/(8s²) ≈ 1.0 [D], confirmed by Monte
    Carlo to 1–5%.
  - The tail is exponential, P(ζ > x) ≈ √(2/π) e^{−sx}: patches that start near χ = 0 inflate longer.
  - These patches are the would-be domain-wall cores (§3.4).
- **Stochastic cross-check** (`waterfall_stoch.py`; phase-space Langevin, first passage, 2×10⁴ realisations).
  - P(δN > 0.67) = 0.199 / 0.206 / 0.229 at Δ = 1 / 0.3 / 0.1.
  - Var = 1.00 / 1.06 / 1.44.
  - The tail slope is 1.14 / 1.10 / 0.95, against s = 1.098 / 1.098 / 1.06.
  - Halving dN changes the tail by 1%.
  - The noise amplitude is uncertain by O(1) for a non-light χ; the tail shape is not.
- **PBHs [E].**
  - The peak scale exits 16–39 e-folds before the end. In instant radiation the horizon mass gives M_PBH ≈ 2×10¹⁶ g
    (fast crossing) up to 10³⁶ g (Δ = 0.01).
  - With the tail β ≈ 0.2–0.27, f_PBH = Ω_PBH/Ω_DM ≈ 10⁶–3×10¹⁶ by radiation-era scaling.
  - That contradicts the measured dark-matter density by many orders, whatever the detailed collapse criterion. Bounds:
    Carr et al. arXiv:2002.12778 [R].
  - The same conclusion is reached by Kawasaki–Tada, arXiv:1512.03515 [R, abstract]: PBHs are overproduced when the mild
    waterfall lasts long enough for the PBH scale to be large.
  - Clesse–García-Bellido, arXiv:1501.07565 [R], tuned a mild waterfall for PBH dark matter. Their parameters are not the
    proposal's, whose v = 1.63 M_P and m_χ = 3H0 fix s = 1.098 and leave no small parameter in the tail.

### 3.4 Domain walls (real χ) [D/Q]

- **Wall formation.** The waterfall breaks χ → −χ on scales that exit near σ_c, 16–40 e-folds before the end: well inside
  today's horizon.
- **Tension.** σ_w = (2√2/3)√λ v³ = (2/3)m_χv² = 5.3 H0 M_P². One wall per horizon reaches ρ_w ≈ σ_wH ≈ 3M_P²H² at
  H ≈ 1.8 H0, so the walls dominate as soon as they are inside the horizon [D].
- **Topological inflation.** v = 3.98×10¹⁸ GeV is 0.99 of the topological-inflation threshold η ≳ 0.33 m_Pl = 4.0×10¹⁸ GeV
  (Sakai–Shinkai–Tachizawa–Maeda gr-qc/9506068; Vilenkin hep-th/9402085; Linde astro-ph/9402031) [R]. The wall cores are
  marginally inflating.
- **One problem, two faces.** The walls and the δN tail of §3.3 are the same physics.
- **The fix is a bias** (§6). A complex χ would give global strings with v ~ M_P.

### 3.5 Holder placement [D/A]

G2's statement "χ = 0 during inflation, so a holder on χ is off" is **not right as stated**.
- **Writing it.** Write the holder as −(g e^{iα}(χ − v)Σ*Y² + h.c.), where the −gv piece is absorbed in μ_bare. Then:
  - u = μ_late + G(χ − v);
  - during inflation u = μ_late − Gv: on, and stationary (χ ≡ 0 classically, δχ ~ H/2π ≪ v);
  - it switches off as χ → v.
- **Its size.** M1r's u = 10i H at μ_late = 7.16 H needs |G|v = 12.3 H, so g = 4.0×10⁻⁷.
- **The hidden-sector tadpole it puts on χ during inflation** (M1r's prepared |Y|, |Σ| as a stand-in) is 2.3×10⁻³ of the
  quantum noise per e-fold [E]. It does not bias the waterfall, so walls still form.
- **The ±v problem [D].**
  - In −v domains the late μ is −7.16 + 20i, i.e. |μ| = 21.2 H instead of 7.16 H, and the CP phase is −0.34 instead of +π/2.
  - So the two vacua carry different hidden sectors and opposite-sign charge.
  - In M1r this does not arise, because the inflaton has a single vacuum φ_I = 0.
- **The release is not M1r's [A].**
  - The holder relaxes over the 17–40 e-fold waterfall; G(χ − v) changes appreciably only in the last ≈ 2 e-folds.
  - It then oscillates with 15% of the full holder amplitude at m = 4.2 H_end.
  - The chain would have to be re-run (mf_ideal / mf_chain machinery).
- **Holder on σ: fails.**
  - After the waterfall m_σ = gv = 99–2964 H_end on the rows above, which is loud.
  - m_σ ≤ 4 H_end needs σ_c ≥ 1.2 M_P, hence Δ ≲ 3×10⁻⁴ and a waterfall far longer than N_*.
- **A spectator holder** (G5) is outside the proposal's field content.

## 4. Does any choice meet the proposal's requirements? (`code/requirements.py` → `results/requirements.json`)

Requirements: m_χ = 3H_i, Γ_eff/H_i ≥ 0.5, H_i ∈ [1.25, 1.55]×10¹¹ GeV, n_s = 0.967–0.973, no defects, a working holder.
With λv⁴/4 = V0 the first one fixes v = 1.633 M_P [D].

| Route to Γ_eff ≥ 0.5 H | Best case | Verdict | Label |
|---|---|---|---|
| perturbative y_N χ N̄N | M_N = yv ≫ m_χ/2 (×4×10⁷). With the bare mass tuned (to 2×10⁻⁸ of yv) the oscillation still blocks it (yΦ/(m/2) = 6×10⁶ at y = 2). The loop caps y ≤ 5.5×10⁻⁴, giving Γ ≤ 7×10⁻⁹ H | fails | C |
| fermion preheating (Greene–Kofman, hep-ph/9807339), bare mass tuned | at y = 2 the first-crossing formula gives 112 × ρ (fully non-perturbative; m_N ≈ 0.48 M_P at the turning points), but y = 2 destroys the potential. At the loop cap: 1.4×10⁻⁷ per crossing, Pauli-limited | fails | C/E |
| bosonic preheating, added scalar ½g_φ²(χ − v)²φ² (needs a Z2-odd trilinear), at the loop cap g_φ² = 4.3×10⁻⁷ (q = 2.4×10⁵) | KLS crossing map (hep-ph/9704452), no back-reaction, instant decay of the products: 10⁻⁷ at the first crossing, ≤ 4×10⁻⁷ by 2 H_end⁻¹. 50% only at 9.7 (maximal growth), 13.4 (random phases) or 23 (WKB phases) H_end⁻¹, i.e. K ≈ 2.1–2.9. At 100× the cap: 6.7–11 H_end⁻¹ | fails Γ_eff ≥ 0.5H; K ≈ 2–3 not excluded (lattice needed) | C on a model / E |
| Higgs portal λ_χH χ²\|H\|² (proposal's sign), EW vacuum tuned | Γ = 0.5H needs λ_χH = 1.0×10⁻⁷ (A = 6.14 H). At χ = 0, m_H² = −λ_χHv² = −9.4×10⁷ H², so the Higgs condenses during inflation and releases ≈ V0 (the tree threshold λ_χH²/λ_χ = 2.09 makes the release exactly V0 for λ_SM → 0): **there is no inflationary valley**. λ_H(full) = 2.1 hits a one-loop Landau pole at ≈ 9×10¹² GeV < V0^{1/4} | fails | D/C |
| Higgs portal with the opposite sign, ∝ (v² − χ²)\|H\|² | heavy during inflation ✓; tachyonic while χ > v, with ∫\|m_H\|dt = 4.3×10³ per window. The full potential is bounded only if λ_SM > 0 [D], which fails above the SM instability scale. Safe \|λ_χH\| ≤ 1–4×10⁻¹⁵ gives Γ ≤ 6×10⁻¹⁶ H | fails | D/C |
| tachyonic preheating of χ itself (Felder et al. hep-ph/0012142, hep-th/0106179; García-Bellido et al. hep-ph/0208228) | the tachyonic rate is ≤ 1.1 H: the "preheating" is the 17-e-fold mild waterfall, with no sub-horizon fragmentation. The energy stays in massive χ, and w → 0 for a quadratic minimum (Lozanov–Amin 1608.01213, 1710.06851) | fails (no SM radiation) | C/Q |
| equation of state after preheating | rapid w → 1/3 (Podolsky et al. hep-ph/0507096) needs efficient preheating into light fields first. Computed w ≈ 0 (−0.03 → 10⁻⁴) | — | C/Q |
| axial coupling (inf_reheat S2) | Γ ≥ 0.5H needs f < H ≪ Φ, outside its EFT [Q inf_reheat skeptic] | fails | Q |

**Other requirements.**
- **H_i range.** Trivial: r = 2.5 / 2.75 / 3.9×10⁻⁷.
- **n_s.** Achievable with a non-convex V_SR (§5).
- **No defects.** Fails for real χ at every g_hyb and V_SR tried (§3.3–3.4). Wrong only in the Clesse regime, N_wf > N_*
  (§6), which is not the proposal.
- **Holder.** Algebraically possible on χ (§3.5).

**Answer [C/A].** No choice of (y_N, λ_χH, g_hyb, V_SR, holder placement) with m_χ = 3H_i and λv⁴/4 = V0 meets all the
requirements. Fast reheating fails on its own, whatever the other choices.
- **The only non-perturbative route not excluded** is an added bosonic channel at the loop cap. It gives K ≈ 2–3 in an
  optimistic linear model, i.e. Γ_eff ≈ ln2/t_½ ≈ 0.03–0.07 H_end, and needs a lattice.
- **A supersymmetric completion** could cancel the N loop [A]. But in F-term form it ties g² = 2λ, so σ_c ≈ v/√2 ~ M_P:
  a different model (Clesse's waterfall-inflation regime). It still needs the bare mass tuned to 2×10⁻⁸. Not pursued.

**G6 (Higgs sign structure) [D].**
- A portal linear in δχ about v gives a sign-oscillating Higgs mass², verified: 2λ_χH v δχ.
- **But G6's "a positive χ²|H|² … keeps the Higgs heavy during inflation" is wrong for this model.** χ = 0 during
  inflation, so χ²|H|² gives no mass there; with the EW vacuum tuned it gives a tachyonic one.
- Only the (v² − χ²)|H|² sign does, and it is tachyonic during the χ > v overshoot.
- The Lebedev–Westphal mechanism (arXiv:1210.6987) [R] needs the field that is large during inflation. Here that is σ:
  λ_σH σ_c² ≥ H² means λ_σH ≥ (H/σ_c)² ≈ 3×10⁻⁹ (σ_c = 10⁻³ M_P) [E].
- The tree threshold λ_χH²/λ_χ is the Elias-Miró–Espinosa–Giudice–Lee / Lebedev effect (arXiv:1203.0237, 1203.0156) [R].
- Higgs fate in preheating: Ema–Mukaida–Nakayama, arXiv:1602.00483 [R]. The brief's "1602.07719" is an unrelated paper.
- With λ_χ ≈ 5×10⁻¹⁵, any portal large enough to reheat dominates the Higgs quartic.

## 5. The n_s claim (`code/ns_hybrid.py` → `results/ns_hybrid.json`) [C]

Inputs:
- N_* = 53.2 (K = 1), 52.8 (K = 2), 52.6 (K = 3) at H0 = 1.3×10¹¹. This is the Planck 2018 X matching form with
  −(1/12) ln g_th.
- N_σ = N_* − N_wf.
- Second-order slow roll; V0 dominance makes ε ≈ 1.7×10⁻⁸.

| V_SR family | n_s | α_s | r (H0 = 1.25 / 1.3 / 1.55×10¹¹) |
|---|---|---|---|
| convex ½m²σ² (η = 0.002 / 0.005 / 0.01 / 0.02) | 1.004 / 1.010 / 1.020 / 1.040: **blue** [D/C] | ≈ 0 | 2.5 / 2.75 / 3.9×10⁻⁷ |
| quadratic hilltop, maximum at σ_m > σ_c | 0.967–0.973 for \|η\| = 0.0136–0.0166, independent of N_σ | −4×10⁻⁹ | same |
| log, V0 A ln σ (CW-like), σ_c ≪ σ_* | 1 − 1/N_σ (second order included): window N_σ = 30.1–36.9; κ = (8π²A)^{1/2} ≈ 0.012–0.013, σ_* ≈ 0.011–0.013 M_P | −1.1×10⁻³ (N_σ = 30) to −7.7×10⁻⁴ (37) | same |
| log, with this waterfall | N_wf = 17.4–40.7 → N_σ = 35.8–12.2 → n_s = 0.972–0.92. **Inside the window only for N_wf ≤ 22.7, i.e. Δ ≳ 0.06–0.1.** At N_wf = 17.4: n_s = 0.972, α_s = −8×10⁻⁴ | | |

Readings:
- **The claim survives only with a non-convex V_SR**, as the brief says [D].
- **The hilltop.** Its direction must be chosen with care. With the proposal's +½g²σ²χ², χ is tachyonic for |σ| < σ_c, so
  the hilltop maximum must sit at σ_m > σ_c, with σ rolling down towards σ_c. G4's example (σ_* = 0.0106 < σ_c = 0.028,
  hilltop at 0) rolls the wrong way for this coupling. It would need an inverted coupling, which is inconsistent with its own
  valley formula [D] (for mf_hybrid).
- **The log term.** It needs an O(10⁻²) coupling to a SUSY-like sector [A]. The χ loop itself, g⁴σ⁴ ln σ²/64π², gives
  ≤ 3% of the required |η| at g ≤ 10⁻³ [E].
- **Retrospective pulls.**
  - n_s = 0.970: +1.2σ against Planck 2018 (0.9649 ± 0.0042, 1807.06211) and −1.3σ against ACT DR6 P-ACT-LB (0.9743 ± 0.0034;
    abstract 0.974 ± 0.003, 2503.14452).
  - The window 0.967–0.973 is at +0.5 to +1.9σ (Planck) and −2.1 to −0.4σ (ACT).
  - ACT DR6 running: dn_s/dlnk = 0.0062 ± 0.0052 (2503.14454). Both α_s values are inside it.
  - r ≈ 3×10⁻⁷ is far below BICEP/Keck's 0.036 (2110.00483) [R].

## 6. Minimal change, and which claims survive

**Variants within the two-field content [C/E/A].**
1. **Bias χ** (a small χ-odd term, a model choice that breaks the reflection).
   - A classical seed χ_b ≥ 5.3 H0 (6.8×10¹¹ GeV) removes the δN tail and leaves < 10⁻² wrong-sign domains in e^{3×36}
     patches.
   - P_ζ at the onset ≈ 8×10⁻⁴ (PBH-safe); the waterfall then lasts ≈ 14 e-folds classically [E].
   - This fixes **defects and PBHs** only. It also selects +v, so the holder on (χ − v) gets a single CP sign set by α and
     the sign of the bias.
2. **Push the waterfall beyond the CMB** (N_wf > N_* ≈ 53; Δ ≲ 3×10⁻³; Clesse arXiv:1006.4522, Kodama–Kohri–Nakayama
   arXiv:1102.5612).
   - Walls and PBHs leave the observable patch, but the CMB then comes from the waterfall, not from V_SR.
   - n_s ≈ 1 − 6Δ(N_wf − N_*) [E] is fine-tuned, and a classical χ seed is needed.
   - Not the proposal; does nothing for reheating.
3. **Nothing in the two-field content fixes reheating while χ is quiet and dominant** (G1, §4).

**The minimal change that works on paper (not built here).** Separate the roles, which is the mf_hybrid architecture:
- (a) **V0 on a heavy waterfall χ** (G3 verified at H0 = 1.3×10¹¹; `g1g2.json` G3_check):
  - for Γ = 0.5 H: λ = 0.013, m_χ = 3.9×10³ H, v = 1.27×10⁻³ M_P, y_max = √(λ/2) = 0.081, y⁴/16π²λ = 2×10⁻⁵;
  - for Γ = H: λ = 0.023, m = 4.4×10³ H, v = 1.1×10⁻³ M_P.
  - The waterfall is then fast (|m_tach| ≈ 2.7×10³ H), with no mild-waterfall tail.
- (b) **The holder on a separate, subdominant, quiet field** (spectator s, G5), with m_s ≈ 3 H_end after inflation.
- (c) **A defect fix.** The cases:
  - a real heavy χ: bias, smooth or shifted hybrid (hep-ph/9506325, hep-ph/0002151);
  - a complex χ: global strings at v ≈ 3×10¹⁵ GeV are a problem.
- (d) **Check σ's share of V0** after a fast waterfall (G4).

This changes m_χ and adds a field. Within the proposal's own field content and requirements, no choice works.

**Claims of the proposal.**

| Claim | Status |
|---|---|
| S/Q < 3.47, Q-matter at Q ~ 10²⁰ | **not supported.** On M1r's own release at e = 10⁻⁶, S/Q = 3.47 would need K ≤ 1.38, i.e. complete reheating within ≈ 2.4 H_end⁻¹ (R_vec ∝ K², dq_entropy closed form) [E]. The proposal's history is not that release (§3.5), and its best channel gives K ≈ 2–3 at most. Q ~ 10²⁰ is conditional on S/Q < 3.47 [Q dq_critic] |
| Higgs rescued by the plasma | **not supported**: no fast plasma (§4). The written portal sign destabilises the Higgs in inflation; the opposite sign gives tachyonic windows |
| leptogenesis kept | **not established**: N with M_N = yv ≈ 3 M_P is not produced perturbatively. With the bare mass tuned, production is non-perturbative and needs y ≫ the loop cap. Not computed |
| V0^{1/4} = 2.4×10¹⁵ GeV | corrected: 7.40×10¹⁴ GeV (reduced M_P) [D] |
| low scale, r ~ 3×10⁻⁷ | **survives**: r = 2.5–3.9×10⁻⁷ [C] |
| n_s 0.967–0.973 from V_SR | **survives with a non-convex V_SR** (§5) |
| fast reheating, Γ_eff/H ≈ 0.5–1 | **fails** for every channel with radiatively acceptable couplings (§4) |
| quiet holder field with m = 3H_i | **survives as a parameter choice** (4.2 H_end at the M1r slice). Its placement on (χ − v) works algebraically, with the ±v caveat |

## 7. Notes for the other tracks

- **mf_hybrid.**
  - G4's hilltop direction (§5).
  - The mild-waterfall floor: any V0-dominated waterfall with m_tach ≲ 2H lasts ≥ 17 e-folds and has an O(1) δN tail, so
    the fast-waterfall condition is |m_tach| ≫ H.
  - In the proposal-as-written regime σ's energy share is ≤ 6×10⁻⁷.
  - The brief's 1602.07719 should read 1602.00483.
- **mf_ideal / mf_chain.**
  - The K table for a sudden decay after the mild waterfall: K(1, 2, 4, 8 H_end⁻¹) = 1.21, 1.31, 1.66, 1.97.
  - K_crit = 1.38 for Q-matter at e = 10⁻⁶ on M1r's release [E].
- **mf_higgs.**
  - With χ = 0 during inflation, only a σ portal or the (v² − χ²) sign keeps the Higgs heavy.
  - Any χ portal strong enough to reheat makes λ_H(full) ≥ 2.

## 8. Limits

- **The model is the brief's summary of the proposal**, not its full text (none is in the repository or on issue #13). V_SR
  near σ_c is unspecified, and I scanned two local shapes.
- **The waterfall.**
  - Linear modes without metric perturbations (exact at linear order for χ̄ = 0). One matching point. A classical
    two-stage background with a WKB-averaged σ (checked a posteriori: m_σ/H = 30 → 3000).
  - The stochastic check uses light-field noise for a non-light χ, so its amplitude is O(1)-uncertain; its tail slope
    matches s.
- **PBHs [E].** ζ threshold 0.5–1, γ = 0.2, instant radiation after the end, no non-Gaussian collapse modelling. With tail
  fractions β ≈ 0.2 the conclusion does not depend on these choices.
- **Preheating.**
  - The bosonic channel is a crossing-map (KLS) model with no back-reaction, rescattering or lattice. It is optimistic
    for the transfer, and assumes instant decay of the products to the SM.
  - Fermions: analytic Greene–Kofman only.
- **The loop caps.** "Residual ≤ 1% of V0 after polynomial counterterms on χ ∈ [0, 1.5v]" is a criterion choice. A looser
  10% moves the caps by 10^{1/4}.
- **The Higgs.**
  - λ_SM(m_χ) is set to 0 for the threshold algebra.
  - The Landau pole is one-loop from 24λ²/16π² only.
  - No SM running was redone.
- **The holder statements are [D/A]**: no release was run on the proposal's history.

## 9. Code and results

| file | content |
|---|---|
| `code/mfg_common.py` | model, V_SR families, background, two-stage (full → WKB-averaged σ) integrator |
| `code/waterfall_dn.py` | linear χ modes, matching, nonlinear N(c), linear and nonlinear δN, PBH estimates, de Sitter regression → `results/waterfall_dn.json` |
| `code/waterfall_stoch.py` | phase-space Langevin first-passage δN → `results/waterfall_stoch.json` |
| `code/waterfall_bg.py` | end of inflation, χ oscillation, w, K(t) → `results/waterfall_bg.json` |
| `code/vscan.py` | m_χ/H_end, durations, overshoot vs v → `results/vscan.json` |
| `code/g1g2_checks.py` | G1, G2, loop caps, Greene–Kofman fractions, Higgs portal, G3 → `results/g1g2.json` |
| `code/preheat_lin.py` | bosonic crossing-map preheating at the loop cap → `results/preheat_lin.json` |
| `code/ns_hybrid.py` | n_s, α_s, r for convex / hilltop / log V_SR with this waterfall → `results/ns_hybrid.json` |
| `code/requirements.py` | Monte Carlo Var(δN), holder algebra, K_crit, biased seed, requirements matrix → `results/requirements.json` |

Reproduce, one core, about 15 minutes in total:

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_gemini/code
python3 -B waterfall_dn.py regress
python3 -B waterfall_dn.py part 0 7; python3 -B waterfall_dn.py part 7 14; python3 -B waterfall_dn.py combine
python3 -B waterfall_stoch.py; python3 -B waterfall_bg.py; python3 -B vscan.py
python3 -B g1g2_checks.py; python3 -B preheat_lin.py; python3 -B ns_hybrid.py; python3 -B requirements.py
```

## Literature (IDs and titles verified against arXiv on 2 Nov 2026; content at abstract level unless stated)

- **Hybrid inflation:**
  - Linde astro-ph/9307002;
  - Copeland–Liddle–Lyth–Stewart–Wands astro-ph/9401011;
  - smooth hybrid, Lazarides–Panagiotakopoulos hep-ph/9506325;
  - shifted hybrid, Jeannerot–Khalil–Lazarides–Shafi hep-ph/0002151;
  - F-term log, Dvali–Shafi–Schaefer hep-ph/9406319.
- **Mild waterfall and PBHs:**
  - Clesse 1006.4522;
  - Kodama–Kohri–Nakayama 1102.5612;
  - Clesse–García-Bellido 1501.07565;
  - Kawasaki–Tada 1512.03515;
  - stochastic δN: Fujita–Kawasaki–Tada–Takesako 1308.4754, Pattison–Vennin–Assadullahi–Wands 1707.00537;
  - PBH constraints, Carr–Kohri–Sendouda–Yokoyama 2002.12778.
- **Topological inflation:** Vilenkin hep-th/9402085; Linde astro-ph/9402031; Sakai et al. gr-qc/9506068 (η_c ≈ 0.33 m_Pl).
- **Preheating:**
  - tachyonic: Felder et al. hep-ph/0012142; Felder–Kofman–Linde hep-th/0106179; García-Bellido–García Pérez–González-Arroyo
    hep-ph/0208228;
  - tachyonic resonance with trilinears: Dufaux et al. hep-ph/0602144;
  - broad and stochastic resonance: Kofman–Linde–Starobinsky hep-ph/9704452; Greene et al. hep-ph/9705347;
  - fermions: Greene–Kofman hep-ph/9807339;
  - instant preheating: Felder–Kofman–Linde hep-ph/9812289.
- **Equation of state:** Podolsky–Felder–Kofman–Peloso hep-ph/0507096; Lozanov–Amin 1608.01213 and 1710.06851.
- **Higgs:**
  - Lebedev–Westphal 1210.6987;
  - Herranen–Markkanen–Nurmi–Rajantie 1407.3141;
  - Ema–Mukaida–Nakayama 1602.00483;
  - threshold stabilisation: Elias-Miró–Espinosa–Giudice–Lee 1203.0237, Lebedev 1203.0156.
- **Data:** Planck 2018 X 1807.06211; BICEP/Keck 2110.00483; ACT DR6 2503.14452 and 2503.14454.
- **Our tracks, read-only:** inf_reheat (width convention, S2, Greene–Kofman check), inf_critic, inf_chain, dq_entropy
  (S/Q closed form), dq_critic, WORKING_THEORY §4b (R_vec ∝ K²).
