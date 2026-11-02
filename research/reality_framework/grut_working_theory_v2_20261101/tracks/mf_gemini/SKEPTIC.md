# SKEPTIC of track mf_gemini (audit of the external two-field waterfall proposal)

2 November 2026 · Claude, adversarial skeptic in the multi-field inflaton pass (BRIEF_MULTIFIELD.md; GitHub issue #13
context). Key: `mf_gemini_skeptic`. Code and results: `tracks/mf_gemini/skeptic/{code,results}`.

This is same-session AI work on one supplied model, not outside replication. Nothing here is a validation, an exclusion
of the real world, a new theory or a GRUT prediction. **Every potential shape, coupling, slice convention, loop
criterion, collapse threshold and ledger used below is a model choice, not a consequence of closure.** Data comparisons
(Planck 2018, BICEP/Keck 2021, ACT DR6) are retrospective. Git was read-only and nothing was committed. All runs were in
the foreground, one process at a time (about 40 min of CPU in total). Scratch:
`$SP/wt2/mf_gemini_skeptic/`.

Labels: [D] derived, [C] computed with this skeptic's own code, [E] estimate, [Q] quoted from a reviewed track,
[A] argued, [R] literature (arXiv IDs and titles checked against the arXiv API on 2 Nov 2026).

Read in full: BRIEF_MULTIFIELD.md; WORKING_THEORY.md §0, §2, §4a, §4b, §7–§9; BRIEF_INFLATON.md; BRIEF_DECISIVE.md; the
track's NOTES.md, all nine scripts and all eight results files. Read in the parts that matter: inf_critic/NOTES.md, the
inf_chain NOTES and SKEPTIC (slice convention, mirror branch), the inf_reheat NOTES and SKEPTIC (width convention,
Greene–Kofman), the inf_potential skeptic (K convention), and the dq_heat and dq_entropy NOTES (K = 1.98581, the S/Q
closed form).

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

Every check below uses **my own code with a different method** (§1). "Stands" means I reproduced the number
independently; "sharpened" means it stands but a stated reason or range changes; "refuted" means a concrete recomputed
number contradicts it.

| Track result | Verdict | Deciding recomputed number | Label |
|---|---|---|---|
| G1 amplitude, blocking, loop | **stands** | Φ = 1.225 / 0.816 / 0.612 M_P; y(Γ = 0.5H) = 2.894; yΦ/(m/2) = 2.95×10⁷; loop/tree = 1.16×10¹³. First rest point after the end: χ = 1.1474 v, H = 0.3166 H0 (the quadratic identity gives 0.295, so the anharmonic correction is 7%) | D/C |
| G2a v, λ, loop, blocking | **stands, one inconsistency** | v = 1.63299 M_P, λ = 4.810×10⁻¹⁵ (solved numerically, not from the closed form). Renormalisation-condition loop caps: y ≤ 2.95×10⁻⁴ (1% change of the tachyonic curvature at χ = 0) to 5.5×10⁻⁴ (max ΔV ≤ 1% V0), matching the track's 3.5–5.5×10⁻⁴. **But at y = 5.5×10⁻⁴, Γ/H = 1.8×10⁻⁸, not the track's 7×10⁻⁹** (that value belongs to the Z2 cap 3.5×10⁻⁴) | C |
| G2b duration of the mild waterfall | **stands** (independent code, < 0.5%) | N_wf (mean) = 17.38 / 18.48 / 21.15 / 28.14 (quad Δ = 1 / 0.3 / 0.1 / 0.03) against the track's 17.45 / 18.48 / 21.15 / 28.14; lin Δ = 0.1 / 0.06: 20.47 / 21.97 against 20.50 / 21.97 | C |
| G2b, concave V_SR (new) | **sharpened** | For a constant slope or a hilltop, a fraction of patches **never completes the waterfall**: σ leaves the band first and χ decays back. That fraction is 3.3×10⁻³ (lin Δ = 0.1), 4.4×10⁻² (hill Δ = 0.1) and 2×10⁻⁵ (hill Δ = 0.06). The track extrapolates these patches with the log law; they actually keep inflating | C |
| G2c P_ζ, δN tail, PBH | **stands; f_PBH overstated ×2.6** | Var(δN) = 1.023, P(ζ > 0.67) = 0.2007, P_ζ,peak = 0.567 (quad Δ = 1; the track gives 1.025 / 0.2005 / 0.561). f_PBH = 3.6×10¹¹–1.1×10¹⁶ in instant radiation, 1.8×10¹⁰–5.6×10¹⁴ with the slow reheating the track itself finds, and ≥ 2×10¹⁰ from the no-waterfall patches alone. The track's (T_f/T_eq)(g_f/g_eq)^{1/3} form omits a factor 0.389 | C/E |
| G2d walls | **stands** | σ_w = 5.333 H0 M_P²; walls dominate once H < 1.78 H0; v/m_Pl = 0.3257, i.e. 0.987 of 0.33 (gr-qc/9506068 abstract: "η ≳ 0.33 m_Pl") | D/R |
| G2e holder on (χ − v) | **stands, strengthened** | In −v domains the late \|μ\| = 21.24 H exceeds the Y_L threshold M_Y²/(√2 h) = 15.91 H, so **m²_YL = −302 H²**. The −v domain ends in a Y-condensed hidden vacuum, not merely the opposite CP sign. Along χ: 0 → −v this happens at \|χ\| = 0.54 v | D |
| REQ (does any choice work?) | **stands** | Every loophole I tried closes (§6): large N_f needs N_f ≥ 7.5×10¹⁴, which puts the species cutoff at 8.9×10¹⁰ GeV = 0.68 H; the Higgs routes of either sign hit a Landau pole at 9×10¹² GeV < V0^{1/4}; the bosonic channel with a consistent decay gives K ≥ 2.2 at best | C/D |
| HIGGS | **stands, reason sharpened** | The threshold 2.094 and the Landau pole at 9.0×10¹² GeV do **not depend on the portal sign**. The opposite-sign valley reaches \|H\| = 3.5×10¹⁴ GeV at the overshoot, 40 times above the pole, so that route fails on perturbativity whatever λ_SM does (the track's stated reason, boundedness for λ_SM < 0, depends on M_t). The Higgs loop on χ is 0.33% of V0, within the 1% criterion | D/C |
| PREHEAT (bosonic, loop cap) | **sharpened** | No-decay regression: 50% in φ at 9.4 H_end⁻¹ (max growth) and 15.5 (random phases), against the track's 9.7 and 13.4. The track's model takes maximal Bose growth **and** instant decay to the SM, which exclude each other. With a decay rate Γ_φ, **radiation** reaches 50% at 11.7 (Γ_φ = 0.3 H_end) or 16.3 (Γ_φ = 1) H_end⁻¹ for max growth; with random phases 21.7 (0.3) and **never within 40 H_end⁻¹** (Γ_φ ≥ 1). For max growth and Γ_φ ≥ 3 it is ≤ 2×10⁻⁶ by 40 H_end⁻¹. K ≈ 2.2–2.9 survives only for Γ_φ ≈ 0.3 H_end | C on a model |
| NS (n_s claim) | **partly refuted: the log-family row** | Exact Mukhanov–Sasaki n_s for the log family with its waterfall attached consistently (finite σ_c): **0.9720 / 0.9721 / 0.9739 / 0.9850** at Δ = 1 / 0.3 / 0.1 / 0.03 (N_wf = 17.5 / 18.2 / 20.3 / 26.8). The track's 1 − 1/N_σ gives 0.9717 / 0.9711 / 0.9692 / 0.9616. Slower crossings make n_s **bluer**, not redder; the family never goes below 0.972 and is in the window only for Δ ≳ 0.2 (N_wf ≲ 19). The track's ns2 also drops the ξ² term (+2.13 ξ² = +8.6×10⁻⁴ here). The convex, hilltop and r statements stand: \|η\| = 0.01356–0.01659, r = 2.75×10⁻⁷ | C |
| G3 heavy waterfall | **stands** | λ = 0.01306, m = 3850 H, v = 1.27×10⁻³ M_P (closed form re-derived) | D |
| SQK (K_crit) | **stands; 1.38 → 1.41** | When the Y deposits are also rescaled as K^{1/2} (inf_potential skeptic, scenario B), K_crit = 1.405. A constant Γ ≥ 0.36–0.39 H_end from ε_H = 1 would reach it; the proposal's own Γ = 0.5–1 H0 would give K = 1.17–1.24. **The only failing link is the decay channel** | C/E |
| BIAS | **stands** | The holder tadpole shifts χ by only 2.4×10⁻⁴ H, against the 5.3 H needed | D/E |

**Bottom line.** The track's main conclusions survive independent recomputation:
- G1 and G2;
- a quantum-seeded waterfall of at least 17 e-folds with an O(1) δN tail;
- PBH overclosure and walls for real χ;
- no channel within the proposal's content that reaches Γ_eff ≥ 0.5 H.

**What I refute.** The log-family n_s mapping is wrong in sign. When the waterfall is attached consistently, slower
crossings move n_s to 0.974–0.998, not 0.92–0.97. The family's n_s is then pinned at ≈ 0.972, the upper edge of the
claimed window.

**What I sharpen.**
- The Higgs failure, the preheating optimism, the f_PBH normalisation and the Γ at the loop cap.
- Concave V_SR shapes, the ones that give the red tilt, leave patches that never finish the waterfall.
- The −v domains of the holder on (χ − v) are Y-tachyonic.

---------------------------------------------------------------------------------------------------------------------

## 1. What is independent

| Track method | This skeptic | File |
|---|---|---|
| closed forms for v, λ | numerical solve of m = 3H, λv⁴/4 = 3M_P²H² | `sk_g1g2.py` |
| least-squares counterterm fit, max residual | renormalisation conditions (V(0) = V0, V(v) = 0, V′(v) = 0; for the shifted form also V′(0) = 0 and V″(v) = m²); measures: change of V″(0), of V″(v), and max ΔV | `sk_g1g2.py` |
| DOP853 in cosmic time; WKB modes at k/aH = 100 | e-fold time; adiabatic vacuum (with the dω/dN term) at k/aH = 50, integrated as v = a^{3/2}u; deeper heavy pre-history | `sk_waterfall.py` |
| σ WKB-averaged once gχ > 30H | σ integrated directly until gχ = 100H, then **dropped** (its energy fraction there is 10⁻²³–10⁻⁹) | `sk_waterfall.py` |
| band exits extrapolated with the log law | band exits detected and counted as no waterfall | `sk_waterfall.py` |
| Planck closed form for N_* | direct horizon matching with entropy conservation | `sk_ns.py` |
| second-order slow roll without ξ² | slow roll with ξ², **plus exact Mukhanov–Sasaki modes** on the numerical σ background | `sk_ns.py` |
| background, K by DOP853 | LSODA with a radiation fluid; **regression: M1r's (3, 0.5) gives K = 1.985805** (dq_heat: 1.98581) | `sk_background.py` |
| KLS map, no decay, "instant decay" | KLS map with a consistent φ decay into a radiation fluid | `sk_preheat.py` |

**Regression of the mode code [C].** In exact de Sitter after 3 e-folds, P_χ is accurate to 3.5×10⁻⁴ (massless) and
3×10⁻⁴ (m² = 2H²); the error is set by the k/aH = 50 start. The m² = −4.5H² growth rate is 1.09795, against the exact
1.09808.

**Two bugs found and fixed in my own code during the run** (no reported number uses the faulty versions):
1. **Mode starts.** An early version started modes deep inside the horizon. The fix starts each mode at
   max(k/aH = 50, background start) in the heavy adiabatic vacuum.
2. **Integration variable.** Integrating u itself lets |u| fall below atol during a long heavy pre-history. The explicit
   integrator then went unstable and gave a spurious "N_wf = 6.7" for log Δ = 0.3. The v = a^{3/2}u variable cures
   this, and all twelve cases were re-run with it.

---------------------------------------------------------------------------------------------------------------------

## 2. G1 and G2 numbers [D/C] (`results/sk_g1g2.json`)

**G1.**
- All identities reproduce (table in §0).
- M1r's oscillator starts at rest. The analogous point on the proposal's background is the first rest point
  (χ = 1.1474 v), where the instantaneous H = 0.3166 H0. The quadratic identity gives 0.295 H0 there; the 7% gap is the
  anharmonicity of the double well [C].

**G2.** For y = 2 (H0 = 1.3×10¹¹ GeV), using renormalisation conditions:

| Form of m_N | ΔV″(0)/V″_tree(0) | ΔV″(v)/V″_tree(v) | max\|ΔV\|/V0 |
|---|---|---|---|
| m_N = yχ (Z2) | 2.1×10¹³ | −2.1×10¹³ | 4.3×10¹² |
| m_N = y(χ − v) (bare mass tuned) | 2.1×10¹³ | 1.1×10⁷ | 1.7×10¹² |

- **The 1% caps** are y ≤ 2.95×10⁻⁴ on the waterfall-rate measure and 4.4–5.5×10⁻⁴ on max|ΔV|, at every H0. That
  confirms the track's 3.5–5.5×10⁻⁴. The waterfall rate, which sets the 17-e-fold floor, is the most sensitive
  measure.
- **The width at the cap.** At y = 5.5×10⁻⁴ the Majorana width is Γ/H = 1.8×10⁻⁸. The track's perturbative row
  ("y ≤ 5.5×10⁻⁴, giving Γ ≤ 7×10⁻⁹ H") mixes the two caps. The factor of 2.5 changes nothing.
- **The odd terms [E].** With the bare mass tuned (m_N = y(χ − v)) the N loop is χ-odd. Its tadpole at χ = 0 is of
  order y⁴v³/8π² × ln ≈ 3×10⁷ H³ at the cap, which would displace χ by about M_P during inflation. Keeping χ = 0 a
  stationary point (V′(0) = 0, imposed above) is therefore a second renormalisation tuning of about 10⁻⁷ on top of the
  2×10⁻⁸ mass tuning. The "bias" of §6 of the track is then a choice of this renormalisation condition, not an extra
  term.

---------------------------------------------------------------------------------------------------------------------

## 3. The mild waterfall, recomputed [C] (`results/sk_waterfall.json`)

The slope at σ_c is the track's p_c = 3.18×10⁻⁴ and H0 = 1.3×10¹¹ GeV. Δ = p_c/σ_c. The matching criterion is the
track's (χ_rms = 10⁻² χ_trap); 10⁶ Monte Carlo draws of χ_L.

| case | N_wf classical / mean | s_eff | Var(δN) | P(ζ > 0.67) | P_ζ,peak (linear δN) | peak exit (e-folds before end) | P(no waterfall) | track N_wf mean / Var / tail / P_ζ |
|---|---|---|---|---|---|---|---|---|
| quad Δ = 1 | 16.80 / 17.38 | 1.098 | 1.023 | 0.2007 | 0.567 | 15.94 | 0 | 17.45 / 1.025 / 0.2005 / 0.561 |
| quad Δ = 0.3 | 17.90 / 18.48 | 1.098 | 1.023 | 0.2007 | 0.358 | 16.90 | 0 | 18.48 / 1.026 / 0.2005 / 0.358 |
| quad Δ = 0.1 | 20.55 / 21.15 | 1.059 | 1.082 | 0.2055 | 0.252 | 19.65 | 0 | 21.15 / 1.084 / 0.2055 / 0.252 |
| quad Δ = 0.03 | 27.43 / 28.14 | 0.868 | 1.558 | 0.2336 | 0.219 | 26.38 | 0 | 28.14 / 1.561 / 0.2335 / 0.218 |
| lin Δ = 0.1 | 19.69 / 20.47 | 0.943 | 2.18 | 0.2161 | 0.315 | 18.69 | **3.3×10⁻³** | 20.50 / 2.32 / 0.2141 / 0.315 |
| lin Δ = 0.06 | 21.39 / 21.97 | 1.096 | 1.047 | 0.2006 | 0.188 | 20.39 | < 10⁻⁶ | 21.97 / 1.048 / 0.2004 / 0.188 |
| hill Δ = 0.1 (η = 0.015) | 19.87 / 20.66 | 0.807 | 1.86 | 0.2413 | 0.430 | 18.76 | **4.4×10⁻²** | — |
| hill Δ = 0.06 | 21.15 / 21.75 | 1.087 | 1.14 | 0.2015 | 0.190 | 20.04 | **2×10⁻⁵** | — |
| log Δ = 1 (s0 = σ_c) | 16.95 / 17.53 | 1.098 | 1.023 | 0.2007 | 0.493 | 15.90 | 0 | — |
| log Δ = 0.3 | 17.66 / 18.24 | 1.098 | 1.023 | 0.2007 | 0.348 | 16.59 | 0 | — |
| log Δ = 0.1 | 19.75 / 20.33 | 1.092 | 1.030 | 0.2014 | 0.235 | 18.79 | 0 | — |
| log Δ = 0.03 | 26.13 / 26.79 | 0.945 | 1.319 | 0.2218 | 0.184 | 25.16 | 0 | — |

**Readings.**
- **Reproduced [C].** The track's durations, variances, tails and peak spectra agree to ≤ 0.5%, with a different time
  variable, mode start, integrator and treatment of σ.
- **Insensitivity [C].** Including the sub-horizon tachyonic band (k/a < |m_χ|) in χ_rms changes nothing at the 10⁻⁴
  level. σ's energy fraction where I drop it is 10⁻²³–10⁻⁹, so G4's "σ holds V0" indeed does not apply.
- **Concave V_SR [C, new].** For lin and hill, patches with |χ_L| below about 4×10⁻³ (lin, Δ = 0.1) or 6×10⁻² (hill,
  Δ = 0.1) χ_rms do not finish the waterfall in this crossing:
  - σ leaves the band, χ turns heavy and decays back (traced explicitly: χ peaks at 3.5×10⁻⁴ M_P, then falls as e^{−N});
  - those regions go on inflating on V0 + V_SR;
  - the track's code labels them as the "band exit" fraction of its c grid (39% of the grid at lin Δ = 0.1). It then
    assigns them a finite N by the log-law extrapolation, which is why its Var(δN) there (2.32) differs from the
    analytic value.

  The hilltop is the family the track uses for the n_s window. On it this fraction is 4.4% at Δ = 0.1. It is a
  separate catastrophe on top of the δN tail: long-inflating regions inside the observable universe.

---------------------------------------------------------------------------------------------------------------------

## 4. PBHs [E] (`results/sk_pbh.json`)

**Normalisation.** I recompute f_PBH from entropy conservation, ρ_PBH,0 = β(3/4)T_f s_0. The track's
f = β(T_f/T_eq)(g_f/g_eq)^{1/3}(Ω_m/Ω_DM) leaves out the change in the radiation density across annihilations,
(g_f/g_eq)(g_sf/g_seq)^{−4/3} = 0.389. Its f is therefore 2.6 times too large [D]. This is immaterial at f ≫ 1.

**Formation epoch.**
- The track forms PBHs in "instant radiation", but its own result is slow reheating.
- At the perturbative cap (Γ ≈ 1.8×10⁻⁸ H_end) the matter era lasts N_MD = (2/3) ln(H_end/Γ) ≈ 11.9 e-folds.
- The peak modes (N_b ≥ 16) then re-enter after reheating, at H_entry = H_end e^{−2N_b − N_MD/2}.
- That raises M_PBH by ×380, to 7×10¹⁸–9×10²⁷ g, and lowers f by ×20, to 1.8×10¹⁰–5.6×10¹⁴ [E].

**No-waterfall patches.** Counting each of them only once, as at least a horizon-mass black hole, gives f ≥ 1.1×10¹³
(lin, Δ = 0.1), 1.4×10¹⁴ (hill, Δ = 0.1) and 2×10¹⁰ (hill, Δ = 0.06) [E].

**Survival.** Every M_PBH is above 5×10¹⁴ g, so the PBHs survive Hawking evaporation to today.

**Verdict.** Overclosure is robust in every variant, by ≥ 10¹⁰. The δN-tail regions are sheet-like, sitting next to
the χ = 0 surface, rather than isolated spherical peaks. So "PBH" here partly means collapsing walls, as the track
says [A]; the conclusion does not depend on which.

---------------------------------------------------------------------------------------------------------------------

## 5. n_s: the log-family row is wrong in sign [C] (`results/sk_ns.json`)

**What the track did.** It attaches V_SR = V0 A ln σ to the waterfall using n_s = 1 − 1/N_σ, which assumes σ_c ≪ σ_*.
But the waterfall's Δ is itself a property of V_SR: Δ = A/σ_c² for pure log, or A/(2σ_c²) with s0 = σ_c. So
σ_c²/(2A) = 1/(2Δ) is **not** small for slow crossings, and the effective e-fold number is N_σ + 1/(2Δ).

**The consistent family.** Take V = V0[1 + (A/2) ln((σ² + σ_c²)/σ_c²)] with σ_c fixed by A_s and N_σ = N_* − N_wf(Δ).
Use N_wf from my own log-family waterfall runs, N_* = 52.82 (K = 2; my direct matching gives 53.17 / 52.82 / 52.62
for K = 1 / 2 / 3, matching the track), and exact Mukhanov–Sasaki modes.

| Δ | N_wf | N_σ | σ_c / σ_* (M_P) | n_s exact (MS) | n_s, 2nd order with ξ² | without ξ² | track formula 1 − 1/N_σ | α_s (MS) |
|---|---|---|---|---|---|---|---|---|
| 1 | 17.53 | 35.29 | 1.09×10⁻³ / 0.0128 | **0.9720** | 0.9725 | 0.9716 | 0.9717 | −7.7×10⁻⁴ |
| 0.3 | 18.24 | 34.59 | 1.98×10⁻³ / 0.0123 | **0.9721** | 0.9724 | 0.9715 | 0.9711 | −7.5×10⁻⁴ |
| 0.1 | 20.33 | 32.49 | 3.42×10⁻³ / 0.0116 | **0.9739** | 0.9740 | 0.9733 | 0.9692 | −6.1×10⁻⁴ |
| 0.03 | 26.79 | 26.03 | 7.15×10⁻³ / 0.0124 | **0.9850** | 0.9850 | 0.9850 | 0.9616 | −2×10⁻⁵ |
| 0.01 (track's N_wf = 40.7) | 40.7 | 12.1 | 0.0187 / 0.0209 | ≈ 0.998 (2nd order) | 0.998 | 0.998 | 0.918 | +1.4×10⁻⁴ |

The MS amplitude matches A_s to 1.5–2.3%.

**Readings.**
- **The direction is reversed [C].** The track's row "N_wf = 17.4–40.7 → n_s = 0.972–0.92" is wrong. Consistently,
  n_s = 0.972 → 0.998: slow crossings exit the window on the **blue** side.
- **Where the window holds.** Only for Δ ≳ 0.2, i.e. N_wf ≲ 19 (the track says ≤ 22.7). There n_s = 0.9720–0.9721, at
  the top of the claimed 0.967–0.973 window. The family cannot give n_s < 0.972.
- **Retrospective pulls of n_s = 0.972:** +1.7σ against Planck 2018 (0.9649 ± 0.0042, 1807.06211) and −0.7σ against
  ACT DR6 P-ACT-LB (0.9743 ± 0.0034; abstract "0.974 ± 0.003", 2503.14452).
- **The ξ² term.** The track's `ns2` multiplies it by 0.0. The Stewart–Lyth term −2(C − 1/3)ξ² = +2.13ξ² shifts n_s by
  +8.6×10⁻⁴ for the log family, because ξ² ≈ 2η² there. The exact MS values confirm that keeping it is right, to
  ≤ 5×10⁻⁴.
- **Running.** α_s = −(6–8)×10⁻⁴ lies inside ACT's dn_s/dln k = 0.0062 ± 0.0052 (2503.14454) [R].
- **What stands [C].** Convex V_SR is blue. The hilltop window is |η| = 0.01356–0.01659, independent of N_σ.
  r = 2.54 / 2.75 / 3.91×10⁻⁷ at H0 = 1.25 / 1.3 / 1.55×10¹¹ GeV.
- **A compatibility point [C/A].** The shapes that give the shortest waterfalls (convex quad or log near σ_c, which
  trap σ) are not the ones that give a red tilt by curvature alone. The red-tilt hilltop is concave in the band, which
  produces the no-waterfall patches of §3. The log family is red at σ_* and convex near 0, but is pinned at
  n_s ≥ 0.972.

---------------------------------------------------------------------------------------------------------------------

## 6. The existence question (item 2): loopholes checked [C/D]

The task: find any (y_N, λ_χH, g_hyb, V_SR, holder) that meets m_χ = 3H_i, λv⁴/4 = V0, Γ_eff ≥ 0.5H, H_i in range, the
n_s window, no defects and a working holder. **None found.** I tried to break the track's "no":

1. **Many N flavours.** At fixed Γ the loop scales as Γ²/N_f, so a large N_f evades it. Γ = 0.5H with
   N_f y⁴ ≤ y_cap⁴ needs N_f ≥ 7.5×10¹⁴ (y = 1.06×10⁻⁷). That is radiatively fine and even kinematically open
   (yΦ/(m/2) = 0.32). But the species cutoff M_P/√N_f is 8.9×10¹⁰ GeV = 0.68 H, below the Hubble scale (Dvali
   0706.2050) [D/R]. **Closed.**
2. **Higgs portal, opposite sign (tachyonic Higgs preheating).** Its windows are an extremely efficient
   non-perturbative transfer into the SM (exponent 4.3×10³), so I checked it hardest.
   - The Higgs loop on χ is 0.33% of V0 after counterterms, which passes the 1% criterion [C].
   - The flat direction χ² − v² = 2|λ_χH||H|²/λ_χ reaches |H| = 3.5×10¹⁴ GeV at the 0.147 v overshoot.
   - The tree threshold λ_χH²/λ_χ = 2.094 holds for **either sign**, so λ_H(UV) ≥ 2.09 and the one-loop Landau pole sits
     at 9.0×10¹² GeV, 40 times below those field values and at 1.2% of V0^{1/4}.
   - The portal also runs λ_χ by 10.6% per e-fold of scale.
   - **Closed by perturbativity**, independently of the SM metastability the track invokes.
3. **The bosonic channel with a consistent decay** (§7): K ≥ 2.2 only for Γ_φ ≈ 0.3 H_end; for faster decays there is no
   transfer. It also adds a field and a Z2-odd trilinear. **Not a two-field solution.**
4. **Granting the claimed Γ.** On the proposal's own background, a constant Γ = 0.5 / 1 H0 switched on at ε_H = 1 would
   give K = 1.245 / 1.166 (a = 1 at ε_H = 1), below K_crit = 1.38–1.41. Γ ≥ 0.36–0.39 H_end suffices [C]. Even then the
   release is not M1r's (§8), so the S/Q claim stays unsupported, but **the channel is the only failing link** in the K²
   estimate.
5. **The proposal's own couplings at the perturbative cap** [E]. K ≈ (H_end/Γ)^{1/3} ≈ 340–520, like TB1/TB2's K ≈ 2000
   on the S1 normalisation. In the inflaton pass that gave free X with S/Q of order 10²–10³ [Q inf_critic], i.e. the
   symmetric-relic outcome, not the proposal's.
6. **Defects** need a χ-odd bias (or the Clesse regime N_wf > N_*). With the bare mass tuned, the bias is a
   renormalisation condition (§2). It removes walls and PBHs but does nothing for reheating.

---------------------------------------------------------------------------------------------------------------------

## 7. Bosonic preheating: the "instant decay" is inconsistent with the Bose growth [C on a model] (`results/sk_preheat.json`)

The model is the track's: ½g²(χ − v)²φ² at the loop cap g² = 4.32×10⁻⁷ (q ≈ 2.4×10⁵), on my own background, no
back-reaction. Times to F = 0.5 are in H_end⁻¹. For Γ_φ = 0, F is the fraction in φ (the track's quantity); otherwise it
is the fraction in radiation.

| Γ_φ/H_end | max growth, cap | max growth, 100×cap | random phases, cap | random, 100×cap |
|---|---|---|---|---|
| 0 (φ fraction) | 9.4 | 6.5 | 15.5 | 11.1 |
| 0.3 | 11.7 | 8.3 | 21.7 | 15.5 |
| 1 | 16.3 | 11.6 | **never by 40** (max 3×10⁻³) | 42.7 |
| 3 | **never by 40** (≤ 1.6×10⁻⁶) | never (5×10⁻⁴) | — | — |
| 10 | never (≤ 4×10⁻⁷) | never (1.3×10⁻⁴) | — | — |

Readings:
- The no-decay row reproduces the track: 9.4 against 9.7, and 15.5 against 13.4 with a different random seed.
- **Decay and Bose enhancement compete.** Decay faster than about one crossing per H_end⁻¹ resets n_k, so the transfer
  collapses to the first-crossing fraction of about 10⁻⁷ per crossing.
- The best case is 50% radiation after 11.7 H_end⁻¹ (Γ_φ = 0.3 H_end, maximal phases), which gives K ≈ 2.2. Random
  phases give 21.7 H_end⁻¹, i.e. K ≈ 2.8.
- The track's "K ≈ 2.1–2.9 not excluded" therefore holds only in a narrow Γ_φ window. It is never below the
  1.38–1.41 that Q-matter needs [C/E]. A lattice would only add back-reaction and rescattering, which slow the transfer
  further [A].

---------------------------------------------------------------------------------------------------------------------

## 8. The holder on (χ − v): the −v domains are Y-tachyonic [D] (`results/sk_g1g2.json`)

**M1r's light-Y mass.** m²_YL = M_Y² − √2|μ|h reproduces ω_YL = 22.249 H at μ = 7.16 (BRIEF_DECISIVE.md).

| Stage | μ | m²_YL |
|---|---|---|
| during inflation | μ_bare = 10i | +334 H² |
| late, +v domain | μ = 7.16 | +495 H² |
| late, −v domain | −7.16 + 20i | **−302 H²** |

- **The threshold** is |μ| = M_Y²/(√2h) = 15.91 H. Along χ: 0 → −v it is crossed at |χ| = 0.54 v, during the last part
  of the waterfall.
- **The −v domain** therefore relaxes to a vacuum with a charge-3 condensate Y. That is a different hidden phase, not
  the mirror of the +v domain.
- **Precedent.** The inf_chain skeptic found the analogous Y condensation for TB's mirror branch at the slice
  (|u| = 17.47).
- **What it changes.** The track's "opposite CP phase (−0.34 vs +π/2)" understates the problem. With a bias, only +v
  exists and this is moot.

**The tadpole.** The holder tadpole 2g|Σ||Y|² ≈ 1.1×10⁻³ H³ against the tachyonic curvature 4.5 H² gives a static
shift of 2.4×10⁻⁴ H. That is negligible against the 5.3 H a bias needs, and confirms the track [D].

---------------------------------------------------------------------------------------------------------------------

## 9. Conventions (slice, H_i, K) [C] (`results/sk_background.json`)

- **Where m_χ/H is measured.** It depends on the slice:

  | Slice | H | m_χ/H |
  |---|---|---|
  | S1, one e-fold before ε_H = 1 (inflaton pass) | 0.9715 H0 | 3.09 |
  | ε_H = 1 (the track) | 0.7126 H0 | 4.21 |
  | first crossing of v | 0.402 H0 | 7.46 |
  | first rest point | 0.3166 H0 | 9.47 |

- **Phase-matched comparison.** M1r's own oscillator (m = 3 at rest) has m/H = 3.34 at its first ε_H = 1 and 6.14 at its
  first zero crossing. So, matched in phase, the proposal is equivalent to **m_I ≈ 3.6–3.8 H_i** in M1r's engine
  units. That is less marginal than the track's "4.21 vs inf_reheat's m_I = 4 marginal" suggests. The release is still
  not M1r's.
- **The dark matter fixes H on the slice.** H_i ≈ 1.29×10¹¹ GeV is the Hubble rate at M1r's slice. On the proposal's
  background that means H0 = H_i/0.7126 = 1.81×10¹¹ GeV, so r = 5.3×10⁻⁷ and V0^{1/4} = 8.7×10¹⁴ GeV. The stated range
  [1.25, 1.55]×10¹¹ on the slice maps to H0 = 1.75–2.18×10¹¹ and r = 5.0–7.7×10⁻⁷. The track's r = 2.5–3.9×10⁻⁷ uses
  H0 as H_i. Both are far below 0.036.
- **K depends on the normalisation.** For a constant Γ = 0.5 H0: K = 1.245 (a = 1 at ε_H = 1), 6.75 (a = 1 at S1),
  1.145 (a = 1 at the first rest point). Only the ε_H = 1 normalisation matches M1r's K = 1.986 (inf_potential
  skeptic); mixing conventions costs a factor of about 5.

---------------------------------------------------------------------------------------------------------------------

## 10. New issues (not in the track)

1. **Y-tachyonic −v domains** of the holder on (χ − v): m²_YL = −302 H². The two vacua differ in hidden phase, not
   only in CP sign [D].
2. **No-waterfall patches for concave V_SR:** 3.3×10⁻³ (lin, Δ = 0.1), 4.4×10⁻² (hill, Δ = 0.1), 2×10⁻⁵ (hill,
   Δ = 0.06). They keep inflating; the track extrapolates them [C].
3. **The log-family n_s runs blue for slow crossings** (0.974–0.998, not 0.92–0.97), is pinned at ≥ 0.972, and sits in
   the window only for N_wf ≲ 19. The ξ² term is dropped in the track's ns2 [C].
4. **The Landau pole is sign-independent.** The opposite-sign Higgs route fails on perturbativity (|H| = 3.5×10¹⁴ GeV
   against a pole at 9×10¹² GeV) whatever M_t is [D/C].
5. **The bosonic channel with a consistent decay** gives radiation at 50% only after 11.7–21.7 H_end⁻¹, for Γ_φ ≈ 0.3
   H_end only; nothing by 40 H_end⁻¹ for faster decays [C].
6. **The large-N_f loophole** is closed by the species bound (cutoff 0.68 H) [D].
7. **f_PBH formula:** a missing g_* radiation factor 0.389 (the track's value is ×2.6 too large). Slow-reheating
   formation gives M ×380 and f ×1/20. Overclosure holds throughout [D/E].
8. **Slice and H_i convention:** with the dark matter fixing H on M1r's slice, H0 = 1.81×10¹¹ GeV and r = 5.3×10⁻⁷.
   K changes by about ×5 between the S1 and ε_H = 1 normalisations [C].
9. **Granting the claimed Γ = 0.5–1 H0 gives K = 1.17–1.24 < K_crit.** The channel, not the history, is the failing
   link. K_crit = 1.405 when the Y deposits also scale as K^{1/2} [C/E].
10. **The bias as a renormalisation condition.** With m_N = y(χ − v), χ = 0 needs a second tuning, of the tadpole, to
    about 10⁻⁷ [E].
11. **Inherited, not computed here [A/Q]:**
    - With a long-lived χ, any induced χ–hidden portal (holder-induced A_χΣ) competes with SM channels capped at
      Γ ≤ 1.8×10⁻⁸ H. That needs the inflaton pass's |A| ≲ 10⁻⁶ H tuning to avoid R_leak.
    - Gravitational ξ_X = −5/6 production runs through the whole long χ-oscillation era.
    - Scalar-induced GWs from P_ζ ~ 0.2–0.6 peak near 1 Hz for the fast crossings. They are unconstrained there and
      immaterial given PBH overclosure [E].

## 11. Not checked

- No lattice: fragmentation, oscillons and back-reaction are not computed. The double-well anharmonicity at the
  overshoot is 0.15, and m/H rises from 4 to 10 over the first swing [A].
- No stochastic δN of my own; the track's Langevin run agrees with both of our nonlinear δN tails.
- No release of the M1r chain on this history (mf_ideal / mf_chain machinery).
- The log-family waterfall runs used the track's p_c (3.18×10⁻⁴), not the CMB-consistent p_c = p_*(r + 1/r)/2
  ≈ (4–7)×10⁻⁴. The track found N_wf changes by < 3% for a change of p at fixed Δ.
- SM running was not redone. The Landau pole is one-loop, from 24λ²/16π².

## 12. Files and reproduction

| file | content |
|---|---|
| `skeptic/code/sk_g1g2.py` → `results/sk_g1g2.json` | G1/G2 by numerical solve, renormalisation-condition loop caps, Higgs (both signs, Landau pole, Higgs loop on χ), large-N_f/species, −v domain Y mass, K_crit with Y deposits |
| `skeptic/code/sk_background.py` → `results/sk_background.json` | LSODA background, M1r K regression (1.985805), slices S1 / ε_H = 1 / v-crossing / rest point, K for sudden and constant-Γ decay in three normalisations, Γ for K_crit |
| `skeptic/code/sk_waterfall.py` → `results/sk_waterfall.json` | independent linear modes + nonlinear δN for 12 cases (quad, lin, hill, log), band-exit detection, de Sitter regression |
| `skeptic/code/sk_ns.py` → `results/sk_ns.json` | direct N_*, hilltop window, consistent log family with exact Mukhanov–Sasaki n_s and α_s |
| `skeptic/code/sk_pbh.py` → `results/sk_pbh.json` | f_PBH by entropy conservation, slow-reheating formation, no-waterfall patches |
| `skeptic/code/sk_preheat.py` → `results/sk_preheat.json` | KLS crossing map with a consistent φ decay into radiation |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_gemini/skeptic/code
python3 -B sk_g1g2.py; python3 -B sk_background.py
python3 -B sk_waterfall.py regress
python3 -B sk_waterfall.py case quad_D1 quad_D0.3 quad_D0.1     # ~1.5-4.5 min per case, one core
python3 -B sk_waterfall.py case quad_D0.03 lin_D0.1 lin_D0.06
python3 -B sk_waterfall.py case hill_D0.1 hill_D0.06 log_D1
python3 -B sk_waterfall.py case log_D0.3 log_D0.1 log_D0.03
python3 -B sk_waterfall.py combine
python3 -B sk_ns.py; python3 -B sk_pbh.py; python3 -B sk_preheat.py
```

## Literature (arXiv IDs and titles verified against the arXiv API, 2 Nov 2026)

- **Topological inflation:** Sakai–Shinkai–Tachizawa–Maeda gr-qc/9506068 (abstract: "η ≳ 0.33 m_Pl").
- **Mild waterfall and PBHs:**
  - Clesse 1006.4522;
  - Kodama–Kohri–Nakayama 1102.5612;
  - Clesse–García-Bellido 1501.07565;
  - Kawasaki–Tada 1512.03515 (abstract: PBHs "rather overproduced" for long waterfalls);
  - PBHs in a matter era: Harada–Yoo–Kohri–Nakao 1609.01588.
- **Preheating:**
  - Greene–Kofman–Linde–Starobinsky hep-ph/9705347;
  - Felder–Kofman–Linde hep-ph/9812289.
- **Higgs:**
  - Elias-Miró–Espinosa–Giudice–Lee 1203.0237;
  - Ema–Mukaida–Nakayama 1602.00483. The brief's 1602.07719 is "Carbon gas in SMC low-metallicity star-forming
    regions", unrelated; the track's correction stands.
- **Species bound:** Dvali 0706.2050.
- **Data:**
  - Planck 2018 X 1807.06211;
  - BICEP/Keck XIII 2110.00483;
  - ACT DR6 2503.14452 (n_s = 0.974 ± 0.003, P-ACT-LB);
  - ACT DR6 2503.14454 (dn_s/dln k = 0.0062 ± 0.0052).
