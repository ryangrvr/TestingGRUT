# Track inf_critic: final verdict of the inflaton-sector pass

3 November 2026 · Claude, inflaton-sector pass (BRIEF_INFLATON.md of this directory; GitHub issue #13 context).

This is same-session AI work on one supplied model, not outside replication. Nothing here is evidence for an inflaton,
the hidden sector, the holder or any dark-matter branch, and nothing is a validation or an exclusion. **Every potential,
slice convention, holder calibration point, decay width, portal, renormalisation condition, ledger and survival
criterion is a model choice; none follows from closure.** Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6) are
retrospective.

Labels: [D] derived, [C] computed (own code in `code/`, output in `results/`, or a reviewed track's computation),
[E] estimate, [Q] quoted from a reviewed track, [A] argued, [R] recalled or literature.

Inputs: the four build/chain tracks and their skeptics, all read in full or in the parts that matter:
- tracks/inf_potential (NOTES, SKEPTIC);
- tracks/inf_reheat (NOTES, SKEPTIC);
- tracks/inf_chain (NOTES, SKEPTIC, results/*.json);
- tracks/inf_higgs (NOTES, SKEPTIC).

Git was read-only and nothing was committed. All runs were in the foreground, one process each (≈ 1 min in total).

## 0. Verdict in brief

1. **No benchmark reproduces M1r's decisive either/or [C].** On every explicit inflaton the declared e = 10⁻⁶ is deep in
   free X (S/Q = 155–987, 45–280× the 3.47 threshold), and no e_crit exists at any e. The relic is a mostly symmetric
   X X̄ gas (p = 8–55 X̄ per net charge, mass per charge M = 17–110 m_X).
2. **Only the quiet T-models with a slow, realisable decay (TB1, α = 0.3; TB2, α = 1) close all four checks
   (CMB, dark matter, reheating, Higgs).** Even they close only conditionally:
   - the Higgs needs an added R-even portal λ_IH ≥ 3×10⁻¹⁰ to survive at all M_t and after inflation;
   - the induced portal |A_IΣ| must be ≲ 10⁻⁶ H_i;
   - the pair-fraction recipe at S/Q ≈ 600 is trusted;
   - ACT DR6 puts their n_s at −3.6σ (retrospective).
   On TB the holder's CP-odd charge sets only 0.9–1.5% of Ω_DM; the hidden entropy sets the rest.
3. **HT6t, T1 (and E1) are not consistent [C].**
   - They miss the dark-matter calibration at their design H_i (ratio 0.21 / 0.076).
   - Their width Γ_I = 0.5 H_i has no realisable channel.
   - HT6t has no φ_I → −φ_I symmetry, so the CP-sign bookkeeping does not apply to it.
4. **Structural reason [C/E, new here].** For single-field α-attractors:
   - a quiet inflaton (m_I ≤ 4 H_i at the slice) needs α ≥ 0.07, so H_i ≥ 4.1×10¹² GeV and r ≥ 3×10⁻⁴;
   - a dark-matter-consistent release with even moderately fast decay (Γ_I ≥ 10⁻³ H_i) needs H_i ≤ 4.5×10¹¹ GeV, and
     with M1r-like fast decay ≈ 1.2–1.7×10¹¹ GeV;
   - the two requirements miss by ≥ ×9 in H_i (≥ ×80 in r).
   M1r's assumed history (m_I = 3, Γ_I = 0.5, H_i ≈ 1.3×10¹¹ GeV) is therefore not realised by any single-field plateau
   tried. The quiet branch is forced into slow reheating, K ≈ 2000, and free X.
5. **Single most important open item:** an inflaton sector, necessarily beyond the single-field plateaus tried, whose
   holder-carrying field is quiet (m ≈ 2–4 H_i) at H_i ≈ 1.2–1.5×10¹¹ GeV and empties into the Standard Model within
   about a Hubble time (K ≈ 2). That would be non-perturbative, since no perturbative channel can do it at an M_P-scale
   amplitude. It must be computed through the chain. It is the only route back to M1r's e-controlled branch point and a
   CP-set relic (§5).

## 1. Contradictions between tracks, and how they settle

| # | Contradiction | Quick check (this track) | Resolution |
|---|---|---|---|
| 1 | inf_potential: "the holder calibration point decides the branch" (S1 → free X 4.6–4.7; ε_H = 1 → Q-matter ≈ 2.0). Skeptic: an omitted homogeneous Y_H remnant (×7116 M1r) puts every case at S/Q ≳ 10² | K1: S/Q re-derived from inf_chain's W for all 7 rows; it reproduces exactly (closed form, rel. 0) | **The skeptic stands** [C]. inf_chain's converged releases give 155 (S1) and 302 (ε_H = 1 holder), both free X. The ε_H = 1 calibration also flips the charge sign (N_X = +3.8×10⁶). The calibration point changes the CP sign and the magnitude, not the branch |
| 2 | Γ_I needed by the dark matter: inf_reheat 0.6–1.1×10⁻⁶ m_I; its skeptic 1.1–3.8×10⁻⁷; inf_chain 1.4–3.4×10⁻⁹ | K3: Γ_cal = (Y_req/(c_Y M))² reproduces inf_chain's 7.512×10⁻⁹ H_i to 10⁻⁴. With M = 1 it would be 3.5×10⁻⁵ H_i (Γ/m = 1.6×10⁻⁵). The rolling release yields 0.25× the static oscillator's charge per Γ^{1/2} | **inf_chain stands** [C]. It is the only track that computes the release, K, vectors and p on one background. The ×10²–10³ shift is the symmetric relic: Γ ∝ M⁻², M² ≈ 4.7×10³. The rolling-release yield pushes back ×14. Earlier values assumed M = 1 |
| 3 | inf_reheat: e_crit ≈ 7–8×10⁻⁵ [E] for B1–B3. inf_chain: no e_crit at any e | inf_chain's Γ scan (read): e_crit = 1.1×10⁻⁴ / 3.2×10⁻⁴ at Γ = 10⁻³ / 10⁻⁴ H_i; none for Γ ≤ 10⁻⁵ | **Both right at their own Γ** [C]. At the dark-matter-consistent Γ = 7.5×10⁻⁹ the matter-era Y_L deposit alone gives S/Q(e → ∞) = 5.1 > 3.47 (margin ×2.4 on R_Y, skeptic) |
| 4 | inf_higgs: "with the S1 plasma, safe to M_t ≥ 174.6". Skeptic: capped at 173.75 (tree) / 174.15 (loop) by inflation; field-grid cut | K7 on TB1 shows the same structure. Without a portal log10⟨N⟩ crosses 0 at 173.10–173.13; with λ_IH = λ_min it is ≤ −3.4×10⁴ to 174.6 | **The skeptic stands** [C]. A plasma cannot act during inflation. Only a portal (λ_IH) lifts the inflationary cap |
| 5 | inf_chain C10: HT6t iteration "diverging". Skeptic: contracting, fixed point near 2×10⁹ | K2: log-mismatch −1.565 → −0.989 (contraction 0.63); secant fixed point 2.04×10⁹ GeV (r = 6.8×10⁻¹¹); T1 1.18×10¹⁰ GeV (r = 2.3×10⁻⁹) | **The skeptic stands** [C on two points]. Neither fixed point is near the design 1.29×10¹¹, and both need the blocked Γ_I = 0.5 |
| 6 | inf_chain: "Q-matter needs slow-roll pre-slice holding". Skeptic: wrong-signed | WORKING_THEORY §4: free X at 10⁻⁶ needs f_eff ≥ 0.47; rolling lowers f_eff | **The skeptic stands** [D]. Rolling helps Q-matter. The missing ingredient is fast reheating with a quiet inflaton at the dark-matter scale (§4) |
| 7 | inf_chain item 6: "the yield moves by K^{3/2}" under the width-switch convention. Skeptic: a re-normalisation | — (skeptic's engine re-run) | **The skeptic stands** [C]: yield +0.5–0.9%, S/Q −1% |
| 8 | inf_potential (late Y_L ledger) vs inf_reheat (prompt radiation) for the φ_I → YY heat: ×10³ apart | inf_chain rows without R_YY: S/Q = 150 / 121 / 615 / 987 | **Immaterial to every verdict** [C]. dq_heat's late ledger is used (M1r's own). It changes H_i,req of HT6t by 2%, of T1 by ×3, and TB not at all |
| 9 | λ_IH window (inf_reheat) vs λ_min (inf_higgs): factor-2 convention | — | Convention m_H² = λ_IH φ_I²/2 (inf_higgs) adopted. TB: λ_min = (20/3)(m_I/M_P)² = 3.1/3.0×10⁻¹⁰, inside inf_reheat's q ≤ 1 window [1.5×10⁻¹⁰, 1.5–2×10⁻⁹] [E] |
| 10 | inf_potential exit factor (half-applied K shift) vs skeptic | — | Superseded by inf_chain's own vectors on own backgrounds (f_eff 0.909 / 0.866, matching the uncorrected numbers; the K effect sits in the deposit) [C] |

## 2. Per-benchmark verdict (e = 10⁻⁶, portals vanishing at M_Y, a = 1 and H = H_i at S1)

Sources: inf_chain results/summary.json [C there]; inf_higgs and both skeptics [Q]; this track's settle_checks.json
(K4–K9) and quiet_scan.json. Full table: `results/verdict.json`.

**Inflation and reheating**

| | HT6t | T1 (E1) | TB1 | TB2 |
|---|---|---|---|---|
| potential | V0(1 − βx − x⁶)²/f_top², v = 0.34694 M_P, β = 1.5529×10⁻⁵ (not reflection-even) | V0 tanh²(φ/√6α), α = 6.0968×10⁻⁵ (E1: E-model, same α) | T-model α = 0.3 | T-model α = 1 |
| m_I/H_i | 42.4 | 128 (256) | 2.21 | 1.54 |
| n_s; α_s | 0.9649; −9.5×10⁻⁴ | 0.9618; −7.3×10⁻⁴ | 0.9619; −7.2×10⁻⁴ | 0.9622; −7.2×10⁻⁴ |
| r | 2.71×10⁻⁷ | 2.71×10⁻⁷ | 1.32×10⁻³ | 4.22×10⁻³ |
| H_* / H_i (GeV) | 1.290×10¹¹ / 1.29×10¹¹ | 1.290×10¹¹ / 1.29×10¹¹ | 8.99×10¹² / 7.47×10¹² | 1.61×10¹³ / 1.06×10¹³ |
| reheating channel | Γ_I = 0.5 H_i engine convention; **no realisable channel** (y φ_I N̄N blocked, y\|φ\|/(m/2) ≈ 1.8×10⁴) | same; blocked (3×10³; E1 2×10²); fragments at 0.35–0.43 H_i⁻¹ | Majorana −(y/2)φ_I N̄N, y = 4.1×10⁻⁴ (inside y_nat 4–9×10⁻⁴), or axial S2 with f ≈ 2.4–9×10¹⁵ GeV (M_N = 1–4×10¹²); Γ_I = 7.5×10⁻⁹ H_i | y = 2.6×10⁻⁴ or S2 f ≈ 0.4–1.4×10¹⁶ GeV; Γ_I = 2.1×10⁻⁹ H_i |
| T_RH (GeV) | 2.1×10¹⁴ (nominal only) | 2.1×10¹⁴ (nominal only) | 2.0×10¹¹ (T_max 1.3×10¹³) | 1.3×10¹¹ (T_max 1.1×10¹³) |
| K = a²H/H_i, a = 1 at S1 (÷ e² for a = 1 at ε_H = 1; M1r: 1.99 at its release) | 12.3 (1.67) | 15.2 (2.06) | 2000 (271); 12.2 matter e-folds | 2589 (350); 13.0 |
| ACT DR6 n_s pull [R] | −2.8σ | −3.7σ | −3.7σ | −3.6σ |

**Chain, relic, calibration, Higgs, by-products**

| | HT6t | T1 (E1) | TB1 | TB2 |
|---|---|---|---|---|
| S/Q at 10⁻⁶ (e → ∞) | 155 (61) | 646 (337); E1 ≈ 285 (148) [Q skeptic] | 615 (5.1) | 987 (7.9) |
| branch; e_crit | free X; none | free X; none | free X; none (×2.4 margin on R_Y) | free X; none (×5.5) |
| dominant heat | Y_H remnant cascade (R = 627) | Born φ_I → Y_LY_L (8.3×10³) | matter-era vectors (1.8×10³) | vectors (3.2×10³) |
| relic: p; M per charge | 8.1; 17 m_X | 36; 72 m_X | 34; 69 m_X (6.3 GeV) | 55; 110 m_X |
| asymmetric share of Ω_DM | 5.8% | 1.4% | 1.5% | 0.9% |
| H_i required / potential | 2.7×10¹⁰ / 1.29×10¹¹ = **0.21** | 9.9×10⁹ / 1.29×10¹¹ = **0.076** | **1.000** (Γ_I absorbs it) | **1.000** |
| fixed point of iteration | ≈ 2.0×10⁹ GeV (r ≈ 7×10⁻¹¹) [C on 2 points] | ≈ 1.2×10¹⁰ (T1_cal1, ratio 1.15; m_I = 1625 H_i; fragments) | — | — |
| Higgs, inflation (ξ_H = 0.1) | safe; fails ≥ 173.75 (tree) / 174.15 (loop) | same | M_t,crit ≈ 173.10–173.13 without a portal [C, K7] | ≈ 173.1 [Q skeptic] |
| Higgs, after inflation, gravitational only | M_t,crit 172.92 (MD) / 172.875 (Γ = 0.5) at 1.29×10¹¹; ≈ 173.4 at its consistent 2.7×10¹⁰ [E] | 172.87–172.92; ≈ 173.6 at 10¹⁰ [E] | ≈ 172.1 [E, K8]: **fails at central M_t** | ≈ 172.1 [E] |
| Higgs with λ_IH ≥ λ_min | safe (λ_min 3.4×10⁻¹¹) [Q] | safe (3.1×10⁻¹⁰; E1 1.2×10⁻⁹) [Q] | inflation log10⟨N⟩ ≤ −3.4×10⁴ up to 174.6 [C]; after inflation transferred [E]. Needs λ_IH ≳ 2.5 λ_min if S1's R-even loop mass is renormalised at M_N | safe at λ_min [E] |
| leptogenesis by-product | not defined (no channel) | not defined | non-thermal φ_I → NN: ε₁ = 1.35×10⁻⁸, i.e. 0.9–4.5×10⁻⁵ of the Davidson–Ibarra bound for M_N = 2–7×10¹²; M_N = 10¹² is washed out [E] | ε₁ = 2.1×10⁻⁸; M_N ≥ 10¹² works [E] |
| remaining open inputs | no channel; not R-even; not DM-consistent | no channel; fragmentation lattice; not DM-consistent | p(S/Q ≈ 600); \|A_IΣ\| ≲ 10⁻⁶ H_i; λ_IH (an added SM portal); Higgs and plasma at H ≈ 7×10¹² not run; the large-e margin; ACT | as TB1 |

Notes on the table:
- HT6t_cal1, HT6t_Send and T1_cal1 are in `results/verdict.json`.
- HT6t_Send: S/Q = 302, opposite charge sign, H_i,req = 1.3×10⁹ GeV.
- The leptogenesis estimate [E] uses Y_L = (3/2)ε₁ T_RH/m_I and Y_B = 8.7×10⁻¹¹, with inverse-decay washout in the
  Buchmüller–Di Bari–Plümacher approximation (hep-ph/0401240 [R]) at K_w = 46, from z = M_N/T_RH.
  - Thermal N are produced at T_max > M_N in the matter era, but the asymmetry that survives comes from the late
    decays near T_RH [A].
  - The needed CP phase is supplied and unrelated to α_I.
- **S1 versus S2 at TB [E, K6].**
  - S1's induced A-term is 140 / 90 TeV, so Aφ_*/H_*² = 0.017 / 0.004. Near zero crossings the most negative m² it gives
    against λ_min is −6×10⁻⁷ H_i², which is harmless.
  - S1's R-even N–L–H loop mass is −18 / −3.3 H_*² per unit of [ln(M²/μ²) − 1]. With μ = M_N that is ≈ −230 / −41 H_*²,
    against the λ_min portal's +190 / +86 H_*². TB1 therefore needs λ_IH ≳ 2.5 λ_min (still q < 1) or the axial S2.
  - S2 has no φ-dependent |m_N|, so it carries neither term.

## 3. Which benchmarks give a consistent M1r, and what that choice fixes

**Answer.** TB1, and TB2 alike, are the only benchmarks that close CMB + dark matter + reheating + Higgs. They do so with
four added or tuned inputs:
1. a realisable perturbative Majorana channel (y ≈ 3–4×10⁻⁴; or S2 axial, f ~ 10¹⁶ GeV);
2. an R-even portal λ_IH ∈ [3×10⁻¹⁰, ~2×10⁻⁹], needed for the Higgs after inflation at H ≈ 10¹³ and above M_t ≈ 173.1;
3. a tuned induced portal |A_IΣ| ≲ 10⁻⁶ H_i;
4. trust in the case-L pair-fraction recipe at S/Q ≈ 600–1000.

TB1 is the primary. Its m_I = 2.21 H_i lies inside the quiet window and its r is lower. TB2 has the larger
no-e_crit margin.

This is a consistent embedding of the M1r **Lagrangian**, not of the M1r **picture**:
- there is no branch point in e;
- the relic is free X at every e;
- the CP-odd holder sets ≈ 1.5% of Ω_DM;
- the abundance is set by the hidden entropy (n_X̄/s_h ≈ 0.05, skeptic [C]) and the reheating dilution.

**Predictions fixed by choosing TB1 (TB2)** [C on the model choices]:
- r = 1.32×10⁻³ (4.22×10⁻³). Single-field consistency gives n_t = −r/8. For the quiet α-attractor family as a whole,
  r ≥ 3×10⁻⁴ for m_I ≤ 4 H_i and r ≥ 6×10⁻⁴ for m_I ≤ 3 H_i [C, quiet_scan]. This would be within the stated reach of
  planned CMB B-mode experiments (σ_r ~ 10⁻³) [R].
- n_s = 0.9619 (0.9622) and running α_s = −7.2×10⁻⁴ (−7.2×10⁻⁴). The long matter era (12–13 e-folds) is already in N_*.
- Relic: free X at every e (no Q-matter), mostly symmetric:
  - p = 34 (55) X̄ per net charge, i.e. 49.3% (49.5%) X̄ by number;
  - mass per net charge 69 (110) m_X;
  - cold, with no annihilation at rest (no lighter hidden states) [Q WORKING_THEORY §3].
- **No objects**: the Q-matter object mass is not applicable.
- T_RH ≈ 2×10¹¹ (1.3×10¹¹) GeV and H_i = 7.5×10¹² (1.06×10¹³) GeV. T_RH inherits p: Γ ∝ M⁻² and T_RH ∝ M⁻¹ at fixed
  H_i [D].

HT6t, T1/E1, the iterates and HT6t_Send are not consistent: dark matter, reheating channel and, for HT6t, the CP
bookkeeping all fail.

## 4. The quiet–scale–reheating triangle (`code/quiet_scan.py` → `results/quiet_scan.json`)

**Method.**
- T- and E-model α-attractors on α = 10⁻⁴–3, with V0 fixed by A_s at N_* = 52.5 (exact Mukhanov–Sasaki, inf_potential's
  code, read-only) and the slice at S1 [C].
- Regression on TB1: H_i = 7.444×10¹² (chain 7.466×10¹² at N_* = 52.34), m_I = 2.2077 H_i (identical), r = 1.307×10⁻³
  (1.315×10⁻³).

**Results.**

| T-model α | 10⁻⁴ | 10⁻³ | 10⁻² | 0.03 | 0.1 | 0.3 | 1 | 3 |
|---|---|---|---|---|---|---|---|---|
| H_i (GeV) | 1.66×10¹¹ | 5.2×10¹¹ | 1.64×10¹² | 2.8×10¹² | 4.8×10¹² | 7.4×10¹² | 1.06×10¹³ | 1.28×10¹³ |
| m_I/H_i | 100 | 32 | 10 | 5.9 | 3.4 | 2.2 | 1.5 | 1.3 |
| r | 4.5×10⁻⁷ | 4.5×10⁻⁶ | 4.5×10⁻⁵ | 1.3×10⁻⁴ | 4.4×10⁻⁴ | 1.3×10⁻³ | 4.2×10⁻³ | 1.2×10⁻² |

- The E-model is twice as fast to oscillate at the same α: it reaches m_I = 4 only at α ≈ 0.35.
- **Quiet threshold** (T-model, at S1): m_I = 4 / 3 / 2 H_i at α = 0.071 / 0.14 / 0.40, with H_i = 4.1 / 5.5 /
  8.2×10¹² GeV and r = 3.1 / 6.1 / 17×10⁻⁴.
- **Dark-matter requirement of a quiet release** [C from inf_chain's scans; E as a family curve].
  - H_i,req(Γ_I) from TB1 and TB2 agree to 3.7–4.0% over Γ = 10⁻¹⁰–10⁻³ H_i, so it is nearly universal.
  - It is 4.5×10¹¹ GeV at Γ = 10⁻³ and rises as Γ^{−0.21}.
  - Extrapolated to Γ = 0.1 / 0.5 it gives 1.7 / 1.2×10¹¹ GeV, consistent with M1r's own (3, 0.5) value of 1.29×10¹¹ [E].
- **Closure for quiet α** (α = 0.1–3): Γ_cal = 1×10⁻⁹–5×10⁻⁸ H_i, K ≈ 1.0–3.8×10³, S/Q(10⁻⁶) ≈ 480–780, always free X [E].
- **Closure with fast decay.** Γ ≥ 10⁻³ needs α ≤ 7×10⁻⁴, i.e. m_I ≥ 37 H_i. That is a noisy release (Y_H remnant,
  φ_I → Y_LY_L), whose H_i,req lies far lower still (T1: 9.9×10⁹).
- **The gap**: ≥ ×9 in H_i and ≥ ×80 in r even at the computed Γ = 10⁻³; ≥ ×24 / ×600 at Γ = 0.1 [C/E].

**Scope [A].** This covers α-attractor plateaus, and the potential track's hilltops and polynomial plateau (all with
m_I ≥ 13 H_i at H_i ≈ 1.3×10¹¹ [Q]).
- **Not covered:** a two-scale single-field potential, i.e. a flat, inflection-like CMB plateau joined to an M_P-wide
  basin (m/H ~ √6 M_P/d for basin width d), and multi-field endings.
- **Even with the shape solved,** fast reheating at an M_P-scale amplitude is non-perturbative. With y|φ_I| < m_I/2
  required at Φ ~ M_P and m_I ~ 4×10¹¹ GeV, y < 2×10⁻⁷, so Γ/m ≲ 10⁻¹⁵ [E].

## 5. The single most important remaining open item

**An inflaton sector that gives the holder field M1r's assumed post-inflation history, computed through the chain.**
It needs:
- a quiet oscillation, m ≈ 2–4 H_i at the slice;
- emptying into the Standard Model within about a Hubble time (K ≈ 2);
- H_i ≈ 1.2–1.5×10¹¹ GeV, with A_s, n_s and r ≈ 3×10⁻⁷.

**What it would take.** A non-perturbative transfer is unavoidable, for example instant preheating into N or H at the
zero crossings. Natural candidates [A]:
- a two-scale single-field potential;
- a hybrid/waterfall ending, with the holder on the waterfall field;
- a light spectator.

**Why it matters.** It is the only remaining route to:
- M1r's e-controlled branch point (free X vs Q-matter at e ≈ 10⁻⁶);
- an asymmetric relic set by α_I and sign(φ_I);
- the Planck-mass objects of the Q-matter branch.

Without it, the consistent embedding is the TB type: free X at every e, mostly symmetric, with r ~ 10⁻³.

**Runner-up open items (ranked).**
1. The pair fraction p at S/Q ~ 10²–10³, checked by a Boltzmann solution beyond the case-L recipe. It sets M, and with
   it Γ_cal ∝ M⁻², T_RH ∝ M⁻¹ and the asymmetric share, on the TB line. A ×2 in p moves H_i,req by ×0.63–1.6.
2. The induced A_IΣ renormalisation condition. It is dominant for slow reheating (R_leak 6×10²–9×10⁵), and must be
   tuned to ≲ 10⁻⁶ H_i or Γ_cal re-solved with it.
3. Higgs survival after inflation at H ≈ 10¹³ GeV with λ_IH and the S1/S2 plasma. It has only been transferred here.
4. ACT DR6's n_s. All T-model benchmarks sit at −3.6σ. A quiet model with n_s ≈ 0.974 would need a different shape [R/A].

## 6. Recommendations to the combiner (WORKING_THEORY)

- **§0, §4, §9.** Replace "the branch depends on the inflaton history" with:
  - on every explicit inflaton tried, e = 10⁻⁶ is deep in free X with a mostly symmetric relic;
  - the decisive either/or survives only for an inflaton sector that is quiet, fast-reheating and at the dark-matter
    scale at once, and none has been constructed.
- **§2 inputs.**
  - (m_I, Γ_I) = (3, 0.5) is not realisable as a single-field plateau with a perturbative channel.
  - The consistent TB embedding adds λ_IH (an R-even inflaton–Higgs portal, so "zero SM portals" no longer holds) and a
    tuned A_IΣ.
- **§5.** H_i on the TB line is fixed by A_s, not fitted. The fit moves to Γ_I (equivalently T_RH). The relic is
  mostly symmetric, M = 69–110 m_X per charge.
- **§7, CP sign.**
  - It applies only to reflection-even potentials (T-models), not to HT6t.
  - Where the holder is calibrated can flip the sign (HT6t_Send).
  - On TB the CP charge controls only ≈ 1.5% of Ω_DM.
- **§8, ranked open items.** Put the triangle of §4 above first, then p(S/Q ≫ 1), then A_IΣ.

## 7. Limits

- Every number inherits its track's model choices.
  - Releases: homogeneous and classical.
  - Vectors: the single-Higgs-direction Proca model.
  - Ledger: dq_heat's.
  - p: the conditional case-L recipe.
  - Higgs: the argued D criterion and the Mantziris–Markkanen–Rajantie bubble count.
- **This track's estimates [E]:**
  - the quiet-release H_i,req curve applied as a family curve and extrapolated to Γ > 10⁻³;
  - the gravitational-only M_t,crit extrapolated in ln H_i from the inf_higgs skeptic's table (made on the HT6t
    background);
  - the post-inflation portal safety at TB, transferred from HT6t/T1;
  - leptogenesis (non-thermal, BDP washout, no flavour effects or ΔL = 2 scattering);
  - S1 loop sizes at φ_* ≈ 4–5 M_P.
- K7 is tree-level, uses the inf_higgs criterion with my TB1 history, and reproduces the inf_chain skeptic's
  log10⟨N⟩ = −35.26 at 172.57 (mine −35.27).
- Not computed here: any release, vector, Higgs mode or kinetic run; a lattice; the two-scale or multi-field
  constructions of §5.
- **Process note.** The inf_chain skeptic's stray copy `/sk_gamma_switch.json` (2.7 KB, filesystem root) still exists.
  The owner should delete it; it is outside the writable tree, so this track did not touch it.

## Code and results

| file | purpose |
|---|---|
| code/quiet_scan.py | α-attractor quiet/scale scan and closure with inf_chain's Γ scans → results/quiet_scan.json |
| code/settle_checks.py | checks K1–K9 (closed form, iteration secant, Γ decomposition, T_RH/T_max, leptogenesis, portal/S1/S2, TB1 Higgs inflation with λ_IH, M_t,crit extrapolation, n_s pulls) → results/settle_checks.json |
| code/verdict.py | assembles the per-benchmark table → results/verdict.json |

Reproduce (one core, about 1 minute):

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_critic/code; export OMP_NUM_THREADS=1
python3 -B quiet_scan.py; python3 -B settle_checks.py; python3 -B verdict.py
```

All imports of other tracks are read-only:
- inf_chain ch_common, which pulls in inf_potential inf_common/calibrate and dq_heat's dq_model_copy;
- inf_higgs ih_inflation/ih_common, with higgs_reheat smpot.

## Credits

- α-attractors: Kallosh–Linde, arXiv:1306.5220 [R].
- Bubble-count criterion: Mantziris–Markkanen–Rajantie, arXiv:2011.03763 [R].
- Leptogenesis washout: Buchmüller–Di Bari–Plümacher, hep-ph/0401240 [R]. Davidson–Ibarra bound: hep-ph/0202239 [R].
  Non-thermal leptogenesis from inflaton decay: Lazarides–Shafi, Phys. Lett. B 258 (1991) 305 [R].
- Instant preheating: Felder–Kofman–Linde, hep-ph/9812289 [R].
- Data, retrospective: Planck 2018, arXiv:1807.06211; BICEP/Keck, arXiv:2110.00483; ACT DR6, arXiv:2503.14452 (via the
  inf_potential skeptic) [R].
- Our tracks, read-only: inf_potential, inf_reheat, inf_chain, inf_higgs (and skeptics), dq_heat, dq_entropy,
  higgs_reheat.
