# Track mf_critic: final verdict of the multi-field inflaton pass

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md of this directory; GitHub issue #13 context).

This is same-session AI work on one supplied model, not outside replication. Nothing here is evidence for a hybrid
inflaton, a spectator, the hidden sector or any dark-matter branch. Nothing is a validation, an exclusion of the real
world, a new theory or a GRUT prediction. **Every architecture, coupling, reheating closure, holder amplitude and phase,
rescaling rule, renormalisation condition, ledger, pair-fraction map and survival criterion is a model choice; none
follows from closure.** Data comparisons (Planck 2018, Planck 2015 XX, BICEP/Keck 2021, ACT DR6, Σm_ν) are
retrospective.

Labels: [D] derived, [C] computed (this track's code in `code/`, output in `results/`, or a reviewed track's
computation), [E] estimate, [Q] quoted from a reviewed track, [A] argued, [R] literature (arXiv numbers checked by web
search in this session or by the quoted track).

**Inputs, read in full or in the parts that matter:** BRIEF_MULTIFIELD.md (all estimates G1–G7), WORKING_THEORY.md §0,
§4a, §4b, §7–§9, BRIEF_INFLATON.md, BRIEF_DECISIVE.md; inf_critic NOTES; inf_chain, inf_reheat NOTES and SKEPTIC (via
inf_critic and the multi-field tracks); mf_gemini, mf_hybrid, mf_chain, mf_higgs NOTES and SKEPTIC; mf_ideal NOTES and
code README; the skeptics' result files (`mf_chain/skeptic/results/sk_summary.json`,
`mf_gemini/skeptic/results/sk_ns.json`).

**Process.** Git read-only, nothing committed. All runs in the foreground, sequential, at most two worker processes
(mf_ideal's vector stage uses two). About 25 min of engine time in total. Engine binaries, tables and run series are in
scratch (`$SP/wt2/mf_critic/`). Python ran with `-B`, so there is no `__pycache__` or binary in the tree. Other tracks'
code (mf_ideal's engine and drivers, mf_chain's `hist_build.py`) was **run read-only, with every output redirected to
this track's scratch**.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict in brief

1. **The external (Gemini) proposal fails as written [D/C, mf_gemini + skeptic, re-derived here].**
   - Its own two conditions, m_χ = 3 H_i and λv⁴/4 = V0, force v = 1.633 M_P and λ = 4.81×10⁻¹⁵. That makes the
     waterfall field a quiet, Planckian holder of V0: the single-field obstruction (G1) moved onto χ.
   - The Yukawa it relies on is kinematically blocked ×4.1×10⁷, and its loop exceeds the tree potential ×2.1×10¹³.
   - The waterfall is mild (≥ 17 e-folds), giving PBH overclosure and walls.
   - No channel tried reaches Γ ≥ 0.5 H.
   - What survives: the low scale (r ≈ 3×10⁻⁷) and a red tilt with a non-convex V_SR. Everything else is unsupported or
     refuted (§2).
2. **The repaired hybrid (mf_hybrid HFa/HFb) passes the inflationary side; its post-inflation side is not settled.**
   - **Passes [C].** CMB (n_s 0.9649 / 0.9743, r = 2.75×10⁻⁷, α_s ≈ +1×10⁻⁴); fast, biased waterfall; no PBH; a
     realisable SM channel (axial couplings to seesaw singlets); Higgs safe given the medium masses [E].
   - **K [C on E].** With Pauli blocking, K is 4–10, not 2.
   - **Walls are not excluded.** The waterfall fragments at t ≈ 0.45 H_i⁻¹, and no lattice has been run.
3. **At the inherited holder calibration (μ_eff(slice) = 10i, "c0 = 1"), every hybrid history is free X at e = 10⁻⁶
   and 3×10⁻⁶, with no e_crit [C, mf_chain, reproduced by its skeptic].**
   - S/Q = 14.6–14.7 (exported) and 27–72 (Pauli-blocked).
   - The self-consistent H_i is ≈ 1.11×10¹¹ GeV (unblocked) or 2.40×10¹¹ GeV (Pauli).
4. **New here [C]: the holder amplitude decides the branch on M1r's own supplied history too.**
   - On M1r's (3, 0.5) static slice, S/Q(10⁻⁶) = 6.07 / 4.955 / 4.25 / 3.51 / 2.12 at c0 = 0.9 / 1 / 1.05 / 1.1 / 1.2.
   - e_crit goes 2.86 → 2.06 → 1.54 → 1.02 → 0.32×10⁻⁶, i.e. e_crit ∝ c0^−4.9 near c0 = 1.
   - **A holder 10.2% stronger than the inherited calibration puts M1r itself at the Q-matter threshold at e = 10⁻⁶.**
   - So the decisive pass's "free X at 10⁻⁶, e_crit = 2.06×10⁻⁶" is a c0 = 1 statement. M1r's branch is a surface in
     (e, c0), and c0 is supplied (it descends from Codex's Option 1A).
5. **M1r's picture is recovered only by re-tuning the holder, and where it lands depends on an uncomputed rate.**
   - **Pauli Γ_th = 0 [C, skeptic; c0 = 1.40 reproduced here to all digits on a rebuilt history].** At c0 ≈ 1.36 and
     H_i ≈ 8.4–8.9×10¹⁰ GeV: free X at 10⁻⁶ (4.16), e_crit = 2.2×10⁻⁶, an asymmetric relic, r ≈ 1.2×10⁻⁷.
   - **Q-matter at 10⁻⁶ on the same history.** At c0 = 1.40 and H_i = 8.40×10¹⁰ GeV, self-consistent to 0.2%:
     S/Q = 2.88, e_crit = 0.60×10⁻⁶, M = 0.807 m_X.
   - **Other end of the thermalisation bracket, Γ_th = 10 H_i (new here) [C].** The same c0 values miss the
     self-consistent H_i by ×0.45–0.54. The M1r-like point moves to c0 ≈ 1.25 at H_i ≈ 5.2×10¹⁰ (S/Q 3.64,
     e_crit 1.3×10⁻⁶).
   - **Q-matter at 10⁻⁶ on Γ_th = 10.** It needs c0 ≈ 1.27–1.30 at H_i ≈ 4.3–5×10¹⁰. There the axial EFT ratio is
     1.05–1.5 under mf_chain's rescaling rule, at or beyond its cutoff.
   - **The tuning.** The window is ±2–3% in |μ_eff| and ±5° in phase. On M1r's own history ±10% suffices.
6. **Single most important open item (§6).** A computed post-waterfall history: a 3D lattice of the biased waterfall
   with the axial singlet channel, its Pauli blocking and the N thermalisation rate. It decides three things:
   - whether walls form;
   - K, and hence H_i;
   - where (and whether inside the EFT) a self-consistent branch point exists.

   Even then the branch at the declared e remains a supplied choice of c0.

---------------------------------------------------------------------------------------------------------------------

## 1. Contradictions between tracks, and how they settle

Quick checks K1–K11 are in `code/settle_checks.py` (→ `results/settle_checks.json`), `code/m1r_c0.py`
(→ `results/m1r_c0.json`) and `code/hyb_c0.py` (→ `results/hyb_c0.json`).

| # | Contradiction | Quick check (this track) | Resolution |
|---|---|---|---|
| 1 | mf_gemini: log-family V_SR gives n_s = 0.92–0.97, redder for slower crossings. Skeptic: 0.972–0.998, bluer | **K3 [D]:** pure log V0(1 + A ln σ) with σ_c finite gives n_s − 1 = −1/(N_σ + 1/2Δ). That is 0.9723 / 0.9727 / 0.9736 / 0.9768 at Δ = 1 / 0.3 / 0.1 / 0.03, against the track's 0.9719 / 0.9714 / 0.9695 / 0.9621 and the skeptic's regulated log 0.9727 / 0.9726 / 0.9742 / 0.9850 | **The skeptic stands on the sign.** The magnitude at small Δ depends on the shape (0.977 pure, 0.985 regulated). Either way n_s ≥ 0.972, and it lies inside 0.967–0.973 only for Δ ≳ 0.2–0.3 (N_wf ≲ 19). The track's 1 − 1/N_σ drops the 1/(2Δ) term |
| 2 | mf_gemini: Higgs-portal failure is "boundedness for λ_SM < 0" (M_t-dependent). Skeptic: Landau pole, sign-independent | **K1 [D]:** λ_χH²/λ_χ = (1.004×10⁻⁷)²/4.81×10⁻¹⁵ = 2.10 for either sign | **The skeptic stands.** The opposite-sign route fails on perturbativity whatever M_t is |
| 3 | mf_gemini: Γ ≤ 7×10⁻⁹ H at the loop cap. Skeptic: 1.8×10⁻⁸ | — | **The skeptic stands** (mixed caps). Immaterial |
| 4 | mf_gemini bosonic preheating K ≈ 2–3. Skeptic: K ≥ 2.2, often never, once the decay is consistent | — | **The skeptic stands** [C on a model]. Both are above K_crit = 1.38–1.41 |
| 5 | mf_hybrid S1: "m_σ ≲ 13 H_i is forced". Skeptic: κ-tuning statement | — | **The skeptic stands**: the track's own (4000, 30) row passes. Qualitative point (σ holds V0 if m_χ ≫ m_σ) stands |
| 6 | mf_hybrid K = 2.05. Skeptic K = 4–7 (≈ 22 tails blocked). mf_chain K = 4.27–9.76; "tails blocked" diverges | My rebuild with mf_chain's hist_build at 8.40×10¹⁰ GeV gives K = 6.499 (Γ_th = 0), identical to the skeptic's; Γ_th = 10 gives 3.19 (8.4×10¹⁰) and 2.49 (5.0×10¹⁰) | **mf_chain and the skeptic stand** [C on the skeptic's E fluid]. K ≈ 2 holds only without Pauli blocking. Realistic K = 2.5–10, depending on the uncomputed Γ_th and on H_i |
| 7 | mf_hybrid: "no walls" (ratio ≥ 33). Skeptic: field fragments at t ≈ 0.45, walls not excluded | — (no lattice allowed here) | **The skeptic stands as "not excluded".** If walls percolate they dominate at t ~ 23 H_i⁻¹, before the bias removes them (t ~ 150). **Open** |
| 8 | Spectator energy fraction at decay: 2.4% (mf_hybrid), 1.3% (skeptic, mf_chain) | — | **1.3%** (unblocked), 0.1–0.2% (Pauli) [Q, two independent derivations agree] |
| 9 | β_iso: mf_hybrid 3.8×10⁻³ p² with p = 4.4–5.1 (0.075–0.098); mf_chain D = 1.36 (6.9×10⁻³); skeptic D = 2.49 (0.022), 0.0385–0.047 at the 3×10⁻⁶ fixed point, 0.12 at the strong holder | **K4 [D on Q]:** β = x/(1 + x), x = (H/2πs_*)²D²/A_s reproduces every quoted value (0.0068, 0.0384, 0.0475, 0.12, 0.071, 0.090; 0.017 for D = 2.49 at 1.11×10¹¹) | **All consistent; D is the variable.** A mostly symmetric relic cancels the charge response (D ≈ 0–2.5). An asymmetric relic (free X with M ≈ 1, or Q-matter) has D ≈ 4.5–9.2 and needs s_* ≳ 1.0–1.4×10¹⁶ GeV, a choice that is always available at fixed \|g_h\|s_* |
| 10 | mf_chain: unblocked fixed point 1.28×10¹¹ GeV, M = 1.58. Skeptic: 1.11×10¹¹, M = 1.77 (gap in the p map) | — | **The skeptic stands** [C, network run inside the gap]. Free X either way |
| 11 | mf_chain: "the Y_L remnant forbids Q-matter even at e → ∞". Skeptic: only at c0 = 1 | **K9 regression [C]:** my rebuilt P0 history at 8.40×10¹⁰, c0 = 1.40 gives S/Q(10⁻⁶) = 2.8754, H_req = 8.3826×10¹⁰, e_crit = 5.99×10⁻⁷, identical to the skeptic. **K10 [C]:** the same c0 freedom moves M1r's own branch (§4a). **K11 [C]:** the self-consistent strong-holder point moves with Γ_th (§4b) | **The skeptic stands**, and the reading is sharpened: c0 is a supplied input that decides the branch on every history, M1r's included. "No Q-matter" and "Q-matter at 10⁻⁶" are both c0 choices |
| 12 | mf_ideal: Q-matter at 10⁻⁶ only in a K = 1, c0 = 1.03–1.05 pocket; "never for fast turn-on or K ≥ 1.6". Skeptic: a second region at c0 ≈ 1.3–1.5 | (K9/K11 add Γ_th = 10 points at c0 = 1.25–1.40) | **The skeptic stands.** The c0 at which a branch point exists depends on the history: ≈ 1.0–1.1 (M1r), 1.03–1.05 (mf_ideal K = 1 corner), ≈ 1.25 (hybrid, Γ_th = 10), ≈ 1.36–1.40 (hybrid, Γ_th = 0) |
| 13 | mf_higgs: pessimistic instantaneous bracket caps M_t at 173.95–174.12 GeV. Skeptic: 173.69 at the exact end, and not a bracket | — | **The skeptic stands** [C]. The constant-H HM edge (ξ_min 0.049 / 0.053, cap > 177 GeV extrapolated) is the inflationary result |
| 14 | mf_higgs gravitational-only M_t,crit = 172.88–172.92 at 1.3×10¹¹. Skeptic: 172.74 at the self-consistent Pauli 2.40×10¹¹ (fails at +1σ) | **K6 [E]:** interpolation in ln H_i of mf_higgs's H7 table: 173.01–173.05 at 8.4×10¹⁰, 173.16–173.21 at 5.2×10¹⁰ | **Both stand**; the gravitational-only margin is history-dependent (172.7–173.2). With the singlet-neutrino gas and the SM hard-gas medium masses, every history is safe to ≥ 174.6 GeV [Q, C on E] |
| 15 | mf_hybrid's Σm_ν and anomalous-coupling issues (skeptic N5, N6) | — | **Unresolved; not addressed by any track.** Four singlets at m̃ = 0.05 eV give Σm_ν = 0.20 eV for a real Casas–Ibarra R, in retrospective tension with the cosmological bounds |

---------------------------------------------------------------------------------------------------------------------

## 2. The external proposal as written: which claims survive, and why (plain terms)

The full original text is not in the repository or on #13. This audits the brief's summary of it.

| Claim | Verdict | Why, in plain terms |
|---|---|---|
| A two-field hybrid at H_i = 1.25–1.55×10¹¹ GeV with m_χ = 3 H_i and λv⁴/4 = 3M_P²H_i² | **fails** [D] | The two conditions together fix v = √(8/3) M_P = 1.63 M_P and λ = 4.8×10⁻¹⁵. The field that must later oscillate quietly also carries the whole vacuum energy, so its swing is Planckian. That is the single-field obstruction in a new place |
| χ reheats fast through y_N ≈ 1.5–2 and a Higgs portal, Γ_eff/H = 0.5–1 | **fails** [D/C] | With y = 2 the decay products would weigh 4×10⁷ times what χ can give them, and the loop of those fermions outweighs χ's own potential ×2×10¹³. Keeping the loop below 1% caps y at 3–5.5×10⁻⁴, giving Γ ≤ 2×10⁻⁸ H. The portal needed (λ_χH ≈ 10⁻⁷) makes the Higgs tachyonic during inflation with the proposal's sign; with the opposite sign it runs into a Landau pole near 9×10¹² GeV. Preheating gives K ≥ 2.2 or never. Even the claimed Γ would need the rest of the history unchanged to give Q-matter, and none tried does |
| V_SR(σ) chosen to give n_s = 0.967–0.973 | **survives only partly** [C] | A convex V_SR gives a blue tilt. A hilltop works with \|η\| = 0.0136–0.0166. A log term gives n_s ≥ 0.972, inside the window only if the waterfall is short (§1 #1). "Removing the ACT tension" is retrospective, and HFb shows a hybrid can sit at 0.974 by choice |
| No defects (implicit) | **fails** for the real χ of the proposal [C] | The quantum-seeded waterfall lasts ≥ 17 e-folds, the δN tail overcloses with PBHs (f_PBH ≥ 10¹⁰), and walls form at 0.99 of the topological-inflation threshold (Sakai et al., gr-qc/9506068 [R]). A bias of ≥ 5 H0 fixes the defects only |
| S/Q < 3.47 | **unsupported as written**; reachable only as a tuned holder choice [C] | No consistent history exists for the proposal. On M1r's own release Q-matter at 10⁻⁶ would need K ≤ 1.38–1.41. On the repaired hybrid at c0 = 1, S/Q = 14.6–72. S/Q < 3.47 appears only with a holder 27–40% above the inherited calibration (or 10% above, on M1r's own history): a choice, not a consequence of two fields |
| Q-matter at Q ~ 10²⁰ | **conditional** [Q] | If a tuned holder gives Q-matter, qm_micro's ripening (with its skeptic's ×0.54–0.70) gives ⟨Q⟩ ≈ 1×10²⁰. That number was computed on M1r's history and is transferred here, not recomputed |
| The Higgs rescued by the plasma | **fails as written; survives on the repaired architecture** [C on E] | The proposal's portal sign makes the Higgs tachyonic during inflation (m² ≈ −9×10⁷ H²), where no plasma can help. Without a portal (HFa/HFb), ξ_H R protects inflation (ξ_min = 0.049–0.053), and afterwards the singlet-neutrino and SM media keep the Higgs heavy to ≥ 174.6 GeV [E] |
| Leptogenesis kept | **unsupported** [A] | Not computed in this pass. The repaired channel uses seesaw singlets (M_N = 1.6–5.2×10¹³ GeV), but Σm_ν (§1 #15) and the anomalous couplings of the axial rotations (mf_hybrid skeptic N6) are open |
| Premises (rolling ⇒ Q-matter; Dirac width; V0^{1/4} = 2.4×10¹⁵; arXiv:1602.07719) | **refuted** [Q/D] | Refuted in the inflaton pass; Majorana Γ = y²mβ³/16π; V0^{1/4} = 7.40×10¹⁴ GeV; 1602.07719 is unrelated (the Higgs-in-preheating paper is 1602.00483) |

**Minimal change that works on the inflationary side [C, mf_gemini G3 + mf_hybrid].** Put V0 on a heavy waterfall field
and the holder on a separate light spectator, with a bias and an axial singlet channel. That is mf_hybrid's
architecture, not the proposal. Even then the "fast" requirement (K ≈ 2) is not met once Pauli blocking is included.

---------------------------------------------------------------------------------------------------------------------

## 3. Benchmark table

Conventions: a = 1 and H = H_i at the slice σ = σ_c (static-slice analogue of M1r's, ε_H ≈ 10⁻⁷, R = 12); K = a²H/H_i in
radiation; M_P = 2.435×10¹⁸ GeV reduced; e is the hidden U(1) coupling; c0 is the holder amplitude in units of M1r's
10i calibration (μ_eff(slice) = μ + c0(10i − μ)) at the 10i phase. Portals vanish at M_Y (c_L = 0) unless stated.
Retrospective pulls: Planck 2018 n_s = 0.9649 ± 0.0042 (arXiv:1807.06211), ACT DR6 P-ACT-LB 0.9743 ± 0.0034
(arXiv:2503.14452).

Columns:
- **P:** the proposal as written.
- **A / B:** HFa / HFb as exported (unblocked axial reheating, c0 = 1).
- **C:** HFa (HFb) Pauli Γ_th = 0, self-consistent, c0 = 1.
- **D:** HFa Pauli Γ_th = 0, strong holder.
- **E:** HFa Pauli Γ_th = 10, strong holder (this track).
- **M:** M1r's supplied (3, 0.5) history (this track).

**Inflation, waterfall, reheating, spectator**

| | P | A (B) | C | D | E | M |
|---|---|---|---|---|---|---|
| architecture | λ/4(χ² − v²)² + ½g²σ²χ² + V_SR, v = 1.633 M_P, λ = 4.81×10⁻¹⁵, m_χ = 3H0 | heavy small-field hybrid: m_χ = 250, m_σ = 600 H_i; λ = 2.32×10⁻⁷, v = 0.0196 M_P, g = 1.63×10⁻³, σ_c = 0.0058 M_P; tilted-hilltop clock with χ loop (κ = 10); bias c = 10⁵ H_i³; spectator m_s = g_sχ (B: m_s² = 9 − ¾R) | as A, rescaled to its fixed point | as C, at 8.4×10¹⁰ | as C, Γ_th = 10 H_i | phenomenological oscillator m_I = 3, Γ_I = 0.5 H_i |
| holder | not working as stated (on (χ − v): −v domains Y-condensed) | s: \|g_h\|s_* = 12.30 H_i (c0 = 1), s_* = 7.3×10¹⁵ GeV | c0 = 1 | **c0 = 1.36 / 1.40** | c0 = 1.25 / 1.30 | inflaton trilinear, c0 = 0.9–1.4 scanned |
| n_s; α_s | ≥ 0.972 (log) or hilltop | 0.9649 (0.9743); +1.1 (+1.4)×10⁻⁴ | 0.9649 retuned (MS not redone) | retuned (MS not redone) | retuned | — |
| r | 2.75×10⁻⁷ | 2.75×10⁻⁷ | 9.4×10⁻⁷ | 1.27 / 1.15×10⁻⁷ | 4.4×10⁻⁸ (at 5.2×10¹⁰) | 2.7×10⁻⁷ (at c0 = 1) |
| H_i (GeV); self-consistent? | 1.25–1.55×10¹¹ asserted; **no** | design 1.30×10¹¹; yes at 1.28×10¹¹ (track map) / ≈ 1.11×10¹¹ (gap-filled) | **2.40×10¹¹, yes** (3 steps); HFa CMB capacity 48.6 < N_* at κ = 10 (needs κ = 15–50); HFb passes | 8.85×10¹⁰ (ratio 1.05) / **8.40×10¹⁰ (0.998)**; capacity 99.8 | ≈ 5.2×10¹⁰ (c0 = 1.25, ratio 1.05); c0 = 1.30: ratio 0.85 | by calibration: 1.52 / 1.29 / 0.96 / 0.72×10¹¹ at c0 = 0.9 / 1 / 1.1 / 1.2 |
| V0^{1/4} (GeV) | 7.40×10¹⁴ | 7.40×10¹⁴ | 1.01×10¹⁵ | 6.0×10¹⁴ | 4.7×10¹⁴ | — |
| waterfall; defects | mild, ≥ 17 e-folds; PBH overclosure; walls | fast (ε_H = 1 at 0.19 H_i⁻¹); no PBH (ζ ≲ 6×10⁻³); **walls not excluded** (fragments at t ≈ 0.45 with Pauli friction) | same | same | same (ε_H = 1 at 0.25 at 5×10¹⁰) | none (single field) |
| reheating channel | y_N χN̄N + portal: blocked / tachyonic | axial −(M/2)N̄e^{iγ5φ/f}N to 4 seesaw singlets → LH; **unblocked (Pauli ignored)** | same, Pauli-blocked fluid (skeptic's [E] model) | same | same | Γ_I = 0.5 H_i (no realisable channel) |
| K | needs ≤ 1.38–1.41; best ≥ 2.2 | 2.05 | 7.84 | 6.50 | 3.19 (8.4×10¹⁰), 2.49 (5.0×10¹⁰) | 1.986 |
| T_RH (GeV) [C, K7]: radiation ≥ 50% / SM ≥ 50% | — | 2.4×10¹⁴ (t = 0.29) / 6.7×10¹³ (t = 7.3) | 1.9×10¹³ (t ≈ 205–211) | 1.1×10¹³ (t ≈ 195–211) | 1.5–1.6×10¹⁴ (t = 0.3–0.5) / 2.0–2.2×10¹³ (t = 29–60) | T_K = 4.3×10¹⁴ (radiation temperature extrapolated to a = 1) |
| spectator s_*; β_iso | — | 7.3×10¹⁵ GeV; 6.9×10⁻³ (track map) / 0.017–0.022 (gap-filled) | ≤ 4.5×10⁻³ (D ≈ 0) | D = 9.2: β = 0.12 at 7.3×10¹⁵ → **s_* ≥ 1.36×10¹⁶ GeV** | not computed (expect D ~ 5–9, s_* ≳ 10¹⁶) [E] | holder on the inflaton (adiabatic) |
| spectator energy; decay | — | 1.35×10⁻⁵ after turn-on; 1.3% at decay; Γ_s = 1.3×10⁻⁶ H_i via N_s (needs soft mass + seesaw, else ΔN_eff ≈ 0.04) | 0.1–0.2% at decay | ×3.5 with the larger s_*, ≲ 1% | — | — |
| ACT DR6 pull (Planck) | −0.7σ (+1.7σ) at 0.972 | −2.8σ (0.0σ); B: 0.0σ (+2.2σ) | as A/B | as A | as A | — |

**Chain, relic, Higgs** (e = 10⁻⁶ / 3×10⁻⁶)

| | P | A (B) | C | D: c0 = 1.36 / 1.40 | E: c0 = 1.25 / 1.30 (5.0×10¹⁰) | M: c0 = 1 / 1.1 / 1.2 |
|---|---|---|---|---|---|---|
| S/Q | not run (no history) | 14.69 / 12.25 (14.58 / 12.15) | 61.5 / 30.2 (HFb 60.0) | 4.16 / 3.28 ; **2.875 / 2.12** | 3.64 / 3.13 ; 2.87 / 2.42 | 4.955 / 2.99 ; 3.51 / 2.40 ; 2.12 / 1.64 |
| e_crit (bulk / Q = 2×10⁵) | — | none | none | 2.18 / 7.7×10⁻⁶ ; **0.60 / 0.95×10⁻⁶** | 1.32 / 6.8×10⁻⁶ ; 0.45 / 0.91×10⁻⁶ | 2.06 / 3.18×10⁻⁶ ; 1.02 / 1.56×10⁻⁶ ; 0.32 / 0.45×10⁻⁶ |
| branch at 10⁻⁶ | — | free X | free X | free X ; **Q-matter** | free X ; Q-matter (H_req/H_i = 0.85) | free X ; free X (×1.01) ; Q-matter |
| relic: mass per charge | — | 1.58 (gap-filled 1.77) m_X, ≈ 57–64% asymmetric | 6.8 m_X, mostly symmetric | 1.000 m_X (asymmetric) ; **0.807 m_X** | 1.000 ; 0.807 | 1.002 ; 1.000 ; 0.798 |
| Q-matter ⟨Q⟩, object mass | — | — | — | c0 = 1.40: ⟨Q⟩ ≈ 1×10²⁰ [Q qm_micro + skeptic, transferred], object ≈ 0.6 M_Pl (≈ 14 μg), vapour few % | as D [Q] | c0 = 1.2: as D [Q] |
| dominant heat (R_vec / R_Y / R_casc at 10⁻⁶) | — | 4.6 / 9.4 / 0.55 | 61.6 / 10.7 / 0.46 | 4.8 / 9.9 / 5.8 ; 4.0 / 5.2 / 6.4 | 1.6 / 5.7 / 2.1 ; 1.3 / 5.8 / 3.9 | 3.6 / 0.38 / 0.02 ; 2.3 / 0.65 / 0.02 ; 1.3 / 0.78 / 0.02 |
| Higgs ξ_H = 0.1: inflation | portal sign makes the Higgs tachyonic: fails | safe, ξ_min = 0.049 / 0.053 at M_t = 172.57 / 172.86 | safe | safe (lower H) [A] | safe [A] | safe (decisive pass) |
| critical M_t, gravitational only (GeV) | — | 172.88 (export), 172.92 (Pauli at 1.3×10¹¹) | 172.74 (fails at +1σ) | ≈ 173.0 [E, K6] | ≈ 173.2 [E] | 172.79–172.89 (decisive pass) |
| with N-gas + SM medium masses | — | safe to ≥ 174.6 | safe to ≥ 174.6 | safe [E, transferred] | safe [E, transferred] | needs a plasma within τ ≈ 8–10 |
| open inputs | — | walls (lattice); Pauli blocking (refuted K); Σm_ν; anomalous couplings; N_s completion; κ; axial UV completion | Γ_th; walls; κ at the fixed point; MS at the new scale | **c0 (±3%) and phase (±5°)**; s_* ≥ 1.36×10¹⁶; c_L ≲ 1/16π²; Γ_th; walls; MS at 8.4×10¹⁰ | c0; **axial EFT at/over its cutoff (1.05–1.5)**; iteration not converged | the inflaton itself (not realisable as a single field; inf_critic) |

Notes on the table:
- Columns D–E use the full chain (mf_ideal's engine, vectors, dq_heat ledger, closed form, p map).
- D's c0 = 1.36 row is the skeptic's [Q]; c0 = 1.40 is reproduced here [C].
- The Q-matter ⟨Q⟩ was computed on M1r's own thermal history and is transferred, not recomputed.
- At the strong holder N_YH = 7–9 (×800 M1r's), so the cascade heat R_casc ≈ 4–6 is a large part of the ledger. The
  branch survives it only because the release is much more efficient (ε ≈ 0.45–0.55, E/E_c ≈ 2.3–2.9).

---------------------------------------------------------------------------------------------------------------------

## 4. Does any benchmark recover the M1r picture? Does any give Q-matter at e = 10⁻⁶?

**M1r's picture** means a branch point in e, a CP-set asymmetric relic, and r ~ 10⁻⁷.

**At the inherited holder calibration (c0 = 1): no [C].** Every hybrid history is free X at every e, with no e_crit.
Two things drive this:
- the fast spectator turn-on at an abrupt end leaves N_YL = 5–11 (M1r: 0.57);
- Pauli blocking raises K and the vector heat.

### 4a. The same freedom on M1r's own history (`code/m1r_c0.py`, K10) [C]

**Method.**
- mf_ideal's bridge table (M1r's recorded (3, 0.5) background and holder c = I/I₀) is copied with its c and ċ columns
  multiplied by c0.
- The preparation is re-solved at u = μ + c0(10i − μ), R = 12.
- Then the full chain.
- **Regression.** c0 = 1 gives S/Q(10⁻⁶) = 4.9550 and e_crit = 2.0578×10⁻⁶ (bulk) / 3.184×10⁻⁶ (Q = 2×10⁵), the
  decisive pass's digits.
- Every preparation is a minimum (lightest physical mass² 8.8–15.8 H_i² > 9/4), with Gauss ≤ 10⁻¹⁰ and N_X frozen to
  ≤ 0.12%.

| c0 | 0.9 | 1.0 | 1.05 | 1.1 | 1.2 | 1.3 | 1.4 |
|---|---|---|---|---|---|---|---|
| ε; N_YL | 0.116; 0.20 | 0.127; 0.57 | 0.137; 0.89 | 0.151; 1.43 | 0.205; 2.76 | 0.302; 3.57 | 0.432; 42.8 |
| S/Q(10⁻⁶) / (3×10⁻⁶) | 6.07 / 3.40 | 4.955 / 2.99 | 4.25 / 2.71 | 3.51 / 2.40 | 2.12 / 1.64 | 1.21 / 0.99 | 2.00 / 1.88 |
| e_crit (bulk) | 2.86×10⁻⁶ | 2.06×10⁻⁶ | 1.54×10⁻⁶ | 1.02×10⁻⁶ | 3.2×10⁻⁷ | < 10⁻⁷ | < 10⁻⁷ |
| H_i the dark matter needs (GeV) | 1.52×10¹¹ | 1.29×10¹¹ | 1.14×10¹¹ | 9.6×10¹⁰ | 7.2×10¹⁰ | 4.5×10¹⁰ | 3.0×10¹⁰ |

**Reading.**
- e_crit ∝ c0^−4.9 near c0 = 1. The bulk threshold reaches e = 10⁻⁶ at **c0 = 1.102**.
- So the decisive pass's branch at the declared e = 10⁻⁶ is set by the supplied holder strength to ±10%, exactly as by
  e itself. That strength was copied from Codex's Option 1A calibration (μ_eff = 10i), not derived.
- The honest statement of M1r is a branch surface in (e, c0): free X for e < e_crit(c0), with
  e_crit = 2.06×10⁻⁶ (c0/1)^−4.9 locally.
- Q-matter at 10⁻⁶ on M1r's own history needs c0 ≳ 1.10, with H_i ≲ 9.6×10¹⁰ GeV and r ≲ 1.5×10⁻⁷.

### 4b. The hybrid at a strong holder, both ends of the Pauli bracket (`code/hyb_c0.py`, K9/K11) [C on the [E] Pauli fluid]

**Γ_th = 0 (the skeptic's corner; K = 6.50 at 8.40×10¹⁰ GeV).** The point is self-consistent within 5–9% for
c0 = 1.36–1.40 at H_i = 8.4–9.2×10¹⁰ GeV:
- c0 = 1.36: free X at 10⁻⁶ (S/Q 4.16, asymmetric, M = 1.000), Q-matter at 3×10⁻⁶, e_crit = 2.2×10⁻⁶. **This is M1r's
  picture.**
- c0 = 1.40: **Q-matter at 10⁻⁶** (S/Q 2.875, objects at Q = 2×10⁵ allowed), e_crit = 0.60×10⁻⁶, H_req/H_i = 0.998.
  It is robust to the cycle-mean R, to rtol and to an independent release solver [Q skeptic]. My rebuilt history
  reproduces it to all printed digits.
- The axial EFT ratio μ5,max/4πf is 0.65–0.76 there.

**Γ_th = 10 H_i (new here).** K = 3.19 at 8.40×10¹⁰ and 2.49 at 5.0×10¹⁰ GeV.

| H_i of the history | c0 | S/Q 10⁻⁶ / 3×10⁻⁶ | e_crit | branch at 10⁻⁶ | H_req/H_i |
|---|---|---|---|---|---|
| 8.4×10¹⁰ | 1.00 | 20.4 / 15.1 | none | free X (M = 2.37) | 1.91 |
| 8.4×10¹⁰ | 1.25 | 6.04 / 5.34 | none | free X | 0.93 |
| 8.4×10¹⁰ | 1.30 | 4.50 / 3.95 | 7.0×10⁻⁵ | free X | 0.67 |
| 8.4×10¹⁰ | 1.36 | 2.35 / 1.89 | 2.8×10⁻⁷ | Q-matter | 0.54 |
| 8.4×10¹⁰ | 1.40 | 1.82 / 1.41 | 1.6×10⁻⁷ | Q-matter | 0.45 |
| 5.0×10¹⁰ | 1.00 | 17.0 / 14.1 | none | free X | 2.45 |
| 5.0×10¹⁰ | **1.25** | **3.64 / 3.13** | **1.32×10⁻⁶** | free X (M = 1.000) | **1.05** |
| 5.0×10¹⁰ | 1.30 | 2.87 / 2.42 | 4.5×10⁻⁷ | Q-matter | 0.85 |
| 5.0×10¹⁰ | 1.33 | 2.38 / 1.94 | 2.5×10⁻⁷ | Q-matter | 0.72 |
| 5.0×10¹⁰ | 1.36 | 1.98 / 1.57 | 1.6×10⁻⁷ | Q-matter | 0.62 |

**Reading.**
- With faster thermalisation, the skeptic's c0 = 1.40 point is no longer self-consistent: the dark matter wants
  H_i = 3.8×10¹⁰, not 8.4×10¹⁰.
- At 8.4×10¹⁰ the self-consistent c0 ≈ 1.24 is free X with no e_crit.
- The M1r-like point moves to c0 ≈ 1.25 at H_i ≈ 5.2×10¹⁰ GeV.
- Q-matter at 10⁻⁶ would need c0 ≈ 1.27–1.30 at H_i ≈ 4.3–5×10¹⁰ GeV.
- There μ5,max/4πf = (0.42–0.49)(1.30×10¹¹/H_i) = 1.05–1.5 under mf_chain's rescaling rule. The axial EFT is at or
  over its cutoff, so a different rule (a different f) and a re-run would be needed.
- The Γ_th = 10 iteration covers two H_i points only and is not converged [E for the fixed points].

### 4c. What that choice fixes and what remains supplied

**Choice D at c0 = 1.40 (Q-matter at e = 10⁻⁶).** It fixes [C on its model choices]:
- **Scale.** H_i = 8.40×10¹⁰ GeV, V0^{1/4} = 6.0×10¹⁴ GeV, r = 1.15×10⁻⁷ (n_t = −r/8).
- **Reheating.** K = 6.50, with the SM bath reaching half the energy at T ≈ 1.1×10¹³ GeV.
- **Branch.** e_crit = 0.60×10⁻⁶ (bulk).
- **Relic.** Q-matter, mass per charge 0.807 m_X. Transferred [Q]: ⟨Q⟩ ≈ 10²⁰, Planck-mass objects, a few-percent
  warm vapour.
- **Isocurvature floor.** s_* ≥ 1.36×10¹⁶ GeV.

The tilt (HFa 0.9649) is put in by hand, not fixed by the choice. It is retuned through the clock, and the
Mukhanov–Sasaki spectrum at the new scale was not recomputed.

**Choice D at c0 = 1.36.** It gives M1r's structure (free X at 10⁻⁶, e_crit 2.2×10⁻⁶, an asymmetric relic) at
H_i ≈ 8.8×10¹⁰ GeV, r ≈ 1.3×10⁻⁷.

**What remains supplied, for both:**
- **The holder.** Amplitude c0 (window ±2–3%) and phase (±5°; purely imaginary μ_eff = 12i–16i stays free X [Q skeptic]).
- **Reheating.** The N thermalisation rate Γ_th, which moves the window (§4b), and the Pauli fluid model itself [E].
- **Defects.** Walls (not excluded).
- **Clock.** The tree quartic κ.
- **Singlet sector.** The axial f_j and M_j and their UV completion (periodic couplings on non-periodic fields;
  anomalous couplings); the seesaw masses (Σm_ν); the N_s completion.
- **Waterfall.** The bias c.
- **Hidden sector.** The sequestering of σ and χ from it; the induced-portal renormalisation (c_L ≲ 1/16π²).
- **Couplings.** e itself.
- **Spectator.** s_*.
- **Higgs.** Survival rests on medium masses [E].

**Summary.** The hybrid can host M1r's picture, but only as a tuned holder choice whose location depends on an
uncomputed rate. M1r's own history hosts the same picture with a looser tune (±10%).

---------------------------------------------------------------------------------------------------------------------

## 5. The brief's estimates G1–G7: final status

| | Status |
|---|---|
| G1 amplitude √6 H M_P/m; blocking; loop | **verified** [D/C, K1: Φ = 1.22 / 0.82 / 0.61 M_P; y = 2.894 for Γ = 0.5H] |
| G2 proposal numbers (v, λ, CW ×2×10¹³, blocking ×4×10⁷, mild waterfall ~17 e-folds, walls with v ~ M_P) | **verified**; "a holder on χ is off during inflation" **refuted as stated** (on (χ − v) it is on, but −v domains end Y-condensed) |
| G3 heavy waterfall decays perturbatively | **algebra verified** (λ = 0.013, m = 3.85×10³ H, v = 1.27×10⁻³ M_P). A general bound follows: holding V0 and decaying in a Hubble time needs m ≥ (384π)^{1/5}(M_P/H)^{2/5}H = 3.4×10³ H at 1.3×10¹¹ [D, K2]. The benchmarks reheat instead through axial couplings |
| G4 σ holds V0 after the waterfall | **confirmed for m_χ ≫ m_σ** (mf_hybrid); the benchmarks avoid it with m_σ > m_χ |
| G5 spectator holder | **coefficients verified** (1.5(s_*/M_P)², ξ_s = −¾). s → −s flips the sign but is not the CP image. The turn-on (≈ 0.2 H_i⁻¹) is too fast for a quiet release |
| G6 Higgs sign structure | **confirmed**; window −2×10⁻⁸ ≲ λ_χH ≲ −7×10⁻¹² for a hypothetical portal; the benchmarks have none |
| G7 R_vec ∝ K², "×4 cut" | R_vec ∝ K^1.8–1.95 **verified**; the "×4 cut" is **refuted**: the realistic K rises to 2.5–10, and the branch is set by Y_L and c0, not by the vectors |

---------------------------------------------------------------------------------------------------------------------

## 6. The single most important remaining open item

**A computed post-waterfall history.** That means a 3D lattice (≤ 64³) of the biased real waterfall together with the
axial singlet channel, including Pauli blocking and the N thermalisation rate Γ_th, run until the χ, σ condensates are
gone. It decides three things at once:
1. **Whether the hybrid survives.** The field fragments at t ≈ 0.45 H_i⁻¹. Percolating walls would dominate at t ~ 23,
   before the bias removes them at t ~ 150.
2. **K.** It is 2.5–10 across the bracket, and with it the self-consistent H_i (5×10¹⁰ to 2.4×10¹¹ GeV), the vector heat
   and the Higgs plasma timing.
3. **Where a branch point can sit.** The c0 window that hosts M1r's picture is c0 ≈ 1.36–1.40 at H_i ≈ 8.6×10¹⁰ inside
   the axial EFT (Γ_th = 0), against c0 ≈ 1.25–1.30 at H_i ≈ 5×10¹⁰ at the EFT edge (Γ_th = 10).

Even with that history in hand, the branch at the declared e = 10⁻⁶ is a supplied choice of the holder amplitude, on the
hybrid and on M1r's own history alike (§4a). The program should state M1r's branch as a surface in (e, c0) unless a
principle fixing |g_h|s_* (or g_I) is supplied.

**Runner-up open items (ranked).**
1. **A principle (or a declared benchmark) for the holder amplitude and phase.** It decides the branch at ±10% (M1r) or
   ±3% (hybrid).
2. **The singlet sector's consistency.**
   - Σm_ν for four or five seesaw singlets (0.20–0.25 eV for a real R, in retrospective tension with 0.064–0.12 eV);
   - the anomalous couplings of the axial rotations (≈ 10² H_i), with SU(2) preheating and possible
     spontaneous lepto/baryogenesis;
   - the N_s completion (otherwise dark radiation).
3. **The induced holder portal's renormalisation condition (c_L).** At the strong holder a threshold-size c_L = 1/4π²
   turns Q-matter back into free X.
4. **The exact Mukhanov–Sasaki spectrum and the κ condition** at the rescaled scales (HFa's clock fails at κ = 10 at
   2.4×10¹¹).
5. **The pair-fraction map in the gap 6.74 < S/Q < 17.1** (now filled by the mf_chain skeptic's network rows; adopt
   them).
6. **Q-matter microphysics (⟨Q⟩) on a hybrid thermal history**, not transferred from M1r's.

---------------------------------------------------------------------------------------------------------------------

## 7. Recommendations to the combiner (WORKING_THEORY)

- **§0, §4b, §9: the branch.** Add that the branch at e = 10⁻⁶ is set jointly by e and the supplied holder amplitude.
  On M1r's own history c0 = 1.10 already gives the Q-matter threshold at 10⁻⁶, and c0 = 0.9 moves e_crit to
  2.9×10⁻⁶ [C here]. "Free X at 10⁻⁶" is a c0 = 1 statement.
- **§4a: add a multi-field row.**
  - The proposal as written fails.
  - The repaired hybrid (HFa/HFb) passes CMB, PBH and Higgs (given the medium), with a realisable axial-singlet channel.
  - It reheats with K = 2.5–10, not 2, once Pauli blocking is included, and walls are not excluded.
  - At c0 = 1 it is free X at every e (S/Q 15–72).
  - With c0 ≈ 1.36–1.40 (±3%, ±5°) on the Γ_th = 0 corner it recovers M1r's picture or Q-matter at 10⁻⁶, at a
    self-consistent H_i ≈ 8.4–8.9×10¹⁰ GeV with r ≈ 1.2×10⁻⁷.
  - That point moves to H_i ≈ 5×10¹⁰ (axial EFT at its cutoff) for Γ_th = 10.
- **§8: re-rank the open items.**
  1. The post-waterfall lattice with the singlet channel.
  2. The holder-amplitude principle.
  3. Σm_ν and the anomalous couplings.
  4. c_L.
- **§7: CP bookkeeping.** On the hybrid the CP sign is α + sign(s_*), and s → −s is not the CP image [Q].

---------------------------------------------------------------------------------------------------------------------

## 8. Limits

- **Inherited layers.** Every number inherits its track's model choices:
  - the homogeneous classical release;
  - the single-Higgs-direction Proca vectors;
  - dq_heat's ledger;
  - the closed-form S/Q;
  - the p(S/Q) recipe (gap noted);
  - qm_micro's M table and ⟨Q⟩;
  - the mf_hybrid skeptic's Pauli fluid (harmonic, constant-amplitude saturation; Γ_th unknown);
  - mf_chain's rescaling rule;
  - the Higgs D criterion and the medium masses.
- **This track's estimates [E]:**
  - the M_t,crit interpolation in ln H_i across histories of different K (K6);
  - the self-consistent locations on the Γ_th = 10 bracket (two H_i points, no converged iteration);
  - the EFT ratio scaling (quoted from the mf_chain skeptic);
  - the expected D on the Γ_th = 10 strong-holder points (not computed).
- **The M1r c0 scan (K10)** uses the bridge table with the inflaton holder scaled; the holder's back-reaction on the
  inflaton is neglected, as in M1r. The phase was not scanned on M1r's history.
- **Not computed here:** any lattice; Higgs modes on the strong-holder or Γ_th = 10 histories; exact MS spectra at the
  rescaled scales; HFb at a strong holder (the skeptic found it approximately similar on the exported history);
  isocurvature on the Γ_th = 10 points; leptogenesis.

## Code and results

| file | purpose |
|---|---|
| `code/crc_common.py` | read-only import of mf_ideal's engine drivers, every output redirected to this track's scratch |
| `code/m1r_c0.py` → `results/m1r_c0.json` | K10: the holder amplitude on M1r's own (3, 0.5) history (bridge table scaled by c0), full chain; c0 = 1 regression |
| `code/hyb_c0.py` → `results/hyb_c0.json` | K9/K11: full chain on mf_chain-built hybrid histories at c0 ≠ 1 (Γ_th = 0 regression at 8.4×10¹⁰; Γ_th = 10 at 8.4×10¹⁰ and 5.0×10¹⁰) |
| `code/settle_checks.py` → `results/settle_checks.json` | K1–K9: G1/G2 algebra, perturbative bound, log-family tilt, β_iso, r/V0/EFT ratio, M_t,crit interpolation, temperatures from every history, e_crit(c0), the self-consistency grid |
| `code/verdict.py` → `results/verdict.json` | the benchmark table of §3 |

Reproduce (one core plus mf_ideal's two vector workers; ≈ 25 min):

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_critic/code; export OMP_NUM_THREADS=1
SPC=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_critic
python3 -B m1r_c0.py 1.0 0.9 1.05 1.1 1.2 1.3 1.4
HB="python3 -B hist_build.py HFa"                               # mf_chain's builder, read-only, scratch redirected
(cd ../../mf_chain/code && MFC_SCRATCH=$SPC $HB crP0_H084  --H 8.4e10 --fric pauli --Gth 0  --tails free --tsw 0.45)
(cd ../../mf_chain/code && MFC_SCRATCH=$SPC $HB crP10_H084 --H 8.4e10 --fric pauli --Gth 10 --tails free --tsw 0.45)
(cd ../../mf_chain/code && MFC_SCRATCH=$SPC $HB crP10_H050 --H 5.0e10 --fric pauli --Gth 10 --tails free --tsw 0.45)
python3 -B hyb_c0.py $SPC/hist/bench_crP0_H084.json crP0H084 1.40
python3 -B hyb_c0.py $SPC/hist/bench_crP10_H084.json crP10H084 1.25 1.30 1.36 1.40
python3 -B hyb_c0.py $SPC/hist/bench_crP10_H050.json crP10H050 1.25 1.30 1.33 1.36
python3 -B settle_checks.py; python3 -B verdict.py
```

Each history build takes ≈ 20 s; each chain run (release, vectors, chain) 1.5–2.5 min.

## Literature (checked by web search in this session unless marked)

- **Inflation and preheating:**
  - Linde, hybrid inflation, arXiv:astro-ph/9307002.
  - Felder, García-Bellido, Greene, Kofman, Linde, Tkachev, tachyonic preheating, arXiv:hep-ph/0012142.
  - Sakai, Shinkai, Tachizawa, Maeda, topological inflation threshold η ≳ 0.33 m_Pl, arXiv:gr-qc/9506068.
- **Isocurvature:** Planck 2015 XX, arXiv:1502.02114 (the β_iso < 0.038 axion-I value, as quoted by the mf_chain
  skeptic); Planck 2018 X, arXiv:1807.06211 [R, via the tracks].
- **Data, retrospective:** ACT DR6, arXiv:2503.14452 (P-ACT n_s = 0.974 ± 0.003); BICEP/Keck 2021, arXiv:2110.00483
  [R, via the tracks].
- **Our tracks, read-only:**
  - mf_gemini, mf_hybrid, mf_ideal, mf_chain, mf_higgs, and their skeptics;
  - inf_critic, inf_chain, inf_reheat;
  - qm_micro and its skeptic;
  - dq_entropy and dq_heat (via mf_ideal's engine).
