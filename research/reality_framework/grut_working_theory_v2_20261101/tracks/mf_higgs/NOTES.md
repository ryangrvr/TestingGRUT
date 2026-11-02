# Track mf_higgs: Higgs survival with ξ_H = 0.1 on the hybrid benchmarks HFa / HFb

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md, track mf_higgs; GitHub issue #13 context).

This is same-session AI work on one supplied model. It is not outside replication. Nothing here is a validation, an
exclusion of the real world, a new theory or a GRUT prediction.

**Every choice below is a model choice, not a consequence of closure.** That covers:
- the benchmark (mf_hybrid's HFa/HFb) and its reheating fluid (unblocked, or the skeptic's Pauli-blocked closure);
- ξ_H = 0.1, read as ξ(M_t);
- the survival criteria (Hawking–Moss bubble count, the scale-resolved danger exponent D, the mean-field edge);
- the hard-gas plasma mass;
- the renormalisation conditions of the induced terms;
- the hypothetical G6 portal.

Data comparisons (PDG top mass, α_s) are retrospective. Git was read-only and nothing was committed. Every run was in
the foreground with one process (the machine showed one core), and no compiled files are in the tree.

Labels:
- [D] derived;
- [C] computed with this track's code (`code/`, output in `results/`);
- [E] estimate, or an argued criterion;
- [A] argued;
- [Q] quoted from a reviewed track;
- [R] literature, arXiv numbers checked by web search in this session.

Units and conventions:
- H_i = 1.30×10¹¹ GeV is mf_hybrid's slice value, with a = 1 and H = H_i at σ = σ_c. τ = H_i t is counted from that
  slice. ε_H = 1 falls at τ = 0.193 (N = 0.193).
- N_obs = N_* + 5.40 is the number of e-folds from the exit of today's horizon to the slice. It is 58.05 at K = 2.05,
  57.70 at K = 4.10, 57.43 at K = 7.05 and 56.89 at K = 20.9. D is quoted at N_obs (D₆₀ in `results/`).
- Inputs: M_t = 172.57 ± 0.29 and M_h = 125.20 ± 0.11 GeV (PDG 2024, [R], pdglive). α_s = 0.1180 ± 0.0009.
- ξ_EW = 0.1 runs to ξ(H_i) = 0.0877.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer on one screen

| # | Result | Label |
|---|---|---|
| 1 | **Inflation, H constant.** Over the whole observable inflation the hybrid has H = H_i to 2×10⁻⁶ and R = 12H² to 10⁻⁷. The plateau benchmarks' end-of-inflation drop of R/H² (inf_higgs R2) is absent. The constant-H Hawking–Moss count gives **ξ_EW,min = 0.049 / 0.053 (tree) and 0.042 / 0.048 (one loop) at M_t = 172.57 / 172.86**, with log₁₀⟨N⟩(ξ = 0.1) = −609 / −427 (tree). At the fit-grid edge, 174.6 GeV, ξ_min is still only 0.073 (tree) / 0.069 (loop), so the inflationary cap of skeptic N2 lies **above 174.6 GeV** (above 176 GeV if extrapolated). This holds at α_s = 0.1171 and 0.1189 and for H_i = 0.6–10×10¹¹ GeV (ξ_min(172.86) = 0.051–0.058). | C [E criterion] |
| 2 | **Inflation, pessimistic brackets.** Taking the instantaneous HM rate literally through the sub-Hubble waterfall gives ξ_min = 0.066 / 0.074 and a cap at **174.12 GeV (HFa) / 173.95 GeV (HFb)**. α_s ∓ 0.0009 moves it to 173.78 / 174.47 (HFa). In the bench, R/H² falls from 12 to 6 within the last ≈ 0.02 H⁻¹, which is much shorter than the Higgs relaxation time of ≈ 3 H⁻¹, so this reading is outside the validity of Hawking–Moss [A]. The non-equilibrium Gaussian count against the instantaneous barrier gives log₁₀(patches) = −974 / −696 / −210 at 172.57 / 172.86 / 174.6. | C / A |
| 3 | **Portal masses during inflation.** The benchmark has no tree portal to σ, χ or s, by choice. The axial N give a field-independent one-loop Higgs mass [D]. The two-loop σ–χ phase mismatch is ≲ 0.02 H_i² [E]. The negative Higgs mass² that inflation tolerates at ξ = 0.1 is 0.72 / 0.66 / 0.39 H_i² at 172.57 / 172.86 / 174.6 (tree). **New:** the spectator singlet N_s, if it gets the seesaw Yukawa that the mf_hybrid skeptic (N8) says it needs, induces δm_h² = −(y_νs²/8π²)M_s²[ln(M_s²/μ²) − 1] with M_s = y_s s_* = 258 H_i during inflation. That is −1.8 H_i² for M_soft = H_i under "vanish at M_soft", which is fatal for ξ = 0.1 unless μ ≳ 25 H_i or M_soft ≲ 0.3 H_i. | D / E |
| 4 | **Waterfall: R windows and lines.** After ε_H = 1, R < 0 for 38% of the time. Each window lasts 0.005 H⁻¹ and reaches R_min/H² = −4.7 (exported) or −5.2 (Pauli), with a tachyonic exponent of 0.0023–0.0026 at ξ(H_i). The 2m_χ and 2m_σ lines have Mathieu q = 2–4×10⁻⁶ and 4×10⁻⁸ (ω/H ≈ 565). Switching the rebuilt lines off changes D by ≤ 0.04. | C |
| 5 | **Plasma timing from N → L H** (this track's fluid with SM number density, reproducing the skeptic's K = 2.03 / 4.10 / 7.05 / ≈ 21). The hard-gas induced mass m_ind² = c ρ_SM/E_SM² (c = 0.0228, E_SM ≈ 4×10³ H_i) reaches **≥ 1 H_i² at τ = 0.191, ≥ 10² at 0.21 and ≥ 10³ at 0.24–0.29** (0.29–0.54 for m̃ = 0.015 eV). That is already at ε_H = 1 (τ = 0.193), before the first R < 0 window (τ = 0.194 Pauli, 0.196 exported), and 10–60 times earlier than every deadline. ρ_SM/ρ ≥ 10⁻³ arrives at τ = 0.26 (unblocked) and 0.30–0.33 (Pauli); ρ_SM ≥ 50% at 7.6 (unblocked) and 155–213 (Pauli). | C on E |
| 6 | **D at ξ_EW = 0.1 after inflation, gravitational only** (no plasma), at M_t = 172.28 / 172.57 / 172.86. Exported (K = 2.05): **−1656 / −299 / −7.8**. Pauli Γ_th = 0 (K = 7.05): **−1658 / −349 / −21.9**. Γ_th = 10 (K = 4.1): −1657 / −320 / −13.7. Tails blocked (K = 21): −1658 / −351 / −22.8. HFb equals HFa to ±0.1. Criterion band at 172.86: [−56, +12] (exported) and [−69, −9] (Pauli). | C [E criterion] |
| 7 | **Critical top mass, gravitational only.** **172.88 GeV (+1.07σ, exported or unblocked) to 172.91–172.92 (+1.2σ, Pauli)**, with criterion band 172.83–173.00. M_h ∓ 1σ gives 172.84 / 172.99. **α_s = 0.1171 / 0.1189 gives 172.56–172.60 / 173.25–173.29.** At α_s = 0.1171, central M_t gives D = +6.8 on the exported history: it fails. N6: the critical M_t falls with H_i, from 173.12–173.16 at 0.6×10¹¹ to 172.74–172.77 at 2.1×10¹¹ and 172.54–172.56 at 5×10¹¹. At mf_chain's interim Pauli calibration (H_i ≈ 2.2–2.5×10¹¹ [Q, unreviewed]) it is ≈ 172.74 (+0.6σ). | C |
| 8 | **D with the N → L H plasma** (field grid to the abundance limit x = 3N + 40, skeptic N1). **No region within the abundance limit rolls at any M_t from 172.86 to 174.6 GeV**, on every history, at α_s ± 1σ, and for H_i = 0.6–5×10¹¹. At 172.28 / 172.57 the finite D ≤ −2100 comes only from h̄ beyond the abundance limit (x ≈ 5000–8000). With the hard-gas coefficient cut by 10⁴, D ≤ −107 at all M_t ≤ 174.6. With a thermal bath only (pessimistic LPM time τ_th = 3.1–8.7) the reach is **173.11–173.44 GeV**. | C on E |
| 9 | **ξ_H window** (ξ_EW; at central M_t, with +1σ in brackets). Lower edge: inflation **0.049 (0.053)** tree, 0.042 (0.048) loop; instantaneous bracket 0.066 (0.074). The post-inflation lower edge is **0.066–0.067 (0.095–0.098) gravitational only**, and none (< 0.02) with the plasma. Upper edge (gravitational mean field 2σ_prod = h_max): **ξ_EW ≈ 113–155 (21–24)**, with ≈ 410 at −1σ. The scale-resolved D is ≤ −1000 up to ξ(H_i) = 1000, with no Hartree runaway. With the plasma barrier included, the produced 2σ stays ≤ 1/10 of the barrier up to ξ(H_i) = 1000 [E]. ξ_H = 0.1 lies inside every version: by about ×2 with the plasma, and by only 2–5% at +1σ without it. | C [E] |
| 10 | **Portal.** HFa/HFb have none (λ_χH = λ_σH = λ_sH = 0), so no portal stabilises the Higgs on the benchmarks. The protection is ξR during inflation and then the decay plasma. **Hypothetical G6 portal** λ_χH(χ² − v²)\|H\|² with λ_χH < 0 [C on E]: m_H² = Λ H_i²(1 − χ²/v²), with \|λ_χH\| = 7.4×10⁻¹² Λ. It is ΛH² during inflation, sign-oscillating at ω ≈ m_χ after the waterfall, and positive on average, Λ f_χ H²/4. On both histories, Λ = 1–3000 gives D ≤ −105 at every M_t ≤ 174.6 with no plasma; Λ = 0.3 fails above ≈ 173.3; Λ ≥ 10⁴ gives a Hartree runaway at τ = 0.22–0.33 (resonant/tachyonic Higgs production). **Window: −2×10⁻⁸ ≲ λ_χH ≲ −7×10⁻¹².** | C on E |

**Bottom line.** On HFa/HFb, ξ_H = 0.1 is decided by inflation and by when the Standard Model (SM) decay plasma arrives.
- **Inflation.** With H constant it is safer than on the single-field plateaus. ξ_min = 0.053 at +1σ, against 0.079 for
  HT6t, and the cap lies above 174.6 GeV, against 173.75 for HT6t.
- **Plasma.** The N → L H gas makes the Higgs heavy (m_ind² ≳ 10–10³ H_i²) by the end of inflation. With the hard-gas
  mass, the post-inflation danger vanishes up to at least 174.6 GeV.
- **Without a plasma.** The gravitational-only verdict is as marginal as M1r's: M_t,crit = 172.88–172.92. It fails at
  central M_t for α_s −1σ, and at +1σ once H_i ≳ 2×10¹¹.
- **What the safety rests on.** It holds only if three things hold: the hard-gas mass [E] (or thermalisation by τ ≈ 3–9,
  which gives 173.1–173.4 GeV); a renormalisation condition for any N_s seesaw that keeps its induced inflationary mass
  ≥ −0.66 H_i²; and the axial reheating channel itself. The skeptic's EFT and UV caveats on that channel stand.

---------------------------------------------------------------------------------------------------------------------

## 1. The estimates of BRIEF_MULTIFIELD.md that touch this track

| | Claude's estimate | Verdict here | Label |
|---|---|---|---|
| G6 | A portal linear in a fluctuation oscillates in sign. Positive χ²\|H\|² or (v² − χ²)\|H\|² keeps the Higgs heavy in inflation | **Confirmed and sharpened.** After the electroweak tuning, λχ²\|H\|² ≡ λ(χ² − v²)\|H\|², so the useful sign is λ_χH < 0. That gives +\|λ\|v² during inflation. After the waterfall, (χ² − v²) oscillates at the *fundamental* m_χ with amplitude 2Λ(A/v). Its mean is +Λ(A/v)² > 0, from the virial ⟨χ(χ² − v²)⟩ = 0 [D], which equals Λ f_χ H²/4. The window is row 10. The broad-resonance edge q = 2Δ/ω² ≈ 1 sits at Λ ~ 10⁴, in agreement with the computed runaway [C] | D / C |
| G6 (mf_hybrid) | "Only Higgs tachyonic window ξR ≥ −0.45 H² for ≲ 0.05 H⁻¹, exponent ≈ 0.03" | **Refined.** It is a train of windows: ξR_min = −0.41 to −0.45 H², each lasting 0.005 H⁻¹ with exponent 0.0023–0.0026, occupying 38% of the time. They do not compound (q ~ 10⁻⁶), and the plasma mass already exceeds them by 20–50× at the first window and by 10³–10⁴× from τ ≈ 0.25 | C |
| G7 | R_vec ∝ K², fast reheating | Not this track's quantity. For the Higgs: with the plasma, the result does not depend on K over 2–21. Gravitational only, a larger K (more matter-like) is slightly *safer* (+0.04 GeV) | C |
| G3 | N with M_N ~ 10³ H_i decays to L H "within the Higgs deadline" | **Confirmed for the benchmark's N** (M = 120–400 H_i, rest width 0.12–1.4 H_i, lab width ≈ 0.07 H_i): the SM appears at τ ≈ 0.19, against deadlines ≥ 2.3 at M_t ≤ 173.75 | C |

## 2. Inflation [C on an argued criterion] (`code/mh_inflation.py` → `results/inflation.json`, `inflation_sens.json`, `inflation_tail.json`)

- **Background.**
  - The bench gives ε_* = 1.8×10⁻⁸, H_* = 1.0000023 H_i and R = 12.00001 H² from t = −6.1 to 0.15. The drop to
    R/H² = 6 at ε_H = 1 (t = 0.193) happens within the last ≈ 0.02 H⁻¹.
  - N0 = N_obs + 0.193 = 58.24 (exported) or 57.62 (Pauli).
- **Criteria.**
  - The constant-H Hawking–Moss count ⟨N⟩ = (4π/9)e^{3N0 − B}. The barrier comes from inf_higgs's `barrier_B` (tree
    RG-improved, or the one-loop de Sitter potential of arXiv:1804.02020 [R]) on higgs_reheat's smpot.
  - This is the criterion of Mantziris–Markkanen–Rajantie, arXiv:2011.03763 [R], here with constant H.

| M_t | 172.28 | 172.57 | 172.86 | 173.15 | 173.45 | 173.75 | 174.05 | 174.35 | 174.6 |
|---|---|---|---|---|---|---|---|---|---|
| ξ_EW,min const-H, tree / loop | 0.044 / 0.034 | **0.049 / 0.042** | **0.053 / 0.048** | 0.057 / 0.052 | 0.061 / 0.056 | 0.064 / 0.060 | 0.067 / 0.064 | 0.070 / 0.067 | 0.073 / 0.069 |
| instantaneous HM (bracket), tree | 0.057 | 0.066 | 0.074 | 0.081 | 0.087 | 0.093 | 0.099 | 0.104 | 0.108 |
| log₁₀⟨N⟩(0.1) const-H tree / inst | | −609 / −157 | −427 / −93 | | −247 / −32 | | | | −111 / +13 |
| tolerated −δm_h² at ξ = 0.1 (H_i², tree / loop) | | 0.72 / 0.82 | 0.66 / 0.74 | | 0.55 / 0.62 | | 0.46 / 0.52 | | 0.39 / 0.44 |

- **HFb** is identical in the constant-H part. Its waterfall bracket is slightly worse (cap 173.95 against 174.12).
- **α_s (N5).** ξ_min (tree) at 172.57 / 172.86:
  - 0.054 / 0.058 at α_s = 0.1171;
  - 0.042 / 0.048 at α_s = 0.1189.
  - The constant-H cap stays above 176 GeV either way.
  - The instantaneous-bracket cap is 173.78 / 174.12 / 174.47 GeV at α_s = 0.1171 / 0.1180 / 0.1189.
- **H_i (N6), constant H.** ξ_min(172.86) (tree) is 0.051 at 0.6×10¹¹, 0.053 at 1.3×10¹¹, 0.055 at 2.1–2.4×10¹¹,
  0.056 at 5×10¹¹ and 0.058 at 10¹² GeV. At 2.4×10¹¹ the instantaneous bracket cap is 174.02 GeV.
- **Why the plateau benchmarks were worse.** On those the drop of R/H² from 12 to 6 took ≈ 1 e-fold, and the literature's
  instantaneous count is dominated by the end (inf_higgs R2, its skeptic N2).
  - In the hybrid the same drop is sub-Hubble. Applying HM there assumes the Higgs distribution re-equilibrates within
    0.02 H⁻¹, when its relaxation time is ≈ 3H/m² ≈ 3 H⁻¹ [A].
  - The baseline is therefore the constant-H count, and the instantaneous one is the pessimistic bracket.
  - The non-equilibrium Gaussian count (mode variance at k < aH against the instantaneous tree barrier, through the
    waterfall) gives log₁₀ = −974 / −696 / −419 / −283 / −210 at 172.57 / 172.86 / 173.45 / 174.05 / 174.6 (HFa; HFb
    within 30).
- **Portal and induced masses** (`results/induced.json`) [D/E]:
  - **Axial N₁–₄.** |m_N| = M_j for every field value, so the N–L–H loop is a constant y_ν²M²/8π² = 2–250 H_i². That
    constant is absorbed in the electroweak tuning (Vissani's seesaw naturalness cost, hep-ph/9709409 [R]). It gives no
    inflationary or post-inflation Higgs mass.
    - N₂,₃,₄ share χ/f, so their relative phases vanish.
    - The σ–χ phase mismatch enters at two loops, at ~(y_ν²)²M₁M_j/(16π²)² = 0.002–0.018 H_i² (one log omitted). It is
      frozen during inflation (≤ 3% of the tolerance) and rotates at ~10⁴ H_i afterwards.
  - **N_s (not in the bench; the skeptic's N8 fix).** If N_s gets a soft mass M_soft and a seesaw Yukawa
    y_νs² = m̃M_soft/v², the field-dependent mass M_s = y_s s_* = 258 H_i during inflation induces the R-even
    (s-even) term of inf_higgs skeptic N4: δm_h² = −c[ln(M_s²/μ²) − 1].

    | M_soft (H_i) | c (H_i²) | vanish at μ = M_soft | at μ = M_s(inflation) | at μ = M_P | μ needed for ≥ −0.66 H_i² |
    |---|---|---|---|---|---|
    | 0.3 | 0.054 | −0.68 | +0.054 | +1.27 | ≥ 0.37 H_i |
    | 1.0 | 0.18 | **−1.8** | +0.18 | +4.2 | ≥ 25 H_i |
    | 1.5 | 0.27 | **−2.5** | +0.27 | +6.4 | ≥ 47 H_i |

    Under the program's "vanish at the heavy mass" style of condition this is fatal for ξ_H = 0.1 during inflation,
    whatever the plasma does afterwards. It is a renormalisation-condition choice, like the hidden-sector portals.

## 3. Post-waterfall history and plasma timing [C on E] (`code/mh_common.py: fluid` → `results/plasma.json`)

**The fluid.** It is my implementation of the mf_hybrid skeptic's per-species fluid (`skeptic/code/sk_pauli.py`). It
starts from the bench state at t = 0.19 and reads the track's exact half-period tables and the skeptic's F_sat(A). I
added one thing: the SM number density. Each N → L H decay gives 2 quanta, so E_SM = ρ_SM/n_SM.

**Validation [C].** It reproduces the skeptic's cases:
- unblocked: K = 2.026 (skeptic 2.025), f_rad = 0.56 / 0.78 / 0.83 / 0.86 at t = 0.3 / 0.6 / 1 / 3;
- Pauli Γ_th = 0: K = 7.054 (7.052);
- Pauli Γ_th = 10: K = 4.099 (4.098);
- tails blocked: K = 20.9 (≈ 22, not converged there).

**The hard-gas mass [E].**
- Formula: m_ind² = c ρ_SM/E_SM².
- Coefficient: c = 0.0228 (inf_higgs/inf_reheat, κ_T-normalised). The inf_higgs skeptic's composition gives 0.031–0.034
  (N11, favourable).
- The SM quanta carry E_SM ≈ 3.9×10³ H_i at birth, from E_N ≈ ¾μ5,max ≈ 7.8×10³ H_i.

| history (HFa) | K | m_ind² ≥ 1 / 10 / 10² / 10³ H_i² at τ = | ρ_SM/ρ ≥ 10⁻⁶ / 10⁻³ / 0.5 | m_ind²/H² at τ = 0.3 / 1 / 10 | pessimistic thermal bath (LPM) |
|---|---|---|---|---|---|
| unblocked (EXP's plasma) | 2.03 | 0.191 / 0.194 / 0.205 / 0.243 | 0.191 / 0.26 / 7.6 | 5×10³ / 4×10⁵ / 10⁸ | τ = 3.1 (T = 8×10¹³ GeV) |
| Pauli Γ_th = 0 | 7.05 | 0.191 / 0.195 / 0.208 / 0.289 | 0.191 / 0.33 / 213 | 1.6×10³ / 4×10⁴ / 3×10⁷ | 5.1 |
| Pauli Γ_th = 10 | 4.10 | same / 0.275 | 0.191 / 0.30 / 155 | 2×10³ / 10⁵ / 7×10⁷ | 3.9 |
| Pauli, tails blocked | 20.9 | same / 0.290 | 0.191 / 0.33 / never | 1.6×10³ / 3×10⁴ / 3×10⁶ | 8.7 |
| unblocked / Pauli, m̃ = 0.015 eV (Σm_ν, skeptic N5) | 2.03 / 7.13 | 0.192 / 0.198–0.200 / 0.22–0.23 / 0.29–0.54 | 0.193 / 0.34–0.61 / 17–241 | 470–1600 / 10⁴–10⁵ | 4.3 / 6.6 |

- HFb is identical to 1%.
- The "plasma ≥ 10⁻³ at τ = 0.64" of mf_hybrid's export is an artefact of the switch at t = 0.6 (the bench's SM column
  is zero before it). The skeptic's 0.31–0.33 is reproduced.
- The first R < 0 window opens at τ = 0.196 (exported) or 0.194 (Pauli), when m_ind² ≈ 10–20 H_i² against ξ|R_min| = 0.41–0.45 H_i².
- **Caveat [E/A].** The fluid starts with the bench's 0.9% radiation already in N₂ at t = 0.19. The real production
  starts during the last χ infall (t ≈ 0.18–0.19). A delay of ~0.01 H⁻¹ changes nothing, because the deadlines are
  ≥ 2.3.

## 4. Higgs modes and D after inflation [C on E] (`code/mh_modes.py`, `code/mh_run.py` → `results/post_D.json`)

**Method.**
- The modes are inf_higgs's equations (canonical normalisation, Hartree, period-mean ⟨R⟩ in the criterion), but my own
  code:
  - a fixed-step vectorised RK4;
  - an exact Bunch–Davies start at t0 = −3, valid because the background is exact de Sitter there;
  - background tables from the bench and the fluid.
- **R(t).** The raw homogeneous R from the bench through t = 0.6 (exported) or 0.19 (fluid histories). After that, the
  fluid mean plus rebuilt lines 9H²f_χ cos(480 t) + 9H²f_σ cos(1200 t), kept to t = 6.
- **D.** higgs_reheat's argued scale-resolved criterion, with the **h̄ grid extended to the abundance limit
  x = 3N + 40** (N1).
- **Numerical checks.**
  - Wronskian error 2×10⁻¹⁵ at the start and 7×10⁻⁴ at t = 80.
  - Controls, as ΔD at 172.86–173.0: Hartree off +2.0; Hartree at 172.92 +1.0; K_max = 15 +0.5; 400 modes ≤ 0.8;
    t0 = −5 −0.08; half step −0.03; R lines off −0.01 (Pauli +0.04); t_end = 300 against 80 ≤ 0.3.
  - So M_t,crit moves by ≤ 0.005 GeV.

| history | K | N_obs | D at 172.28 / 172.57 / 172.86 (grav) | M_t,crit grav (band) [N = 60] | M_h ∓1σ | α_s 0.1171 / 0.1189 | deadline τ at 173.0 / 173.15 / 173.45 / 173.75 | with plasma |
|---|---|---|---|---|---|---|---|---|
| HFa exported (EXP) | 2.05 | 58.05 | −1656 / −299 / **−7.8** | **172.880** (172.826–172.979) [172.885] | 172.840 / 172.987 | **172.562** / 173.246 | 7.0 / 4.8 / 3.0 / 2.3 | no roll ≤ 174.6 |
| HFa unblocked fluid | 2.03 | 58.05 | −1657 / −302 / −7.5 | 172.880 | 172.840 / 172.989 | 172.561 / 173.246 | 6.9 / 4.8 / 3.0 / 2.3 | no roll |
| HFa Pauli Γ_th = 10 | 4.10 | 57.70 | −1657 / −320 / −13.7 | 172.896 (172.854–172.990) | 172.847 / 173.008 | 172.570 / 173.269 | 8.8 / 5.8 / 3.7 / 2.8 | no roll |
| HFa Pauli Γ_th = 0 | 7.05 | 57.43 | −1658 / −349 / **−21.9** | **172.914** (172.885–173.002) [172.937] | 172.856 / 173.026 | **172.601** / 173.290 | 11.5 / 7.3 / 4.5 / 3.3 | no roll |
| HFa Pauli, tails blocked | 20.9 | 56.89 | −1658 / −351 / −22.8 | 172.918 (172.891–173.005) | 172.859 / 173.032 | 172.606 / 173.296 | 11.9 / 7.6 / 4.7 / 3.4 | no roll |
| HFb EXP / Pauli Γ_th = 0 | 2.05 / 7.05 | | as HFa to ±0.1 | 172.880 / 172.914 | | | | no roll |

- **Reading the table.** These are the same two regimes as inf_higgs:
  - fast radiation (⟨R⟩ → 0) is like G05: 172.88;
  - a long χ condensate, ⟨R⟩ = 3H²f_m, is like MD: 172.91–172.92.

  The hybrid's histories sit between them, with a slightly larger N_obs (58 against 54.7–57.1).
- **With the plasma** (the hard-gas mass in m_c², the field grid to the abundance limit):
  - D = −inf ("no roll within x ≤ 3N + 40") at every M_t ∈ [172.86, 174.6] on every history.
  - The same holds at α_s = 0.1171 and 0.1189 (`results/alphas.json`) and at H_i = 0.6–5×10¹¹ (`results/Hi_*.json`).
  - With c × 10⁻⁴, D ≤ −107 at 174.6.
- **Thermal-only reading.** A sharp onset at the pessimistic thermal time, compared with the deadline curve
  (`results/thermal_reach.json`), gives M_t reach 173.44 (exported), 173.39 (Pauli 0), 173.43 (Pauli 10) and
  173.11 (tails blocked).
- **H_i dependence (N6), gravitational only**, fixed dimensionless history (`results/Hi_HFa_*.json`):

| H_i (GeV) | 0.6×10¹¹ | 1.0×10¹¹ | 1.3×10¹¹ | 2.1×10¹¹ | 2.4×10¹¹ | 3×10¹¹ | 5×10¹¹ |
|---|---|---|---|---|---|---|---|
| M_t,crit EXP / Pauli 0 | 173.12 / 173.16 | 172.96 / 172.99 | 172.88 / 172.91 | 172.74 / 172.77 | — / 172.74 | 172.64 / 172.68 | 172.54 / 172.56 |

  - In this convention the plasma rescales as Γ_rest/H_i ∝ H_i and m_ind²/H_i² ∝ (M_P/H_i)². It stays prompt
    (m_ind² ≥ 1 at τ = 0.191–0.193) at every H_i.

## 5. The ξ_H window [C on E] (`results/xiscan_*.json`, `upper_*.json`, `upper_plasma.json`)

**Lower edge after inflation, gravitational only** (D = 0; criterion band in brackets):

| | 172.28 | 172.57 | 172.86 | 173.15 |
|---|---|---|---|---|
| EXP | 0.025 | 0.067 [0.063–0.071] | 0.098 [0.091–0.103] | 0.118 |
| Pauli Γ_th = 0 | 0.025 | 0.066 [0.062–0.069] | 0.095 [0.090–0.098] | 0.114 |
| with plasma | < 0.02 | < 0.02 | < 0.02 | < 0.02 |

**Upper edge.**
- **Mean field, 2σ_prod = h_max (inf_higgs convention; produced variance mapped back adiabatically).**

  | | −1σ | central | +1σ |
  |---|---|---|---|
  | ξ(H_i), EXP | 483 | 134 | 29 |
  | ξ(H_i), Pauli | 492 | 183 | 25 |
  | ξ_EW, EXP | 408 | 113 | 24 |
  | ξ_EW, Pauli | 416 | 155 | 21 |

  For comparison: HT6t 29 (11), T1 70 (37) [Q inf_higgs].
- **Scale-resolved D.** It is ≤ −1000 at every ξ(H_i) ≤ 1000, with no Hartree runaway. The produced variance sits at
  k/a ~ 10²–10³ H_i (n_k ≤ 0.4–0.7), which drives the R lines and the first negative swing.
- **With the plasma.** The barrier of λ_eff φ⁴/4 + ½(m_ind² + ξ⟨R⟩)φ² sits at ≥ 10× the produced 2σ for t ≥ 0.24 at
  every ξ(H_i) ≤ 1000 [E]. Before t ≈ 0.24, production and plasma build-up overlap, and the modes did not include the
  plasma mass (conservative).
- **Comparison.** Figueroa–Rajantie–Torrentí (arXiv:1709.00398 [R]) find a lattice instability for ξ ≳ 4–5 at
  M_t ≈ 173.3 in a different (chaotic, gravitational-only) history; Herranen et al. (1506.04065 [R]) report
  significant production for ξ ≳ 1. The hybrid's ω/H ≈ 500 makes the lines weak, which is why its edge is higher [A].

**ξ_H = 0.1** is inside every version. With the plasma the margin is ≈ ×2 above the inflationary lower edge. Gravitational
only, it is 2% (exported) or 5% (Pauli) above the post-inflation edge at +1σ.

## 6. The G6 portal (hypothetical; the benchmark has none) [C on E] (`code/mh_extra.py g6` → `results/g6_portal_*.json`)

**Setup.**
- The term: λ_χH(χ² − v²)|H|² with λ_χH < 0, giving m_H² = Λ H_i²(1 − χ²/v²), with Λ = |λ_χH|v²/H_i² and
  v = 3.67×10⁵ H_i.
- χ(t): exact from the bench up to the join; after it, χ = v(1 + δ), with 1 − χ²/v² = −2(A/v)cos(m_χ t) + (A/v)² and
  A = √(2ρ_χ)/m_χ.
- The modes carry the instantaneous mass. The criterion uses its mean over exactly two periods. The run is gravitational
  only, with no plasma.

| Λ (λ_χH) | 0.3 (−2×10⁻¹²) | 1 (−7×10⁻¹²) | 1.5 | 3 | 30 | 300 | 3000 (−2×10⁻⁸) | 10⁴ (−7×10⁻⁸) | 3×10⁴ |
|---|---|---|---|---|---|---|---|---|---|
| EXP: D at 172.86 / 173.45 / 174.6 | −219 / +44 / +70 | −1106 / −580 / −274 | −1244 / −640 / −166 | −1844 / −707 / −105 | ≤ −275 | ≤ −532 | ≤ −264 | runaway τ = 0.33 | runaway τ = 0.23 |
| Pauli 0: same | −291 / +33 / +67 | −1056 / −575 / −269 | ≤ −276 | ≤ −326 | ≤ −2590 | ≤ −449 | ≤ −254 | runaway, D ≈ +78 | runaway |
| inflation log₁₀⟨N⟩ at 174.6 (tree) | −230 | −618 | −989 | −2567 | −1.5×10⁵ | −1.3×10⁷ | −1.2×10⁹ | | |

- **Lower edge, Λ ≈ 1.** The Higgs must be heavy enough during inflation (m² ≈ 2–4 H²) to suppress the super-horizon
  variance. The post-inflation mean Λf_χ/4 is only 0.1–0.6 H² at Λ ≤ 3.
- **Upper edge, Λ ≈ 3×10³–10⁴.** It is set by tachyonic and resonant production from the sign-oscillating part, at
  amplitude 2Λ(A/v) ~ 0.8Λ H² at τ ≈ 0.2–0.3. This is the preheating-stage upper bound of the type of Ema–Mukaida–Nakayama
  (arXiv:1602.00483 [R]).
- **Other effects.** The portal's own χ → HH* width is 4.7×10⁻¹⁵Λ² H_i (negligible). Its loop correction to λ_χ is ≪ λ.
- **Caveats.**
  - On the Pauli history the rebuilt oscillation joins the exact χ at t = 0.19 with a phase jump. This affects the
    Λ ≥ 3000 rows (σ² burst at the join), so the upper edge is [E].
  - The plasma would only help.

## 7. What is not established

- **The criteria.** D is argued, not simulated. HM with constant H is the literature's criterion; its use through a
  sub-Hubble waterfall is bracketed, not settled. There is no lattice: the waterfall fragmentation the mf_hybrid
  skeptic found (N2) would make R(t) and the plasma inhomogeneous.
- **The plasma.**
  - The hard-gas mass is [E]: an expectation of h̄-dependent quantum masses in a non-thermal L + H gas. Its sign rests on
    the longitudinal W/Z dominating the negative radial λ (skeptic N11).
  - The result survives a ×10⁻⁴ cut in c. A gas that thermalises only by the pessimistic LPM time gives 173.1–173.4 GeV.
- **The reheating channel** (axial N at μ5/4πf ≈ 0.4–0.5, periodic couplings on non-periodic fields; skeptic N6
  anomalies) is taken as given. If it failed, the gravitational-only rows apply, as do H_i at mf_chain's calibration.
- **Not done.**
  - The plasma mass inside the mode equations (only in the criterion; conservative).
  - α_s for the upper edge.
  - A computed N_s sector; a two-loop calculation of the phase mismatch.
  - The SM above the fit grid (M_t > 174.6 is extrapolation).
  - HFb's G6 and upper-edge runs (its post-waterfall history equals HFa's to 0.1%).
- **Input editions (inf_higgs skeptic N7).** At 172.60 ± 0.27 / 125.13 the σ language shifts by about −0.2σ. It was not
  recomputed.

## 8. Code, results and reproduction (one core; ≈ 45 min in total)

| file | content → results |
|---|---|
| `code/mh_common.py` | bench loader, N_* bookkeeping, post-waterfall fluid with SM number (Pauli or unblocked), background tables, hard-gas plasma |
| `code/mh_modes.py` | Higgs modes (RK4, exact BD start), danger D with the abundance-limit grid |
| `code/mh_inflation.py` | constant-H and instantaneous HM, ξ_min, M_t caps, α_s and H_i → `inflation.json`, `inflation_sens.json` |
| `code/mh_run.py` | bench D, M_t,crit, variants, deadlines (`post_D.json`); plasma (`plasma.json`); ξ scans (`xiscan_*`); upper edge (`upper_*`); H_i (`Hi_*`) |
| `code/mh_extra.py` | controls, Gaussian tail, R windows and lines, G6 portal, plasma barrier, induced masses, thermal reach |
| `code/mh_alphas.py` | α_s for the plasma D and the instantaneous bracket → `alphas.json` |
| `code/mh_summary.py` | → `summary.json` |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_higgs/code; export OMP_NUM_THREADS=1
python3 -B mh_inflation.py main; python3 -B mh_inflation.py sens
python3 -B mh_run.py bench HFa EXP; python3 -B mh_run.py bench HFa PAULI0,UNBL; python3 -B mh_run.py bench HFa PAULI10,PAULIb
python3 -B mh_run.py bench HFb EXP,PAULI0; python3 -B mh_run.py plasma
python3 -B mh_run.py xiscan HFa EXP; python3 -B mh_run.py xiscan HFa PAULI0; python3 -B mh_run.py upper HFa EXP; python3 -B mh_run.py upper HFa PAULI0
python3 -B mh_run.py Hi HFa EXP; python3 -B mh_run.py Hi HFa PAULI0
for c in controls induced tail thermal windows upperpl; do python3 -B mh_extra.py $c; done
python3 -B mh_extra.py g6 EXP; python3 -B mh_extra.py g6 PAULI0; python3 -B mh_alphas.py; python3 -B mh_summary.py
```

Fluids, SM tables and mode runs are cached in
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_higgs/`. These read-only
imports are used:
- tracks/inf_higgs/code (`ih_common.SMTab`, `ih_inflation.barrier_B`);
- tracks/higgs_reheat/code (`smpot`);
- tracks/inf_potential/code (through ih_common);
- tracks/mf_hybrid/results and skeptic/results (bench, axial tables, F_sat).

## Literature (checked by web search in this session unless marked)

- **Inflation and its end.**
  - Mantziris–Markkanen–Rajantie arXiv:2011.03763 (bubble count in the past light cone; decay most likely a few
    e-folds before the end) and arXiv:2207.00696.
  - Markkanen–Nurmi–Rajantie–Stopyra arXiv:1804.02020 (one-loop curved-space SM potential).
  - Herranen–Markkanen–Nurmi–Rajantie arXiv:1407.3141 (ξ_EW ≳ 0.06 stabilises during inflation) and arXiv:1506.04065
    (production for ξ ≳ 1 at the end of inflation).
- **After inflation.**
  - Figueroa–Rajantie–Torrentí arXiv:1709.00398 (lattice: instability for ξ ≳ 4–5 at M_t ≈ 173.3).
  - Ema–Mukaida–Nakayama arXiv:1602.00483 (preheating upper bounds on Higgs–inflaton and Higgs–curvature couplings).
  - Lebedev–Westphal arXiv:1210.6987 (a positive Higgs–inflaton coupling keeps the Higgs heavy).
- **Thermalisation, SM input, seesaw.**
  - Harigaya–Mukaida arXiv:1312.3097 (near-instant thermalisation).
  - Buttazzo et al. arXiv:1307.3536 (SM running).
  - Vissani hep-ph/9709409 (seesaw naturalness).
  - PDG 2024 t-quark listing: 172.57 ± 0.29 GeV (S = 1.5).
- **Not confirmed:** arXiv:1602.07719 (in the brief's list). The search did not find it, and it is not used.
- **Quoted through the tracks** [Q]: Greene–Kofman hep-ph/9807339, Planck 2018 (1807.06211), DESI DR2 (2503.14744).
- **Our tracks, read-only:** inf_higgs (+ skeptic), higgs_reheat, inf_critic, inf_chain (+ skeptic), inf_reheat
  (+ skeptic), mf_hybrid (+ skeptic). mf_chain's interim iteration files were glanced at for H_i only, are unreviewed,
  and are not used for any number here.
