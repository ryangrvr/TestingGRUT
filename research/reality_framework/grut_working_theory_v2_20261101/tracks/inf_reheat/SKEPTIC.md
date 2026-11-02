# SKEPTIC of track inf_reheat (key inf_reheat_skeptic)

3 November 2026 · Claude, adversarial pass on `tracks/inf_reheat` (BRIEF_INFLATON.md; GitHub issue #13 context).
- This is same-session AI work on one supplied model, not outside replication.
- Nothing here is a validation, an exclusion of the real world or a new theory.
- Every choice in the track and here is a model choice, not a consequence of closure.
- Git was read-only and nothing was committed.

Labels: [D] derived, [C] computed with my own code (`skeptic/code/`, output in `skeptic/results/`), [E] estimate,
[A] argued, [Q] quoted (abstract level), [R] recalled.

Units: H_i = 1 at the slice, M_P = 2.435×10¹⁸ GeV (reduced). "b" means m_I/H_i.

**Method policy.** Every recomputation below uses a different method or code path from the track's. The method used is
stated each time. I did not re-run the track's engine (rprh.c). Its N_X, E_X and N_YL are taken as given and are marked
"not rechecked".

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

| Track result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| R1 spurion, mirror μ_eff = 2μ − 10i | **stands** | algebra redone: 14.32 − 10i, \|u\| = 17.47 | D |
| R2 Option B (Aφ\|H\|²) | **stands** | A = 6.140 H_i, m_H² = 9.456×10⁷, growth exponent 7767 (∫√cos = 2.3963), δλ = −2.094 | C |
| R3 Z4 and S1/S2 realisations | **partly refuted** | The Z4 assignment and the Z2-odd failure stand. The S1 claim does not: "M_N soft breaking, induced odd terms ≤ 10⁻¹² of the tree slope, A_ind ~ 10⁻⁹ H_i" fails at the benchmark couplings. The R-odd one-loop slope is 0.05–34 times the tree slope, and A_ind = 0.35–17×10⁻⁶ H_i (§2.2). The lepton-only Z4 is also SU(2)_L-anomalous; it can be repaired (§4, N11) | C/D |
| R4 Majorana width | **stands** | explicit spinor sum gives 13.5 = 2m²β² and Γ/(y²m) = β³/16π exactly | C |
| R5 Greene–Kofman | **stands** | exact 2×2 matrix-exponential propagator: max\|n_k − GK\| = 1.08×10⁻³ (yΦ/m = 20) and 2.2×10⁻⁴ (100); N/N_GK = 1.007 and 1.001 | C |
| R6 K = c_K(b)Γ^{−1/3}, D = K^{3/2} | **stands for the static slice** | my integrator agrees to ≤ 1.4×10⁻¹¹; analytic limit c_K = C x₀^{2/3} with universal C = 0.50128 (to 4×10⁻⁵). **But see R11: the T-model's own rolling slice gives K 30–48% lower** | C |
| R7 yield and H_i scalings | **stands** (arithmetic) | 2.206 Γ^{1/2}, 0.587 Γ^{−1/3} and Γ ≥ 4.2×10⁻⁹ for r < 0.036 reproduced by hand. N_X not rechecked | D |
| R8 N_YL(b), "only m_I ≲ 3 quiet" | not rechecked | engine output; the resonance at b ≈ m_YL is consistent with a linear Y drive through the Y condensate | — |
| R9 holder leakage | not rechecked; **extended** | The induced A_IΣφ\|Σ\|² that the track leaves unevaluated gives R_inf = 0.06–0.09 if the portals vanish at M_Y, and 0.76–1.36 if they vanish at M_P. That is at the track's own R_inf ≤ 1 edge (§4, N8) | E |
| R10 y_nat | stands as a criterion | at the benchmark y, \|V_CW\|/V_tree = 110–260 at φ_*: the one-loop term exceeds the T-model by two orders (§2.2) | C |
| R11 B2 numbers | **refuted in part** | Γ_I/m_I, T_RH and the couplings are not self-consistent with the T-model background the benchmark itself uses. On the rolling slice: Γ/m = 1.30×10⁻⁷ (not 6.2×10⁻⁷), T_RH = 1.23×10¹² GeV (not 2.7×10¹²), y_S1 = 0.0032, y_p = 0.0028, f = 1.5×10¹⁵ GeV, N_MD = 10.4. n_s = 0.9627 (MS, corrected N_*), r = 4.0×10⁻³. H_i and K are essentially unchanged (§2.1, §2.3) | C |
| R12 B1, B3 numbers | **refuted in part** | same correction: B1 Γ/m = 3.8×10⁻⁷, T_RH = 2.1×10¹²; B3 Γ/m = 1.08×10⁻⁷, T_RH = 1.1×10¹² GeV | C |
| R13 plasma times | not rechecked; one input changes | S2's τ ≈ 0.02–0.03 rests on μ5,max = 1.1–1.8 × 4πf, outside its EFT at the track's Γ; on the rolling slice μ5,max/4πf = 0.2–0.6. A new early SM source exists through λ_IH itself (§4, N9) | E |
| R14 λ_IH window | stands as a criterion; **convention factor 2** | (λ/2)φ²\|H\|² gives m_H² = λφ²/2, not λφ². S1's own N loop induces λ_ind ≈ 2×10⁻⁹–10⁻⁷, at or above the window (§4, N5, N10) | D/E |
| R15 leptogenesis | **stands** | κ_f = 5.24×10⁻³, ε₁ = 1.20×10⁻⁵ = 3.0% of Davidson–Ibarra | C |
| R16 S/Q ≈ 190–280 | branch direction stands; **value halves** | with the rolling-slice x₀ in the track's own rescaling, S/Q(10⁻⁶) = 114–129 (≈ 60 with the decisive pass's rolling horizon-exit factor ≈ 0.4); e_crit ≈ 2–5×10⁻⁵ | E |
| R17 ξ_A for S3 | **stands** | evaluated at CMB exit, correctly (6.75 at B2 by hand). Caveat: the ξ ≲ 2.5 bound is for an abelian field; S3 couples gluons | D |
| Headline: no perturbative, radiatively acceptable channel gives Γ/m = 1/6 | **stands** | Every route I tried fails: S1 needs y ≈ 2.9 (CW, blocking); S2 needs f < H_i ≪ Φ (outside the EFT); λ_IH needs q ~ 10⁶; S3 has ξ ≫ 2.5; Option B fails as in R2; large N_f with y ≤ 10⁻⁷ needs N_f ~ 10¹⁵ (species cutoff below H_i) | C/A |
| Headline: "realisable channels give Γ/m ≈ 10⁻⁶" | **misframed** | The 10⁻⁶ is the dark-matter calibration's demand in the T-model stand-in on the static slice, not a property of the channels. It is 1.1–3.8×10⁻⁷ on the rolling slice, and (1+2p)⁻² lower again with a symmetric relic (§4, N7, N13) | A |

---------------------------------------------------------------------------------------------------------------------

## 1. What I recomputed, and how the method differs

| Check | Track's method | Mine | File |
|---|---|---|---|
| n_s, r, A_s | first-order slow roll at φ_* | numerical Mukhanov–Sasaki integration of ζ_k and h_k (Bunch–Davies at k/aH = 300, read at 10⁻³), five k's around k_* | `sk_tmodel.py` |
| N_* | Planck 2018 X eq. 47 constant with −(1/12) ln(g/106.75) and a sudden-decay ρ_th | explicit e-fold bookkeeping from K = a²H (a = 1 at ε_H = 1), with T₀, g_s0 = 3.931, H₀ and k_* = 0.05/Mpc written out | `sk_tmodel.py` |
| K(b, Γ) | static slice: quadratic oscillator released at rest, then a fluid | (i) the same static slice with my own integrator; (ii) the analytic limit C x₀^{2/3}; (iii) **the rolling T-model slice**: the full tanh² oscillator continued from ε_H = 1 with the decay term Γφ̇, then a fluid | `sk_tmodel.py` |
| Greene–Kofman | RK on the two-component equations | exact 2×2 matrix exponentials with piecewise-constant mass, vectorised over k | `sk_channels.py` |
| Majorana width | Dirac trace plus Wick counting | explicit Dirac-basis spinors, summed over spins | `sk_channels.py` |
| S1 radiative terms | Δη criterion for the even part; odd part estimated in words | full (M + yφ)⁴ ln one-loop, with the odd part's slope at φ_* computed numerically; induced A and λ from the loop structure | `sk_channels.py` |
| Leakage | holder Y channels | adds the induced A_IΣ → σ mixing → A_L A_L (equivalence theorem) | `sk_channels.py` |

## 2. Corrections with numbers

### 2.1 The benchmarks used the wrong post-inflation background [C]

- **The inconsistency.** The track sets H_i := H_end at ε_H = 1 of the T-model. It takes Φ₁ from that rolling background, but
  it takes K from the engine's static slice: a quadratic oscillator released **at rest** with ρ = 3M_P²H_i².
- **Why it matters.** A released-at-rest field first sits nearly frozen (w ≈ −1), so ρa³ grows: x₀ = ρa³ in the matter era is
  3.96 / 4.82 / 5.69 at b = 3.21 / 2.48 / 2.17. On the rolling slice the field already carries kinetic energy at ε_H = 1, and
  ρa³ falls through the first zero crossing: x₀ = 2.32 / 2.21 / 2.16.
- **Consequence for K.** Since K → C x₀^{2/3} Γ^{−1/3} (C = 0.50128, a universal constant I computed), K on the rolling slice
  is 0.699 / 0.595 / 0.525 of the track's value at the same Γ.
- **Size.** This background part of the "static slice" systematic is −30% to −48% in H_req at fixed Γ. The track's quoted
  ±20% came from at-rest plateau variants only.
- **Self-consistent shift.** H_req ∝ K/N_X^{2/3} and K ∝ Γ^{−1/3}, so the calibration fixes K and moves Γ by
  (K_roll/K_static)³:

| | B1 (α = 0.3) | B2 (α = 1) | B3 (α = 3) |
|---|---|---|---|
| Γ shift | ×0.341 | ×0.211 | ×0.144 |
| Γ_I/m_I | 3.8×10⁻⁷ (track 1.12×10⁻⁶) | 1.30×10⁻⁷ (6.2×10⁻⁷) | 1.08×10⁻⁷ (7.5×10⁻⁷) |
| T_RH (GeV) | 2.1×10¹² (3.6×10¹²) | 1.23×10¹² (2.7×10¹²) | 1.11×10¹² (2.9×10¹²) |
| K, H_i, r | unchanged (82, 5.1×10¹², 1.3×10⁻³) | unchanged (124, 6.6×10¹², 4.2×10⁻³) | unchanged (136, 7.4×10¹², 1.2×10⁻²) |
| N_MD | 9.4 (8.7) | 10.4 (9.3) | 10.7 (9.4) |
| S1 y (y/y_nat) | 0.0054 (8.3) | 0.0032 (4.4) | 0.0029 (3.6) |
| S2 y_p, f (GeV) | 0.0047, 8.7×10¹⁴ | 0.0028, 1.5×10¹⁵ | 0.0025, 1.6×10¹⁵ |
| S2 μ5,max/4πf | 0.62 (track's Γ: 1.82) | 0.23 (1.09) | 0.20 (1.37) |
| z_RH = M_N/T_RH | 1.9 | 3.3 | 3.6 |

- **Caveat.** The hidden-sector N_X is still the static-slice engine value. A release on the rolling slice (the chain
  track's job) changes N_X too; compensating this K shift would need N_X lower by a factor 0.38–0.58.
- **Verdict.** The quoted Γ_I/m_I, T_RH and couplings of B1–B3 are not self-consistent with the T-model background the
  benchmarks use. H_i, r and K survive.

### 2.2 S1's "soft" R breaking and its induced couplings are not small [C/E]

All numbers are at the track's benchmark y. With m_N = M_N + yφ and a Majorana N (2 dof):
V_CW = −(1/32π²) m_N⁴ [ln(m_N²/μ²) − 3/2].

| | B1 | B2 | B3 | Track's statement |
|---|---|---|---|---|
| \|V_CW,even\|/V_tree at φ_* (μ = yφ_*) | 255 | 114 | 206 | "10⁴ flatness tuning" |
| R-odd slope / tree slope at φ_*, μ = M_N | −34 | −5.2 | −2.8 | "≤ 10⁻¹²" |
| same, μ = yφ_* (most favourable) | 0.57 | 0.087 | 0.046 | "≤ 10⁻¹²" |
| induced A φ\|H\|², y y_ν² M_N/8π² × [1, ln(M_P²/M_N²)] | 0.64–17×10⁻⁶ H_i | 0.36–9.7×10⁻⁶ | 0.35–9.2×10⁻⁶ | "~10⁻⁹ H_i" |
| inflation-branch bound without λ_IH (\|A\|φ_slice < 0.5 H_i²) | 1.55×10⁻⁶ | 1.61×10⁻⁶ | 1.62×10⁻⁶ | — |
| λ_IH,min φ_slice (portal mass that dominates A) | 4.9×10⁻⁵ | 4.6×10⁻⁵ | 4.5×10⁻⁵ | — |
| induced R-even portal λ_ind ≈ y² y_ν²/16π² × [1, ln] | 3.7×10⁻⁹–10⁻⁷ | 2.0×10⁻⁹–5.4×10⁻⁸ | 2.4×10⁻⁹–6.4×10⁻⁸ | window [1.5×10⁻¹⁰, 1.6–2.1×10⁻⁹] |

- **The cubic term.** The odd part is ≈ −M_N y³φ³ ln(...)/8π², suppressed relative to the even part only by
  4M_N/(yφ_*) = 1–1.8×10⁻⁴. The even part is 10² × V_tree, so the odd part is O(10⁻²) of V_tree and O(0.05–30) of its slope.
  NOTES §1 wrote it as "yM_Nφ³", which drops the y² and the loop factor.
- **Consequences.**
  - **Bookkeeping.** In S1 the inflaton potential is not R-even unless the tadpole and cubic counterterms are tuned.
    "α_I ≡ α_I + π" then holds only for (G, c_odd) ≡ (−G, −c_odd): the two sides of V_I inflate differently [D].
  - **Without λ_IH, S1 is an Option-B model.** A_ind exceeds the inflation-branch bound by up to ×10, so sign(A_ind)
    would select the CP sign. **With** λ_IH ≥ λ_min, λφ_slice ≈ 4.5×10⁻⁵ H_i > A_ind, and the bookkeeping is protected [E].
  - **The portal.** S1's own N–L–H loop generates λ_IH at or above the window. The window then needs a tuned counterterm,
    or the q ≤ 350 reading.
  - **The flatness tuning.** It is not a one-parameter 10⁴ tuning. With V_CW ≈ 10² V_tree, the T-model must be the
    loop-corrected potential, which is a functional tuning [A].
- **Verdict.** In the track's three channels, **S2 is the only radiatively consistent one.** S1 at y ≤ y_nat would need
  Γ/m ≲ 7×10⁻⁹. That gives H_req ≈ 1.5–3×10¹³ GeV (rolling to static K), above every r-allowed T-model H_end (≤ 8×10¹² GeV), so
  there is no self-consistent S1 member at all unless p ≫ 1 (N7) [C/E].

### 2.3 CMB numbers: Mukhanov–Sasaki and the e-fold bookkeeping [C]

| | B1 | B2 | B3 |
|---|---|---|---|
| n_s track (slow roll, first order) | 0.96181 | 0.96218 | 0.96263 |
| n_s Mukhanov–Sasaki (at the track's N_*) | 0.96280 | 0.96313 | 0.96349 |
| N_* track → mine (from K) | 53.69 → 53.19 | 54.01 → 53.44 | 54.47 → 53.84 |
| n_s corrected (MS, my N_*) | 0.9624 | 0.9627 | 0.9631 |
| r track → MS | 1.30 → 1.25×10⁻³ | 4.17 → 4.01×10⁻³ | 11.6 → 11.2×10⁻³ |
| H_end track → MS normalisation (GeV) | 5.10 → 5.00×10¹² | 6.56 → 6.43×10¹² | 7.38 → 7.24×10¹² |
| running dn_s/d ln k (MS) | −6.9×10⁻⁴ | −6.8×10⁻⁴ | −6.6×10⁻⁴ |

- **The N_* offset.** I derive the constant explicitly:
  66.894 = −½ ln 3 − ln H₀ + ln T₀ + ⅓ ln g_s0 − ¼ ln(30/π²).
  So the Planck form is 66.89 − ln(k/a₀H₀) + ¼ ln(V²/ρ_end) + (1/12) ln(ρ_th/ρ_end) **− (1/12) ln g_th**. The track wrote
  −(1/12) ln(g/106.75), dropping −0.39 e-folds. A further −0.18 comes from the sudden ρ_th = 3Γ²M_P² against the exact K.
- **Size.** These are small (Δn_s ≈ −4×10⁻⁴ and +1×10⁻³), but they are recomputed corrections to R11/R12's quoted n_s and r.
- **H_end.** It falls by 2% with the MS normalisation, so the self-consistent Γ rises by ×1.06.

## 3. What I confirm (different method, same number)

- **R1 [D].** The holder is the only R-odd term in the hidden sector, so theory(G) = theory(−G). The mirror history gives
  u = 2μ − 10i.
- **R2 [C].** The Option-B numbers are reproduced to all printed digits.
- **R4 [C].** For −(y/2)φN̄N the spinor sum is Σ|ūv|² = 2m_φ²β² (13.5 at m_φ = 3, m = 0.75). With the identical-particle ½,
  Γ = y²mβ³/16π exactly.
- **R5 [C].** Exact matrix-exponential propagation through one crossing gives n_k(0) = 1 to 10⁻¹⁴ and agrees with GK to
  1.1×10⁻³ / 2.2×10⁻⁴.
- **R6 [C].** Static-slice K is reproduced to ≤ 1.4×10⁻¹¹ at b = 2, 2.5, 3 and Γ = 10⁻⁴, 10⁻⁶. The asymptotic law is
  K Γ^{1/3} = 0.50128 x₀^{2/3}.
- **R7 [D].**
  - Y/Y_Q = 3.446×10⁻⁴ × 9316/c_K^{3/2} Γ^{1/2} = 2.206 Γ^{1/2}.
  - H_max(r = 0.036) = πM_P(A_s r/2)^{1/2} = 4.70×10¹³ GeV, so Γ ≥ 4.2×10⁻⁹ H_i.
- **R15 [C].** K_N = 46.3, z_B = 8.2, κ_f = 5.24×10⁻³, n_B/s = 7.25×10⁻⁶ ε₁, ε₁ = 1.20×10⁻⁵ = 3.0% of ε_DI = 4.0×10⁻⁴.
- **R17 [D].** ξ_A = g√(2ε_*)M_P/2 at CMB exit gives 6.75 at B2. The track evaluated it at the right time.
- **The CMB inputs [R/Q].**
  - A_s = 2.10×10⁻⁹: ln(10¹⁰A_s) = 3.044 (Planck 2018, 1807.06209 / 1807.06211).
  - r_0.05 < 0.036 at 95% (BICEP/Keck 2021, 2110.00483).
  - The pivot is 0.05/Mpc.
  - B1–B3 are within 1σ of Planck 2018's n_s = 0.9649 ± 0.0042, and below the r bound.

## 4. New issues

- **N1 [C] Rolling-slice K** (§2.1). The background part of the static-slice systematic is a factor 0.52–0.70 in K. It moves
  the self-consistent Γ by ×0.14–0.34 and T_RH by ×0.38–0.58. It is the largest correction to the benchmark table I found.

- **N2 [C] e-fold bookkeeping** (§2.3). N_* is too large by 0.51–0.63: −0.39 from the dropped −(1/12) ln g_th, and −0.18
  from the sudden ρ_th.

- **N3 [C] Higher-order spectra** (§2.3). MS against first-order slow roll gives n_s +0.0010, r −4% and H_end −2%.

- **N4 [Q/C] Retrospective comparison with post-Planck CMB.** With n_s ≈ 0.9624–0.9631, the pulls are:
  - −0.4 to −0.6σ against Planck 2018 (0.9649 ± 0.0042);
  - −1.8 to −2.0σ against CMB-SPA (SPT-3G D1 + Planck + ACT: 0.9684 ± 0.0030, as quoted in 2510.18656 from 2506.20707);
  - −3.3 to −3.5σ against P-ACT-LB (ACT DR6: 0.9743 ± 0.0034; Louis et al. 2503.14452).

  The long matter era (9–11 e-folds, w = 0) lowers N_* by about 2.5 and n_s by about 0.0017, against the direction the
  newer data prefer. This is a statement about the stand-in family and a retrospective data comparison, not an exclusion.
  The inflaton-potential track should report against these combinations too.

- **N5 [C/E] S1 induced couplings** (§2.2). CW ≈ 10² V_tree; R-odd slope 0.05–34 × tree; A_ind up to 1.7×10⁻⁵ H_i;
  λ_ind ≈ 2×10⁻⁹–10⁻⁷. S1 is not a consistent channel for the T-model stand-in at the benchmark y, and S2 remains.

- **N6 [C/A] S2 EFT and side effects.**
  - **At the track's Γ.** μ5,max/4πf = 1.1–1.8. φ̇_end/(4πf)² = 0.31–0.43. ρ_end^{1/4}/4πf = 0.62–0.73. The phase winds
    2600–4700 times between CMB exit and the end. ξ_ψ = φ̇/2fH = 31–76 at CMB exit and about 3×10³ at ε_H = 1. Any UV
    completion with a radial mode m_ρ ≪ 4πf is excited [A].
  - **On the rolling slice** (N1), f rises to 0.9–1.6×10¹⁵ GeV and μ5,max/4πf falls to 0.2–0.6, inside the EFT.
  - **Gauge production [A].** In the basis where the lepton doublets carry the phase, there is an explicit (φ/f_A)WW̃
    (1/f_A = 3α₂/4πf). Its ξ_W is 0.2–0.45 at CMB exit (the track's "≈ 0.3") but 15–20 at the end of inflation. In the
    basis without the rotation, the coupling to the SM vanishes as y_ν → 0. Physical W production therefore needs the
    y_ν-mediated chirality flips and is suppressed. Not computed.
  - **Spontaneous leptogenesis [A], not considered by the track.** After rotation the axial coupling is (∂φ/2f)J_L, a
    lepton chemical potential φ̇/2f ~ 10³ H_i during the oscillation. With ΔL = 2 processes active (T_RH ~ M_N), a
    Cohen–Kaplan-type asymmetry can form. Its sign is set by sign(f) and the oscillation phase, not by α_I. Its size
    (suppressed by Γ_ΔL/m_I) is not estimated reliably here.

- **N7 [D/E] The symmetric-relic feedback is misplaced.** The track says p ≳ 1 at S/Q ~ 10² would lower H_i by
  (1+2p)^{−2/3}. In the self-consistent T-model, though, H_i = H_end is fixed by the inflaton. The dark-matter calibration
  therefore moves **Γ** instead, by (1+2p)⁻² (D ∝ 1 + 2p at fixed H_i). At B2:

| p | 0.1 | 1 | 10 | 30 |
|---|---|---|---|---|
| Γ/m | 4.3×10⁻⁷ | 6.9×10⁻⁸ | 1.4×10⁻⁹ | 1.7×10⁻¹⁰ |
| T_RH (GeV) | 2.2×10¹² | 8.9×10¹¹ | 1.3×10¹¹ | 4.4×10¹⁰ |
| K | 140 | 258 | 943 | 1920 |
| S1 y / y_nat | 7.9 | 3.2 | 0.45 | 0.16 |
| z_RH | 1.8 | 4.6 | 32 | 93 |
| S/Q(10⁻⁶), track's rescaling | 213 | 268 | 436 | 569 |

  (Static-slice K; multiply Γ by a further 0.21 for N1.)
  - **No runaway [A].** With p ∝ S/Q ∝ K^{3/8}, the loop has a fixed point (p ∝ (1+2p)^{1/4}).
  - **The real gap.** p(S/Q ≈ 10²) is unknown: the map ends at 6.7. With a hidden sector whose pairs can only leave
    through number-changing processes, O(1–10) is plausible [A].
  - **Consequence.** The benchmark Γ, T_RH, leptogenesis regime (thermal at p ≲ 1, lost at p ≳ 10) and S1 naturalness are
    hostage to p.

- **N8 [E] Induced A_IΣ leakage.**
  - **Size.** A natural-size A_IΣ ≈ 2|Re(Gμ*)|/16π² × ln: 0.9–1.1×10⁻⁵ H_i if the portals vanish at M_Y, 3.9×10⁻⁵ at M_P.
    At the new g_I = 3.8–4.0×10⁻⁵, it mixes φ with σ (m_σ = 40 H_i).
  - **Rate.** Γ(φ → A_L A_L) = A²m_I³/(32π m_σ⁴) = 5×10⁻¹⁸–2×10⁻¹⁶ H_i.
  - **Heat.** R_inf = 0.06–0.09 (M_Y) or 0.76–1.36 (M_P), against the track's requirement R_inf ≤ 1 (Γ_hid ≤ 7×10⁻¹⁷).
  - **Weight.** Subdominant to R_vec ≈ 10², but no longer negligible, and it scales as g_I² ∝ (M_P/H_i)⁻².

- **N9 [E] λ_IH is itself an early Higgs-coupled plasma source.**
  - **Channel.** The oscillating portal mass drives condensate annihilation φφ → HH*. It is open for q < 1/2, i.e. for the
    lower half of the window.
  - **Rate.** Γ_eff = λ²Φ₁²β/(16π m_I) = 2–30×10⁻¹² H_i, 10⁻⁶–10⁻⁵ of Γ_I.
  - **Effect.** It injects 0.06–3 H_i⁴ of SM energy in the first Hubble time: T ≈ 0.2–0.5 H_i and a hard-gas
    m_ind² ≈ 2×10⁻⁴–5×10⁻³ H_i². It gives an early Higgs-coupled plasma independent of the N channel, in the helpful
    direction. There is no narrow resonance (q²m/H ≈ 0.02).

- **N10 [D] λ_IH convention.** For a canonical doublet, (λ/2)φ²|H|² gives m_H² = λφ²/2. The quoted window
  [1.5×10⁻¹⁰, 2×10⁻⁹] is therefore the window for λφ²|H|²; for (λ/2)φ²|H|² it is twice as large. Literature comparisons
  (Enqvist et al., λ_nat) must use one convention.

- **N11 [D/A] The Z4 as stated is anomalous.**
  - **The anomaly.** With L → iL, e^c → −ie^c and quarks neutral, the SU(2)_L 't Hooft vertex carries Z4 charge 3 (mod 4).
    So R is broken by electroweak instantons and is not a discrete gauge symmetry.
  - **The repair.** Q → iQ, u^c → −iu^c, d^c → −id^c (a Z4 of 3B + L) keeps every SM Yukawa. It makes the SU(2)² vertex
    charge 12 ≡ 0, the SU(3)² sum 0, and the gravitational sum 0.
  - **Why it matters here [A].** φ_* = 4–8 M_P, so any R-odd Planck-suppressed operator (φ^{2n+1}|H|²/M_P^{2n−1}, …) is
    unsuppressed during inflation. The CP bookkeeping needs R respected by gravity, i.e. a gauged, anomaly-free version.
    The same concern applies to the plateau and the holder sequestering.

- **N12 [E] Pauli blocking of φ → NN.** The N thermalise through y_ν (rate ≈ 10⁻⁴ H_i ≫ H near reheating). The decay rate
  then carries 1 − 2n_F(m_I/2) = tanh(m_I/4T): 0.81–0.91 at T_RH, 0.65–0.78 one e-fold earlier, 0.49–0.62 two e-folds
  earlier. This is a thermal correction the track omitted. It needs about a 1.2–1.3× larger vacuum Γ for the same K, and
  applies to S1 and S2.

- **N13 [A] Framing of the headline.** "Realisable channels give Γ_I/m_I ≈ 10⁻⁶" should read: the dark-matter calibration
  inside the T-model stand-in demands Γ_I/m_I ≈ 10⁻⁶ on the static slice (1–4×10⁻⁷ on the rolling slice, lower still with
  p). S2 can deliver that Γ; S1 cannot do so naturally (N5).

- **N14 [E] Branch, re-estimated.**
  - **Rolling slice.** The matter-era vector crossing scales with x₀, so R_vec falls by 0.38–0.58 and S/Q(10⁻⁶) = 114–129.
    With the decisive pass's rolling horizon-exit factor (~0.4) it is ≈ 60.
  - **e_crit.** Since R_vec ∝ e^{−3/2}, e_crit ≈ 4–5×10⁻⁵, or ≈ 2–3×10⁻⁵ with that factor.
  - **The branch.** Free X, far from the 3.47 threshold, stands in direction. The magnitude is uncertain by about 3×, and
    p(S/Q) is uncomputed there.

## 5. Omitted-effects checklist (from the brief)

| Effect | Status |
|---|---|
| Preheating | GK (S1) confirmed; S2 Landau–Zener outside its EFT at the track's Γ (N6); λ_IH q < 1, no resonance (N9); no inflaton self-resonance for α ≥ 0.3 [Q, Lozanov–Amin] |
| Fragmentation, oscillons | not expected for √(6α) ≥ 1.3 M_P [Q]; not recomputed |
| Thermal corrections | Pauli blocking (N12, ×0.8–0.9 near T_RH). Inflaton thermal mass y²T²/12 ≈ 10⁻⁴ H_i² at T_max ≈ 5 H_i: negligible. N thermal mass ≈ 0.15 H_i: decay stays open |
| Induced couplings | S1: CW, R-odd cubic and tadpole, A φ\|H\|², λ φ²\|H\|² (N5). Holder: A_IΣ (N8). The holder's loop correction to m_I² rises ×2.3×10³ with g_I but stays ≤ 2×10⁻⁶ H_i² [E, scaled from fx_holder] |
| Kinematic and perturbativity conditions | S1 blocking reproduced (yΦ₁/(m/2) = 660–860). S2 at Γ/m = 1/6 needs y_p ≈ 1.7–2.9, i.e. f ≲ H_i ≪ Φ ≈ 10⁷ H_i, outside its EFT |
| CP-sign bookkeeping | S2: preserved (no potential is generated, and sign(f) does not enter the hidden sector). S1: needs λ_IH and tuned odd counterterms (§2.2). Z4 needs the quark extension (N11) |
| M1r assumptions | static-slice N_X (not rechecked); g_I ×48; induced-portal renormalisation condition now matters for R_inf (N8) |

## 6. What this skeptic did not establish

- **The engine.** I did not recompute N_X, E_X, N_YL or the b-scan (R8). The rolling-slice correction applies to the
  background only.
- **The rescaled pieces.** S/Q, the vector crossing and R_inf use the track's rescaling [E]. I did not run dq_heat's mode
  code on any slow background.
- **The Higgs at the new H_i.** Not redone; it is the track's open item.
- **Spontaneous leptogenesis and W production in S2.** Argued, not computed.
- **The brief's Peloso–Sorbo reference correction.** Not checked.

## 7. Files and reproduction

| File | Content |
|---|---|
| `skeptic/code/sk_tmodel.py` → `skeptic/results/sk_tmodel.json` | static-K check and universal C; Mukhanov–Sasaki n_s, r, running and normalisation; N_* from K; Φ₁; rolling-slice K and the self-consistent Γ shift |
| `skeptic/code/sk_channels.py` → `skeptic/results/sk_channels.json` | GK by matrix exponentials; spinor-sum width; Option B; S1 odd and even one-loop terms, A_ind and λ_ind; portal annihilation; S2 EFT ratios; CMB pulls; p feedback; S/Q on the rolling slice; A_IΣ leakage; Pauli blocking; rolling self-consistent table |
| `skeptic/code/sk_summary.py` → `skeptic/results/sk_summary.json` | headline numbers |

Reproduction, in `skeptic/code/`, one thread: `python3 -B sk_tmodel.py` (about 2 min), then `python3 -B sk_channels.py` (about
1 min), then `python3 -B sk_summary.py`. The scripts read the track's `results/benchmarks.json` and `results/fluid.json`
read-only.

**Literature** (abstract level unless stated):
- Planck 2018 VI/X: 1807.06209, 1807.06211 [R].
- BICEP/Keck 2021: 2110.00483 [R].
- ACT DR6: Louis et al. 2503.14452 (P-ACT-LB n_s = 0.9743 ± 0.0034) [Q, via secondary sources].
- SPT-3G D1: Camphuis et al. 2506.20707 (SPT alone n_s = 0.951 ± 0.011) [Q]; the CMB-SPA 0.9684 ± 0.0030 is as quoted
  in 2510.18656 [Q].
- Greene–Kofman hep-ph/9807339 [R].
- Graham–Mardon–Rajendran 1504.02102 [R].
- Cohen–Kaplan spontaneous baryogenesis (Phys. Lett. B 199 (1987) 251) [R].
- Ibáñez–Ross discrete-anomaly conditions (Phys. Lett. B 260 (1991) 291) [R].
- Buchmüller–Di Bari–Plümacher hep-ph/0401240 [R].
- Davidson–Ibarra hep-ph/0202239 [R].
