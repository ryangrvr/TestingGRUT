# Track inf_reheat: realisable reheating channels from the inflaton to the SM for M1r

3 November 2026 · Claude, inflaton-sector pass (BRIEF_INFLATON.md; GitHub issue #13 context). This is same-session AI
work on one supplied model, not outside replication. Nothing here is a validation, an exclusion of the real world or
a new theory. **Every choice below is a model choice, not a consequence of closure.** That covers the channel, the
couplings, the symmetry assignments, the stand-in inflaton families, M_N, m̃ and λ_IH. Git was read-only and nothing
was committed.

Labels: [D] derived, [C] computed with our own code (`code/`, output in `results/`), [E] estimate or criterion choice,
[A] argued, [Q] quoted (abstract level unless stated), [R] recalled.

Units: H_i is the Hubble rate at the M1r slice, and τ = H_i t. M_P = 2.435×10¹⁸ GeV (reduced). "b" means m_I/H_i.

**What I did not have.** The sibling inflaton-potential track had produced no benchmarks when this track ran.
- I therefore pair the channels with stand-ins:
  - F0, the engine's phenomenological oscillator;
  - the α-attractor T-model V₀ tanh²(φ/√(6α)M_P), integrated exactly to ε_H = 1 (the slice: H_i := H_end).
- The engine runs use the static slice M1r specifies (released at rest, R = 12 H²). The real rolling slice is the chain
  track's job.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer on one screen

| # | Result | Label |
|---|---|---|
| 1 | **Symmetry.** The holder −(g_I e^{iα_I} φ_I Σ*Y² + h.c.) is the only R-odd term (R: φ_I → −φ_I). With an R-even inflaton sector, R is a spurion relation G ↔ −G. Physics depends on Gφ_I, so the CP sign is fixed by α_I together with sign(φ_I) during inflation. The mirror history (φ_I → −φ_I) is not the CP image: it gives μ_eff = 2μ − 10i = 14.32 − 10i H_i instead of 10i, so the opposite sign and a different magnitude. **Recommended: keep R.** Breaking it with Aφ_I\|H\|² is untenable at decay strength (row 2) | D |
| 2 | **Option B (break R, Aφ_I\|H\|²).** Γ/m = 1/6 at m_I = 3 needs A = 6.14 H_i = 7.9×10¹¹ GeV. Then: Higgs mass² ±Aφ_I = ±9.5×10⁷ H_i² at the slice (tachyonic on one branch); a tachyonic growth exponent of 7.8×10³ per negative half-cycle; a tree quartic shift −A²/2m_I² = −2.09; one-loop Higgs correction to V_I of 5% of V, with a slope 7×10² times the tree slope. A safe A is ≤ 4 TeV (inflation-branch stability) or ≤ 75 TeV (no tachyonic resonance). That gives Γ_hh ≤ 4.5×10⁻¹⁵ H_i: it cannot reheat, but it would give a marginal Higgs plasma (T ≈ 0.26 H_i at τ = 8 for instant thermalisation). Then sign(A) selects the surviving CP sign | C/E |
| 3 | **Keeping R needs N → ±iγ5N (a Z4 of lepton type: L → iL, e^c → −ie^c).** It forbids a bare Majorana M_N, the Weinberg operator and φ_I\|H\|², and allows yφ_I N̄N, y_ν L̄H̃N and the R-even portal (λ_IH/2)φ_I²\|H\|². Variants: **S1** = scalar y with M_N as the only (soft) R breaking (type-I seesaw; induced odd terms negligible); **S2** = axial M_N N̄e^{iγ5φ_I/f}N (R includes parity on N; \|m_N\| independent of φ_I); a Z2-odd partner N′ gives a stable N′ that overcloses (Ωh² ≈ 4×10⁷ M′/GeV) or, if lighter than 3 eV, gives ΔN_eff ≈ 2.8 | D/E |
| 4 | **Width convention.** With −(y/2)φN̄N (so m_N = M + yφ): Γ = y² m β³/16π. The brief's y²m/8π is the Dirac normalisation (trace computed). Γ/m = 1/6 needs N_f y² = 8.38 (Majorana) or 4.19 (Dirac) | D/C |
| 5 | **Blocking at the engine amplitude** (Φ = 0.816 M_P, m_I = 3). The (3, 0.5) coupling y = 2.89 has yΦ/(m/2) = 3.0×10⁷, an open fraction of 2×10⁻⁸ and a Greene–Kofman energy transfer of 5×10² ρ_φ per crossing: fully non-perturbative. **Greene–Kofman checked numerically**: n_k after one crossing equals exp(−πk²/yΦm) to ≤ 1.1×10⁻³, number and energy within 2%, and the next crossing stays Pauli-bounded (n_k ≤ 0.70) | C |
| 6 | **Radiative flatness.** The fermion one-loop potential −y⁴φ⁴ln φ/32π² must not spoil n_s: \|Δη\| ≤ 0.002 with one tuned counterterm gives y_nat = (0.008π² V_*/φ_*²)^{1/4} ≈ 4–9×10⁻⁴ in the T-model family. **The (3, 0.5) history (y ≈ 3) is out by 10⁴ in y.** S2 has no such potential | E |
| 7 | **Slow-reheating background** (engine-validated to ≤ 1.6×10⁻³). K = c_K(b)(Γ_I/H_i)^{−1/3}, with c_K = 2.721, 1.712, 1.417, **1.284**, 1.170, 1.121 for b = 1.5, 2, 2.5, 3, 4, 5. The dilution D = K^{3/2}. The matter era lasts N_eq = (2/3)ln(H_i/Γ_I) + 0.35 e-folds (b = 3) | C |
| 8 | **Engine (rprh.c, the copy of rpdq.c), slow decay.** At b = 3 the charge is \|N_X\| = 9316 (×1.158 of (3, 0.5)) with E_X ×1.098. **Y/Y_Q(H_M1) = 2.206 (Γ_I/H_i)^{1/2}**, so **H_i = 1.291×10¹¹ GeV × 0.587 (Γ_I/H_i)^{−1/3}**. r < 0.036 (with H_CMB = H_i) needs Γ_I ≥ 4.2×10⁻⁹ H_i, i.e. y ≥ 3.3×10⁻⁴. **New:** with the trilinear holder, the oscillating μ_eff drives Y_L: N_YL = 0.57 (b = 3) → 2.2 (4.5) → 41 (7.2) → 589 (11.4) → 5.5×10³ (22 ≈ m_YL) → 102 (33) → 57 (102). Only **m_I ≲ 3 H_i** is quiet | C |
| 9 | **Holder leakage into the hidden sector.** Γ(φ_I → hidden) is 10⁻⁴⁰–10⁻³² H_i below m_YL = 22.25 H_i (six-body), 10⁻²⁵–10⁻²⁰ for 25–40 H_i (Y_L + 3X off-shell), and 2–3×10⁻¹⁰ above 2m_YL = 44.5 H_i. Heat R_inf = 3M_P²K² BR/E_X (= 5.7×10⁹ K² BR at H_i = 1.29×10¹¹). Small-α attractors (b ≈ 100) give R_inf ≈ 64 even with fast decay | C/E |
| 10 | **Self-consistent benchmarks** (T-model, DM calibration H_req = H_end; table in §4). Quiet members (α = 0.3–10, b = 3.2–2.0) need **Γ_I/m_I ≈ 0.6–1.1×10⁻⁶**, giving H_i = 5.1–7.8×10¹² GeV, **r = 1.3×10⁻³–3.2×10⁻²**, n_s = 0.962–0.963, T_RH ≈ 3×10¹² GeV, **K = 82–138**, a matter era of 8.7–9.4 e-folds, and y = 0.007–0.009 for S1 (y/y_nat ≈ 9–14, a 10⁴ flatness tuning) or an axial y_p = M_N/f = 0.006–0.008 for S2 (f ≈ 5–7×10¹⁴ GeV). The S3 gauge channel φGG̃ fails the inflationary gauge-production bound (ξ_A = 5.0–23 > 2.5) | C/E |
| 11 | **Higgs-coupled plasma.** S1 is perturbatively blocked until τ ≈ 440–570, but instant preheating (heavy Greene–Kofman N → LH) gives an induced Higgs mass² ≥ 10⁻³ H_i² by τ ≈ 0.3–0.8. S2 reaches it by τ ≈ 0.02–0.03. A fully thermal bath (pessimistic LPM time) comes only at τ ≈ 2×10³. **The deadline becomes moot with the R-even portal λ_IH**: 1.5×10⁻¹⁰ ≤ λ_IH ≤ 1.6–2×10⁻⁹ (no broad resonance; up to ~6×10⁻⁷ by a transfer of Enqvist et al.). This keeps λ⟨φ_I²⟩ ≥ 10H² through the whole matter era and makes the Higgs heavy during inflation; it is radiatively natural (λ_nat ≈ 3–6×10⁻⁷) | E |
| 12 | **Consequence for the branch [E, rescaled budget].** K ≈ 100 and the vector relic set in the matter era (Γ_I < m_A) raise R_vec(10⁻⁶) from 3.6 to ≈ 7×10², so **S/Q(e = 10⁻⁶) ≈ 2–3×10²**: deep free X, outside the computed p(S/Q) map (which ends at 6.7). e_crit (bulk) ≈ 7–8×10⁻⁵. Leptogenesis (T_RH ≈ M_N, strong washout, K = 46) needs ε₁ ≈ 1.2×10⁻⁵ = 3% of the Davidson–Ibarra bound. It is a by-product with its own CP phase, unrelated to α_I; a near-maximal phase overproduces n_B by ≈ 30 | E |

**Bottom line.** A realisable, perturbative, CP-bookkeeping-preserving reheating channel exists:
- keep R; reheat through singlet N (S1 tuned, or S2 axial and radiatively flat), decaying N → LH; add the R-even Higgs
  portal λ_IH ≈ 10⁻¹⁰–10⁻⁹;
- but it cannot be fast. Γ_I/m_I ~ 10⁻⁶ replaces the supplied 1/6.

That changes the M1r numbers it was built on, in these model-choice stand-ins:
- the calibrated H_i rises ~50× (to ~6×10¹² GeV) and r to 10⁻³–10⁻²;
- the hidden plasma is ~40× hotter (S/Q ~ 10²);
- only m_I ≲ 3 H_i keeps the holder quiet.

---------------------------------------------------------------------------------------------------------------------

## 1. Item 1: the reflection and the CP-sign bookkeeping [D]

**The spurion statement.** Write the holder as u Σ*Y² + h.c. with u = μ + Gφ̂_I, G = g_I e^{iα_I}.
- μ and κ are made real by rephasing X, Y and Σ (two non-gauge phases remove both).
- The late hidden sector is therefore CP-conserving. The CP-odd content is arg(u_slice/μ_late), and with M1r's numbers
  u_slice = 10i H_i.
- Suppose every other term is even under R: φ_I → −φ_I (V_I(φ) even, decay couplings even). Then the field
  redefinition φ_I → −φ_I maps the theory with G onto the theory with −G. So α_I ≡ α_I + π once sign(φ_I) during
  inflation is counted as part of the input.
- The opposite-side history (inflation at −φ_I) gives u = μ − Gφ̂ = 2μ − 10i = 14.32 − 10i. That is Im < 0, so the
  **opposite charge sign**, but with |u| = 17.4 and a different held phase, so **not the same magnitude** [D].
- Which side our patch inflated on is an initial condition. For a hilltop with ⟨φ_I⟩ = ±v_I, it is which vacuum our
  patch rolled to (the walls lie outside the horizon) [A].

**Radiative closure of the R assignment [A].**
- Inside the hidden sector, G is the only R-odd coupling, so every induced R-odd inflaton term (tadpole, φ³,
  A_IΣ φ|Σ|², A_IY φ|Y|²) is ∝ Re(Gμ*). These are the "induced portals" whose renormalisation condition M1r supplies.
- φ_I|H|² is not generated at any order by the hidden sector, because M1r has no SM portal. A = 0 is technically
  natural.

**Option A: keep R.** A singlet inflaton has no renormalizable Yukawa to chiral SM fermions, so the decay goes to
singlets.
- **N → ±iγ5N.** (In two-component form ν → ∓iν. These are the only phase actions compatible with the Majorana
  condition.) Then N̄N and N̄iγ5N are both R-odd, so yφN̄N and y′φN̄iγ5N are allowed and a bare M_N N̄N is forbidden.
- The seesaw Yukawa y_ν L̄H̃N is R-even only if L → iL and e^c → −ie^c. R is then a **Z4 subgroup of lepton number**. It
  also forbids the Weinberg operator (LH)² unless it is multiplied by φ_I.
- **Consequences.**
  - With ⟨φ_I⟩ = 0 (any even plateau or oscillator) the N are massless in vacuum: no seesaw, Dirac-like neutrinos.
    Reheating products would be N's coupled only through y_ν ~ 10⁻¹³: dark radiation, which is untenable [A].
  - Three realisations avoid this:
    - **S1, soft breaking by M_N.** The R-odd mass is dimension 3. The induced odd terms are ~yM_N³φ/16π² (a tadpole)
      and yM_Nφ³: ≤ 10⁻¹² of the tree slope [E], and A_ind ~ y y_ν² M_N/16π² ~ 10⁻⁹ H_i for φ|H|². So "α_I ≡ α_I + π" holds
      to that accuracy. Type-I seesaw works with m̃ = y_ν²v²/M_N.
    - **S2, axial.** −(M_N/2)N̄e^{iγ5φ_I/f}N with R = φ_I → −φ_I combined with parity on N. |m_N| = M_N for every φ_I,
      so the N loop gives no φ_I potential and there is no kinematic blocking. R is broken only by the chiral SM
      couplings, which after rephasing leave derivative/anomaly couplings (φ/f)WW̃ ∝ α₂/8π (ξ_A ≈ 0.3 there) [A].
    - **Spontaneous breaking, ⟨φ_I⟩ = ±v_I (hilltop).** M_N = yv_I. The holder then needs φ_I(slice) − v_I, because
      μ_late = μ_bare + Gv_I is rephased real. But hilltops have m_I ≫ 44 H_i and fail on items 8–9.
  - **Z2-odd partner N′** (yφ N N′ with N even and N′ odd). N′ is stable and every decay makes one:
    - Y_N′ ≈ (3/4)T_RH/m_I ≈ 0.14 at the benchmark T_RH, so Ωh² ≈ 2.75×10⁸ (M′/GeV) Y = 4×10⁷ (M′/GeV);
    - for M′ < 3 eV it is instead dark radiation with half the energy: ΔN_eff ≈ (106.75/1.75)(10.75/106.75)^{4/3} = 2.8.
    - Untenable unless R is softly broken so that N′ decays [E].

**Option B: break R by Aφ_I|H|²** (`code/channels.py: optionB`, at m_I = 3 H_i, Φ = φ_slice = 0.816 M_P) [C/E].
- The A needed for Γ_I/m_I = 1/6 (Γ(φ → HH*) = A²/8πm_I, four real components) is A = 6.14 H_i = 7.9×10¹¹ GeV.
- Consequences at that A:
  - Higgs mass² ±Aφ_I = ±9.5×10⁷ H_i² at the slice (m_H = 1.3×10¹⁵ GeV), tachyonic on the branch Aφ_I < 0.
  - Even on the stable branch, each negative half-oscillation gives ∫√(Aφ)dt = 7.8×10³: explosive tachyonic Higgs
    growth into λ < 0.
  - Integrating out φ_I shifts the low-energy quartic by −A²/2m_I² = −2.09.
  - One-loop Higgs correction to V_I: 5% of V, with a slope 706× the tree slope at r ~ 3×10⁻⁷. This destroys the
    plateau.
- **Stability bounds.**
  - |A|φ_slice < 0.5 H_i², keeping half of ξ_H R on the negative branch: A ≤ 3.2×10⁻⁸ H_i = 4.2 TeV.
  - q_A = |A|Φ/m_I² < 1, no tachyonic resonance (Enqvist et al. 1608.08848 find |σ_hφ| ≲ 10⁸ GeV, i.e. q ~ 1, in
    their chaotic model): A ≤ 5.8×10⁻⁷ H_i = 75 TeV.
- **What such an A can do.** Γ_hh ≤ 4.5×10⁻¹⁵ H_i, so it cannot reheat. It gives T_SM(τ = 8) = 0.26 H_i (instant
  thermalisation) or 0.06 H_i (at the inflation-branch bound). That is a marginal plasma.
- **CP bookkeeping in Option B.** Patches with Aφ_I < 0 during inflation destabilise the Higgs once |A|φ_* > 0.5 H². The
  surviving asymmetry sign is then α_I with sign(φ_I) = sign(A). The CP sign becomes (α_I, sign A) [D, given the
  stability criterion].

**Decision (model choice).** Keep R (Option A). Use S1 or S2 for the decay, and add the R-even (λ_IH/2)φ_I²|H|² for
the Higgs (§2.5).

## 2. Item 2: the singlet-fermion channel

### 2.1 Width and kinematics [D/C]

- `channels.py: widths` computes Σ|ū(p₁)v(p₂)|² = Tr[(p̸₁ + m)(p̸₂ − m)] = 2m_φ²β² with explicit Dirac matrices
  (16.72 = 16.72), so Γ_Dirac = y²mβ³/8π.
- Majorana with −(y/2)φN̄N: two Wick contractions × ½ give the vertex y, and identical particles give ½, so
  Γ = y²mβ³/16π. Writing −yφN̄N instead gives y²mβ³/4π. The pseudoscalar (S2) has β¹ instead of β³.
- Γ/m = 1/6 needs N_f y² = 8.38 (Majorana, adopted), 4.19 (Dirac; WORKING_THEORY's "≈ 4") or 2.09 (−yφN̄N).
- Blocking needs |m_N(φ)| = |M_N + yφ| > m_I/2. The open fraction is (2/π)arcsin(m_I/2yΦ). In the matter era
  Φ ∝ a^{−3/2}, so it unblocks at a_unb = (2yΦ_1/m_I)^{2/3}.

| y (engine amplitude, m_I = 3) | Γ₀/m_I | yΦ/(m_I/2) | open fraction | a_unb | GK energy fraction per crossing |
|---|---|---|---|---|---|
| 2.89 (Γ/m = 1/6) | 0.167 | 3.0×10⁷ | 2.1×10⁻⁸ | 9.6×10⁴ | 5.2×10² (fully non-perturbative) |
| 1 | 0.020 | 1.0×10⁷ | 6.2×10⁻⁸ | 4.7×10⁴ | 37 |
| 10⁻² | 2.0×10⁻⁶ | 1.0×10⁵ | 6.2×10⁻⁶ | 2.2×10³ | 3.7×10⁻⁴ |
| 10⁻³ | 2.0×10⁻⁸ | 1.0×10⁴ | 6.2×10⁻⁵ | 4.7×10² | 1.2×10⁻⁶ |
| 10⁻⁴ | 2.0×10⁻¹⁰ | 1.0×10³ | 6.2×10⁻⁴ | 1.0×10² | 3.7×10⁻⁹ |

**Non-perturbative production** (Greene–Kofman hep-ph/9807339, hep-ph/0003018; Peloso–Sorbo hep-ph/0003045) [C].
- Per zero crossing, n_k = exp(−πk²/yΦm_I), so n = (yΦm_I)^{3/2}/4π³ Majorana N, and the energy at the next turning
  point is ≈ n yΦ. The energy fraction is f₁ = y^{5/2}(Φ/m_I)^{1/2}/2π³.
- `code/fermion_modes.py` solves the two-component mode equations through two crossings at yΦ/m = 20, 100, 500:
  - |n_k − GK| ≤ 1.1×10⁻³, 2.2×10⁻⁴, 4.3×10⁻⁵;
  - N/N_GK = 1.007, 1.001, 1.000;
  - E/(N yΦ) = 1.019, 1.003, 1.000;
  - after the second crossing n_k ≤ 0.64, 0.70, 0.52 (Pauli: N₂/N₁ = 1.49, 1.39, 1.13).
- Production saturates unless the N decay between crossings.

**Instant preheating (S1)** [E]. The GK-produced N are heavy at the turning points (m_N ~ yΦ ~ 10³ H_i in the
benchmarks). They decay at y_ν²m_N/8π ≈ 0.1–0.3 H_i into LH (Felder–Kofman–Linde hep-ph/9812289), so a
Higgs-coupled population appears before perturbative unblocking (§2.4).

**Axial production (S2)** [E].
- The mass phase rotates at φ̇/f. In the rotated basis this is an axial chemical potential μ5 = φ̇/2f, and each helicity
  has k_eff = hk + μ5.
- Landau–Zener at k_eff = 0 with gap M_N, twice per half-period, gives n_k ≈ 2πM_N²/|μ̇5| ≈ 6×10⁻⁴ (helicity
  suppressed; cf. Adshead–Sfakianakis 1508.00891) for k < μ5,max ≈ 1.4–2.2×10³ H_i.
- The resulting Γ_eff ≈ 10³ × the perturbative width initially, but Γ_eff/H ≈ 10⁻³ stays constant. So ≈ 1% of the
  energy moves over the matter era. It provides an early plasma without changing D or K.
- μ5,max is ~ the EFT cutoff 4πf. This is strained; a mode calculation is needed.

### 2.2 Radiative flatness [E]

- With m_N = yφ, V_CW = −(N_f y⁴φ⁴/32π²)[ln(y²φ²/μ²) − 3/2]. The quartic counterterm is chosen so that the CW slope
  vanishes at φ_* (one tuned parameter). Then Δη = −N_f y⁴φ_*²/(4π²V_*).
- Requiring |Δη| ≤ 0.002 (Δn_s ≤ 0.004, 1σ) gives y_nat = (0.008π²V_*/N_fφ_*²)^{1/4}.
- T-model values: 4.4×10⁻⁴ (α = 10⁻⁴) to 9×10⁻⁴ (α = 10). The flatness tuning is (y/y_nat)⁴.
- The same criterion for the bosonic λ_IH gives λ_nat = (0.004π²V_*/φ_*²)^{1/2} = 1.4–5.8×10⁻⁷. Enqvist et al. quote
  10⁻⁶ from the same requirement [Q].
- S2 is exempt (|m_N| constant).

### 2.3 Seesaw N [C]

- y_ν² = m̃M_N/v² and Γ_N = m̃M_N²/(8πv²), with v = 174.1 GeV and m̃ = 0.05 eV.
- At M_N = 1 H_i (H_i = 1.29×10¹¹): y_ν² = 2.1×10⁻⁴ and τ_N = 1.2×10⁵ H_i⁻¹.
- τ_N = 8 H_i⁻¹ needs M_N = 121 H_i, above m_I/2 for every quiet m_I. So perturbative φ → NN → LH alone cannot meet
  the τ ≈ 8–10 deadline; this confirms the higgs_reheat skeptic's 5×10⁴ H_i⁻¹.
- In the benchmarks, M_N = m_I/4 ≈ 4×10¹² GeV, so y_ν² ≈ 6.7×10⁻³, τ_N,lab ≈ 1.2×10⁴ H_i⁻¹ ≪ Γ_I⁻¹, and the N decay
  promptly on the reheating scale.

### 2.4 Plasma time [E] (`benchmarks.py: plasma_solve`)

Model:
- matter era H = 1/(1 + 1.5τ), with N and SM energy fed at Γ_I (after unblocking for S1) and at Γ_N/γ;
- the instant-preheating source of §2.1 for S1;
- a hard-gas induced Higgs mass m² = 0.018 ρ_SM/E², the coefficient matched to κ_T = 0.11 for a thermal SM gas;
- a pessimistic thermalisation time t_th = α^{−2}(E/T_*)^{1/2}/T_* (LPM-type, α = 1/40). By contrast, Harigaya–Mukaida
  (1312.3097) find near-instantaneous thermalisation; their condition α^{8/5} ≫ (m/M_P)(M_P²Γ/m³)^{1/5} holds with a
  margin of 10²–4×10³ here [Q/E].

| Channel (α = 0.3 / 1 / 3) | τ with m_ind² ≥ 10⁻³ H_i² | without instant preheating | thermal bath (pessimistic) |
|---|---|---|---|
| S1 | 0.24 / 0.66 / 0.75 | ≈ 540 / 464 / 550 | ≈ 1.9 / 2.7 / 2.8 ×10³ |
| S2 | 0.018 / 0.031 / 0.032 | — | ≈ 1.7 / 2.6 / 2.7 ×10³ |

**Caveat.** The τ ≈ 8–10 deadline was computed at H_i = 1.29–1.51×10¹¹ and M_t + 1σ. At the H_i these benchmarks need
(5–8×10¹² GeV, H_i/h_max ≈ 20–30) the Higgs analysis must be redone (higgs_reheat code). That redo is not done here.

### 2.5 The R-even Higgs portal removes the deadline [E]

- (λ_IH/2)φ_I²|H|² gives a Higgs mass² λ_IHφ_I² ≥ 0 on both sides: there is no tachyonic branch.
- In the matter era λ⟨φ²⟩/H² = 3λM_P²/m_I² is constant. λ ≥ λ_min = (10/3)(m_I/M_P)² keeps it ≥ 10 H² until the decay
  plasma takes over, and makes the Higgs heavy during inflation (λφ_*² ~ 10³ H_i² at α = 1).
- The upper bound is no broad resonance, q = λΦ_1²/4m_I² ≤ 1. Transferring the q of Enqvist et al.'s λ_hφ < 3×10⁻⁸
  (q ≈ 350) gives a looser bound [Q/E].
- Window: λ_IH ∈ [1.4–1.6×10⁻¹⁰, 1.5–2.7×10⁻⁹] (q ≤ 1), up to 5–9×10⁻⁷ (q ≤ 350), for α = 0.1–10. It lies below
  λ_nat and inside Enqvist et al.'s favoured 10⁻¹⁰–10⁻⁸ [Q].
- Ema–Mukaida–Nakayama 1602.00483 give the preheating-stage upper bounds in general [Q abstract].
- At the engine oscillator (H_i = 1.29×10¹¹): [8.4×10⁻¹⁴, 1.5×10⁻¹³] (q ≤ 1) or up to 5×10⁻¹¹.

### 2.6 T_RH, K, the matter era, leptogenesis [C/E]

- T_RH = (90/π²g_*)^{1/4}(Γ_I M_P)^{1/2} = 0.546(Γ_I M_P)^{1/2} (H = Γ).
- K and the matter era: item 7 of §0, and the table in §4.
- **Leptogenesis** (Fukugita–Yanagida [R]).
  - Non-thermal N from φ decay give n_B/s = (28/79)(3/2)ε₁T_RH/m_I. But T_RH ≈ M_N (z_RH = M_N/T_RH = 1.2–1.5) and
    K_N = m̃/m_* = 46, so the inverse-decay washout integral is ≈ 50 and the non-thermal asymmetry is erased.
  - Strong-washout thermal leptogenesis (Buchmüller–Di Bari–Plümacher hep-ph/0401240: κ_f = 2/(K z_B) = 5.2×10⁻³) gives
    n_B/s = 7.3×10⁻⁶ ε₁. The observed 8.7×10⁻¹¹ needs ε₁ = 1.2×10⁻⁵ = 3% of ε_DI = (3/16π)M_N m_atm/v² = 4.0×10⁻⁴
    (Davidson–Ibarra hep-ph/0202239 [R]).
  - **Verdict: a bonus that needs one more CP phase** (in y_ν, not related to α_I or sign(φ_I)). A maximal phase
    overproduces by ≈ 30, so it is a mild problem only if that phase is O(1). ΔL = 2 scattering at T ~ 3×10¹² adds
    some washout, which is not computed.

## 3. Item 3: consequence for M1r: dilution, H_i, r [C]

**Engine** (`code/rprh.c`, copy of rpdq.c + RH1, an optional T-model plateau, identical arithmetic for the quadratic):
- reproduction of (3, 0.5): N_X(400–1000) = −8045.485140302082, N_YL = 0.566640911258934, E_X = 188706.29275778163 (all
  digits), N_YH to 8×10⁻¹⁶; D = 2.79837, K = 1.985805;
- checks: Gauss ≤ 2.5×10⁻¹², H ≤ 3.2×10⁻¹², work ≤ 3.2×10⁻¹³.

| (b = 3, Γ_I) | \|N_X\| | E_X | N_YL | D = K^{3/2} | K | Y/Y_Q at H_M1 |
|---|---|---|---|---|---|---|
| 0.5 | 8045.5 | 188706 | 0.567 | 2.798 | 1.986 | 0.9907 |
| 0.2 | 8693.9 | 199486 | 0.577 | 3.715 | 2.399 | 0.8065 |
| 0.1 | 8978.9 | 203280 | 0.585 | 4.922 | 2.894 | 0.6286 |
| 0.05 | 9140.2 | 205261 | 0.589 | 6.726 | 3.563 | 0.4683 |
| 0.01 | 9277.9 | 206879 | 0.593 | 14.62 | 5.978 | 0.2187 |
| 0.005 | 9295.6 | 207108 | 0.593 | 20.60 | 7.514 | 0.1555 |
| → 0 (τ ≤ 400) | 9316.3 | 207141 | 0.594 | — | — | — |

- The charge is fixed early. Γ_I ≤ 10⁻³ changes |N_X| by ≤ 4×10⁻⁴, so the slow-decay yield factorises as
  Y ∝ H_i^{3/2}|N_X|/D.
- The fluid model (`code/fluid.py`) extends D to Γ_I = 10⁻¹² H_i. It agrees with the engine to ≤ 1.6×10⁻³ and moves by
  ≤ 0.1% when the match time τ is doubled.

**Scalings for the chain track [C].**
- Y/Y_Q(H_M1; b = 3) = 2.206 (Γ_I/H_i)^{1/2}.
- H_i = 1.291×10¹¹ GeV × 0.5866 (Γ_I/H_i)^{−1/3}.
- r = 2.70×10⁻⁷ (H_i/1.291×10¹¹)², a floor because H_CMB ≥ H_i.

| y (F0, b = 3, M_N = m_I/4) | Γ_I/m_I | K | H_i needed (GeV) | r (floor) | T_RH (GeV) |
|---|---|---|---|---|---|
| 10⁻⁴ | 1.3×10⁻¹⁰ | 1761 | 1.0×10¹⁴ | 0.18 | 1.7×10¹¹ |
| 3.3×10⁻⁴ | 1.4×10⁻⁹ | — | 4.7×10¹³ | **0.036** | — |
| 10⁻³ | 1.3×10⁻⁸ | 379 | 2.2×10¹³ | 8.2×10⁻³ | 7.9×10¹¹ |
| 10⁻² | 1.3×10⁻⁶ | 82 | 4.8×10¹² | 3.8×10⁻⁴ | 3.7×10¹² |
| 10⁻¹ | 1.3×10⁻⁴ | 17.6 | 1.0×10¹² | 1.8×10⁻⁵ | 1.7×10¹³ |
| 2.9 ((3, 0.5)) | 0.167 | 1.99 | 1.29×10¹¹ | 2.7×10⁻⁷ | 2×10¹⁴ |

**Mass dependence** (slow decay, τ ≤ 400) [C]:

| b | 1.5 | 2 | 2.5 | 3 | 4 | 5 | 7.15 | 11.4 | 22 | 33.2 | 102 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| \|N_X\| | 6696 | 9726 | 10930 | 9312 | 7183 | 7113 | 6163 | 8335 | 5993 | 10933 | 5897 |
| N_YL | 0.025 | 0.095 | 0.29 | 0.59 | 2.34 | 2.04 | 41 | 589 | 5497 | 102 | 57 |

- The holder's u = μ + Gφ̂_I oscillates at m_I and drives Y_L directly (resonant at m_I ≈ m_YL = 22.25).
- REVIEW_REHEATING's "m_I ≳ 30 is quiet" was for a model without the trilinear holder. **In M1r only m_I ≲ 3 H_i is
  quiet** [C].

**Direct runs at the benchmark points** (`run_engine.py bench`, τ ≤ 400 at their own Γ_I ≈ 1.5–3.4×10⁻⁶) [C]:
- |N_X| = 8586 / 10958 / 10523 / 9964 at b = 3.21 / 2.48 / 2.17 / 2.04, against 8856 / 10875 / 10144 / 9829 from the
  b interpolation: −3.1% to +3.7%, so ≤ 2.5% in H_i;
- N_YL = 0.63 / 0.28 / 0.15 / 0.11;
- the benchmark table uses the direct values.

**The released-at-rest plateau variant (RH1)** at b = 3, α = 0.17 / 0.38 / 1.04: |N_X| = 13517 / 9824 / 9488 against
9312, and K(Γ = 10⁻⁶) = 345 / 159 / 137 against 128. A plateau slice is a ±20% systematic on H_i [C/E]. The rolling
slice is for the chain track.

**The vector heat and S/Q [E].**
- The longitudinal relic number is N_A ∝ a_*³m_A, with a_* at H = m_A.
  - When Γ_I > m_A this is in radiation, a_* = (K/m_A)^{1/2}, which reproduces dq_heat's K^{3/2}.
  - When Γ_I < m_A (= 2.4×10⁻⁴ H_i at e = 10⁻⁶, every benchmark here) it is in the matter era:
    a_*³ = x₀/(3m_A²), with x₀ = ρ_I a³ = 3.3–12.4 from the background.
- The deposit per vector scales as K^{1/2}.
- This form is consistent with the literature statement that Ω_A ∝ T_RH when reheating ends after H = m [A]. It still
  needs dq_heat's mode code on a recorded slow-reheating background.
- Result: R_vec(10⁻⁶) ≈ 690–900 (vs 3.61), R_Y ≈ 0.4–2.9 (b ≤ 3.2), S/Q(10⁻⁶) ≈ 190–280.
- R_vec ∝ e^{−3/2} for a matter-era crossing, so e_crit (bulk) ≈ 7.0–7.8×10⁻⁵ for α ≥ 1. For b ≥ 3.2, R_Y alone exceeds
  R_crit = 2.14, so there is no Q-matter at any e.
- The relic would be free X far outside the computed p(S/Q) map (≤ 6.7, p = 0.017). A large symmetric part, p ≳ 1, is
  likely; the mass-per-charge correction 1 + 2p would then lower H_i by up to (1 + 2p)^{−2/3} [E].

## 4. Item 4: benchmarks (`results/benchmarks.json`, `results/summary.json`)

- The T-model member is solved self-consistently: H_req(b, Γ_I) = H_end, with N_* from the reheating history (Planck
  2018 X eq. 47 form, w = 0).
- m_I = 1.6×10¹³ GeV for every α.
- Columns: S1 values in "y" (y/y_nat), S2 values in "y_p, f".

| Family (stand-in) | b | H_i (GeV) | r | n_s | Γ_I/m_I | T_RH (GeV) | K | N_MD | N_YL | S1: y (y/y_nat) | S2: y_p, f (GeV) | plasma τ (S1 / S2) | S/Q(10⁻⁶) [E] | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T α = 10⁻⁴ (r ~ 10⁻⁷ type) | 102 | 1.64×10¹¹ | 4.5×10⁻⁷ | 0.961 | 2.2×10⁻³ | 1.6×10¹⁴ | 2.1 | 3.4 | 57 | 0.41 (930) | 0.35, 1.2×10¹³ | — | 550 | fails: holder leak R_inf = 64, N_YL = 57, y ≫ y_nat |
| T α = 10⁻³ | 33 | 5.0×10¹¹ | 4.6×10⁻⁶ | 0.961 | 4.9×10⁻⁵ | 2.4×10¹³ | 9.5 | 4.6 | 102 | 0.062 (129) | 0.054, 7.8×10¹³ | — | 2.2×10³ | fails: N_YL = 102 |
| T α = 0.01 | 11.4 | 1.46×10¹² | 4.5×10⁻⁵ | 0.961 | 1.0×10⁻⁵ | 1.1×10¹³ | 23 | 6.3 | 589 | 0.028 (53) | 0.024, 1.7×10¹⁴ | — | 500 | fails: Y resonance |
| T α = 0.03 | 7.2 | 2.3×10¹² | 1.3×10⁻⁴ | 0.961 | 7.4×10⁻⁶ | 9.4×10¹² | 30 | 6.8 | 41 | 0.024 (42) | 0.021, 2.0×10¹⁴ | — | 370 | fails: N_YL = 41 |
| T α = 0.1 | 4.5 | 3.6×10¹² | 4.4×10⁻⁴ | 0.962 | 2.4×10⁻⁶ | 5.3×10¹² | 52 | 7.9 | 2.2 | 0.014 (22) | 0.012, 3.5×10¹⁴ | 0.08 / 0.009 | 178 | marginal (R_Y = 13) |
| **T α = 0.3** | 3.2 | 5.1×10¹² | 1.3×10⁻³ | 0.962 | 1.12×10⁻⁶ | 3.6×10¹² | 82 | 8.7 | 0.63 | 0.0093 (14) | 0.0081, 5.1×10¹⁴ | 0.24 / 0.018 | 193 | **benchmark B1** |
| **T α = 1** | 2.5 | 6.6×10¹² | 4.2×10⁻³ | 0.962 | 6.2×10⁻⁷ | 2.7×10¹² | 124 | 9.3 | 0.28 | 0.0069 (9.5) | 0.0060, 6.8×10¹⁴ | 0.66 / 0.031 | 204 | **benchmark B2** |
| **T α = 3** | 2.2 | 7.4×10¹² | 1.2×10⁻² | 0.963 | 7.5×10⁻⁷ | 2.9×10¹² | 136 | 9.4 | 0.15 | 0.0076 (9.5) | 0.0066, 6.1×10¹⁴ | 0.75 / 0.032 | 249 | **benchmark B3** |
| T α = 10 | 2.0 | 7.8×10¹² | 3.2×10⁻² | 0.963 | 8.9×10⁻⁷ | 3.2×10¹² | 138 | 9.3 | 0.11 | 0.0083 (9.2) | 0.0072, 5.5×10¹⁴ | 0.79 / 0.031 | 280 | near r bound |
| F0 engine oscillator | 3 | 2.2×10¹³ at y = 10⁻³ | 8×10⁻³ (floor) | — | 1.3×10⁻⁸ | 7.9×10¹¹ | 379 | 11.7 | 0.59 | 10⁻³ (y_nat 2.3×10⁻³ at φ_* = 0.82 M_P, lenient) | — | — | 324 | no potential; r < 0.036 needs y ≥ 3.3×10⁻⁴ |

**Notes on the table.**
- For B1–B3, M_N = m_I/4 ≈ 4×10¹² GeV, m̃ = 0.05 eV, τ_N,lab ≈ 1.2×10⁴ H_i⁻¹, z_RH = 1.2–1.5, thermal leptogenesis
  ε₁ ≈ 1.2×10⁻⁵.
- Higgs portal window: λ_IH ∈ [1.5×10⁻¹⁰, 2×10⁻⁹].
- The holder is re-fixed at φ_end = 0.67 / 0.84 / 0.94 M_P: g_I ≈ 3.8–4.0×10⁻⁵.
- S1 blocking: yΦ_1/(m_I/2) = 666–861, unblocking at τ ≈ 440–570. The Greene–Kofman first-crossing energy fraction is
  1.4–2.6×10⁻⁵.
- S3 (φGG̃) needs 1/g = 3–4×10¹⁵ GeV and gives ξ_A = 5.0–12.5 > 2.5 for B1–B3, so it fails.
- S2 is the radiatively natural choice. S1 needs a 10⁴ flatness tuning.

## 5. New issues for the chain and the working theory

1. **[C] The holder drives the mediator for m_I ≳ 4 H_i** (N_YL ≥ 2, up to 5.5×10³ at m_I ≈ m_YL). With the
   trilinear holder, the quiet inflaton masses are m_I ≲ 3 H_i. REVIEW_REHEATING's "≳ 30 quiet" does not carry over to
   M1r.
2. **[C/E] Inflaton → hidden leakage through the holder** (§0 item 9). Above 2m_YL = 44.5 H_i it is fatal for any
   Γ_I ≲ m_I. The "induced portals" (A_IΣ φ|Σ|² → Goldstone pairs, open at any m_I) obey R_inf = 3M_P²K²BR/E_X ≤ O(1).
   At benchmark B2 (K = 123, M_P/H_i = 3.7×10⁵) the coefficient is 2.2×10¹⁰. So R_inf ≤ 1 needs BR ≤ 4.5×10⁻¹¹, i.e.
   Γ_hid ≤ 7×10⁻¹⁷ H_i [C]. This tightens the renormalisation condition the decisive pass left open.
3. **[E] "Γ_I = 0.5 H_i" is not realisable**: perturbative channels give Γ_I/m_I ~ 10⁻⁶ in the self-consistent
   stand-ins. The supplied (3, 0.5) numbers (H_i = 1.29×10¹¹, r = 2.7×10⁻⁷, S/Q = 4.955, K = 1.986) should be read as
   belonging to an unrealised history.
4. **[E] The branch.** With realistic slow reheating the vector heat is ~10² times larger (matter-era crossing, √K
   deposit), so S/Q ~ 10². This is opposite in direction to the rolling-pre-slice effect (which lowers S/Q). It needs a
   p(S/Q) map beyond 6.7, and dq_heat's vector modes on the recorded slow background.
5. **[A] Higgs at H_i ~ 6×10¹² GeV.** The higgs_reheat windows were computed at 1.3–1.5×10¹¹. The λ_IH portal makes the
   Higgs heavy during inflation and through reheating, which is the cleanest fix. Without it the inflation-stage ξ_min
   must be recomputed at the new H_i.
6. **[Q/E] S3** (shift-symmetric φGG̃) fails the ξ_A ≲ 2.5 gauge-production bound in every quiet member.

## 6. What is not established

- **The stand-in families.** The T-model is not the potential track's output. Its N_X and D come from the static-slice
  engine at b = m_I/H_end; the plateau-at-rest variant moves H_i by ±20%.
- **The plasma model.** The hard-gas mass coefficient, the LPM thermalisation time and the instant-preheating refill
  (Pauli and decay balance) are estimates.
- **S2's axial non-perturbative production.** Landau–Zener estimate only, near its EFT cutoff.
- **The S/Q rescaling.** The matter-era vector crossing (a_*³m_A ansatz), the ledger's K^{1/2} decay factor (valid while
  Γ_Y_L, Γ_A ≪ Γ_I) and the absence of a p map at S/Q ~ 10² are open. These are [E], not [C].
- **Flatness criterion.** One tuned counterterm, |Δη| ≤ 0.002. A different criterion moves y_nat by O(1).
- **Not computed.** Leptogenesis ΔL = 2 washout; the Higgs at the new H_i; the induced-portal values;
  fragmentation/oscillons (Lozanov–Amin 1608.01213 and 1710.06851 find fragmentation only for plateau scales M ≪ M_P;
  the quiet members have √(6α) ≥ 1.3 M_P [Q abstract]).

## 7. Files and reproduction

| File | Content |
|---|---|
| `code/rprh.c` | engine: copy of `tracks/dq_entropy/code/rpdq.c` + RH1 (T-model plateau option; quadratic path identical) |
| `code/rh_common.py` | constants and driver (`run` is a copy of dq_common.run with wT; fx_common imported read-only) |
| `code/run_engine.py` → `results/engine_runs.json` | reproduction; Γ scan at b = 3; Γ → 0; b scan 1.5–102; benchmark points; plateau-at-rest variants |
| `code/fluid.py` → `results/fluid.json` | K, D, N_eq, N_RH for b = 1.5–5, Γ_I = 0.5–10⁻¹²; validation against the engine |
| `code/channels.py` → `results/channels.json` | width conventions (trace check), Option B, blocking and GK, seesaw, holder leakage, λ_IH window, thermalisation |
| `code/fermion_modes.py` → `results/fermion_modes.json` | numerical Greene–Kofman check through two crossings |
| `code/benchmarks.py` → `results/benchmarks.json` | T-model self-consistent benchmarks (S1/S2/S3), F0 rows, plasma time, leptogenesis, S/Q [E] |
| `code/summary.py` → `results/summary.json` | headline numbers and the benchmark table |

Reproduction, in `code/`, one thread:

```
python3 -B run_engine.py repro gscan early bscan tmodel   # ~45 s
python3 -B run_engine.py bhigh                            # ~10 s
python3 -B benchmarks.py; python3 -B run_engine.py bench; python3 -B benchmarks.py   # bench points need a first pass
python3 -B fluid.py                                       # ~100 s
python3 -B channels.py                                    # ~2 s
python3 -B fermion_modes.py                               # ~2 min
python3 -B summary.py
```

The scratch build and run series are in `$SP/wt2/inf_reheat/` (`build/rprh`, `runs/*.npz`).

**Literature** (abstract-level unless stated):
- Greene–Kofman hep-ph/9807339 and hep-ph/0003018; Peloso–Sorbo hep-ph/0003045 (the brief's hep-ph/0001263 is an
  unrelated paper; corrected here); Felder–Kofman–Linde hep-ph/9812289 [R].
- Davidson–Sarkar hep-ph/0009078; Harigaya–Mukaida 1312.3097.
- Enqvist–Karčiauskas–Lebedev–Rusak–Zatta 1608.08848 (λ_hφ ≲ 3×10⁻⁸ and |σ_hφ| ≲ 10⁸ GeV at m = 1.6×10¹³ GeV;
  favoured 10⁻¹⁰–10⁻⁸; read from the ar5iv text); Ema–Mukaida–Nakayama 1602.00483.
- Lozanov–Amin 1608.01213 and 1710.06851; Adshead–Sfakianakis 1508.00891 [R].
- Davidson–Ibarra hep-ph/0202239 [R]; Buchmüller–Di Bari–Plümacher hep-ph/0401240 [R]; Fukugita–Yanagida (1986) [R];
  Barnaby–Peloso 1011.1500 [R]; BICEP/Keck 2110.00483 [R]; Planck 2018 X (1807.06211) [R].
