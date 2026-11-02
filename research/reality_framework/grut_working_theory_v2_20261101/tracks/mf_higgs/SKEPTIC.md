# SKEPTIC of track mf_higgs (key mf_higgs_skeptic)

2 November 2026 · Claude, adversarial pass on `tracks/mf_higgs` (BRIEF_MULTIFIELD.md, track mf_higgs; GitHub issue #13
context).

This is same-session AI work on one supplied model. It is not outside replication. Nothing here is a validation, an
exclusion of the real world, a new theory or a GRUT prediction.

**Every choice in the track and here is a model choice.** That covers:
- the benchmark (HFa/HFb) and its post-waterfall closure (unblocked, or Pauli-blocked);
- ξ_H = 0.1, read as ξ(M_t);
- the survival criteria (Hawking–Moss, the argued danger exponent D, the mean-field edge);
- the medium masses;
- the renormalisation conditions;
- the hypothetical G6 portal.

None of these follows from closure. Data comparisons (PDG, Planck, DESI) are retrospective. Git was read-only and nothing
was committed. Every run was in the foreground with one process, and there is no compiled file or `__pycache__` in the
tree.

Labels:
- [D] derived;
- [C] computed with my own code (`mf_higgs/skeptic/code/`, output in `mf_higgs/skeptic/results/`);
- [E] estimate, or an argued criterion;
- [A] argued;
- [Q] quoted from a track;
- [R] literature, with arXiv numbers checked by web search in this session.

**Method policy.**
- **Shared on purpose.** The SM input is `higgs_reheat/code/smpot.py`, imported read-only: 3-loop running, the flat 2-loop
  λ_eff, ξ running and the one-loop de Sitter potential. The brief says not to redo the SM.
- **Not imported.** The track's code (`mf_higgs/code`, and its use of `inf_higgs/code`) is not imported. Everything else is
  re-written:
  - Barrier: four potential and RG schemes (§1).
  - Higgs modes: the friction-free variable u_k = a^{3/2}h_k, with Ω_k² = k²/a² + (ξ − ¼)R + ¾H². This is an exact
    rewrite using R = 6(2H² + Ḣ).
  - Integrator: adaptive DOP853 per stage on cubic-spline backgrounds, with a Bunch–Davies start from Hankel functions
    (or the adiabatic vacuum when m > 3H/2).
  - D: my own implementation of the argued criterion.
- **A second background code path.** Besides mf_hybrid's export, I use mf_chain's rebuilt histories
  (`mf_chain/results/histories/bench_itHFa_track_it1.json`, unblocked, H_i = 1.284×10¹¹, K = 2.04; and
  `bench_itHFa_P0_it2.json`, Pauli Γ_th = 0, H_i = 2.40×10¹¹, K = 7.84). Their stage 1 carries the Pauli friction and the
  seesaw decays from the slice onward, so they are a different construction from the track's fluid, which restarts from
  the unblocked bench at t = 0.19 [Q, mf_chain unreviewed].
- **Plasma.** The medium masses are formed from the histories' own energy fractions, with the quantum energies bounded
  rather than evolved (§2.4).

Units and conventions, checked:
- H_i is the Hubble rate at σ = σ_c, with a = 1 there, as in the track and the bench.
- ε_H = 1 falls at τ = 0.1927 (bench) or 0.178 (mf_chain's Pauli rebuild).
- M_P = 2.435×10¹⁸ GeV (reduced).
- N_obs = N_* + ln(k_*/a₀H₀) = 52.644 + 5.405 = 58.049 at K = 2.048, with N_* = const − ½ ln K + ½ ln H_i. This is
  reproduced, including mf_chain's N_* = 52.280 at (2.40×10¹¹, K = 7.84).
- The track's CMB inputs enter only through N_*. CMB pivot quantities are mf_hybrid's (and its skeptic's) and are not
  rechecked here.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

| Track result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| H1 constant-H Hawking–Moss, ξ_EW,min = 0.049 / 0.053 (tree), cap > 174.6 GeV | **stands** | **Own barrier, scheme A (the track's tree, μ² = h² + R):** ξ_min = 0.0489 / 0.0534 / 0.0726 at 172.57 / 172.86 / 174.6, identical to the track.<br>**Other schemes are safer:** flat 2-loop λ_eff + ξR gives 0.041 / 0.046 / 0.068; the one-loop de Sitter potential (1804.02020) gives 0.042 / 0.048 / 0.069 at μ² = h² + R or μ = h.<br>**Past light cone** (η₀ = 3.2/a₀H₀): +0.0003. **Cap:** above 177 GeV in every scheme and at α_s ± 1σ (extrapolated beyond the SMDR grid end at 174.6). | C |
| H2 instantaneous-HM bracket: cap 174.12 (HFa) / 173.95 (HFb); 173.78 / 174.47 at α_s ∓ 1σ; ξ_min 0.066 / 0.074 | **refuted (numbers); the bracket is not a bracket** | The track's stride-4 table stops at t = 0.19207, where R/H² = 6.74, not at ε_H = 1 (t = 0.19266, R/H² = 6); there I reproduce its −156.6 / −93.1 / −31.9 / −2.6 exactly.<br>**With the ε_H = 1 row:** log₁₀⟨N⟩ = −112.1 / −60.8 / −11.5 / +12.0 at 172.57 / 172.86 / 173.45 / 174.05. The cap is **173.69 GeV for both HFa and HFb**, 173.35 / 174.03 at α_s = 0.1171 / 0.1189, and ξ_min = **0.071 / 0.080**.<br>**Past ε_H = 1:** B(172.57) falls to 115 / 15 / 1.7 at R/H² = 3 / 1 / 0.3, which the bench reaches within 0.003 H⁻¹. Read literally, the instantaneous count fails at every M_t, so it bounds nothing (§3.1). | C |
| H3 induced masses; N_s seesaw −0.68 / −1.8 / −2.5 H_i², "needs μ ≳ 0.37 / 25 / 47 H_i" | **numbers stand; the "fatal" reading is conditional** | **Reproduced:** −0.68 / −1.85 / −2.56 H_i² at M_soft = 0.3 / 1 / 1.5 H_i, at m̃_s = 0.05 eV.<br>**Smaller m̃_s:** the tolerance holds for m̃_s ≤ 0.048 / 0.018 / 0.013 eV under "vanish at M_soft". N_s needs only m̃_s ≳ 3×10⁻⁷–7×10⁻⁶ eV to decay before it dominates, or no seesaw at all (dark radiation).<br>**Σm_ν:** 0.05 eV per singlet is already disfavoured by Σm_ν. The s-odd piece makes the inflationary Higgs mass CP-branch dependent by 0.003–0.065 H_i² (§4, N5). | D/C |
| H4 R < 0 windows (38%, 0.005 H⁻¹, −4.7 / −5.2 H², exponent 0.0023–0.0026) and the R lines | **stands** | Bench export: R_min/H² = −4.67; R < 0 for 36% of [t(ε_H = 1), 0.6]; largest window exponent 0.0023, summed 0.035.<br>mf_chain Pauli rebuild: −5.52, 32%, 0.0024.<br>Mathieu q of the R line is 0.018 (export) / 0.049 (Pauli) even at ξ = 1000. The line sits at ω = 481 H_i on the export but at 413 H_i on the Pauli rebuild (anharmonic χ); the track uses 480 for all histories, which is immaterial. | C |
| H5 plasma timing: m_ind² ≥ 1 / 10² / 10³ H_i² at τ = 0.191 / 0.21 / 0.24–0.29; ρ_SM/ρ ≥ 10⁻³ at 0.26–0.33 | **stands, and strengthened** | On mf_chain's rebuilds the SM hard gas reaches ≥ 1 / 10² / 10³ H_i² at τ = 0.187 / 0.208 / 0.260 (unblocked) and 0.172 / 0.194 / 0.294 (Pauli).<br>The export's f_SM = 0 before 0.6 is confirmed to be an artefact.<br>**New:** the singlet-neutrino gas itself gives m² = y_ν²ρ_N/E_N² ≥ 1 / 10² / 10⁴ H_i² at τ = 0.13 / 0.17 / 0.18–0.19, before ε_H = 1 (§2.4, N2). | C on E |
| H6 D at ξ_EW = 0.1, gravitational only | **stands** | Own modes and D:<br>• export: −1819 / −298.5 / **−5.9** at 172.28 / 172.57 / 172.86 (track −1656 / −299 / −7.8);<br>• mf_chain Pauli rebuild at 1.3×10¹¹: −2012 / −353 / **−22.6** (track Pauli Γ_th = 0: −1658 / −349 / −21.9).<br>The 172.86 difference comes mostly from the shared SMTab's h_max grid (Δ ln h = 1.15%). | C |
| H7 M_t,crit gravitational only: 172.88–172.92, band 172.83–173.00; α_s; M_h; H_i | **stands** | **M_t,crit:** export 172.877 (band 172.824–172.976); unblocked rebuild 172.882; Pauli rebuild at 1.3×10¹¹: 172.915 (172.887–173.002).<br>**At mf_chain's self-consistent 2.40×10¹¹:** 172.737 (172.712–172.818), with **D(172.86) = +34.5**.<br>**α_s 0.1171 / 0.1189:** 172.552 / 173.205 (export); 172.410 / 173.066 at 2.4×10¹¹, where central M_t gives D = +40.8.<br>**M_h ∓ 1σ:** 172.818 / 172.939 (track 172.840 / 172.987).<br>**Deadlines** at 173.0 / 173.15 / 173.45 / 173.75: 7.0 / 4.8 / 3.05 / 2.3 (export), 11.8 / 7.6 / 4.75 / 3.5 (Pauli).<br>The track's bracketed "[N = 60]" values (172.885, 172.937) are coarse-grid artefacts. At N = 60 the crossings are **172.870 / 172.906**, below the N_obs values as they must be (N6). | C |
| H8 D with the plasma: no roll ≤ 174.6 GeV; c × 10⁻⁴ still gives D ≤ −107 | **stands, but the margin depends on the energy law** | **N gas + SM gas** in the criterion: safe at 174.6 down to ×1.8×10⁻⁴ (export) / ×5.6×10⁻⁴ (Pauli) of the medium mass with the injection energies held fixed (conservative), and down to ×10⁻⁶ with redshifting energies.<br>**N gas alone (×1):** D(174.6) = −756 / −4581.<br>**SM hard gas alone at ×10⁻⁴, energies fixed:** D(174.05 / 174.6) = +2.6 / +25 (export). So the track's "×10⁻⁴" relies on its redshifting E_SM (N4). | C on E |
| H9 ξ window: lower 0.049 (0.053); upper ξ_EW ≈ 113–155 (21–24) | **stands** | **Lower edge, inflation:** H1.<br>**Lower edge after inflation, gravitational only (own ξ scan):** 0.0675 / 0.0985 (export) and 0.0665 / 0.0956 (Pauli rebuild) at 172.57 / 172.86; the track has 0.067 / 0.098 and 0.066 / 0.095.<br>**Upper edge** (2σ_prod = φ_max, own modes on the export, ξ(H_i) = 10–300): ξ_EW ≈ 111 at 172.57 and ≈ 21 at 172.86. | C [E] |
| H10 G6 window −2×10⁻⁸ ≲ λ_χH ≲ −7×10⁻¹² | **stands (the edge is uncertain to ×3)** | **Virial mean:** confirmed on the exact anharmonic valley orbit. The track's harmonic (A/v)² underestimates it by 1–19% at A/v = 0.1–0.3 and by 72% at 0.4.<br>**Rebuilt oscillation:** its negative swing is about 20% too shallow (−0.57 against −0.69 at A/v ≈ 0.3).<br>**Floquet per-mode growth through the first band:** ≤ 0.21 at Λ = 3×10³; 1.7–2.5 at 10⁴; 15–21 at 3×10⁴, with a k = 0 tachyonic instability. This is consistent with the runaway at Λ ≥ 10⁴. | D/C |
| O1 open items | — | the lattice and the N_s sector remain open; see §7 | — |

**Headline.**
- **The track survives an independent recomputation almost everywhere.** The constant-H Hawking–Moss edge, every
  gravitational-only D and M_t,crit, the deadlines, the R windows, the upper edge and the G6 window are reproduced to
  ±0.005 GeV in M_t,crit (±0.05 GeV for the M_h variant) and ±2 in D, on the export and on mf_chain's independently
  built histories.
- **What does not hold** is the "pessimistic instantaneous bracket" (cap 173.95–174.12 GeV). Its value is a sampling
  artefact (173.69 GeV at the exact end of inflation). More importantly, it is not a bracket: the same reading carried
  0.003 H⁻¹ further gives decay at every M_t.
- **The post-inflation rescue is stronger than the track says.** Beyond the SM hard gas, the singlet-neutrino gas
  produced by the waterfall gives the Higgs a positive mass y_ν²ρ_N/E_N² ≳ 10³–10⁴ H_i² already before ε_H = 1. Its sign
  is fixed by the seesaw mixing, with no gauge or composition ambiguity.
- **The gravitational-only verdict is marginal and history-dependent.** At mf_chain's self-consistent Pauli H_i
  (2.4×10¹¹ GeV) it fails at +1σ (D = +34.5), and at central M_t for α_s = 0.1171 (D = +40.8). With the medium masses it
  does not matter.

---------------------------------------------------------------------------------------------------------------------

## 1. What was recomputed, and how

| Check | Track | Skeptic | File |
|---|---|---|---|
| Hawking–Moss barrier | `inf_higgs.barrier_B`, tree μ² = h² + R; one-loop dS | own scan and refine; schemes A (track tree), B (flat 2-loop λ_eff(h) + ξ(h)R e^{2Γ}h²/2), C (one-loop dS, μ² = h² + R), D (one-loop dS, μ = h); light-cone factor; N₀ for the Pauli K | `sk_infl.py` |
| Instantaneous count | stride-4 bench rows to t ≤ t(ε_H = 1) | all rows + the exact ε_H = 1 row (R = 6H²); B past ε_H = 1 | `sk_infl.py inst` |
| Higgs modes | h_k, ḣ_k, fixed-step RK4, track Background (raw R to 0.6, rebuilt lines) | u_k = a^{3/2}h_k, DOP853 (rtol 10⁻⁸) per stage, cubic splines, no rebuilt lines; Wronskian error 5×10⁻⁷ | `sk_modes.py` |
| Backgrounds | export; own fluid from the bench state at t = 0.19 | export; mf_chain's unblocked and Pauli rebuilds (Pauli friction and seesaw decays from the slice) | `sk_modes.Hist` |
| D | `mh_modes.danger` (SMTab h_max on a 1.15% grid) | own implementation; h_max refined by a parabola; same argued definition and abundance limit x ≤ 3N + 40 | `sk_modes.danger` |
| Medium masses | own fluid with the SM number density; c ρ_SM/E_SM² | the histories' own fractions; SM hard gas and the new N-gas term; energies bounded (fixed) or redshifted | `sk_plasma.py`, `sk_medium.py` |
| G6 | harmonic rebuild, mode runs | exact anharmonic valley orbit (virial); Floquet exponents, exact and harmonic | `sk_g6.py` |
| Upper edge | modes to ξ(H_i) = 1000 | own modes at ξ(H_i) = 10–300 (adiabatic start) | `sk_upper.py` |
| N_s | one-loop formula at m̃ = 0.05 eV | m̃_s dependence, decay-before-domination bound, s-odd piece | `sk_ns.py` |

**Numerical controls of my mode code** (export, ΔD at 172.86 / 172.88):
- rtol 10⁻⁹: 0.000; t_end = 150: 0.000;
- n = 400 modes: −0.6 / +0.3;
- t₀ = −5: −0.1;
- K_max = 15: −1.1 (this is the Hartree subtraction convention);
- Hartree off: +1.96 / +1.91 (the track found +2.0).

Together these move M_t,crit by at most 0.003 GeV.

## 2. Reproductions (all [C])

### 2.1 Inflation, constant H (`results/sk_infl.json`)

| M_t | 172.57 | 172.86 | 173.45 | 174.05 | 174.6 |
|---|---|---|---|---|---|
| ξ_min, A (track's tree) | 0.0489 | 0.0534 | 0.0609 | 0.0674 | 0.0726 |
| ξ_min, B (flat 2-loop λ_eff + ξR) | 0.0406 | 0.0464 | 0.0553 | 0.0625 | 0.0682 |
| ξ_min, C / D (one-loop dS) | 0.0419 / 0.0418 | 0.0475 / 0.0475 | 0.0564 / 0.0563 | 0.0635 / 0.0634 | 0.0691 / 0.0691 |
| log₁₀⟨N⟩(ξ_EW = 0.1), A / B / C | −609 / −1186 / −1098 | −427 / −713 / −672 | −247 / −358 / −342 | −159 / −217 / −208 | −111 / −149 / −142 |

- **Light cone.** Replacing 1/(a₀H₀) by the past-light-cone radius η₀ ≈ 3.2/(a₀H₀) adds 3 ln 3.2 = 3.5 to ln⟨N⟩, i.e.
  +0.0002 to +0.0004 in ξ_min.
- **Pauli K.** Using N₀ for the Pauli K changes ξ_min by −0.0001 to −0.0002.
- **Cap.** log₁₀⟨N⟩(ξ = 0.1) stays below zero up to 177 GeV in every scheme (−21 / −32 / −29 at 177; extrapolated) and at
  α_s = 0.1171 / 0.1189.
- **Verdict.** The track's constant-H edge is the most pessimistic of the four schemes.

### 2.2 After inflation, gravitational only (`results/sk_postD.json`)

| history (H_i) | N_obs | D at 172.28 / 172.57 / 172.86 | M_t,crit (band) | N = 60 | α_s 0.1171 / 0.1189 | M_h 125.09 / 125.31 |
|---|---|---|---|---|---|---|
| export (1.30×10¹¹) | 58.049 | −1819 / −298.5 / −5.9 | **172.877** (172.824–172.976) | 172.870 | 172.552 / 173.205 | 172.818 / 172.939 |
| mf_chain unblocked (1.284×10¹¹) | 58.045 | −1843 / −306 / −7.7 | **172.882** (172.828–172.979) | 172.875 | 172.556 / 173.210 | 172.821 / 172.942 |
| mf_chain Pauli, SM at 1.30×10¹¹ | 57.378 | −2012 / −353 / −22.6 | **172.915** (172.887–173.002) | 172.906 | 172.589 / 173.241 | 172.854 / 172.977 |
| mf_chain Pauli at its own 2.40×10¹¹ | 57.684 | −925 / −108 / **+34.5** | **172.737** (172.712–172.818) | 172.730 | **172.410** / 173.066 | 172.679 / 172.798 |

- **Track for comparison:** 172.880 (exported), 172.914 (Pauli Γ_th = 0), ≈ 172.74 at 2.4×10¹¹, and 172.562 / 173.246 for
  α_s.
- **Where the D values differ** (172.28, where D is irrelevant; and ±2 at 172.86), the cause is the region grid and the
  shared SMTab's h_max quantisation (1.2394 against 1.2445×10¹¹ GeV at 172.86).
- **α_s = 0.1171, D at central M_t:** +6.6 (export; track +6.8), −6.6 (Pauli at 1.3×10¹¹), **+40.8** (Pauli at 2.4×10¹¹).
- **Deadlines** D(τ_th) = 0 at 173.0 / 173.15 / 173.45 / 173.75:
  - export: 6.99 / 4.82 / 3.05 / 2.32 (track 7.0 / 4.8 / 3.0 / 2.3);
  - Pauli at 1.3×10¹¹: 11.8 / 7.6 / 4.75 / 3.47 (track 11.5 / 7.3 / 4.5 / 3.3);
  - Pauli at 2.4×10¹¹: 7.2 / 5.7 / 3.9 / 3.0.

### 2.3 Windows, upper edge, G6 (`results/sk_windows.json`, `sk_upper.json`, `sk_g6.json`)

- **R < 0 windows.** See §0, row H4.
- **Upper edge (2σ_prod/φ_max).**

  | ξ(H_i) | 10 | 30 | 100 | 300 |
  |---|---|---|---|---|
  | at 172.57 | 0.24 | 0.38 | 0.80 | 1.98 |
  | at 172.86 | 0.69 | 1.07 | 2.28 | 5.62 |

  Interpolated in ln ξ, the edge is at ξ(H_i) ≈ 131 (ξ_EW ≈ 111) and 25 (ξ_EW ≈ 21). The track has 134 / 29 (export).
- **G6, virial mean.** On the exact orbit, ⟨1 − χ²/v²⟩ / (2ρ_χ/m_χ²v²) = 1.003 / 1.013 / 1.06 / 1.19 / 1.72 at turning
  amplitude A/v = 0.05 / 0.1 / 0.2 / 0.3 / 0.4. Orbits with A/v ≥ √2 − 1 cross χ = 0.
  - The track's statement "mean = +Λf_χH²/4 > 0" holds and is conservative.
  - The rebuilt −2(A/v)cos(m_χt) + (A/v)² has a negative swing about 20% too shallow at A/v ≈ 0.3, and a frequency
    10–20% too high at large amplitude.
- **G6, Floquet.** Largest exponent, per-mode growth μ × q/2H through the first band:

  | Λ | at τ ≈ 0.25–0.3 (export / Pauli) | at τ = 1–3 |
  |---|---|---|
  | 3×10³ | 0.16 / 0.21 | 0.01–0.03 |
  | 10⁴ | 1.7 / 2.5 | 0.18–0.50 |
  | 3×10⁴ | 15 / 21 (k = 0 tachyonic, μ = 31–47 H_i) | 1.7–4.6 |

  The upper edge Λ ≈ 3×10³–10⁴ is confirmed to within the factor-3 grid.

### 2.4 Medium masses after the waterfall (`results/sk_plasma.json`, `sk_medium.json`)

**Model** (for a region field φ̄):
- **Formula.** m² = Σ_i κ_i n_i⟨1/E_i⟩, with m_i²(h) = m_i²(0) + κ_i h². This is the forward-scattering (medium) mass.
  The same formula gives the textbook thermal masses, e.g. y²T²/6 for a Dirac fermion [D].
- **SM hard gas** from N → L H: c ρ_SM/E_SM², with c = 0.0228 (the track's value) and E_SM = E_N/2.
- **N gas (new).** The seesaw mixing in a Higgs background raises the heavy mass, M_h² = M² + y_ν²h² (m_D = y_νh/√2), so
  κ_N = y_ν² > 0 and m² = y_ν²ρ_N/E_N². I take the smallest Yukawa, y_ν² = 0.026 (N₄; 0.086 for N₂), at m̃ = 0.05 eV.
- **Energies.** E_N is held at its injection value ¾μ5,max at ε_H = 1 (8.2–8.3×10³ H_i at 1.3×10¹¹; 4.5×10³ at 2.4×10¹¹).
  This is conservative, since m² ∝ 1/E².

| history | N gas m² ≥ 1 / 10² / 10³ / 10⁴ H_i² at τ = | SM hard gas ≥ 1 / 10² / 10³ at τ = | t(ε_H = 1) |
|---|---|---|---|
| export (f_SM = 0 before 0.6: an artefact) | 0.138 / 0.176 / 0.186 / 0.193 | (only after 0.6) | 0.193 |
| mf_chain unblocked | 0.138 / 0.176 / 0.186 / 0.193 | 0.187 / 0.208 / 0.260 | 0.193 |
| mf_chain Pauli (2.4×10¹¹) | 0.126 / 0.163 / 0.172 / 0.180 | 0.172 / 0.194 / 0.294 | 0.178 |

- **Size against the barrier.** Inside the abundance limit (φ̄ ≤ 3.4 H_i), −V″ ≤ 0.48–0.50 H_i² even at 174.6 GeV. So
  either medium term exceeds it by 10²–10⁵ from the end of inflation onward.
- **Placed in the criterion** (own modes):
  - every combination gives D < 0 to 174.6 GeV;
  - the smallest safe factor on N + SM is 1.8×10⁻⁴ (export) or 5.6×10⁻⁴ (Pauli) with fixed energies, and ≈ 10⁻⁶ with
    redshifting energies.

---------------------------------------------------------------------------------------------------------------------

## 3. Refutation

### 3.1 H2: the instantaneous-HM "pessimistic bracket" [C] (`sk_infl.py inst`, `results/sk_infl.json`)

**What the track computed.**
- `mh_inflation.waterfall_table` keeps bench rows with t ≤ t(ε_H = 1), and `log10N_inst` uses every fourth row.
- The last row used is t = 0.19207, with R/H² = 6.74 and B(172.57) = 529.7.
- ε_H = 1 is at t = 0.19266, with R/H² = 6.00 and B = 425.1.

**Why the endpoint matters.** The integrand e^{3(N₀−N_b) − B(N_b)} is dominated by the very last instant: dB/dN ≈ 10⁵
there. So the 6×10⁻⁴ H⁻¹ that the stride skips changes log₁₀⟨N⟩ by 20–45.

| M_t | 172.57 | 172.86 | 173.45 | 174.05 |
|---|---|---|---|---|
| log₁₀⟨N⟩, track's sampling (reproduced exactly) | −156.6 | −93.1 | −31.9 | −2.6 |
| with the ε_H = 1 row | **−112.1** | **−60.8** | **−11.5** | **+12.0** |

**Corrected values.**
- The cap at ξ_EW = 0.1 is **173.69 GeV for HFa and 173.69 for HFb**. The track's HFa/HFb difference (174.12 against 173.95)
  is the same sampling artefact.
- At α_s = 0.1171 / 0.1189 the cap is 173.35 / 174.03.
- ξ_min(instantaneous) is 0.071 / 0.080 at 172.57 / 172.86.
- The Gaussian non-equilibrium tail at the exact ε_H = 1 row (own modes, σ²(k < aH) = 0.0284) is log₁₀ = −799 / −564 / −158
  at 172.57 / 172.86 / 174.6, against the track's −974 / −696 / −210, which were sampled at R/H² = 7.27. These are still far
  below zero.

**Why it is not a bracket.**
- The bench reaches R/H² = 3 / 1 / 0.3 within the next 0.003 H⁻¹. There B(172.57) = 115 / 15 / 1.7.
- Continued even that far, the literal instantaneous count gives ⟨N⟩ ≫ 1 at every M_t.
- So the instantaneous count is not a pessimistic bracket of the constant-H baseline. It is undefined for a sub-Hubble end
  of inflation, as the track argued for its validity [A]. It cannot set an "overall reach".
- **What remains.** The constant-H count (cap > 177 GeV, extrapolated) and the non-equilibrium Gaussian count (≤ 10⁻¹⁵⁸
  at 174.6) are the meaningful inflationary criteria.
- **Correction to the track.** The digest's "overall reach … 173.95–174.12 GeV on the pessimistic bracket" should read:
  "173.7 GeV only if the instantaneous count is stopped exactly at ε_H = 1, an arbitrary point; otherwise undefined".

## 4. New issues (ranked)

1. **N1 [C] Instantaneous-bracket sampling and ill-definedness** (§3.1).
   - With the exact end the cap is 173.69, not 174.12 / 173.95, and HFa = HFb.
   - Past ε_H = 1 the criterion diverges. It should not be quoted as a reach.

2. **N2 [D/E, favourable] The singlet-neutrino gas is itself a Higgs-coupled medium.**
   - **The term.** m² = Σ_j y_νj² ρ_Nj/E_Nj², positive by the seesaw level repulsion (M_h = M + m_D²/M).
   - **Timing.** It exists as soon as the N are produced, during the last χ infall: ≥ 10³ H_i² at τ = 0.17–0.19, before
     ε_H = 1, on all three histories, with the smallest Yukawa and unredshifted energies.
   - **What it removes.** It does not need N decays or any SM gas. Its sign carries no gauge or composition
     ambiguity (compare N3).
   - **What it adds to the track.** The track's "[D] the axial N give no field-dependent Higgs mass" is true of the
     vacuum loop only.
   - **Scaling.** It scales ∝ m̃: ×0.3 at 0.015 eV.

3. **N3 [A] The SM hard-gas sign is gauge-sensitive; N2 makes it moot.**
   - c > 0 (the track's 0.0228, or 0.026–0.032 from the composition count) counts hard Goldstone quanta as longitudinal
     W/Z, with κ = g²/4 or (g² + g′²)/4.
   - In Landau gauge, a Goldstone off a stationary point carries κ = λ < 0, which gives c ∝ 6λ/8 < 0 for λ(μ ~ E) < 0.
   - The physical (unitary-gauge) counting is the more plausible one. But the forward-scattering potential off a
     non-stationary h̄ is gauge-dependent, and the track's safety statement rested on this term alone.
   - The N gas (N2) and, after thermalisation, the standard thermal mass (κ_T ≈ 0.11 > 0) do not have this issue.

4. **N4 [C] The medium-mass margin depends on the energy law** (`sk_medium.json`).
   - With injection energies held fixed, the SM term alone at ×10⁻⁴ fails at ≥ 174.05 GeV (D = +2.6 / +25 on the
     export; +44 / +53 Pauli).
   - N + SM is safe to ×(2–6)×10⁻⁴. With redshifting energies it is safe to ×10⁻⁶.
   - The track's "c × 10⁻⁴: D ≤ −107" holds only with its redshifting E_SM. The conclusion stands, with a margin of
     ≳ 10³ rather than 10⁴.

5. **N5 [D/C] N_s is not forced to be fatal** (`sk_ns.json`).
   - **Threshold in m̃_s.** At M_soft = 0.3 / 1 / 1.5 H_i the inflationary tolerance (−0.66 H_i² at 172.86) holds for
     m̃_s ≤ 0.048 / 0.018 / 0.013 eV, even under "vanish at M_soft".
   - **What N_s needs.** It only needs m̃_s ≳ 7×10⁻⁶ / 6×10⁻⁷ / 3×10⁻⁷ eV to decay before it would dominate (it is
     produced at a ≈ 1.8×10³ with 1.3% of ρ, matter-like afterwards).
   - **No seesaw at all.** It can also stay dark radiation, ΔN_eff ≈ 0.04 [Q mf_hybrid skeptic], which is retrospectively
     inside Planck's N_eff. The Planck number is recalled, not re-searched here.
   - **Σm_ν.** Five singlets at 0.05 eV each give Σm_ν = 0.25 eV for a real Casas–Ibarra R. That is already above the
     quoted Planck+BAO 0.12 eV and DESI DR2 0.064 eV [Q mf_hybrid skeptic]. So 0.05 eV for N_s is the least natural choice.
   - **CP bookkeeping.** M(s) = |M_soft + y_s s| contains an s-odd piece. The inflationary Higgs mass then differs between
     the two CP branches (sign s_*) by 0.003 / 0.031 / 0.065 H_i² at M_soft = 0.3 / 1 / 1.5 H_i (≤ 10% of the tolerance).
     It is small, but it is a third CP-correlated input, as in inf_higgs skeptic N10.

6. **N6 [C] The "[N = 60]" bracketed crossings are a coarse-grid artefact.**
   - `Mt_crit_grav_N60` is interpolated on the 0.29 GeV MT_EXT grid between D(172.86) and D(173.15).
   - That puts the N = 60 crossing above the N_obs one, the wrong way round, since N = 60 > N_obs means more regions.
   - On a 0.02 GeV grid: 172.870 (export) and 172.906 (Pauli), against the track's 172.885 / 172.937. This is immaterial.

7. **N7 [C, with Q] At the self-consistent Pauli H_i the gravitational-only verdict fails at +1σ.**
   - On mf_chain's own Pauli history at 2.40×10¹¹ GeV: D(172.86) = +34.5 and M_t,crit = 172.74 (band 172.71–172.82).
   - At α_s = 0.1171, central M_t fails: D = +40.8, M_t,crit = 172.41.
   - mf_chain also finds HFa's clock is not CMB-flat there at κ = 10 (capacity 48.6 < N_* = 52.3) [Q, unreviewed].
   - So, on the more realistic closure, survival above central M_t rests entirely on the medium masses (N2, H8).

8. **N8 [D/C] G6 rebuild.**
   - The harmonic rebuild underestimates the period-mean positive mass (by up to 72% at A/v = 0.4). It underestimates
     the depth of the negative swing (about 20%) and overestimates the line frequency by 10–20%.
   - On the Pauli history the condensate (f_χ ≈ 0.85 to τ ~ 30) keeps the first band active: per-mode growth 0.5 at
     τ = 3 for Λ = 10⁴.
   - The upper edge is uncertain to the factor-3 grid. The lower edge is unaffected.

9. **N9 [C] h_max quantisation.**
   - The shared SMTab finds h_max on a Δ ln h = 1.15% grid. That gives a D jitter of up to ±2 and ±0.003 GeV in M_t,crit.
   - It is immaterial, but it is the main difference between my D and the track's at 172.86.

10. **N10 [R/A] Thermal tunnelling at T ~ 10¹³–10¹⁴ GeV is not computed by the track.**
    - Salvio–Strumia–Tetradis–Urbano (arXiv:1608.02555 [R]) find only minor thermal corrections to SM vacuum decay. The
      higgs_reheat skeptic found S₃/T ≥ 478 [Q].
    - With m_ind² ≫ T²-scale terms and φ/T ≪ 1, I see no route to thermal decay here. It is not computed.

11. **N11 [A, favourable] The SU(2) anomaly couplings of the rotating axial phases** (mf_hybrid skeptic N6) could add a
    gauge-boson population, which would add a positive Higgs mass. It is not computed.

## 5. Omitted-effects checklist

| effect | status |
|---|---|
| preheating of the Higgs by ξR, and the R lines | none: q ≤ 0.05 even at ξ = 1000; window exponents ≤ 0.0024 [C] |
| tachyonic/broad resonance from a portal | G6 only (hypothetical); edge Λ ≈ 3×10³–10⁴ [C] |
| waterfall fragmentation (mf_hybrid skeptic N2) | not computed (no lattice). The media are produced on scales 1/m_χ ≪ 1/H, so super-Hubble regions see the mean; walls would raise ⟨R⟩ (w → −2/3), which helps the Higgs [A] |
| thermal / medium corrections | the N gas (new, positive), the SM hard gas (positive in the physical counting), the thermal mass (positive); thermal tunnelling minor [R] |
| induced couplings | axial N: the vacuum loop is constant [D], two-loop mismatch ≲ 0.02 H_i² [Q]; N_s: conditional (N5); no hidden–Higgs portal (sequestered by choice) |
| gravitational production | inside the mode equations; upper edge reproduced [C] |
| isocurvature | the Higgs has m² = 1.05 H² during inflation, so its super-Hubble modes decay ∝ a^{−0.4} and carry no frozen isocurvature [D] |
| decays into the hidden sector | none for H, N₁–N₄ or the medium (no portal, by choice) [A] |
| kinematics, perturbativity | N → L H open (M_N ≫ m_H(h̄)); y_ν ≤ 0.29; G6 \|λ_χH\| ≤ 2×10⁻⁸ ≪ √(λ_χλ_H) [D] |
| CP-sign bookkeeping | the Higgs criteria are branch-blind except the N_s s-odd piece (N5) |
| M1r assumptions | ξ_H = 0.1 at M_t; H_i design value vs mf_chain calibration 0.73–3.9×10¹¹ [Q]; zero SM portals (true of the bench) |
| CMB numbers | enter only through N_*, reproduced. Planck 2018 / BICEP-Keck values are not used by the Higgs criterion |

## 6. Corrected verdict for ξ_EW = 0.1 (model choices as in the track)

- **During inflation.** Safe. The constant-H edge is ξ_min = 0.049 / 0.053 (tree; 0.041–0.048 in the other three
  schemes), and the cap is above 177 GeV (extrapolated). The non-equilibrium Gaussian count gives ≤ 10⁻¹⁵⁸ at 174.6. The
  instantaneous reading gives no usable bound (§3.1).
- **After inflation, gravitational only.** Reproduced and marginal:
  - M_t,crit = 172.88 (unblocked) to 172.92 (Pauli) at 1.3×10¹¹ GeV, i.e. +1.1σ to +1.2σ;
  - at the self-consistent Pauli H_i = 2.4×10¹¹ it is 172.74, failing at +1σ;
  - α_s ∓ 1σ moves it by ∓0.33 GeV.
- **With the medium.** Safe to at least 174.6 GeV on every history, with a margin of 10³–10⁶ in the medium-mass
  normalisation. The rescue now rests on two independent positive terms: the N gas from the end of inflation, and the SM
  hard gas from τ ≈ 0.2. The N-gas term is [E] (hard, non-thermal, Pauli-blocked quanta) but sign-robust.
- **ξ window at central M_t (+1σ).**
  - Lower edge: 0.049 (0.053), from inflation. Gravitational-only after inflation it is 0.067 (0.098) on the export
    and 0.066 (0.096) on the Pauli rebuild (recomputed).
  - Upper edge (mean-field): ξ_EW ≈ 111 (21).
  - ξ_H = 0.1 lies inside; with the medium, by about ×2 above the lower edge.

## 7. Not established here

- No lattice. The waterfall fragmentation and walls of mf_hybrid skeptic N2 are not computed.
- The N-gas and SM medium masses are hard-gas estimates without back-reaction. Thermalisation would only raise them.
- The SM input is shared (smpot). Above 174.6 GeV it is a quadratic extrapolation of the SMDR grid.
- The input-edition shift (inf_higgs skeptic N7) and the MC-versus-pole top mass are not recomputed.
- The G6 runaway was not rerun with modes; only Floquet and virial checks were done.
- The N_s sector is not built.
- The ξ upper edge was checked only on the export and only to ξ(H_i) = 300.

## 8. Files and reproduction (one core, ≈ 25 min in total)

All under `research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_higgs/skeptic/`. Mode runs and SM tables
are cached in `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_higgs_skeptic/`.

| code | content → results |
|---|---|
| `sk_common.py` | paths, smpot access (read-only), own flat-potential table (h_max refined, V″ canonical), history loader |
| `sk_infl.py main / cap / inst` | Hawking–Moss in four schemes, light cone, caps, instantaneous count with and without the endpoint → `sk_infl.json` |
| `sk_modes.py runs` | own u_k modes on the export and mf_chain's rebuilds; own D (scratch npz) |
| `sk_postD.py` | D, M_t,crit, band, N = 60, α_s, M_h, deadlines, medium terms in the criterion → `sk_postD.json` (+ controls) |
| `sk_medium.py` | the medium-mass factor needed for safety at 174.6, two energy laws → `sk_medium.json` |
| `sk_plasma.py` | N-gas and SM hard-gas masses on three histories → `sk_plasma.json` |
| `sk_windows.py` | R < 0 windows, R line frequency and q → `sk_windows.json` |
| `sk_upper.py` | upper ξ edge → `sk_upper.json` |
| `sk_g6.py` | G6 virial and Floquet → `sk_g6.json` |
| `sk_ns.py` | N_s versus m̃_s → `sk_ns.json` |
| `sk_summary.py` | → `sk_summary.json` |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_higgs/skeptic/code; export OMP_NUM_THREADS=1
python3 -B sk_infl.py main; python3 -B sk_infl.py cap; python3 -B sk_infl.py inst
python3 -B sk_modes.py runs; python3 -B sk_postD.py; python3 -B sk_medium.py; python3 -B sk_plasma.py
python3 -B sk_windows.py; python3 -B sk_upper.py; python3 -B sk_g6.py; python3 -B sk_ns.py; python3 -B sk_summary.py
```

The ξ_min-with-endpoint, B-past-ε_H = 1, Gaussian-tail and mode-control entries were added by short inline drivers. They
use only the functions above (`sk_infl.inst_log10N`/`Bexp`, `sk_modes.run`/`danger`).

## Literature (checked by web search in this session)

- **Inflation-stage criterion.** Mantziris–Markkanen–Rajantie, arXiv:2011.03763 (bubble count; ξ_EW ≳ 0.06 in field
  models, ≳ 0.1 for Starobinsky).
- **Curved-space potential.** Markkanen–Nurmi–Rajantie–Stopyra, arXiv:1804.02020 (one-loop SM potential in curved space;
  scheme C/D here).
- **Production at the end of inflation.** Herranen–Markkanen–Nurmi–Rajantie, arXiv:1506.04065 (significant production for
  ξ ≳ 1).
- **Lattice upper bound.** Figueroa–Rajantie–Torrentí, arXiv:1709.00398 (instability for ξ ≳ 4–5 at their M_t; gauge
  fields have minimal impact).
- **Preheating bounds.** Ema–Mukaida–Nakayama, arXiv:1602.00483 (upper bounds on Higgs–inflaton and Higgs–curvature
  couplings).
- **Thermal decay.** Salvio–Strumia–Tetradis–Urbano, arXiv:1608.02555 (minor thermal corrections to SM vacuum decay).
- **SM inputs.** PDG 2024 top quark listing: M_t = 172.57 ± 0.29 GeV (S = 1.5). α_s = 0.1180 ± 0.0009 and
  M_h = 125.20 ± 0.11 GeV are confirmed through secondary sources quoting PDG 2024.
- **Not confirmed:** arXiv:1602.07719 (the brief's list), as in the track. It is not used.
- **Quoted through the tracks** [Q]: Planck 2018 X (1807.06211), Planck 2018 VI (1807.06209), DESI DR2 (2503.14744), via the
  mf_hybrid skeptic.

Our tracks, read-only: mf_higgs (code, results), inf_higgs and its skeptic, higgs_reheat (smpot), mf_hybrid and its
skeptic (bench, tables), and mf_chain (rebuilt histories, the H_i iteration and CMB rescale; unreviewed).
