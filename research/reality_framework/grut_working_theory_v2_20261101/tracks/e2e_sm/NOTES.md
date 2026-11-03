# e2e_sm: the Standard Model side of M2

Track of the end-to-end pass ([BRIEF_E2E.md](../../BRIEF_E2E.md), section e2e_sm). Every choice below is a model choice. Closure fixes
none of these numbers. Nothing here is a GRUT prediction, a validation or an exclusion. Data comparisons are
retrospective. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Inputs reused read-only: the mf_hybrid singlet sector (`mf_hybrid/results/bench_HFa.json`); the e2e_kinetic history
`e2e_kinetic/results/fluid_HFa.json` / `fluid_HFb.json` (in progress at the time of writing: K = 16.3, T_RH = 3.06×10¹² GeV,
H_i = 1.3×10¹¹ GeV, mf_hybrid's m̃ = 0.05 eV); the mf_chain spectator numbers; the mf_higgs skeptic's N-gas medium
mass. No code was copied. All new code is in `code/`.

## 0. Verdict

| item | result | label |
|---|---|---|
| Light ν | NO, m1 = 0 (choice), NuFIT 6.0 Δm² fitted. **Σm_ν = 0.0588 eV**, below DESI DR2+CMB 0.0642 eV and Planck+BAO 0.12 eV | D/Q |
| Seesaw (S1) | complex Casas–Ibarra R, 3×4. m̃ = (0.050, 0.091, 0.050, 0.050) eV keeps every N rate and the N-gas Higgs mass at or above the mf_hybrid / mf_higgs-skeptic values. Cost: a ×4.1 cancellation in m_ν (max\|R\| = 1.35). One-loop non-universal shift ≤ 1.4% of m3 | C |
| New constraint | a CP-conserving R exists only if Σ_j η_j m̃_j = Σ_i m_i, with η_j = ±1 the N CP parities. "Four singlets at 0.05 eV" is CP-conserving only with N2's m̃ raised to 0.091 eV | D |
| N_s | **massless dark radiation** (no soft mass, no seesaw Yukawa; N_s → iγ5N_s exact). **ΔN_eff = 0.001–0.004** on the e2e_kinetic history | C |
| UV of the axial couplings | Φ_a = (ρ_a/√2)e^{iθ_a}, with σ = f₁θ₁ and χ = f₂θ₂ as angles, and the hybrid potential as a monodromy potential. One-loop flat in θ; the CP bookkeeping is unchanged; no new light states; no strings if U(1)_a is never restored (T_max/f₁ = 0.35). **But a weakly coupled completion is not frozen at the slice amplitude:** θ̇/m_ρ = 1.5–2.2 at λ = 4π; it needs λ ≥ 110 (σ) and 248 (χ). At strong coupling (m_ρ = 4πf) the ratio is 0.6–0.9 | C/A |
| Anomalous SU(2) couplings | In the mass-phase basis there is no θWW̃. The (α₂N_g/8π)θWW̃ of the mf_hybrid skeptic (N6) comes with an equal derivative coupling to J_L. The physical θ dependence is only through the Majorana mass. No tree-level W preheating | A |
| η_B | **overproduced at a generic CP phase:** η_B = 1.7×10⁻⁴ δ (δ the explicit CP phase of R). The observed 6.1×10⁻¹⁰ needs **δ ≈ 3.5×10⁻⁶** (window 1.6–7.8×10⁻⁶), i.e. \|ε₄\| ≈ 2×10⁻⁸ | C |
| Source | non-thermal N4 from the late χ-remnant decays (χ → N4N4) at T ≈ T_RH ≈ M₄/5. Thermal leptogenesis alone would give 3.2×10⁻⁷ δ (δ ≈ 1.9×10⁻³) | C |
| Spontaneous lepto/baryogenesis | erased: an asymmetry present at t ≤ 30 H_i⁻¹ survives the washout by ≤ 10⁻⁴⁵ | C |
| BBN; long-lived states | every new massive state decays by ~10⁻²⁹ s; N_s is the only new light state | C |

## 1. Neutrinos (`code/sm_seesaw.py` → `results/seesaw_H*.json`)

**Light sector.**
- NuFIT 6.0, IC24 with SK-atm, normal ordering (arXiv:2410.05380, Table 1, verified):
  - Δm²₂₁ = 7.49 ± 0.19 ×10⁻⁵ eV²; Δm²₃₁ = 2.513 (+0.021/−0.019) ×10⁻³ eV²;
  - θ12 = 33.68°, θ13 = 8.56°, θ23 = 43.3°, δ = 212°.
- m1 = 0 is a choice. Then (m1, m2, m3) = (0, 8.65, 50.13) meV and **Σm_ν = 0.0588 eV**. The window m1 ≲ 4.2 meV keeps
  Σ < 0.0642 eV.
- Retrospective bounds:
  - DESI DR2 BAO + CMB, ΛCDM: Σm_ν < 0.0642 eV (95%) (arXiv:2503.14744, verified). Later DESI combinations quote 0.0565–0.059 eV (arXiv:2607.27410, seen in search), at the NO floor.
  - Planck 2018+BAO: < 0.12 eV (arXiv:1807.06209).
  - Inverted ordering needs Σ ≥ 0.0995 eV, so IO is not used.

**Heavy sector.** mf_hybrid's four axial Majorana singlets: M_j/H_i = 250 (N1, from σ) and 400, 200, 120 (N2–N4, from χ).
In GeV this is 1.56–5.2×10¹³ at H_i = 1.3×10¹¹. Type-I seesaw with Casas–Ibarra (hep-ph/0103065):
- m_D = i U* √m R √M, with R a complex 3×4 matrix and R Rᵀ = 1₃.
- Unflavoured, h†h = √M R† m R √M / v² does not depend on U [D].
- m̃_j = Σ_i m_i |R_ij|²; Γ_j = (h†h)_jj M_j/8π; K_j = m̃_j/m*, with m* = 1.08×10⁻³ eV (hep-ph/0401240).

**[D] A constraint the brief's "complex R" route meets.**
- If R is CP-conserving (each column real or imaginary, with CP parity η_j = ±1), then R Rᵀ = 1 gives
  **Σ_j η_j m̃_j = Σ_i m_i**.
- With four equal m̃ = 0.05 eV this has no solution (0.05 × integer ≠ 0.0588). So "four singlets at 0.05 eV with
  Σm_ν = 0.0588" is unavoidably CP-violating, with ε_j ~ 10⁻³.
- I instead choose η = (+, −, +, +) and **m̃ = (0.050, 0.0912, 0.050, 0.050) eV** (scenario S1).
  - R0 is the least-norm solution, max|R| = 1.35.
  - CP is then broken by one explicit phase: R(δ) = R0 exp(δA), with A complex antisymmetric and CP-odd. R Rᵀ = 1
    holds exactly, so the light spectrum is exact for every δ.
  - Results: ε_j = δ × (2.07, 0.32, 0.96, 0.67)×10⁻² (N1..N4) at H_i = 1.3×10¹¹. m̃ moves only at O(δ²): 0.03% at
    δ = 10⁻².

**Costs and checks** [C/E]:
- **Tuning.** Σm̃/Σm = 4.1, i.e. a ×4 cancellation in m_ν.
- **One-loop (Z, h) correction.** Grimus–Lavoura type, heavy-M limit [E]. It has two parts:
  - a universal 22%, absorbed in the fit;
  - a non-universal part ≤ 6.9×10⁻⁴ eV (1.4% of m3, 0.2% of m2), which R must be re-tuned to absorb.
- **Reheating rates** (Γ_rest/H_i = 0.53, 2.49, 0.34, 0.12): every one is at or above mf_hybrid's. N2's rate is ×1.82.
  **e2e_kinetic should rerun with m̃₂ = 0.0912 eV**; as quoted, its history used 0.05 eV for all four.
- **Higgs medium.** The N-gas Higgs medium mass m² = Σ y_j²ρ_j/E_j² uses y_j² = (0.054, 0.156, 0.043, 0.026). That is at
  or above the mf_higgs skeptic's smallest-Yukawa value (0.026), so the skeptic's rescue stands [D].
- **Alternative S2** (m̃ = 0.03, 0.03, 0.03, 0.0288 eV; tuning ×2.0): rates ×0.6, medium mass ×0.6, within the skeptic's
  ≳ 10³ margin. Its leptogenesis result (δ ≈ 3.3×10⁻⁶) uses the S1 history, so it is indicative only.

**N_s.** Fixed as massless dark radiation. It has no soft Majorana mass and no LH Yukawa, so N_s → iγ5N_s is exact and the
holder stays the only s-odd term [D].
- It has no SM coupling, hence no BBN or Σm_ν effect.
- **ΔN_eff (computed on the e2e_kinetic history).** s is matter-like through the χ-remnant era and grows ∝ a only after
  reheating. It decays at t = 1/Γ_s = 7.7×10⁵ H_i⁻¹ (a ≈ 5×10³), so
  ΔN_eff = 2.87 ρ_DR/ρ_SM = **9×10⁻⁴ (s_* = 3×10⁻³ M_P) to 3.4×10⁻³ (s_* = 1.4×10¹⁶ GeV), 4.3×10⁻³ at H_i = 8.4×10¹⁰**.
- That is below the retrospective N_eff = 2.99 ± 0.17 (Planck 2018, 1807.06209 [R]) and 2.86 ± 0.13 (ACT DR6,
  arXiv:2503.14454, verified).
- The skeptics' 0.04 assumed a 1.3% s fraction at decay on the unblocked history. On mf_chain's Pauli history (0.1–0.2%)
  the same formula gives 0.002–0.005.
- **Variant (not adopted):** M_soft = 1 H_i with m̃_s = 10⁻⁵ eV, inside the mf_higgs-skeptic window. N_s then decays at
  1.7×10⁻⁹ H_i, out of equilibrium, with η_B,max ≈ 4×10⁻⁸. That needs a second CP phase ≲ 0.015 and adds an s-odd piece [E].

## 2. UV completion of the axial couplings (`code/sm_checks.py` → `results/checks.json` "uv")

**Construction [A].**
- Two complex scalars Φ_a = (ρ_a/√2) e^{iθ_a}, with ⟨ρ_a⟩ = f_a (f₁ = 1200, f₂ = 2000 H_i). Yukawas y_j Φ_a N_jN_j give
  M_j = y_j f_a/√2.
- **The phase fields are σ ≡ f₁θ₁ and χ ≡ f₂θ₂ themselves.** The fields that couple to N are the angular modes, with no
  separate Goldstone.
- The non-periodic hybrid potential V(σ, χ) is then an explicit breaking of U(1)₁×U(1)₂ by a multi-branched (monodromy)
  potential. Examples are the Kaloper–Sorbo 4-form coupling (arXiv:0811.1989) and axion monodromy (arXiv:0808.0706). Both
  arXiv numbers were confirmed through secondary references.
- The phases wind 10²–10³ times (v/f₂ = 183.5 rad).

**Checks.**
- **One-loop flatness [D].** |M_j| = y_j ρ_a/√2 does not depend on θ_a, so N loops give no θ potential, and radial
  loops couple to θ only through ρ²(∂θ)².
  - The Yukawas break U(1)₁×U(1)₂ to the diagonal (lepton number). The relative angle σ/f₁ − χ/f₂ gets mf_hybrid's
    phase-mismatch potential (multi-loop, ≲ 10⁻⁵ of V′) [Q].
- **CP bookkeeping [D].** Nothing new couples to s. The holder remains the only s-odd term, and the DM sign is still
  α + sign(s_*).
  - The vacuum phase of N2–4 (v/f₂ mod 2π = 1.31 rad in the +v vacuum, its negative in −v) is part of the leptogenesis
    phase δ.
  - **So the baryon sign is a separate supplied input, uncorrelated with the dark-matter sign** [D/A].
- **New states [C/A].**
  - The radial modes are heavy (m_ρ = √(2λ) f ~ 10¹⁴–10¹⁵ GeV) and decay to NN at once.
  - There is no new light Goldstone.
  - No strings form: U(1)_a is broken before inflation (f ≫ H_i) and not thermally restored (T_max = 5.5×10¹³ GeV, against
    f₁ = 1.56×10¹⁴; restoration needs T ≳ √3 f).
  - The hybrid's own ±v walls map onto θ₂ = ±v/f₂ branches. They are e2e_lattice's question, and the completion adds none.
- **[C] Radial test: the decisive negative.**
  - The fast rotation stretches ρ centrifugally: ρ² = f² + θ̇²/λ.
  - On the e2e_kinetic history θ̇/f peaks at **7.4 (σ) and 11.1 (χ)**.
  - Freezing ρ (m_ρ ≥ 2θ̇) needs **λ ≥ 110 (σ) and 248 (χ)**, far beyond perturbative λ ≲ 4π.
  - At λ = 4π: ρ²/f² = 2.7–3.7 at the peak, so μ5 is cut ×1.6–1.9, and θ̇ < m_ρ/2 is reached only at t ≈ 1.4–1.9 H_i⁻¹,
    after the main non-perturbative transfer (f_rad = 0.36 by t ≈ 0.45).
  - At strong coupling (m_ρ = 4πf) θ̇/m_ρ = 0.59–0.89: the EFT edge mf_hybrid already flagged, restated.
  - **What fixes it:** raise f by ×3 (σ) and ×4.4 (χ) for a weakly coupled completion. That cuts μ5,max, hence the
    non-perturbative rates, and raises K, so e2e_kinetic must rerun. Otherwise accept a strongly coupled completion at its
    cutoff.
- **Anomalous SU(2) couplings [A].**
  - In the Lagrangian as written, θ enters only via M e^{iθ}NN with N a gauge singlet. Rotating N alone (anomaly-free)
    moves θ into the Yukawa.
  - The skeptic's (θ/2)N_g(α₂/8π)WW̃ appears only after an anomalous rotation of L. That rotation also generates
    (∂θ/2)J_L, and by parts the two cancel up to the M-proportional L violation.
  - So the physical θ–SM coupling needs both y and M, and there is no unsuppressed tachyonic W growth (κ ≈ 10² H_i is a
    basis artefact). The physical remnant is spontaneous leptogenesis through the ΔL processes (Domcke–Ema–Mukaida–Yamada,
    arXiv:2006.03148, verified), which is washed out (§3).
  - A skeptic should check this basis argument against the invariant-coupling treatment of 2006.03148.

## 3. Baryogenesis (`code/sm_lepto.py`, `code/run_lepto.py` → `results/lepto.json`)

**Equations.** Unflavoured, in cosmic time, comoving densities.
- N_j^th relaxes to equilibrium at Γ_j⟨K₁/K₂⟩.
- Source: Σ_j ε_j[Γ_j⟨K₁/K₂⟩(N^th − N^eq) + Γ_j(M_j/E_j)N_j^NT].
- Washout: inverse decays Γ_W = Γ_j z²K₁(z)/4, plus ΔL = 2 at 2.62×10⁻²⁹ T³ Σm² (BDP hep-ph/0401240 forms, verified)
  [Q].
- N_j^NT, the axially produced and χ-remnant N, comes from e2e_kinetic's f_Nj with E_j = max(¾A_j m_par, M_j, m_par/2)
  [E]. T is e2e_kinetic's thermalised bath.
- After the history ends: radiation. Y_B = (28/79)Y_{B−L} [R, Harvey–Turner]; η_B = 7.04 Y_B.
- The integrator is exponential (exact per step) on 6×10⁴ log-steps.
- **Check:** with the non-thermal source switched off, it gives the textbook strong-washout result η_B ≈ 0.96×10⁻² ε₄κ_f
  (K = 46): 3.2×10⁻⁷ δ.

**Results** [C on Q/E rates], η_obs = 6.12×10⁻¹⁰ (Planck 2018 ω_b = 0.02237 [R]):

| case | η_B / δ | thermal only / δ | δ needed |
|---|---|---|---|
| S1, HFa (and HFb: same history to 10⁻⁴), H_i = 1.3×10¹¹ | 1.74×10⁻⁴ | 3.2×10⁻⁷ | **3.5×10⁻⁶** |
| washout ×0.5 / ×2 | 3.8 / 0.78 ×10⁻⁴ | 7.2 / 1.4 ×10⁻⁷ | 1.6 / 7.8 ×10⁻⁶ |
| no non-thermal N | 3.2×10⁻⁷ | same | 1.9×10⁻³ |
| rescaled to H_i = 8.4×10¹⁰ (same history in H_i units) | 2.9×10⁻⁴ | 3.6×10⁻⁷ | 2.1×10⁻⁶ |
| S2 (history not consistent) | 1.9×10⁻⁴ | 6.1×10⁻⁷ | 3.3×10⁻⁶ |

**Reading.**
- On this history reheating ends late, at t ≈ 4400 H_i⁻¹ and T_RH = 3.1×10¹² GeV ≈ M₄/5. The last of the χ remnant decays
  to N4 pairs near rest while the inverse-decay washout is switching off.
- Survival of an asymmetry present at t: 1.3×10⁻⁹ at t = 3000, 0.048 at 10⁴ and 0.90 at 3×10⁴. So non-thermal N4
  leptogenesis dominates over the thermal yield by ×540.
- **η_B is overproduced at any natural phase:** ×3×10⁴ at δ = 0.1.
- **It is reached only for δ ≈ (1.6–7.8)×10⁻⁶,** i.e. |ε₄| ≈ 2×10⁻⁸.
- δ includes the vacuum phase v/f₂ (1.31 rad), so the Yukawa phases must cancel it to ~10⁻⁶. That is a tuning [A].
- What fixes η_B otherwise:
  - (i) the supplied CP phase;
  - (ii) a higher T_RH/M₄, with full washout, which takes it to the thermal value (needs δ ~ 2×10⁻³);
  - (iii) a χ-remnant channel into an N with tiny ε, or flavour alignment;
  - (iv) a late entropy release.
- The result depends exponentially on the washout near t ~ 5×10³ and on e2e_kinetic's final history. ΔL = 1 scatterings,
  thermal masses, spectator processes and flavour (T ≈ 2–3×10¹² GeV is at the τ-equilibration edge) are not included:
  O(1)–×10 in δ [E].
- **Spontaneous lepto/baryogenesis** from the rotating phases. θ̇ oscillates in sign, at ω = m_σ, m_χ ≫ Γ_ΔL, and
  sphalerons are out of equilibrium at T ≳ 10¹² GeV (Γ_sph/H = 0.2 at 2.6×10¹² GeV [E]). Whatever is generated before
  t ≈ 30 H_i⁻¹ survives by ≤ 10⁻⁴⁵ [C]. It is irrelevant.

## 4. Other checks (`results/checks.json`)

- **Lifetimes.** N: 2–4×10⁻³⁵ s. s: 3.9×10⁻³⁰ s. Radial modes: shorter. Nothing new is present at BBN except N_s
  (ΔN_eff ≤ 0.004) [C].
- **Sphaleron conversion.** B−L is made at T ≈ 2–3×10¹² GeV and conserved by sphalerons. Conversion (28/79) applies once
  they equilibrate, at T ≲ 10¹² GeV [R/E].
- **Washout.** Computed above. Strong for all N: K = 46–85 (S1) [C].
- **Gravitational anomaly** of the N rotations: θRR̃ vanishes on FRW [D].

## 5. Open, and dependencies for other tracks

- **e2e_kinetic:** rerun with m̃ = (0.05, 0.0912, 0.05, 0.05) eV. If the weakly coupled UV is wanted, rerun with f ×3/×4.4.
  Its final T_RH and remnant-decay timing set δ.
- **e2e_chain:** δ is one more supplied CP input, independent of the holder phase. N_s gives ΔN_eff ≈ 10⁻³.
- **Not done:** flavoured Boltzmann equations; ΔL = 1 and thermal corrections; the strong-coupling UV (composite N); the
  parametric radial excitation and possible non-thermal restoration; the soft-N_s variant beyond an estimate.

## Files and reproduction (one core, < 1 min)

```
cd code
python3 -B sm_seesaw.py 1.3e11 ../results/seesaw_H1.3e11.json
python3 -B sm_seesaw.py 8.4e10 ../results/seesaw_H8.4e10.json
python3 -B run_lepto.py          # -> results/lepto.json (reads ../../e2e_kinetic/results/fluid_HF*.json)
python3 -B sm_checks.py          # -> results/checks.json
python3 -B make_card.py          # -> results/sm_card.json (every SM-side parameter with its status)
```

## References (arXiv numbers verified by search this pass unless marked [R])

- NuFIT 6.0, arXiv:2410.05380 (Table 1 read).
- DESI DR2 neutrino constraints, arXiv:2503.14744.
- Planck 2018 VI, arXiv:1807.06209 (ω_b, N_eff values [R]).
- ACT DR6 extended models, arXiv:2503.14454.
- Casas–Ibarra, hep-ph/0103065 [R].
- Covi–Roulet–Vissani, hep-ph/9605319.
- Davidson–Ibarra, hep-ph/0202239.
- Buchmüller–Di Bari–Plümacher, hep-ph/0401240.
- Domcke–Ema–Mukaida–Yamada, arXiv:2006.03148.
- Kaloper–Sorbo, arXiv:0811.1989.
- McAllister–Silverstein–Westphal, arXiv:0808.0706.
- Fields–Olive–Yeh–Young, arXiv:1912.01132 (BBN η, context only).
- Harvey–Turner 1990 (28/79), no arXiv [R].
