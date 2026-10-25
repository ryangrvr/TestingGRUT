# Track sp_yledger: the Y-debris ledger of Codex 373d37e, and what it leaves open

25 October 2026 · Claude (independent review track, issue #13) · packet `373d37e`
(`origin/codex/grut-inflation-scale-scan-20261025`, `research/reality_framework/grut_inflation_scan_20261025/`),
read only through `git archive 373d37e … | tar -x` into the scratchpad. Our own code is in `code/`, outputs in
`results/`. No Codex code was run, and nothing was committed. Our earlier gp_leftovers code is reused for the
thermal benchmark, as permitted by the brief.

**This is one supplied tree-level working model. Nothing here is evidence that X, Y, Σ, φ, the vector or the
inflationary inputs exist. It is not a validation, an exclusion or a GRUT prediction.** All agents involved are
same-session AI, so this is not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Inputs [Q]

- Couplings (unchanged): M_Y = 1e12, h = 1e11, μ = 5e11 GeV, κ = 1e-3, e = 1e-7, λ_Y = λ_Σ = 1, λ_X = 1e-3.
- Light target: g_X = 0.1, λ_φ = 0.03, v = 0.924089637 GeV, m_X = 0.092409 GeV, Q = 2e5.
- Calibration: Y_Q = 4.8294e-9, B = 4.4e-10 GeV per visible entropy.
- Leading point: H_i = 1.1211167304108865e13 GeV, c_H = 1. Aligned inflationary minimum with
  Σ = −|Σ| e^{iπ/2}, zero velocities, instant removal of the Hubble mass, radiation, g* = 106.75.
- Codex targets (CLAUDE_HANDOFF.md, Y_SOURCE_RESULTS.json):
  - release: ε_X = −0.24066124, N_X = −238.64868, E_X = 741.64988, J_Y = 13.26903, Y_X = −65.9823 Y_Q;
  - debris: Γ_Y = 377.95 GeV, a_dec = 1.722e5, E_Y/E_X = 274.85 at decay;
  - pairs: 0.08340 per net charge against an allowance of 0.0049296, a ratio of 16.92;
  - energy-equivalent fraction 0.003555.

## 2. Method (own code)

- `code/yrelease.c` [C]. Homogeneous X, Y, Σ release in cosmic time with physical fields, integrated with an
  adaptive Dormand–Prince 5(4) scheme. Codex instead uses conformal fields with DOP853.
  - Backgrounds: radiation (bg 0); an inflaton matter era with decay width Γ_φ set by T_R (bg 1); and bg 1 with
    the Hubble mass kept as −c_H(ρ_φ/3)|X|² (bg 2).
  - Diagnostics: charges, E_X, ε, and an energy ledger W = a³ρ integrated as an ODE.
  - Y excitation, measured about the static forced solution M²y_p + 2μΣy_p* = −κX³ (with its exact time
    derivative). It is reported two ways: Codex's single-mass number J_Y = E_exc/M, and our projection onto the two
    real Majorana eigenmodes, m²_{H,L} = M² ± 2μ|Σ|, in the frame that co-rotates with the Σ phase,
    N_Y = Σ_j E_j/ω_j.
  - Σ radial number N_σ = a³(ṙ² + λ_Σ(r² − h²/2)²)/m_σ.
  - Window means are weighted by da, as in Codex's windows.
- `code/ymodel.py`: parameters, the multi-start minimum polished by a Newton root of the gradient (|∇V| =
  1.2e-15), and the run wrapper.
- `code/point_ledger.py` [C, D, E]: the leading point to a = 3000, with tolerance and preparation controls; tree
  widths; decay epochs; pairs; and deposited energies.
- `code/rates.py` [E]: the Γ/H table for every channel in item (b).
- `code/thermal_relic.py` [C on the gp_leftovers benchmark]: the maximal-equilibrium relic with the debris energy
  added to the single onset budget.
- `code/md_release.py`, `code/md_leading_TR.py` [C]: the inflaton-matter release scan, and its T_R dependence at
  the leading point, for item (c).
- `code/phase_spread.py` [C, side check]: the inflationary spread of the source phase, and an ensemble run at
  that spread.

## 3. Results

### 3a. Reproduction of the leading point [C, own code]

| Quantity (window a = 600–1000 unless stated) | Codex 373d37e | This track | Verdict |
|---|---:|---:|---|
| Inflationary minimum \|X\|, \|Y\|, \|Σ\| / H_i | 23.826064299, 1.893761922, 0.430868093 | 23.826064299, 1.893761922, 0.430868093 (\|∇V\| = 1.2e-15) | CONFIRMED |
| ε_X | −0.24066124 | −0.24066124 (rtol 1e-11/1e-12/1e-13 agree to 2.5e-8) | CONFIRMED |
| N_X, E_X | −238.64868, 741.64988 | −238.64867, 741.64988 | CONFIRMED |
| J_Y (single-mass convention) | 13.26903 | 13.26903 | CONFIRMED |
| Y_X / Y_Q | −65.98 | −65.98 (same N_X, same instant-reheating entropy) | CONFIRMED |
| δ → δ + 1e-9 rad | Δε ≈ 2.6e-7 | Δε = 2.57e-7 | CONFIRMED |
| Late drift (a = 2000–3000 mean versus 600–1000) | ≈ 0.00046 | 0.00046 (ε = −0.2402016) | CONFIRMED |

- **State agreement.** Final states agree with Codex's published canonical states:
  - at a = 1000 (Codex rtol 2e-14, ours 1e-13), to 4e-8 (X), 3e-7 (Y) and 5e-7 (Σ);
  - at a = 3000, to 6e-6 to 3e-5.
- **New control: radial preparation.** Multiplying |Σ| by (1 + 2.08e-7) shifts ε by 1.28e-4, about 500 times the
  shift from the 1e-9 phase seed. The supplied preparation has to be specified radially to about 1e-7 to give ε to
  1e-4. This does not affect the ensemble mean that Codex reports.
- The energy ledger closes to 1.2e-8 and Gauss's law to 2.7e-9 at rtol 1e-12.

### 3b. Codex's Y excitation and its extrapolation

- **Forced-solution subtraction: CONFIRMED.** Our independent physical-field version of the static forced solution
  and its derivative reproduces J_Y to 1e-8.
- **Neglected dynamic response is small** [D]. The static solution omits the response at the drive frequency 3ω_X.
  The relative correction is (3ω_X)²/m_Y² ≈ 2e3/a², so 0.5% in amplitude at a = 600 and 1e-5 in energy.
- **The excitation is not an artefact.** E_exc ≈ 950 (conformal units) exceeds E_X at a = 800. The number is
  stable across windows (J_Y = 13.29, 13.27, 13.26, 13.26 from 300–600 to 2000–3000).
- **Our number is 3.2% larger: CORRECTED.** The eigenmode count is N_Y = 13.69 at a = 1000–3000, not 13.26.
  - The two real Majorana eigenstates are split, m_H/m_L = 1.0348/0.9640 M_Y, ΔM = 7.08e10 GeV.
  - The single-mass convention undercounts the light state, which holds 94% of the quanta late on (N_YL = 12.84,
    N_YH = 0.85).
  - The heavy population is 2.92 at a = 100–300. It is depleted by Σ-rotation mixing: the Gauss-law Σ phase rotates
    at θ̇/2 ≈ 2.4e5/a³, which equals ΔM/2 near a ≈ 425.
- **Width: CONFIRMED, with a mass-eigenstate correction** [D].
  - Γ(Y → 3X) = 3κ²M/(256π³) = 377.95 GeV, with |M|² = (6κ)² and a 1/3! symmetry factor.
  - The state that actually decays is the light real eigenstate. It decays to XXX or X̄X̄X̄ with half the rate
    each, so Γ_L = 3κ²m_L/(256π³) = 364.3 GeV.
- **Decay energy: CORRECTED from 275 to 349 E_X** [D, C].
  - Codex sets H = Γ, giving a_dec = 1.722e5. The energy actually deposited by an exponential decay in radiation is
    weighted by ⟨a⟩ = Γ(3/2)√(2H_i/Γ) = 2.199e5.
  - With N_Y = 13.69 and m_L, the deposited energy is 349 E_X, not 275 E_X (×1.27).
  - The energy-equivalent fraction for Y alone is then 0.0030, not 0.00356.

### 3c. Item (a): decay products and the symmetric population [D, C, E]

| Channel | Width (GeV) | Γ = H at | Products |
|---|---:|---:|---|
| Y_H → Y_L + A_L (Goldstone equivalence: coupling √2μ G Y₁Y₂) | 1.27e9 | a ≈ 94 | one vector of 6.8e10 GeV per heavy quantum; Y number unchanged (**new channel**) |
| σ (Σ radial, m_σ = 1.414e11 = 2ΔM) → A_L A_L | 2.81e9 | a ≈ 63 | two vectors of 7.07e10 GeV per σ quantum (**new source**) |
| Y_L → XXX or X̄X̄X̄ | 364.3 | ⟨a⟩ = 2.20e5 (T_vis ≈ 1.3e10 GeV) | three X with dN/dE = 8E/m_L² on [0, m_L/2], mean m_L/3 = 3.2e11 GeV |
| Y → XX through one condensate leg | Γ/H = 3.5e-3 (constant in radiation) | — | 2.7% of the Y decays before the three-body time; neglected |
| A → X X̄ (m_A = 6eh = 6e4 GeV) | 3.98e-12 | ⟨a⟩ = 2.10e12 (T_vis ≈ 1.3 TeV) | one non-relativistic pair of 3e4 GeV quanta per vector |

- **Charge symmetry: CONFIRMED (Codex assumed it).**
  - Y_L is a real (CP-even) eigenmode, and the Y–Ȳ oscillation frequency ΔM is 2e8 times Γ_L.
  - There is no physical CP phase at dimension ≤ 4 with Z2_φ (gp_phase [Q]).
  - The residual Y-condensate charge is 3N_Y/N_X = −8.2e-4 at a = 2000–3000 and still falling: N_Y = 9.77, 6.59,
    2.04, 0.36, 0.065 across the windows. N_eff → N_X as a result.
  - The decays therefore add no net X charge at the 1e-3 level, and the final net charge is N_X.
- **The Σ excitation is real.** N_σ = 4.895, 4.889, 4.893 in the windows 600–1000, 1000–2000 and 2000–3000.
  - It is a 10% radial oscillation at a = 1000, consistent with Codex's endpoint |B|/A = 0.0638 against the vacuum
    value 0.0707.
  - The centrifugal and μY² forced shifts contribute ≤ 1e-3 quanta by a = 1000.
  - Codex does not count it.
- **The vectors are sterile until they decay.** The coupling is e = 1e-7, and the potential has no X–Σ portal, so the
  Goldstone reaches X only through heavy Y exchange.
- **Pairs per net charge, if never thermalised** (N_X = −238.20, late):
  - Y daughters: 1.5 N_Y/|N_X| = **0.0862** (Codex 0.0834);
  - σ vectors: 2N_σ/|N_X| = **0.0411**;
  - Y_H vectors: 0.0036 (late heavy population) to 0.0123 (a = 100–300 population);
  - **total 0.131–0.140**, which is **26.5–28.3 times** the 1% allowance of 0.0049296. Codex has 16.9 times.
- **Mass fraction: CORRECTED from 0.169 to 0.27–0.28.** This is the mass of these pairs relative to calibrated dark
  matter once the net yield equals Y_Q.
- **The frozen X remainder is larger again.** At ε = 0.2402 it alone gives f = 1.0143(1/ε − 1) = 3.21 (gp_leftovers
  number-count law). **The frozen total is f ≈ 3.5, which is 350 times the 1% allowance.**
- **Deposited energies, all disjoint** [C, E]:
  - Y daughters: 349 E_X;
  - vectors: 15.2 E_X per vector quantum (they are non-relativistic for about four decades before decaying), so
    161–193 E_X in total;
  - **total injection 510–542 E_X.** The energy-equivalent fraction is 0.0021–0.0022.
- **How double counting is avoided.**
  - E_X, E_Y,exc and E_σ are separate terms of the Hamiltonian.
  - The vector energy is taken out of the Y (ΔM per heavy quantum) and σ budgets.
  - The X onset budget ρ_c + V0 remains the single budget for remainder, binding pairs and quench debris (gp_critic).

### 3d. Item (b): thermalisation and number-changing processes [E, rates.py]

All rates are evaluated with H = H_i/a² (instant visible reheating). If the hidden sector thermalises,
T_h/T_vis = 0.023 before Y decay and 0.10–0.11 after it.

- **Prompt decays.** Y_H → Y_L A and σ → AA have Γ/H > 1 from a ≈ 60–100 onward.
- **Condensate fragmentation.** The quartic self-resonance has a Floquet exponent ≈ 0.036 per unit of the Lamé
  variable x = √(2λ_X)|X|a, which is 1.31 per unit a (Greene–Kofman–Linde–Starobinsky [R]). That gives growth of
  e^30 by a ≈ 600, so the homogeneous condensate is expected to fragment inside the release window. Lattice
  calculations were not run. N_X and E_X are conserved, so the ledgers above are unaffected. The late homogeneous ε
  drift of 5e-4 is not certified past fragmentation.
- **Energy sharing of the hard daughters.** E = 3.2e11 GeV against soft quanta of 5e7 GeV, with
  σ = λ_X²/(π s): Γ/H = 1.1e-2 at ⟨a_dec⟩, growing ∝ a. It reaches 1 at a ≈ 3e7 (T_vis ≈ 9e7 GeV). The cascade
  then finishes quickly, since Γ ∝ 1/E. Using the X–φ coupling g_X² = 0.01 once φ is populated makes this 100
  times faster.
- **Number changing.** Thermal 2→4 (∝ λ_X³, or g_X⁶ through φ) reaches Γ = H at T_vis ≈ 0.9 TeV (λ_X only) or about 900 TeV
  (with g_X²). Both are far above the light-scale onset (T_vis ≲ 54 GeV here), and condensate-assisted 2↔3 is
  faster still.
- **Vector decay products** (3e4 GeV quanta at T_vis ≈ 1.3 TeV) share their energy at Γ/H ≈ 1e5.
- **Channels at the light scale** [Q, gp_leftovers]:
  - 3→2 is forbidden (2m_X < m_s, and D is Z₂-odd), so 4→2 is the cannibal channel;
  - X X̄ → φφ and s → X X̄ only exchange energy within the hidden sector;
  - X X̄ → A' is closed (m_A/2m_X = 3.2e5);
  - the portal is 20 orders of magnitude too slow.
- **Conclusion.** **No channel removes hidden energy.** Kinetic and chemical equilibrium are expected well before
  the onset, so the maximal-equilibrium benchmark is the relevant regime, not the frozen one.

**Relic in that benchmark** (our gp_leftovers code with E0 = ρ_c(1 + R) + V0; `thermal_relic.py`) [C on the benchmark]:

| Hidden energy added, R = E_inj/E_X | σ (entropy per net charge) | T_h0/m_X | p_f (X̄ per net charge) | f (symmetric mass / calibrated DM) | f / 1% |
|---|---:|---:|---:|---:|---:|
| 0 (X only, control) | 0.99 | 3.5 | 0.0024 | 0.0049 (case L: erased) | 0.49 (passes) |
| 274.9 (Codex's Y energy) | 48.6 | 13.0 | 4.58 | 9.28 | 928 |
| 349 (corrected Y) | 58.1 | 13.8 | 5.54 | 11.2 | 1124 |
| 510–542 (Y + vectors) | 77–81 | 15.1–15.4 | 7.5–7.8 | **15.2–15.9** | **1520–1590** |
| 510, with the hidden charge 66 times denser during freeze-out (dilution after freeze-out) | 77 | 15.1 | 6.75 | 13.7 | 1370 |

- Case L (lumps plus conversion) equals case H in every injected row: with σ ≈ 50–80 no condensed phase ever forms,
  so lumps never form.
- The X-only control reproduces the review's 0.19 benchmark: f = 0.49% at ε = 0.240.
- **Thermalising the Y debris makes the leftover problem about 4 times worse than the frozen limit, not better.**
  Thermalisation turns the 510 E_X of injected energy into about 15 extra X X̄ per net charge.

### 3e. Item (c): can the 66-fold yield be diluted while keeping the window? [C, A, E]

- **Required dilution** [D]. If N_X were unchanged, a matter-dominated inflaton era dilutes by T_R/T_RI.
  T_RI(H_i) = 2.82e15 GeV, so 1/66 needs **T_R = 4.28e13 GeV** (Γ_φ = 2.6e9 GeV). The matter era then ends at
  H ≈ 2.3e-4 H_i, which is inside the charge-generation epoch (the RD-equivalent a ≈ 66). The release has to be
  recomputed.
- **Inflaton matter era at the leading H_i, Hubble mass removed at a = 1** (`md_release.py`; 11 phases):
  - at δ = π/2, N_X falls from −238.6 to −28.4 and ε_X from −0.241 to −0.0253, so Y_X = −0.11 Y_Q;
  - the best phase gives |ε| = 0.092 at Y = 0.41 Y_Q.
  - T_R dependence at δ = π/2 (`md_leading_TR.py`):

    | T_R (GeV) | T_R/T_RI | ε | Y/Y_Q |
    |---:|---:|---:|---:|
    | 1e13 | 0.0035 | −0.025 | −0.026 |
    | 4.28e13 | 0.0152 | −0.025 | −0.11 |
    | 1e14 | 0.0354 | −0.020 | −0.21 |
    | 3e14 | 0.106 | +0.025 | +0.77 |
    | 1e15 | 0.354 | +0.201 | +20.2 |

    A short matter era (dilution ×3.3) can keep |ε| ≈ 0.2, though with the sign flipped, which shows how sensitive
    the window is to the background. It still leaves Y = 20 Y_Q.
- **Hubble mass persisting until reheating** (bg 2): |ε| ≤ 0.064 and |Y| ≤ 0.35 Y_Q over 11 phases.
  - Argued [A]: with a persistent −c_H H²|X|², the scale-invariant quartic tracking solution only releases when the
    inflaton decays. The release then resembles Codex's instant release at H ≈ Γ_φ, and the yield relation is
    unchanged rather than diluted.
- **Coarse matter-era scan:** 132 runs. These cover H_i = 3e12–4.7e13 GeV with 11 phases each at T_R = 4.28e13 and
  1.5e14 GeV, and include the 22 leading-point runs.
  - **None** has |ε| ≥ 0.19 together with |Y| ≥ Y_Q.
  - The largest |ε| is 0.135 (H_i = 6e12, T_R = 1.5e14, δ = 3π/4, Y = −1.14 Y_Q). The median |ε| is 0.063, and
    14% of runs have |ε| > 0.1.
  - The yield reaches Y_Q only at |ε| ≈ 0.1 (up to 7.1 Y_Q at H_i = 4.7e13).
  - Controls: Gauss ≤ 1.2e-4 and ledger ≤ 1e-5 at rtol 1e-10.
  - This is not an exclusion: the RD window itself is only ±2% in H and ±0.06 rad wide.
- **Late visible-only entropy injection** (a new field, not in the model) [E]:
  - It leaves the release, ε and the pairs and entropy per net charge unchanged, so it can fix the yield by
    construction. It needs T_dom/T_d ≈ 66 with T_d above BBN, about 5 MeV.
  - Its hidden branching must stay below about r₀/D^{4/3} ≈ 2.7e-6/267 ≈ **1e-8**, or the hidden energy doubles.
    Here r₀ = ρ_h/ρ_vis ≈ 2.7e-6 after all deposits.
  - It does not help the leftovers: f = 15.2 if it acts before hidden freeze-out, 13.7 if after (§3d).
- **Longitudinal vector from inflation** (OPEN) [E].
  - Graham–Mardon–Rajendran (GMR) with the late mass m_A = 6e4 GeV scales gp_critic's ρ/s = 0.307 GeV by H_i²
    to 0.0175 GeV. That is **0.91 pairs per net charge** (185 times the allowance).
  - The GMR premise of a constant mass fails here: F_g(inflation) = 3.9e14 GeV against 6e11 GeV late. A frozen
    gauge angle (A_L = kθ/e is conserved superhorizon) rescales this by (F_late/F_inf)² = 2.4e-6, giving 2e-6 pairs.
  - The light physical mode (ψ₂, 0.52 rad rms) also feeds the late Σ phase, which is eaten. That transfer has not
    been computed.
  - The vector source therefore spans about 1e-6 to 1 pair per net charge. It is not needed for the verdict.

### 3f. Side check: the inflationary spread of the source phase [C, `phase_spread.py`]

- **Spread.** The light angular mode at this point has m² = 0.3716 H_i², which reproduces Codex's angular masses. Its
  stochastic equilibrium, ⟨φ²⟩ = 3H⁴/(8π²m²) [R: Starobinsky–Yokoyama], gives the source phase
  ψ₂ = 2θ_Y − θ_Σ an rms of **0.52 rad**. The heavy combination ψ₁ has 0.003 rad.
- **Ensemble.** Gaussian noise on the Σ phase about π/2, N = 128, seed 2026102501:
  - σ = 0.03 gives ε = −0.2462 ± 0.0034, consistent with Codex's −0.2491 ± 0.0033;
  - σ = 0.52 gives a net of **−0.0527**, mean −0.052 ± 0.010, with 64% of members negative and a 5–95% range of
    −0.25 to +0.12.
- **Reading.** Without a phase selector, the natural inflationary spread at this H_i reduces the net charge
  fraction to about 0.05. That is below 0.19 even for X alone (case H gives f ≈ 0.44 at ε = 0.05). The scope is
  phase noise only: radial fluctuations are not added, and this is not a computed spatial spectrum.

### 3g. Item (d): combined leftover verdict at this point

| Leftover limit | Symmetric relic f (relative to calibrated DM, after the yield is set to Y_Q) | Against 1% |
|---|---:|---:|
| X only, maximal equilibrium (Codex's "clears 0.19" test) | 0.0049 | passes (×0.49) |
| Frozen: X remainder + Y daughters + vectors, no number changing | 3.21 + 0.27–0.28 ≈ **3.5** | ×350 |
| Maximal equilibrium including the debris energy (the regime the rates point to) | **15.2–15.9** | ×1520–1590 |
| Same, with dilution only after hidden freeze-out | 13.7 | ×1370 |
| Lumps plus conversion (case L) with the debris energy | same as case H (no condensed phase forms) | ×1520 |
| Intermediate: kinetic but not chemical equilibrium (4.4 quanta per charge, T ≈ 66 m_X, μ/T ≈ −5.8), then 4→2 | σ ≥ 43, so f ≳ 8 [E] | ≳ ×800 |
| Natural phase spread (net ε ≈ 0.05), any of the above | worse | — |

**Verdict: leftovers FAIL at this point** on homogeneous calculations plus estimates. The excess is at least about
350 times the allowance in every limit computed. This is independent of how the 66-fold yield is diluted, because
visible-only dilution leaves the pairs and the entropy per net charge unchanged. The X-only pass does not survive
the release's own Y, Σ and vector debris.

- **Not included** (each would only add to the symmetric population):
  - longitudinal vectors from inflation (§3e);
  - non-perturbative vector production in the changing Higgs background during the release;
  - string and wall debris;
  - the 2.7% condensate-assisted Y → XX channel.

## 4. Verdict table

| # | Codex claim (373d37e) | Codex value | This track | Verdict |
|---|---|---|---|---|
| 1 | Leading point: ε_X, N_X, E_X, J_Y, Y_X | −0.24066124, −238.64868, 741.64988, 13.26903, −65.98 Y_Q | −0.24066124, −238.64867, 741.64988, 13.26903, −65.98 Y_Q (independent formulation and integrator) [C] | CONFIRMED |
| 2 | 1e-9 phase shift → Δε ≈ 2.6e-7; late drift 0.00046 | as stated | 2.57e-7; 0.00046 [C] | CONFIRMED (new: radial sensitivity, 2e-7 in \|Σ\| → 1.3e-4) |
| 3 | Late forced-solution subtraction gives the Y excitation | J_Y = 13.269 | same; the dynamic correction is 1e-5 in energy [C, D] | CONFIRMED |
| 4 | Y number | 13.27 (single mass) | 13.69 (Majorana eigenmodes; 94% light) [C] | CORRECTED (+3.2%) |
| 5 | Γ_Y = 3κ²M_Y/(256π³) | 377.95 GeV | 377.95; the decaying light state has 364.3 [D] | CONFIRMED (eigenstate value CORRECTED) |
| 6 | Y energy at decay | 275 E_X (at H = Γ) | 349 E_X (exponential-decay average, N_Y, m_L) [D, C] | CORRECTED |
| 7 | Energy-equivalent fraction (diagnostic only) | 0.00356 | 0.0030 (Y only); 0.0021–0.0022 (Y + vectors) [C] | CORRECTED |
| 8 | Extra pairs per net charge, unthermalised | 0.0834 (16.9× the allowance) | 0.131–0.140 (26.5–28.3×), adding the σ → AA and Y_H → Y_L A vectors (each A → X X̄) [C, E] | CORRECTED |
| 9 | Unthermalised Y mass fraction | 0.169 | 0.27–0.28; frozen total with the X remainder f ≈ 3.5 [C] | CORRECTED |
| 10 | Zero mean Y charge (assumed) | assumed | Y_L is real; no CP phase; 3N_Y/N_X = −8e-4 and falling [C, D] | CONFIRMED |
| 11 | "Do not apply 0.19 after the injection; not a solved relic" | caveat | Solved in the maximal-equilibrium benchmark: f = 15.2–15.9 (σ = 77–81) [C] | CONFIRMED caveat; relic now computed: FAIL |
| 12 | Raw 66-fold yield, no dilution supplied | open | An inflaton matter era with T_R = 4.28e13 GeV (the ×66 dilution) cuts ε to 0.025 at the point; over the 132 scan runs \|ε\| ≤ 0.135, and none has \|ε\| ≥ 0.19 with Y ≥ Y_Q. A short era (T_R = 1e15) keeps \|ε\| = 0.20 but leaves Y = 20 Y_Q. Late visible-only injection needs a new field with hidden branching ≲ 1e-8 and leaves the leftovers unchanged [C, E] | NOT RESOLVED; the dilution routes checked fail or need new physics |
| 13 | A blanket "no scale works" claim is premature | stated | Agreed for the X-only homogeneous test. At this point the leftover requirement fails regardless; with the natural ψ₂ spread (0.52 rad), net ε ≈ 0.05 [C] | CONFIRMED (scope); the point itself fails on leftovers |

## 5. Limits

- Homogeneous k = 0 releases; no lattice. Fragmentation is expected near a ≈ 600 [E].
- Decay widths are tree level with vacuum masses. The Y_H → Y_L A count during the Σ-rotation era (a < 425) is
  bracketed, not solved (0.85–2.9 quanta).
- σ → AA uses the late classical N_σ. Broad-resonance Goldstone production would give more and softer vectors, hence
  more pairs.
- The thermal relic is our gp_leftovers maximal-equilibrium benchmark: X and X̄ degrees of freedom only, 4→2 at
  threshold without Bose enhancement, and coarse g*. s and D add ≤ 19% to σ [Q].
- Thermalisation rates are order-of-magnitude estimates.
- The inflaton-matter scan is coarse: 5 further H values × 11 phases × 2 T_R, plus 22 runs at the leading point (132
  in all), and 5 T_R values at δ = π/2. It is not an exclusion.
- The vector from inflation is OPEN.
- Codex's code was not run. The comparisons in §3a use its published JSON.

## 6. Credit

- Affleck–Dine; Dine–Randall–Thomas.
- Q-ball formation: Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki [R].
- λφ⁴ self-resonance: Greene, Kofman, Linde & Starobinsky, PRD 56 (1997) 6175 [R]. Turbulent thermalisation:
  Micha & Tkachev [R].
- Goldstone equivalence (Cornwall–Levin–Tiktopoulos; Lee–Quigg–Thacker) [R].
- Cannibal dark matter: Carlson–Machacek–Hall [R].
- Graham–Mardon–Rajendran (vector from inflation) [R].
- Starobinsky–Yokoyama (stochastic spectator equilibrium) [R].
- Non-instant reheating and entropy dilution: Scherrer–Turner; Giudice–Kolb–Riotto [R].
- Mixed or decaying dark matter structure constraints are not used here.

## 7. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/sp_yledger/code
python3 point_ledger.py ../results/point_ledger.json                              # ~2 min (5 releases to a=3000)
python3 rates.py ../results/point_ledger.json ../results/rates.json               # 1 s
python3 thermal_relic.py ../results/point_ledger.json ../results/thermal_relic.json  # ~10 min (uses ../../gp_leftovers/code)
python3 md_release.py ../results/md_release.json                                  # ~20 min, 4 cores
python3 md_leading_TR.py ../results/md_leading_TR.json                            # ~2 min
python3 phase_spread.py ../results/phase_spread.json                              # ~4 min
```
yrelease.c is compiled automatically (gcc -O3) into the scratchpad build directory.
