# SKEPTIC of track mf_hybrid (key mf_hybrid_skeptic)

2 November 2026 · Claude, adversarial pass on `tracks/mf_hybrid` (BRIEF_MULTIFIELD.md; GitHub issue #13 context).
- This is same-session AI work on one supplied model, not outside replication.
- Nothing here is a validation, an exclusion of the real world, a new theory or a GRUT prediction.
- Every potential, coupling, renormalisation condition, slice, friction model and criterion, in the track and here, is a
  model choice. None follows from closure.
- Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6, DESI DR2) are retrospective.
- Git was read-only and nothing was committed. Every run was in the foreground with one process. No compiled files and
  no `__pycache__` are in the tree (python `-B`). No lattice was run.

Labels: [D] derived, [C] computed with my own code (`skeptic/code/`, output in `skeptic/results/`), [E] estimate,
[A] argued, [Q] quoted from a reviewed track, [R] literature (arXiv numbers checked by web search in this session).

Units: H_ref = 1.30×10¹¹ GeV (the track's design H; H_i/H_ref = 1.0000187), M_P = 2.435×10¹⁸ GeV (reduced),
t in H_i⁻¹, a = 1 at σ = σ_c (the track's slice).

**Method policy.** Every number below was recomputed with a different method or code path from the track's:
- my own re-implementation of the potential, taken from the exported parameters;
- a different Mukhanov–Sasaki formulation;
- an exact 2×2 exponential propagator in place of RK4;
- a coupled (δσ, δχ) mode system in place of δχ alone;
- a per-species fluid.

The track's production tables `axial_tables.json` are used as input only after I reproduced them: at 8 points to
10⁻⁶, and at 5 more to the 3 digits printed.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

| Track result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| B1, B2: CMB of HFa and HFb | **stand** | Independent MS (v-equation, finite-difference V′, V″, exact tensor modes). HFa: n_s = 0.964901, A_s = 2.0999×10⁻⁹, α_s = +1.069×10⁻⁴, r = 2.7506×10⁻⁷. HFb: 0.974299, 2.0998×10⁻⁹, +1.434×10⁻⁴, 2.7507×10⁻⁷. N_* = 52.6445 from my own K bookkeeping. Pulls 0.0/−2.76σ and +2.24/0.0σ. n_s barely depends on K: K = 3.3 → 8 moves it by only +2.5×10⁻⁵ to +7.4×10⁻⁵ | C |
| A1: exact half-period axial production | **stands as a single half-period** | My exponential propagator reproduces ε(A, M) to ≤ 10⁻⁶ at all 8 check points, and the tables at 5 more to 3 digits | C |
| A1 used as a sustained rate Γ = (m/π)ε/(S²/2) | **refuted** | Both helicities evolved exactly through 16 half-periods: the band **saturates (Pauli) after 3–6 half-periods**. F_sat = 0.13 (χ→N₂), 0.085 (N₃), 0.055 (N₄) of ρ_χ and 0.14 (σ→N₁) of ρ_σ. The track's rate keeps adding 0.091 / 0.038 / 0.015 / 0.050 every half-period (§2.1) | C |
| B4: radiation 54/77/85% at t = 0.3/0.6/3; K = 2.05 | **refuted** | A per-species fluid with that saturation, unblocked only by N decays: f_rad = 0.21 / 0.20 / 0.15 at t = 0.3 / 0.6 / 3. **K = 7.1** if the perturbative tails are free; 6.4 / 4.1 with an extra N-thermalisation rate of 1 / 10 H_i; **≈ 22** if the tails are blocked too. My unblocked control reproduces the track (K = 2.03) | C on E |
| B5: SM plasma ≥ 10⁻³ at τ = 0.64, ≥ 50% at τ = 7.3; T_eq ≈ 2.2×10¹⁴ GeV | **refuted in part** | ρ_SM/ρ ≥ 10⁻³ at τ = 0.31–0.33 (0.60 if m̃ = 0.015 eV) still holds. But ρ_SM ≥ 50% arrives only at τ ≈ 160–240, and there is no early radiation–matter equality | C on E |
| B3: fast waterfall, ε_H = 1 at 0.19, H = 0.997 H_i | **stands** | t = 0.1926, H = 0.99665, R_min = −4.54 | C |
| B3: "no walls; classical/quantum ≥ 33; ln P ≈ −550" | **refuted as stated** | With coupled modes and the track's own friction the ratio is 20.5 (still safe). With the Pauli-limited friction, δχ_rms reaches 0.3 v at t = 0.43, and ratio(t) falls below 2 at t = 0.46, while χ's inner turning points sit at 0.37–0.59 v. The waterfall fragments, and walls are not excluded without a lattice. If walls percolate they dominate (t ≈ 23) before the bias removes them (t ≈ 150) (§2.2) | C on E |
| S1: one-loop flatness forces m_σ ≲ 13 H_i (X ≲ 0.5–4) for a heavy waterfall | **refuted (the number)**; the valley argument stands | The track's own row (4000, 30), with X = 3.7, passes at all three κ (N_max = 108/85/54). My scan (21 κ values, V0-dominated, no near-inflection) at m_χ = 4000 H_i passes for m_σ = 13 (21/21 κ), 30 (20/21), 60 (11/21), 90 (4/21) and 150 (1/21), i.e. X up to 18. The bound is a κ (tree-quartic) tuning statement, not m_σ ≲ 13. "σ holds V0 when m_χ ≫ m_σ" stands (§2.3) | C |
| S2: perturbative bound m ≥ (384π)^{1/5}(M_P/H)^{2/5} H = 3.4×10³ H_i | **stands** | Re-derived. Dirac gives 2.9×10³, N_f = 10 gives 2.1×10³, a trilinear to a complex scalar gives a similar value | D |
| B6: spectator | **refuted in part** | s_* and β_iso = 3.82×10⁻³ p² are reproduced; s(t) was reproduced by mf_ideal to 5×10⁻³ [Q]. **The energy fraction at s decay is 1.3% (Γ_s t = 1) or 0.94% (Γ_s = H), not 2.4%**: the formula drops a_on³/K². "No hidden heat" holds only if the induced holder portal A_sΣ s\|Σ\|² vanishes (§2.4). The Planck bound is β_iso < 0.038, so p ≤ 3.2; mf_ideal's p ≈ 4.4–5.1 would give β = 0.075–0.098 | C/D |
| B7: gravitational X ≈ 5×10⁻⁴ of E_X; no line near 22.25/44.5 | **stands (×2)** | The rate is A²ω/64π per real d.o.f. with ω = 2m_χ, not m_χ, so 1.1×10⁻³ of E_X. Still negligible. No line near the Y resonance | D |
| O1 (open) | — | not addressed | — |

---------------------------------------------------------------------------------------------------------------------

## 1. What I recomputed, and how the method differs

| Check | Track's method | Mine | File |
|---|---|---|---|
| A_s, n_s, α_s, r | ζ-equation, analytic V′, V″, r from 2H²/π²M_P² | v = zζ equation with z″/z from the exact background (V″ by Richardson differences), 7 k values, cubic fit, exact tensor modes. N_* from my own bookkeeping (g_s0 = 3.909 and 3.931) | `sk_cmb.py` |
| e-fold capacity | 12 pivots, 3 κ | my own loop derivatives, 100 pivots ≤ 0.2 M_P, 21 κ in 0.1–10³, a V0-dominance cut (\|V_loop\|/V0 < 0.05), a regularity cut (V′ ≥ 0.1 V′_*, i.e. no hidden near-inflection) | `sk_structure.py` |
| axial production | RK4, one half-period, one helicity | exact 2×2 exponential (unitary); both helicities through 16 half-periods at fixed amplitude; occupations at every turning point | `sk_axial.py` |
| reheating fluid | one Γ per condensate from the tables at every half-period; mean N width | per-species N energy and number; in-band blocking energy with redshift and band-edge shrinkage; Pauli saturation from my mode runs; perturbative tails free or blocked; N-thermalisation bracket Γ_th = 0, 1, 10 H_i | `sk_pauli.py` |
| waterfall and defects | δχ_k only, mass M_χ²; vacuum spread added in full; ratio sampled at 7 times | coupled (δσ, δχ) with the full tree Hessian (V_σχ = 2g²σχ ~ m_χ²); two vacuum solutions; amplified spread only; ratio minimised over the whole series, including χ's inner turning points; three friction models | `sk_waterfall.py` |
| bench JSON | — | an audit of every column against its definition (N, a, ε_H, R, K, fractions, s, f_s) | `sk_bench_check.py` |
| spectator, gravitational, ν masses, portal leak | estimates | re-derived (§2.4–2.6) | `sk_misc.py` |

---------------------------------------------------------------------------------------------------------------------

## 2. Corrections with numbers

### 2.1 Pauli blocking of the axial channels [C, then C on E] (`sk_axial.json`, `sk_pauli.json`)

**What the exact modes show.**
- In the rotated frame each (k, h) mode is a two-level system. n_k ≤ 1 holds automatically, and repeated crossings do
  not keep producing: the occupation oscillates (Greene–Kofman hep-ph/9807339, hep-ph/0003018 [R]).
- At the slice amplitude (A = 42.3 for χ with f/m = 8; A = 12.3 for σ with f/m = 2), the cumulative energy fraction is:

| channel | track's F per half-period | track's fluid after 16 half-periods | exact, after 2 / 4 / 16 half-periods | saturated F_sat | band n̄ |
|---|---|---|---|---|---|
| χ → N₂ (M = 400) | 0.091 | 1.46 | 0.18 / 0.13 / 0.13 | **0.13** | 0.27 |
| χ → N₃ (M = 200) | 0.038 | 0.60 | 0.076 / 0.11 / 0.085 | **0.085** | 0.16 |
| χ → N₄ (M = 120) | 0.015 | 0.24 | 0.030 / 0.054 / 0.051 | **0.055** | 0.11 |
| σ → N₁ (M = 250) | 0.050 | 0.80 | 0.10 / 0.16 / 0.145 | **0.14** | 0.17 |

- So the three χ species together can take at most ≈ 27% of ρ_χ, and σ's species ≈ 15% of ρ_σ.
  - F_sat scales roughly as A^{1.7–2}: it is computed at A = 10, 20, 42 (3, 6, 12 for σ) in `sk_axial_Fsat_vs_A.json`.
  - This is the "≈ 30–40% cap" the track names in its caveats, now computed. It is incompatible with the 54% and 77%
    radiation the same track exports at t = 0.3 and 0.6.
- The N's leave the band only by decaying.
  - Their lab-frame rates are Γ_rest M/E ≈ 0.07 (N₂), 0.009 (N₃), 0.002 (N₄) and 0.03 (N₁) H_i at E ≈ ¾μ5.
  - Redshift does not unblock: the band edge μ5 ∝ ρ_c^{1/2} shrinks faster than the N momenta.
  - N–N scattering through χ exchange (y_p = M/f ≤ 0.2) is ≈ 10⁻² H_i [E].

**The fluid** (bench state at t = 0.19, HFa; HFb differs by ≤ 0.1%):

| case | f_rad at t = 0.3 / 0.6 / 1 / 3 | K | ρ_SM ≥ 10⁻³ | ρ_SM ≥ 50% |
|---|---|---|---|---|
| unblocked (my control) | 0.56 / 0.78 / 0.83 / 0.86 | 2.03 | τ = 0.27 | 7.7 |
| track (bench_HFa) | 0.54 / 0.77 / 0.83 / 0.86 | 2.05 | 0.64 | 7.3 |
| Pauli, tails free, Γ_th = 0 | 0.21 / 0.20 / 0.18 / 0.15 | **7.05** | 0.33 | 214 |
| Pauli, tails free, Γ_th = 1 H_i | 0.22 / 0.23 / 0.23 / 0.21 | **6.41** | 0.33 | 210 |
| Pauli, tails free, Γ_th = 10 H_i | 0.28 / 0.42 / 0.48 / 0.51 | **4.10** | 0.31 | 157 |
| Pauli, tails blocked too | 0.21 / 0.19 / 0.17 / 0.12 | **≈ 22** (not converged by 3×10⁵) | 0.33 | not reached |
| Pauli, tails free, m̃ = 0.015 eV | 0.21 / 0.20 / 0.18 / 0.14 | 7.13 | 0.60 | 240 |

- "Tails free" treats the perturbative products (σ → N₁ at 0.29 H_i, χ → N₄ at 5×10⁻³ H_i) as unblocked. A
  thermalised N bath would allow this; the products sit inside the filled band, so it is the optimistic case.
- K therefore lies in **4–7 for any N-thermalisation rate up to 10 H_i, and near 20 without it**, not 2.05 [C on E].
  Since R_vec ∝ K², the vector heat rises by a factor of 4–12.
- The track's own bracket (K = 3.26 if the non-perturbative transfer stops at t = 0.6) is too optimistic, because it
  starts from the unblocked 77% radiation.
- N_* moves by −0.35 to −0.6 e-folds. That leaves n_s unchanged to 10⁻⁴ and needs A_s retuned by 1–2% (§0).

### 2.2 Waterfall, defects and fragmentation [C on E] (`sk_waterfall_HFa.json`, `sk_waterfall_HFb.json`)

| friction model | f_rad(0.3) | min ratio \|χ_cl\|/δχ_rms | δχ_rms = 0.3 v at | ratio < 2 at | first inner turning points (χ/v) |
|---|---|---|---|---|---|
| track (unblocked) | 0.54 | 20.5 (t = 0.19) | never (≤ 0.037 v) | never | 0.44, 0.56, 0.63 |
| Pauli-limited (§2.1) | 0.16 | 0.15 (linear, t = 0.6) | t = 0.43 | t = 0.46 | 0.37, 0.43, 0.47, 0.50 |
| Hubble only | 0 | 0.004 | t = 0.34 | t = 0.34 | 0.25, 0.34, 0.40 |

- **The coupled modes grow faster** than the track's δχ-only modes, even with the track's friction: the k = 0 growth is
  ×7.7×10³ against ×1.35×10³, and δχ_rms/v at t = 0.6 is 0.019 against 0.0053.
- **The track's "ratio ≥ 33" comes from its earliest sample (t = 0.05).** At that time its δχ_rms is dominated by the
  unamplified vacuum. That choice is conservative, but the sampling misses the minimum.
- **Where the danger lies.**
  - With only the friction that the Pauli cap allows, χ keeps ~75% of V0.
  - Its inner turning points stay inside the spinodal region |χ| < v/√3, so tachyonic self-resonance continues
    (Felder et al. hep-ph/0012142, hep-th/0106179 [R]).
  - The linear spread reaches |χ_cl| at the turning points by t ≈ 0.45.
  - The Gaussian wrong-side probability there is O(10⁻²–0.4) per correlation volume. The bias (ΔV/V0 = 7×10⁻⁵) is far
    too weak to stop a non-linear field from populating −v.
- **Walls, if they form and percolate** [E, scaling estimates]:
  - tension σ_w = (2√2/3)√λ v³ = 2.2×10¹³ H_i³;
  - they would dominate at H ≈ σ_w/3M_P² = 0.021 H_i (t ≈ 23);
  - the bias would only remove them at H ≈ 2cv/σ_w = 3.3×10⁻³ H_i (t ≈ 150).
- **Consequence.** "No walls" is not established for the HF sector. The waterfall's outcome needs the 64³ lattice the
  brief allows, run with the fermion back-reaction (or a Pauli-limited friction).
- **HFa's spectator.** Its mass, m_s = g_s|χ|, would then be inhomogeneous at O(1). That breaks the homogeneous release
  the chain assumes. HFb's ξ_s R mechanism does not have this problem.

### 2.3 The structural bound S1 [C] (`sk_structure.json`)

- **Same definition, more pivots and κ values.** I used the track's definition of the capacity (tree c1, m2 solved at
  each pivot for A_s and n_s; N = ∫V/M_P²V′) on 100 pivots and 21 κ, keeping only V0-dominated solutions with no
  near-inflection.
- **Regression.** (250, 600) gives N_max = 55.8 / 74.3 at κ = 1 / 10 (track: 55.8 / 73.0).
- **The scan.** The largest m_σ with N_max ≥ 52.6 for some κ is:
  - ≈ 150 H_i at m_χ = 4000 (X = 18.5; 1 of 21 κ);
  - 90 H_i passes for 4 κ, and 60 H_i for 11;
  - at m_χ = 3000, 150 H_i passes for 1 κ.
- **So the "bound" is a statement about κ tuning, not m_σ ≲ 13.** The track's own (4000, 30) row already passes at all
  three of its κ.
- **Unchanged:**
  - σ, not χ, holds V0 whenever m_χ ≫ m_σ;
  - S2 then forbids fast perturbative reheating of σ.
- **What changes.** The heavy-waterfall class is not forced to be quiet (2–4 H_i). m_σ = 30–150 H_i clocks with
  amplitude σ_c = 2√3 M_P/m_σ = 0.02–0.12 M_P exist. They could be emptied by the same axial route at the EFT edge,
  subject to the same Pauli cap. The track leaves this open; I did not compute it.

### 2.4 The spectator s [C/D/E] (`sk_bench_check.json`, `sk_misc.json`)

- **Energy fraction at decay.**
  - From the exported history, f_s/a = 7.35×10⁻⁶ in radiation.
  - At Γ_s t = 1 (a = 1.80×10³) this gives **1.32%**, and at Γ_s = H (a = 1.27×10³) it gives **0.94%**.
  - The track's 2.43% uses 1.5(s_*/M_P)² a_dec, which drops the factor a_on³/K² = 0.545. The entropy injection into the
    SM scales down with it.
- **Isocurvature.**
  - β_iso = (H/2πs_*)² p²/A_s = 3.82×10⁻³ p² is reproduced.
  - Planck 2018 X gives β_iso < 0.038 for uncorrelated, n_II = 1 CDI ("axion I") at k_mid = 0.05 Mpc⁻¹
    (arXiv:1807.06211 [R]). So p ≤ 3.2.
  - mf_ideal's p = d ln Ω/d ln s_* on fast-turn-on histories is 4.4–5.1 [Q mf_ideal isocurvature.json]. At s_* =
    3×10⁻³ M_P that gives β = 0.075–0.098, a retrospective conflict.
  - The bench would need s_* ≳ 4.1–4.7×10⁻³ M_P. That raises f_s by a factor of 1.9–2.5, and mf_chain must compute p on
    the bench.
- **The decay channel is under-specified.**
  - y_s s N̄_sN_s with N_s → iγ5N_s leaves N_s massless in the vacuum. To reach the SM, N_s needs a soft Majorana mass
    (< m_s/2 = 1.5 H_i) and a seesaw Yukawa: a fifth singlet, which adds to Σm_ν (§2.6).
  - The induced s-odd terms are then small: a drift of ~10² H_i over 50 e-folds against s_* = 5.6×10⁴ H_i [E].
  - Without that mass and Yukawa, N_s is dark radiation, with ΔN_eff ≈ 0.04 at f_s = 1.3% [E].
- **The N_s loop sets the s mass during inflation** [D]: m_s² = −(y_s⁴s_*²/32π²)(12L − 4) = −0.0045(12L − 4) H_i², with
  L = ln(y_s²s²/μ²).
  - "≤ 0.01 H²" holds only for |L − 1/3| ≲ 0.2, which is a renormalisation choice.
  - At L = 1 the mass² is −0.036 H², so s drifts by ~1% per e-fold.
- **"No hidden heat"** omits the induced holder portal A_sΣ s|Σ|², generated by the Y loop through g_hμ*.
  - A threshold-size coefficient gives A = 2–8×10⁻⁵ H_i and BR(s → A_LA_L) = 3×10⁻¹¹–5×10⁻¹⁰ [E].
  - mf_ideal's ledger turns this into R_leak ≈ 0.36–5.8 at HFa's Γ_s, comparable to the whole Q-matter budget [Q].
  - It vanishes only under the "portals vanish at M_Y" condition.

### 2.5 Gravitational production and the Y resonance [D] (`sk_misc.json`)

- For V = ½ A cos(ωt) x² the energy rate is A²ω/64π per real d.o.f. (from Γ(φ → xx) = σ²/32πm). The R oscillation has
  ω = 2m_χ.
- For complex X with A = (5/6) × 3ρ_χ/M_P² = 6.4 H_i², the rate is 202 H_i⁵, against the track's 101: **1.1×10⁻³ of
  M1r's E_X**. Negligible.
- With Pauli blocking the χ remnant persists, but A ∝ H² keeps the integral dominated by t ≲ 1 [E].
- No background line lies near 22.25 or 44.5 H_i [C, reproduced].

### 2.6 CMB, one-loop potential, CP bookkeeping, kinematics: what I checked and found standing

- **CMB.** Reproduced to 3×10⁻⁵ in A_s, 10⁻⁶ in n_s and 10⁻³ in α_s (relative). r ≪ 0.036 (BICEP/Keck 2021,
  arXiv:2110.00483 [R]). The pulls use Planck 2018 0.9649 ± 0.0042 and ACT DR6 0.9743 ± 0.0034 (arXiv:2503.14452
  [R]; v2 abstract 0.974 ± 0.003).
- **One-loop potential** [D/C].
  - The χ loop normalisation (one real d.o.f., MS-bar) and its derivatives are re-derived. The loop's share of η at the
    pivot is 58% (HFa) and 74% (HFb) with my finite differences, as stated.
  - The other fields add no σ dependence: the N masses are constant, and s, Higgs and the hidden sector have no portal.
  - Labelling: the scale is μ² = 10 g²(2σ_c)² = 2×10 m_χ², not "10 M_χ²(2σ_c)". M_χ at the pivot is 622 (HFa) and 556 (HFb) H_i, not ≈ 800.
- **Phase-mismatch potential** [A]. It needs two L-violating mass insertions on different N's (three loops, ∝ y⁴M M′Λ²).
  The track's two-loop 400 H_i⁴ is therefore conservative.
- **CP bookkeeping** [D].
  - The holder is the only s-odd term. g_s²χ²s², ξ_sRs², m_s0²s² and y_s sN̄_sN_s (with N_s → iγ5N_s) are even.
  - The soft N_s mass of §2.4 breaks this only at the drift level.
  - The CP sign is α + sign(s_*), as stated.
- **Kinematics.**
  - Open: σ → N₁N₁ (500 < 600) and χ → N₄N₄ (240 < 250).
  - Closed: χ → σσ, and perturbative χ → N₂,₃.
  - Γ(χ → ss) = g_s⁴v²/8πm_χ = 9.6×10⁻¹⁴ H_i.
  - These are confirmed, but every one of the axial channels is Pauli-limited (§2.1).

### 2.7 Are the exported histories correct and usable by the chain? [C] (`sk_bench_check.json`)

**Correct:**
- N = ∫H dt to 3×10⁻⁴;
- a = e^N and K_running = a²H to 10⁻⁶;
- R = 6H²(2 − ε_H) in the homogeneous part to 10⁻⁵;
- R = 3H²(f_σ + f_χ) in the fluid part to 10⁻⁶ (except the switch row);
- the fractions sum to 1 to 10⁻⁷;
- the s envelope joins the exact s at t = 60 within 3%;
- the size is 0.57 MB;
- mf_ideal's reader reproduces s(t) to 5×10⁻³ [Q].

**Defects of the export:**
1. **f_s in 0.6 < t < 60 is only the potential energy** ½m_s²s_*²s², oscillating down to 5×10⁻¹¹. The full energy is
   4.3–12×10⁻⁵.
2. **R switches from the raw value to the cycle mean at t = 0.6.** There the 2m_χ amplitude R_osc,χ = 3ρ_χ/M_P² ≈ 0.69
   H_i² is three times the mean (0.23). mf_ideal found a 7% shift in S/Q between raw and cycle-mean R. With Pauli
   blocking the χ remnant (~75%) makes this larger.
3. **The whole post-0.19 history (fractions, K, R, m_s(t) for HFa) is the unblocked one of §2.1.** It should be treated
   as the K ≈ 2 corner, not as the HF model's history.

---------------------------------------------------------------------------------------------------------------------

## 3. New issues (ranked)

1. **N1. Pauli blocking caps the axial reheating** [C/E].
   - It invalidates K = 2.05, the radiation fractions and T_eq.
   - K ≈ 4–7 with the perturbative tails free, ≈ 22 with them blocked.
   - Since R_vec ∝ K², the vector heat rises by a factor of 4–12 (§2.1).
2. **N2. The waterfall's defect-freedom rests on the [E] unblocked friction** [C/E].
   - With the Pauli-limited friction the coupled linear modes reach |χ_cl| at χ's inner turning points by t ≈ 0.45.
   - Walls cannot be excluded without a lattice. If they percolate, they would dominate (t ≈ 23) before the bias
     removes them (t ≈ 150) (§2.2).
3. **N3. The coupled (δσ, δχ) modes grow ×3–6 faster than δχ alone** [C]. Even with the track's friction the minimum
   ratio is 20.5, not 33.
4. **N4. Isocurvature.** β_iso < 0.038 (Planck 2018) gives p ≤ 3.2. With mf_ideal's p = 4.4–5.1, s_* = 3×10⁻³ M_P gives
   β = 0.075–0.098 (retrospective). s_* ≳ 4.1–4.7×10⁻³ M_P is needed unless mf_chain finds a smaller p on the bench
   [C/Q].
5. **N5. Neutrino-mass sum** [D, retrospective].
   - Four (or five) seesaw singlets with m̃ = 0.05 eV each give Σm_ν = Σm̃ = 0.20 (0.25) eV for a real Casas–Ibarra R
     (hep-ph/0103065 [R]).
   - Planck 2018+BAO gives < 0.12 eV (arXiv:1807.06209 [R]); DESI DR2+CMB gives < 0.064 eV (arXiv:2503.14744 [R]).
   - A complex R with Σ_j|R_ji|² ≈ 3 is needed: cancellations, and large CP phases in N decays, whose baryon yield must
     then be checked. Alternatively m̃ must be smaller, which leaves K at 7.1 but delays ρ_SM ≥ 10⁻³ to τ ≈ 0.6.
6. **N6. Anomalous couplings of the axial rotations are not considered** [E].
   - Moving the N phases into the lepton doublets generates (θ/2)N_g(α₂/8π)WW̃, or equivalently derivative couplings to
     anomalous currents.
   - The effective rate is κ ≈ (α₂N_g/4π)θ̇ ≈ 10²H_i at θ̇ = χ̇/f ≈ 2×10⁴ H_i, comparable to m_χ. This compares with
     ξ_A ≈ 0.3 for inf_reheat's S2.
   - SU(2) gauge preheating (cf. Adshead et al. arXiv:1502.06506 [R] for the U(1) case) and spontaneous lepto- or
     baryogenesis may therefore not be negligible. They could add an SM channel that relieves N1 and the Higgs, or
     overproduce the baryon asymmetry. Not computed here or in the track.
7. **N7. Hidden heat from s** [E/Q]. The induced holder portal A_sΣ s|Σ|² gives R_leak ≈ 0.4–6 unless it vanishes at the
   stated renormalisation scale (§2.4).
8. **N8. N_s is under-specified** [D/E]. It needs a soft Majorana mass and a seesaw Yukawa; otherwise it is dark
   radiation (ΔN_eff ≈ 0.04). The N_s loop sets m_s² ≈ −0.0045(12L − 4) H_i² during inflation, which is ≤ 0.01 H² only
   for a tuned μ (§2.4).
9. **N9. The structural bound S1 is κ-dependent.** Heavy-waterfall clocks with m_σ ≈ 30–150 H_i pass, so G3/G4's quiet
   conclusion is not forced (§2.3).
10. **N10. Export defects** (§2.7): f_s is potential-only for 0.6 < t < 60; R switches from raw to cycle mean at t = 0.6;
    the spectator fraction at decay is 1.3%, not 2.4%.
11. **N11. If χ fragments, HFa's waterfall-mass spectator is inhomogeneous**, and the homogeneous release does not apply.
    HFb (ξ_s R) is the more robust variant.
12. **N12. Gravitational X production is ×2 low** (1.1×10⁻³ of E_X). This changes no conclusion.

---------------------------------------------------------------------------------------------------------------------

## 4. What stands

- **CMB.** The CMB numbers of HFa and HFb stand. n_s, A_s, r and α_s are insensitive to K over 2–8, apart from a 1–2%
  retune of A_s.
- **Waterfall timing.** ε_H = 1 at t = 0.19 with H = 0.997 H_i (an abrupt end); R reaches −4.5 to −4.75. There are no
  PBH: the tachyonic phase lasts ~0.2 H_i⁻¹, and ζ is ~4.6×10⁻⁵ on super-horizon scales.
- **Structure.** The valley argument (σ holds V0 when m_χ ≫ m_σ) and the perturbative bound S2 stand.
- **Axial production.** The single half-period axial production stands, including LZ within ×1.15–2.3. It is exact.
- **CP sign and decays.** The CP-sign bookkeeping (α + sign(s_*)) and the open or closed status of the decay channels
  stand.
- **Spectator.** The spectator's freezing, t_on ≈ 0.2 H_i⁻¹ and the β_iso formula stand.
- **SM plasma onset.** A Higgs-coupled SM component (ρ_SM/ρ ≥ 10⁻³) appears early, at τ ≈ 0.3–0.6, even with Pauli
  blocking. Whether that is enough for the Higgs is mf_higgs's question.

## 5. Limits of this pass

- The Pauli saturation is computed for a harmonic, constant-amplitude oscillation. The real χ oscillation about v is
  anharmonic and redshifting.
- The fluid's blocking model (in-band energy, band-edge shrinkage, unblocking by N decays plus Γ_th) is my [E] model on
  that [C] input.
- The waterfall is homogeneous plus linear modes. I ran no lattice and included no metric perturbations, so the wall
  outcome is not computed, only shown to be undecided.
- Γ_th (N thermalisation) is bracketed, not computed.
- mf_ideal's p and leak ledger are quoted, not recomputed.
- The anomaly item (N6) is an order-of-magnitude [E].
- I did not run the release or the chain. The branch is not prejudged.

## 6. Code and results (`skeptic/code`, `skeptic/results`; one core, ≈ 12 min in total)

| file | content |
|---|---|
| `sk_common.py` | independent potential from the exported parameters; vacuum; JSON helpers |
| `sk_cmb.py` → `sk_cmb.json` | Mukhanov–Sasaki (v-equation) + tensors; N_* bookkeeping; K sensitivity |
| `sk_structure.py` → `sk_structure.json` | e-fold capacity (track rows + 21-κ scans with V0-dominance and regularity cuts); S2 variants |
| `sk_axial.py` → `sk_axial.json`, `sk_axial_Fsat_vs_A.json` | exponential propagator; half-period regression; multi-period Pauli saturation |
| `sk_pauli.py` → `sk_pauli.json` | per-species reheating fluid with and without blocking; K, fractions, SM times |
| `sk_waterfall.py` → `sk_waterfall_HFa.json`, `sk_waterfall_HFb.json` | two-field waterfall with coupled modes; three friction models; wall ratios |
| `sk_bench_check.py` → `sk_bench_check.json` | audit of the exported histories |
| `sk_misc.py` → `sk_misc.json` | gravitational X, wall timescales, Σm_ν, isocurvature, N_s loop, A_sΣ leak, m̃ variant |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_hybrid/skeptic/code
python3 -B sk_cmb.py; python3 -B sk_structure.py; python3 -B sk_axial.py; python3 -B sk_axial.py fsat
python3 -B sk_pauli.py; python3 -B sk_waterfall.py HFa; python3 -B sk_waterfall.py HFb
python3 -B sk_bench_check.py; python3 -B sk_misc.py
```

## Literature (checked by web search in this session)

- Planck 2018 X, arXiv:1807.06211: β_iso < 0.038 for axion I (uncorrelated, n_II = 1), k_mid = 0.05 Mpc⁻¹.
- Planck 2018 VI, arXiv:1807.06209: Σm_ν < 0.12 eV (with BAO).
- ACT DR6, arXiv:2503.14452: n_s = 0.974 ± 0.003 (P-ACT-LB).
- BICEP/Keck 2021, arXiv:2110.00483: r_0.05 < 0.036.
- DESI DR2 neutrino constraints, arXiv:2503.14744: Σm_ν < 0.0642 eV.
- Greene–Kofman, hep-ph/9807339 and hep-ph/0003018: fermionic preheating and Pauli blocking.
- Casas–Ibarra, hep-ph/0103065.
- Adshead–Giblin–Scully–Sfakianakis, arXiv:1502.06506: gauge preheating.
- Recalled, already checked by the track: Felder et al. hep-ph/0012142 and hep-th/0106179.
