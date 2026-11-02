# Track inf_chain: the inflaton benchmarks pushed through the full M1r chain

3 November 2026 · Claude, inflaton-sector pass (BRIEF_INFLATON.md of this directory; GitHub issue #13 context). This is
same-session AI work on one supplied model, not outside replication. Nothing here is evidence for an inflaton, the hidden
sector, the holder or any dark-matter branch, and nothing is a validation or an exclusion. **Every potential, slice
convention, holder calibration point, decay width, renormalisation condition and ledger is a model choice; none follows
from closure.** Data comparisons (Planck 2018 A_s and n_s, BICEP/Keck r) are retrospective.

Labels: [D] derived, [C] computed with our own code (`code/`, output in `results/`), [E] estimate, [Q] quoted from a
reviewed track, [A] argued, [R] recalled or literature.

Git was read-only and nothing was committed. All runs were in the foreground, with at most two worker processes. Two
process notes:
- One shell wait loop (its `pgrep` matched its own command line) exceeded the call limit and was moved to the background
  by the harness. The engine job it was waiting for had already finished and was written to disk. The loop was later
  stopped at its time limit, and no work was lost.
- Two parallel kinetic runs overwrote each other's entries in one JSON file. They were rerun one after the other.

## 0. Result in brief

**Which benchmarks.** The potential track recommended HT6t (primary) and T1 (stress test); E1 was optional (same inflation
as T1) and is not run. The reheating track built quiet T-model members; I take its B1 (α = 0.3) and B2 (α = 1) as TB1 and
TB2, re-cast on this chain's slice convention.
- **Slice and holder.** For every benchmark the slice and the holder sit at S1, one e-fold before ε_H = 1, so a = 1 and
  H = H_i there. The preparation is re-solved at the slice, with R = 6(2 − ε_H)H_i², and g_I is re-fixed so that
  μ_eff(S1) = 10i.
- **The Send bracket.** HT6t_Send uses the holder calibrated at ε_H = 1 (g_I × 3.29).
- **The one-step iteration.** HT6t_cal1 and T1_cal1 re-calibrate the potential at the H_i that the dark-matter calibration
  demanded.

**The per-benchmark table** (e = 10⁻⁶, portals vanishing at M_Y, induced-leak channel excluded; see §5):

| | HT6t | HT6t_cal1 | HT6t_Send | T1 | T1_cal1 | TB1 | TB2 |
|---|---|---|---|---|---|---|---|
| potential | tilted p = 6 hilltop | same, re-tuned | HT6t | T-model α = 6.10×10⁻⁵ | α = 3.79×10⁻⁷ | T-model α = 0.3 | T-model α = 1 |
| H_i = H(S1) of the potential (GeV) | 1.290×10¹¹ | 2.805×10¹⁰ | 1.290×10¹¹ | 1.290×10¹¹ | 1.043×10¹⁰ | 7.47×10¹² | 1.060×10¹³ |
| m_I/H_i; Γ_I/H_i | 42.4; 0.5 | 119.9; 0.5 | 42.4; 0.5 | 128; 0.5 | 1625; 0.5 | 2.21; **7.5×10⁻⁹** (cal.) | 1.54; **2.1×10⁻⁹** (cal.) |
| **S/Q at e = 10⁻⁶** | 155 | 618 | 302 | 646 | 447 | 615 | 987 |
| S/Q without the Born φ_I → YY heat | 150 | 73 | 302 | 121 | 85 | 615 | 987 |
| S/Q at e → ∞ | 61 | 325 | 134 | 337 | 228 | 5.1 | 7.9 |
| **branch at 10⁻⁶** | free X | free X | free X | free X | free X | free X | free X |
| e_crit (bulk 3.47; Q = 2×10⁵) | none at any e | none | none | none | none | none | none |
| relic | X + X̄ gas, p = 8.1 X̄/charge, M = 17 m_X | p = 34, M = 69 | p = 16, M = 34 | p = 36, M = 72 | p = 24, M = 50 | p = 34, M = 69 | p = 55, M = 110 |
| H_i required by Ω_DM (GeV) | **2.70×10¹⁰** | 1.04×10¹⁰ | 1.26×10⁹ | **9.9×10⁹** | 1.19×10¹⁰ | 7.47×10¹² (fixes Γ_I) | 1.060×10¹³ (fixes Γ_I) |
| consistency: required/potential | **0.21 → iterate** | 0.37 (still falling) | 0.010 | **0.076 → iterate** | 1.15 (near a fixed point) | 1.000 | 1.000 |
| H_* at pivot exit (GeV) | 1.290×10¹¹ | 2.805×10¹⁰ | 1.290×10¹¹ | 1.290×10¹¹ | 1.043×10¹⁰ | 8.99×10¹² | 1.611×10¹³ |
| r (with H_*) | 2.71×10⁻⁷ | 1.28×10⁻⁸ | 2.71×10⁻⁷ | 2.71×10⁻⁷ | 1.77×10⁻⁹ | 1.32×10⁻³ | 4.22×10⁻³ |
| n_s (Planck 2018: 0.9649 ± 0.0042) | 0.9649 | 0.9649 | 0.9649 | 0.9618 | 0.9608 | 0.9619 | 0.9622 |

Labels:
- Release, vectors, ledger arithmetic, S/Q and the CMB numbers are [C] given the model choices.
- The ledger itself is dq_heat's [Q/E].
- p comes from a conditional kinetic recipe [C on a recipe; E for its validity], and so does every H_i that depends on it.

**Five findings**

1. **No benchmark gives an M1r-like branch point [C].** On every benchmark the declared e = 10⁻⁶ is free X, with S/Q 45–280
   times the threshold 3.47.
   - S/Q stays above 3.47 at every e. With vectors and A_L deposits removed (e → ∞), the Y_L deposits alone give S/Q of
     5.1–337.
   - No benchmark therefore has a Q-matter branch or an e_crit. This confirms the inf_potential skeptic's result on our own
     converged releases, own vectors and own K. We find S/Q = 155 (HT6t) and 646 (T1); the skeptic had 149 and 561.
2. **Two different mechanisms do it.**
   - **HT6t and T1 (fast oscillation, Γ_I = 0.5) [C].** The release is non-adiabatic, because m_I is comparable to the Y
     masses. It leaves a homogeneous Y_H remnant of N_YH = 1483 / 1649 (a = 1 at S1), i.e. 74 / 82 at ε_H = 1, against
     M1r's 0.0104. The cascade A_L from that remnant dominates the heat (R_casc = 627 / 692 E_X).
     - For T1, and for HT6t once it is re-calibrated, the open fundamental φ_I → Y_LY_L adds R_YY ≈ 5.6–8.3×10³ E_X.
   - **TB1 and TB2 (quiet: m_I = 2.2 and 1.5 H_i) [C].** The release is M1r-like: N_YL = 0.30, N_YH = 0.008, ε = 0.083–0.086.
     The realisable perturbative width makes the difference.
     - The dark-matter calibration fixes Γ_I ≈ 2–8×10⁻⁹ H_i, i.e. Γ/m_I ≈ 1.4–3.4×10⁻⁹. That gives T_RH = 1.3–2.0×10¹¹ GeV,
       a matter era of 12–13 e-folds and K ≈ 2000–2600 (a = 1 at S1).
     - The relic vectors cross H = m_A in the matter era (R_vec ∝ e^{−1.42}, computed). Their deposit carries K^{1/2}, so
       R_vec(10⁻⁶) = 1.8–3.2×10³ E_X.
3. **The relic is a mostly symmetric free-X gas [C on a conditional recipe].** I extended the pair-fraction map with the
   reviewed case-L moment network to S/Q = 17–1308. It gives p ≈ 0.05 S/Q, linear above S/Q ≈ 500.
   - Running the network's visible clock with the calibrated charge yield (Y_Q × 0.9374/M) raises p by 6–10%.
   - The mass per charge M = 1 + 2p is 17–110 m_X. The asymmetric-relic calibration is then a mostly symmetric calibration:
     H_i ∝ M^{−2/3}.
4. **Consistency and the iteration [C].**
   - **HT6t** needs H_i = 2.7×10¹⁰ GeV, not 1.29×10¹¹. Re-calibrated there (v = 0.1226 M_P, β = 1.235×10⁻⁶, still
     n_s = 0.9649), m_I rises to 120 H_i. φ_I → Y_LY_L opens and the requirement falls again, to 1.04×10¹⁰. The single step
     moves away from consistency, and HT6t has no consistent point at or near its design H_i.
   - **T1** needs 9.9×10⁹; re-calibrated there (α = 3.79×10⁻⁷, m_I = 1625 H_i) it needs 1.19×10¹⁰. That is a nearly
     self-consistent point at H_i ≈ 1.1×10¹⁰ GeV, with r ≈ 2×10⁻⁹. T1 fragments (Lozanov–Amin arXiv:1608.01213), so its
     homogeneous release is indicative only.
   - **TB1 and TB2** are consistent by construction, because the width absorbs the calibration. One N_* step with the actual
     reheating history converges: N_* = 52.35 / 52.61, n_s = 0.9619 / 0.9622, r = 1.32 / 4.22×10⁻³.
   - Their needed coupling y ≈ 2.6–4.1×10⁻⁴ (Majorana φN̄N, N_f = 1) lies inside the reheating track's radiative window,
     y_nat ≈ 4–9×10⁻⁴ [Q inf_reheat].
5. **A new heat source for slow reheating dominates the T-model ledger [E].** Inflaton → A_LA_L proceeds through the induced
   φ_I–σ mixing (the A_IΣ portal), and here it is counted with the long-lived-vector ledger.
   - A slowly decaying inflaton holds about 10¹² quanta per H_i³, so a branching of 10⁻¹¹ already matters.
   - The renormalisation-condition bracket gives R_leak = 6×10², 1×10⁴ and 2×10⁶ E_X for TB1, against R_vec = 1.8×10³.
   - It is negligible for every Γ_I = 0.5 benchmark, except HT6t with M_P portals (m_I ≈ m_σ, near-resonant mixing):
     R_leak = 1.3×10³.
   - It only pushes further into free X, but it lowers the T-models' H_i,req at fixed Γ by up to ×30. The induced A_IΣ must
     be tuned below about 10⁻⁶ H_i for the TB numbers above to stand.

**Bottom line [A].** With an explicit inflaton the M1r chain no longer has its decisive either/or at e = 10⁻⁶. Every
benchmark is deep in free X, with a relic that is mostly X X̄ pairs. The calibrated H_i lands either:
- far from the potential's design value (HT6t, T1 at 1.29×10¹¹), or
- at H_i ~ 10¹³ GeV with r ~ 10⁻³ (the quiet T-models with their realisable slow decay).

The Q-matter branch would need, simultaneously, a quiet inflaton (m_I ≈ 2–4 H_i), fast reheating (K ≈ 2) and slow-roll
pre-slice holding. No benchmark here, and no realisable perturbative channel of the reheating track, supplies all three.

## 1. Benchmarks and conventions (`code/ch_common.py`, `code/benchmarks_reg.py`)

| key | potential | origin | slice / holder | Γ_I |
|---|---|---|---|---|
| HT6t | V0(1 − βx − x⁶)²/f_top², v = 0.34694 M_P, β = 1.5529×10⁻⁵ | inf_potential bench_HT6t.json [Q] | S1 / S1 | 0.5 H_i (the engine convention the potential track used) |
| HT6t_Send | same | — | S1 / ε_H = 1 (μ_eff(S1) = μ + 3.29(10i − μ) = −16.4 + 32.9i) | 0.5 |
| HT6t_cal1 | v = 0.122614 M_P, β = 1.23528×10⁻⁶ at H_i = 2.805×10¹⁰ (A_s = 2.100×10⁻⁹, n_s = 0.9649) | `recal.py` [C] (inf_potential calibrate.py, read-only) | S1 / S1 | 0.5 |
| T1 | V0 tanh²(φ/√(6α)), α = 6.0968×10⁻⁵ | inf_potential [Q] | S1 / S1 | 0.5 |
| T1_cal1 | α = 3.7855×10⁻⁷ at H_i = 1.043×10¹⁰ | `recal.py` [C] | S1 / S1 | 0.5 |
| TB1 | T-model α = 0.3 (inf_reheat B1) | V0 from A_s at N_* (iterated, §4) | S1 / S1: ε_H = 0.135, R = 11.19 | calibrated, 7.5×10⁻⁹ H_i |
| TB2 | T-model α = 1 (inf_reheat B2) | as TB1 | S1 / S1: ε_H = 0.227, R = 10.64 | calibrated, 2.1×10⁻⁹ H_i |

How to read this table:
- The hidden sector is M1r's at λ_Y = λ_Σ = 0.5, in units of H_i = H(S1) [Q]. Its dimensionless history therefore depends
  only on the inflaton shape in H_i units.
- **Holder couplings** (g_I; α_I): HT6t 2.56×10⁻⁶ (−0.949); HT6t_cal1 1.38×10⁻⁶; T1 6.79×10⁻⁶ (2.192); T1_cal1 4.63×10⁻⁶;
  TB1 2.46×10⁻⁵ (2.192); TB2 2.95×10⁻⁵ (2.192).
- **Prepared states at the slice** (|X|, |Y|, |Σ|) in H_i, all heavy:
  - HT6t / T1: (90.78, 6.884, 28.576), lightest m² = 10.19;
  - TB1: (85.17, 5.864, 28.497), 10.17;
  - TB2: (81.50, 5.226, 28.453), 10.15;
  - HT6t_Send: Y-condensed, (508.2, 96.1, 73.5), lightest m² = 163.
- TB slice choice [A]. inf_reheat put its slice at ε_H = 1 (H_i := H_end). There the static preparation is not adiabatic,
  and R = 6H² halves the X hold. Starting at S1 lets the engine evolve the last e-fold itself. In this convention TB1 has
  m_I = 2.21 H_i, inside the inf_potential skeptic's quiet window (2–4 H_i); TB2 has 1.54.

## 2. Item 1: the release to freeze-out (`code/rpchain.c`, `code/release.py`) [C]

**Engine.** `rpchain.c` is a copy of inf_potential's `rpinf.c`, which is itself a copy of dq_entropy's `rpdq.c` plus an
explicit-potential MODE 4. One addition, CH1: an optional step-cap relaxation after the inflaton has decayed, used only in
long runs from τ = 60.
- With CH1 unset, the output is byte-identical to rpinf.c on HT6t and T1.
- The slice setup reproduces inf_potential's engine inputs to every digit (G, I0).

**Runs.**
- A dense run to τ = 100 (dt = 0.005) gives the vector background.
- A long run to τ = 3×10⁴ (Γ = 0.5) or 2×10⁴ (TB, Γ = 0 in the engine: N_X does not depend on Γ for Γ ≲ 10⁻³ [Q inf_reheat])
  uses conformal bins.
- The X ↔ Σ exchange is fitted as N_inf + c₂a⁻² over windows of 4 in η, with dq_entropy's alternatives.
- The fluctuation decays as a^{−2.00…−2.21} (TB1: −1.60), i.e. the a⁻² law holds.

| a = 1 at S1 | HT6t | HT6t_cal1 | HT6t_Send | T1 | T1_cal1 | TB1 | TB2 | M1r [Q] |
|---|---|---|---|---|---|---|---|---|
| N_X (frozen) | −75913 ± 2 | −91760 ± 2 | **+3.836×10⁶** ± 5 | −112285 ± 0.3 | −145345 ± 0.3 | −8554 ± 0.8 | −4635 ± 0.3 | −8047.8 |
| E_X (conformal) | 1.981×10⁶ | 2.439×10⁶ | 1.579×10⁸ | 2.228×10⁶ | 2.849×10⁶ | 3.477×10⁵ | 1.607×10⁵ | 1.887×10⁵ |
| ε = N_X/N_c | 0.2061 | 0.2131 | 0.3904 | 0.2790 | 0.3004 | 0.0856 | 0.0827 | 0.1274 |
| E/E_c | 8.216 | 7.857 | 3.505 | 5.485 | 4.970 | 26.50 | 27.73 | 15.60 |
| K = a²H (late) | 12.341 | 13.645 | 12.341 | 15.200 | 17.108 | fluid: 2000 / 2589 at Γ_cal | | 1.986 |
| N_YL / N_YH / N_σ | 105 / 1483 / 0.22 | 66 / 676 / 0.79 | 5.85×10⁵ / 3.21×10⁵ / 1.68×10⁵ | 234 / 1649 / 5.0 | 23 / 1401 / 6.1 | 0.301 / 0.0079 / 4×10⁻⁶ | 0.307 / 0.0072 / 6×10⁻⁶ | 0.567 / 0.0104 / 2×10⁻⁵ |
| Y/Y_Q at H_M1 | 0.6034 | 0.6274 | 30.50 | 0.6529 | 0.7078 | 3.3×10⁻⁵ (Γ_cal) | 1.2×10⁻⁵ | 0.9910 |
| Gauss / hidden work | 5×10⁻⁹ / 1×10⁻¹⁰ | 3×10⁻⁹ / 1×10⁻¹⁰ | 1×10⁻⁸ / 9×10⁻¹¹ | 6×10⁻⁹ / 7×10⁻¹¹ | 4×10⁻¹¹ / 2×10⁻¹⁰ | 4×10⁻¹¹ / 3×10⁻¹¹ | 6×10⁻¹¹ / 4×10⁻¹¹ | — |

Notes on the table:
- **Converting to the potential track's normalisation** (a = 1 at ε_H = 1): divide comoving numbers by e³, E_X by e⁴ and K
  by e².
  - HT6t: N_X = −3778, E_X = 3.63×10⁴, N_YL = 5.23, N_YH = 73.8. This reproduces the skeptic's smoke numbers (−3768, 5.23,
    73.84). T1: N_YH = 82.1.
- **K for HT6t.** It is 1.670 in that normalisation, against inf_potential's 1.516. The engine applies the Γ term from the
  slice, through the last e-fold, which the potential track's Python history did not. I use the engine's K throughout.
- **H-constraint check.** The relative H check grows to 3×10⁻⁵ – 9×10⁻³ by τ = 3×10⁴. This is a constant absolute offset of
  the auxiliary H, ≤ 1.5×10⁻⁷ H_i, divided by the falling H. It does not affect K, which is computed from the constraint.
- **HT6t_Send** flips the sign of the charge relative to S1 at the same α_I. Its Y-condensed preparation makes a 50 times
  larger, opposite-sign release [C].

## 3. Item 2: vector heat with the benchmark's H(N) and F(N) (`code/vmodes.c`, `code/vectors.py`) [C]

**Method.**
- **Mode equation.** The longitudinal Proca modes use the Stueckelberg angle in cosmic time (the dq_heat skeptic's action),
  with p = a³F²k²/(k² + a²e²F²) and q = aF²k².
- **Integration.** Each mode is started as an exact massless de Sitter Goldstone at k/aH = 50, with the slice's prepared F.
  It is integrated in C (Dormand–Prince) through the recorded engine background (dense τ ≤ 100). The continuation is exact
  radiation (Γ = 0.5) or my matter-plus-radiation fluid (TB).
- **Fluid check.** The fluid reproduces the inf_reheat skeptic's constant K Γ^{1/3} = 0.50128 x₀^{2/3} to 3×10⁻⁶ at small Γ.
- **Read-out.** The adiabatic number is read once Ω/H > 100, and checked at 200. Weighted modes (k ≤ 1) converge to ≤ 1%; the
  k ≥ 10 tail is noisy but carries < 0.3% of N_A.
- **Exit factor.** The pre-slice history enters through the dq_heat skeptic's frozen-angle factor with the benchmark's own
  exit-time F and H: f(N) = [h_r F_s/F(s, h_r, ε)]². F is re-solved along the benchmark's path (s = φ_I/φ_I(S1),
  h_r = H/H_i, R = 6(2 − ε)h_r²), and the mode k exits at k = e^{−N}h_r(N). Exit states stay heavy (m²/H² ≥ 8.7).

**Validation on M1r (MODE 0, (3, 0.5)).** Against dq_heat's spectrum:
- N_A(k ≤ 1) agrees to 0.05% / 0.13% / 0.23% / 0.33% at e = 10⁻⁷ / 10⁻⁶ / 3×10⁻⁶ / 10⁻⁵;
- the per-mode n_k agrees within [0.987, 1.005] for k ≤ 0.3;
- the release band (k ≥ 1) gives 0.0527, against 0.0512.

| e = 10⁻⁶, a = 1 at S1 | HT6t | HT6t_cal1 | HT6t_Send | T1 | T1_cal1 | TB1 (Γ_cal) | TB2 (Γ_cal) |
|---|---|---|---|---|---|---|---|
| N_A | 28.1 | 34.4 | 1124 (k ≥ 1: 1122) | 36.5 | 46.4 | 58.5 | 42.7 |
| f_eff (exit factor, k ≤ 1) | 0.909 | 0.959 | 0.780 | 0.866 | 0.925 | **0.512** | **0.468** |
| deposit per vector / E_X | 0.42 | 0.36 | 0.0053 | 0.42 | 0.35 | 30.7 | 75.5 |
| R_vec (own E_X) | 11.9 | 12.4 | 6.0 | 15.2 | 16.1 | **1796** | **3226** |
| R_vec scaling in e | e^{−1.00} | e^{−1.00} | — | e^{−1.00} | e^{−1.00} | **e^{−1.42}** (matter-era crossing) | e^{−1.42} |

Reading:
- The HT6t and T1 f_eff match inf_potential's (0.910 / 0.867).
- The T-model members roll more before S1: s = 1.4–1.7 at 2–4 e-folds. That halves f_eff, the "rolling lowers the vector
  heat" effect of the dq_heat skeptic.
- Their K ≈ 2000 matter-era history raises the deposit by K^{1/2}.
- For TB, N_A is independent of Γ_I within 2% for Γ_I ≤ 10⁻⁴ H_i, because the crossing happens in the matter era.
- In HT6t_Send the violent release of the Y-condensed state produces sub-horizon Goldstones, N_A(k ≥ 1) = 1122 (in the
  single-Higgs-direction model).

**Mediators and cascade.**
- **Y_L deposit.** Each Y_L deposits m_L⟨a_L⟩, with ⟨a⟩ integrated on the benchmark's own a(t). For TB, Γ_I < Γ_L = 10⁻⁷
  H_i, so Y_L decays inside the matter era and R_Y saturates at about 2–4.
- **Y_H.** Each Y_H gives one Y_L plus one A_L, and σ gives 2 A_L (dq_heat ledger [Q]).
- **Born φ_I → YY.** I recount the holder-induced pairs from S1 (`yy_born.py`; same Born formulas as inf_potential,
  verified by its skeptic).
  - N_YL / N_YH = 67 / 46 (HT6t), against inf_potential's e³-converted 59 / 40;
  - 1.23 / 1.12×10⁴ (T1), against 1.05 / 0.96×10⁴;
  - 9.7 / 9.0×10³ (HT6t_cal1), 9.3 / 9.3×10³ (T1_cal1);
  - zero for TB (all harmonics below 44.5 H_i).

## 4. Item 3: S/Q, branch, relic, H_i and r (`code/chain.py`, `code/calib_tb.py`, `code/pmap_kin.py`, `code/recal.py`)

**S/Q [C].** I use dq_entropy's closed form S/Q = 0.190255(W − 0.97494)^0.75013 − 0.00428, with W = (E/E_c)(1 + R).
- It was fitted on W = 13–140. I checked it here up to W = 1.3×10⁵ (`sq_check.py`): max |rel| = 6×10⁻⁴.
- The check uses the copied relaxed-Hartree solver on the condensed branch up to W = 5400, where the circular condensate
  evaporates (condensed fraction 0.045 at W = 5000). Above that it uses our own no-condensate Hartree branch
  (`hartree_nocond.py`; μ/M = 0.99 → 0.09 over W = 5400 → 1.3×10⁵).

**S/Q(e)** (own E_X units; the e → ∞ column keeps the Y_L deposits only):

| | 10⁻⁷ | 10⁻⁶ | 10⁻⁵ | 10⁻⁴ | 10⁻³ | → ∞ |
|---|---|---|---|---|---|---|
| HT6t | 320 | 155 | 94 | 72 | — | 61 |
| HT6t_cal1 | 1138 | 618 | 425 | 358 | — | 325 |
| HT6t_Send | 588 | 302 | 192 | 153 | — | 134 |
| T1 | 1192 | 646 | 442 | 372 | — | 337 |
| T1_cal1 | 835 | 447 | 303 | 252 | — | 228 |
| TB1 (Γ_cal) | 7416 | 615 | 55 | 8.5 | 5.3 | 5.1 |
| TB2 (Γ_cal) | 1.2×10⁴ | 987 | 88 | 13.7 | 8.2 | 7.9 |

There is no e_crit (bulk or Q = 2×10⁵) in any row.
- The only Q-matter seen anywhere was TB1 at a fast width (Γ_I ≥ 10⁻⁴ H_i, K ≤ 85): S/Q(e → ∞) = 3.28, with e_crit = 3×10⁻⁴.
- That width misses the dark-matter calibration by a factor 10 in H_i [C].

**Relic: pair fraction at large heat [C on a conditional recipe, E].** wt_critic's gas-only map ends at S/Q = 6.74. I ran the
reviewed case-L moment network (fx_chain's M1-patched copy, imported read-only; gas only; forced-equilibrium a_ref start)
at R_spot = 30–10⁴.

| S/Q | 17.1 | 21.8 | 41.6 | 94.4 | 233 | 530 | 1308 |
|---|---|---|---|---|---|---|---|
| p (X̄ per charge) | 0.503 | 0.744 | 1.78 | 4.54 | 11.7 | 26.9 | 65.8 |

- The anchor at S/Q = 21.8 reproduces REVIEW_KINETIC's 0.736 (our M1 inputs: 0.744).
- The slope tends to 1 (d ln p/d ln(S/Q) = 0.99 over 530–1308), i.e. symmetric freeze-out at fixed T_h/m_X.
- **Self-consistent yield.** With M ≫ 1 the calibrated charge yield is Y_Q × 0.9374/M, which makes the hidden sector colder
  relative to the visible bath. Running the network with that yield raises p by 6.4 / 7.5 / 8.2 / 8.7 / 8.8 / 9.6% at
  S/Q = 155 / 302 / 447 / 600 / 646 / 958. This is used in every row above.
- The recipe's limits are unchanged: forced equilibrium start, no objects, MB fugacities. A factor 2 in p moves H_i by
  ×0.63–1.6 (bracket in `summary.json`).

**H_i calibration [D on C].** H_i,req = H_M1(0.9374/(M·Y/Y_Q))^{2/3}, with Y/Y_Q = |N_X|H_M1³/(a³s) at H_M1. This is
fx_holder's convention: Y ∝ H_i^{3/2} at fixed dimensionless history.
- For HT6t/T1 the history depends on H_i through the potential (m_I/H_i, φ_I(t)/φ_I(S1)), so a mismatch requires a
  re-calibrated potential and a new release. That is the iteration.

| | potential H_i | required H_i | after one step |
|---|---|---|---|
| HT6t | 1.290×10¹¹ | 2.70×10¹⁰ (p × 2: 1.7–4.1×10¹⁰) | HT6t_cal1 at 2.805×10¹⁰ → 1.04×10¹⁰ (not converged; m_I 42 → 120 H_i opens φ_I → Y_LY_L) |
| T1 | 1.290×10¹¹ | 9.9×10⁹ | T1_cal1 at 1.043×10¹⁰ → 1.19×10¹⁰ (near a fixed point at ≈ 1.1×10¹⁰) |
| TB1 / TB2 | fixed by A_s at α | equal by the choice of Γ_I | N_* step: 52.38 → 52.35 / 52.65 → 52.61; Γ_cal 7.5 / 2.1×10⁻⁹ H_i converged |

**r and n_s [C].** r = 2H_*²/(π²A_sM_P²) uses H_* at pivot exit from each potential.
- **HT6t and T1.** H_* = H_i to 3×10⁻⁴ [Q inf_potential], so r = 2.71×10⁻⁷ at the design H_i. After the iteration r falls to
  1.3×10⁻⁸ (HT6t) and 1.8×10⁻⁹ (T1); at T1's near-fixed point r ≈ 2×10⁻⁹.
- **TB1 / TB2.** The plateau rolls, so H_* = 1.20 / 1.52 H_i: r = 1.32 / 4.22×10⁻³ from exact Mukhanov–Sasaki spectra at
  the iterated N_*. n_s = 0.9619 / 0.9622, with α_s ≈ −7×10⁻⁴.
- **Retrospective n_s comparisons [R, via the inf_potential skeptic].** Against Planck 2018 every benchmark is within 1σ.
  Against ACT DR6 P-ACT-LB (n_s = 0.9743 ± 0.0034, arXiv:2503.14452), T1/TB are at −3.6 to −3.7σ and HT6t at −2.8σ. Not an exclusion.

**Reheating of TB [C, with Q for the channel].**
- Γ_I/m_I = 3.4 / 1.4×10⁻⁹, T_RH = 2.0 / 1.3×10¹¹ GeV, with 12.2 / 13.0 matter e-folds after ε_H = 1.
- A singlet Majorana channel needs y = (16πΓ/m_I)^{1/2} = 4.1 / 2.6×10⁻⁴, inside inf_reheat's radiatively natural window
  (4–9×10⁻⁴; S1 or S2 [Q]).
- Kinematic blocking lifts long before the decay (y|φ_I| ≈ 10¹⁵ GeV against m_I/2 ≈ 10¹³ GeV at the slice; unblocked after
  ≈ 3 e-folds of redshift, while the decay comes after ≈ 12) [E].
- The Higgs-plasma timing at these H_i was not recomputed. inf_reheat's R-even portal λ_IH is the stated route [Q].

## 5. New issues

1. **Induced φ_I–σ mixing leaks inflaton quanta into long-lived A_L [E]** (`leak_estimate.py`).
   - **Mechanism.** The portal condition's A_IΣ φ_I|Σ|², with A = 2μRe(G)c, mixes φ_I with σ. σ(∂π)²/h then gives
     Γ = A²m_I³/[32π(m_σ² − m_I²)²], the inf_reheat skeptic's form. Each A_L deposits m_A⟨a_A⟩, about 10³ times its production
     energy.
   - **Bracket** c = 1/16π², 1/4π² (threshold size, M_Y condition) and 13.31/4π² (M_P condition):

   | benchmark | R_leak at 10⁻⁶ (bracket) | S/Q at 10⁻⁶ (bracket) |
   |---|---|---|
   | TB1 | 6.3×10², 1.0×10⁴, 1.8×10⁶ | 771, 2539, 1.1×10⁵ |
   | TB2 | 1.5×10³, 2.4×10⁴, 4.2×10⁶ | — |
   | HT6t | 0.5, 7, 1.3×10³ (near-resonant m_I ≈ m_σ) | — |
   | T1 | ≤ 21 | — |

   - At fixed Γ_I it lowers the TB H_i,req to 6.4×10¹², 2.9×10¹² and 2.3×10¹¹ GeV. Re-calibrating would raise Γ_I.
   - The quiet slow-reheating benchmarks therefore require |A_IΣ| ≲ 10⁻⁶ H_i, i.e. a tuned induced-portal condition, for
     their computed numbers to hold. This extends WORKING_THEORY's open input 3 (the portal renormalisation condition) from a
     17% effect to the dominant one.
2. **HT6t has no dark-matter-consistent point near its design H_i [C on one iteration step].** Lowering H_i raises
   m_I/H_i ∝ v⁻¹ and opens φ_I → Y_LY_L above 44.5 H_i.
   - The inf_potential "clean window" (30 < m_I < 44.5 H_i at H_i = 1.2–2.2×10¹¹) is incompatible with the free-X calibration
     once the relic is mostly symmetric.
3. **The Send calibration flips the sign of the charge [C].** HT6t_Send has N_X = +3.84×10⁶, against −7.6×10⁴ at S1, with
   the same α_I. With a held Y condensate at the slice, the CP sign is not fixed by α_I and sign(φ_I) alone; it also depends
   on where the holder is calibrated.
4. **The asymmetric-relic calibration becomes a symmetric one [C/E].** With S/Q ≳ 10², M = 1 + 2p ≈ 17–110. Two
   consequences:
   - The mass test of WORKING_THEORY §5 and the "H_i absorbs the mass per charge" statement now hinge on the gas-only
     freeze-out recipe at S/Q ~ 10²–10³, which has never been checked against a Boltzmann solution there.
   - The free-X relic is then roughly 49% X̄ by number. This is CDM-like and has no annihilation channel at rest (no lighter
     hidden states) [Q WORKING_THEORY §3].
5. **The S/Q closed form holds far beyond its fit range [C],** through the condensate evaporation (W ≈ 5400) up to
   W = 1.3×10⁵, to 6×10⁻⁴. Later tracks can use it there.
6. **The engine K includes the decay during the last e-fold [C].** For HT6t and T1 it is 10% (HT6t) and 11% (T1) above the
   potential track's Python history value, which switched the width on at ε_H = 1. Both are conventions; the yield moves
   by K^{3/2}.

## 6. Limits

- **Homogeneous release.** The release is homogeneous and classical, with linear vector modes in the single-Higgs-direction
  Proca model. T1 and T1_cal1 fragment within ≈ 0.4 H_i⁻¹ [Q inf_potential], although the Y_H remnant forms before that
  [Q its skeptic]. Their Born R_YY uses the homogeneous background and late-vacuum thresholds, uncertain by ≳ ×2 [Q skeptic].
- **The width Γ_I = 0.5 H_i** of the HT6t/T1 rows is kinematically blocked for the singlet-fermion channel [Q
  inf_potential skeptic]. A realisable slow width raises K, as in TB, and only increases S/Q.
- **Ledger.** It is dq_heat's: decay at rest, full deposit into the X plasma, branching 1 [Q/E]. Finite-k Y production
  with the new inflatons is not computed beyond the Born count.
- **Pair fraction.** p(S/Q > 6.7) comes from a conditional moment recipe (forced equilibrium start, gas only, MB
  fugacities), run far above the heat where it was reviewed.
- **Hilltop vacuum.** The HT6t rows inherit the potential track's choice of the φ = +v vacuum (no domain walls in the
  observable patch [A]).
- **Not recomputed here:** Higgs survival at the new H_i, and the TB slice at ε_H = 1 (inf_reheat's own convention).
- **Exit factor.** It assumes the hidden sector tracks its conditional minimum adiabatically along the pre-slice path
  (lightest m²/H² ≥ 8.7 on every path).

## Code (`code/`, own; read-only imports stated in each header)

| file | purpose |
|---|---|
| rpchain.c | engine: copy of inf_potential/code/rpinf.c (← dq_entropy rpdq.c) + CH1 step-cap relaxation |
| ch_common.py | benchmarks, slice preparation, holder, engine driver (imports inf_common, calibrate, dq_model_copy read-only) |
| benchmarks_reg.py | the benchmark registry |
| release.py | item 1 → results/release_*.json |
| vmodes.c, vectors.py | item 2: mode integrator; backgrounds, exit factor, N_A, R_vec → results/vectors_*.json |
| yy_born.py | Born φ_I → YY from the slice → results/yy_born_*.json |
| chain.py | ledger, S/Q(e), e_crit, relic, H_i,req → results/chain_*.json |
| hartree_nocond.py, sq_check.py | S/Q closed-form check above the condensate → results/sq_check.json |
| pmap_kin.py | case-L network at large heat, self-consistent yield → results/pmap_kin.json |
| recal.py | one-step re-calibration → results/recal_*.json |
| calib_tb.py | T-model width calibration and N_* step → results/calib_TB*.json |
| leak_estimate.py | induced-mixing leakage [E] → results/leak_estimate.json |
| summary.py | → results/summary.json |

## Reproduce (one or two cores; ≈ 30 min)

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_chain/code; export OMP_NUM_THREADS=1
gcc -O2 -o $SP/wt2/inf_chain/build/vmodes vmodes.c -lm          # rpchain is built by ch_common
for k in HT6t T1 HT6t_Send; do python3 -B release.py $k 100 30000; done     # ~2 min each
for k in TB1 TB2; do python3 -B release.py $k 100 20000; done
python3 -B vectors.py M1r                                                   # validation vs dq_heat
for k in HT6t T1 HT6t_Send; do python3 -B vectors.py $k; python3 -B yy_born.py $k; done
for g in 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3; do for k in TB1 TB2; do CH_ES=1e-7,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3 python3 -B vectors.py $k $g; done; done
python3 -B pmap_kin.py 30 41.707 100 300 1000 3000 10000; python3 -B pmap_kin.py sc 155 302 447 600 646 958   # sequentially
python3 -B sq_check.py
for k in HT6t T1 HT6t_Send; do python3 -B chain.py $k; done; for k in TB1 TB2; do python3 -B calib_tb.py $k; done
python3 -B recal.py HT6t_cal1 HT6t 2.805e10 0.121 1.26e-6; python3 -B recal.py T1_cal1 T1 1.043e10
for k in HT6t_cal1 T1_cal1; do python3 -B release.py $k 100 30000; python3 -B yy_born.py $k; python3 -B vectors.py $k; python3 -B chain.py $k; done
python3 -B vectors.py TB1 7.512e-9; python3 -B chain.py TB1 7.512e-9; python3 -B vectors.py TB2 2.113e-9; python3 -B chain.py TB2 2.113e-9
python3 -B leak_estimate.py; python3 -B summary.py
```

Scratch: `$SP/wt2/inf_chain/` (binaries in `build/`, run series in `runs/`, background tables).

## Credits

- Longitudinal-vector production: Graham–Mardon–Rajendran, arXiv:1504.02102 [R].
- α-attractors: Kallosh–Linde, arXiv:1306.5220 [R]. Squared hilltop: Boubekeur–Lyth, hep-ph/0502047 [R].
- Fragmentation: Lozanov–Amin, arXiv:1608.01213 [R].
- Data: Planck 2018 X, arXiv:1807.06211; BICEP/Keck, arXiv:2110.00483; ACT DR6, arXiv:2503.14452, quoted via the
  inf_potential skeptic [R].
- Greene–Kofman fermion production: hep-ph/9807339 [R].
- Our own tracks, reused read-only: inf_potential (rpinf.c, inf_common, calibrate, bench files, the Born formulas),
  inf_reheat (fluid constant, channel windows), dq_entropy (rpdq.c lineage, the closed form, the Hartree solver), dq_heat and
  its skeptic (ledger, mode action, frozen-angle exit factor, dq_model_copy), fx_chain (the kinetic network copy), wt_critic
  (p map) and dq_critic (Q-matter M table).
