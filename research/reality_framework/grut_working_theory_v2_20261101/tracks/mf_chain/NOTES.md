# Track mf_chain: the full M1r chain on the mf_hybrid benchmarks (HFa, HFb)

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md of this directory; GitHub issue #13 context).

This is same-session AI work on one supplied model, not outside replication. Nothing here is evidence for the hybrid,
the spectator, the hidden sector or any dark-matter branch. Nothing is a validation, an exclusion of the real world, a
new theory or a GRUT prediction. **Every history, rescaling rule, ledger entry, renormalisation condition and criterion
below is a model choice, and none follows from closure.** Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6) are
retrospective.

Labels: [D] derived, [C] computed with this track's runs (`code/`, output in `results/`), [E] estimate, [Q] quoted from a
reviewed track, [A] argued, [R] literature (arXiv numbers checked by web search in this session).

Git was read-only and nothing was committed. All runs were in the foreground. Binaries, run series and the corrected
history files live in scratch (`$SP/wt2/mf_chain/`); no compiled file or `__pycache__` is in the tree. Process note:
one early call ran two chain runs side by side, so for a few seconds the two vector stages could have used up to four
worker processes; every later run was sequential (at most two processes). One history build was stopped by its own
time-out (foreground) and replaced (§2.3).

## 0. Result in brief

1. **On every history tried, both benchmarks are free X at e = 10⁻⁶ and at 3×10⁻⁶, and there is no e_crit at any e [C].**
   - As exported (mf_hybrid's unblocked axial reheating, K = 2.05): S/Q(10⁻⁶) = 14.69 (HFa) and 14.58 (HFb); 12.25 /
     12.15 at 3×10⁻⁶; 10.86 / 10.77 at e → ∞.
   - With the mf_hybrid skeptic's Pauli blocking applied (K = 4.3–9.8 at the design H_i): S/Q(10⁻⁶) = 27–72, 3×10⁻⁶:
     18–34, e → ∞: 13–18.
   - **What forbids Q-matter even at e → ∞ is the light mediator Y_L**: N_YL = 5.4–11 (M1r: 0.57), R_Y = 9–15 E_X, while
     the bulk threshold allows a total heat R < 0.6–1.3. The vectors then add R_vec = 4.6 (K = 2) to 80 (K = 9.8).
2. **Self-consistent H_i (V0 rescaled, rule in §6) [C].**
   - Exported/unblocked: the dark matter needs 0.986 × the design value; the iteration converges at **H_i = 1.278×10¹¹
     (HFa) / 1.279×10¹¹ GeV (HFb)** in two steps (K = 2.04, r = 2.68×10⁻⁷, CMB capacity 75 ≥ N_* = 52.6).
   - Pauli-blocked (Γ_th = 0, perturbative tails free): 1.30 → 2.30 → 2.40×10¹¹ GeV, **H_i = 2.40×10¹¹ GeV** in three steps
     (K = 7.84, S/Q(10⁻⁶) = 61.5, r = 9.4×10⁻⁷). Brackets: Γ_th = 10 H_i gives 2.21×10¹¹ (K = 5.81); tails blocked while
     inside the Pauli band gives 2.51×10¹¹ (K = 9.64).
   - **At the Pauli fixed point HFa's clock is no longer CMB-flat at the benchmark's κ = 10** (X = g²σ_c grows ∝ H_i to
     0.53; e-fold capacity 48.6 < N_* = 52.3). κ = 15–50 restores it (61–102 with the skeptic's regularity cut). HFb
     (n_s = 0.9743) passes at κ = 10 (53.7).
   - Calibrated at e = 3×10⁻⁶ instead: 1.75×10¹¹ GeV (unblocked) and 3.92×10¹¹ GeV (Pauli). At the latter no κ gives a
     flat clock for (m_χ, m_σ) = (250, 600) H_i (N_max ≤ 50.2).
3. **Relic [C on the inherited p(S/Q) map].** Free X X̄, M = 1.56–1.58 m_X per charge (unblocked; 63–64% asymmetric) and
   4.6–8.8 m_X (Pauli; 11–22% asymmetric) at 10⁻⁶. Q-matter is never reached, so qm_micro's ⟨Q⟩ and M are not needed.
4. **Hidden heat from s decay [E]** (induced holder portal A_sΣ s|Σ|² → s mixes with the Σ radial mode → A_L A_L).
   R_leak = 0.5–0.7 at c_L = 1/16π², 8.5–12 at 1/4π², 1.5–1.9×10³ if the portals vanish at M_P. It only raises S/Q (to
   15–20 / 62–86 at threshold size, 460–720 with M_P portals). With c_L = 1/4π² included in the calibration, the unblocked
   HFa fixed point moves to H_i = 7.3×10¹⁰ GeV (three steps). σ and χ have no hidden channel (sequestered by choice). The
   s decay into the SM dilutes the yield by 0.1–0.9% only.
5. **Gravitational production at the abrupt end [C/E]:** finite-k X carries 2.8×10⁻⁴ (unblocked) to 5.6×10⁻³ (Pauli) of
   E_X; gravitational Y_L, Y_H, σ_Σ quanta ≤ 0.015 per H_i³. ΔR ≤ 0.25 against R_total = 15–80: negligible.
6. **Isocurvature [C].** D = d ln Ω_DM/d ln s_* = 1.36–1.38 (unblocked) and −0.07 to −0.01 (Pauli) at 10⁻⁶; 3.0 and
   0.46–0.59 at 3×10⁻⁶. The charge alone responds as N_X ∝ s_*^4.54–4.59, but the symmetric part (M) cancels most of it.
   At the exported s_* = 3×10⁻³ M_P = 7.3×10¹⁵ GeV, β_iso = 6.9×10⁻³ (unblocked, 10⁻⁶), 0.033 (3×10⁻⁶) and ≤ 4.5×10⁻³
   (Pauli), against Planck 2018's 0.038 (retrospective).
7. **Against mf_ideal's map** at the same K and m_s = 3: the benchmarks sit at 0.5–0.8 of the map's S/Q (map 29.8 / 18.4 for
   t_on = 0.1 / 1 at K = 2.05; extrapolated 104 / 93 at K = 7.8). The hybrid's end is gentler than an instantaneous R step
   (larger ε), but its fast spectator turn-on still excites Y_L hard.

**Plain answer.** The repaired hybrid does not recover M1r's e-controlled branch point: its fast end and fast spectator
turn-on leave a Y_L remnant that makes the relic free X at every e, and the correct (Pauli-blocked) reheating is slower,
which raises the vector heat and the calibrated H_i to ≈ 2.4×10¹¹ GeV, where HFa's clock needs a different
renormalisation condition to stay CMB-flat.

## 1. Engine and validation

**Reuse, not rewrite.** The engine is mf_ideal's (`tracks/mf_ideal/code`: `rpmfi.c` MODE 5, `vmodes.c`, `mfi_common.py`,
`release_mfi.py`, `vectors_mfi.py`, `chain_mfi.py`, `bench_run.py`). It is **imported read-only and unchanged**:
- `code/mfc_common.py` only redirects the module-global output paths to `results/engine/` and the scratch area here;
- the binaries are compiled from mf_ideal's sources into this track's scratch.
No fix was needed, so nothing was copied.

| check | result |
|---|---|
| bench_HFa as exported, raw R | S/Q(10⁻⁶) = **14.6933**, identical to mf_ideal's smoke test (14.693300) [C] |
| my history builder (`hist_build.py`, unblocked friction, switch at 0.6) vs the exported bench_HFa | K = 2.0481 (bench 2.0482); f_rad = 0.538 / 0.769 / 0.855 at t = 0.3 / 0.6 / 3 (bench identical to 3 digits); ε_H = 1 at t = 0.1926 (0.1927); H to 3.5×10⁻⁵ and R to 0.08 through the waterfall; chain S/Q(10⁻⁶) = 14.63 (0.4%) [C] |
| HFb rebuilt the same way | K = 2.0471 (bench 2.0473), S/Q = 14.54 (bench run 14.58) |
| every release | Gauss ≤ 2×10⁻¹⁰, hidden work ≤ 5×10⁻¹¹, N_X frozen to ≤ 0.2% |
| pre-slice exit factor f = [h_r F_s/F(exit)]² from the benchmark's own H(N), F(N) at 1–4 e-folds | **1 + 2.0×10⁻⁸** (ε_H ≈ 8×10⁻⁸ before the slice), so the vectors use the exact static-slice start [C, bench_run] |

## 2. The histories

### 2.1 As exported

mf_hybrid's bench_HFa / bench_HFb (§5 of its NOTES): homogeneous two-field run with unblocked axial decay friction to
t = 0.6 H_i⁻¹, then its fluid. K = 2.048 / 2.047. Run with raw R (resolved at Δt = 4×10⁻⁴) and with the cycle-mean R
(`--R trace`).

### 2.2 Skeptic corrections applied (`code/hist_build.py`)

The mf_hybrid skeptic's concrete corrections, all applied:
- **Pauli blocking (N1).** Stage 1, the homogeneous σ, χ run from the slice to t_sw = 0.45 H_i⁻¹ (the skeptic's
  fragmentation time, ratio < 2 at 0.46), uses the skeptic's Pauli-limited friction. Stage 2 is a copy of the skeptic's
  `sk_pauli.run` (production tables of mf_hybrid, F_sat(A) of the skeptic's exact multi-period modes, all read-only).
  - Changed in the copy: the per-species N energy, number and band-blocking energy are carried over from stage 1;
    M_P/H_i and the seesaw rate scale are parameters; the late tail uses an adaptive step.
  - The skeptic's Γ_th brackets (0, 1, 10 H_i) are run.
- **R(t).** Raw R through the coherent stage, then the cycle mean ρ_m/M_P² (w = 0 condensates) after fragmentation; the
  2m_χ and 2m_σ amplitudes are exported as R_osc. Switching at 0.6 instead of 0.45 moves S/Q(10⁻⁶) by +5% (53.9 vs 51.4).
- **f_s (N10).** Full spectator energy (kinetic plus potential) in every row.
- **m_s(t).** As mf_hybrid's bench.py: HFa g_s²χ² in stage 1 and g_s²(v² + ρ_χ/m_χ²) after; HFb 9 − (3/4)R.

### 2.3 Tails blocked: the skeptic's bracket does not converge

With the perturbative tails always blocked, the condensate can only empty through the band. But F_sat(A) → 0 as the
amplitude decays, so K grows without limit (matter forever; the skeptic's "≈ 22, not converged" was a snapshot) [C].
Physically, perturbative products at k = (m/2)β lie inside the Pauli band only while μ5,max > k. **"Tails blocked while
in the band"** is this track's refinement [A]. It converges, at K = 9.76.

### 2.4 Post-inflation histories at the design H_i = 1.30×10¹¹ GeV [C on the [E] rates]

| history (HFa; HFb within 0.1% in K) | K | f_rad at t = 0.3 / 0.6 / 3 / 100 | ρ_SM ≥ 10⁻³ at τ | ρ_SM ≥ 50% at τ |
|---|---|---|---|---|
| exported (unblocked) | 2.048 | 0.54 / 0.77 / 0.86 / 0.81 | 0.27 | 7.6 |
| Pauli, Γ_th = 0, tails free (**P0**) | 7.17 | 0.17 / 0.17 / 0.15 / 0.31 | 0.36 | ≈ 2×10² |
| Pauli, Γ_th = 1 H_i | 6.57 | 0.18 / 0.20 / 0.21 / 0.32 | 0.35 | ≈ 2×10² |
| Pauli, Γ_th = 10 H_i | 4.27 | 0.23 / 0.38 / 0.49 / 0.44 | 0.32 | ≈ 1.6×10² |
| Pauli, tails blocked in band | 9.76 | 0.17 / 0.16 / 0.10 / 0.08 | 0.36 | ≈ 3.6×10² |

These reproduce the skeptic's fluid (K = 7.05 / 6.41 / 4.10) to 2–4%. The difference is the coherent stage to 0.45
instead of a fluid from 0.19. The Higgs-plasma timing is mf_higgs's question; it is quoted only.

## 3. The chain on each history [C]

The release is mf_ideal's MODE 5 on the table: s evolves with the bench's m_s²(t), Γ_s = 1.3×10⁻⁶ H_i, the preparation
is the held point at μ_eff = 10i, R = 12. Then vectors on the history's own background, dq_heat's ledger, the closed form
S/Q(W), and the inherited p(S/Q) map.

| history (design H_i) | N_X | ε | E/E_c | N_YL | N_YH | N_σ | N_A(10⁻⁶) | R_vec / R_Y / R_casc (10⁻⁶) | S/Q 10⁻⁶ / 3×10⁻⁶ / ∞ | M (10⁻⁶) |
|---|---|---|---|---|---|---|---|---|---|---|
| HFa exported, raw R | −5404 | 0.101 | 21.2 | 11.07 | 0.20 | 0.021 | 2.03 | 4.59 / 9.42 / 0.55 | **14.69 / 12.25 / 10.86** | 1.58 |
| HFa exported, cycle-mean R | −5018 | 0.098 | 22.2 | 11.03 | 0.19 | 0.021 | 2.03 | 4.83 / 9.86 / 0.56 | 15.74 / 13.11 / 11.62 | 1.75 |
| HFb exported, raw R | −5412 | 0.102 | 21.1 | 10.98 | 0.19 | 0.021 | 2.03 | 4.59 / 9.34 / 0.52 | **14.58 / 12.15 / 10.77** | 1.57 |
| HFa P0 | −4069 | 0.082 | 28.3 | 6.82 | 0.11 | 0.010 | 10.4 | 48.2 / 11.8 / 0.62 | **51.4 / 27.1 / 15.8** | 5.72 |
| HFa P0, cycle-mean R | −3942 | 0.080 | 28.9 | 7.73 | 0.13 | 0.010 | 10.4 | 49.3 / 13.7 / 0.70 | 54.1 / 29.2 / 17.8 | 6.01 |
| HFa P1 | −4123 | 0.082 | 27.9 | 6.89 | 0.11 | 0.011 | 9.22 | 40.8 / 11.4 / 0.60 | 45.9 / 25.0 / 15.3 | 5.14 |
| HFa P10 | −4431 | 0.087 | 26.0 | 7.54 | 0.12 | 0.012 | 5.12 | 17.8 / 9.82 / 0.51 | 27.4 / 17.7 / 13.0 | 3.09 |
| HFa tails blocked in band | −4060 | 0.081 | 28.3 | 6.82 | 0.11 | 0.010 | 14.9 | 80.4 / 13.8 / 0.73 | 71.7 / 34.0 / 17.6 | 7.95 |
| HFb P0 / P10 / band | −4129 / −4476 / −4120 | 0.083–0.088 | 25.6–27.8 | 7.2–7.8 | 0.08–0.10 | 0.010–0.012 | 10.4 / 5.12 / 14.9 | — | 51.0 / 27.3 / 71.0 at 10⁻⁶ | 5.68 / 3.08 / 7.87 |

**Readings.**
- **No e_crit anywhere.** The Y_L deposit alone (R_Y = 9–15) keeps S/Q(e → ∞) ≥ 10.8. The bulk threshold S/Q = 3.47
  corresponds to W = 49.1, i.e. R_total < 1.3 at E/E_c = 21 and < 0.75 at E/E_c = 28.
- **Why N_YL is large.** The spectator mass turns on over ≈ 0.2 H_i⁻¹ (HFa: with the waterfall; HFb: with the drop of R)
  at an abrupt end. This is mf_ideal's fast-turn-on corner (its §4c: N_YL = 2–19 for t_on ≤ 0.1). The ξ_s mechanism
  (HFb) is no better than the waterfall-mass mechanism (HFa).
- **Pauli blocking cuts N_YL** (11 → 7) but lowers ε (0.10 → 0.08) and raises the vector number ∝ K^1.3. The net effect
  is S/Q ×3.5 (K = 7.2).
- **G7 check [C]:** N_A(10⁻⁶) = 2.03 / 5.12 / 9.22 / 10.4 / 14.9 at K = 2.05 / 4.27 / 6.57 / 7.17 / 9.76, so N_A ∝ K^1.28.
  R_vec (own E_X) ∝ K^1.83. "R_vec ∝ K²" holds approximately; there is no "×4 cut", because K rises.
- **Relic.** Free X X̄ with p = 0.28–0.29 X̄ per charge (exported) and 1.0–3.9 (Pauli) at 10⁻⁶. It is mostly symmetric
  on the Pauli histories, as on inf_chain's single-field benchmarks. Q-matter is never reached, so qm_micro's ⟨Q⟩ and M
  (and its skeptic's ×0.54–0.70 correction) are not invoked.

## 4. Heavy remnants and gravitational production at the abrupt end

**Homogeneous remnants [C].** N_YH = 0.19–0.20 (exported) and 0.05–0.12 (Pauli), against M1r's 0.0104. N_σ = 0.006–0.021
(M1r: 2×10⁻⁵). The cascade heat R_casc = 0.5–0.7 at 10⁻⁶ is in the ledger. N_YL is in §3.

**Finite-k gravitational production [C, `code/grav.py`].**
- **Method.** Free fields with mass m and coupling ξ on the history's own a(t) and raw R(t), mode by mode:
  u″ + [k² + a²(m² + (ξ − 1/6)R)]u = 0.
  - The modes start adiabatic in the de Sitter stage and are read out after the switch to the cycle-mean R.
  - A pure de Sitter control with identical initial data is subtracted mode by mode (start ambiguity plus integrator).
    The Wronskian holds to 5×10⁻⁷.
  - The post-switch tail is the perturbative rate (per real d.o.f.) ρ̇ = A²ω/64π, with A = |ξ − 1/6|R_osc and ω = 2m_χ
    or 2m_σ [E]. The conformal coefficient is (ξ − 1/6), cf. Ema et al. arXiv:1502.02475 [R].
- **Fields.** X-like: ξ = −5/6, m² = 20.2 H_i², so that m² + ξR = 10.2 in inflation, as the prepared state's lightest
  mode; bracket m² = 15–30. Y_L (22.25), Y_H (36.1) and σ_Σ (40 H_i): minimal.

| history | X: E/E_X (modes + tail) | N(Y_L) | N(Y_H) | N(σ_Σ) | ΔR in the ledger |
|---|---|---|---|---|---|
| HFa exported | 2.8×10⁻⁴ (38 + 4.7 H_i⁴; bracket 2.7–2.9×10⁻⁴) | 1.6×10⁻⁴ | 1.6×10⁻⁴ | 8×10⁻⁵ | 0.0012 |
| HFa P0 (design) | 4.2×10⁻³ (91 + 484) | 0.013 | 0.012 | 0.006 | ≈ 0.2 |
| HFa P0 (H_i = 2.40×10¹¹) | 5.6×10⁻³ (89 + 624) | 0.015 | 0.015 | 0.008 | 0.22 |

- The R step (12 → −4.5 → mean) produces X at k ≲ 5 H_i with n_k ≈ 0.04–0.09. Heavy Y and σ are adiabatic through the
  step: after the control subtraction, ≤ 6×10⁻⁴ quanta.
- The Pauli tail is larger because χ keeps ~80% of the energy for t ≲ 10².
- **Conclusion:** negligible on every history (≤ 0.4% of S/Q). The mf_hybrid skeptic's 1.1×10⁻³ of E_X for the unblocked
  remnant is consistent with the tail here once the remnant's actual fraction is used.

## 5. Hidden heat from s, σ and χ decays [E, `code/post.py`]

**The s leak.** The holder u = μ + G s/s_*, with G = 10i − μ, induces (dq_entropy's |u|² − μ² structure) the term
A(s/s_*)|Σ|², with A = 2μRe(G)c_L = −102.5 c_L H_i². Mixing with the Σ radial mode (m_σ = 40 H_i) and σ(∂π)²/h give
[D]

  Γ(s → A_L A_L) = (A/s_*)² m_s³ / [32π (m_σ² − m_s²)²],

independent of h; this is inf_chain's form with the canonical normalisation explicit. Then:
- s quanta per comoving volume: n_s a³ = 3M_P²K²(f_s/a)/m_s, from each history's late f_s/a;
- the branching BR = Γ_leak/Γ_s, with Γ_s = 1.3×10⁻⁶ H_i × (H_i/1.3×10¹¹ GeV) (mf_hybrid's light-s Yukawa; y_s² ∝ H/M_P
  under the rescaling rule) [model choice];
- each A_L deposits the vector deposit m_A⟨a_A⟩ of the history's own vector run.

| c_L (renormalisation condition) | 0 (vanish at M_Y, no threshold piece) | 1/16π² | 1/4π² (threshold size) | ln(M_P/M_Y)/4π² (vanish at M_P) |
|---|---|---|---|---|
| R_leak, exported / Pauli-sc | 0 | 0.53 / 0.68 | 8.5 / 11 | 1.5×10³ / 1.8×10³ |
| S/Q(10⁻⁶), exported / Pauli-sc | 14.7 / 61.6 | 15.1 / 62.1 | 20.4 / 68.3 | 462 / 686 |
| H_req/H_i at fixed history | 0.99 / 1.00 | 0.97 / 0.99 | 0.76 / 0.93 | 0.10 / 0.20 |

- **With c_L = 1/4π² inside the calibration** (unblocked HFa, three steps, `results/iterate_itHFa_trackleak.json`):
  H_i = 9.77 → 8.41 → 7.25×10¹⁰ GeV. The ratio goes 0.861 → 0.928 → 1.006, with K = 1.75, R_leak = 13.7,
  S/Q(10⁻⁶) = 22.4 (free X) and M = 2.58.
  - R_leak ∝ 1/Γ_s ∝ 1/H_i, so the leak strengthens as H_i falls.
  - The clock there is comfortably flat (N_max = 110 at κ = 10), and r = 8.6×10⁻⁸.
- **The s decay into the SM** (N_s → LH, a model choice; otherwise N_s is dark radiation, skeptic N8):
  - f_s at Γ_s t = 1 is 1.30–1.33% (unblocked) and 0.09–0.21% (Pauli; s stays a fixed fraction of the long matter era);
  - the entropy ratio is D = (1 + 0.886 f_s)^{3/4} = 1.0087 / 1.0009, so H_i,req rises by D^{2/3}: +0.6% / +0.06%.
  - There is no early s domination.
- **σ and χ decays.** Their only channels are the axial N's → LH; there are no σ- or χ-hidden portals (sequestered by
  choice) [A]. HFa's χ → ss (Γ = 9.6×10⁻¹⁴ H_i) adds n_s ~ 10³ against the condensate's ~10¹⁰ [E]. So their hidden heat
  is zero, up to the gravitational production of §4.
- **Direction.** Every hidden-heat item raises S/Q. None can open a Q-matter branch here.

## 6. The self-consistent H_i calibration (`code/iterate.py`, `code/cmb_rescale.py`)

**Rescaling rule (model choice).**
- Fixed in H_i units: m_χ, m_σ, the clock's tree c1 and m2 (so ε_c ∝ H² at fixed A_s), μ_R² (κ = 10), the singlet masses
  M_j and decay constants f_j, m_s = 3 and |g_h|s_* = 12.3.
- Fixed in M_P units: v, σ_c, s_*.
- Fixed in classical units: c/(λv³).
- Hence λ ∝ H², the axial amplitude A = Σ/2f ∝ M_P/H_i, and the seesaw N widths ∝ H_i (m̃ = 0.05 eV).
- The whole post-slice history is rebuilt and rerun at each step. The update is a secant in ln H_i.

| case | iteration: H_i (10¹¹ GeV) → H_req/H_i | converged H_i (GeV) | K | S/Q(10⁻⁶) | M | N_* | CMB capacity N_max (κ = 10) | r |
|---|---|---|---|---|---|---|---|---|
| HFa unblocked (e = 10⁻⁶) | 1.300 → 0.988; 1.284 → 0.995 | 1.28×10¹¹ | 2.04 | 14.6 | 1.57 | 52.64 | 75.0 OK | 2.68×10⁻⁷ |
| HFb unblocked | 1.278 → 1.001 | 1.28×10¹¹ | 2.03 | 14.5 | 1.56 | 52.64 | 84.6 OK | 2.66×10⁻⁷ |
| HFa Pauli P0 | 1.300 → 1.771; 2.302 → 1.039; 2.400 → 0.999 | **2.40×10¹¹** | 7.84 | 61.5 | 6.82 | 52.28 | **48.6 FAILS** | 9.37×10⁻⁷ |
| HFb Pauli P0 | 2.400 → 0.998 | 2.40×10¹¹ | 7.86 | 60.0 | 6.65 | 52.28 | 53.7 OK (n_s 0.9743) | 9.37×10⁻⁷ |
| HFa Pauli Γ_th = 10 | 1.300 → 1.503; 1.954 → 1.113; 2.260 → 0.982; 2.214 → 1.000 | 2.21×10¹¹ | 5.81 | 40.6 | 4.57 | 52.39 | 51.2 FAILS | 7.98×10⁻⁷ |
| HFa tails blocked in band | 1.300 → 1.938; 2.519 → 0.996 | 2.51×10¹¹ | 9.64 | 78.9 | 8.76 | 52.20 | 47.0 FAILS | 1.03×10⁻⁶ |
| HFa unblocked, calibrated at 3×10⁻⁶ | 1.453 → 1.074; 1.561 → 1.045; 1.752 → 1.003 | 1.75×10¹¹ | 2.37 | 12.6 (at 3×10⁻⁶) | 1.33 | 52.72 | 60.2 OK | 5.0×10⁻⁷ |
| HFa Pauli, calibrated at 3×10⁻⁶ | 3.825 → 1.025; 3.919 → 1.001 | 3.92×10¹¹ | 7.75 | 30.2 (at 3×10⁻⁶) | 3.38 | 52.53 | 38.4; **≤ 50.2 for every κ** | 2.5×10⁻⁶ |

- **Why the Pauli histories need a larger H_i.** The yield ∝ |N_X| K^{−3/2} falls, and the relic is more symmetric
  (M ≈ 7). The response is weak: d ln H_req/d ln H_i ≈ 0.03–0.07, so the iteration contracts fast.
- **K rises with H_i on the unblocked history** (2.03 → 2.37 from 1.28 to 1.75×10¹¹ GeV), because A ∝ 1/H_i weakens the
  non-perturbative transfer.
- **CMB [C on the inherited first-order capacity definition]** (mf_hybrid skeptic's `sk_structure`, read-only; regression
  74.30 at the design H_i, skeptic 74.3).
  - X = g²σ_c ∝ H_i rises to 0.49–0.56 at the Pauli fixed points.
  - At κ = 10 the e-fold capacity then falls below N_* for HFa (n_s = 0.9649).
  - With the skeptic's regularity cut (V′ ≥ 0.1 V′_*), κ = 15–50 gives N_max = 61–102 at 2.2–2.5×10¹¹ GeV
    (`results/cmb_kappa_scan.json`). So the self-consistent HFa needs a different tree quartic (renormalisation condition).
  - At 3.92×10¹¹ GeV no κ works for (250, 600) H_i.
  - n_s and A_s are retuned to target by the capacity construction. The exact Mukhanov–Sasaki spectrum at the new scales
    was not recomputed. ΔN_* = −0.36 shifts n_s by ≈ α_s ΔN ≈ 4×10⁻⁵ at fixed shape [E].
- **Retrospective comparison.** r = 0.9–2.5×10⁻⁶ is far below BICEP/Keck's 0.036 (arXiv:2110.00483 [R]). ACT DR6
  (arXiv:2503.14452 [R], n_s = 0.974 ± 0.003) pulls are unchanged from mf_hybrid's (HFa −2.8σ, HFb 0.0σ), up to the
  unrecomputed MS shift.

## 7. Isocurvature (`code/iso.py`) [C, with R for the bound]

- **Method.** D = d ln Ω_DM/d ln s_* at fixed g_s and H_i, from c0 = 0.99 / 1.01 releases (preparation re-solved, same
  table) with their own vectors and chain; Ω ∝ M|N_X|K^{−3/2} (mf_ideal's recipe). The bound is Planck 2018's
  β_iso < 0.038 (uncorrelated CDI, n = 1, k = 0.05 Mpc⁻¹; arXiv:1807.06211 [R], value as checked by mf_ideal and the
  mf_hybrid skeptic). β_iso = (H_i/2πs_*)²D²/A_s, with s_* = 3×10⁻³ M_P = 7.3×10¹⁵ GeV.

| converged history | D (10⁻⁶) | charge only | β_iso (10⁻⁶) | s_* needed (10⁻⁶) | D (3×10⁻⁶) | β_iso (3×10⁻⁶) |
|---|---|---|---|---|---|---|
| HFa unblocked (1.28×10¹¹) | 1.36 | 4.54 | 6.9×10⁻³ | > 3.1×10¹⁵ GeV | 3.03 | 0.033 (s_* > 6.8×10¹⁵) |
| HFb unblocked | 1.38 | 4.54 | 7.0×10⁻³ | > 3.1×10¹⁵ | 3.03 | 0.033 |
| HFa Pauli (2.40×10¹¹) | −0.07 | 4.54 | 6×10⁻⁵ | > 2.8×10¹⁴ | 0.46 | 2.7×10⁻³ |
| HFb Pauli | −0.01 | 4.59 | 2×10⁻⁶ | > 5×10¹³ | 0.59 | 4.5×10⁻³ |

- **The skeptic's N4 worry does not materialise on these histories.** With the charge-only response (p = 4.5) β would be
  0.07–0.21. But a larger s_* gives a larger kick, hence a lower S/Q and fewer pairs: M falls and cancels most of it.
- The cancellation is a property of the mostly symmetric relic. In an asymmetric relic the full 4.5 would apply.
- The unblocked history at 3×10⁻⁶ sits within 7% of the bound in s_*.

## 8. Against mf_ideal's map

mf_ideal's primary map at m_s = 3, interpolated (extrapolated linearly beyond K = 3, flagged):

| history | K | benchmark S/Q(10⁻⁶) | map, t_on = 0.1 | map, t_on = 1 |
|---|---|---|---|---|
| exported / unblocked sc | 2.03–2.05 | 14.5–14.7 | 29.6–29.8 | 18.2–18.4 |
| Pauli Γ_th = 10, sc | 5.81 | 40.6 | 77.7 (extrap.) | 67.0 (extrap.) |
| Pauli P0, sc | 7.84 | 61.5 / 60.0 | 103.6 (extrap.) | 93.3 (extrap.) |
| tails blocked in band, sc | 9.64 | 78.9 | 126.6 (extrap.) | 116.5 (extrap.) |

**Readings.**
- The benchmarks do better than the map at the same K: S/Q is ×0.5–0.8. The release is more efficient (ε = 0.10 against
  0.075 for the map's K = 2, t_on = 0.1 corner), because the hybrid's R falls over the 0.19 H_i⁻¹ waterfall instead of
  stepping from 12 to ρ_m.
- But N_YL is higher (11 against 6.4), the fast-turn-on mark.
- Both agree that fast turn-on plus K ≥ 2 leaves no e_crit.
- mf_ideal's best corners (t_on ≈ 0.5–1 with K ≈ 1) are not what HFa or HFb provide: t_on ≈ 0.2 and K = 2–10.

## 9. The brief's estimates, as far as this track tests them

| | verdict here |
|---|---|
| G5 | energy fraction after the turn-on 1.35×10⁻⁵ = 1.5(s_*/M_P)² **verified** on the benchmark [C]. The fraction at decay is 1.3% (unblocked) or 0.1–0.2% (Pauli), consistent with the skeptic's 1.3% and not mf_hybrid's 2.4%. Isocurvature: s_* ≳ 2.27×10¹⁵ GeV × \|D\| holds; D on the benchmark is 1.4 or ~0, not the charge-only 4.5 [C]. The ξ_s mechanism (HFb) gives the same fast turn-on as the waterfall mechanism (HFa) [C] |
| G7 | R_vec ∝ K^1.83 (own E_X), N_A ∝ K^1.28 [C]: approximately K². The fast benchmarks do not reach K → 1–1.3. Their branch is set by Y_L, not by the vectors [C] |
| G1–G4, G6 | not tested here (mf_gemini, mf_hybrid, mf_higgs) |

## 10. Limits

- **Inherited model.**
  - Homogeneous classical release; the single-Higgs-direction Proca vectors; dq_heat's ledger (decay at rest, full
    deposit).
  - The closed form S/Q(W) and the p(S/Q) recipe at S/Q = 15–80 (between wt_critic's gas map and inf_chain's case-L
    runs) [Q].
- **Pauli-blocked histories [E].** The rates and the blocking model are the skeptic's (harmonic, constant-amplitude
  saturation, in-band energy, Γ_th unknown). The "tails blocked while in band" rule is this track's.
- **Fragmentation.** The skeptic found that the waterfall fragments at t ≈ 0.45 and that walls are not excluded. No
  lattice was run here. After 0.45 the release sees only the cycle-mean R; HFa's m_s = g_s|χ| would be inhomogeneous at
  the ~5% level [E], which the homogeneous chain cannot represent.
- **Rescaling rule.** It is one choice, and others (e.g. a fixed EFT margin μ5,max/4πf) would give other K(H_i). The CMB
  side uses the first-order capacity, not exact MS spectra.
- **Hidden heat.** The s leak, Γ_s scaling and N_s completion are [E] brackets and model choices. The leak-included
  iteration was done only for the unblocked HFa at c_L = 1/4π².
- **Not addressed:** Higgs survival (mf_higgs); the skeptic's N5 (Σm_ν), N6 (anomalous couplings) and N8 (N_s loop
  tuning); finite-k Y driven by the holder modulation (adiabatic: m_s = 3 ≪ m_L); X ↔ Σ dephasing at fragmentation.

## 11. Code and results

| file | content |
|---|---|
| `code/mfc_common.py` | read-only import of mf_ideal's engine and drivers, output redirection |
| `code/run_bench.py` | whole chain on one bench-format file (mf_ideal's bench_run.py, unchanged) + isocurvature releases |
| `code/hist_build.py` | skeptic-corrected histories (stage 1 homogeneous with Pauli-limited friction, stage 2 Pauli fluid), rescaling rule |
| `code/iterate.py` | self-consistent H_i iteration |
| `code/cmb_rescale.py` | e-fold capacity (skeptic's sk_structure, read-only), N_* and r at the converged scales |
| `code/grav.py` | finite-k gravitational production with de Sitter control subtraction, plus the perturbative tail |
| `code/post.py` | final ledger: hidden heat brackets, s-decay entropy, gravitational add-ons, isocurvature β, map comparison |
| `code/iso.py`, `code/table.py`, `code/summary.py` | isocurvature driver, run table, summary |
| `results/summary.json` | headline numbers, every run, iterations, CMB checks, gravitational outputs |
| `results/final_ledger.json` | the eight final cases with every ledger item |
| `results/iterate_*.json`, `cmb_rescale.json`, `cmb_kappa_scan.json`, `grav/` | iterations, CMB, gravitational production |
| `results/engine/{release,vectors,chain,bench}/` | mf_ideal-format outputs of every run |
| `results/histories/` | the converged corrected histories (HFa and HFb Pauli P0 at 2.40×10¹¹ GeV, HFa unblocked at 1.28×10¹¹ GeV), bench format |

Reproduce (one or two cores, ≈ 45 min in total):

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_chain/code; export OMP_NUM_THREADS=1
python3 -B run_bench.py ../../mf_hybrid/results/bench_HFa.json HFa_raw          # likewise HFb_raw, --R trace
python3 -B hist_build.py HFa HFa_P0 --Gth 0          # P1, P10, --tails band; HFb; --fric track --tsw 0.6 (validation)
python3 -B run_bench.py $SP/wt2/mf_chain/hist/bench_HFa_P0.json HFa_P0
python3 -B iterate.py HFa itHFa_P0 --fric pauli --Gth 0   # itHFa_track (--fric track --tsw 0.6), itHFa_P10, itHFa_PBb, HFb, --e 3e-6
python3 -B iso.py itHFa_track_it1 itHFa_P0_it2 itHFb_track_it0 itHFb_P0_it0
python3 -B grav.py ../../mf_hybrid/results/bench_HFa.json HFa_raw               # and on the corrected histories
python3 -B cmb_rescale.py 2.3997e11:7.840:0.9649 ...; python3 -B post.py; python3 -B summary.py
```

## Literature (arXiv numbers checked by web search in this session)

- Planck 2018 X, arXiv:1807.06211: isocurvature bounds (β_iso < 0.038 as quoted by mf_ideal and the mf_hybrid skeptic).
- ACT DR6, arXiv:2503.14452: n_s = 0.974 ± 0.003.
- BICEP/Keck 2021, arXiv:2110.00483: r_0.05 < 0.036.
- Graham–Mardon–Rajendran, arXiv:1504.02102: longitudinal vectors from inflation.
- Greene–Kofman, hep-ph/9807339: fermionic preheating and Pauli blocking.
- Ema–Jinno–Mukaida–Nakayama, arXiv:1502.02475: gravitational couplings of an oscillating inflaton to matter.
- Our tracks, read-only:
  - mf_ideal (engine, drivers, map);
  - mf_hybrid (bench files, axial tables, mfh_common) and its skeptic (sk_pauli, sk_structure, sk_common, F_sat tables);
  - inf_chain (leak formula), dq_heat skeptic (dq_model_copy), wt_critic and inf_chain (p map), qm_micro and its skeptic
    (not needed: no Q-matter).
