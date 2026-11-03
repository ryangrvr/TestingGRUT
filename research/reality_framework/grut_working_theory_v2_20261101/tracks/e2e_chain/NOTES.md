# Track e2e_chain: the final self-consistent chain of M2

3 November 2026 · Claude, end-to-end pass (BRIEF_E2E.md of this directory; GitHub issue #13 context).

Same-session AI work on one supplied model. **Every history, splice and hand-over rule, closure, holder amplitude and
phase, rescaling rule, ledger entry, renormalisation condition and criterion is a model choice; none follows from
closure.** Nothing here is a validation, an exclusion of the real world, a new theory or a GRUT prediction. Data
comparisons are retrospective. Labels: [D] derived, [C] computed here, [Q] quoted from a reviewed track, [E] estimate,
[A] argued, [R] literature.

**Inputs read:** BRIEF_E2E.md (this track), WORKING_THEORY.md §0, §2, §4b, §7–§9, mf_critic NOTES, BRIEF_MULTIFIELD
(relevant parts), mf_chain NOTES (leak, iterate, post), mf_chain skeptic (phase scan), mf_hybrid NOTES (spectator
masses), e2e_lattice / e2e_kinetic / e2e_relic digests, results and interfaces. No history skeptic was run (budget);
the splice checks below are this track's own.

**Reuse (all read-only; outputs redirected to this track's scratch `$SP = .../scratchpad/wt2/e2e_chain`):**
mf_ideal's engine and drivers (rpmfi.c MODE 5, vmodes.c, release_mfi, vectors_mfi, chain_mfi, bench_run), compiled into
scratch; e2e_kinetic's `ek_run.run_history` / `export_fluid` and `ek_scan.occ_from_moments` (in-process, `EK_SCRATCH` →
`$SP/kin`, its p-tables copied there); e2e_relic's `pmap.p_selfconsistent`, `relic.relic`, `relic.observables`;
mf_hybrid's `inflation.calibrate/ms`, `mfh_common.derived/Nstar_from_K` and its skeptic's `sk_structure.nefolds`
(capacity, as mf_chain cmb_rescale); the leak formula and grav numbers of mf_chain post.py / grav (re-typed / quoted).
What changed is stated in each file header. Python ran with `-B`; no binary or `__pycache__` is in the tree.

---------------------------------------------------------------------------------------------------------------------

## 0. Result in brief

1. **The final histories (§1).** Lattice coherent stage to t_sw = 3 H_i⁻¹, then the kinetic network restarted from
   the lattice's own moments (the kinetic track's preferred interface), energy-conserving after two fixes. K [C]:

   | history | HFa 1.3×10¹¹ | HFb 1.3×10¹¹ | HFa 8.4×10¹⁰ |
   |---|---|---|---|
   | GM (box-mean Pauli saturation) | **15.28** | 15.26 | 16.42 |
   | GM + fragment spread σ_lnk = 0.3 | 12.56 | — | — |
   | P0 (local per-site blocking) | **8.77** | — | 9.21 |

   T_RH (SM at half the energy) = 3.1–4.3×10¹² GeV (1.3×10¹¹) and 2.0×10¹² GeV (8.4×10¹⁰). The blocking closure is
   the dominant uncertainty: it moves K by ×1.7 and decides whether M2-Q exists.
2. **M2-X exists on the primary (GM) history [C at nodes, E between].** c0 = 1.55, δ = 0 (μ_eff = −2.38 + 15.81i at
   the slice), self-consistent H_i = 1.12×10¹¹ GeV, K = 15.7, S/Q(10⁻⁶) = 4.03 (free X, M = 1.0003 m_X, p = 1.6×10⁻⁴),
   S/Q(3×10⁻⁶) = 3.05 (Q-matter), **e_crit = 1.72×10⁻⁶**. r = 2.05×10⁻⁷, N_* = 51.55, capacity 82 ≥ N_*.
3. **M2-Q does not exist on GM.** Over c0 = 1.45–1.8 and δ = −8°…+8°, S/Q(10⁻⁶) never falls below 3.78 at the
   1.3×10¹¹ node (3.49 at 8.4×10¹⁰, where the point is not self-consistent; 3.49 with the fragment spread). The floor is
   set by the vector heat (R_vec ≈ 5–8, ∝ K^1.9) plus the Y_H cascade (R_casc ≈ 15–25, which grows with c0).
4. **M2-Q exists on the P0 closure.** c0 = 1.45, δ = 0, H_i = 9.89×10¹⁰ GeV, K = 9.05, S/Q(10⁻⁶) = 2.43,
   e_crit = 0.43×10⁻⁶, Q-matter with M = 0.802 m_X, ⟨Q⟩ = 7.5×10¹⁹, objects of 9.6 μg (0.44 M_Pl), R = 0.73 nm,
   0.55 crossings m⁻² yr⁻¹. r = 1.59×10⁻⁷, N_* = 51.76, capacity 89.
5. **Windows.** M2-X: free X at 10⁻⁶ with an e_crit for self-consistent c0 ≈ 1.45–1.60 (H_i 1.5→1.0×10¹¹), δ ≈ −2°…+6°;
   e_crit ≈ 2×10⁻⁶ only near c0 ≈ 1.55. M2-Q (P0 only): c0 ≈ 1.40 (marginal) – 1.50, δ ≈ 0…+4°.
6. **Isocurvature.** D = d ln Ω/d ln s_* = 5.77 (M2-X) / 7.87 (M2-Q), so s_* ≥ 1.13×10¹⁶ / 1.36×10¹⁶ GeV
   (β_iso < 0.038). At the bench's s_* = 7.3×10¹⁵ GeV, β = 0.087 / 0.12.
7. **The spectator spread does not break the release, but the release is fragile.**
   - The physical spread δs/s_* ≤ 6×10⁻⁴ shifts S/Q by 0.002 / 0.006.
   - A separate-universe ±1σ m_s² patch would break it (S/Q 4.0 → 18–20; 2.4 → 5.9–6.5). That regime does not apply,
     because the m_s² inhomogeneity sits at k ~ m_χ ≫ m_s [A].
   - A ±3% shift of the mean lattice-stage m_s² moves M2-X's e_crit from 1.9 to 9.6 / 2.3×10⁻⁶. At the lattice's
     resolution, M2-X's e_crit is good to a factor of a few.
8. **c0 is not natural under de Sitter equilibrium.**
   - Equilibrium for s's residual inflationary mass (m² = 4.5×10⁻⁸ – 3×10⁻⁵ H²) gives s_rms = 36–920 H (4×10¹² –
     1×10¹⁴ GeV), reached only after 10⁵–7×10⁷ e-folds.
   - The isocurvature floor needs s_* ≈ 10⁵ H, 110–3900 × s_rms.
   - With g_h fixed by the floor, P(window) ~ exp(−5×10³) or smaller. With g_h free, the best case is P ≈ 3–5%, but
     then s ~ s_rms and β_iso ≈ 1.
   - So s_* is an initial condition from before the observable e-folds. It is frozen there: the tadpole drift is
     10 H in 50 e-folds [Q mf_hybrid].

---------------------------------------------------------------------------------------------------------------------

## 1. The final history and its sanity checks (`code/splice.py`, `code/handover.py`, `results/splice_checks.json`)

**Construction [model choices].**
- t < 0: the mf_hybrid bench rows (exit factor only).
- 0 ≤ t ≤ 3: e2e_lattice box means. Columns: a, H, the cycle-mean R (the raw R is aliased, with three samples per
  2m_χ period), and m_s² as a running mean of the box mean over two χ² periods.
- At t_sw = 3 (the lattice end, inside both tracks' advised windows): a moments hand-over to e2e_kinetic.
  - Condensates: ρ_σ = f_σ ρ; ρ_χ = (f_χ + f_grad) ρ. The fragments go to the homogeneous χ; the fragment spread is
    bracketed by the σ = 0.3 run.
  - SM: the energy split hard/thermal in the kinetic run's own ratio.
  - N: the energy split over the four species in the kinetic run's own ratio.
- The kinetic dephased network runs from t_sw to 10⁸. HFb uses its own lattice (B_P0GM) and kinetic split. The 8.4×10¹⁰
  histories use the lattice H84_* runs and the design kinetic split [A].

**Checks [C].**
- **Energy at the hand-over.**
  - The interface as shipped lost 19% (GM) and 75% (P0) of the handed-over N energy. Two causes: the n̄ ≤ ½ cap in
    `occ_from_moments`, and the coherent-grid k_max when the coherent stage is re-entered.
  - Fixes: a degenerate sphere (n̄ = ½, k_max = (8π²E)^{1/4}) where needed (every species; n̄ needed 0.5–3.5), and the
    dephased network started at t_sw. The kinetic reference itself switches there; t_coh 3 → 10 moves its K by 0.2%.
  - After the fixes, H is continuous to 5×10⁻⁴. K was insensitive to the lost energy: 15.26 → 15.28, 8.59 → 8.77.
- **Splice method.**
  - A naive row splice (homogeneous kinetic rows after t_sw, a rescaled; `sHFa_*`) gives K = 17.35 / 14.30. That is
    14% above the moments hand-over, because the lattice GM state has less χ (0.53 against 0.72 of ρ) and more σ
    (0.22 against 0.06).
  - A χ-weighted estimate K ∝ (comoving χ energy)^{2/3} [E] gives 15.43, consistent with the hand-over.
  - The naive splices are kept as K-sensitivity runs only.
- **Interface.** At t_sw the GM lattice and the homogeneous kinetic run agree on f_rad (0.248 against 0.219) but not
  on the σ:χ split. P0 has f_rad = 0.60, so a row splice is invalid there, and only the hand-over is used.
- **Lattice energy budget.** |budget| ≤ 7.4×10⁻³ (GM) [Q e2e_lattice].
- **Friedmann consistency of the table.** d ln a/dt against H matches to 4% inside the lattice block (aliased
  sampling) and to better than 10⁻³ in the kinetic block. The release engine's own H-consistency check reads
  1.2–1.3×10⁻².
- **K convergence.** K drifts by 3.5×10⁻⁵ over the last decade in t.
- **Lattice convergence** [Q e2e_lattice]: 2–5% in fractions at the base spacing. Its consequence through m_s² is §4.

---------------------------------------------------------------------------------------------------------------------

## 2. The holder scan (`code/chain.py`, `code/scan.py`, `results/scan_*.json`)

**Chain per point.**
- The mf_ideal MODE 5 release runs on the history's table, with the spectator evolved on the table's m_s²(t). The holder
  at the slice is u₀ = μ + c0 e^{iδ}(10i − μ); δ rotates the holder coupling's phase.
- Then the vectors, the dq_heat ledger and the closed-form S/Q.
- **Corrected relic** (e2e_relic): free X → M = 1 + 2p with the self-consistent p; Q-matter → `relic()`.
- **Leak.** A = 2μ Re(G_eff) c_L with G_eff = c0e^{iδ}(10i − μ). c_L = 0 is the declared condition: portals vanish at
  M_Y with no threshold piece. The variant c_L = 1/16π² adds R_leak = 0.36–0.51 and moves S/Q by ≤ +0.05.
- **Gravitational production** at the abrupt end: R ≈ 0.03 [Q mf_chain grav, scaled to E_X].
- **s decay entropy:** H_req rises by 0.07% / 0.26%.
- H_req = chain_mfi.H_of(M, Y/Y_Q / D_s).

**Selected points at e = 10⁻⁶ (c_L = 0)** [C].

| history (K) | c0, δ | S/Q 10⁻⁶ / 3×10⁻⁶ | e_crit | H_req/H_i |
|---|---|---|---|---|
| GM 1.3e11 (15.28) | 1.00 | 83.2 / 37.1 | — | 2.69 |
| | 1.45 | 4.02 / 2.87 | 1.51e-6 | 1.169 |
| | 1.50 | 3.86 / 2.82 | 1.38e-6 | 0.980 |
| | 1.50, +4° | **3.78** / 2.56 | 1.23e-6 | 1.068 |
| | 1.55 | 4.11 / 3.13 | 1.88e-6 | 0.854 |
| | 1.60 | 4.05 / 3.09 | 1.79e-6 | 0.761 |
| GM 8.4e10 (16.42) | 1.50 | 3.49 / 2.49 | 1.01e-6 | 1.561 |
| | 1.55 | 3.86 / 2.90 | 1.44e-6 | 1.368 |
| GM+σ0.3 1.3e11 (12.56) | 1.50, +4° | 3.49 / 2.37 | 1.01e-6 | 0.878 |
| P0 1.3e11 (8.77) | 1.35 | 6.42 / 5.65 | — | 0.981 |
| | 1.40 | **3.29 Q** / 2.58 | 0.84e-6 | 0.890 |
| | 1.45 | **2.45 Q** / 1.73 | 0.45e-6 | 0.735 |
| P0 8.4e10 (9.21) | 1.40 | 3.81 / 3.16 | 1.59e-6 | 1.269 |
| | 1.45 | **2.42 Q** / 1.76 | 0.43e-6 | 1.201 |
| HFb GM 1.3e11 (15.26) | 1.50 / 1.55 | 3.86 / 4.11 | 1.39 / 1.89e-6 | 0.980 / 0.854 |
| K = 17.35 splice (`sHFa_GM`) | 1.5 ± 8°, 1.6 + 8°, 1.7–1.8 | ≥ 4.0 everywhere | — or ≤ 2.3e-6 | |

- **HFb.** HFb reproduces HFa to < 1% in S/Q, e_crit and H_req.
- **What c0 trades.**
  - With c0, ε rises (0.58 → 0.79) and E/E_c falls (2.06 → 1.36).
  - N_YH and the cascade heat grow (R_casc 11 → 23).
  - N_YL has an interference minimum near arg μ_eff ≈ 108° (N_YL = 2.7–5.6).
  - S/Q therefore has a floor in (c0, δ) on every history. The floor sits below 3.47 only for K ≲ 12.

---------------------------------------------------------------------------------------------------------------------

## 3. Self-consistent points (`code/fixpoint.py` → `results/fixpoints.json`; `code/cmb.py` → `results/cmb.json`)

f = ln(H_req/H_i) is interpolated linearly in ln H between the two final histories (8.4×10¹⁰, 1.3×10¹¹). That is the
V0 rescaling of mf_chain's iterate.py rule, with K and the early stage recomputed at each node. The slope is
d ln(H_req/H_i)/d ln H ≈ −1.05 (GM) and −0.8 to −1.1 (P0). The root is the self-consistent H_i [C at nodes, E
between]. N_* is then re-solved with the root's K (mf_hybrid's Nstar_from_K). The clock is re-calibrated at the rescaled
V0 and the exact Mukhanov–Sasaki spectrum is computed. n_s is the retuned target, not an output.

| | **M2-X** (GM, HFa) | **M2-Q** (P0, HFa) |
|---|---|---|
| c0, δ; μ_eff(slice) | 1.55, 0; −2.38 + 15.81i (\|μ\| 15.99, 104.3°) | 1.45, 0; −1.66 + 14.85i (102.5°) |
| \|g_h\|s_* | 19.1 H_i | 17.8 H_i |
| H_i (self-consistent); V0^{1/4} | **1.123×10¹¹ GeV** (w = 0.66); 6.85×10¹⁴ GeV | **9.89×10¹⁰ GeV** (w = 0.37); 6.43×10¹⁴ GeV |
| K; T_RH | 15.7; ≈ 2.5–3×10¹² GeV [E between nodes] | 9.05; ≈ 2–3×10¹² GeV [E] |
| ε, E/E_c, N_YL, N_YH, N_σ | 0.75, 1.48, 35, 23, 4.0 | 0.61, 1.93, 6.4, 9.2, 1.5 |
| ledger at 10⁻⁶: R_vec, R_Y, R_casc, R_grav, R_leak | 5.6, 13.7, 19.6, 0.03, 0 (c_L = 0) | 3.2, 3.7, 8.2, 0.03, 0 |
| S/Q at 10⁻⁶ / 3×10⁻⁶ / ∞ | 4.03 / 3.05 / 1.86 | 2.43 / 1.75 / 0.91 |
| branch at 10⁻⁶ | free X | Q-matter (Q = 2×10⁵ objects allowed) |
| e_crit (bulk / Q = 2×10⁵) | **1.72×10⁻⁶** / 3.8×10⁻⁶ | **0.43×10⁻⁶** / 0.63×10⁻⁶ |
| relic (e2e_relic) | free X: p = 1.6×10⁻⁴, M = 1.0003 m_X; n = 3.2 cm⁻³; σ/m = 6.0×10⁻⁸ (XX), 7.1×10⁻⁷ (XX̄) cm²/g; no annihilation channel | Q-matter: M = 0.802 m_X; ⟨Q⟩ = 7.2–7.5×10¹⁹; f_v = 0.094 (perfect) / 0.138 (Beliaev); 9.6 μg = 0.44 M_Pl; R = 0.73 nm; 0.55 m⁻² yr⁻¹; impulse 3.6×10⁻²⁷ kg m/s on 1 mg at 1 mm |
| m_A at e = 10⁻⁶ | 2.7×10⁷ GeV (Higgsed; the relic is e-independent at fixed S/Q) | 2.4×10⁷ GeV |
| N_*; r; capacity | 51.55; **2.05×10⁻⁷**; N_max = 82 ✓ | 51.76; **1.59×10⁻⁷**; N_max = 89 ✓ |
| MS at the final scale | n_s 0.96490 (target), A_s 2.1000×10⁻⁹, **α_s = +9.5×10⁻⁵**, r 2.05×10⁻⁷ | n_s 0.96490, **α_s = +7.4×10⁻⁵**, r 1.59×10⁻⁷ |
| HFb at the M2-X scale | N_* 51.55, n_s 0.97430 (target), α_s +1.1×10⁻⁴, N_max 93 | not run |
| isocurvature D (charge-only) | 5.77 (5.77) | 7.87 (8.02) |
| s_* floor (β_iso < 0.038) | **1.13×10¹⁶ GeV** (4.65×10⁻³ M_P) | **1.36×10¹⁶ GeV** (5.6×10⁻³ M_P) |
| s energy at decay (s_* at floor); entropy ratio | 0.16%; 1.0010 (H_req +0.07%) | 0.60%; 1.0040 (+0.26%) |
| axial EFT μ5,max/4πf (mf_critic scaling (0.42–0.49)·1.3×10¹¹/H_i) [E] | 0.49–0.57 | 0.55–0.64 |

- The CMB comparisons are retrospective [R, arXiv numbers as verified by mf_critic / mf_hybrid]. n_s is put in (Planck
  2018 0.9649 ± 0.0042, arXiv:1807.06211; ACT DR6 0.9743 ± 0.0034, arXiv:2503.14452). r ≈ 2×10⁻⁷ is far below
  BICEP/Keck's 0.036 (arXiv:2110.00483). α_s ≈ +10⁻⁴ is unobservable.
- **Heavy remnants.** Y_L and Y_H decay deep in radiation (Γ_L = 1.0×10⁻⁷ H_i). Their heat is R_Y in the ledger, and
  nothing heavy survives [Q ledger]. The σ_Σ quanta are in R_casc.
- **Gravitational production.** The X heat is 0.6% of E_X at c0 = 1 [Q mf_chain], 0.03 of the ledger here. Y and Σ add
  1.5×10⁻² and 8×10⁻³ quanta: negligible.
- **The s decay.** N_s → LH is a model choice; otherwise N_s is dark radiation (skeptic N8). Γ_s = 1.3×10⁻⁶ H_i
  (H_i/1.3×10¹¹), so it decays at a ≈ 4.7–4.8×10³ in radiation.

---------------------------------------------------------------------------------------------------------------------

## 4. Does the lattice's spectator spread break the release? (`scan_*_sp±1`, `scan_*_ms*`) [C/A]

- **Physical spread.** δs/s_* = 6×10⁻⁴ → 2×10⁻⁴ [Q e2e_lattice]. The engine's m_s² enters only the s equation, and the
  holder only through s/s_*, so δs/s_* acts like δc0/c0. dS/Q/d ln c0 = +2.8 (M2-X) and −10.2 (M2-Q) give
  ΔS/Q = 0.002 / 0.006. **No break.**
- **Separate-universe bound** (not physical here). Patches with the lattice's ±1σ m_s² for the whole coherent stage
  change c(t) by up to 0.15 at t = 1–3:
  - M2-X: S/Q 4.11 → 17.9 / 20.2, N_YL 38 → 1090 / 785;
  - M2-Q: S/Q 2.42 → 5.9 / 6.5.

  The inhomogeneity is at k ~ m_χ ≫ m_s, and the box is sub-Hubble. Super-Hubble modulation of the waterfall comes
  from δσ, which is adiabatic (a time shift). So s sees the box mean [A]. Patches larger than 1/m_s with ±1σ m_s²
  would break the release.
- **Mean-profile sensitivity** (an uncertainty of the history, not a spread). Lattice-stage m_s² × 1.03 / × 0.97:
  - M2-X: S/Q 4.11 → 4.88 / 4.27; e_crit 1.88 → 9.6 / 2.3×10⁻⁶.
  - M2-Q: × 1.03 gives 2.42 → 2.59, still Q-matter.

  M2-X's e_crit is fragile at the lattice's few-% resolution. M2-Q's branch is not.

---------------------------------------------------------------------------------------------------------------------

## 5. Is c0 natural? (`code/points.py`, `naturalness` in `results/m2_points.json`) [D/E]

**What equilibrium gives.** A light s in de Sitter relaxes to the Starobinsky–Yokoyama distribution
(arXiv:astro-ph/9407016, verified this session): P(s) ∝ exp(−8π²V/3H⁴), so ⟨s²⟩ = 3H⁴/(8π²m²), after
N_relax ≈ 3H²/m² e-folds. The residual inflationary masses are [Q mf_hybrid]:
- the holder loop, 2.5×10⁻⁷ H²;
- HFa's waterfall-induced g_s²χ_min², ≤ 3×10⁻⁵ H²;
- HFb's (m_s0² + ξ_sR), which equals 4.5 ε_H H² ≈ 4.5×10⁻⁸ – 4.5×10⁻⁷ H².

| m²/H² | s_rms/H | s_rms (GeV) | N_relax | s_*,min/s_rms (M2-X / M2-Q) |
|---|---|---|---|---|
| 3×10⁻⁵ | 36 | 4×10¹² | 10⁵ | 2830 / 3870 |
| 4.5×10⁻⁷ | 290 | 3.3×10¹³ | 6.7×10⁶ | 347 / 473 |
| 2.5×10⁻⁷ | 390 | 4.4×10¹³ | 1.2×10⁷ | 259 / 353 |
| 4.5×10⁻⁸ | 920 | 1.0×10¹⁴ | 6.7×10⁷ | 110 / 150 |

**Reading.**
- **Equilibrium cannot reach s_* ~ 10¹⁶ GeV.** It falls short by ×110–3900. Reaching it by free diffusion
  (Δs = H√N/2π) would need about 4×10¹¹ e-folds, longer than every relaxation time above. So equilibrium caps s first.
- **With g_h fixed by the isocurvature floor,** c0 = c0,centre s/s_*,min. Then P(c0 ∈ window) ~ exp(−5×10³) to
  exp(−7×10⁶).
- **With g_h free,** the best case is when the equilibrium width matches the window (σ_c0 ≈ 1.5). P = 4.8% (M2-X,
  1.45–1.60) or 3.3% (M2-Q, 1.40–1.50). But then s ~ s_rms ~ 10² H and β_iso ≈ 0.93–1.0 (D = 1 already; worse at
  D = 6–8), far above Planck's bound (retrospective).
- **What sets s_*.** It is an initial condition fixed before the observable e-folds and frozen through them. The
  tadpole drift is 10 H in 50 e-folds, and the classical roll at m² ≤ 3×10⁻⁵ H² over 60 e-folds is < 10⁻³. Examples:
  the end of an earlier stage, such as a negative Hubble-induced mass in a pre-inflationary era, or inflation that
  lasts less than N_relax [A]. The phase δ (equivalently α) is not distributed by equilibrium at all: it is the
  supplied CP phase. **c0 and δ remain supplied numbers.**

---------------------------------------------------------------------------------------------------------------------

## 6. Parameter-card entries this track fixes (all supplied or derived as marked)

- **c0 (supplied).** M2-X 1.55 (window ≈ 1.45–1.60 self-consistent, e_crit ≈ 2×10⁻⁶ only near 1.55). M2-Q 1.45
  (window ≈ 1.40–1.50, P0 closure only).
- **δ (supplied).** 0 (window ≈ −2°…+6° for M2-X; 0…+4° for M2-Q).
- **e (supplied).** 10⁻⁶.
- **s_* (supplied).** ≥ 1.13×10¹⁶ (M2-X) / ≥ 1.36×10¹⁶ GeV (M2-Q); |g_h| = 19.1 / 17.8 H_i/s_*.
- **c_L (supplied).** 0.
- **Blocking closure (supplied, decisive).** GM for M2-X, P0 for M2-Q.
- **t_sw = 3 and the hand-over rule (supplied).**
- **H_i (fitted to Ω_DM through the V0 rescaling).** 1.12×10¹¹ / 9.9×10¹⁰ GeV.
- **Derived / computed.** K, S/Q, e_crit, M, ⟨Q⟩, r, α_s, N_*.

---------------------------------------------------------------------------------------------------------------------

## 7. Limitations [A]

- **No history skeptic.** The splice checks are self-checks.
- **The blocking closure is open.** It decides whether M2-Q exists (§0.3–4). Which closure is physical was not
  computed here or by the lattice.
- **The self-consistent H_i interpolates between two final histories.** Their early stages are separate lattice runs
  at 8.4×10¹⁰ and 1.3×10¹¹. The 8.4×10¹⁰ hand-overs use the design kinetic N split. The c0 = 1.45 GM root is
  extrapolated (w = 1.34). No history was built at the roots themselves.
- **Lattice-stage inputs.** The lattice-stage R is the cycle mean, which mf_chain's skeptic found the release robust
  to; the raw-R effect was not re-checked. The m_s² profile is the running box mean, and the release is sensitive to
  it at the few-% level (§4).
- **The fragment spread** was bracketed (σ_lnk = 0.3 → K 12.6) only at the 1.3×10¹¹ GM node, not at 8.4×10¹⁰.
- **Phase scans are coarse** (±4–8° steps). The phase window edges for M2-X come partly from the K = 17.35 row-splice
  history.
- **Isocurvature D** comes from one node per point (1.3×10¹¹ for M2-X, 8.4×10¹⁰ for M2-Q), not at the root.
- **Q-matter microphysics** (⟨Q⟩, f_v) is e2e_relic's transfer of qm_micro with the Y_Q exponent from S/Q = 3, not
  recomputed on this history.
- **Not covered here:** the Higgs (e2e_higgs) and the SM side (e2e_sm).
- **Not run (budget):** histories at 2.4×10¹¹ GeV, more phase points, D at the roots, a P0 fragment-spread run, a
  history skeptic.

## 8. Files

- `code/e2c_common.py` (redirects), `splice.py` (final histories), `handover.py` (lattice → kinetic moments),
  `chain.py` (release c0/δ, corrected relic, leak, grav), `scan.py` (driver), `fixpoint.py`, `cmb.py`, `points.py`.
- `results/m2_points.json` (every number of §3–§5), `fixpoints.json`, `cmb.json`, `splice_checks.json`,
  `handover_*.json` (kinetic summaries), `scan_*.json` (every chain run; `sHFa_*` = naive row splices kept for K
  sensitivity; `*_sp±1`, `*_ms*` = §4 checks).
- Scratch (not in tree): `$SP/hist/bench_*.json` (final histories), `$SP/kin/fluid_*.json` (kinetic exports),
  `$SP/engine/` (release/vector/chain files).

**Literature verified this session:** Starobinsky & Yokoyama, arXiv:astro-ph/9407016. Quoted through reviewed tracks
[R]: Planck 2018 arXiv:1807.06211; ACT DR6 arXiv:2503.14452; BICEP/Keck arXiv:2110.00483.
