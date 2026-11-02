# SKEPTIC of tracks mf_chain and mf_ideal (key mf_chain_skeptic)

2 November 2026 · Claude, adversarial pass in the multi-field inflaton pass (BRIEF_MULTIFIELD.md; GitHub issue #13 context).
- This is same-session AI work on one supplied model, not outside replication.
- Nothing here is a validation, an exclusion of the real world, a new theory or a GRUT prediction.
- **Every history, holder amplitude and phase, pair-fraction map, rescaling rule, ledger entry and criterion, in the tracks and
  here, is a model choice.** None follows from closure.
- Data comparisons (Planck 2018, BICEP/Keck 2021, ACT DR6) are retrospective.
- Git was read-only and I committed nothing. (The orchestrator's interim commit 45cfc2a picked up early versions of my files;
  that was not my action.)
- All runs were in the foreground. The Python driver plus mf_ideal's two vector-mode subprocesses is the largest load used.
- No lattice was run. No compiled file or `__pycache__` was left in the tree: Python ran with `-B`, and the binaries are in
  scratch (`$SP/wt2/mf_chain_skeptic/`).

Labels: [D] derived, [C] computed with my code (`skeptic/code/`, output in `skeptic/results/`), [E] estimate, [A] argued,
[Q] quoted from a reviewed track, [R] literature (arXiv numbers checked by web search in this session).

Units: H_i = 1 at the slice, a = 1 there, M_P = 2.435×10¹⁸ GeV (reduced), K = a²H in radiation.
"c0" is the holder amplitude in units of M1r's calibration. It sets μ_eff(slice) = μ + c0(10i − μ), so c0 = 1 is 10i and
|g_h|s_* = 12.30 c0 H_i, at the benchmark's phase α. Changing c0 is the same as rescaling s_* at fixed g_h.

---------------------------------------------------------------------------------------------------------------------

## 0. Verdict on one screen

**What stands.**
- The engine and the vector heat survive two independent re-implementations: release to ≤ 0.14% in N_X, vectors to ≤ 5×10⁻⁴
  in n_k.
- The exported and Pauli-blocked branch numbers at the benchmark holder (c0 = 1) stand. So do the leak, the s-decay
  entropy, the gravitational production and the CMB capacity.

**What is refuted.**
- **The plain answer.** mf_chain says the repaired hybrid "does not recover M1r's e-controlled branch point", with Y_L
  "forbidding" Q-matter. That holds only at the copied M1r holder calibration, μ_eff(slice) = 10i.
  - With the holder 36–40% stronger at the same phase, the Pauli-corrected HFa reaches its own self-consistent H_i.
  - At c0 = 1.40: H_i = 8.40×10¹⁰ GeV, H_req/H_i = 0.998. The relic is Q-matter at e = 10⁻⁶ (S/Q = 2.88), with
    e_crit = 0.60×10⁻⁶.
  - At c0 = 1.36 (H_req/H_i = 1.05) it gives M1r's own picture: free X at 10⁻⁶ (4.16), Q-matter at 3×10⁻⁶,
    e_crit = 2.2×10⁻⁶.
- **The relic, H_i and isocurvature on the unblocked histories.** These rest on an unsampled gap of the inherited
  pair-fraction map.
- **β_iso at e = 3×10⁻⁶.** It was evaluated at the wrong H_i.
- **mf_ideal's "pocket only" and "never for fast turn-on or K ≥ 1.6".** Both fail beyond c0 = 1.2.

**Verdict table, mf_chain**

| Track result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| C1 regression | **stands** | My own Python release solver (DOP853, own equations, splines, preparation and N_YL read-out; `sk_release.py`) on bench_HFa: late (75–100) means N_X −5306.5 vs engine −5299.2 (1.4×10⁻³), E_X 7×10⁻⁴, ε 8×10⁻⁴, N_YL 11.072 vs 11.074 (1.3×10⁻⁴). On the Pauli fixed-point history: N_X 5×10⁻⁴, N_YL 7×10⁻³. Gauss ≤ 2×10⁻¹¹. Full-chain regression in my scratch: S/Q(10⁻⁶) = 14.6933 (§1.1) | C |
| C2 exported branch (c0 = 1) | **branch stands; relic mass refuted** | Free X at all e is reproduced. M = 1.58 rests on the map gap (C6/C8): with the case-L network run inside the gap, M = 1.77 at 10⁻⁶ and 1.52 at 3×10⁻⁶ (not 1.58 / 1.30) | C |
| C3 Pauli branch (c0 = 1) | **stands** (K not rechecked) | My release reproduces itHFa_P0_it2 (N_YL 5.445 vs 5.483). The Pauli relics sit above the map gap, so M is unchanged | C |
| C4 "the Y_L remnant exceeds the whole Q-matter budget by itself" | **refuted as a general statement** | The budget is not fixed. R_total < W_thr/(E/E_c) − 1, and E/E_c falls steeply with the holder amplitude. On the same exported history at c0 = 1.30: N_YL = 25 and R_Y = 6.2, yet E/E_c = 3.14, the budget is R < 14.6, and S/Q(10⁻⁶) = 2.75 (Q-matter). At c0 = 1 the statement is reproduced (§2.1) | C |
| C5 tails-blocked bracket | not rechecked | — | — |
| C6 self-consistent H_i | **Pauli numbers stand; unblocked fixed point refuted** | With the gap-filled p map, the unblocked HFa history rebuilt at 1.12×10¹¹ (track's own hist_build, rescaling rule) needs 1.112×10¹¹: the **fixed point is ≈ 1.11×10¹¹ GeV, not 1.28×10¹¹**. At 3×10⁻⁶ the H_req at the 1.75×10¹¹ history falls to 1.58×10¹¹, so that fixed point also moves down [C step / E fixed point] | C |
| C7 CMB at fixed points | **stands** | e-fold capacity reproduced: X = 0.5335, N_max = 48.57 < N_* = 52.28 at 2.40×10¹¹ (κ = 10); 74.30 at the design | C |
| C8 isocurvature | **refuted (values)** | (i) The unblocked D comes from the map's log-linear bridge (slope 3.68). The network slope is 1.93, so D(10⁻⁶) = **2.49, not 1.36**, and β_iso = 0.022, not 6.9×10⁻³. (ii) β at 3×10⁻⁶ was evaluated at the 10⁻⁶ H_i. At the track's own 3×10⁻⁶ fixed point (1.752×10¹¹ GeV, its own map) D = 2.40 and **β_iso = 0.0385**, not 0.033: at or above Planck's 0.038, and s_* ≥ 7.36×10¹⁵ GeV is needed (bench 7.305). With the gap-filled map on that history, D = 2.68 and β = 0.047 | C |
| C9 gravitational production | **stands** (tail ×1.25–1.5) | My direct integral of 2A²ω/64π·a⁴ over the history's R_osc columns gives 776 H_i⁴ (P0 at 2.40×10¹¹; track 624) and 7.2 (exported; track 4.7). That is still < 1% of E_X | C |
| C10 vs mf_ideal's map | **stands** | mf_ideal's dust fluid at K = 7.17 gives S/Q(10⁻⁶) = 96.1 against the track's extrapolated 104 / 93 | C |
| E1 s → A_LA_L leak | **stands** (reproduced) | The formula was re-derived in the linear basis: the direct s ππ vertex plus σ exchange gives the amplitude A_s m_s²/(m_s² − m_σ²). R_leak = 0.534 / 8.54 / 1521 (exported) and 0.680 / 10.87 / 1762 (Pauli), identical to the track | D/C |
| E2 s-decay entropy | **stands** | f_s/a = 7.35×10⁻⁶; at Γ_s t = 1 (a = 1.78×10³), f_s = 1.31% | C |

**Verdict table, mf_ideal**

| Track result | Verdict | Recomputed / reason | Label |
|---|---|---|---|
| R1, R2 (c0 = 1) | **stand** | My independent solver reproduces the engine on two hybrid histories. K = 2, m_s = 3, t_on = 0.1 at c0 = 1 re-run: 29.15 (map 29.1) | C |
| R3 "e_crit never for t_on ≤ 0.1 or K ≥ 1.6" | **refuted beyond c0 = 1.2** | The same history (K = 2, m_s = 3, t_on = 0.1) at c0 = 1.4: ε = 0.487, S/Q(10⁻⁶) = **2.26 (Q-matter)**, e_crit = 2.5×10⁻⁷. Not universal: K = 3, m_s = 4, t_on = 0.1 gives 19.2 / 9.1 at c0 = 1.4 / 1.5, and K = 1, m_s = 3, t_on = 1 gives 70 at c0 = 1.4 (N_YL = 1280) | C |
| R4 "Q-matter at 10⁻⁶ only in a narrow pocket (K = 1, c0 = 1.03–1.05)" | **refuted** | A second, separate region exists at c0 ≈ 1.3–1.5: the K = 2 fast corner above, every hybrid history tried (§2), HFb | C |
| R5 release mechanism | stands | Hybrid, holder held: N_YL = 1.41 (R history alone) | C |
| R6 "N_A depends only on K (~K^1.45)" | **approximately stands** | At K = 7.17 the dust fluid gives N_A(10⁻⁶) = 9.34 against 10.40 on the Pauli hybrid (11%). The K^1.45 law extrapolated from K = 3 gives 12.4, so it steepens too fast | C |
| R8 isocurvature | stands at c0 = 1, with the caveat of R9 | — | — |
| R9 relic and H_i | **refuted in part** | 47 of the c0 = 1 histories have S/Q(10⁻⁶) inside the map gap (6.74–17.1). There M rises by 2–17% and H_req falls by 2–10% at fixed history (`sk_mfideal_gap_rows.json`) | C |
| R12 leak constraints | stands (same formula, reproduced in E1) | — | D |

---------------------------------------------------------------------------------------------------------------------

## 1. Independent recomputations

### 1.1 The release (`sk_release.py`) [C]

**What is independent.** None of rpmfi.c is used:
- equations re-derived from the M1r potential (BRIEF_DECISIVE.md, λ = 0.5);
- scipy DOP853 in place of the engine's step-capped DP5;
- bench background as two-piece not-a-knot splines (H, R, m_s²) joined at the bench's switch, with ln a the exact integral
  of the H spline;
- the preparation solved by my own least squares, with continuation in the holder amplitude;
- N_YL, N_YH read with my own forced-point subtraction (finite-difference forced velocity).

**Results.**

| history | quantity (late mean, τ 75–100) | mine | engine | rel. diff |
|---|---|---|---|---|
| bench_HFa (exported) | N_X / E_X / ε | −5306.5 / 1.5054×10⁵ / 0.09955 | −5299.2 / 1.5043×10⁵ / 0.09947 | 1.4×10⁻³ / 7×10⁻⁴ / 8×10⁻⁴ |
| | N_YL / N_YH / N_σ | 11.072 / 0.200 / 0.0209 | 11.074 / 0.2025 / 0.0207 | 1×10⁻⁴ / 1×10⁻² / 5×10⁻³ |
| HFa Pauli P0 at 2.40×10¹¹ | N_X / N_YL | −3705.3 / 5.445 | −3703.3 / 5.483 | 5×10⁻⁴ / 7×10⁻³ |

- **The large N_YL is real, not an engine artefact.** The c0 = 1.40 point of §2 agrees between the two codes too
  (ε 0.5336 vs 0.5339, N_YL 13.65 vs 13.74), and is identical at rtol 10⁻¹⁰ and 10⁻¹².
- **A sensitivity worth knowing [C].** My first version used one global spline across the bench's switch at t = 0.6, where
  the sampling jumps from Δt = 4×10⁻⁴ to 0.06. The spline overshoots there in R and m_s² for 0.06 H_i⁻¹, and that moves N_X
  by −9% and N_YL by −2.5%. The engine's finite-difference Hermite slopes are benign. Future exports should carry a dense
  overlap across the switch.

### 1.2 The vectors (`sk_vec.py`) [C]

**What is independent.**
- F(t) from my own release;
- e-folds as the independent variable, scipy DOP853;
- my own de Sitter Goldstone start and read-out.

Same action as vmodes.c.

**Result.** On the Pauli fixed-point history at e = 10⁻⁶:

| k | 0.01 | 0.063 | 0.158 | 0.398 | 1 |
|---|---|---|---|---|---|
| n_k mine / track | 1.00004 | 1.000005 | 1.00005 | 0.99966 | 0.99945 |

The Wronskian is 1 to 3×10⁻¹¹. The vector heat on the new histories stands.

### 1.3 The ledger, closed form and calibration (`sk_ledger.py`) [C]

- My own abundance arithmetic (ρ_r a⁴ = 3K²M_P²H², s from g = 106.75, ρ_DM/s = M m_X Y = 4.4×10⁻¹⁰ GeV) reproduces every
  H_req of the track to 10⁻⁴.
- The ledger sums and S/Q reproduce exactly.
- The leak (`sk_leak.py`) and the CMB capacity (via the mf_hybrid skeptic's sk_structure, read-only) also reproduce.

### 1.4 The gap of the pair-fraction map (`sk_pmap.py`) [C]

**What the track does.** chain_mfi's p(S/Q) joins wt_critic's gas map (S/Q ≤ 6.74) to inf_chain's case-L runs (≥ 17.1) by a
log-linear bridge.

**What I did.** I ran the same case-L network (fx_chain's m1_kin, read-only, gas rule, a_ref start, exactly inf_chain's
pmap_kin recipe) inside the gap.

| S/Q | 6.742 | 8 | 10 | 12 | 14 | 15 | 16 | 17.13 | 21.8 |
|---|---|---|---|---|---|---|---|---|---|
| p, track's bridge | 0.0166 | 0.031 | 0.071 | 0.138 | 0.244 | 0.315 | 0.399 | 0.513 | 0.762 |
| p, case-L network (mine) | **0.0308** | 0.071 | 0.154 | 0.247 | 0.345 | 0.395 | 0.446 | 0.503 (inf_chain 0.5034) | 0.744 (inf_chain 0.7439) |

- **The two recipes disagree by ×1.85 at their junction.** The bridge slope is 3.68; the network's is 1.93 at 14–16.
- **The unblocked hybrid relics (S/Q = 12.2–15.9) sit in the gap.** So do 47 of mf_ideal's c0 = 1 histories.

---------------------------------------------------------------------------------------------------------------------

## 2. The main refutation: the holder amplitude decides the branch on the hybrid histories [C]

**Starting point.** mf_hybrid copied M1r's calibration, |g_h|s_* = |10i − μ| = 12.30 H_i. mf_chain ran every history at that
value only. But the holder strength is a free model choice, exactly like g_I in M1r.

**The preparation stays healthy up to c0 = 2 [C, `sk_release.prep_hessian`].** The held stationary point:
- exists and is a minimum (lightest mass² 10.1 / 9.5 / 12.2 / 15.7 / 20.2 H_i² at c0 = 1 / 1.2 / 1.3 / 1.4 / 1.5, all above
  9/4);
- grows in amplitude: |X| = 91 → 157 at c0 = 1.4, |Y| = 6.9 → 24, with Y partly condensed;
- is matched to 10⁻⁶ by the engine's own preparation.

### 2.1 Exported (unblocked) HFa, design H_i = 1.30×10¹¹

Approximate S/Q uses my release (late means) with the base history's vector deposit; the full chain uses mf_ideal's engine,
vectors and ledger.

| c0 | ε | E/E_c | N_YL | S/Q(10⁻⁶) | S/Q(3×10⁻⁶) | S/Q(∞) | e_crit (bulk) | H_req (GeV) |
|---|---|---|---|---|---|---|---|---|
| 1.00 (track) | 0.101 | 21.2 | 11.1 | 14.69 | 12.25 | 10.86 | none | 1.28×10¹¹ |
| 1.10 (approx.) | 0.145 | 13.1 | 12.0 | 8.1 | 6.9 | 6.2 | none | — |
| 1.20 (approx.) | 0.225 | 7.3 | 14.8 | 4.2 | 3.6 | 3.25 | 5.8×10⁻⁶ | — |
| **1.25 (full)** | 0.302 | 4.93 | 23.8 | **3.74** | **3.31** | 2.88 | **1.81×10⁻⁶** | 4.7×10¹⁰ |
| **1.30 (full)** | 0.424 | 3.14 | 25.0 | **2.75** | 2.36 | 1.90 | 3.5×10⁻⁷ | 3.85×10¹⁰ |
| **1.40 (full)** | 0.642 | 1.81 | 6.3 | **1.43** | 1.10 | 0.63 | < 10⁻⁷ | 2.37×10¹⁰ |
| 1.70 (approx.) | 0.79 | 1.37 | 2788 | 20 | 16.8 | 12.2 | none | — |

- **Q-matter window.** It is broad on this history: c0 ≈ 1.27–1.5 computed, with the upper edge between 1.5 and 1.7. HFb (ξ_s mechanism) behaves the same: approximately 3.10
  at c0 = 1.3 and 1.33 at c0 = 1.4.
- **But the dark matter then wants H_i ≈ 2–5×10¹⁰ GeV.** The charge yield grows ×7–25.
- **The axial EFT breaks there [D].** Under the track's own rescaling rule (f fixed in H_i units), μ5,max/4πf =
  0.42–0.49 × (1.30×10¹¹/H_i) reaches 1.4–1.7 at ~4×10¹⁰. So the unblocked strong-holder branch is not consistent.

### 2.2 Pauli-blocked HFa (the corrected history), iterated to self-consistency

**Method.**
- The history is rebuilt with mf_chain's own hist_build.py (Pauli friction to t = 0.45, then the skeptic's fluid, Γ_th = 0,
  tails free, the track's rescaling rule).
- I ran it read-only with its scratch redirected to mine.
- c0 was scanned at each H_i.

| history H_i (K) | c0 | ε | E/E_c | N_YL | S/Q(10⁻⁶) | S/Q(3×10⁻⁶) | e_crit bulk / Q = 2×10⁵ | H_req (GeV) | H_req/H_i |
|---|---|---|---|---|---|---|---|---|---|
| 2.40×10¹¹ (7.84) | 1.00 (track) | 0.078 | 30.1 | 5.5 | 61.5 | 30.2 | none | 2.40×10¹¹ | 1.00 |
| 2.40×10¹¹ | 1.38 | 0.426 | 3.12 | 19.7 | 4.23 | 3.12 | 1.88 / 4.2×10⁻⁶ | 1.12×10¹¹ | 0.47 |
| 2.40×10¹¹ | 1.42 | 0.510 | 2.46 | 13.4 | 3.31 | 2.36 | 0.88 / 1.4×10⁻⁶ | 1.05×10¹¹ | 0.44 |
| 1.10×10¹¹ (6.94) | 1.38 | 0.476 | 2.69 | 18.5 | 3.55 | 2.66 | 1.07 / 1.9×10⁻⁶ | 8.86×10¹⁰ | 0.81 |
| 1.10×10¹¹ | 1.40 | 0.518 | 2.41 | 12.5 | 3.02 | 2.19 | 0.69 / 1.1×10⁻⁶ | 9.27×10¹⁰ | 0.84 |
| **8.40×10¹⁰ (6.50)** | **1.36** | 0.450 | 2.90 | 28.4 | **4.16 (free X)** | **3.28 (Q-matter)** | **2.18 / 7.7×10⁻⁶** | 8.85×10¹⁰ | **1.05** |
| 8.40×10¹⁰ | 1.38 | 0.493 | 2.57 | 17.7 | 3.30 | 2.48 | 0.86 / 1.5×10⁻⁶ | 9.17×10¹⁰ | 1.09 |
| **8.40×10¹⁰** | **1.40** | 0.534 | 2.31 | 13.7 | **2.875 (Q-matter, Q = 2×10⁵ allowed)** | 2.12 | **0.60 / 0.95×10⁻⁶** | **8.38×10¹⁰** | **0.998** |

**Ledger at the converged c0 = 1.40 point (e = 10⁻⁶).** R_vec = 4.03, R_Y = 5.24, R_casc = 6.36, W = 38.4 (inside the closed
form's fitted range 13–140). The relic mass is M = 0.807 m_X [Q qm_micro table].

**Robustness of that point [C].**
- Cycle-mean R instead of raw R: 2.86 (c0 = 1.36: 3.83, e_crit 1.41×10⁻⁶).
- My independent release: ε, N_YL within 0.7%.
- rtol 10⁻¹⁰ vs 10⁻¹²: identical to 10⁻⁸.

**CMB [C].**
- At 8.40×10¹⁰ GeV, X = g²σ_c = 0.187, so the clock is flatter: N_max = 99.8 ≥ N_* = 51.85 at κ = 10.
- r = 1.15×10⁻⁷ (BICEP/Keck 2021, arXiv:2110.00483, r < 0.036 [R]).
- n_s is retuned by the capacity construction, as in the track. The exact Mukhanov–Sasaki spectrum was not recomputed.

**Axial EFT [D].** μ5,max/4πf = 0.65–0.76 at 8.4×10¹⁰ by the track's scaling: inside the cutoff.

**What the point needs (each a model condition, computed):**

1. **A tuned holder.** On this history the Q-matter window is c0 ≈ 1.37–1.45 (an interference minimum of N_YL). It needs the
   phase as well:
   - purely imaginary μ_eff = 12i–16i stays free X at all e (S/Q(∞) ≥ 4.1; `sk_scan_P0H084_imag.json`);
   - at |μ_eff| = 14.29 the window is arg μ_eff ≈ 99°–109° for Q-matter at 10⁻⁶ and ≈ 95°–112° for an e_crit
     (`sk_scan_P0H084_phase.json`). This scan is approximate and reads S/Q about 8% high against the full chain at
     c0 = 1.40, so the true window is somewhat wider.

   So it is a patch of about ±3% in modulus by ±5° in phase. That is comparable to mf_ideal's pocket, at a different place.
2. **Isocurvature [C].** In the Q-matter branch M is fixed, so D = d ln Ω/d ln s_* = 9.2 (charge only, 9.5; c0 = 1.386 /
   1.414 full chains).
   - At the bench's s_* = 7.3×10¹⁵ GeV, β_iso = 0.12, above Planck 2018's 0.038 [R] (retrospective).
   - s_* ≥ 1.36×10¹⁶ GeV (5.6×10⁻³ M_P) is needed. The release does not depend on s_* at fixed |g_h|s_* (mf_ideal §5).
   - f_s rises ×3.5, still ≲ 1% at s decay.
3. **The induced holder portal [C].**
   - c_L = 1/16π²: R_leak = 0.49, S/Q = 2.94 (bulk Q-matter only).
   - c_L = 1/4π²: 7.9, giving 3.88 (free X).
   - Vanish at M_P: free X.

   M1r's benchmark condition (vanish at M_Y, c_L = 0) is needed, or a coefficient ≲ 1/16π².

**Reading.** On the corrected hybrid, the "free X at every e" outcome is a property of the copied holder calibration, not of
the hybrid history. A tuned holder restores the e-controlled branch point at a self-consistent H_i ≈ 8.4–8.9×10¹⁰ GeV. This
is [C] on inherited [E]/[Q] layers:
- the Pauli fluid;
- the homogeneous release;
- the single-Higgs-direction vectors;
- the dq_heat ledger;
- the closed form;
- the qm_micro M.

### 2.3 What else moves with c0 (on the exported history, approximate) [C]

- **N_YL is monotonic in m_s at c0 = 1.** At m_s = 1 / 1.5 / 2 / 2.5 / 3 / 3.5 / 4 / 5 it is 1.9 / 2.9 / 4.9 / 8.3 / 11.1 /
  18.2 / 22.5 / 21.9, and S/Q(∞) ≥ 5.2. So at c0 = 1 the "no e_crit" holds for every m_s tried. This supports C2/C4 at the
  benchmark holder.
- **With the holder held** (s frozen, no turn-on), the background alone gives N_YL = 1.41.

---------------------------------------------------------------------------------------------------------------------

## 3. New issues (ranked)

1. **N1. The holder amplitude and phase are free and decide the branch** [C, §2].
   - mf_chain's and mf_ideal's "no Q-matter / no e_crit" verdicts are c0 = 1 (μ_eff(slice) = 10i) statements.
   - On the Pauli-corrected HFa, a self-consistent point exists with Q-matter at e = 10⁻⁶ (c0 = 1.40, H_i = 8.40×10¹⁰ GeV,
     S/Q = 2.88, e_crit = 0.60×10⁻⁶).
   - Another (c0 = 1.36, H_req/H_i = 1.05) reproduces M1r's branch structure (free X at 10⁻⁶, e_crit = 2.2×10⁻⁶).
   - The window is tuned to about ±3% in |μ_eff| and ±5° in phase.
2. **N2. Isocurvature at the strong-holder point.** D ≈ 9.2 needs s_* ≥ 1.4×10¹⁶ GeV (β_iso = 0.12 at the bench's
   7.3×10¹⁵) [C, retrospective].
3. **N3. The pair-fraction map has an unsampled gap (6.74 < S/Q < 17.1), bridging two recipes that differ ×1.85 at the
   junction** [C].
   - Corrected on the unblocked hybrid: M = 1.77, not 1.58 (57% asymmetric, not 64%).
   - H_i fixed point ≈ 1.11×10¹¹, not 1.28×10¹¹ GeV.
   - D(10⁻⁶) = 2.5, not 1.36.
   - 47 mf_ideal histories shift (H_req −2 to −10%).
4. **N4. β_iso at e = 3×10⁻⁶ was evaluated at the 10⁻⁶ calibration.** At its own fixed point it is 0.0385 (track map) to
   0.047 (gap-filled map, same history): at or above the bound [C].
5. **N5. The unblocked strong-holder branch calibrates to H_i ≈ 2–5×10¹⁰ GeV**, where the axial couplings exceed their EFT
   cutoff (μ5,max/4πf ≈ 1.4–1.7) under the track's rescaling rule [D/C]. Only the Pauli branch gives a consistent strong
   holder.
6. **N6. The leak condition is binding at the strong-holder point.** A threshold-size c_L = 1/4π² gives free X; c_L ≲ 1/16π²
   or vanish-at-M_Y is required [C].
7. **N7. Release sensitivity to the export's sampling switch** at t = 0.6: a ringing interpolation moves N_X by 9% [C]. The
   engine is fine, but exports should overlap densely across the switch.
8. **N8. Gravitational tail** ×1.25–1.5 of the track's value; immaterial (< 1% of E_X) [C].
9. **N9. N_A universality in K holds to about 11%**, and mf_ideal's K^1.45 overestimates N_A at K ≈ 7 by 33% [C].
10. **N10. Checked and negligible [E]: medium-enhanced Y_L → 3X.** Bose enhancement from the X gas could make Y_L decay early.
    But the X quanta sit at comoving k ~ 3, while the Y_L products need k ≈ 7.4a. Before fragmentation (τ ≈ 60) the X is
    coherent, and the 3ω_X = m_L resonance needs a < 1. So the vacuum width and the late deposit stand.

---------------------------------------------------------------------------------------------------------------------

## 4. What stands

- **The engine** (rpmfi.c MODE 5) and the vector code, both independently reproduced.
- **The c0 = 1 branch results** on the exported and Pauli histories, including "no e_crit at any e" for every m_s in 1–5.
- **The leak formula and its numbers; the s-decay entropy; the CMB capacity failure** of HFa at the Pauli fixed point
  (κ = 10).
- **The gravitational-production conclusion.**
- **mf_ideal's c0 = 1 map** and its regression.

## 5. Limits of this pass

- **Approximate scans.** Rows marked "approx." use late means and the base history's vector deposit (−2% in ε against the
  frozen value). Every branch-changing claim was re-run through the full chain.
- **Inherited layers.** The Pauli history (the [E] rates and blocking model of the mf_hybrid skeptic), the homogeneous
  release, the vector model, the ledger, the closed form, the case-L recipe and qm_micro's M are inherited. A strong holder
  makes the release more violent (ε ≈ 0.5, Y partly condensed); fragmentation of that state is not computed.
- **Iteration coverage.** The self-consistent iteration at c0 = 1.36 was not taken to |ratio − 1| < 1%. The c0 = 1.40 point
  is at 0.998. HFb, Γ_th > 0 and band-tail histories at the strong holder were not run (HFb only approximately, on the
  exported history).
- **Isocurvature derivative.** D at the strong holder is a ±1% finite difference inside the Q-matter branch.
- **Not addressed:** Higgs survival at 8.4×10¹⁰ GeV, walls and fragmentation of the waterfall, the exact Mukhanov–Sasaki
  spectrum at the rescaled scales, and C5 (tails blocked).

## 6. Code and results (`skeptic/code`, `skeptic/results`; ≈ 40 min on one or two cores in total)

| file | content |
|---|---|
| `sk_release.py` | independent homogeneous release (own equations, DOP853, two-piece splines, own preparation with continuation, own N_YL read-out; `prep_hessian`) |
| `sk_vec.py` | independent longitudinal-vector modes in e-folds on the skeptic's own F(t) |
| `sk_pmap.py` → `sk_pmap.json` | case-L network inside the p-map gap |
| `sk_ledger.py` → `sk_ledger.json` | independent ledger, S/Q, abundance and H_req, gap-filled relic, D and β_iso |
| `sk_leak.py` → `sk_leak.json` | leak re-derivation and its effect at the strong-holder points |
| `sk_scan.py` → `sk_scan_*.json` | m_s, held-holder and c0 scans with the independent release |
| `sk_common.py`, `sk_chainrun.py`, `sk_c0chain.py`, `sk_idealc0.py` | full chain (mf_ideal's engine, read-only, outputs in scratch) on bench histories, at c0 ≠ 1, and on mf_ideal's fluid corners |
| `sk_summary.py` → `sk_summary.json` | every number above (regressions, vector check, full-chain summaries, isocurvature, CMB) |
| `sk_scan_P0H084_imag.json`, `sk_scan_P0H084_phase.json`, `sk_mfideal_gap_rows.json` | phase dependence; mf_ideal gap rows |

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/mf_chain/skeptic/code; export OMP_NUM_THREADS=1
python3 -B sk_release.py ../../../mf_hybrid/results/bench_HFa.json HFa_raw
python3 -B sk_pmap.py 6.742 8 10 12 14 15 16 17.131 21.797
python3 -B sk_chainrun.py ../../../mf_hybrid/results/bench_HFa.json skHFa; python3 -B sk_c0chain.py skHFa 1.25 1.3 1.4
# Pauli history at 8.4e10 (mf_chain's hist_build, read-only, scratch redirected):
(cd ../../code && MFC_SCRATCH=$SP python3 -B hist_build.py HFa skP0_H084 --H 8.4e10 --fric pauli --Gth 0 --tails free --tsw 0.45)
python3 -B sk_chainrun.py $SP/hist/bench_skP0_H084.json skP0H084; python3 -B sk_c0chain.py skP0H084 1.36 1.38 1.40 1.386 1.414
python3 -B sk_idealc0.py 2 3 0.1 1.0 1.3 1.4; python3 -B sk_ledger.py; python3 -B sk_leak.py; python3 -B sk_summary.py
```

## Literature (checked by web search in this session)

- **Planck 2018 X, arXiv:1807.06211:** n_s = 0.9649 ± 0.0042; isocurvature constraints.
  - The β_iso < 0.038 value for uncorrelated n_II = 1 CDI ("axion I", k = 0.05 Mpc⁻¹) is stated explicitly in Planck 2015
    XX (arXiv:1502.02114).
  - For 2018 I rely on later literature quoting the same number, and on mf_ideal and the mf_hybrid skeptic.
- **BICEP/Keck 2021, arXiv:2110.00483:** r₀.₀₅ < 0.036.
- **ACT DR6, arXiv:2503.14452:** n_s = 0.974 ± 0.003 (P-ACT-LB).
- **Ema–Jinno–Mukaida–Nakayama, arXiv:1502.02475:** gravitational production by an oscillating background.
- **Graham–Mardon–Rajendran, arXiv:1504.02102:** longitudinal vectors from inflation.
- **Our tracks, read-only:**
  - mf_chain (hist_build, cmb_rescale, results);
  - mf_ideal (engine, drivers, map);
  - mf_hybrid and its skeptic (bench files, sk_structure, sk_pauli tables);
  - fx_chain (m1_kin case-L network);
  - inf_chain (pmap_kin rows), wt_critic (gas map), dq_heat skeptic (dq_model_copy), qm_micro (M table).
