# Track mf_ideal: the M1r chain on idealised fast histories with a spectator holder

2 November 2026 · Claude, multi-field inflaton pass (BRIEF_MULTIFIELD.md of this directory; GitHub issue #13 context).

This is same-session AI work on one supplied model, not outside replication. Nothing here is evidence for an inflaton
sector, a spectator, the hidden sector or any dark-matter branch. Nothing is a validation, an exclusion, a new theory or a
GRUT prediction. **Every history below is a model choice, and none follows from closure.** That covers the background, K,
the abrupt or smooth end, the R(t) treatment, the ramp of m_s, Γ_s, s_*, the holder calibration μ_eff(slice) = 10i and
the ledger. Data comparisons (Planck 2018 isocurvature) are retrospective.

Labels: [D] derived, [C] computed with this track's code (`code/`, output in `results/`), [E] estimate, [Q] quoted from
a reviewed track, [A] argued, [R] literature (arXiv numbers checked by web search in this session).

Git was read-only and nothing was committed. All runs were in the foreground, with at most two worker processes (≈ 1 h
of CPU). Binaries and series are in scratch (`$SP/wt2/mf_ideal/`). No compiled file or `__pycache__` is in the tree.

## 0. Result in brief

**The map** (`results/sq_map.json`, 152 histories; primary grid in §3):
- **At c0 = 1, the brief's μ_eff(slice) = 10i calibration, no fast spectator history gives Q-matter at e = 10⁻⁶ [C].**
  This covers all 30 primary histories (K = 1.0–3.0, m_s = 2–4, t_on = 0.1, 1) and the 81 at c0 = 1 in the other groups:
  - the fine K = 1 scan;
  - ξ_s tracking;
  - small or fast s decay;
  - oscillating-matter backgrounds;
  - smooth ends;
  - s_* variants.
  Primary S/Q(10⁻⁶) = 4.97–53.9; over all c0 = 1 histories the minimum is 4.38.
- **The best corners sit at M1r's own level, S/Q(10⁻⁶) = 4.4–5.0:**

  | corner | S/Q(10⁻⁶) | e_crit (bulk) |
  |---|---|---|
  | K = 1, m_s = 3, t_on = 0.5 | 4.49 | 1.9×10⁻⁶ |
  | K = 1, m_s = 2.5, t_on = 0.5 | 4.52 | 2.4×10⁻⁶ |
  | t_end = 0.03 smooth end, t_on = 0.5 | 4.38 | 1.7×10⁻⁶ |
  | m_bg = 3 oscillator, K = 2, t_on = 1 | 4.82 | 1.8×10⁻⁶ |
  | K = 1, m_s = 4, t_on = 1 | 4.97 | 2.2×10⁻⁶ |

  M1r's static slice gives 4.955 and 2.06×10⁻⁶. These corners are Q-matter at e = 3×10⁻⁶.
- **e_crit (bulk) in the primary grid exists only in three histories:**
  - K = 1, m_s = 4, t_on = 1: 2.19×10⁻⁶ (Q = 2×10⁵: 3.5×10⁻⁶);
  - K = 1.3, m_s = 3, t_on = 1: 8.2×10⁻⁶;
  - K = 1, m_s = 3, t_on = 1: 6.3×10⁻⁴.
  The other 27 have none at any e. Fast turn-on (t_on ≤ 0.1) never has one, and neither does any K ≥ 1.6.
- **One narrow pocket gives Q-matter at e = 10⁻⁶ [C], and only with a holder 3–5% stronger than M1r's calibration.**
  - The pocket: K = 1 (instant radiation), m_s = 3, t_on = 0.5, c0 = 1.03–1.05, i.e. |g_s s_*| 3–5% above the value for
    μ_eff = 10i.
  - S/Q(10⁻⁶) = 3.36 / 3.31, e_crit = 0.94 / 0.90×10⁻⁶.
  - The neighbours are all free X: c0 = 1.00 (4.49), 1.01 (4.03), 1.07 (4.20), 1.10 (7.4); and the other corners at
    c0 = 1.03–1.07 (3.6–7.6).
  - It is an interference minimum of the Y_L excitation (N_YL = 0.018 at c0 = 1.03). It is converged numerically (rtol
    10⁻¹⁰ → 10⁻¹² moves N_X by 10⁻¹² and N_YL by 10⁻⁶), but it is fine-tuned in (K, t_on, c0).
- **Why the fast histories do not reach Q-matter [C]** (§4):
  1. Without M1r's oscillating inflaton, R(t) no longer pumps X through ξ_X R|X|². E_X stays at 5.7–7.1×10⁴, against
     M1r's 1.89×10⁵. ε falls to 0.05–0.15, and the vector heat per unit E_X rises.
  2. The abrupt end itself (R stepping from 12 to ρ_m) excites Y_L: N_YL ≈ 1.9 with the holder held.
  3. A fast mass turn-on (t_on ≤ 0.1, including the ξ_s R mechanism) adds to (2) constructively: N_YL = 2–19.
- **G7's "×4 cut of the vector heat at K → 1" does not carry over to S/Q [C].**
  - N_A × deposit does scale as K^1.95.
  - But in units of the release's own E_X, R_vec(10⁻⁶) falls only from 3.61 (M1r) to 2.6–2.8 at K = 1, because E_X
    falls ×2.7.
- **Isocurvature [C].**
  - On the map, D = d ln Ω_DM/d ln s_* = −1.6 to +5.0 at e = 10⁻⁶; the charge alone gives 4.4–5.2.
  - In the Q-matter pocket D ≈ 6.2 (charge only).
  - Planck 2018 then needs s_* > 1.75×10⁴ |D| H_i [R: β_iso < 0.038, arXiv:1807.06211], i.e. s_* ≳ 1.1×10¹⁶ GeV in the
    pocket.
  - Near e_crit a 1% shift in s_* crosses the branch, so the response there is not smooth.
- **Calibration [C].**
  - The near-threshold free-X histories need H_i = 1.0–1.3×10¹¹ GeV; the pocket's Q-matter needs 1.01–1.09×10¹¹ GeV.
    That gives r = 1.7–1.9×10⁻⁷ at H_* = H_i.
  - The fast-turn-on histories have a partly symmetric relic, M = 1.1–6 m_X per charge.
- **Bench interface [C]** (§9).
  - `bench_run.py` runs the whole chain on mf_hybrid's `bench_HFa.json`. The evolved s matches the bench's s(t) to
    5×10⁻³, and the pre-slice exit factor is f = 1.00000002.
  - As a smoke test only (mf_chain owns the analysis), HFa lands in the fast-turn-on corner: N_YL = 11, S/Q(10⁻⁶) = 14.7,
    free X, no e_crit.

**Plain answer to the brief.**
- With the holder calibrated as M1r does it, no idealised fast spectator history gives Q-matter at e = 10⁻⁶.
- e_crit sits at 1.7–2.4×10⁻⁶ in the best corners. Those need K ≈ 1, a turn-on time t_on ≈ 0.5–1 H_i⁻¹ (not fast), and
  an end of inflation that is either abrupt with an averaged R or M1r-like.
- Elsewhere e_crit is ≥ 7×10⁻⁶ or absent.
- Q-matter at 10⁻⁶ needs, in addition, a holder 3–5% stronger and a tuned turn-on, inside a pocket a few percent wide.

## 1. Engine and regression (`code/rpmfi.c`, `code/regress.py` → `results/regression.json`)

**The engine.** `rpmfi.c` is a copy of inf_chain's `rpchain.c` (lineage: rpinf.c ← dq_entropy rpdq.c ← fx_holder
rptri.c). The additions are marked MF0–MF2:
- **MF0, the stage time.**
- **MF1, MODE 5:** an external background (a built-in decaying-dust fluid, or a table `t lna H R [c cdot [ms2]]`, with
  Hermite interpolation and a radiation tail) and the holder u = μ + G_eff c on a spectator c = s/s_*.
  - The spectator either evolves by c'' + (3H + Γ_s)c' + m_s²(t)c + (H_i/s_*)² ∂V_h/∂c = 0, with the exact back-reaction
    of the holder on s and the hidden-work ledger, or follows the table.
  - m_s²(t) is a smootherstep ramp over t_on, or the table column, plus an optional ξ_s R.
- **MF2:** a state dump.

**Checks.**

| check | result |
|---|---|
| R1: MODE 0, (3, 0.5), τ ≤ 1000 | stdout **byte-identical** to rpdq.c and rpchain.c. N_X(400–1000) = −8045.485140302082, N_YL = 0.566640911258934, E_X = 188706.29275778163 (the dq_entropy / fx_holder digits) [C] |
| R2: the table bridge (the bench path). M1r's own background and holder c = I/I_0 written as a 7-column table, fed to MODE 5 | late means agree with MODE 0 to 8×10⁻⁸ (Δt = 0.002) and 1×10⁻⁸ (Δt = 0.001): cubic convergence [C] |
| Chain on the bridge: own vectors, ledger, closed form | S/Q(10⁻⁶) = 4.955, e_crit = 2.06×10⁻⁶ (bulk) / 3.18×10⁻⁶ (Q = 2×10⁵), R_vec = 3.614, R_Y = 0.379. This is the decisive pass's 4.955 / 2.06 / 3.19 [C] |
| Vectors | vmodes.c on the bridge reproduces inf_chain's M1r validation exactly (N_A = 2.0302), within 0.13% of dq_heat [C] |
| Every history | Gauss ≤ 2×10⁻⁹, hidden work ≤ 5×10⁻¹¹, H/ln a consistency ≤ 2×10⁻⁶; N_X frozen to ≤ 0.6% (a⁻² fits to τ = 10⁴; ≤ 0.3% except the slow-turn-on t_on ≥ 2 rows) [C] |
| Hartree vs closed form (`results/hartree_check.json`) | max \|rel\| = 1.3×10⁻⁴ for W = 16–3832 on two releases. W = 6866 (the K = 21 smooth end only) exceeds the condensed solver, where inf_chain's no-condensate check applies [Q] |

## 2. The idealised histories (model choices)

**Preparation.** Constant-H inflation up to the slice, with s frozen at s_*. The preparation is the M1r held stationary
point at u = μ + (10i − μ)c0 and R = 12, which for c0 = 1 is (90.800, 6.888, 28.576) H_i. With constant H the
frozen-angle exit factor is exactly f = 1 (the static slice).

**Background (abrupt end).**
- At the slice ρ = 3 sits in dust (R = ρ_m), which decays into radiation at Γ_bg. The R step from 12 to 3 (or to 0) is
  the heavy-waterfall limit, m_χ ≫ H, where the hidden sector sees the cycle-mean R.
- K = a²H in radiation, tuned by Γ_bg [C, `mfi_common.gam_for_K`]:

  | K | 1.0 | 1.3 | 1.6 | 2.0 | 3.0 |
  |---|---|---|---|---|---|
  | Γ_bg (H_i) | ∞ (radiation from the slice) | 1.0829 | 0.39651 | 0.16946 | 0.044288 |

**Spectator.**
- m_s² = m_s² S(t/t_on), with S the smootherstep 6x⁵ − 15x⁴ + 10x³.
- m_s ∈ {2, 3, 4} H_i, t_on ∈ {0.1, 1} H_i⁻¹ (the fine scan adds 0, 0.3, 0.5, 2, 3).
- Γ_s = 0 (primary), 0.05 or 0.5.
- s_* = 10¹⁶ GeV (back-reaction (H_i/s_*)² = 1.8×10⁻¹⁰). 10¹⁴ and 10¹⁸ are checked in §5.

## 3. The primary map S/Q(K, m_s, t_on, e) [C]

Columns: S/Q at e = 10⁻⁷ / 10⁻⁶ / 3×10⁻⁶ / 10⁻⁵ / 10⁻⁴ / ∞, then e_crit (bulk; Q = 2×10⁵), N_YL and ε. Full rows, with
R_vec, R_Y, W, the relic, H_i and the checks, are in `results/sq_map.json` (group "primary").

| K | m_s | t_on | S/Q (10⁻⁷ … ∞) | e_crit | N_YL | ε |
|---|---|---|---|---|---|---|
| 1 | 2 | 0.1 | 21.0 / 8.73 / 7.59 / 7.16 / 6.97 / 6.94 | – | 5.28 | 0.127 |
| 1 | 2 | 1 | 24.9 / 7.58 / 5.82 / 5.15 / 4.86 / 4.81 | – | 1.82 | 0.098 |
| 1 | 3 | 0.1 | 28.3 / 15.1 / 13.8 / 13.3 / 13.0 / 12.9 | – | 11.3 | 0.116 |
| 1 | 3 | 1 | 16.4 / 5.24 / 4.12 / 3.69 / 3.50 / 3.46 | 6.3×10⁻⁴; – | 2.22 | 0.151 |
| 1 | 4 | 0.1 | 46.2 / 28.0 / 25.8 / 24.6 / 23.9 / 23.6 | – | 19.3 | 0.094 |
| 1 | 4 | 1 | 21.3 / **4.97** / 3.09 / 2.28 / 1.89 / 1.81 | **2.19×10⁻⁶; 3.53×10⁻⁶** | 0.010 | 0.109 |
| 1.3 | 2 | 0.1 | 36.0 / 10.0 / 7.35 / 6.34 / 5.92 / 5.87 | – | 2.29 | 0.102 |
| 1.3 | 2 | 1 | 39.7 / 9.13 / 5.80 / 4.48 / 3.93 / 3.86 | – | 0.80 | 0.092 |
| 1.3 | 3 | 0.1 | 49.8 / 19.2 / 16.2 / 15.0 / 14.5 / 14.3 | – | 6.51 | 0.082 |
| 1.3 | 3 | 1 | 39.0 / 8.24 / 4.76 / 3.32 / 2.70 / 2.61 | 8.2×10⁻⁶; 2.4×10⁻⁵ | 0.24 | 0.092 |
| 1.3 | 4 | 0.1 | 61.2 / 31.4 / 28.3 / 26.9 / 26.0 / 25.7 | – | 15.2 | 0.083 |
| 1.3 | 4 | 1 | 53.3 / 11.7 / 7.11 / 5.23 / 4.42 / 4.31 | – | 0.47 | 0.068 |
| 1.6 | 2 | 0.1 | 53.0 / 13.2 / 9.06 / 7.45 / 6.81 / 6.74 | – | 2.21 | 0.094 |
| 1.6 | 2 | 1 | 57.4 / 12.3 / 7.30 / 5.29 / 4.46 / 4.36 | – | 0.83 | 0.088 |
| 1.6 | 3 | 0.1 | 70.2 / 23.9 / 19.3 / 17.5 / 16.7 / 16.6 | – | 6.40 | 0.077 |
| 1.6 | 3 | 1 | 62.1 / 12.7 / 7.05 / 4.73 / 3.74 / 3.62 | – | 0.41 | 0.080 |
| 1.6 | 4 | 0.1 | 78.9 / 36.4 / 32.0 / 30.2 / 29.1 / 28.8 | – | 15.1 | 0.081 |
| 1.6 | 4 | 1 | 85.8 / 18.0 / 10.4 / 7.25 / 5.92 / 5.76 | – | 0.59 | 0.058 |
| 2 | 2 | 0.1 | 75.2 / 17.1 / 10.9 / 8.48 / 7.52 / 7.43 | – | 2.19 | 0.092 |
| 2 | 2 | 1 | 81.1 / 16.3 / 9.05 / 6.07 / 4.85 / 4.73 | – | 0.84 | 0.086 |
| 2 | 3 | 0.1 | 96.9 / 29.1 / 22.3 / 19.7 / 18.6 / 18.4 | – | 6.37 | 0.075 |
| 2 | 3 | 1 | 91.3 / 17.8 / 9.39 / 5.88 / 4.39 / 4.24 | – | 0.48 | 0.076 |
| 2 | 4 | 0.1 | 102 / 41.6 / 35.5 / 33.0 / 31.7 / 31.3 | – | 15.1 | 0.081 |
| 2 | 4 | 1 | 128 / 25.4 / 13.7 / 8.92 / 6.90 / 6.69 | – | 0.63 | 0.054 |
| 3 | 2 | 0.1 | 137 / 27.2 / 15.2 / 10.4 / 8.70 / 8.57 | – | 2.18 | 0.090 |
| 3 | 2 | 1 | 147 / 27.1 / 13.4 / 7.70 / 5.49 / 5.32 | – | 0.84 | 0.086 |
| 3 | 3 | 0.1 | 170 / 41.9 / 28.7 / 23.7 / 21.8 / 21.6 | – | 6.35 | 0.074 |
| 3 | 3 | 1 | 170 / 30.7 / 14.7 / 7.89 / 5.16 / 4.95 | – | 0.52 | 0.074 |
| 3 | 4 | 0.1 | 165 / 53.9 / 42.8 / 38.5 / 36.6 / 36.2 | – | 15.1 | 0.082 |
| 3 | 4 | 1 | 240 / 43.7 / 21.2 / 11.8 / 8.07 / 7.77 | – | 0.66 | 0.052 |

**Branches.**
- At e = 10⁻⁶: free X in 30 of 30.
- At e = 3×10⁻⁶: Q-matter only for K = 1, m_s = 4, t_on = 1 (S/Q = 3.09, bulk only, since > 2.93).

**Ledger at K = 1** (E_X units of each release): R_vec(10⁻⁶) = 2.65–2.78, R_Y = 1.31 N_YL, R_casc ≤ 4.2. The vector number
depends only on K [C]:

| K | 1 | 1.3 | 1.6 | 2 | 3 |
|---|---|---|---|---|---|
| N_A(10⁻⁶) | 0.78 | 1.10 | 1.46 | 2.00 | 3.45 |

That is N_A(k ≤ 1) ∝ K^1.45, plus a release band of 0.04–0.08 at k ≥ 1. With the deposit ∝ K^½ the total is ∝ K^1.95.

## 4. What sets S/Q on fast histories [C, `code/` inline diagnostics]

**(a) The background's R(t) controls the release.**
- M1r's own background, with the spectator holder (m_s = 3, Γ_s = 0.5, t_on = 0), reproduces M1r to all printed digits.
  With Γ_s = 0 (switch-off by redshift only) it gives N_X = −7978 against −8079, and N_YL = 0.58 against 0.57.
- The same holder on a dust fluid with the same Γ = 0.5 (K = 1.52) gives N_X = −2340, E_X = 6.3×10⁴ and N_YL = 6.6.
- So it is not the decay-versus-redshift switch-off that changes the release (G7's worry). It is the absence of the
  oscillating R. In M1r, R = ρ(1 + 3cos 2m_I t), between 12 and −6 at frequency 6, which pumps X through ξ_X R|X|²
  (E_X 6.8 → 18.9×10⁴).

**(b) The abrupt end excites Y_L.**
- With the holder held (m_s = 0), K = 1: N_YL rises to 8 by t = 0.2 and settles at 1.9.
- On M1r's continuous-R background the held value is 0.35.

**(c) The turn-on interferes with (b).** At K = 1, m_s = 3:

| t_on | 0 | 0.1 | 0.3 | 0.5 | 1 | 2 | 3 |
|---|---|---|---|---|---|---|---|
| N_YL | 11.6 | 11.3 | 2.56 | 0.58 | 2.22 | 1.51 | 1.88 |

Slow turn-on also lowers ε (0.075 at t_on = 3), because the held phase persists while X relaxes.

**(d) The oscillator background (`osc_bg.py`) shows the crossover.** Matter of mass m_bg released at rest, K = 2, m_s = 3:

| m_bg | ε (t_on = 0.1 / 1) | N_YL | S/Q(10⁻⁶) |
|---|---|---|---|
| 3 (M1r-like) | 0.121 / 0.124 | 0.57 / 0.075 | 5.15 / 4.82 (e_crit 2.3 / 1.8×10⁻⁶) |
| 6 | 0.10 / 0.11 | 3.7 / 6.5 | 12.8 / 15.5 |
| 10 (2m_bg ≈ m_YL) | 0.081 / 0.085 | 27 / 14 | 52 / 34 |
| 30 | 0.076 / 0.077 | 7.3 / 0.77 | 30.3 / 18.1 (= the fluid limit) |

**K_min of an oscillator released at rest [C].** Such an oscillator cannot reach small K at all:

| m_bg | 3 | 6 | 10 | 30 |
|---|---|---|---|---|
| K_min (at Γ_bg) | 1.86 (1.1) | 1.31 (2.7) | 1.17 (5) | 1.09 (5) |

So "quiet (m ≈ 3) and K → 1" is not available from a single oscillating field. K ≈ 1 needs heavy matter or instant
conversion, and that is the R-step regime.

**(e) Smooth end (`smooth_end.py`).** Vacuum converts to radiation over t_end, so R is continuous:

| t_end | K | S/Q(10⁻⁶) at t_on = 0.1 / 0.5 / 1 |
|---|---|---|
| 0.03 | 1.06 | 14.0 / **4.38** / 5.35 (e_crit 1.67×10⁻⁶ at t_on = 0.5) |
| 0.1 | 1.24 | 9.4 / 5.96 / 5.04 |
| 0.3 | 2.11 | 12.3 / 17.7 / 23.7 |
| 1 | 21.1 | 93–144 (inflation continues; not a fast history) |

## 5. Other sensitivities [C]

| variant | result |
|---|---|
| ξ_s tracking (G5 b: m_s² = m_s²(1 − R/12), all K, m_s) | m_s turns on at once at an abrupt end, so this equals t_on = 0: N_YL = 2–19, S/Q(10⁻⁶) = 8.9–43.8, no e_crit |
| small decay Γ_s = 0.05 (m_s = 3, all K) | within 4% of Γ_s = 0 (K = 1, t_on = 1: 5.04 vs 5.24) |
| fast decay Γ_s = 0.5 (K = 1, m_s = 3, t_on = 0.5) | S/Q 4.49 → 5.95: worse |
| s_* = 10¹⁴ / 10¹⁸ GeV | S/Q(10⁻⁶) = 4.45 / 4.49 (vs 4.49 at 10¹⁶): the release is s_*-independent above ~10¹⁵ GeV |
| mirror s_* < 0 (c0 = −1: u_slice = 2μ − 10i) | Y-condensed preparation; N_X = **+35654** (opposite sign, ×9), N_YL = 857, S/Q = 38. The reflection s → −s flips the sign but **is not the CP image** [C; cf. inf_reheat R1, inf_chain's Send flip] |
| holder amplitude c0 at K = 1, m_s = 3, t_on = 0.5 | S/Q(10⁻⁶) = 8.35 / 4.97 / 4.49 / 4.03 / **3.36** / **3.31** / 4.20 / 7.37 / 22.6 at c0 = 0.90 / 0.99 / 1.00 / 1.01 / 1.03 / 1.05 / 1.07 / 1.10 / 1.20. Above 1.07 Y_L is driven hard (N_YL = 12 at 1.10, 168 at 1.20) |
| c0 = 1.03–1.07 at other corners | K = 1, m_s = 2.5, t_on = 0.5: 4.06 / 3.81 / 3.61 (free X, e_crit 1.1–1.6×10⁻⁶). K = 1, m_s = 4, t_on = 1: 5.2–7.6. K = 1.3, m_s = 3, t_on = 1: 5.4–6.8 |
| rtol 10⁻¹² vs 10⁻¹⁰ in the pocket (`results/rtol_check_*.json`) | N_X 10⁻¹², E_X 10⁻¹², N_YL 10⁻⁶, N_YH 2×10⁻⁵ |

## 6. Branch, relic and H_i calibration [C on the inherited recipes]

These use dq_heat's ledger, the closed form, wt_critic's / inf_chain's p(S/Q) map with the self-consistent yield
correction [Q], qm_micro's M for Q-matter [Q], and Y ∝ H_i^{3/2}.

| history | e = 10⁻⁶: branch, M (m_X), H_i,req (GeV), r | e = 3×10⁻⁶ |
|---|---|---|
| M1r bridge | free X, 1.00, 1.29×10¹¹, 2.7×10⁻⁷ | Q-matter |
| K = 1, m_s = 3, t_on = 0.5 | free X (×1.29), 1.00, 1.06×10¹¹, 1.8×10⁻⁷ | Q-matter (bulk) |
| K = 1, m_s = 4, t_on = 1 | free X (×1.43), 1.00, 1.21×10¹¹ | Q-matter (bulk) |
| pocket, c0 = 1.03 / 1.05 | **Q-matter**, 0.811, 1.09 / 1.01×10¹¹, 1.9 / 1.7×10⁻⁷ | Q-matter (Q = 2×10⁵ allowed) |
| fast turn-on (t_on ≤ 0.1) | free X, partly symmetric: M = 1.1–6.0, H_i,req = 0.6–2.1×10¹¹ | free X |

In the free-X corners with M ≈ 1 the relic is asymmetric (p ≲ 10⁻²), so α + sign(s_*) sets it. With fast turn-on,
1/M = 17–90% of Ω_DM is the asymmetric part.

## 7. Isocurvature (`results/isocurvature.json`) [C, with R for the bound]

**Method.** D = d ln Ω_DM/d ln s_* at fixed g_s and H_i, from c0 = 1 ± 0.01 releases (preparation re-solved) with their
own vectors. Ω_DM ∝ M|N_X|K^{−3/2}.

| (K, t_on), m_s = 3 | (1, 0.1) | (1, 0.5) | (1, 1) | (1.3, 0.1) | (1.3, 1) | (1.6, 0.1) | (1.6, 1) | (2, 0.1) | (2, 1) | (3, 0.1) | (3, 1) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D at 10⁻⁶ | 4.43 | 4.99 | 4.42 | 2.77 | 3.35 | 2.56 | −0.99 | 2.10 | −1.55 | 1.46 | −0.72 |
| charge-only | 4.91 | 5.07 | 4.42 | 4.41 | 4.91 | 4.39 | 4.88 | 4.54 | 5.02 | 4.67 | 5.19 |

- The charge responds as N_X ∝ s_*^{4.4–5.2}. A larger μ_eff gives a larger kick.
- In the partly symmetric relics, M falls as the charge grows and partly cancels it.
- In the pocket, D ≈ 6.2 (charge only, from c0 = 1.03 → 1.05, with M fixed in the Q-matter branch).

**The bound.**
- Planck 2018 (arXiv:1807.06211 [R]): β_iso < 0.038 at 95% CL for uncorrelated n_II = 1 CDM isocurvature, at
  k = 0.05 Mpc⁻¹. A spectator that is massless during inflation gives n_II ≈ 1 and no correlation.
- Hence s_* > |D| H_i/(2π√(βA_s/(1 − β))) = 1.75×10⁴ |D| H_i.
  - In the near-threshold corners (D ≈ 5): s_* ≳ 9×10¹⁵ GeV.
  - In the pocket: s_* ≳ 1.1×10¹⁶ GeV.
  - The fiducial 10¹⁶ GeV is at the edge.
- **Near e_crit the response is not smooth.** At K = 1, t_on = 0.5, e = 3×10⁻⁶, c0 = 0.99 is free X (S/Q 3.50) and
  c0 = 1.01 is Q-matter (2.65). With δs/s_* = H/2πs_* ≈ 2×10⁻⁶ this needs S/Q within ~10⁻⁵ of the threshold to matter,
  but such a model sits on a knife edge in s_* [A].

## 8. The brief's estimates G1–G7: verified or refuted here

| | verdict | basis |
|---|---|---|
| G1 | arithmetic **verified**: √6 H M_P/m = 0.816 M_P at m = 3H [D]. Blocking and loop consequences not computed here (mf_gemini) | `summary.py` G_checks |
| G2 | arithmetic **verified** [D]: v_χ = 1.633 M_P; λ_χ = 4.8×10⁻¹⁵ at 1.3×10¹¹ GeV (5.1×10⁻¹⁵ at 1.343×10¹¹); y⁴/16π²λ = 2.0×10¹³; yv/(m/2) = 3.9×10⁷; m²(χ = 0) = −4.5 H². The waterfall duration, walls and holder placement are not computed here | same |
| G3, G4, G6 | not addressed (mf_hybrid / mf_gemini) | — |
| G5 | coefficients **verified** [D]: s_* ≳ 2.27×10¹⁵ GeV × \|D\| at H_i = 1.3×10¹¹ (β_iso 0.038 checked [R]); ρ_s/ρ = (m_s²/6)(s_*/M_P)² = 1.5(s_*/M_P)²; ξ_s = −3/4. **D computed** [C]: \|D\| = 0.7–5.0 (charge-only 4.4–5.2; pocket 6.2), so s_* ≳ 10¹⁶ GeV. "The CP sign is α + sign(s_*)" **holds in sign only**: s_* → −s_* is not the CP image (§5). Option (b) ξ_s turns m_s on instantly at an abrupt end, which is the worst corner (§5) | §5, §7 |
| G7 | "R_vec ∝ K²" **verified** for N_A × deposit (∝ K^1.95). "K ≥ 1 from the slice" verified; an oscillator released at rest has K ≥ 1.86 at m = 3 [C]. "×4 cut at K → 1–1.3" **refuted for S/Q**: per unit of the release's own E_X, R_vec falls only ×1.3, and ε drops. "The release changes when the holder switches off by turn-on and redshift": **the switch-off mode barely matters** (−1.3% in N_X on M1r's background); **the background's R(t) and the abruptness of the end change it** | §3, §4 |

## 9. Bench interface (`code/bench_run.py`, `code/README.md`) and a smoke test on bench_HFa [C, interface test only]

- **What bench_run.py does.** It reads `history.{t, N, H, R, m_s, s_over_sstar, s_is_envelope, R_osc_sig}` and
  `params.{H_i_GeV, s_star_GeV}`, writes a 7-column MODE 5 table for t ≥ 0, and evolves s with the bench's m_s²(t)
  (`MFI_SM2TAB`) and Γ_s. It then runs release, vectors and chain.
- **Checks on HFa.**
  - The evolved s agrees with the bench's exact s(t ≤ 60) to 5.2×10⁻³.
  - The pre-slice frozen-angle factor is 1.00000002 at 1–4 e-folds (ε_H ≈ 8×10⁻⁸).
  - The bench resolves the 2m_χ oscillation of R (Δt = 4×10⁻⁴), so raw R is used. The cycle-mean variant (`--R trace`)
    moves S/Q by 7%.
- **Smoke result.** K = 2.048, ε = 0.101, N_YL = 11.1, S/Q(10⁻⁶) = 14.7 (cycle-mean R: 15.7). Free X, no e_crit.
  - This is consistent with the map: HFa's spectator mass comes from the waterfall over ~0.2 H_i⁻¹ (fast turn-on), and its
    end is abrupt.
  - This is not mf_chain's result. The H_i iteration, s and σ decay heat, gravitational production and remnants remain
    with mf_chain.

## 10. For mf_chain and the critic: what a fast spectator history needs [C/A]

1. **A slow mass turn-on.** t_on ≈ 0.5–1 H_i⁻¹, against the Y_L period of ~0.3 H_i⁻¹. The waterfall-mass (m_s = g_s χ)
   and ξ_s R mechanisms both switch on with the waterfall or the R drop, which is fast in a heavy hybrid. A separately
   slow channel would be needed.
2. **K ≈ 1, with E_X not reduced.** Either M1r's oscillating R (m_bg ≈ 3, but then K ≥ 1.86), or a smooth end over
   t_end ≲ 0.03 H_i⁻¹ with K ≈ 1.
3. **A holder amplitude tuned to within a few percent.** This is needed for Q-matter at e = 10⁻⁶; otherwise free X with
   e_crit ≈ 2×10⁻⁶.
4. **s_* ≳ 10¹⁶ GeV** for isocurvature.
5. **s must decay before it dominates** [E]: at a ≈ 6K²/[m_s²(s_*/M_P)²] = 4×10⁴ (K = 1, s_* = 10¹⁶ GeV), i.e.
   Γ_s ≳ 6×10⁻¹⁰ H_i.
6. **The induced holder portal A_sΣ s|Σ|² must be tiny** [E, inf_chain's leak formula Γ = A²m_s³/32π(m_σ² − m_s²)²].
   - With a threshold-size coefficient (c = 1/16π² to 1/4π²) it sends s quanta into long-lived A_L, with
     R_leak ≈ 0.36–5.8 × (1.3×10⁻⁶ H_i/Γ_s). This is independent of s_*.
   - With HFa's Γ_s = 1.3×10⁻⁶ this is comparable to or larger than the whole Q-matter budget (R_crit ≈ 2).
   - Either Γ_s ≳ 10⁻⁴ H_i or the "vanish at M_Y" condition (A = 0) is needed.

## 11. Limits

- **The model.** Homogeneous, classical release; the single-Higgs-direction Proca vector model; dq_heat's ledger (decay
  at rest, full deposit); the inherited closed form, p(S/Q) recipe and qm_micro M table [Q].
- **R treatment.**
  - The fluid uses the cycle-mean R (heavy-waterfall limit). Gravitational particle production at the R step
    (ξ_X = −5/6) is not computed. My estimate [E]: n ~ 0.05 H³, ≲ 10⁻⁶ of E_X.
  - Parametric effects of an R oscillation at 2m_χ ≫ 40 H_i are neglected [A].
- **Ramp shape.** The smootherstep ramp is one choice. N_YL depends on the interference phase, so other ramp shapes move
  the t_on optimum.
- **Grid coverage.** The pocket is resolved only on a coarse c0 grid, at one (m_s, t_on) corner and at K = 1 exactly.
- **Not computed here:** hidden heat from s decays (only the [E] leak bracket); the H_i iteration of a real potential;
  Higgs survival; fragmentation; finite-k Y production.

## 12. Code and results

| file | content |
|---|---|
| `code/rpmfi.c`, `code/vmodes.c` | engine (MODE 5) and vector modes (copy) |
| `code/mfi_common.py` | driver, preparation, fluid K(Γ) |
| `code/release_mfi.py`, `vectors_mfi.py`, `chain_mfi.py` | release, vectors, ledger/S/Q/e_crit/relic/H_i/isocurvature/Hartree |
| `code/grid.py`, `osc_bg.py`, `smooth_end.py`, `extras.py` | history grids and sensitivity runs |
| `code/bench_run.py`, `code/README.md` | **bench_*.json interface and usage** |
| `code/regress.py`, `code/summary.py` | regression; map and summary |
| `results/release/`, `results/vectors/`, `results/chain/` | per-history JSON (152 histories) |
| `results/sq_map.json` | the map S/Q(K, m_s, t_on, e) by group, with the relic, H_i and checks |
| `results/summary.json` | headline numbers, regression, isocurvature, G checks, bench smoke |
| `results/regression.json`, `hartree_check.json`, `isocurvature.json`, `rtol_check_*.json`, `bench/` | checks |

To reproduce, see `code/README.md`. The full set is ≈ 1 h on two cores.

## Credits

- **Literature [R]** (arXiv numbers checked by web search in this session):
  - Planck 2018 X, arXiv:1807.06211 (β_iso < 0.038, uncorrelated n_II = 1, k = 0.05 Mpc⁻¹, TT,TE,EE+lowE+lensing);
  - Linde, hybrid inflation, astro-ph/9307002;
  - Felder, García-Bellido, Greene, Kofman, Linde, Tkachev, tachyonic preheating, hep-ph/0012142 (in hybrid inflation the
    oscillation stage ends after about one oscillation: the abrupt-end regime);
  - Graham–Mardon–Rajendran, longitudinal vector production, arXiv:1504.02102.
- **Our tracks, read-only:**
  - inf_chain: rpchain.c, vmodes.c, release/vectors/chain methods, pmap_kin.json;
  - dq_heat skeptic: dq_model_copy.py;
  - dq_entropy: rpdq.c reference, closed form, hartree_dq.py;
  - fx_holder: rptri.c lineage;
  - wt_critic: chain_consistency.json;
  - dq_critic / qm_micro: the Q-matter M table;
  - mf_hybrid: bench_HFa.json, used for the interface test only.
