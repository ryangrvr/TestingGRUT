# Track rp_critic: consolidation of rp_reheating, rp_preparation and rp_portals (Codex e160356)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. The consolidated review, verdict
tables, corrections, overreach flags, status change and the decisive next calculation are in
[CONSOLIDATED.md](CONSOLIDATED.md). This file records method, own numbers and limits.

## Method

- Read the three tracks' NOTES.md in full and the packet's CLAUDE_HANDOFF.md, REPORT.md, MATHEMATICS.md, BOUNDARY.md,
  PORTALS.md, REHEATING_RESULTS.json, INDEPENDENT_REHEATING.json, REHEATING_VERIFICATION.json, VERIFICATION_RESULTS.json,
  BOUNDARY_RESULTS.json, PORTAL_RESULTS.json (git archive e160356 into the scratchpad; the head of code/reheat.c was read
  for the state-vector layout only; no Codex code was run).
- Recomputed with own code: (1) `rc_reheat.py` — the homogeneous Jordan-frame system written from MATHEMATICS.md in SciPy
  DOP853 with H evolved, own preparation solve (fsolve on the conditional stationarity with the actual R(0)), yield computed
  directly in GeV units, window diagnostics; (2) `rc_source.py` — momentum-resolved inverse-decay source (per-q LSODA system)
  and the MB moment closure; (3) `rc_boundary.py` — third build of the three-mode neutral model (neutral block enumerated,
  ladder-rule matrix elements); (4) `rc_arith.py` — dilution bookkeeping from Codex's recorded final states, T_RH, Higgs
  source, widths; (5) `rc_summary.py` — collects everything into `results/summary.json`.
- Reuse: none of Codex's or the three tracks' code. The Y free-mode diagnostic convention is re-implemented from its
  formula [Q hw_verify/Codex]. Thresholds, R_med, R_A, onset coordinate quoted from hw_critic, hw_gates, kin_source.

## Own numbers (results/rc_*.json, summary.json)

| Run (m_I/H_i, Γ_I/H_i, τ_end) | ε | N_X | N_YL | yield/Y_Q | D | notes |
|---|---|---|---|---|---|---|
| 2, 0.1, 250 | −0.198769 | −12263.84 | 3.5065e-3 | 0.55956 | 7.861 | Codex −0.1987691 / 0.5595552 / 3.5065e-3 |
| 2, 0.3, 83.3 | −0.218682 | −12041.88 | 3.5116e-3 | 0.78427 | 5.507 | Codex identical |
| 10, 0.1, 250 | −0.181528 | −3314.59 | 1.70713 | 0.32505 | 3.657 | Codex identical |
| 12, 0.1, 250 | −0.178595 | −3156.67 | 2.22325 | 0.31209 | 3.628 | resonance peak (rp_reheating 2.22) |
| 3, 0.5, 400 | −0.19929 | −8120.99 | 6.8156e-3 | 1.0408 | 2.798 | windows 40–50 / 64–80 / 160–200 / 320–400: 1.037 / 1.045 / 1.042 / 1.041 |
| 5, 1, 400 | −0.18399 | −5205.79 | 1.7768e-2 | 1.0414 | 1.793 | 1.036 / 1.047 / 1.042 / 1.041 |
| 10, 3, 400 | −0.17172 | −3670.25 | 0.76425 | 1.0136 | 1.299 | 1.005 / 1.019 / 1.013 / 1.014; Γt = 25 (τ ≤ 8.3): 1.070 |
| 3, 0.9, 80 | −0.20488 | −8176.63 | 6.6018e-3 | 1.1482 | 2.554 | Γ/m = 0.3 edge |

Preparation: own solve x = 76.39978265923669, y = −2.1115098541668 − 0.9393290405025 i, σ = 28.29314212586876
(Codex 76.39978266076902, −2.111509854283515, −0.9393290404896636, 28.293142125868524), R(0) = 12.000, residual 4e-12.
Charge timing: 50/90/99% of the late |N_X| at τ = 0.41/2.51/3.08 (m_I = 2) and 0.16/0.40/2.36 (m_I = 10), when radiation is
0.6–24% of ρ. Vector source: R_thermal momentum-resolved / closure = 717.9/724.1 (ε = 1e-9), 68902/69501 (1e-8),
968381/983804 (1e-7), 485009/502714 (3e-7), 162892/176058 (1e-6), 44326/46597 (1e-5); Y_pair identity to 1e-6; J ≥ 5Y_pair
everywhere. Boundary: dims 31/58/88; E(0, ±π/2, π) = 14.311195/14.084921/13.858647 (cutoff 2); T(+π/2) = −1.35626541,
−1.35764469, −1.35764502; 6 dE/dδ = T to 2e-9; ⟨cb²⟩ → 2√2 i; eigenstate torques 0. Higgs: Γ_hss 4.81568e-20 GeV,
Y_pair 2.23656e-8, f_rest 9.394, R_h 1942.8. T_RH = 0.5406 √(Γ_I M_r): 5.3e13 (Γ = 0.03 H_i) … 5.3e14 GeV (3 H_i).

## Limits

- Own checks: Gauss ≤ 7e-8 absolute (N ~ 1e4); relative Friedmann residual ≤ 1e-8 (τ_end ≤ 83), 7e-6 (τ = 250), 2e-5
  (τ = 400) from the late-time drift of the evolved H at atol 1e-13 while ρ falls by 1e-6–1e-8 — a diagnostic of this code;
  the physical outputs agree with the constraint-solved engine of rp_reheating and with Codex to ≤ 1e-9.
- ξR extremes are read on a 20001-point grid; peaks between samples are not resolved (range [−10, +5] H² holds analytically).
- Debris (R_central), thresholds, R_A, the onset coordinate and the clock/leakage/matching items are quoted, not recomputed.
- Same-session AI work; one supplied model; no observational statement; Candidate 1 unchanged; vacuum density underived.

Reproduction commands are at the end of CONSOLIDATED.md (≈ 12 min on 4 cores).
