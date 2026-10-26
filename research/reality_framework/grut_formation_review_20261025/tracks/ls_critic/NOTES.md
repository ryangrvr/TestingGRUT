# Track ls_critic: notes on the cross-track recomputations (Codex b448f25)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. The consolidated verdicts are in
[CONSOLIDATED.md](CONSOLIDATED.md); this file records only this track's own runs. Nothing was committed. No Codex
code was run.

**Scope.** One supplied tree-level working model. This is not a validation or an exclusion.

All code below is ours. It imports `ls_threshold/code/lsrelic.py`, `sp_yledger/code/ymodel.py` (yrelease.c) and
`ls_ledger/code/ls_project.py` unchanged, with no bytecode written there.

## Runs

| Script | Question | Result |
|---|---|---|
| `heavy_threshold.py` | At the ls_scaling heavy point (ε = 0.14571744, normalised yield): ke f(R), ke + s, μ_c = 0.923, and the proxy R_L by bisection | Proxy R_L = 81.875–82.1875 (the table gave 87.9). ke at μ_c = 0.861 erases for R ≤ 5 and gives f = 0.045 at R = 8.9, 0.383 at 15 and 1.36 at 34.2. ke + s gives 8e-187 at 8.9, 4.6e-14 at 15 and 0.511 at 34.2. ke at μ_c = 0.923 erases at R ≤ 3 and gives 0.249 at 8.9 [C] |
| `ke_refine.py` | Edge of the ke crossing and N control | μ_c = 0.861: f = 4e-162 (R = 6), 1.6e-26 (7), 0.0057 (8), 0.045 (8.9). At N = 1000: 1.7e-173 (5), 1.4e-36 (7), 0.041 (8.9). ke + s: 0.093 (20), 0.241 (25). μ_c = 0.923: 3.7e-4 at R = 5 [C] |
| `relaxed_check.py` | Re-check of ls_scaling's Y-relaxed preparation, with our own stationarity solve, the sp_yledger release to a = 1000 and our own closed-form ledger | Candidate 10: ε = 0.164869, R_lower = 148.855, R_central = 178.953, yield 5.970 Y_Q. Candidate 17: 0.157832, 192.685, 224.613, 5.701. The ls_scaling values are 148.92 / 192.76. The window and snapshot projections agree to 6e-4 [C] |
| `stack_check.py` | Stacked favourable choices in the proxy at the Y-relaxed candidates | Candidate 10: R_L = 110.0–110.3 (×1), 158.0–158.4 (σv ×72) and 158.9–159.4 (×72 + s). At R_lower, f = 6.58 (×1) and 7.8e-7 (×72 + s, passes); at R_central, f = 2.26 (×72 + s). Candidate 17 (×72 + s): R_L = 149.5–150.0; f = 5.23 at R_lower [C] |
| `stack_ke.py` | That passing corner in ke (×72, s on) | f = 2.64 (μ_c = 0.861) and 3.11 (0.923) [C] |

**Derived from these runs [C, E].**
- **ke crossing at the heavy point.** Log-interpolated, it is R ≈ 8.2 (μ_c = 0.861) and 7.0 (0.923). With the s route
  it lies between 15 and 20.
- **Vector headroom at the base point.** Defined as (R_L − 1.1445)/(19.642 × 2788.18) vectors per net charge:

  | Threshold | Headroom |
  |---|---:|
  | Table | 1.58e-3 |
  | Direct proxy | 1.48e-3 |
  | ke | 1.1e-4 to 1.3e-4 |
  | ke + s | 2.5e-4 to 3.4e-4 |

- **Window edge in m = M_Y/H** [E]. Interpolated with R ∝ m^−4.7 between m = 15 and 20, it is about 20–21 under ke
  and about 17–18 under ke + s. With R ∝ m^−3.6 between m = 10 and 15, it is about 12 under the proxy.

## Limits

- **ke is an equilibrium proxy.** It assumes instant condensation and hidden kinetic equilibrium.
- **Grid.** The cumulative energy-work residual is 1–2% at N = 500.
- **Root fallbacks.** Runs with about 200 fallbacks sit at the p floor (1e-250) and mean only "erased".
- **ε at the heavy point.** It is held at the base-point value for every R. The m = 20 point actually has ε = −0.142.
- **The ×72 conversion factor** is ls_threshold's sensitivity value, not a derived rate.

## Reproduction

See CONSOLIDATED.md, "Reproduction". Large logs are in the scratchpad under `review/ls_critic/`.
