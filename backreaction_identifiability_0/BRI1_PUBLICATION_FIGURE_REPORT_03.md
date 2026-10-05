# BRI1 PUBLICATION FIGURE REPORT 03 — V3-R2
## CORRECTED TENSOR QUADRATURE (remediation of V3-R1 broadcasting defect)
## Executed under `BRI1_PUBLICATION_VERIFICATION_CHARTER_01.md` (frozen `8d765ea`), V3 lane only
## Updated after post-execution reporting/provenance audit and remote verification

**Status:** V3-R2 — illustrative / evidence-grade numerical demonstration. **Not part of the analytic proof.** V1 (the proof), V2 (KNOWN-RESULT-NEW-FRAMING), and the original V3 / V3-R1 reports are unchanged.

**Authoritative reproduction script:** `publication_verification/v3_r2/v3_r2_authoritative.py`
**Authoritative machine-readable output:** `publication_verification/v3_r2/v3_r2_results.json`
**CSV table:** `publication_verification/v3_r2/v3_r2_table.csv` (generated from JSON)
**Figures:** `publication_verification/v3_r2/v3_r2_figures.png` (generated from JSON)

**Independent code-path numerical reproduction:** a separately written implementation (DOP853, rtol = 1e-12, 160-point tensor GL grid on [-6,6]^2) reproduced the repository values to all printed digits for checked cases. Recorded below. Not external peer review, not independent external verification, not independent human replication.

---

## 1. Frozen scientific choices (identical to V3-R1; no retuning)

- Primary time: **t★ = 0.5**
- Diagnostic times: **t ∈ {0.25, 0.75, 1.0}**
- Bath-size grid: **N_B ∈ {4, 8, 16, 32, 64, 128}**
- Quadrature domain: **[-6, 6]^2** (correct mapping via `xs = L * xg`)
- Quadrature resolutions: **nx ∈ {120, 240, 480}**
- RK4 refinement: **dt ∈ {10^-3, 5×10^-4, 2.5×10^-4}**
- Protocols: **P0 and P1** (frozen X1)
- Scaling fit: OLS on **all six** N_B values
- The **only** permitted change: correcting the implementation.

## 2. Corrected implementation

Explicit 2-D tensor grids with **every `(x_i, p_j)` pair evolved independently** (see `v3_r2_authoritative.py`). Assertions enforce correct shapes.

## 3. Mandatory P0 stationarity firewall

At the reference resolution (nx=240, dt=5×10^-4):

| t | E[x(t)] | Var[x(t)] | κ₃[x(t)] |
|---|---|---|---|
| 0.0 | 0.000e+00 | 0.4679199170 | 5.551e-17 |
| 0.25 | 5.551e-17 | 0.4679199128 | -7.792e-17 |
| 0.5 | 0.000e+00 | 0.4679199051 | 0.000e+00 |
| 0.75 | 0.000e+00 | 0.4679199058 | 0.000e+00 |
| 1.0 | -1.388e-17 | 0.4679199137 | 2.642e-17 |

**All P0 stationarity conditions satisfied.** Variance matches the independent reference to relative discrepancy 2.373e-16. **The primary implementation firewall passes.**

### Exact virial identity (Gibbs)

**⟨x²⟩ + ⟨x⁴⟩ = 1: residual = 6.661e-16** (machine epsilon level).

## 4. Independent x-marginal check (separate 1-D pipeline)

| Moment | 1-D result | Reference | Agreement |
|---|---|---|---|
| ⟨x²⟩ | 0.467919916973665 | 0.467919916974 | ✓ |
| ⟨x⁴⟩ | 0.532080083026335 | 0.532080083026 | ✓ |
| ⟨x⁶⟩ | 0.871679667894661 | 0.871679667895 | ✓ |
| ⟨x²⟩+⟨x⁴⟩ | 1.000000000000000 | 1 (exact) | ✓ |

## 5. P1 finite-ε results (primary witness t★ = 0.5)

| N_B | ε | κ₃(X^ε) | Var(X^ε) | κ₃(F) | γ₁(F) | N_B·γ₁ |
|---|---|---|---|---|---|---|
| 4 | 0.5000 | -3.417347366832e-06 | 0.4679199048 | -1.708673683416e-06 | -5.338286160021e-06 | -2.135314464008e-05 |
| 8 | 0.3536 | -2.416429497459e-06 | 0.4679199050 | -8.543368419563e-07 | -2.669143079768e-06 | -2.135314463814e-05 |
| 16 | 0.2500 | -1.708673684248e-06 | 0.4679199050 | -4.271684210620e-07 | -1.334571539891e-06 | -2.135314463826e-05 |
| 32 | 0.1768 | -1.208214749028e-06 | 0.4679199051 | -2.135842105419e-07 | -6.672857699159e-07 | -2.135314463731e-05 |
| 64 | 0.1250 | -8.543368422763e-07 | 0.4679199051 | -1.067921052845e-07 | -3.336428849845e-07 | -2.135314463901e-05 |
| 128 | 0.0884 | -6.041073745303e-07 | 0.4679199051 | -5.339605263690e-08 | -1.668214424715e-07 | -2.135314463635e-05 |

## 6. Diagnostic times (all reported)

| t | γ₁(F) at N_B=4 | γ₁(F) at N_B=128 | N_B·γ₁ (N_B=4) |
|---|---|---|---|
| 0.25 | -4.74716774e-08 | -1.48348995e-09 | -1.89886710e-07 |
| 0.5 | -5.33828616e-06 | -1.66821442e-07 | -2.13531446e-05 |
| 0.75 | -7.66741491e-05 | -2.39606714e-06 | -3.06696596e-04 |
| 1.0 | -4.64190824e-04 | -1.45059613e-05 | -1.85676330e-03 |

All four times: nonzero P1 skewness, clean N_B^-1 scaling, stable N_B·γ₁.

## 7. Odd-ε check (implementation/symmetry control)

| N_B | κ₃(+ε) | κ₃(-ε) | sum | ratio |
|---|---|---|---|---|
| 4 | -3.417347366832e-06 | 3.417347366754e-06 | -7.78999e-17 | -0.999999999977205 |
| 16 | -1.708673684248e-06 | 1.708673684192e-06 | -5.55383e-17 | -0.999999999967496 |

κ₃(-ε) = -κ₃(+ε) to 10 significant figures. Implementation/symmetry control — **not** independent confirmation of V1 physics.

## 8. Scaling fit

**Global fitted p = 1.0000000000328206**

Local effective exponents: ['1.0000000001313', '0.9999999999921', '1.0000000000641', '0.9999999998851', '1.0000000001795']

Deviations from p=1 are far below displayed physical significance. Reported as: **p ≈ 1.000000**.

## 9. Finite-N corrections (explicitly reported)

At t=0.5, N_B·γ₁ varies across the grid by relative amount ~1.7×10^-10 — no physically meaningful correction is resolved.

At t=1.0:

| N_B | N_B·γ₁ |
|---|---|
| 4 | -1.8567632957546506e-03 |
| 8 | -1.8567631648600351e-03 |
| 16 | -1.8567630994129301e-03 |
| 32 | -1.8567630666878328e-03 |
| 64 | -1.8567630503264713e-03 |
| 128 | -1.8567630421445214e-03 |

Relative drift (t=1.0, N_B 4→128): 1.748e-10. A genuine monotone correction is visible.

At t=1.0, κ₃(X^ε)/ε shows somewhat larger relative variation (~10^-6) because the simultaneous variance correction partly cancels it in standardized skewness.

> The small finite-N corrections are consistent with the weak smoothstep forcing and near-linear-response regime.

This is an interpretation, not a proved cause.

## 10. Independent code-path numerical reproduction

A separately written implementation (DOP853, rtol = 1e-12, 160-point tensor GL grid on [-6,6]^2) reproduced the following repository values to all printed digits:

- t=0.5: N_B=4 γ₁ = -5.338286e-06; N_B=16 γ₁ = -1.334572e-06
- t=1.0: N_B=4 γ₁ = -4.641908e-04; N_B=16 γ₁ = -1.160477e-04
- P0 Var = 0.4679199

Recorded as: **independent code-path numerical reproduction**. Not: external peer review, independent external verification, independent human replication.

## 11. Post-execution reporting audit

The following provenance defects were discovered in sequence:

1. V3-R1 tensor broadcasting defect (pairwise instead of tensor evolution).
2. First R2 report populated cumulant columns incorrectly.
3. Odd-ε table inherited the same packaging error.
4. Gibbs residual was misreported as 8.28×10^-8.
5. Prose was written around a number not generated by the operative pipeline.
6. An abandoned [0,6] draft remained in the reproduction directory.
7. Publication-facing values were manually hard-coded in `save_results.py`.
8. **Commit `ea89f68` described a corrected Report 03 and removal of hard-coded publication data, but remote verification showed that those changes were absent from the pushed commit. The commit message therefore overstated the contents of the commit.**

The earlier residual (8.28×10^-8) and its explanation were not produced by the operative numerical pipeline. Both were removed during provenance repair.

**None of these reporting/provenance failures has survived into the authoritative result pipeline. The corrected tensor computation itself agrees with a separate independent code-path reproduction.**

**Adopted rules:**

> A task is not complete when the working tree looks correct. Completion requires verification of the pushed remote commit itself.

> Prose may interpret only values present in the authoritative machine-readable output or an explicitly identified independent check.

> Resolution convergence demonstrates stability of the implemented computation, not correctness of the implementation. Every numerical illustration must include at least one exact or independently known control that is not merely the same symmetry as the target observable.

## 12. Interpretation

1. P0 Gibbs stationarity: **Yes**. 2. Independent marginal check: **Yes**. 3. P1 nonzero at t★=0.5: **Yes**. 4. Sign agrees with V1: **Yes**. 5. Witness decreases with N_B: **Yes**. 6. Scaling compatible with 1/N_B: **Yes** (p ≈ 1.000000). 7. Finite-N corrections visible: **Yes** (~10^-10 at t=0.5, ~10^-7 at t=1.0). 8. Odd-ε behavior: **Yes** (to 10 significant figures; implementation control). 9. Quadrature/integration errors below effect: **Yes**. 10. Anything contradicts V1: **No**.

## 13. Comparison with invalid V3-R1

| Property | V3-R1 (invalid) | V3-R2 (corrected) | Different? |
|---|---|---|---|
| P0 Var(t★) | 0.96005 (wrong) | 0.46792 (correct) | **Yes** |
| P0 stationarity | Failed | Passed | **Yes** |
| Sign of γ₁(P1) | Negative | Negative | No |
| Fitted p | 1.0000 | 1.000000 | Similar |
| Exceeds Gibbs bound? | Yes | No | **Yes** |

## 14. Disposition

# **V3-R2-ILLUSTRATION-PASS**

All mandatory controls pass on the authoritative pipeline. Independent code-path reproduction confirms the values. Evidence-grade only.