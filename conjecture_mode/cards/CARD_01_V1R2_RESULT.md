# CARD #1 v1R2 — RESULT — LOWER-MINIMUM REPRODUCIBILITY AUDIT (POST-DATA NUMERICAL VERIFICATION)

## **TERMINAL: CARD-01-v1R2 — COMPUTATIONALLY UNRESOLVED AT AVAILABLE MINIMIZATION METHOD**

**Why.** Stage A failed on one of the four core minima: **fixed ε = −0.08 was not reproduced** under the frozen lower-minimum criterion within its maximum of 4 genuinely distinct starts. Under the hard stop:
- **no Stage B** was run (no continuation profile, no crossing verification);
- **no S5 state** is assigned;
- no extensions, control, MCMC, other optimizer or v1R3 were run.

**Status of the campaign series.** v1R2 was the final local-optimizer campaign. v1 (`826f0d1`) and v1R (`1c02f7d`) stay unresolved and unaltered. Card #1 remains **POSTULATED**.

| Item | Commit / value |
|---|---|
| Charter, frozen before any evaluation | `af647ac` |
| Starts | `code/v1r2_inputs/stageA_starts.json` |
| Driver log | Round 1 ran routes 2 and 3 for all four targets. Round 2 ran route 4 for ΛCDM and ε = −0.08, which were not yet reproduced. All 10 v1R2 runs exited rc = 0 and all are cycle-converged |
| Analyzer output | `code/card01_v1r2_summary.json` |

**Criterion** (charter §2). L = the lowest χ² over all executed starts (v1R arms plus v1R2 routes). A target is reproduced iff **≥ 2 cycle-converged starts have χ² ≤ L + 0.2**. Higher local minima are recorded but do not fail a target.

## 1. Stage-A record: every executed start (all cycle-converged)

### ΛCDM (ε = 0): **REPRODUCED**

| χ² | Start | H0 |
|---|---|---|
| **10976.084** | v1R2 route 2: best official ΛCDM sample (chain 1) | 68.11 |
| **10976.230** | v1R2 route 4: best sample of official ΛCDM chain 3 | 68.19 |
| 10976.390 | v1R2 route 3: best sample of official ΛCDM chain 2 | 68.07 |
| 10977.178 | v1R arm A (from v1 best) | 68.06 |
| 10978.132 | v1R arm B (random start) | 68.01 |

→ 2 within 0.2 of L = 10976.084. The ΛCDM minimum is **1.09 below** the v1R value.

### CPL w0waCDM: **REPRODUCED**

| χ² | Start | w0 | wa |
|---|---|---|---|
| **10963.805** | v1R2 route 3: best sample of official w0wa chain 1 | −0.393 | −1.806 |
| **10963.820** | v1R arm B: best sample of official w0wa chain 2 | −0.425 | −1.759 |
| 10964.061 | v1R2 route 2: best sample of official w0wa chain 4 | −0.402 | −1.820 |
| 10970.178 | v1R arm A (from v1 best), a higher local basin | −0.957 | −0.294 |

→ 2 within 0.2 of L = 10963.805.

### Card, fixed ε = −0.08: **NOT REPRODUCED** (the terminal failure)

| χ² | Start | H0 |
|---|---|---|
| **10972.769** | v1R arm A (from v1 best) | 69.27 |
| 10973.056 | v1R2 route 2: free-ε low solution fixed to −0.08 | 69.22 |
| 10973.376 | v1R2 route 3: ε = −0.06 solution continued to −0.08 | 69.38 |
| 10974.039 | v1R2 route 4: ε = −0.10 solution continued to −0.08 | 69.11 |
| 10983.994 | v1R arm B (random start), a higher local basin | 68.99 |

→ only 1 start within 0.2 of L = 10972.769. The closest independent route is **0.287** above it. All 4 genuinely distinct starts are used up (route 1 plus routes 2–4).

### Card, free ε: **REPRODUCED**

| χ² | Start | ε |
|---|---|---|
| **10972.720** | v1R2 route 2: ε = −0.08 solution, ε released | −0.0835 |
| **10972.906** | v1R arm A (from v1 best free) | −0.0975 |
| 10973.110 | v1R arm B (start ε = −0.20) | −0.0861 |
| 10973.660 | v1R2 route 3: ε = −0.10 solution released | −0.1007 |

→ 2 within 0.2 of L = 10972.720.

## 2. Numerical facts established (no S5 state; orientation only)

| Quantity | Reproduced lowest basin |
|---|---|
| χ²_Λ | **10976.08** |
| χ²_CPL | **10963.81** at (w0, wa) ≈ (−0.39, −1.81) |
| χ²_card, free-ε fit | **10972.72** at ε ≈ −0.084 |
| χ²_card, fixed ε = −0.08 | not reproduced (lowest 10972.77; next distinct route 10973.06) |

**Orientation, explicitly NOT an S5 evaluation.** Because Stage A failed, the profile set was never resolved and `decide()` is not applied. Combining the three reproduced core minima gives:

| Quantity | Value |
|---|---|
| I_card = χ²_Λ − χ²_card(free) | 10976.08 − 10972.72 = **3.36** |
| I_CPL = χ²_Λ − χ²_CPL | 10976.08 − 10963.81 = **12.28** |
| F | ≈ **0.27** |

Relative to v1R's best-found values (I_card 4.41, I_CPL 13.36, F 0.33):
- the **ΛCDM minimum dropped by 1.09** once it was started from the official DESI ΛCDM samples, which lowered both I_card and I_CPL;
- the Card improvement over ΛCDM is now **below the S3 cut of 4**;
- the Card recovers about a quarter of the CPL improvement.

In every campaign, the CPL comparator remains the far stronger fit.

## 3. Execution note

The container restarted at about 15:56 UTC, **after** the driver had already printed its terminal and "STOP FOR OWNER REVIEW" (log complete; all 10 result files present). No run was affected and nothing was rerun.

## 4. Firewall

- This is a **POST-DATA NUMERICAL VERIFICATION RESULT — NOT PREREGISTERED v1**.
- No S5 state is assigned. Card #1 is **not** EXPLANATORY-SURVIVES, EXPLANATORY-KILLED or MODEL-KILLED by any accepted run.
- Card #1 remains POSTULATED, with no PREDICTED or VALIDATED status.
- Any further work needs a new owner ruling after review of this terminal, e.g. a fundamentally different, independent likelihood / minimization implementation.
