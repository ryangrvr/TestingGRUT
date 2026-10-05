# CARD #1 v1R — RESULT (PRIMARY) — POST-DATA VERIFICATION-REPAIR

## **CARD-01-v1R — STILL NUMERICALLY UNRESOLVED**

**Mechanical S5 state: NOT ASSIGNED.** Under the v1R charter, the frozen `decide()` is applied only if every load-bearing target passes the convergence check. Five of seven do not.

**This is a post-data verification-repair result, not the preregistered v1 result.** v1 remains as recorded at `826f0d1`: **CARD-01-v1 — NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT**. Card #1 remains **POSTULATED**.

| Item | Commit / value |
|---|---|
| v1R protocol, frozen before execution | `60e5a18` |
| Execution repair ER-G (per-cycle subprocesses after OOM) | `bca8f8c` |
| Analyzer output | `code/card01_v1r_primary_summary.json` (`code/card01_v1r_analyze.py`) |
| Procedure (unchanged from the charter) | BOBYQA `rhoend` 0.005; official-DESI-chain proposal covariances; cycles until a full cycle improves χ² by < 0.1 (cap 6); independent arm B for every load-bearing target; acceptance `|A − B| ≤ 0.2` |
| Run | 22 arm-A targets (all cycle-converged) and 7 arm-B targets (all cycle-converged). The load-bearing set was stable at the second identification. The driver stopped at "STOP FOR OWNER REVIEW" |

## 1. Convergence table (χ²_eff per refinement cycle; LB = load-bearing)

| Target | LB | Arm A cycles | Arm B cycles | Best | Check |
|---|---|---|---|---|---|
| ε = 0 (ΛCDM) | **LB** | 10977.20 → 10977.18 | 10979.04 → 10978.59 → 10978.17 → 10978.13 | 10977.18 | **FAIL**, \|A−B\| = 0.955 |
| ε = −0.005 | **LB** | 10976.07 → 10976.00 | 10977.44 → 10977.35 | 10976.00 | **FAIL**, 1.352 |
| ε = −0.010 | | 10975.69 → 10975.65 | — | 10975.65 | (not LB) |
| ε = −0.015 | | 10975.77 → 10975.76 | — | 10975.76 | (not LB) |
| ε = −0.020 | | 10974.96 → 10974.95 | — | 10974.95 | (not LB) |
| ε = −0.030 | | 10974.21 → 10974.20 | — | 10974.20 | (not LB) |
| ε = −0.040 | | 10974.91 → 10974.87 | — | 10974.87 | (not LB) |
| ε = −0.050 | | 10974.95 → 10973.69 → 10973.62 | — | 10973.62 | (not LB) |
| ε = −0.060 | | 10973.64 → 10973.57 | — | 10973.57 | (not LB) |
| ε = −0.080 | **LB** | 10972.93 → 10972.79 → 10972.77 | 10985.44 → 10984.38 → 10984.05 → 10983.99 | 10972.77 | **FAIL**, 11.225 |
| ε = −0.100 | | 10973.70 → 10973.68 | — | 10973.68 | (not LB) |
| ε = −0.120 | | 10973.19 → 10973.13 | — | 10973.13 | (not LB) |
| ε = −0.150 | **LB** | 10974.32 → 10974.30 | 10974.31 → 10974.27 | 10974.27 | **PASS**, 0.036 |
| ε = −0.200 | **LB** | 10977.05 → 10976.94 → 10976.88 | 10977.64 → 10976.87 → 10976.84 | 10976.84 | **PASS**, 0.041 |
| ε = −0.250 | | 10980.00 → 10979.99 | — | 10979.99 | (not LB) |
| ε = −0.300 | | 10985.32 → 10985.30 | — | 10985.30 | (not LB) |
| ε = −0.400 | | 10997.51 → 10997.41 | — | 10997.41 | (not LB) |
| ε = −0.500 | | 11013.95 → 11013.68 → 11013.68 | — | 11013.68 | (not LB) |
| ε = −0.700 | | 11050.20 → 11050.19 | — | 11050.19 | (not LB) |
| ε = −1.000 | | 11113.91 → 11113.46 → 11113.25 → 11113.00 → 11112.94 | — | 11112.94 | (not LB) |
| free-ε | **LB** | 10973.08 → 10972.93 → 10972.91 (ε = −0.0975) | 10973.18 → 10973.11 | 10972.91 | **FAIL**, 0.204 |
| CPL w0waCDM | **LB** | 10970.19 → 10970.18 (w0 −0.96, wa −0.29) | 10963.84 → 10963.82 (w0 −0.425, wa −1.759) | **10963.82** | **FAIL**, 6.358 |

**Load-bearing failures:**

| Target | \|A−B\| | Which arm is lower |
|---|---|---|
| ε = −0.08 | 11.225 | A |
| CPL | 6.358 | **B** |
| ε = −0.005 | 1.352 | A |
| ε = 0 | 0.955 | A |
| free-ε | 0.204 | A |

**Passes:** ε = −0.15 and ε = −0.20, the lower 95 % crossing bracket.

## 2. Best values found (best-of-arms; **not resolved** minima)

| Quantity | Value |
|---|---|
| χ²_Λ (ε = 0) | 10977.178 (arms disagree by 0.96) |
| χ²_card, branch best | 10972.769 at ε̂ = −0.080 (grid). Free-ε: 10972.91 at ε = −0.0975 |
| χ²_w0wa | **10963.820** at (w0, wa) = (−0.425, −1.759), from arm B |
| I_card = χ²_Λ − χ²_card | 4.41 |
| I_CPL = χ²_Λ − χ²_w0wa | 13.36 |
| F = I_card / I_CPL | 0.33 |
| Descriptive 95 % profile set | [−0.196, −0.002] |

These numbers are reported because the owner ruling requires them. Under the charter they carry **no S5 state**, because the load-bearing minima are not numerically resolved.

## 3. What v1R did establish (numerical facts only)

1. **The v1 CPL comparator was badly under-converged.** Started from the official DESI highest-likelihood w0wa sample, arm B reached χ² = 10963.82 at (w0, wa) = (−0.425, −1.759), **6.36 below** arm A, which started from the v1 best (10970.18 at (−0.96, −0.29)).
   - The arm-B minimum sits at the official DESI posterior location for this combination (D1 mean ≈ (−0.42, −1.75)).
   - v1's I_CPL (7.05) and F (0.61) were therefore distorted by comparator non-convergence.
   - CPL is still formally unresolved, because no second independent start has reached 10963.8.
2. **Arm-A cycle convergence does not guarantee the global minimum.** Every arm met the cycle criterion (improvement < 0.1), yet independent starts landed up to 11 χ² apart, at ε = −0.08 where arm B stalled near 10984.
   - In four of the five failures the v1-seeded arm A is the lower one; for CPL it is arm B.
   - The likelihood surface traps BOBYQA in distinct local stalls even at `rhoend` 0.005.
3. **The arm-A profile is smoother than v1's but not monotone.** Residual jitter of up to ~0.9 χ² between neighbouring points remains (e.g. ε = −0.03 → −0.04: +0.67; ε = −0.10 → −0.12: −0.55). That is consistent with the remaining minimizer stalls. Its lowest values lie near ε ≈ −0.08 to −0.12. Under §2 this profile still carries no S5 state.

## 4. Not run (owner ruling)

- the SN extensions;
- the ε > 0 control;
- the MCMC posterior.

No further repair was attempted. Any additional numerical protocol (e.g. more independent starts per load-bearing target, or a different optimizer) would need a new owner decision, frozen before it runs.

## Provenance

| Item | Status |
|---|---|
| v1 raw outputs, archived mis-seeded runs, `826f0d1`, the v1 run configuration and the owner threshold ruling | unchanged |
| v1R files | all additive |
| v1R authorization | given only after the v1 convergence failure and preliminary values were observed (see `CARD_01_V1R_CHARTER.md`) |
| `vsctestinggrut`, `grut-backreaction-identifiability-0`, frozen scientific and governance branches | not touched |
