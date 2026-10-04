# CARD #1 — RESULT (v1, PRIMARY)

## **CARD-01-v1 — NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT**

**Mechanical S5 output: `EXPLANATORY-SURVIVES`.**
**MECHANICAL OUTPUT — NOT SCIENTIFICALLY ACCEPTED DUE TO FROZEN CONVERGENCE FAILURE.**

This is **not** EXPLANATORY-SURVIVES, EXPLANATORY-KILLED or MODEL-KILLED. Its basis is the frozen `CHI2_START_DISAGREEMENT_FLAG = 0.2` (run configuration `c9899abc`), applied under the owner convergence ruling (`CARD_01_OWNER_CONVERGENCE_RULING.md`).

| Item | Commit / value |
|---|---|
| Pre-data spec | `8c634145` |
| Owner threshold ruling | `86ccbb2e` |
| Run configuration | `c9899abc` |
| Pre-evaluation corrections CR-1 to CR-4 | `b50f5224` |
| Convergence ruling and execution repairs | `baa6c053` |
| Primary combination | DESI DR2 BAO + DESI DR2 baseline CMB (Planck PR3 low-ℓ TT Commander and EE SimAll via clik; Planck PR4 NPIPE CamSpec TTTEEE; Planck + ACT DR6 lensing v1.2) |
| Minimizer (frozen) | BOBYQA (Cobaya default `rhoend` 0.05), 3 seeded starts per point; 66/66 fits completed |
| Analyzer output | `code/card01_d3_primary_summary.json` (`code/card01_d3_analyze.py`) |

*History:* this file previously recorded CARD-01-ACCESS-BLOCKED (commit `ae07fa1c`). That state was superseded when the owner enabled full network access and the card resumed under the unchanged ruling.

---

## 1. Frozen convergence diagnostic

Starts are reported as computed: not averaged, not refined. χ²_eff = −2 ln L + nuisance Gaussian-prior terms, with flat-prior constants removed (CR-3).

| ε | start χ² (sorted) | best | spread | > 0.2? | H0 | ω_c |
|---|---|---|---|---|---|---|
| −1.000 | 11114.44, 11116.36, 11119.45 | 11114.44 | 5.01 | **YES** | 81.86 | 0.1243 |
| −0.700 | 11051.69, 11051.93, 11053.88 | 11051.69 | 2.19 | **YES** | 77.89 | 0.1228 |
| −0.500 | 11014.15, 11014.16, 11018.13 | 11014.15 | 3.98 | **YES** | 75.19 | 0.1216 |
| −0.400 | 10997.59, 10998.81, 11001.07 | 10997.59 | 3.48 | **YES** | 73.71 | 0.1211 |
| −0.300 | 10986.20, 10986.29, 10986.46 | 10986.20 | 0.26 | **YES** | 72.33 | 0.1204 |
| −0.250 | 10980.08, 10980.78, 10981.61 | 10980.08 | 1.53 | **YES** | 71.58 | 0.1201 |
| −0.200 | 10977.25, 10977.71, 10977.80 | 10977.25 | 0.56 | **YES** | 70.91 | 0.1197 |
| −0.150 | 10974.47, 10974.70, 10975.04 | 10974.47 | 0.57 | **YES** | 70.19 | 0.1193 |
| −0.120 | 10973.48, 10973.91, 10975.17 | 10973.48 | 1.68 | **YES** | 69.77 | 0.1190 |
| −0.100 | 10973.83, 10974.33, 10975.77 | 10973.83 | 1.94 | **YES** | 69.50 | 0.1188 |
| **−0.080** | 10972.99, 10973.42, 10974.61 | **10972.99** | 1.62 | **YES** | 69.28 | 0.1185 |
| −0.060 | 10973.70, 10973.73, 10977.79 | 10973.70 | 4.09 | **YES** | 69.01 | 0.1183 |
| −0.050 | 10975.31, 10976.83, 10977.64 | 10975.31 | 2.33 | **YES** | 68.69 | 0.1185 |
| −0.040 | 10975.17, 10975.70, 10978.42 | 10975.17 | 3.25 | **YES** | 68.62 | 0.1182 |
| −0.030 | 10974.34, 10974.94, 10978.56 | 10974.34 | 4.22 | **YES** | 68.47 | 0.1182 |
| −0.020 | 10975.27, 10975.59, 10976.02 | 10975.27 | 0.75 | **YES** | 68.36 | 0.1181 |
| −0.015 | 10977.19, 10977.51, 10978.48 | 10977.19 | 1.29 | **YES** | 68.30 | 0.1180 |
| −0.010 | 10976.57, 10977.48, 10979.85 | 10976.57 | 3.28 | **YES** | 68.29 | 0.1178 |
| −0.005 | 10976.47, 10976.85, 10976.98 | 10976.47 | 0.51 | **YES** | 68.21 | 0.1177 |
| **0.000 (ΛCDM)** | 10977.28, 10977.47, 10977.79 | **10977.28** | 0.51 | **YES** | 68.06 | 0.1179 |

| Other load-bearing fit | start χ² (sorted) | best | spread | > 0.2? |
|---|---|---|---|---|
| free-ε, seeded at best grid ε = −0.080 | 10973.32 (ε = −0.0965), 10974.14 (ε = −0.0768), 10976.08 (ε = −0.0589) | 10973.32 | 2.77 | **YES** |
| CPL w0waCDM comparator | 10970.23, 10970.42, 10973.72 | 10970.23 (w0 = −0.957, wa = −0.292) | 3.49 | **YES** |

**All 20 grid points, the free-ε fit and the CPL comparator exceed the frozen 0.2 flag.** So does every load-bearing quantity named in the owner ruling:

| Load-bearing quantity | Spread |
|---|---|
| ΛCDM null | 0.51 |
| CPL comparator | 3.49 |
| Branch best grid point (ε = −0.08) | 1.62 |
| Free-ε fit | 2.77 |
| Profile-set brackets (ε = −0.20/−0.15 and −0.005/0) | 0.56 / 0.57 and 0.51 / 0.51 |

The best-of-3 profile is visibly non-smooth: neighbouring grid points jump by 1–2 in χ² (e.g. ε = −0.030 → −0.020 → −0.015). That is the signature of under-converged minima, not of the likelihood.

## 2. Mechanical S5 output (frozen `decide()`, applied exactly as committed)

| Quantity | Value |
|---|---|
| χ²_Λ (ε = 0, best of 3) | 10977.277 |
| χ²_card, branch best | 10972.995 at ε̂ = −0.080 (grid; the free-ε fit reached 10973.317 at ε = −0.0965) |
| χ²_w0wa (best of 3) | 10970.226 |
| **I_CPL** = χ²_Λ − χ²_w0wa | **7.05** (≥ 4, so a targeted signal exists) |
| **I_card** = χ²_Λ − χ²_card | **4.28** (≥ 4 by **0.28**) |
| **F** = I_card / I_CPL | **0.607** (≥ 0.5) |
| ε̂ | −0.080 (< −0.03; not at the ε = −1 boundary) |
| Descriptive 95 % profile set (Δχ² ≤ 3.84 from branch best) | [−0.193, −0.003] |
| Boundary-null Chernoff p (secondary) | 0.019 |
| Branch vs CPL tension, χ²_card − χ²_w0wa | 2.77 |

→ `decide()` returns **EXPLANATORY-SURVIVES**.

**Why it is not accepted.** The margin over the S3 cut is **0.28** in χ². The start-to-start spread on the two quantities defining I_card is 0.51 (ΛCDM) and 1.62 (branch best), and the spread on the CPL comparator that defines F is 3.49. Re-picking among the computed starts alone moves I_card across 4 in either direction:
- ΛCDM start 10977.79 vs branch start 10974.61 gives **3.18**;
- ΛCDM start 10977.28 vs free-ε best 10973.32 gives **3.96**.

**F's denominator is also not converged.** The best CPL point (w0 = −0.957, wa = −0.292) lies far from the official DESI DR2 posterior mean for this same combination (w0 ≈ −0.42, wa ≈ −1.75, per D1; diagnostic only). This is consistent with the minimizer stopping in a flat region of the w0–wa degeneracy rather than at the CPL minimum. If so, I_CPL is understated and F overstated.

Under the frozen 0.2 flag, the S5 state is therefore **not resolved** by this run.

## 3. D1/D2 diagnostics (official DESI chains; DIAGNOSTIC — NOT THE KILL TEST)

From `code/card01_d1d2_diagnostic.json`; all 24 official chain files verified against DESI's published sha256sum.

| Combination | w0wa posterior mean | P(w0 > −1, wa < 0) | Card #1 primary ray (ε < 0) reaches HPD level | Inside 95 %? |
|---|---|---|---|---|
| DESI + CMB (primary) | (−0.42, −1.75) | 0.998 | 0.991 | no |
| + Pantheon+ | (−0.84, −0.62) | 0.999 | 0.997 | no |
| + Union3 | (−0.67, −1.09) | 1.000 | 1.000 | no |
| + DESY5 | (−0.75, −0.86) | 1.000 | 1.000 | no |

These diagnostics cannot alter the D3 state (S8).

## 4. Not run (owner ruling §6)

- the SN extensions (Pantheon+, Union3, DESY5);
- the ε > 0 control;
- the secondary MCMC posterior.

The driver stopped at "STOP FOR OWNER REVIEW".

## 5. Execution events during the run

See `CARD_01_EXECUTION_REPAIRS.md`:
- ER-A: OOM;
- ER-B: mis-seeded free-ε fits, archived and excluded;
- ER-C: container restart;
- ER-F: a second container restart at ≈ 23:12 UTC killed the last in-flight fit (free-ε start 2). It was rerun under the unchanged protocol from the same best-grid seed (−0.080). All other 65 fits were already complete.

No minimizer setting, likelihood, model, threshold or result logic was changed.

## 6. What is open for the owner (no action taken)

- **CR-5** remains a **POST-DATA REPAIR PROPOSAL — NOT PART OF CARD #1 v1 — NOT YET AUTHORIZED**.
- Any repaired computational analysis would need its convergence protocol frozen before it runs. It would be labeled post-data / verification-repair, not the preregistered v1 test.
- The mechanical numbers above (I_card ≈ 4.3, F ≈ 0.6) **must not be cited as a result**. They sit inside the minimizer's own scatter.

## Firewall

- Card #1 remains **POSTULATED**.
- Card #1 v1 has **no accepted scientific verdict**.
- Even a future accepted EXPLANATORY-SURVIVES would not be support for GRUT. ε = 0 is ΛCDM, and the comparison is retrospective (DR2 already exists).
- The historical C2-R reconnaissance grade is unchanged.
