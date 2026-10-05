# CARD #1 v1R — POST-DATA CONVERGENCE VERIFICATION — CHARTER and FROZEN NUMERICAL PROTOCOL

**Label.** Every v1R output is **POST-DATA / VERIFICATION-REPAIR**. It is **never** the preregistered v1 result.

**Authorization.**
- The owner authorized v1R (owner ruling, 2026-10-05) **only after observing**:
  - the v1 convergence failure;
  - the preliminary v1 values (mechanical I_card = 4.28, F = 0.61).
- **v1 is unchanged and permanently recorded** at `826f0d1` as **CARD-01-v1 — NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT**. Its mechanical `EXPLANATORY-SURVIVES` is recorded but **not scientifically accepted**.
- v1R does not overwrite or retroactively repair v1. All v1R files are additive.

**Sole purpose.** Determine the actual PRIMARY likelihood profile and minima closely enough to know what the unresolved v1 computation was trying to estimate.

**Commit order.** This charter and the files listed below are committed **before any v1R optimization is executed**.

## 1. Inherited unchanged from v1 (frozen)

All of the following are inherited exactly from the frozen v1 configuration (`card01_run_config.py` at `c9899abc`, plus the pre-evaluation corrections CR-1 to CR-4 at `b50f5224`) and the owner threshold ruling (`86ccbb2e`):

- **Model:** the Card #1 physical postulate; the ε < 0 primary branch; ε = 0 as the ΛCDM null.
- **Profile:** the ε grid (20 points, 0 to −1).
- **Data and physics:** the PRIMARY dataset (DESI DR2 BAO + DESI DR2 baseline CMB); cosmological parameter ranges and priors; CAMB physics and accuracy settings; likelihood versions and data, hashes unchanged.
- **Comparison and decision:** the CPL comparator definition (including the frozen `w0 + wa < 0` constraint); thresholds S1–S8; the S5 `decide()` function; the χ²_eff definition (CR-3).

None of these is altered in response to the observed v1 values.

## 2. New: the numerical procedure (CR-5, owner-authorized for v1R only)

**One uniform procedure** applies to every target:
- all 20 fixed-ε Card profile points (ε = 0, the ΛCDM null, included);
- the free-ε Card minimization;
- the CPL comparator.

Each target is refined in a series of **cycles**:

| Element | Frozen value |
|---|---|
| Optimizer | Cobaya BOBYQA |
| `rhoend` | **0.005** (v1: 0.05) |
| Starting point of each cycle | exactly the current best point |
| Proposal covariance | from the official, checksum-verified DESI DR2 chains for the primary combination. `v1r_covmat_lcdm.txt` comes from `base/` and is used for Card / ΛCDM / free-ε. `v1r_covmat_w0wa.txt` comes from `base_w_wa/` and is used for CPL. Both use the post-burn-in (0.3) weighted samples |
| Convergence | one complete cycle improves χ² by **< 0.1** (so at least 2 cycles) |
| Cap | **6 cycles**. A target that has not converged is **NOT CONVERGED**, and no state is forced |
| Early stop | **none** based on which side of S3/S5 a value lies |

**Arm A (all 22 targets)** starts from the **best valid v1 solution** for that target, taken from the v1 raw outputs. The archived mis-seeded free-ε fits are excluded.

**Neighbor warm starts** (permitted as optional by the ruling) are **not used**. Every fixed-ε point starts from its own v1 best.

**Smoothing:** no smoothing, interpolating away or replacing of computed points. Interpolation is used **only** for the registered descriptive 95 % crossing, after the profile values exist.

## 3. Independent convergence check (arm B)

Every **load-bearing** target is minimized again, with the same CR-5 cycles, from a **distinct admissible start**:

| Target type | Arm-B start |
|---|---|
| Card fixed-ε / ΛCDM | a fresh seeded draw from the v1 reference distributions (seed base 90017, distinct from the v1 seeds 17 / 1017 / 2017) |
| free-ε | the same kind of draw, with ε starting at **−0.20** |
| CPL | the **highest-likelihood sample of the official DESI DR2 base_w_wa primary chains** (`v1r_cpl_startB.json`: w0 = −0.4248, wa = −1.7630). Ruling §5: a numerical starting location only; it is never used as a fit result, and no posterior mean is imported |

**Acceptance.** A target passes only if all of the following hold:
- both arms reached cycle convergence;
- `|χ²_best,A − χ²_best,B| ≤ 0.2`, retaining the frozen v1 0.2 scale.

A failure marks that target numerically unresolved.

**Load-bearing set** (identified mechanically from the best-of-arms profile):
- the ΛCDM null;
- the CPL comparator;
- the free-ε fit;
- the best grid point;
- the grid points bracketing each descriptive 95 % crossing (Δχ² = 3.84).

The set is re-identified after the arm-B results, up to 3 times, so that any newly load-bearing target also receives its check.

## 4. Result reporting

- **A. Numerical status:**
  - `CARD-01-v1R — NUMERICALLY RESOLVED` if every load-bearing target passes;
  - otherwise `CARD-01-v1R — STILL NUMERICALLY UNRESOLVED`.
- **B. Mechanical S5 state:** reported **if and only if** the status is RESOLVED, using the unchanged `decide()`. It is labeled **POST-DATA VERIFICATION-REPAIR RESULT — NOT THE PREREGISTERED v1 RESULT**. It is never upgraded to PREDICTED or VALIDATED. Card #1 remains POSTULATED.

## 5. Scope

**PRIMARY only.** The following are **not run**:
- the Pantheon+, Union3 and DESY5 extensions;
- the ε > 0 control;
- the posterior MCMC.

The driver stops for owner review.

**Not touched:**
- any other repository;
- `grut-backreaction-identifiability-0`;
- every frozen scientific branch and every governance branch.

## 6. Files (committed with this charter, before execution)

| File | Content |
|---|---|
| `code/card01_v1r_config.py` | frozen numerical constants |
| `code/card01_v1r_refine.py` | one target × one arm, CR-5 cycles |
| `code/card01_v1r_driver.py` | arm A for all targets, then arm B for the load-bearing set, then STOP |
| `code/card01_v1r_analyze.py` | convergence table, numerical status, gated `decide()` |
| `code/card01_v1r_inputs.py` and `code/v1r_inputs/` | covariances, CPL arm-B start and nuisance-parameter definitions, all derived from the official chains and the Cobaya likelihood defaults |
| `CARD_01_V1R_OFFICIAL_CHAINS.sha256sum` | DESI-published checksums of the official ΛCDM chains used for the covariance (6/6 OK). The w0wa chain checksums are already in `CARD_01_D1_OFFICIAL_CHAINS.sha256sum` |

**Provenance preserved:**
- all v1 raw outputs;
- the archived mis-seeded runs;
- `826f0d1`;
- the original run configuration;
- the owner threshold ruling.
