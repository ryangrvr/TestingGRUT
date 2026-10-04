# CARD #1 v1 — OWNER CONVERGENCE RULING

**Issued by** the scientific owner, after the executor reported a minimizer start-to-start scatter of 2.5–5 in χ², which is comparable to the S3 cut `I_card ≥ 4`.

**Transcribed by** the executor. The owner threshold ruling (`86ccbb2e`) and the run configuration (`c9899abc`, plus the pre-evaluation corrections CR-1 to CR-4 at `b50f5224`) are **not amended**.

## Decision: KEEP THE FROZEN PROTOCOL

- **CR-5 is not added** to Card #1 v1.
- The number of starts is **not increased** from 3 to 5.

**Reason.** The frozen run configuration states that nothing in it may be changed after data access. It freezes:
- BOBYQA;
- 3 starts;
- the minimizer protocol;
- `CHI2_START_DISAGREEMENT_FLAG = 0.2`.

Likelihood values have now been inspected, including values near the owner-locked `I_card = 4` threshold. Changing the convergence tolerance, proposal covariance or refinement logic now would be a post-data modification of the frozen scientific run configuration.

## 1. Allowed execution repairs

These do not alter the model, likelihood, minimizer definition or decision thresholds:
- reduce the worker / thread count to avoid OOM;
- restart crashed jobs;
- preserve completed valid jobs;
- add a driver guard that prevents the free-ε minimization until the full frozen grid is complete;
- archive and exclude the two prematurely started, wrongly seeded free-ε fits;
- rerun those free-ε fits later under the frozen protocol, with the correct best-grid seed.

All of these are documented in `CARD_01_EXECUTION_REPAIRS.md`. The archived erroneous runs are **not deleted**.

## 2. Finish the primary frozen three-start profile first

The archived wrongly seeded free-ε fits are used for **no statistic**.

## 3. Preserve the frozen convergence diagnostic

For every load-bearing minimization, report:
- all three start χ² values;
- the best χ²;
- the start-to-start spread;
- whether `spread > 0.2`.

Starts are **not averaged** and are **not replaced** by a tighter optimization.

## 4. Mechanical result vs accepted scientific verdict

Apply the frozen S5 `decide()` exactly as committed, and report its raw output.

If unresolved minimizer disagreement greater than 0.2 affects any quantity needed to distinguish among the S5 result states, the run is marked:

> **CARD-01-v1 — NUMERICALLY UNRESOLVED / NO ACCEPTED SCIENTIFIC VERDICT**

The quantities that count are:
- the ΛCDM null;
- the Card branch best fit / profile;
- the CPL comparator;
- the profile-set crossings.

In that case the mechanical S5 output is still shown, but labeled:

> **MECHANICAL OUTPUT — NOT SCIENTIFICALLY ACCEPTED DUE TO FROZEN CONVERGENCE FAILURE.**

This label is not EXPLANATORY-SURVIVES, EXPLANATORY-KILLED or MODEL-KILLED. Its basis is the already-frozen `CHI2_START_DISAGREEMENT_FLAG = 0.2`.

If the completed frozen run resolves the relevant start disagreements below 0.2, S1–S8 are applied normally.

**Implementation** (`code/card01_d3_analyze.py`). The **load-bearing minimizations** are:
- the ΛCDM null (ε = 0);
- the CPL comparator;
- the branch best grid point;
- the free-ε fit;
- the grid points bracketing each 95 % profile-set crossing.

If any of them has a spread above 0.2, the run is NUMERICALLY UNRESOLVED.

## 5. CR-5 — recorded verbatim, not executed

> **POST-DATA REPAIR PROPOSAL — NOT PART OF CARD #1 v1 — NOT YET AUTHORIZED.**
>
> CR-5 (executor proposal): after the three-start stage, refine every point the same way — every ε, the ΛCDM null, the CPL comparator and the free-ε fit. Each refinement restarts from that point's best start with a 10× tighter convergence tolerance (BOBYQA `rhoend` 0.05 → 0.005) and a proposal covariance taken from DESI's official chains, repeating until χ² improves by less than 0.1. Because it's the same rule for every model, it can't favour one over another.

If v1 is numerically unresolved, the owner will separately decide whether to open a repaired computational analysis. Its convergence protocol would be frozen before that repair runs. It would be labeled post-data / verification-repair and would **not** be described as the preregistered v1 test.

## 6. Extensions

The three SN extensions and the positive-ε control are **not run** until the primary frozen run is complete and its numerical convergence status is known. Stopping after the primary for owner review changes nothing in the primary.

## 7. Record integrity

The execution-repair log and the driver guards are committed separately.

**Owner note:** the preliminary numbers **must not be interpreted scientifically**. The tempting figure 10977.47 − 10974.07 ≈ 3.4 comes from wrongly seeded free-ε fits, and the minimizer scatter is several χ². "The useful result right now is methodological: the preregistered optimizer may be too loose to resolve the preregistered scientific threshold."
