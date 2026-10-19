# Execution and provenance record

This additive packet continues the work of 18 October 2026 (America/Chicago); final integration and fresh replay were executed on 19 October UTC. The folder's 20261018 identifier is retained. The parent repository state is `c4ea0604e183f3dd62d11516c2536359774489e7`, the already posted classical-stability follow-up. Candidate 1's construction is at `c2e2cb6d7be1b68ce2b4370d5fb87b08412d5e8e`.

## Executed calculations

| Calculation | Implementation and frozen protocol | Execution evidence |
|---|---|---|
| Exact scalar beta functions and RG identity | `rg_symbolic.py`, `rg_protocol.json` | `rg_results.json`, `rg_polynomials.json`, `rg_execution.log`; 290 exact checks, zero polynomial residual |
| Constant-background potential, analytic derivatives, vacuum and bulk minima, negative-mode maps | `effective_potential.py`, `potential_protocol.json` | `potential_results.json`, `potential_execution.jsonl`, two raw CSV scans; final recorded run about 1.58 s |
| Numerical RG flow and scale/order tests | `flow_check.py`, `flow_protocol.json` | `flow_results.json`, `flow_scale_values.csv`, `flow_execution.jsonl`; six declared acceptance conditions, about 0.024 s |
| Cross-implementation and boundary checks | `verify.py` | `VERIFICATION.json`; 39/39 checks |
| Scientific figure | `plot_results.py` | `quantum_diagnostics.png`; generated from the recorded scans and visually inspected |
| End-to-end fresh-directory replay | `reproduce.py` | `REPRODUCTION.json`; all five scripts return zero and all 16 substantive comparisons pass |

All calculations used the existing local Python runtime. Version strings are recorded in the outputs. The symbolic calculation needs only the standard library; numerical calculations use NumPy/SciPy, and the optional figure uses Matplotlib. No installation, account access, paid computation, new dataset, external experiment or indefinitely running job was used. Repository posting uses the owner's existing authorization.

## Errors and retained corrections

- The first potential result reused the optimized bulk mass array in a field labeled as the tree-background masses. The corrected source has distinct variables and output labels. The original JSON/log and `potential_correction_record.json` preserve this event. Tadpole, CW value, bulk correction and vacuum-energy results were unchanged.
- The integrating verification script initially applied Python's scalar `abs` to a JSON list. It raised `TypeError` before completing verification. Replacing that operation with NumPy's array absolute value corrected the check; the physical calculations and thresholds were unchanged. Final verification passes 39/39.
- The first fresh-replay wrapper copied `rg_results.json` into its new destination before running a symbolic script that intentionally refuses to overwrite that file. The replay failed at that guard, with no change to the posted calculation. The wrapper now omits that generated file when copying to a new destination. A second fresh directory completed all five scripts and all 16 comparisons. The failed destination was not reused or treated as a successful replay.

The three frozen physical protocols were not altered in response to an unfavorable physical output. No physical parameter or matching constant was tuned.

## Division of work and review limit

The constructor derived/implemented the exact RG sector; the calculator implemented the scalar potential and derivatives; the challenger inspected their actual source/results and normalization; the integrating worker implemented the independent numerical Hessian/flow, cross-checks, replay and report. These workers share the model/framework and project context. This is internal adversarial review and numerical cross-implementation, not independent peer review or experimental validation.

`NEXT_CALCULATION.md` derives a starting fluctuation operator and specifies a subsequent test. Its spectrum, determinant, quantum energy, pole mass and lifetime were **not executed here**. Nothing in the execution count refers to a universe simulation or empirical test.

## Preserved inputs and scope

The only numerical profile imported is the original Candidate 1 `code/results_profile.npz`, SHA-256 `0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3`. The replay copies this one input into the same relative layout. Candidate 1 and classical-stability files are unchanged. The packet manifest records its own delivered source, protocols, outputs and reports; generated caches and the manifest itself are excluded from that list.
