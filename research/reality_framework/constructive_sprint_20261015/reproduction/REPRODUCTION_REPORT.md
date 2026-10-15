# Fresh reproduction of the charged-accessibility packet

The previously pending packaged reproduction at commit
`1801afee27096c03c2a3b64cc811ac0be0d2a4a7` has now been executed successfully.
Both the exact verifier and the complete numerical audit passed without changing
their published source. This resolves an execution gap; it does not establish
macroscopic geometry formation.

All 15 immutable packet files were restored through GitHub blob reads. Their
Git blob hashes and byte lengths matched before execution. The exact verifier
also checked its packet manifest and six frozen parent dependency hashes.
The earlier unpublished local one-target files were preserved separately under
`unpublished_local_backup/`; they are not substituted for the published packet.

## Execution

| Calculation | Result | Wall time |
|---|---|---:|
| Published `node verify.mjs` | Pass | 190.54 s |
| Published `audit.py`, complete and without cache | Pass | 531.86 s |
| Fresh Python one-target integer recurrence | Pass | 0.124 s |

The numerical audit built both models afresh with `CP12_FORMATION_CACHE` unset.
OpenBLAS, OMP and MKL thread limits were one. Its internal calculation timer
reported 527.37 s; the wrapper wall time includes process startup and imports.
The two primary runs were launched concurrently. Exact shell commands,
environment, exit status and timestamps are retained in the execution JSONs.

The exact verifier reproduced all 960 minimum-order forward amplitudes, the
eight independently reversed representative amplitudes, the total integer
norm, and the rational order-44 probability coefficient. It also executed
20 collective-bound evaluations and 680 independent-spin remainder checks.
The additional Python control returned the published signed amplitude
`-49583226880` for the same labelled order-22 target. This is a fresh
cross-language implementation control, not an additional physical assumption.

## Numerical agreement

The complete N5 results matched the saved packet values exactly at the printed
precision. N4 probability differences were at most approximately
`7e-18`; the projected time-average difference was approximately `8.4e-17`.
The comparison uses explicit mixed tolerances and records small residuals
without requiring bitwise agreement at the eigensolver rounding scale.
The analytic table differs at one last bit: the N=1024, T=0.75 probability
upper bound is `7.129462761222235e-165` in the fresh Node execution versus
`7.129462761222236e-165` in the saved table. The new comparator initially
required bitwise table equality; it was corrected to use a `2e-13` relative
floating tolerance and preserves this discrepancy in `COMPARISON.json`.
The approximately `1.4e-16` relative change has no physical significance.

| N5 charged full-circulation target, Q=10 | Retained probability | Full-law probability enclosure |
|---|---:|---:|
| T=1 | 2.1256306851140608e-19 | [0, 7.299721998951689e-13] |
| T=2 | 2.2586580273902814e-9 | [0, 1.0313583399591289e-5] |
| T=3 | 3.744554469770909e-7 | [0, 0.0025380112417202954] |

The N5 closure has 39,969 retained signed permutation orbits, 17,418 outside
orbits, and 66 target orbits. The outgoing residual includes omitted electric
configurations. All three lower bounds remain zero.

The N4 spectral benchmark reproduced 577 retained modes, 86 target orbits,
initial-weighted target Gram sum `0.14752894104361644`, and weighted full-H
trial residual `1.4352906680026263e-13`. For example, its T=5 target probability
is `0.05502210766596675`. This Q=4 benchmark is not the Q=2N extensive family.
The full-H residual bound retains its time-dependent defect; it is not an
all-time spectral bound after discarding that defect.

## Interpretation and limitations

No scientific discrepancy or required packaging fix was found. The original
packet's static verification output still says that the packaged Python rerun
had not happened at its publication time. This new report and execution data
record that it has now happened; the immutable historical record is preserved.

This reproduction confirms the reported finite calculations. Exact charged
accessibility and early-time collective suppression remain compatible: a
nonzero coefficient does not imply appreciable formation probability. N5 at
Q=10 has the forced graph K5, so it cannot discriminate dynamically selected
large-scale host geometries. The numerical enclosures evaluate analytical
electric-truncation inequalities in floating point; they do not provide formal
outward-rounded machine arithmetic or external mathematical certification.
Persistence, dimension selection and spatially constrained propagation are
not established by rerunning this packet.

## Reproduce

From the repository root:

```sh
python research/reality_framework/constructive_sprint_20261015/reproduction/run_reproduction.py node
python research/reality_framework/constructive_sprint_20261015/reproduction/run_reproduction.py python
python research/reality_framework/constructive_sprint_20261015/reproduction/python_onset_check.py
node research/reality_framework/cp12_charged_accessibility_gate/analytic_bounds.mjs > research/reality_framework/constructive_sprint_20261015/reproduction/ANALYTIC_BOUNDS_RERUN.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/constructive_sprint_20261015/reproduction/compare_results.py
```

The immutable source report is
[`../../cp12_charged_accessibility_gate/REPORT.md`](../../cp12_charged_accessibility_gate/REPORT.md).
Fresh results are in `python.stdout.log`, `node.stdout.log`,
`PYTHON_ONSET_CHECK.json`, `ANALYTIC_BOUNDS_RERUN.json` and `COMPARISON.json`.
