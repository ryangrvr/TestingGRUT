# Execution disclosure

The prospective physical scaling protocol was not executed. The only new
executions are `controls.py` and its identical-source replay. They perform
exact Fraction calculations and L=4,6,8 matrix checks, without any trajectory.
All historical code and records remain unchanged.

The source and DERIVATION.md hashes were frozen before both executions.
ROUTE_AND_DECISION.md contains an explicitly later analytic velocity upper
bound, not a new control or revised protocol. It has no new numeric run.

| Metric | Primary | Replay |
|---|---:|---:|
| Program seconds | 0.496398895 | 0.469613960 |
| Whole worker seconds | 1.127843921 | 0.831278327 |
| Peak RSS MiB | 65.1328125 | 64.0468750 |

Combined whole-worker time: 1.959122248 seconds, below the 60-second ceiling.
Single-thread BLAS variables were set. Workers used a 20-second timeout;
RSS/storage ceilings were monitored, not enforced by OS resource limits.
No paid allocation, installation, persistent service or background work.
Python 3.12.14, NumPy 2.3.5; runtime details are archived in EXECUTION.json.
CONTROLS.json is SHA256-identical on replay. Metrics naturally differ.

The controls reuse B1's `model.py` and B2's `Colored` builder. They do not run
B2's main trajectory suite or its Propagator. Their matrix checks use
ordinary binary64 eigh and sparse arithmetic with tolerances. These are
diagnostics, not conditional propagation enclosures or interval certificates.
Exact arithmetic applies only to explicitly recorded rational expressions.

One assistant performed all research roles. No subagent, independent
implementation, external reviewer, empirical data, or mathematical
validation by another party is claimed. The bounded cycle stops with its
mathematical status unresolved for the operating task.
