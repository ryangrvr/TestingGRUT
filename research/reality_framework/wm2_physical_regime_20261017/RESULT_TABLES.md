# Frozen finite recovery outcomes

Model target time is 100 in receiver-hopping units. Worker wall/CPU
seconds are computational resources. A stopped run earns no target-time
probability. The ledger column describes accumulated costs only.

| Remaining counts / original weight | Status | Completed model time | Wall s | CPU s | Peak MiB | Ledger |
|---|---|---:|---:|---:|---:|---:|
| [2, 2] / 3/5 | UNRESOLVED | 65.032680 | 1202.430 | 1200.406 | 68.305 | 0.143917 |
| [3, 1] / 2/5 | UNRESOLVED | 76.746725 | 1208.667 | 1206.954 | 59.840 | 0.171175 |

| Counts | Orthogonal loss | Linear in-range defect | Strang | Field product | Stop |
|---|---:|---:|---:|---:|---|
| [2, 2] | 0.05979118 | 4.19e-08 | 0.07306207 | 0.01096802 | whole-worker wall-time ceiling |
| [3, 1] | 0.07188512 | 3.18e-08 | 0.08610870 | 0.01306111 | whole-worker wall-time ceiling |

The worker limit is 1200 seconds, checked at complete physical steps;
overshoot is retained in each raw JSON. Rank ceiling 96, RSS ceiling
2048 MiB, and all error gates remain frozen. Retained tensor entries
are not total memory. No equal-error/budget speed claim is made.


Full original-weight recovery: **FAIL**.
No combined T=100 family interval is claimed.
No fallback or six-event target belongs to S9.
