# Recorded primary outcomes

Times are completed whole physical steps. A mid-sweep stop does not
provide a stored state at that time or a target-time probability.
Retained tensors/coefficients are not total resident memory.

| Protocol | Remaining counts | Status | Completed / required | Wall s | CPU s | Peak MiB | Accumulated ledger |
|---|---|---|---:|---:|---:|---:|---:|
| Original S8 | [2, 2] | UNRESOLVED (mid-sweep) | 30.065359 / 100 | 53.022 | 53.020 | 56.305 | 0.388869 |
| Original S8 | [3, 1] | UNRESOLVED (mid-sweep) | 57.314410 / 100 | 185.399 | 185.371 | 58.746 | 0.424129 |
| S8-G1 | [2, 2] | UNRESOLVED (mid-sweep) | 88.344227 / 100 | 839.938 | 839.900 | 90.102 | 0.464346 |
| S8-G1 | [3, 1] | COMPLETE | 100.000000 / 100 | 778.287 | 778.238 | 72.445 | 0.335272 |
| Fresh S7 | [2, 2] | COMPLETE | 100.000000 / 100 | 118.904 | 118.890 | 86.496 | 0.549138 |
| Fresh S7 | [3, 1] | COMPLETE | 100.000000 / 100 | 101.578 | 101.569 | 66.062 | 0.547993 |

All limits, allocation choices and error components remain in the raw JSON.
Original S8 has a 16-per-charge cap and effective total <=80 at N=5;
G1 has a global cap 96. S8/G1 have 900 seconds; S7 has 600 seconds.
These are different numerical protocols, not an equal-rank/equal-error
speed comparison. A stopped run does not demonstrate a cheaper target.
A mid-sweep ledger does not certify its partially updated state against
physical H evolution at the last completed-step time.

| Complete N=5 family | Macro lower | Full-assembly lower | Recovery gate |
|---|---:|---:|---|
| Original S8 | No target-time bound | No target-time bound | Not passed |
| S8-G1 | No target-time bound | No target-time bound | Not passed |
| Fresh S7 | 0.698713 | 0.579462 | reference |

Displayed probability lower endpoints are rounded downward; exact
endpoints and any upper bounds are in summary.json.
