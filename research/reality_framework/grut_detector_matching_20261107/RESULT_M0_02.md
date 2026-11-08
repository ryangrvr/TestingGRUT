# M0-02: endpoint improvement is substantial, but the frozen gate still misses

Source and protocol were frozen in commit `c5e107a9920c57b13482d63f290fa9f7bc1fc7b5` before execution. This was a separately declared hypothesis following the known M0-01 miss. [Protocol](PROTOCOL_M0_02.md), [raw analysis](raw_m0_02/analysis.json), [source provenance](raw_m0_02/provenance.json).

The analytically motivated ghost-endpoint convention reduces the maximum N=20 history discrepancy from 0.01001135 to **0.00167772**. The fixed acceptance window is 0.001; the gate therefore remains a **miss**. Two of four finest-point protocols pass individually, but the complete comparison does not.

| Sites | TV at (0) | TV at (0,1) | TV at (0,2) | TV at (0,1,2) |
|---:|---:|---:|---:|---:|
| 10 | 0.00141892 | 0.00208079 | 0.00464575 | 0.00464385 |
| 14 | 0.00081338 | 0.00127063 | 0.00262637 | 0.00284452 |
| 18 | 0.00054522 | 0.00089966 | 0.00172750 | 0.00195740 |
| 20 | 0.00046283 | 0.00077995 | 0.00145567 | 0.00167772 |

The last three points decrease monotonically for every protocol. The largest 18-to-20 change, 0.00027968, remains above the frozen 0.00025 refinement threshold. No continuum extrapolation replaces these failed finite-resolution gates.

The N=20 mean differs from the continuum mean by 0.00015344, compared with -0.01581027 in M0-01. Its gap error decreases from -0.04670902 to -0.00295126. The independent free-chain identity underlying the change was verified before production; it does not prove an interacting boundary theorem. Bulk mass, detector quadrature and endpoint errors have not been independently separated.

The declared alternating-constant subtraction control reduces the N=20 maximum TV to 0.00125273; it also misses 0.001 and does not replace the primary verdict. The frozen generic analyzer reports the raw mean for both controls; their physical subtracted means are the `continuum_observable_mean` fields in each raw point. A descriptive string in the M0-02 JSON still mentions the older bare-zero control, but the registered tasks and executed code use the specified smoothed controls. That stale string is preserved and disclosed here; it did not change the primary mass or calculation.

Reference and numerical diagnostics pass. Tighter eigensolving and split propagation change N=14 histories by at most 1.94e-15. The analytic continuum tail bound remains 5.50e-7. These checks bound different contributions and do not turn a finite-step trend into a continuum error theorem.

**Conclusion:** the boundary dictionary matters materially, and the selected massless histories approach the independent continuum target. Matching is still unresolved at the registered accuracy and resource boundary. The improved dictionary can be investigated in a separately specified massive reference, but this result supplies no massive error bound. Causal implementation of the ideal extended detector, uniform continuum matching and distinctive GRUT physics remain open. Zenodo remains paused.
