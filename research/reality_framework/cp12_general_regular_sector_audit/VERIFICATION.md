# Executed verification snapshot (2026-10-13)

## Parent `0463ec1` verbatim independent run

The parent repository's `cp12_regularity_fragmentation_review/check.py` Git blob SHA was verified as `798e922ffb591a3ecfcb375a102b461f1da4005d` in the independent runtime. The unchanged script printed `CP12_REGULARITY_FRAGMENTATION_INDEPENDENT_PASS`.

Regular graph and one-defect counts for N=4..7: (3,12), (12,120), (70,1080), (465,10500).

At N=1000, the exact recurrence gives single-cycle probability for component penalties betaDelta = 0,1,2,4 respectively: 0.05936593571247547; 0.3413342945351178; 0.6695270403884722; 0.946752189196841.

## New independent check.py

**Git blob SHA `abc2744b6d51d222d3f71a9fd3e9f8faadefc5f4`** exactly matches the Python file run independently before publication. Output marker: `CP12_GENERAL_REGULAR_SECTOR_CHECK_PASS`.

| (N,k) | k-regular graphs | distinct one-pair defect graphs | maximum predecessor multiplicity |
|---|---:|---:|---:|
| (6,2) | 70 | 1080 | 3 |
| (6,3) | 70 | 1080 | 3 |
| (6,4) | 15 | 180 | 2 |
| (7,2) | 465 | 10500 | 3 |
| (7,4) | 465 | 10500 | 3 |

The script verifies counting and switching, **not** the infinite-dimensional physical quantum operator norm.

## Optional geometry control

**Git blob SHA `690b69cceaa7fe7ce00d31f7d507046a9058d5bc`** exactly matches locally executed `geometry_control.py`; requires networkx (tested 3.6.1). Output marker: `CP12_SPATIAL_METRIC_ADVERSARY_PASS`.

At N=343, Q=1029 and degree=6, the periodic 3D torus ball sizes for radius 0..6 are [1,7,25,63,123,195,263], diameter 9. Seeded random regular graphs 1..5 are connected, have radius-three balls [140,148,153,152,146], diameter 5. These are graph-metric, not operational quantum detector, results.

## Review boundary

Independent execution and Git blob identity are established for the listed scripts. The generalized full CP-12 quantum cut proof and the limiting regular-graph enumeration **still require** Claude Code's independent analytical and regulator audit. No primary GRUT model or branch was changed.
