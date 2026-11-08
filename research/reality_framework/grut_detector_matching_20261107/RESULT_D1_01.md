# D1-01: energy populations miss future histories; 32 retained coherent levels also miss

Frozen before the new timings were evaluated: `90745ef0d7aaf4c6d8fe3a11edc7ee0d12f7d339`. [Protocol and bound](PROTOCOL_D1_01.md), [full raw result](raw_d1_01/result.json). Execution used a clean checkout and took 36.1 seconds. This is a test of the existing 1280-state interacting parent, not a reopening of its failed continuum gate.

The dephased comparator matches the single-time calibration, retains all 1280 energy populations, and preserves mean field energy at each dephasing operation. It nevertheless gives different later detector histories. Coherences discarded after an instrument are operationally consequential in this finite interacting model.

| New time protocol | Energy-population comparator TV | 32 coherent levels TV | Bound for 32 levels |
|---|---:|---:|---:|
| (0,0.7,1.9) | 0.02214397 | 0.00240298 | 0.17763962 |
| (0,0.6,1.5,2.4) | 0.03745807 | 0.00427073 | 0.32145093 |

The registered contrast threshold, 0.001, is exceeded. The registered 32-level error window, 0.00001, and bound window, 0.001, both miss. Retaining coherence helps, but this chosen energy reduction is not adequate at its requested accuracy. No larger K is substituted for the frozen primary.

| Retained levels | Three-time TV | Four-time TV |
|---:|---:|---:|
| 4 | 0.04031218 | 0.06310531 |
| 8 | 0.01410558 | 0.02407306 |
| 16 | 0.00706668 | 0.01190588 |
| 32 | 0.00240298 | 0.00427073 |

The constructive bound dominates every measured discrepancy, and overflow accounting closes. Its looseness is visible in the table; it is not a 0.1% certificate. The parent itself already has overflow 0.00134507 and 0.00183825 on these new protocols. The reduced calculation records additional projection losses separately, approximately 0.00249920 and 0.00442731 at K=32. Those losses occur on the reduced paths and must not be interpreted as an independent additive estimate of the parent's continuum error.

The dephasing comparison is a declared change of instrument dynamics, implementable by coupling to an additional record/environment. It is not a claim that field coherence vanishes without a physical interaction. Its first-outcome calibration and immediate energy populations match; later energy exchanges need not match after subsequent measurements. It rejects the particular energy-population-only description for these histories, not all possible reduced descriptions, classical effective regimes, or decoherence models.

The result connects the integrated framework's sufficiency criterion to an interacting field realization. The claimed mechanism is standard quantum interference and measurement backaction. No restriction beyond the adopted quantum and operational data has been derived. M1-01 remains unresolved; Zenodo remains paused.
