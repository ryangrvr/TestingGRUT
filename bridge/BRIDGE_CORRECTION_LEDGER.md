# BRIDGE CORRECTION LEDGER (Bridge-1)

**Ordering.** Corrections are additive and owner-directed. Historical script logs (`bridge/b1/*.log`) are kept exactly as
emitted. Where a log line's wording is broader than the corrected scope, the correction below takes precedence.

## BRIDGE REPAIR 01 — B1 scope wording (owner, after `add9f79`; none of these changes the verdict)

| ID | location | original wording | corrected reading | effect on the verdict |
|---|---|---|---|---|
| **BR1-01** | `B1_SIGMA_LOCAL_NET_RESULT.md` (chain row, memory-kernel row); `BRIDGE_ZOOM_OUT_01.md` §1; ledger BI-05; log line "Lanczos from ANY unit start vector" | "every q gives a connected tridiagonal chain frame" | **Every cyclic start vector q gives a full connected Jacobi / tridiagonal representation.** For a simple spectrum the cyclic vectors form an open, dense, full-measure set. Exceptional non-cyclic q (orthogonal to some eigenvector) terminate Lanczos early | none: the continuum / non-uniqueness result is unchanged |
| **BR1-02** | `B1_SIGMA_LOCAL_NET_RESULT.md` (sparsity row); `BRIDGE_ZOOM_OUT_01.md` §1 | "optimum of every factorization criterion" | **"optimum of the tested entry-sparsity / autonomy / response-factorization criteria."** No claim is made about arbitrary possible criteria | none |
| **BR1-03** | `B1_SIGMA_LOCAL_NET_RESULT.md` (C1-a row); `BRIDGE_ZOOM_OUT_01.md` §1; ledger BI-02 | "known theorem: path DLS" | the fact that the unit-weight path is determined by its Laplacian spectrum is a **KNOWN-RESULT IMPORT**. Bridge-1 gives no independent proof | none: still **CONDITIONAL COMPRESSION / CRITERION-PRICED** |

**Status:** B0 and B1 are **accepted provisionally** by the owner, with BR1-01 … BR1-03 applied.
