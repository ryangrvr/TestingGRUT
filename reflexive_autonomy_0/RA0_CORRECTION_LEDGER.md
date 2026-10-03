# RA0 CORRECTION LEDGER (additive)

**Rule:** logs are kept as emitted. Where document wording differs, this ledger takes precedence.

## RA0 REPAIR 01 (owner review of `0870843`)

G1 is accepted after this repair. G2 is accepted provisionally after it.

| ID | item | correction |
|---|---|---|
| **RA1-01** | G1 PROP 6 codimension (k − 1)(n − k) | Regraded to **KNOWN / REDERIVED IN RA0**. It matches the asymptotic chi-square degrees of freedom of classical lumpability tests for a fixed k-block lumping (owner-cited; Statistics & Probability Letters, ScienceDirect S0167715203001263; not re-read here). The internal proof is retained; the measure-zero corollary stands. No novelty claimed |
| **RA1-02** | G1.6 autonomy wording | A factor partition is lumpable iff its own transition law is free of hidden coordinates: P(A_{t+1} \| A_t, B_t) = P(A_{t+1} \| A_t), i.e. **no back-action onto the coarse variable**. Full factorization P_A ⊗ P_B is **not** required; one-way / triangular interacting dynamics may qualify. "Exactly non-interacting factor" → **"exact autonomous factor / no back-action onto the coarse variable"**. Finite reflexive autonomy detects autonomous factors, not necessarily independent ones. The frozen-witness conclusion is unchanged |
| **RA1-03** | G2 B2 mean-density defect | The original code maximized each source class independently, ignoring Σ_i n_i = N. Its output is regraded **MEASURE-FREE UPPER BOUND ON THE CONSTRAINED MEAN-DENSITY DEFECT**. **Repaired:** `g2/g2_b2_constrained.py` computes the exact constrained supremum by dynamic programming. It **agrees with the upper bound to 4 decimals** at L = 64 … 4096 for all three decompositions. The L^(−1/2) trend for √L cells is unchanged. The old log is preserved. B2 is not load-bearing |
| **RA1-04** | G2 B3 slow eigenspace | Preserved: at L = 8, 10, 12 the lowest nonzero eigenspace lies in the linear-occupation span, and the lowest gap equals the one-particle gap (consistent with Caputo–Liggett–Richthammer). Grade: **CANONICAL FINITE-SIZE LOWEST EIGENSPACE**. Density-mode identification: **NUMERICALLY EXACT AT TESTED SIZES / CONSISTENT WITH KNOWN STRUCTURE**. The gap theorem is **not** promoted to "the whole asymptotic slow sector is density modes" |
| **RA1-05** | G2 "canonical bounded-ratio slow subspace" | **Withdrawn.** "All modes with λ/λ₁ bounded as N → ∞" needs a cross-N mode correspondence; "λ ≤ Cλ₁" needs a supplied C. Replaced by **INTRINSIC SPECTRAL SLOWING / SCALE HIERARCHY PRESENT; A UNIQUE ASYMPTOTIC SLOW SUBSPACE IS NOT YET DEFINED**. The finite-N lowest eigenspace remains canonical |
| **RA1-06** | G2 twin diagnostic | Contiguous vs interleaved: **DYNAMICAL DIAGNOSTIC DISTINGUISHES THE TWO SUPPLIED CANDIDATES**. It does not define or select a partition. No Σ selector |

**Applied to:** `G1_FINITE_REFLEXIVE_AUTONOMY.md`, `G2_ASYMPTOTIC_AUTONOMY.md`, `RA0_ZOOM_OUT_01.md`,
`RA0_ZOOM_OUT_02.md`.

## G3 note (no repair)

- **Numerical defect fixed before the logged G3 run.** At L = 10⁴ the one-particle levels (~10⁻¹¹) fell below an absolute
  zero / merge tolerance, producing a spurious "ratio 1.78 at rank 8". Fixed with relative tolerances and 1 − cos = 2 sin²
  (no cancellation). The corrected value is ratio 4.000 at rank 2.
- **Metastable printout:** a crash (only one cut in a degenerate spectrum) was fixed.
- **D3:** a non-symmetric metastable control was added after the first run, to test whether algebra closure in D1 is an
  artefact of basin symmetry. It is: closure is exact only for symmetric basins, and asymptotic otherwise.
