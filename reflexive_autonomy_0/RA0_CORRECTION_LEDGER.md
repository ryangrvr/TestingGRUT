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

## RA0 REPAIR 02 (owner review of `e263656`)

G3 is accepted provisionally after this repair. All numerical values are preserved; logs are kept as emitted.

| ID | item | correction |
|---|---|---|
| **RA2-01** | D3 algebra closure | "The D3 slow subspace becomes an algebra asymptotically" → **"In the tested non-symmetric metastable sequence, the product-closure residual decreases strongly with M. This is numerical evidence for asymptotic algebraization, not a proof of convergence to zero."** Reason: each M used a newly generated random matrix, so the tested objects are not a nested realization of one family with a proved limit. Additionally (self-reported): the G3 residual was measured on eigenbasis pairs in the Euclidean norm, which is basis-dependent; G4 replaces it with the basis-free L²(π) defect. This ledger row supersedes the G3-note sentence "closure is exact only for symmetric basins, and asymptotic otherwise" |
| **RA2-02** | finite-N canonicality vs cross-N hierarchy | **Earned, every finite N:** ordered levels intrinsic; spectral projectors at a specified gap basis-free; degeneracies need no eigenbasis choice. **Conditional:** a cross-N filtration is canonical only if the sequence of gap ranks is itself specified by spectral data without matching modes. Established in the constructed fixed-rank controls D1, D2 and the ladder β > 2. **No general theorem for arbitrary families** |
| **RA2-03** | SSEP counting profile | "2√x" → **for fixed x and large L, λ_q/λ₁ → q² implies the ±q counting profile approaches 2⌊√x⌋, away from finite-size / Nyquist qualifications.** At the tested squares x = 1, 4, 9, 16 this is exactly 2, 4, 6, 8. No G3 terminal changes |
| **RA2-04** | positive-control interpretation | D1 – D3 were **built** with metastable scale separation. They establish that the no-ε spectral criterion detects true asymptotic scale separation and returns the correct basis-free projector / filtration in positive controls. They do **not** show that generic dynamics generate metastability |
| **RA2-05** | final G3 grade | Preserved: **CONDITIONAL DYNAMICAL HIERARCHY DERIVATION**; no frozen witness distinguished; no residual input eliminated; family / scaling priced; no consciousness claim |

**Applied to:** `G3_SPECTRAL_HIERARCHY.md`, `RA0_ZOOM_OUT_03.md`.

## G4 note (no repair)

- **Numerical defect fixed before the logged G4 run.** The subspace gap metric was first computed as sqrt(1 − σ_min²),
  which reads ~1e-6 for exactly coincident subspaces (cancellation). It was replaced by the direct operator norm
  ‖(I − E_A)E_V‖. The logged run uses the direct form.
- **Control added before the logged run.** The G3 ladder is translation-invariant, so its lane indicators are **exactly**
  lumpable for every β (an exact G1-type factor, not an emergent one). G4 therefore adds a site-modulated, non-identical-lane
  ladder that is not lumpable. Both are logged.
- **Floating-point floor.** For the exactly lumpable ladder, Δ_alg reads 1e-12 … 1e-9 (eigenvector round-off), slightly
  above the proved bound s·(4/√p_min + C_V) ~ 1e-14 … 1e-12. This is round-off on an exact zero, not a violation; every
  other logged bound holds.

## G5 note (owner ruling on `573898e`; no repair)

G4 is accepted provisionally as CONDITIONAL DYNAMICAL PARTITION DERIVATION. The owner's three tightenings are recorded
here as standing readings:

| ID | item | standing reading |
|---|---|---|
| **G5-R1** | G4 rank ≥ 3 rounding | Numerical only in G4. G5 targets the metastable recovery theorem directly, **not** CONJECTURE G4-C, which remains open |
| **G5-R2** | G4-T H3 | Load-bearing and **priced**. G4 did not show "diverging gap alone ⇒ partition algebra". It showed "nearly decoupled + reversible + regular block masses / eigenfunctions ⇒ partition algebra asymptotically". G5 weakens H3 to H3′ (p_min ≥ p_*, s²·C_V → 0) but does not remove it |
| **G5-R3** | "blind" | Separated into **mathematical canonicality** (the idempotent and primitive sets are fixed by V and π) and **algorithmic discovery** (seeded Newton or tensor-power starts are numerical solvers only) |

**Process notes, before the logged G5 run:**
- The first draft of the tensor bound (Lemma G5-2) mis-stated the third-order term as O(s³) and missed the first-order
  cancellation. The corrected bound, with first-order terms ≤ Λs² because A is an algebra and third order ≤ 2s²(C_V + Λ),
  is the one logged.
- The rounding-lemma constant was tightened from 4 to 2 (|d_i| + |d_j| ≥ 1 ⇒ d_i² + d_j² ≥ 1/2). Both versions are valid.
