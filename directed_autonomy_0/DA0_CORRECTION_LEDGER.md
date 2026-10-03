# DA0 CORRECTION LEDGER (additive)

**Rule:** logs are kept as emitted. Where document wording differs, this ledger takes precedence.

## DA0 CHARTER REPAIR 01 (owner review of `bc9d3d2`; applied before any computation)

| ID | item | correction |
|---|---|---|
| **DC1-01** | cycle decomposition is not canonical | J, the cycle space and edge affinities are canonical. Any particular fundamental-cycle decomposition is **BASIS-PRICED / DIAGNOSTIC ONLY**, unless a basis-independent selector is proved. A "canonical cycle" must be invariantly identified from D; otherwise the terminal is **FLUX STRUCTURE ONLY — NO CANONICAL CYCLE**. EP = ½ Σ J log[(π_x G_xy)/(π_y G_yx)] ≥ 0 only for mutually supported edges; one-way edges need an explicit extended convention |
| **DC1-02** | real-part gap is not a sufficient non-normal cut certificate | Use r(λ) = −Re λ. Spectral sets must be conjugation-closed. Report E_N, ‖E_N‖, K_N(Γ) and conditioning. Perturbative stability needs the contour-resolvent product to vanish. Exploding conditioning ⇒ **SPECTRAL SEPARATION PRESENT — NONNORMAL RECOVERY UNSTABLE**, not a recovered macrostructure |
| **DC1-03** | orientation bit vs supplied non-reversible structure | P vs P* is one supplied bit. A and J are part of supplied D and carry more. The positive grade is **CONDITIONAL DERIVATION OF DIRECTED MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D**. Test: covariance under P → P*; matched families with the same π and S but different A may differ (**A-PRICED**). R1 becomes a matched pair |
| **DC1-04** | partition recovery and directedness are separate | C1-P (partition without supplied inputs) and C1-D (directed macro-generator, conditional on C1-P) are separate questions. C1-D without C1-P earns no directed endogenous partition. For Q vs Q*: covariant partition, sign-reversed current and affinities, no supplied cycle basis. Final information accounting: derived / supplied / not earned |

## C1 note (process; no repair)

These changes were made **before the logged C1 run**:
- `macro_cycle` crashed on a round-off-level macro current in R1b. A declared computational zero (ZERO_J = 1e-13) and a
  guard requiring a cycle of length ≥ 3 were added.
- The Kato bound was first coded in the crude form ρR₀·δ/(1 − δ), which is O(√(η/g)). It was replaced by the
  residue-expanded form 2eaη + (e + ρa)δ²/(1 − δ), which is O(η/g) and is the one stated in THEOREM C1-T. Both are valid.
- The macro-current error constant was corrected from 2·q_max·m to **4**·q_max·m, as in PROP C1-D. All logged m are 0.
- An audit mapping blind macro-cycle labels to hidden labels, the R1 matched-pair check, and the DC1-02 defect-ring probe
  were added.
- The ad-hoc defect-ring probe that preceded this was run once outside the log. Its values are reproduced in the logged
  probe.

## DA0 C1 REPAIR 01 (owner review of `feb6bc3`)

The C1 terminal is preserved: **NONREVERSIBLE PARTITION DERIVATION**, graded **CONDITIONAL DERIVATION OF DIRECTED
MACROSTRUCTURE FROM SUPPLIED NONREVERSIBLE D — A-PRICED**. The original log is preserved.

| ID | item | correction |
|---|---|---|
| **C1R-01** | rank-K proof in THEOREM C1-T | The ‖E − E₀‖ < 1 argument is withdrawn; the hypotheses do not make the explicit bound < 1 at finite N. Replaced by a **homotopy** G_t = G₀ + tG₁ on the fixed circle |z| = ρ. For all t, ℓ_t ≤ ℓ < ρ and δ_t ≤ δ < 1, so the circle stays in the resolvent set and the rank is constant: rank E₁ = rank E₀ = K. The projector bound is kept for convergence only |
| **C1R-02** | macro-current recovery condition | N4 is replaced by **N4a** (j_min,N = min nonzero \|J_𝔅\|) and **N4b** (q_max,N·m_N / j_min,N → 0). Sufficient specialisation when j_min ≥ c_J η: q_max·η/g² → 0. N3 alone (η/g → 0) does not suffice when g → 0, e.g. η = g^{3/2}. **Audit:** R1a, R2 and R5 satisfy N4a / N4b asymptotically (g = O(1), η ∝ 1/M); **R4 fails both** (j_min/η → 0; q_max·η/g² grows ∝ M) |
| **C1R-03** | zero currents | "Support and sign pattern are recovered" is withdrawn. Under N4b, non-degenerate nonzero entries have their signs and relative magnitudes recovered. Hidden zeros are shown only to converge to zero, and an exact zero is never inferred from a vanishing sequence. Structural zeros by construction (R2's non-adjacent basins; R5's pair (1, 3)) are family-specific facts |
| **C1R-04** | macro EP | General macro-EP recovery is removed from PROP C1-D (log traffic ratios can be ill-conditioned). Kept: blind EP_Q, the numerical values, EP(G*) = EP(G), EP_Q ≤ EP |
| **C1R-05** | time-reversal certification | "G* satisfies N1 – N4 with the same constants" is withdrawn: the decoupled operator of G* need not be G₀*. Kept: J(G*) = −J(G) exactly, and J_Q(G*) = −J_Q(G) for any fixed partition. If G and G* are each certified for the same asymptotic partition, the recovered partitions agree and the currents reverse. The numerical covariance is evidence |
| **C1R-06** | K2 comparator scope | The theorem-like wording is replaced by a comparator observation: in the ideal first-order Markov description with predictively distinct macro-rows, causal states identify the macrostates. Verdict unchanged: **NOT EXCEEDED** |
| **C1R-07** | multi-cycle control | **R5 added**: K = 4, complete macro graph (cycle-space dim 3), supplied superposition of two independent circulations. No spanning tree, cycle basis or decomposition is used. Partition recovered blindly (misclassified 0 at M = 10 … 80); canonical J_Q reported as a matrix; J_Q is not a multiple of any single cycle (5 support edges on 4 vertices); time-reversal covariance holds to ≤ 1.4e-16. Terminal: **NONREVERSIBLE PARTITION DERIVATION + CANONICAL MACRO FLUX; NO CANONICAL CYCLE DECOMPOSITION** |
| **C1R-08** | non-normal instability | Remains an **untriggered possible terminal**; its absence is not claimed in general. N2 stays priced. No pathological family is required before C2 |

**Applied to:** `C1_NONREVERSIBLE_MACROSTRUCTURE.md`, `DA0_ZOOM_OUT_01.md`, `DA0_STATUS.md`. New numerics:
`c1/c1_repair.py`, `c1/c1_repair.log`.

## Owner ruling on `1c3e9ad` (no repair)

- **C1 ACCEPTED AFTER REPAIR 01. C2 APPROVED TO CHARTER.**
- **Wording discipline (standing).** At the logged finite sizes, the N4b sufficient quantity q_max·η/g² is still > 1 for
  R1a and R5 and ≈ 1 for R2 at the largest M. Their certification is **asymptotic, from the analytic 1/M scaling**, not
  from having entered a numerically small-error regime.
