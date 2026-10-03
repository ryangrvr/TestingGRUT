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
