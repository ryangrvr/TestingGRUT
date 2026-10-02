# INDEPENDENT REVIEW LEDGER

## Theorem-level grades (owner's independent hostile review, first pass)

| item | grade |
|---|---|
| ΣH-0 Prop. 1 (state productness selects no TPS) | **PASS** |
| ΣH-0 dimension count (REPAIR 06: codim 2N − 2 − 2Σ(d_i − 1)) | **PASS** |
| ΣH-0 Prop. 2, generic exact incompatibility | **CONDITIONAL PASS** (IR-01) |
| ΣH-0 Prop. 2, Pareto dichotomy | **RECLASSIFIED: NUMERICAL / SCOUT CLASSIFICATION, not a theorem** (IR-02) |
| D0 Lemma 1 (recurrence) | **PASS** |
| D0 Theorem 1 | **PASS** |
| D0 Corollary 1 | **PASS** |
| D0 Corollary 1′ (finite-horizon concatenation) | **PASS** |
| D0 Corollary 2 (Σ-compatible time reversal) | **PASS** under its stated assumptions |
| D0 Proposition 2 (H-only covariant states are stationary) | **PASS** |
| D0 overall | **THEOREM CORE: PASS; EXAMPLE / SCOPE REPAIR NEEDED** (IR-03, IR-04) |

## Findings

| ID | target | finding | consequence | residual boundary affected? |
|---|---|---|---|---|
| **IR-01** | ΣH-0 Prop. 2; S2-Σ / ΣH CPR citations | CPR's generic local-TPS uniqueness is **conditional**: the general algebraic argument assumes at least one Hamiltonian with a unique local TPS in the locality class. The explicit numerical establishment is for a restricted translation-invariant nearest-neighbour qubit class. Extension to arbitrary finite (k, d, n) is presented cautiously. **Also:** if H has several genuine dual local TPS classes, ψ-productness *can* distinguish among them | Do not attribute an unconditional all-(k, d, n) theorem to CPR. "For a compatible pair the joint principle adds no selection power beyond H-locality" holds **only when H's physical local TPS class is already unique** | **no** |
| **IR-02** | ΣH-0 Prop. 2 dichotomy (b) | F_H ∩ F_ψ = ∅ shows only that no frame achieves exact locality *and* exact productness. It does **not** imply several Pareto-incomparable solutions: a single dominating compromise is mathematically possible | The Pareto nonuniqueness of incompatible pairs is a **numerical SCOUT result** (8-case S2-ΣH plus S2-Σ), not a theorem. Prop. 2 earns "exact joint compatibility is non-generic under stated conditions"; the numerics earn "the remaining choices form competing Pareto alternatives in the models tested" | **no** (it sharpens the theorem / numerics division) |
| **IR-03** | D0 "covered examples" list | "Rényi entropies" is too broad. S₀(ρ) = log rank ρ is discontinuous: ρ_ε = diag(1−ε, ε) → diag(1, 0) jumps from log 2 to 0 | Restrict to **continuous Rényi entropies (finite-dimensional, α > 0)**; exclude rank entropy and other discontinuous functionals explicitly | **no** |
| **IR-04** | D0 covered list + D0 Prop. 2 wording | (a) The **thresholded** S2-8 / D7 redundancy `I(S:F) >= (1-δ) H(S)` is an integer count, generally **discontinuous** at threshold crossings, so D0 does not directly cover it. D0 does cover the continuous mutual informations it is built from. Record reversal was shown separately and numerically (D7: redundancy 6 → 0). (b) "The minimal case is Σ" is not proved. Σ is a **demonstrated sufficient** additional covariant structure, not a proved minimal one | Remove the thresholded redundancy from the covered list (keep the continuous MI); demote "minimal" to "sufficient (minimality not proved)" | **no** |

| **IR-05** | S2-ΣH H1 "conflict" case (found by numerical reproduction) | The original conflict state is exactly product across the id-frame 2-cut (01235)\|(4). Its non-dominance came from an intra-frame cut trade-off, not a clean H-local-frame vs ψ-product-frame conflict. Weak-conflict instances can be coarse-compatible and then show dominance (correctly) | Erratum E-06: describe H1 as coarse-grouping-compatible with an inner-cut trade-off. The intended conflict claim is **reproduced 3/3 with certified-incompatible pairs** (NR-ΣH-2-certified). Compatibility must be checked at every grouping | **no** (it strengthens the Σ-relativity of compatibility) |

| **IR-06** | frozen intermediate ledgers (owner's second-reader pass) | Several frozen ledgers deliberately keep historical classifications that were later repaired: (a) S2-Σ commutant copies read as "degeneracy-protected" (superseded by Y-07: H-relative gauge); (b) the H2 row ends at H_env (superseded by D-arrow: H_corr\|Σ); (c) "G2 / G3 / G4 not yet probed" (historical; the S2-G campaign was run); (d) JOINT_BOUNDARY H1 "incompatible" (superseded by IR-05) | **Reviewed precedence:** (1) explicit repair / correction, (2) final frozen handoff, (3) independent-review errata, (4) historical intermediate ledger row. The review handoff uses `review/REVIEW_ACCOUNTING_LEDGER.md`, never a frozen ledger row directly. Frozen ledgers are not rewritten | **no** |
| **IR-07** | `ledgers/ARROW_ORIGIN_LEDGER.md` AR-01 | AR-01 lists the record arrow as "NO" under D0's Theorem 1 / Cor. 1, which reads as extending D0 to the thresholded integer redundancy R_δ. D0 covers continuous subsystem entropies, continuous mutual informations and continuous fixed-POVM coarse observables. It does **not** cover rank entropy, thresholded R_δ or other discontinuous counts | AR-01's record entry should read "continuous MI: no every-state arrow (D0); R_δ: numerically reversed (D7; NR-D7), not theorem-covered". Companion to IR-04 | **no** |

**Net effect of IR-01 … IR-07:** the final residual `C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res` is **unchanged**. The
theorem package is review-clean only after the errata in `review/ERRATA_PROPOSED.md`.

## Numerical reproduction queue (owner-specified order)

| ID | target | status (`review/NR_SIGMAH_RESULT.md`) |
|---|---|---|
| NR-ΣH-1 | positive compatible case | **REPRODUCED** 3/3 |
| NR-ΣH-2 | H1 conflict | as constructed: **not reproduced → IR-05**; certified-incompatible: **REPRODUCED** 3/3 |
| NR-ΣH-3 | Haar case | **REPRODUCED** 3/3 |
| NR-ΣH-4 | local-dimension tie | **REPRODUCED** 3/3 |
| NR-ΣH-5 | epoch covariance | **REPRODUCED** 3/3 |

**IR-05 lesson:** compatibility must be tested over the complete candidate grouping class, not merely at the finest TPS
or the nominal frame. Applied to every SCOUT-2 case classified "incompatible" (`review/ir06_check.log`):
- H5a: min C_min over all 202 id-frame groupings = 1.46 bits;
- H7: 1.64 bits.

Both are certified. **No IR-06.**

**S2-Σ (NR-Σ-1 … 7):** REVIEW-CONFIRMED-WITH-SCOPE (`review/NR_SIGMA_RESULT.md`).

Next: the D-arrow key numbers, then H2.
