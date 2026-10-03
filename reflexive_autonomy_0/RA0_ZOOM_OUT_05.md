# RA0 ZOOM-OUT 05 (after G5) — HARD STOP for owner review

> **Repaired by RA0 REPAIR 03** (RA3-01 … 05). Where wording differs, the ledger takes precedence.

**Branch:** `grut-reflexive-autonomy-0`. Not frozen. No consciousness, quantum, collapse or Born-rule work.

All "proved" items below are proved in `G5_METASTABLE_RECOVERY_THEOREM.md`. They are **not externally reviewed**.

## 1. Was tensor near-odeco structure proved?

**Yes (Lemmas G5-1 and G5-2).**
- The partition algebra's cubic moment tensor is exactly odeco: T_A = Σ p_i^{−1/2} u_i^{⊗3}.
- The slow space's tensor, transported by the canonical polar isometry W, satisfies ‖T̃ − T_A‖ ≤ s²(11Λ + 2C_V).
- The perturbation is **second order** in s, because the first-order terms cancel: A is an algebra.
- F1 confirms ε ∝ s².

## 2. Was blind component recovery proved?

**Yes, by a self-contained identifiability lemma (G5-3)**, not by upgrading an algorithmic result. The lemma uses
Newton–Kantorovich plus a global exclusion argument.

For ‖E‖ < 1/(36·K^{3/2}·Λ²):
- the compressed product has exactly 2^K idempotents;
- exactly K of them are primitive, by the canonical spectral count (one eigenvalue > 1/2);
- each is within (8/3)ε of a block indicator.

**Source audit.**
- AGHKT 2014 is algorithmic and randomized, so it is not used.
- Mu–Hsu–Goldfarb 2015 and Auddy–Yuan 2023 (deterministic) corroborate stability but are not used.
- Robeva 2016 is rederived.

## 3. Was a new margin assumption needed?

**Not for vanishing misclassified π-mass:** Chebyshev rounding, Lemma G5-4.

**Yes for exact eventual recovery:** a new supplied pointwise margin, H4 (sup-norm eigenfunction convergence). In F1 the
sup error falls 0.18 → 0.02, which is consistent with H4 but is not a proof.

## 4. What recovery notion was earned?

**Earned:**
- projector convergence;
- L² convergence of the canonical primitive idempotents;
- misclassified-or-tied π-mass ≤ (256/9)Kε² + 8s² → 0.

**Not earned:** exact block equality.

## 5. Did rank 2 emerge as a special case?

**Yes, exactly (PROP G5-R2).** For every 2-dimensional unital V, the idempotents are {0, 1, e₊, e₋} in closed form, and the
G5 argmax partition is identical to G4-P's {f > s/2}. Verified to 1e-14.

## 6. Did any result require hidden labels computationally?

**No.**
- The recovery map R uses only V and π.
- The hidden partition appears only in proofs and in the post-hoc audit, which runs after R returns.
- Seeded Newton starts are a solver for a fixed polynomial system. The theory does not depend on them, and the audit
  shows that the numerical solver found 2^k distinct idempotents. That is evidence, not proof of completeness [RA3-02].

## 7. Is the G4-T partition derivation now theorem-grade?

**Yes, in the vanishing-mass sense: CONDITIONAL DYNAMICAL PARTITION DERIVATION — THEOREM-GRADE IN G4-T CLASS.**

Conditions (priced):
- a supplied family that is reversible and irreducible;
- η/g → 0;
- p_min bounded below;
- s²C_V → 0 (H3 weakened);
- a divergence certificate that is a proof about the family.

**The bounds are loose.** [RA3-01] On F1 the numerical lower estimate ε_num falls below ε_glob from M = 80, but that does
not certify the lemma's condition, and no logged M is rigorously certified. The theorem's own sufficient
condition (via the proved ε bound) is only met at about M ≈ 400 by extrapolation, and much later if the analytic
Davis–Kahan bound on s is used.

## 8. Is the general G4-C conjecture still open?

**Yes.** G5 never uses Δ_alg; it bypasses the conjecture.

R alone is not a metastability detector: forced onto a diffusive ladder it returns a partition. Refusal of non-metastable
families still rests on the certified-cut rule.

## 9. Is a consciousness hypothesis sector now mathematically eligible?

**Under the owner's entry ruling, yes as a NEW HYPOTHESIS sector, with a narrow starting point:**
- Within class G4-T, an endogenous macro-algebra and macro-partition exist (up to vanishing π-mass) **before** any
  observer or subject is specified.
- The partition is derived from P_N with no supplied Π, K, ε, objective or labels.

**What any such sector must still respect:**
- The partition is defined up to vanishing π-mass (exact only under H4).
- It is a **macrostate** partition, not a subsystem split.
- The family, reversibility and regularity remain supplied.
- G5 itself establishes no awareness, experience, self-reference, memory or access.

The new sector must not take over A_partition's job. G5 gives the macro-partition, so consciousness does not need to
supply one.

## Scorecard

| gate | terminal |
|---|---|
| G1 | FINITE REFLEXIVE AUTONOMY NO-GO |
| G2 | CONDITIONAL MACROVARIABLE DERIVATION FROM A SUPPLIED FAMILY |
| G3 | CONDITIONAL DYNAMICAL HIERARCHY DERIVATION |
| G4 | CONDITIONAL DYNAMICAL PARTITION DERIVATION (provisionally accepted) |
| G5 | **METASTABLE BLIND RECOVERY PROVED** (vanishing π-mass; exact recovery needs margin H4) → **CONDITIONAL DYNAMICAL PARTITION DERIVATION — THEOREM-GRADE IN G4-T CLASS** |
| CONJECTURE G4-C | open |
| frozen witnesses distinguished | 0 |
| residual inputs eliminated (globally) | 0 |

**HARD STOP.** Returning for owner review. Not frozen. Consciousness not begun.
