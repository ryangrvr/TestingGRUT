# QFT-SCOUT-1 — MODULAR RESIDUAL COMPRESSION (charter)

**Status:** non-canonical playground campaign, authorized by owner ruling (this thread, after the Bridge-1 freeze).

| field | value |
|---|---|
| Branch | `qft-scout-1`, created from `grut-bridge-1-frozen @ 753b90a3e3976ed286f385ca1b455e04c92a2d26` |
| Frozen inputs (read-only) | `scout-2-reviewed @ 1f08ae1`; `grut-bridge-1-frozen @ 753b90a`; canonical GRUT `GRUT-RAI/master-w25bu9 @ b935099` |
| Not modified | scout-0, scout-1, scout-2, scout-2-review, scout-2-reviewed, grut-bridge-1, grut-bridge-1-frozen, GRUT-RAI |

## 0. Central question

The frozen Bridge-1 residual is

    D_dyn ⊕ [Σ ⊗ (H_marginals, H_cross, H_epoch)] ⊕ [A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time],
    A_closure = f(D, A_seed; R_closure)

Does a materially enlarged premise envelope **narrow, condition, or compress** any frozen component, rather than merely
re-expressing it? The enlargement is: algebraic QFT, infinite-dimensional local algebras, and modular theory.

## 1. Premise change (what is new, priced)

| ID | new premise | price / status |
|---|---|---|
| **QP-1** | Local algebras are von Neumann algebras on an infinite-dimensional Hilbert space, possibly non-type-I | supplied |
| **QP-2** | A Haag–Kastler-type net O ↦ 𝒜(O) (isotony, locality / Einstein causality), **where used** | supplied. Using it for Σ questions is RELOCATION by definition |
| **QP-3** | A vacuum / standard state (cyclic and separating), with positive energy and covariance **where used** | supplied |
| **QP-4** | Specific structural properties: Haag duality, the split property / nuclearity, type III₁ | supplied **or** imported as theorems for named models (e.g. the free field). Each use is recorded |
| **QP-5** | A physical quantum lift (ħ, the complex structure) | supplied. Canonical GRUT's lift is NONUNIQUE-LIFT, so every result here is **auxiliary to canonical GRUT** unless the lift itself is selected |
| **QP-6** | Modular-position data (half-sided modular inclusions, geometric modular action), **G3 only** | supplied unless derived |

## 2. Gates, order and hard stops

| gate | target component | question | hard stop |
|---|---|---|---|
| **G1** | **H_cross \| Σ** | Is H_cross = 0 (a normal product state) **admissible** for the chosen localization? Compare sharp complementary localization with split-buffer localization | **HARD STOP after G1 → `QFT_ZOOM_OUT_01.md`** |
| G2 | orientation / H_epoch | Does modular flow (𝒜, ω) → σ_t^ω select time or orientation, beyond (a) the canonical flow and (b) the interpretive identification? | after G2 |
| G3 | Σ + Lorentz | Do algebras + state + modular-position relations reconstruct net / geometry / symmetry without supplying an equivalent amount? | after G3 |
| G4 | payoff | the B5 seven-criterion bar, unchanged | final |

## 3. Preregistered success hierarchy (stricter than Bridge-1)

    FORBIDDEN  <  CONSTRAINED-NONUNIQUE  <  CONDITIONALLY SELECTED  <  TRUE COMPRESSION

| level | meaning |
|---|---|
| **FORBIDDEN** | a member of a frozen component's class is shown inadmissible in a stated class (e.g. "H_cross = 0 inadmissible in the sharp normal-state class") |
| **CONSTRAINED-NONUNIQUE** | the admissible set shrinks, but many members remain |
| **CONDITIONALLY SELECTED** | a unique member is selected given stated supplied inputs (e.g. a net recovered from a state plus modular-position assumptions) |
| **TRUE COMPRESSION** | the target information is reconstructed **without supplying an equivalent amount of it** |

The Bridge-1 compression vocabulary (RELOCATION, RENAMING, GAUGE, NO BRIDGE, BLOCKED) still applies. Proving that zero
cross-correlation is impossible is **CONSTRAINED-NONUNIQUE (or FORBIDDEN-in-class), not a derivation of H_cross.**

## 4. Mandatory finite / type-I controls (one per gate)

| gate | control |
|---|---|
| G1 | ordinary tensor-factor product states (B(H₁) ⊗ B(H₂), always admissible and normal) **vs** sharp AQFT localization **vs** split-buffer localization. Plus a lattice free-field cutoff sequence a → 0, labelled **illustration, not proof** |
| G2 | finite-dimensional modular flow from (ρ, B(H)) (exists trivially: σ_t(x) = ρ^{it} x ρ^{−it}) **vs** the geometric Bisognano–Wichmann identification |
| G3 | what modular data can and cannot reconstruct in a finite matrix-algebra setup |

The campaign must state exactly **what infinite / type-III structure buys**. "Infinite-dimensional" is never itself a
selector.

## 5. Smuggling firewalls

1. **Supplied net.** If the local net or localization structure is supplied, every Σ result is RELOCATION.
2. **Priced inputs.** The vacuum, positive energy, Haag duality, split property / nuclearity, type, modular position and
   the KMS sign are supplied inputs unless derived. They are counted like R_closure and the Gibbs postulate in Bridge-1.
3. **Theorem vs interpretation.** (𝒜, ω) → σ_t^ω is theorem-grade (Tomita–Takesaki). "σ_t^ω = physical time / orientation"
   (the thermal-time hypothesis) is an **interpretive hypothesis**, never a derivation.
4. **Separate imports.** Reeh–Schlieder (the vacuum is cyclic and separating for local algebras under the standard AQFT
   hypotheses) is imported **separately** from any bounded-energy-state entanglement extension. Each is sourced and graded
   on its own.
5. **Not type III alone.** "Type III ⇒ no product state" is **not** chartered. The chartered question is whether a
   **normal product state exists for the chosen sharp complementary localization**. That requires stating Haag duality,
   factoriality and normality explicitly. Split-buffer localization is the contrasting case.
6. **Lift.** Every result is auxiliary to canonical GRUT (NONUNIQUE-LIFT) unless the lift is selected.
7. **Lattice controls.** Lattice numerics are illustrations of continuum statements. They never adjudicate theorem-level
   claims.

## 6. Source grading

| grade | meaning |
|---|---|
| **KNOWN-RESULT IMPORT — PRIMARY-SOURCE VERIFIED** | statement checked against the primary text |
| **KNOWN-RESULT IMPORT — SOURCE LOCATED, TEXT NOT RE-READ** | bibliographic record located (web search / owner pointer). Direct fetches from this environment are blocked by the egress proxy (arxiv.org and link.springer.com returned EGRESS_BLOCKED on 2026-10-02) |
| **BRIDGE / QFT-SCOUT PROPOSITION** | proved or sketched here, with its review status stated |

## 7. Payoff bar (carried unchanged from Bridge-1 B5-0)

A distinctive payoff requires all of: forced · residual-input invariant · convention invariant · non-trivial · distinctive
vs the ordinary comparison class · observable · falsifiable.

Baseline: **ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS.**

## 8. Rules carried

- No PR unless asked. No canonical edit.
- Commit trailers as in Bridge-1. No model identifiers in commits or artifacts.
- Numerics are an **independent code path, not an independent reviewer.**
