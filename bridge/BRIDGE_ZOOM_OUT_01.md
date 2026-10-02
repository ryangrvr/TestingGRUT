# BRIDGE ZOOM-OUT 01 (after B0 + B1) — HARD STOP for owner review

**Branch:** `grut-bridge-1` (from `scout-2-reviewed @ 1f08ae1`).
**Canonical source:** GRUT-RAI `master-w25bu9 @ b935099`, read-only and untouched.

## 1. Does GRUT independently reconstruct Σ?

**No.**
- With the local net deleted, the earned linear generator K does **not** reconstruct it:
  - The sparsest / most factorized frame is the decoupled normal-mode frame, which has no edges.
  - The connected-chain frames form a continuum.
  - Passivity leaves a continuum of inequivalent nets that all pass the geometry predicate.
  - A 2D-grid operator admits a passive 1D chain reading.
  - Inside the L0-1b class itself there is an explicit isospectral, inequivalent net at N = 12 whose geometry predicate fails.
- Unique reconstruction occurs only in the narrow C1-a class (unit springs; DLS). That class *is* a supplied net.

## 2. Is the canonical local net redundant, consistency-checked or still supplied?

> **Σ REMAINS SUPPLIED** in the earned linear core.
> - It is **CONSISTENCY-CHECKED / OVERDETERMINED** where on-site quartic drift (L0-1c) or non-uniform site noise (L0-1e)
>   is also supplied: the net is exactly recoverable from those layers.
> - That is **RELOCATION / REDUNDANT SUPPLY**, not elimination. Those layers were written site-locally, and a control
>   quartic recovers whatever frame it was written in.

Proposed: **CANONICAL UPDATE CANDIDATE B1-CUC-1** (`B1_SIGMA_LOCAL_NET_RESULT.md`). It sharpens EA0_OWNER_RULING_02 §3.
Not applied.

## 3. Did any GRUT-specific structure do work that generic SCOUT dynamics did not?

**Narrowly, and only as pricing or relocation:**
- **Passivity (pins ≥ 0) and the uniform-pin L0-1b form** shrink the family of admissible nets from the whole sphere to a
  neighbourhood. They do not reach uniqueness (counterexamples in both classes).
- **Unit springs (C1-a)** give uniqueness. This is CRITERION-PRICED.
- **Site-local nonlinear drift / non-uniform noise** select the net. This is RELOCATION.
- **Negative:** GRUT's classical direct-sum locality is a *weaker* object than SCOUT's tensor factorization (a frame plus
  weights, not a TPS), and it is still not reconstructed.
- **No GRUT-specific structure produced TRUE COMPRESSION.**

## 4. Did the reviewed residual shrink?

**No.**

    C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res

The residual is unchanged.
- **Gained:** an identification, not a reduction. Σ ↔ GRUT local net (OVERLAPPING, category difference), and D_dyn ↔ GRUT
  generator (IDENTICAL on "not derived").
- **Gained:** a bookkeeping redundancy. In canonical GRUT the net is supplied **more than once** (substrate + drift +
  possibly noise).
- **Not gained:** the time orientation is still not selected, and Lorentz structure is still not derived. Bridge-1 has not
  touched either.

Ledger tally: TRUE COMPRESSION 0.

## 5. Which bridge next: A_res (B3) or H_corr (B2)?

**Recommend B3 (A_res vs EA-0) first**, in line with the owner's order. Reasons:
1. GRUT's richest relevant structure is on the access side: EA-0 L7 (non-trivial access only via exact local invariant
   structure), Theorem S, GS1 (site-resolved access selects geometry), and the open `L0_ACCESS_BRIDGE_01` Route A/B/C
   taxonomy. H_corr has a thinner GRUT counterpart: the "preparation-relative" S6, which is already SCOUT STRICTLY SHARPER
   at B0.
2. B1 suggests a concrete B3 hypothesis. EA-0 L7 ties non-trivial access to *local* invariant structure, and B1 shows the
   net is relocated into drift and noise. So GRUT access is likely **RELOCATION of Σ**, not compression of A_res. That is a
   sharp test.
3. B2 risks pure **RENAMING** (preparation ≡ H_corr\|Σ) unless the GRUT environment layer can be shown to carry it. It is
   better run after B3 fixes how the access is factored.

**B3 guardrails:**
- No silent lift (EA-0 is UNFORMULABLE at earned scope; quantum results are auxiliary only).
- No EA-1.
- Do not optimize for a unique lift.
- Record a Route C (category boundary) outcome as a result.

## Status

- **B0:** DONE. **B1:** DONE.
- **TRUE COMPRESSION:** 0. Residual unchanged. One CANONICAL UPDATE CANDIDATE proposed (B1-CUC-1).
- **Empirical status:** zero distinctive GRUT predictions; B5 not reached.
- **HARD STOP. Awaiting owner review before B3.**
