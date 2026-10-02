# SCOUT-2 REVIEW STATUS (live; `scout-2-review`)

**Codes:**
- REVIEW-CONFIRMED;
- REVIEW-CONFIRMED-WITH-SCOPE;
- ERRATUM-REQUIRED;
- NOT-REPRODUCED;
- UNREVIEWED.

The frozen `scout-2` (`af0042f`) is never edited. Errata live in `ERRATA_PROPOSED.md`.

**Numerical reproductions are an *independent code path, not an independent reviewer*.**

## 1. Theorem review (owner's independent review)

| frozen item | status | refs |
|---|---|---|
| ΣH-0 Prop. 1 (productness selects no TPS) | REVIEW-CONFIRMED | — |
| ΣH-0 dimension count (REPAIR 06) | REVIEW-CONFIRMED | Y-13 |
| ΣH-0 Prop. 2 (compatibility non-generic) | ERRATUM-REQUIRED (conditional on CPR-type uniqueness; dual-class caveat) | IR-01, E-01 |
| ΣH-0 Prop. 2 Pareto dichotomy | ERRATUM-REQUIRED (reclassify as numerical) | IR-02, E-02 |
| D0 Lemma 1 / Thm 1 / Cor 1 / Cor 1′ / Cor 2 / Prop 2 | REVIEW-CONFIRMED (Cor 2 under its stated assumptions) | — |
| D0 covered-examples list | ERRATUM-REQUIRED (Rényi-0; thresholded redundancy) | IR-03, IR-04, E-03, E-04 |
| D0 Prop. 2 "minimal case is Σ" | ERRATUM-REQUIRED (sufficient, not proved minimal) | IR-04, E-04 |

## 2. Numerical reproduction

| frozen claim set | status | file |
|---|---|---|
| S2-ΣH: positive compatible case | REVIEW-CONFIRMED | NR_SIGMAH_RESULT.md |
| S2-ΣH: H1 conflict as characterized | ERRATUM-REQUIRED (construction flaw: coarse-compatible) | IR-05, E-06 |
| S2-ΣH: conflict claim (certified incompatible) | REVIEW-CONFIRMED (re-established 3/3) | NR_SIGMAH_RESULT.md |
| S2-ΣH: Haar / local-d tie / epoch covariance | REVIEW-CONFIRMED | NR_SIGMAH_RESULT.md |
| S2-ΣH: H5a, H7 incompatibility (complete grouping test) | REVIEW-CONFIRMED (1.46 / 1.64 bits; no IR-06) | ir06_check.log |
| S2-Σ: the 8 load-bearing claims | REVIEW-CONFIRMED-WITH-SCOPE (δ-threshold stable window noted) | NR_SIGMA_RESULT.md |
| S2-D-arrow key numbers (D1, D3/D6, D4, D5, D7, D8, D2) | REVIEW-CONFIRMED (numerics; theorem errata IR-03 / IR-04 separate) | NR_DARROW_RESULT.md |
| S2-H2 key numbers (structural, relocation, alignment, stripes, invertibility) | REVIEW-CONFIRMED (H4 every-state claim rests on the majority argument, scope noted) | NR_H2_RESULT.md |
| S2-1 individuation | REVIEW-CONFIRMED-WITH-SCOPE (local ≠ global; IR-01) | NR_LOWER_RESULT.md |
| S2-1b locality criteria | REVIEW-CONFIRMED | NR_LOWER_RESULT.md |
| S2-2 composition | REVIEW-CONFIRMED-WITH-SCOPE (family discriminator) | NR_LOWER_RESULT.md |
| S2-3 / 3b probability and measure | REVIEW-CONFIRMED-WITH-SCOPE (compact theorem imported) | NR_LOWER_RESULT.md |
| S2-4 convexity | REVIEW-CONFIRMED | NR_LOWER_RESULT.md |
| S2-5 / 6 shared reference | REVIEW-CONFIRMED-WITH-SCOPE (Y-05, Y-06) | NR_LOWER_RESULT.md |
| S2-7 histories | REVIEW-CONFIRMED | NR_LOWER_RESULT.md |
| S2-8 Darwinism | REVIEW-CONFIRMED (R_δ numerical, IR-07) | NR_LOWER_RESULT.md |
| S2-G2 / G3 / G4 | REVIEW-CONFIRMED (G3 with REPAIR 05 scope) | NR_LOWER_RESULT.md |

## 3. Information-accounting review

| item | status |
|---|---|
| ledgers (accounting, measure, arrow, dimension, joint) | AUDITED → `REVIEW_ACCOUNTING_LEDGER.md` (IR-06 precedence, IR-07) |
| further decomposition of A_res | UNREVIEWED |
| IR-05 lesson | **compatibility must be tested over the complete candidate grouping class**, not just the finest TPS or the nominal frame. Applied to every SCOUT-2 case classified "incompatible": only H1 failed |

## 4. Unresolved findings

| ID | status |
|---|---|
| IR-01 … IR-07 | recorded; errata E-01 … E-08 proposed; none affects the residual boundary |
| open | none threatening `C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res` |
