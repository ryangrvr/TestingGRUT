# SCOUT-2 REVIEW HANDOFF — the hostile audit and reviewed interpretation

> **Branch:** `scout-2-review`, built on the frozen SCOUT-2 head `af0042f`. The frozen record (`scout-2`) is unchanged.
>
> **This document, not the frozen branch, is the reviewed interpretation of SCOUT-2.**
>
> **Review layers:**
> - **Theorem review:** the owner's independent hostile review (IR-01 … IR-04).
> - **Ledger audit:** the owner's second-reader pass (IR-06, IR-07).
> - **Numerical reproductions:** an **independent code path, not an independent reviewer** (same agent; no shared code;
>   different routines; fresh frames / seeds). This produced IR-05.
>
> **Interpretation precedence (IR-06):**
> 1. explicit repair / correction;
> 2. final frozen handoff;
> 3. independent-review errata (`review/ERRATA_PROPOSED.md`);
> 4. historical intermediate ledger rows.
>
> The reviewed accounting is `review/REVIEW_ACCOUNTING_LEDGER.md`.
>
> **No `scout-2-reviewed` branch has been created.** That decision is brought to the owner (Q9).

---

## Q1. Which frozen claims are REVIEW-CONFIRMED?

| claim set | evidence |
|---|---|
| D0 Lemma 1, Theorem 1, Corollaries 1 and 1′, Proposition 2 (theorem core); Corollary 2 under its stated Σ-compatible time-reversal assumptions | owner theorem review |
| ΣH-0 Proposition 1 and the repaired dimension count (Y-13) | owner theorem review |
| S2-ΣH: positive compatible case; Haar no-dominance with a trivial correlation-first winner; local-dimension tie; epoch covariance (C₀ loses dominance; C_min class unchanged, t* shifted) | NR-ΣH, 3/3 seeds |
| S2-ΣH: the **conflict claim**, re-established with certified-incompatible pairs | NR-ΣH-2-certified, 3/3 |
| S2-D-arrow numerics (D1 mirror / anti-arrow; D3 / D6 correlations vs marginals; D4 maximally mixed fresh bath; D5 reuse; D7 record reversal; D8 Σ-dependence; D2 typicality) | NR-D |
| S2-H2 numerics (structural controls; relocation; all-to-all vs ring; frozen stripes; exact invertibility) | NR-H |
| S2-1b, S2-4, S2-7, S2-8, S2-G2, S2-G4 | NR-L / NR-G |

## Q2. Which are CONFIRMED-WITH-SCOPE?

| claim | scope carried |
|---|---|
| **S2-Σ** (all 8 load-bearing claims) | objective dependence shown for the *tested natural criteria*, not a proof that no fundamental objective exists; CPR-type selection needs (n, d, k) supplied; epoch statement only for exp(−iHs) (a nonlinear commutant element is *not* a time shift in one exploratory instance; the full commutant is not adjudicated); scale-pricing holds for τ and ε, with a **stable window for the record threshold δ** |
| **S2-1** | local ≠ global uniqueness; CPR conditional (IR-01) |
| **S2-2** | FAMILY COMPOSITION DISCRIMINATOR; GMCD scoped to boxworld |
| **S2-3 / 3b** | the compact-quotient unique ergodicity is imported (Furstenberg), not simulated |
| **S2-5 / 6** | one J⊗J obstruction (not full LT); asymptotic coherence (not "iff") |
| **S2-G3** | MM within its calibration domain only (REPAIR 05) |
| **H2 complete-graph Ising** | the every-state claim rests on the odd-N majority argument (checked), not on simulating every state |

## Q3. Which require errata? (`review/ERRATA_PROPOSED.md`; frozen files untouched)

| erratum | finding | correction |
|---|---|---|
| **E-01** | IR-01 | ΣH-0 Prop. 2 is conditional on CPR-type uniqueness of H's local class; dual local classes can be distinguished by ψ |
| **E-02** | IR-02 | ΣH-0 Prop. 2's Pareto dichotomy is a numerical classification, not a theorem |
| **E-03** | IR-03 | D0 covers continuous Rényi entropies (α > 0), not rank entropy |
| **E-04** | IR-04 | D0 covers continuous MI, not thresholded R_δ; Σ is a *sufficient*, not proved minimal, extra structure |
| **E-05** | — | handoff wording updates following E-01 … E-04 |
| **E-06** | IR-05 | the S2-ΣH H1 case was coarse-compatible (product across (01235)\|(4)); its non-dominance was an intra-frame cut trade-off |
| **E-07** | IR-07 | ARROW_ORIGIN_LEDGER AR-01 record cell |
| **E-08** | IR-06 | precedence rule for superseded ledger rows (no file edits) |

## Q4. Which were not reproduced?

**Only the S2-ΣH H1 conflict case as originally characterized** (IR-05). The intended incompatibility result was
**re-established on a separate code path** with certified-incompatible pairs, 3/3. Every other targeted load-bearing
conclusion reproduced.

## Q5. Did any failure affect the residual boundary?

**No.** IR-01 … IR-07 are theorem-scope, example-list, bookkeeping or construction issues. IR-05 if anything
*strengthens* the Σ ⊗ H_corr coupling: compatibility is itself grouping-relative, so it must be tested over the complete
grouping class. That lesson was applied to every "incompatible" case; only H1 failed, so no further incompatibility misclassification was found.

## Q6. Reviewed final residual

    C5  →  D_dyn  ⊕  [Σ ⊗ H_corr]_coupled  ⊕  A_res

- **Time orientation:** NOT selected.
- **Lorentz structure:** NOT derived.

| component | reviewed reading |
|---|---|
| **D_dyn** | the Hamiltonian / transition law, connectivity and propagation rule. Graph / spectral / causal-order statistics, rigid-D invariant measures and structural attractors are **derived from** D (conditional derivations). D itself is not derived |
| **[Σ ⊗ H_corr]_coupled** | the state price is **H_corr\|Σ**: independence / freshness / low inter-subsystem correlation relative to a factorization, at a special moment. It is **not H_env and not low entropy**. Σ's physical non-uniqueness rests on **reproduced objective / dimension / scale / state conflicts**, not on commutant copies (which are H-relative gauge) |
| **A_res** | coupled to Σ but not reducible to Σ or D in any tested case: fragment grouping, record threshold, readout basis, history times, resolution, effect restrictions |

## Q7. What remains genuinely unreviewed?

- an **independent human or second-model reviewer** for the numerics (the reproductions are an independent code path by
  the same agent);
- **further decomposition of A_res**;
- **re-fetching the primary sources from this environment** (arXiv / APS blocked; literature upgrades rest on the owner's
  checks; `review/LITERATURE_REVIEW_OVERLAY.md`);
- **exhaustive historical-number reproduction** (by design only load-bearing conclusions were graded);
- **the full Stoica theorem** and **the full-commutant epoch question** (not adjudicated);
- **out-of-envelope questions** (S2-D-arrow-∞; QFT / algebraic subsystems), out of scope by construction.

## Q8. Is SCOUT-2 still saturated after hostile review?

**Yes.** The hostile review found:
- one construction flaw (IR-05);
- four theorem / example scope issues (IR-01 … IR-04);
- two bookkeeping issues (IR-06, IR-07).

Each was repaired or re-established, and **the residual boundary survived every correction**. No reproduction suggested
a hostile test that would move the conceptual boundary inside the finite-dimensional premise envelope.

> **SCOUT-2 SCIENTIFICALLY SATURATED AT CURRENT PREMISE ENVELOPE — HOSTILE REVIEW COMPLETED (independent code path +
> owner theorem / ledger review); RESIDUAL BOUNDARY CONFIRMED.**

## Q9. Should a `scout-2-reviewed` immutable branch be created?

**Recommended: yes, after the owner reads this handoff.** Create it from the current `scout-2-review` head, as a
citation target. Preserve three objects:

| branch | role |
|---|---|
| `scout-2` (`af0042f`) | the frozen experiment record (never edited) |
| `scout-2-review` | the hostile audit (may receive later review additions) |
| `scout-2-reviewed` | an immutable snapshot of the reviewed interpretation |

**Do not merge review into `scout-2`.** Creation is left to the owner.

---

**Index** (`review/`):
- `README.md`;
- `REVIEW_STATUS.md`;
- `INDEPENDENT_REVIEW_LEDGER.md` (IR-01 … IR-07, grades);
- `ERRATA_PROPOSED.md` (E-01 … E-08);
- `REVIEW_ACCOUNTING_LEDGER.md`;
- `LITERATURE_REVIEW_OVERLAY.md`;
- results: `NR_SIGMAH_RESULT.md`, `NR_SIGMA_RESULT.md`, `NR_DARROW_RESULT.md`, `NR_H2_RESULT.md`, `NR_LOWER_RESULT.md`;
- scripts + logs: `nr_*.py`, `ir06_check.py`, `nr2_*.py`.

**Empirical status (unchanged):** zero distinctive GRUT empirical predictions.
