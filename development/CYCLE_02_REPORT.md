# Issue #2 — engineering cycle 2, review package

**ENGINEERING BLOCKED. SCIENTIFIC REVIEW_REQUIRED. CARD 1 CLOSED.**

Seven concrete repairs are demonstrated in a disposable checkout. They are supplied as an exact review-only patch, not applied to protected inputs. The four undeclared pointer cases and two orphaned OPEN passes remain active obligations. The 25 bank flags (24 flagged changes and the deletion of `rung1_inin_action`) are unchanged. No foundations lane, candidate law, scoring, pricing, card or IP-12 draft was opened.

This work follows the owner’s stop-foundations direction in issue #2. PR #1 remains draft, targeting `scout-0`; no canonical merge or independent approval is claimed.

## Demonstrated repair

| Defect | Repair in the disposable checkout | Validation |
|---|---|---|
| Two “valid” classifier fixture tests fail on the live invalid repository | Detach a copied declaration table, exact synthetic enumerators and a minimal OPEN-pass inventory; restore original globals afterward | Valid fixtures pass even when live OPEN passes are contaminated; unmodified live classifier still rejects all four cases/two orphans |
| CLOSED/MOOT mutations target uncited P1A and never exercise their stated condition | Mutate actually cited P4 in the isolated fixture | Both mutations are refused for NON-OPEN PASS; full mutation unit battery passes |
| Calculation manifest drift | Regenerate strictly from tracked `calc/*.py` | 149 tracked calculation files; no mutant or untracked file enters the manifest |
| Public generated counts drift | Regenerate PUBLIC_NUMBERS and its dependent public document with an explicit proposed wave stamp | 22 test files, 242 tests; generated output matches register and source |
| Stamped test count 240 versus 242 | Change the stamped constant to the exact collection count | All 242 original node identities remain; exactly seven FAIL→PASS changes |

The full overlay is **239 PASS / 2 FAIL / 1 SKIP**. Every other outcome equals the frozen manifest. The two remaining failing test nodes are the original sealed-numeric and sealed-pointer guards. Both declaration bytes and every preregistration/result/seal manifest remain identical.

The public update proposes fixed STAMP_DATE/WAVE_DATE `2026-10-10`; historical CORRECTION_DATE `2026-08-14` is preserved. These are proposed public wave changes, requiring review rather than silent issuance. The exact patch has SHA256 `4c46efcb263d738d1909474cecb5905ecd1c94db5a6da11885ee08c0c2bf8ba1`.

Files: [repair.patch](review_cycle_02/repair.patch), [overlay_results.json](review_cycle_02/overlay_results.json), [remaining_owner_items.json](review_cycle_02/remaining_owner_items.json).

Reproduce without editing protected inputs:

```bash
python development/reconciliation_review.py --output development/results/review-overlay --stamp 2026-10-10
```

The script verifies the locked protected input digest before cloning, reproduces the seven failures, builds the patch, compares all 242 outcomes, runs the unchanged live classifier, and checks declaration/pass/seal invariance. It refuses source drift, unexpected changed paths, missing tests or any other outcome change. A dedicated hosted repair-overlay job retains the patch, raw JUnit and classifier rejection alongside the unchanged default job.

The first local overlay correctly rejected an unexpected outcome: same-length stamp edits in the same timestamp tick reused pre-repair bytecode. The review runner now disables bytecode writes in the disposable checkout. The successful rerun proves the repaired source rather than its stale cache. No guard was loosened.

## Owner dispositions still required

1. **Four pointer cases:** both TERMINATION V3 and V4 point to RESULT_TERMINATION_events, where line 312 reports “6.24 -> 4.87 sigma”. Each seal triggers two existing numeric-leak patterns. These are actual reported figures, not a frozen threshold false positive. P2’s earlier architectural defense was recorded for identifiers; a ruling must explicitly address both seals and these numeric cases. Keep all four undeclared reds until that ruling. P4’s separate KAPPA-results question is not silently folded into P2.
2. **P1A:** zero tier cases does not establish its stated MOOT condition. The old rung1 ID is absent, but `rung1_inin_formalism`, `rung1_ontology_finite_memory` and `rung2_kms_gate` still depend on the booked, assumed `background_time_translation_flow`. The stated discharge did not happen. Do not manufacture CLOSED/MOOT or add a symptomless exception. Record an owner disposition with its actual cause.
3. **P6:** the stale-net enumerator is empty, but no owner closure follows from that test. Its record says “CLOSES WHEN: the owner rules on whether to amend the four documents, and the amendment lands.” The current document hashes are provided; record the amendment’s authority or another explicit disposition.
4. **Bank inventory:** retain all 24 changed-node flags plus the deletion. The record provides no independent approval for those changes.

These are applications of the existing OPEN_PASSES rules, not new gates. That file states “Closing a pass is a HUMAN ACT” and distinguishes a ruling (CLOSED) from dissolution (MOOT). It forbids treating a green guard as the ruling itself.

## CR-5 separate implementation check

The private bundle is supplied out of band; no private Stage-3 source is copied into TestingGRUT. The harness pins the complete bundle, verifies its member manifest and all 20 payload hashes before loading code, checks the issued ruling hash, and cross-checks the CR-5 changed-file manifest. It writes the separate oracle/control checkpoint before importing the supplied comparison implementation.

**183 exact comparisons agree**: interval boundaries and ambiguity, exact/oscillation/dyadic error arithmetic, context/pipeline/point binding changes, missing receipts/current images, NR-17 all-removed fallback, Q3’s independent two chains, conservative linear transport, the shared full-K NR-4 mask, missing ablation images retained in the denominator, responsibility, and VOID fractions.

A hostile invalid-input fixture supplies reference weights `[-1, 2]` and a hit only at the second point. The supplied appearance routine accepts them and returns fraction `2`; the separate implementation rejects negative measure weights. This is an **input-contract validation gap**, not a counterexample to certified arithmetic with valid nonnegative weights. It must be addressed in the reviewed caller/kit contract. No automatic reconciliation or private kit edit was made.

This implementation separation is a builder check, not independent external authorship. Synthetic registry entries stand for already-proved enclosures; they do not prove a physical limit, certify a carrier or establish a lock. The bundle reports canonical publication/re-freeze false and external approval false. The broader evaluator source is not included and is not reviewed by this harness.

Files: [cr5_comparison.json](review_cycle_02/cr5_comparison.json), [cr5_independent_ready.json](review_cycle_02/cr5_independent_ready.json).

```bash
python -m development.cr5_review.check --bundle /absolute/path/GRUT_G2_14R_CR5_REVIEW_BUNDLE.zip --output development/results/cr5-review
```

## Reporter and environment repairs

The development checker now refuses duplicated or contradictory classifier verdicts, even when the expected substring and diagnostics are present. Hostile regression tests cover both red and green transcript contradictions. No production classifier declaration was altered.

The first full-mutation execution failed on an unmutated control because SymPy was absent. Its failed artifact is retained. Separate expensive-profile numerical dependencies are now pinned (mpmath 1.3.0, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0), checked for exact installed versions and reported. The five-package fast-loop requirements and frozen observational manifest remain untouched.

## Execution evidence and readiness

| Run | Raw outcome | Contract result | Scope |
|---|---|---|---|
| Development regression tests | 34 PASS | PASS | Reporter, source/fixture isolation, hash checks and exact arithmetic |
| Protected default suite + actual classifier | 232 PASS / 9 FAIL / 1 SKIP | Exact baseline reproduced, classifier agreement true | Engineering BLOCKED; scientific REVIEW_REQUIRED |
| Repair-only overlay | 239 PASS / 2 FAIL / 1 SKIP | Exactly seven repaired nodes, no other changes | Proposed patch; four cases/two orphans still rejected |
| Full mutation, initial environment | 8 PASS / 1 FAIL | FAIL | Missing SymPy; unmutated control refused; failed evidence retained |
| Full mutation, pinned numerical environment | 9 PASS | PASS, zero deltas/errors | All 14 controls and 55 registered mutants executed |
| Slow falsifiers, pinned numerical environment | 233 PASS / 9 FAIL | PASS, zero deltas/errors | All 15 registered cited calculations exit zero; former skip enabled |
| Separate CR-5 implementation | 183 agreeing comparisons | Valid-input arithmetic PASS | Signed-weight validation finding remains; no external approval |

Successful numerical guard profiles do not replace the default classifier, resolve a bank flag or certify new physics. Default and both expensive profiles verify unchanged protected inputs and retain the existing 7/4/2 adjudication inventory. The default case digest remains `fc6c0495455622e495d16dcb48bbbdd1ba2fb3622ffd5ca102e8ab96013cfc42`; protected source digest remains `efc48249a1657d5cd0e6bd3bde1aae0a653ceb1e6e78f26b23cf95184da851e1` across 1,757 files. The observational manifest and its hash lock are unchanged.

The machine records identify the actual execution checkout `b5d0493a23e5b741397f4dd3867bd8a09da08858` and its dirty development inputs honestly; they are not relabelled as the later publication commit. Hosted verification of the publication is recorded on PR #1.

## Required next act

Review the concrete five-file repair patch and proposed public wave update, then record the existing pointer/pass dispositions and independently check the evaluator/CR-5 implementation. A visible, reviewed observational-manifest transition follows any accepted protected change; this work never performs one automatically.

There is no engineering clearance, scientific approval or Card-1 authorization in this package. It creates no further foundations prerequisite lane. Card budget consumed: **0**.

