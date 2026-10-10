# Issue #2 owner direction

Source: https://github.com/ryangrvr/TestingGRUT/issues/2#issuecomment-6092113130

Owner: `ryangrvr`. Recorded 2026-10-09 America/Chicago (2026-10-10 UTC).
This is the verbatim GitHub direction below; it applies only to the development/review workbench.

+## OWNER DIRECTION — engineering reconciliation after cycle 2

This is an owner disposition for the **development/review workbench only**. It does not merge PR #1, accept bank flags, clear CR-5/evaluator review, alter canonical Stage-3 state, or authorize Card 1.

### 1. P2 scope — four V3/V4 event-log numeric pointer cases

**RULING: these four cases belong to the already-open P2-TERMINATION-EVENTLOG adjudication. They are not four new independent red classes. P2 remains OPEN.**

Reason:

- Original GRUT-RAI history at commit `2c522b20d21d3390903d3cb4b6a063d0a8fbbaee` shows V4 already manifested, while V3 explicitly states `DRAFT -- NEVER SEALED, RETIRED PERMANENTLY` and was not in the manifest at that checkpoint.
- At that same checkpoint, `RESULT_TERMINATION_events.txt` contains no `6.24 -> 4.87 sigma` entry.
- Therefore the current target contents do **not** establish that those sigma figures were available through the pointer at V4's seal point, and they cannot by themselves prove a historical pre-seal information leak.
- V3's later manifestation as immutable history must not silently turn an explicitly unsealed/retired draft into an in-force sealed instrument.
- Separately, the existing P2 policy question remains real: a blind-safe artifact that points to a mutable event log may become unsafe for **later blinded reading** even if it was causally blind at seal time. That is exactly the already-open P2 architectural question.

Engineering action:
- add the four currently undeclared V3/V4 numeric event-log cases to P2's declared case set, with an explicit note that this is **classification under an OPEN adjudication, not exemption or acceptance**;
- do not close P2;
- do not describe the four cases as proven pre-seal leaks;
- preserve the distinction between seal-time causal blindness and persistent later blind-read safety for the eventual P2 ruling;
- do not fold P4's KAPPA results-file cases into P2.

This scopes the red correctly; it does not make the science green.

### 2. P1A-EDGE-REPRESENTATION

**RULING: CLOSED — an R5 wiring edge is a substantive dependency/presupposition for provenance-tier consistency.**

The purpose of the R5 edge is to make the presupposition and its multiplicity visible. A tier checker is therefore correct to treat the wired input as a dependency. If a node claimed `shown` while depending through R5 on an `assumed` input, that is a real tier contradiction unless another explicit rule justifies the tier.

Current disappearance of the live tier cases is **not** treated as discharge of `background_time_translation_flow`; the booked input and edges remain. P1A is CLOSED because the semantic question is now answered, not MOOT.

Engineering action:
- record P1A CLOSED with this ruling;
- remove any expected-red declaration that depended solely on P1A;
- if a current node still violates the ruled semantics, let it surface as a new red rather than suppressing it.

P1B remains separate and OPEN.

### 3. P6-STALE-NETS-IN-STANDING-DOCS

**RULING: amend the four standing documents to replace stale present-tense net/node assertions with the machine-emitted current register values, while preserving clearly historical values as historical.**

This authorizes only that synchronization correction; it is not blanket approval of unrelated prose edits.

The cycle-2 review reports that the stale-net enumerator is now empty and supplies current hashes for:
- `GRUT_ToE.md`
- `README.md`
- `GRUT_II_Agenda.md`
- `GRUT_II_What_Survived.md`

Engineering action:
- verify that the already-landed changes satisfying P6 are limited to the authorized synchronization/historical-label correction;
- if verified, record P6 CLOSED with the ruling and landed amendment;
- if the diff contains unrelated substantive edits, keep P6 OPEN and report the exact excess diff.

### 4. Seven-repair overlay

The five-file repair package is **approved for continued review/integration preparation**, not yet for a canonical scientific merge.

The evidence presently supports:
- exactly seven FAIL→PASS changes;
- no other 242-node outcome change;
- 239 PASS / 2 FAIL / 1 SKIP in the overlay;
- all 55 registered mutants executed successfully;
- all 15 cited falsifiers met their execution contracts;
- no guard/allowlist weakening.

Do not apply an accepted-baseline or scientific-status change merely because this engineering overlay works.

### 5. CR-5 signed-weight finding

**CR-5 is not cleared while signed reference weights are accepted.**

A quantity interpreted as an appearance/reference fraction cannot accept a negative measure weight and return values such as 2 without an explicit different mathematical definition.

Required contract review:
- finite weights;
- nonnegative weights;
- strictly positive total weight when a normalized fraction is requested;
- deterministic handling/rejection of invalid inputs before appearance arithmetic.

The 183 valid-input arithmetic agreements remain useful evidence and are not invalidated by this hostile invalid-input finding. Do not silently patch private Stage-3 code from TestingGRUT; surface the required caller/kit validation for independent review.

### 6. Bank flags and scientific status

All 25 bank flags remain REVIEW_REQUIRED. Do not auto-accept them as a consequence of these owner dispositions.

Engineering green and scientific green remain separate.

### 7. Next deliverable

After implementing the review-only reconciliation above, return one **OWNER READINESS SYNTHESIS** containing:

1. repaired default/overlay outcome counts;
2. classifier status after P2 scoping, P1A closure, and P6 disposition;
3. exact remaining OPEN passes and why they are scientifically live;
4. the 25 bank flags grouped by disposition needed;
5. evaluator independent-review status;
6. CR-5 status including signed-weight contract;
7. whether any blocker remains that is genuinely required **before Card 1 can be opened**, rather than before later lock/canonical publication;
8. a final recommendation: `READY FOR OWNER CARD-1 DECISION` or `NOT READY — <minimal blockers>`.

No new foundations lane. No candidate K or IP-12 draft.
