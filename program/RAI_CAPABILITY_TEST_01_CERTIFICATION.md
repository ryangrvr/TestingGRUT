# RAI_CAPABILITY_TEST_01_CERTIFICATION — CONDITIONAL certification of |N| = 0 (Phase-0 gate re-run)

**Certification statement (the whole result, one sentence):**

> **|N| = 0 — CERTIFIED, CONDITIONAL on the rung1/rung2 edge representation, which resolved
> the foundation tier contradiction by demotion rather than by adjudication
> (P1A-EDGE-REPRESENTATION: MOOT, unruled).**

Not a clean certification, by design. The conditional is on the face of the certificate
because the green this certification rests on embeds one unadjudicated foundational
assumption, named below. A clean "|N| = 0 CERTIFIED" would launder that assumption into a
settled state.

Authorized by owner directive 2026-09-12: *"yes, authorize it — but authorize it with the
P1A residual named on the cert, not clean."*

---

## 1. What is being certified, and what is not

**Certified:** the R/G/N classification frozen in `RAI_CAPABILITY_TEST_01.{md,json}` —
eight claims binned, zero survivors in (N) — now stands on a gate-passing instrument
rather than a red one. The PROVISIONAL status those frozen artifacts carry is superseded
by this file; the artifacts themselves are untouched (§4).

**Not certified:** correspondence to truth (forbidden by the test's own success criterion);
any closure of the open passes P2/P3/P4/P5/P6; any ruling on P1A; any authorization of the
Born Phase-2A hostile audit or Phase 2B. |N| = 0 remains a statement about what the
instrument can generate from the corpus it holds — it bounds generative reach, not
physics' content (frozen caveat 1 — the PROVISIONAL/re-run clause — is discharged by the
§2 re-run this file records; caveats 2 and 3 carry forward unchanged).

## 2. Phase-0 gate re-run (the precondition the frozen artifact failed)

Run 2026-09-12 at HEAD `161a883` ("Provenance suite restored: failing set == declared
set"), working tree clean of tracked modifications (sole untracked file: this
certification draft itself):

| Instrument | Result |
|---|---|
| `pytest provenance/` | **3 failed / 236 passed / 1 skipped** |
| `expected_red.py` | **exit 0 — failing set == declared set** |
| Register invariants | 74 claims total, 53 GRUT-scope, net **+16** |
| rung1_inin_formalism / rung2_kms_gate | `derived-pending` (both) |

The three failures are, by name, the three DECLARED adjudications, each citing an OPEN
pass: `test_no_standing_doc_asserts_a_stale_net` (9 enumerated prose cases,
P6-STALE-NETS-IN-STANDING-DOCS) and the two `TestBlindSafe` prereg tests
(P2-TERMINATION-EVENTLOG / P3-SEALED-HISTORY-POLICY / P4-TERMINATION-KAPPA-RESULT;
sealed preregs are immutable, so these reds are permanent until ruled).

**What "green" means here — and does not mean:** green is *failing set == declared set*,
classification not suppression. It is NOT zero failures. A reader who wants zero failures
is asking for the laundering this suite exists to prevent. The frozen artifact's Phase-0
failure was a red state the two records count differently — 14 failures (artifact record)
/ 15 (restoration commit); the discrepancy is noted here, not resolved. What was
undiagnosed was the SET (failing set != declared set), not every member: three were
long-declared adjudications. The restoration (commit `161a883`) repaired the rest at
their causes and left exactly those three declared, open-pass-cited reds.

## 3. THE CONDITIONAL — what this green rests on

**The mechanics, exactly.** The register's foundation contradiction (nodes at tier
`shown` resting on the `assumed` input `background_time_translation_flow` through their
R5 edges) was dissolved by the honest tier propagation: `rung1_inin_formalism` and
`rung2_kms_gate` moved `shown` → `derived-pending`. That is a resolution **by demotion**.
The same propagation (one commit, `67059a8`) cascaded six further nodes through the
rung1/rung2 dependency cone, `shown` → `derived-pending`: `rung4_love_kk`,
`kr_contract_retarded_tier4` (banked at `shown` on owner relay eleven days earlier,
`d5e9a99`), `info_i1_renorm_as_information`, `passivity_channel_diagonal`,
`rung7_w1_wz_map`, `u1_form_universality` — eight demotions in all; the rung2/rung4
doc-pin re-acceptances recorded in `161a883` travel with them. The owner's condition
names rung1/rung2 because they carry the contested edge; an adverse edge ruling puts all
eight demotions, and the pin re-acceptances resting on them, up for review. (The in-tree
`expected_red.py` comment names only four of the eight;
`git diff d5e9a99..67059a8 -- provenance/claims.json` is the authority.)

The open pass P1A-EDGE-REPRESENTATION existed to adjudicate the prior question — *does an
R5 wiring edge, which exists to make a presupposition's multiplicity visible, mean what
the tier-contradiction check reads it as meaning?* — and it was never ruled. It went
**MOOT** on 2026-09-12: the symptom dissolved, by a cause its stanza did not anticipate,
without the question being answered. The MOOT status was set by the builder in `161a883`
— the stanza's only pre-authorized moot route (the static-patch discharge) did not occur
— and is encompassed by the owner's 2026-09-12 authorization, which post-dates it.

**Why that conditions this certificate.** The demotion *presupposes* the
tier-contradiction check's reading of the edge — it demoted the tiers rather than
re-typing the edge. So the register this test read (and this gate certified) is the
register *as booked under that reading*. The frozen artifact's Phase-1 foundation line
("rung1/rung2 = `derived-pending` (honest cascade from `background_time_translation_flow`
= assumed)") states the demoted tiers as fact; they are fact only conditional on the edge
reading.

**What happens if the conditional fails.** If the owner ever rules the edge
representation wrong, the demotions are up for review (recorded in the P1A stanza), the
register as booked is revised, and **this certification is VOID pending re-run** — it
does not survive the re-ruling by default.

**Where the condition attaches.** The gate did not require the propagation: before
`67059a8` the resident tier-contradiction was a DECLARED red citing OPEN P1A, so a lawful
green (failing set == declared set) was reachable with the tiers still `shown`. The
condition therefore attaches to the certified content — the frozen artifact's Phase-1
line and the register as booked — not to the gate-pass; the counterfactual route would
have carried the same P1A residual as an open declared red.

**Scope note (bounded, not exculpatory).** The R/G/N reductions in the frozen artifact
are content-based: no row's reduction argument invokes the tier labels of any of the
eight demoted nodes as a premise (the only mentions anywhere in the artifact are the
Phase-1 foundation line and row 7's title parenthetical) (row 7 rests on the
REALITY_CHECK_03 / Experiment P effective-layer results and the ontological-surplus
argument; rows 1–6 and 8 reduce to recovered physics or gate consequences independent of
the foundation tiers). So an adverse edge ruling would not, by itself, be expected to
manufacture an (N). But that expectation is an *argument*, not a certification — a re-run
after any edge re-ruling must re-establish it against the revised register, not inherit
it from this paragraph.

**Residuals under the same green, disclosed for completeness (audit-surfaced):**
- P1B-SHOWN-ON-LEDGER-INPUTS remains OPEN with `MAY_HAVE_NO_SYMPTOM: yes` — no instrument
  can see it; the runner prints it on every run under its standing-questions banner
  ("A green suite is not an answer to these"); the stanza's own words: "A GREEN SUITE
  WOULD THEN BE THE ANSWER, AND NOBODY WOULD HAVE HAD TO GIVE ONE." Coupled to the
  conditional above: the demotion that mooted P1A also removed P1B's live subject — no
  `shown` node now carries prose-priced inputs (the four persist on
  `rung1_inin_formalism` at `derived-pending`, ledger_delta 4) — so an adverse edge
  ruling re-promoting the node would make P1B acutely live again.
- P5-RESIDENT-SCOPE remains OPEN with `MAY_HAVE_NO_SYMPTOM: yes`, printed by this same
  run beside P1B: `resident.check_change()` is not scope-aware and judges the
  vacuum-cluster nodes against GRUT's tier vocabulary; which vocabulary the resident
  enforces per scope awaits the owner. Its failure mode is SAFE — it blocks, it does not
  wave through — but it is live under this green and invisible to it.
- The two `test_resident` ORPHANED-RESULT guards are exercised on synthetic in-memory
  tier promotions (dated fixture comments in `test_resident.py`), because the propagation
  left no live result-tier node with non-empty `depends_on`. The guard code itself is
  unchanged and still runs on the live register; what the suite no longer demonstrates is
  the flag firing on a *live* node — there is none for it to fire on.

## 4. Provenance and freeze discipline

- Frozen artifacts untouched, verified byte-identical at certification time:
  `program/RAI_CAPABILITY_TEST_01.md` blob `1555a863cbb7e81b8b88e88d89d74d04f83b5f70`,
  `program/RAI_CAPABILITY_TEST_01.json` blob `fa9dbd43b546bc9369c48a4c5254283de3462f28`
  (identical to their `89445a8` archive commit blobs).
- This file is the separate citing result, per the program's standing discipline (the
  certifier must not sit inside the thing being certified) and per precedent
  (`EC_BORN_01_PHASE2A_ADJUDICATION`, `HISA01_REPAIRED_43_ADJUDICATION`).
- `REALITY_STATUS_01.md` self-declares read-only/frozen and is untouched; its component
  table predates the capability test and carries no |N| row; this file ADDS the
  conditional |N| = 0 certification to the program record and supersedes nothing in that
  table. Its "Genuine GRUT prediction — PREDICTED = EMPTY" row is REALITY_CHECK_01's
  distinct result (prediction extraction — a different instrument and measurand from
  generative reach) and stands untouched.
- The 9 stale-net prose cases remain declared-red under the open P6 pass, deliberately
  unfixed pending the owner's ruling (owner directive 2026-09-12: leave them exactly
  there — this certificate is the first repository record of that directive).
  Disclosure: the declared set was TEN cases until `161a883` — the certified commit
  itself — inserted an "as of the 2026-07-02 posing" historical cue into P6-governed
  `GRUT_II_Agenda.md:7`, whereupon the enumerator stopped emitting the case and its
  declaration was removed (10 → 9). That false-positive call was the builder's alone,
  made without the P6 ruling; it is owner-ratifiable under P6, and the P6 stanza still
  reads "ten lines" against the live nine pending that ratification.

## 5. Two-audit record

Per the standing protocol (provenance verifier + hostile refuter run independently +
mandatory referee), executed 2026-09-12 against this document's draft. All three roles
ran as independent agents with repository access; the referee received both reports and
the repository.

- **Verifier: PASS_WITH_CORRECTIONS.** All gate numbers, failing-test identities and
  pass citations, register invariants, P1A/P1B stanza content, blob hashes, freeze
  claims, and directive fidelity verified against the repository, with a row-by-row
  finding that none of the eight reduction arguments invokes a demoted tier label as a
  premise. Five corrections: the pre-narrated §5; the phantom REALITY_STATUS |N|-row
  supersession; the working-tree qualification; the condensed foundation-line quotation
  presented as verbatim; the P1B banner misattribution. All incorporated.
- **Refuter: NOT_REFUTED_WITH_AMENDMENTS.** Nine findings sustained — one
  blocking-curable (the pre-narrated §5 itself); five material (phantom |N|-row
  supersession; P5-RESIDENT-SCOPE omitted from the residuals; the demotion cone
  understated at two nodes against git's eight; the undisclosed 10 → 9 P6
  declared-case reduction inside the certified commit; the caveat-1
  supersede-and-carry-forward contradiction); three minor (the MOOT-flip actor
  undisclosed; the P1B attribution and P1A–P1B coupling; the 14-vs-15 wording). Seven
  attacks pressed and NOT sustained: joint placement of the conditional; tier-label
  independence of the eight reductions; asterisk laundering of "CERTIFIED, CONDITIONAL";
  scope-note defusal of the conditional; the synthetic ORPHANED-RESULT fixtures; the
  dynamized expected_red self-tests; the legitimacy of the frozen-artifact supersession.
- **Referee: ISSUE_AMENDED.** Every verifier/refuter contradiction ruled with repository
  evidence; two verifier verifications overruled (the "mechanics matches" finding — the
  cone was understated; the "unfixed/honest range" finding on P6 and 14-vs-15); the
  demotion cone corrected from both roles' counts (two and four) to eight from
  `git diff d5e9a99..67059a8`; the scope note retained ONLY as fenced (the
  no-inheritance clause is load-bearing; weakening it in any future edit reopens the
  laundering attack); the conditional's placement affirmed as certified-content, with
  the placement paragraph ordered onto the cert. Final: issue in the amended form —
  this document.

The draft of this certificate carried a pre-written projection of this section; it was
struck as a sustained BLOCKING defect — a certificate must not pre-narrate its own audit
— and replaced by this record.

---

*Certification issued 2026-09-12. Owner-authorized, conditional as stated. Q3 remains
OPEN; ledger delta 0; no dispatch sent.*
