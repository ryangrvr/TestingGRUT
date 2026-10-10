# Current review lane — Issue #9

The current branch is `codex/final-precard-clearance`, starting exactly at
`b89a0e607e71450fce07d5cc55daae8224acad18`. See
[FINAL PRE-CARD REVIEW PACKET](final_precard/REVIEW_PACKET.md) and its pinned
`REVIEW_MANIFEST.json`. Status: **READY FOR EXTERNAL CHECK**; same-author evidence,
not external clearance. P6 is CLOSED by the exact owner-directed clause removal.
Card 1 remains closed and all 25 bank flags remain unchanged.

```bash
bash development/final_precard/review.sh /absolute/path/GRUT_FINAL_PRECARD_REVIEW_INPUTS.zip /fresh/output/directory
python development/checks.py --reconciliation development/final_precard/P6_TRANSITION.json
python development/reconciliation_review.py --output development/results/fresh-overlay --stamp 2026-10-10 --reconciliation development/final_precard/P6_TRANSITION.json
```

The first command requires Python 3.12.14 and the separately supplied private
archive; no private Stage-3 source is committed. The record commands additionally
use the previously pinned test dependencies. The actual record is 234/7/1; the
unintegrated review-only repair is 239/2/1 with a clean classifier. Five repair
failures remain in actual source. The two axes continue to report engineering
BLOCKED / science REVIEW_REQUIRED. No allowlist expansion follows from readiness.

The material below is the **historical b89a0e6 workbench checkpoint**, retained for
continuity. Its P6 OPEN status and original transition describe that checkpoint.
Current invocations use `final_precard/P6_TRANSITION.json`, never refresh the
original observational manifest, and admit only the exact P6 source/outcome effects.

# GRUT engineering workbench

The owner requested completion after the Issue #9 packet. The exact reviewed
five-file repair is now applied on `codex/engineering-finish`. Run:

```bash
python development/finish_checks.py
```

This checks the actual applied tree against the fixed c32 review snapshot and
the already frozen repair hash. Only five named FAIL→PASS effects are permitted;
the original observation manifest, declarations, OPEN passes, seals and bank
inventory remain unchanged. The full current provenance suite and live
classifier are run. Historical owner/packet tests are executed at their exact
input snapshot; additional hostile tests cover the applied transition.

The full check reproduced a live-registry mutation by the hostile pin fixture.
Its additional, separately frozen repair runs the unchanged checker on copied
inputs. The copied baseline must pass and the changed tier must fail with the
target named. The ordinary pin guard still checks the actual register. No test
identity or expected outcome changes under this isolation repair.

The expected result is **engineering PASS**, **scientific REVIEW_REQUIRED**.
The immutable evaluator/CR-5 packet remains at
`c32f672288fcad2a9d08772b8ef2bf46ddee70a2`; replay it from that exact checkout.
Its status is still **READY FOR EXTERNAL CHECK**, not external approval.
Card 1 is CLOSED. The older sections below are historical engineering records.

This branch implements [issue #2](https://github.com/ryangrvr/TestingGRUT/issues/2)
and remains the development-only [draft PR #1](https://github.com/ryangrvr/TestingGRUT/pull/1).
It supplies no scientific approval. Card 1 remains closed.

```bash
python3 -m venv development/venv
development/venv/bin/python -m pip install -r development/requirements-test.txt
development/venv/bin/python development/checks.py --reconciliation development/owner_reconciliation/transition.json
```

The frozen observation was independently reproduced in a pristine worktree at
`114af03f5a9851ccc7ae5d681e2457203d192e80`: 242 tests, 232 PASS, 9 FAIL,
1 SKIP. The manifest records every exact node/outcome, declared and live case
set, OPEN_PASS state, full structured bank inventory and held flags, protected
source digest, and measured tool versions. The prior hosted run and its artifact
digest are retained as provenance. The new reproducible environment pins Python
3.12.14 and all five test packages; the prior host was reported as Python 3.12.15.

`expected_red_manifest.json` is an **observational snapshot**, not an expansion
of `provenance/expected_red.py`'s adjudication authority. Its adjacent SHA256 lock
must match. `freeze_manifest.py` refuses to overwrite either file. The runner
never regenerates them. Future changes require explicit review of both files;
the issue's source pin is enforced separately in code.

The default run validates source syntax and selected safe imports, runs both
validators and the bank gate, runs the reporting regression tests, executes raw
pytest, and **runs the unmodified `provenance/expected_red.py`**. It compares the
classifier's independently executed failure inventory and diagnostics with raw
pytest and an independent case/adjudication audit. It verifies that inputs and
source identity stayed unchanged throughout the run.

The two report axes and GitHub checks are distinct:

| Axis | Meaning |
|---|---|
| Engineering integrity PASS | Infrastructure works, every baseline observation matches, and all observed failing tests/cases are validly declared against open adjudications |
| Engineering integrity BLOCKED | Baseline and instruments agree, but the required adjudication condition is unsatisfied |
| Engineering integrity FAIL | Infrastructure error, classifier disagreement, or an unapproved observation/input/version delta |
| Scientific status REVIEW_REQUIRED | Bank flags, open passes or provenance failures remain; engineering success does not resolve them |

After the Issue #2 owner direction, P2 owns the four current event-log numeric
cases (classification while OPEN, not acceptance), and P1A is CLOSED by semantic
ruling. P6 remains OPEN: original document diffs contain an extra consolidation/
retirement clause, so its conditional closure requirement is not certified.
The actual default retains seven unrepaired machinery failures; the review-only
overlay repairs those seven and leaves P6 as its only classifier obligation.
The 25 bank flags (24 flagged changes plus one deletion) remain unaccepted.
See [OWNER_READINESS_SYNTHESIS.md](OWNER_READINESS_SYNTHESIS.md).

The original manifest/lock remain historical and unchanged. The explicit locked
`owner_reconciliation/transition.json` records only the authorized P2/P1A changes
and the exact two protected source paths; it cannot refresh test outcomes, bank
flags or any other adjudication. Omitting `--reconciliation` still compares with
the original snapshot and correctly rejects the now-authorized source delta.

Every invocation uses a new results directory. Requested JUnit files are removed
before a subprocess runs; missing/empty/malformed reports, duplicate node IDs,
conflicting outcomes, dishonest counts, unexpected pytest exits and missing
assets fail. JSON/Markdown summaries record actual checkout SHA, dirty state,
all deltas, every adjudication blocker, and hashes of the run's report assets.
Outputs live under `development/results/` and remain ignored by Git.

Separate expensive guard profiles are available:

```bash
development/venv/bin/python development/checks.py --profile full-mutation --reconciliation development/owner_reconciliation/transition.json
development/venv/bin/python development/checks.py --profile slow-falsifiers --reconciliation development/owner_reconciliation/transition.json
```

The first enables `GRUT_FULL_MUTATION=1` and runs `test_mutation_battery.py`.
The second enables `GRUT_RUN_SLOW=1` for the complete suite. Their outcomes are
compared with explicitly defined profile expectations: the mutation subset must
retain its passing baseline, and the formerly skipped falsifier execution guard
must pass when enabled. These are pre-existing guard contracts, not observed
claims that the expensive runs succeeded. A profile run reports its execution
integrity and leaves the default engineering axis NOT_ESTABLISHED; it cannot
replace the default classifier run. Baseline input/adjudication/bank drift is
still rejected. Cycle 1 did not execute the expensive jobs. Cycle 2 executes both profiles successfully with pinned numerical dependencies; see CYCLE_02_REPORT.md.

To activate the separate hosted jobs on this development branch, manually edit
draft PR #1's title to include `[full-mutation]` or `[slow-falsifiers]`. Only a
title-edit event starts the selected expensive job; normal pushes keep the fast
loop. Restore the ordinary title afterward. Both jobs retain logs even on failure.
`workflow_dispatch` inputs are also configured. GitHub requires a workflow on the
default branch for that dispatch route; the PR-title route works with this draft
branch and leaves the frozen default untouched. See
[GitHub's manual-run documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow).

Actions use pinned release commits and read-only tokens. Scientific status has
its own failing check whenever review is outstanding. The PR stays draft.
CR-5's private inputs and independent review lane are outside this cycle.

## Cycle 2 review-only repair

See [CYCLE_02_REPORT.md](CYCLE_02_REPORT.md) for the concrete five-file patch, exact 242-case overlay comparison, owner dispositions, separate CR-5 review and expensive-profile execution evidence. The default protected record and locked manifest remain unchanged.

```bash
python -m pip install -r development/requirements-numerical.txt
python development/reconciliation_review.py --output development/results/review-overlay --stamp 2026-10-10 --reconciliation development/owner_reconciliation/transition.json
python -m development.cr5_review.check --bundle /absolute/path/GRUT_G2_14R_CR5_REVIEW_BUNDLE.zip --output development/results/cr5-review
```

The CR-5 bundle is kept outside this repository. Its exact hash is pinned in the harness. A separate arithmetic implementation is not external authorship; synthetic receipts do not certify physical pipelines. Invalid signed reference weights remain an explicit review finding.
