# GRUT engineering workbench

This branch implements [issue #2](https://github.com/ryangrvr/TestingGRUT/issues/2)
and remains the development-only [draft PR #1](https://github.com/ryangrvr/TestingGRUT/pull/1).
It supplies no scientific approval. Card 1 remains closed.

```bash
python3 -m venv development/venv
development/venv/bin/python -m pip install -r development/requirements-test.txt
development/venv/bin/python development/checks.py
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

The present blocker is explicit: nine observed failing tests versus two declared
tests, four additional live pointer-leak cases, and uncited nonsymptomless OPEN
passes P1A and P6. Keeping these unchanged does not authorize engineering green.
The 25 bank flags (24 flagged changes plus one deletion) retain their original
review state. No accepted baseline, seal, ruling or adjudication is rewritten.

Every invocation uses a new results directory. Requested JUnit files are removed
before a subprocess runs; missing/empty/malformed reports, duplicate node IDs,
conflicting outcomes, dishonest counts, unexpected pytest exits and missing
assets fail. JSON/Markdown summaries record actual checkout SHA, dirty state,
all deltas, every adjudication blocker, and hashes of the run's report assets.
Outputs live under `development/results/` and remain ignored by Git.

Separate expensive guard profiles are available:

```bash
development/venv/bin/python development/checks.py --profile full-mutation
development/venv/bin/python development/checks.py --profile slow-falsifiers
```

The first enables `GRUT_FULL_MUTATION=1` and runs `test_mutation_battery.py`.
The second enables `GRUT_RUN_SLOW=1` for the complete suite. Their outcomes are
compared with explicitly defined profile expectations: the mutation subset must
retain its passing baseline, and the formerly skipped falsifier execution guard
must pass when enabled. These are pre-existing guard contracts, not observed
claims that the expensive runs succeeded. A profile run reports its execution
integrity and leaves the default engineering axis NOT_ESTABLISHED; it cannot
replace the default classifier run. Baseline input/adjudication/bank drift is
still rejected. This cycle does not claim expensive-job execution.

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
