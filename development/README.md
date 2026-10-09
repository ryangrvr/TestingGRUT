# Development checks

This branch adds a repeatable build/check loop to `ryangrvr/TestingGRUT`, starting
from `scout-0` at `ab2da47407bb670d94e2a52c87599fa13fd8ab99`.
It changes development infrastructure, not the deposited scientific record.

```bash
python3 -m venv development/venv
development/venv/bin/python -m pip install -r development/requirements-test.txt
development/venv/bin/python development/checks.py
```

The command parses the active provenance Python files, runs both provenance
validators, checks the bank gate, tests its own reporting, and runs the complete
default provenance test suite. Logs, JUnit results, and `summary.json` go in
`development/results/`, which is ignored by Git. The source commit and whether
the working tree was dirty are recorded. Python 3.12 is the CI target.

The bank gate deliberately returns exit zero for `FLAG-FOR-FIREWALL`.
This runner reports that outcome as `REVIEW_REQUIRED` and exits nonzero overall.
A raw test failure also remains a failure, including declared adjudications.
Nothing accepts a register baseline, seals a preregistration, issues a ruling,
or supplies external scientific approval. For the separate adjudication process,
see the existing `HOW_TO_VERIFY.md` and `provenance/expected_red.py`.

The default profile excludes slow mutants and long falsifier runs, even when
their environment variables are set in the caller. Run the explicit commands
in `HOW_TO_VERIFY.md` when those checks are required.

The prepared GitHub Actions workflow runs this same command on branch pushes,
pull requests, and manual dispatch, retaining reports even when checks fail.
Actions are pinned to release commit hashes. Its token has read-only repository
access. A GitHub-hosted execution is not established until the branch is actually
pushed and a workflow run is observed.

TestingGRUT is public. This branch contains no private Site content or current
private Stage-3 research files. The latest owner-approved CR-5 review target is
maintained separately; this repository's older default snapshot is not substituted
for the current Stage-3 canonical record. Card 1 remains closed.

At the base commit, the full default suite collected 242 tests: 232 passed,
9 failed, and 1 skipped. Both provenance validators passed. The bank gate
reported 25 new flags and `FLAG-FOR-FIREWALL`. The new runner reproduced all
242 baseline case outcomes exactly, and its four reporting tests passed.
These are local builder checks, not an independent review or a clean baseline.
