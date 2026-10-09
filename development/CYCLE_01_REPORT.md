# Issue #2 — engineering cycle 1

Engineering integrity is **BLOCKED**; scientific status is **REVIEW_REQUIRED**.
The workbench reproduces the frozen observation but cannot satisfy the queue's
condition that every failing test/case be validly declared against an open pass.

The pristine reference at `114af03f5a9851ccc7ae5d681e2457203d192e80`
reproduced 242 outcomes: 232 PASS, 9 FAIL, 1 SKIP. The final local instrument
reproduced all of them with **zero baseline deltas** and zero infrastructure
errors. Both validators, selected imports, syntax and all **22** reporting
tests pass. The unchanged legacy classifier was actually executed; its failure
inventory, verdict and diagnostics agree with the independent audit.

The exact blockers are seven failing tests without declarations, four live
pointer-leak cases outside the declared case set, and uncited nonsymptomless
OPEN passes **P1A-EDGE-REPRESENTATION** and **P6-STALE-NETS-IN-STANDING-DOCS**.
They are individually listed in `CYCLE_01_LOCAL_RESULTS.json`. Both declared
tests remain red. No declaration has been expanded to cover the other seven.
Bank inventory remains **25 unreviewed flags**: 24 flagged changes and one deletion.

Implemented a development-only observational manifest/hash lock, case-level
comparison, source and mid-run code stability, actual classifier execution,
independent adjudication auditing, fresh report directories, strict JUnit/count
validation, artifact completeness and source/dirty metadata. GitHub gets two
separate status checks. Missing reports and conflicting status lines cannot
certify green. An actual small pytest fixture demonstrates detection of an
injected new red and a removed declared red.

The manifest's protected-input digest covers **1,757** scientific/source inputs.
Its lock is `b46f65a0ad10193f9443dffee4227b041949dda8541d7f7e98a109ead70a396d`.
The default case digest is
`fc6c0495455622e495d16dcb48bbbdd1ba2fb3622ffd5ca102e8ab96013cfc42`.

Separate full-mutation and slow-falsifier jobs and profile comparisons are
configured, with manual PR-title activation and retained artifacts. Those
expensive jobs were **not executed in this cycle**. Their successful execution
is not inferred from their configuration or from default tests.

Evidence grade: **builder-executed, same-author cross-check**. This is not
external approval. Local measurements record source `114af03` with a dirty
development tree honestly; no post-commit local run is invented. Hosted evidence
and the published commit are recorded in draft PR #1 after upload.

No accepted register baseline, owner ruling, charter, preregistration, bank
decision, candidate law, scoring, pricing, or Card-1 gate was changed. CR-5's
private review lane remains separate. Card 1 is closed; cards consumed: **0**.

Next narrow task: owner/reviewer reconciliation of the classifier declarations,
sealed-pointer scope and P1A/P6 state, with explicit review of any subsequent
observational-manifest refresh. This queue cannot adjudicate those questions.
