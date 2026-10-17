# GRUT-RAI bounded autonomous discovery

This experimental workstream starts at P1 commit
`643ec7f9af41ec5d8d7816a5dcddbe24de2de7a7`. It does not reopen the canonical
GRUT-RAI SYN campaign or modify previous candidates. The owner requested a
new TestingGRUT branch and an executed discovery cycle.

The runner uses the actual GRUT-RAI discipline engine,
`provenance/auditor.py`, with a separate claim ledger. It is a deterministic,
bounded experiment runner, not a second AI service or a consciousness claim.
Builder, Challenger, Integrator and Recorder are functions performed in this
session, with shared model/context ancestry. There is no independent reviewer.

## Run and preserve

From the repository root, with Python, NumPy and SciPy available:

```sh
python research/reality_framework/grut_rai_discovery_20261017/runner.py verify
python research/reality_framework/grut_rai_discovery_20261017/runner.py run --id replay-001
python -m unittest discover -s research/reality_framework/grut_rai_discovery_20261017 -p 'test_*.py' -v
```

`run` requires a committed `FREEZE.json` containing the protocol, implementation,
claim ledger, auditor and P1 implementation hashes. A changed file is rejected.
Run directories are created exclusively and never overwritten. A stopped or
failed run remains visible and a rerun must use a new ID. The process does not
loop, call paid APIs, change parameters, publish, or continue in the background.
New hypotheses require a separately reviewed and frozen successor packet.

Read `PROTOCOL.json` before execution, `MATHEMATICS.md` for the exact
two-state diagnostic, and `REPORT.md` for the executed outcome. Every small
matrix here is a diagnostic, not an all-size geometry theorem.
