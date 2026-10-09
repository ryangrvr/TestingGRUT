# Meta-law lab — Issue #5

Foundations-only audit on `codex/meta-law-lab-05`, based on PR #4 head
`c670710bd217ba01bdf8f5a58d84f072dafb5984`. Draft review surface, not canonical.

- [Decision report](REPORT.md): REDUCES TO FIXED-LAW FORMALISM, with scope.
- [Proofs and counterexamples](PROOFS.md): intervention-preserving lift,
  finite-state caveat, rank exclusion, selection and urn/do contrast.
- [Architecture and premise audit](AUDIT.md): Issue #5 sections A–F.
- [Operational discriminators](EXPERIMENTS.md): concrete comparator boundaries.
- [Generated exact tables](TABLES.md), [machine results](RESULTS.json),
  [source ledger](SOURCES.json), [local verification](VALIDATION.json).

Reproduce without dependencies on Python 3.12.14:

```sh
python development/meta_law_lab/run.py --check
```

`run.py` without `--check` explicitly regenerates machine results and tables.
CI uses check mode and fails on missing/different/stale evidence or code.
Inherited Issue #2 engineering/scientific checks run separately and keep
their declared red outcomes. Do not use this control job to bank a result
or clear scientific flags. Card 1 CLOSED; zero cards or IP-12 drafts.
