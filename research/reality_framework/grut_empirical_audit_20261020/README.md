# Candidate 1 empirical decision audit

Continue from `479f111117ea7162bd7888ff8f3792b2f6c45fdb`. This is an audit of the
unchanged Candidate-1 action, not a new empirical model or a retuned benchmark.
The owner explicitly permits a negative or underdetermined scientific outcome.

Read `REPORT.md` when the audit is complete. The pre-computation documents are
`DERIVATION.md`, `ASSUMPTIONS.md`, `PARAMETER_PROVENANCE.md`, `FEASIBILITY.md`
and `PROTOCOL.json`. `FREEZE.json` records their identities and local time;
`PUBLIC_FREEZE.json` will identify the public source commit. The maser table was
already inspected during source auditing. This is a retrospective constraint
check, not a blinded or prospective validation.

The baseline packet is `../grut_3d_empirical_20261020/`; its files are read only.
No unavailable Claude results are used. No new nonlinear production campaign,
population simulation, action extension, PR, merge or external message is part
of this audit. Physical finite-time durability and prior failed symmetry
controls retain their original scope. Same-author checks are not independent
scientific replication.

From this directory after the public freeze:

```bash
python audit.py audit
python audit.py predict
python audit.py compare
python audit.py verify
```

`predict` reads only the frozen protocol and its named calibration file.
`compare` then reads the published disk-distance table. `verify` checks frozen
identities and recomputes results without rewriting them. `audit` streams all
old raw archives to verify the manifest; it does not extract or alter them.
Output files in this new directory are generated reproducibly. Do not rerun
`freeze.py` over an existing freeze or run old preflight/freezing commands over
the baseline. A report and an independent-review checklist accompany results.
