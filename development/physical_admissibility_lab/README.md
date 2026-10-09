# Physical admissibility lab — Issue #7

**Foundations-only standard controls; DRAFT, external review pending.**
Card 1 CLOSED; zero cards. This branch does not modify canonical science,
the expected-red manifest or the inherited engineering adjudications.

- [REPORT.md](REPORT.md): verdict and direct answers.
- [PROOFS.md](PROOFS.md): explicit embeddings, obstructions and resource proofs.
- [AUDIT.md](AUDIT.md): comparator partial order and standard physical exclusions.
- [EXPERIMENTS.md](EXPERIMENTS.md): operational witnesses and scoped experimental evidence.
- [RESULTS.json](RESULTS.json) / [TABLES.md](TABLES.md): machine results / generated table.
- [SOURCES.json](SOURCES.json) / [INPUTS.json](INPUTS.json): inspection depth / pinned provenance.
- [VALIDATION.json](VALIDATION.json) / `evidence/`: same-author checks and raw baseline evidence.

With Python **3.12.14**, no scientific packages required:

```bash
python development/physical_admissibility_lab/run.py --check
python development/physical_admissibility_lab/verify_fixtures.py
```

`--check` cannot regenerate the expected evidence: changed code, results,
source/proof text or generated table fails. `run.py` without `--check`
explicitly regenerates this lane's results only; it never changes the
engineering manifest. Fractions use `{"rational": [numerator,denominator]}`.
Numerical transcendental illustrations are rounded to 12 decimal places
for portable serialization and tested with declared tolerance; they are
not exact-theorem evidence.

The inherited two-axis judge runs separately:

```bash
python development/checks.py
```

It is expected to exit 1 with engineering BLOCKED / science REVIEW_REQUIRED
until existing adjudications are genuinely resolved. A green new control
workflow is not permission to weaken that judge or merge canonical science.
