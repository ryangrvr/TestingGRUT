# Reproduce the exact architecture checks

From the repository root, Python and NumPy are sufficient:

```bash
python -B research/reality_framework/grut_exact_architecture_20261026/code/history.py
python -B research/reality_framework/grut_exact_architecture_20261026/code/operators.py
python -B research/reality_framework/grut_exact_architecture_20261026/code/package.py
```

The first script constructs the 64 by 64 constraint, checks its spectrum and
physical POVMs, computes five example conditional histories, checks the
maximum-entropy state, and tests an additional state filter. The p-family
and generic mixed-state result are derived analytically in RECORDS.md.

The second script enumerates every exponent tuple through dimension four,
uses exact rational polynomial derivatives for the scalar Hessian identity,
runs 40 independent gradient-difference comparisons, evaluates three
vector-potential derivative resolutions and four analytic e sensitivities,
and checks the curvature-mass representative on three ideal backgrounds.
The large gauge-potential derivatives are explicitly not physical pole masses.

Eight groups pass. package.py aggregates their actual results, produces
THEORY_DATA.json, pins inspected internal sources and checks JSON, Python
syntax, relative links and the file manifest. It does not execute the pinned
earlier research or papers. No root artifact, kinetic solver, original data
or reviewer branch is modified by these scripts.

An initial attempt to use SymPy found it unavailable in the current runtime.
The calculation instead uses a small exact Fraction polynomial algebra plus
independent numerical derivatives; no new dependency was installed. Some
documentation patches with stale context stopped before modifying files and
were applied with corrected context. No scientific input or tolerance was
changed to manufacture passing checks.

The finite model is an engineered standard-quantum history construction.
It is not the scalar field theory, its ground energy is not the universe's
measured total energy, its finite clock does not imply cyclic cosmological
time, and its record interval proves no macroscopic lifetime. The field
action remains an EFT with supplied matching and state data. Full quantum
gravity, a cosmological state law, a physical phase holder, observed vacuum
density, the Born rule and a distinctive empirical signal are not derived.
