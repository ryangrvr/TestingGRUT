# Reproduce

From this directory with Python, numpy, scipy and matplotlib installed:

```sh
OPENBLAS_NUM_THREADS=1 python code/calculate.py
python code/plot.py
python code/verify.py
python code/verify.py --replay
```

The calculator reads only INPUTS.json; it does not import old packet solvers
or fetch observations. It writes 34 endpoints, five sampled histories,
RESULTS.json and ENVIRONMENT.json, and runs 41 internal numerical conditions.
The plot uses these recorded outputs. No new dependency installation was needed.

For a new draft, run `python code/make_manifest.py` only after finishing
all documentation and figure changes. The published manifest hashes all
packet files except itself and Python cache directories. The read-only
verifier does not update outputs or regenerate hashes.

The default verifier checks the full manifest and file set, exact endpoint
CSV/JSON agreement, raw Friedmann budgets, separate energy and scalar-momentum
quadratures and analytic attractor eigenvalues. `--replay` copies inputs and
source to a temporary directory, executes a fresh calculation, requires
identical scientific key sets, compares every scientific JSON leaf and all
six raw CSV arrays, and deletes the temporary output. Solver evaluation
counts are excluded from physical comparisons; all scientific result and
check fields remain included.

Same-code replay, a second numerical clock and separately evaluated
identities are internal verification by one author. They do not constitute
outside replication or a test of the underlying action in nature. To
compare with the archived run without rewriting it, use `--replay` or a
separate checkout.
