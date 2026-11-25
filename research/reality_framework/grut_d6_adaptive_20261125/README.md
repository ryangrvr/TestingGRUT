# Reproduce GRUT-D6

This packet is rooted at draft [PR #18](https://github.com/ryangrvr/TestingGRUT/pull/18).
The [frozen model](MODEL.md) specifies each link and preparation.
[PROTOCOL.json](PROTOCOL.json) supplies the full grid, readouts,
and acceptance thresholds. [FREEZE.json](FREEZE.json) records the
pre-production commit and content hashes.

| File | Content |
|---|---|
| [REPORT.md](REPORT.md) | Results and separate acceptance verdicts |
| [PROOF.md](PROOF.md) | Exact solution, optimal information bound and finite late window |
| [CONTROLS.md](CONTROLS.md) | Intervention definitions, comparison failures and timing diagnostics |
| [RESOURCES.md](RESOURCES.md) | Preparation, dimensions, clock, energy, work and entropy |
| [CRITIC.md](CRITIC.md) | Independent internal review and replication scope |
| [RESULTS.json](RESULTS.json) | Full set of 864 conditional rows and numeric checks |
| [ORACLE_RESULTS.json](ORACLE_RESULTS.json) | Independent physical propagation for the registered small-system rows |

## Reduced exact model

Run from this directory with Python 3.11 or newer. The main
calculation uses only the standard library.

```bash
python3 run_d6.py --verify RESULTS.json
python3 run_d6.py --output /tmp/grut-d6-results.json
```

`columns` in RESULTS.json gives the fields in each conditional
row. Coins are averaged in every row; source, challenge and arm
remain separate. `null` in the endpoint's log-tail field means
an exactly zero unfinished probability. Finite-time scores
include arrival and completion as defined in the model.

## Independent physical propagation

The reviewer builds the reachable matrix directly from the
local basis transitions and the witness term. The implementation
uses NumPy and SciPy. The recorded run used Python 3.12.14,
NumPy 2.3.5, and SciPy 1.17.0.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 independent_oracle.py \
  --output /tmp/grut-d6-oracle.json
python3 cross_check.py --oracle /tmp/grut-d6-oracle.json
python3 cross_check.py
```

The comparison uses the frozen grid and checks both solvers
against the registered tolerance. The direct physical matrices
are tested at n=2. The larger-size results use the exact
one-head theorem. External RAI review remains open.
