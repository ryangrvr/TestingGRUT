# Reproduce GRUT-D7

This additive packet extends [draft PR #19](https://github.com/ryangrvr/TestingGRUT/pull/19).
The local model and exact comparison are in [MODEL.md](MODEL.md)
and [PROOF.md](PROOF.md). [REPORT.md](REPORT.md) gives the
outcomes. [CRITIC.md](CRITIC.md) specifies the recovery access
boundary and the limits of the result. [FREEZE.json](FREEZE.json)
records the commit and hashes for the original model and
[protocol](PROTOCOL.json).

The recorded environment used Python 3.12.14, NumPy 2.3.5 and
SciPy 1.17.0. In this directory run

~~~bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 run_d7.py --verify RESULTS.json
python3 physical_oracle.py --verify ORACLE_RESULTS.json
~~~

Each script also accepts --output and an output path to write
a fresh replay. The physical implementation checks each
registered row against RESULTS.json. It constructs the reachable
Hamiltonian from the local forward and reverse exchanges and
uses a sparse matrix exponential. It also compares the initial
energy law and verifies full recurrence. L=32 is checked by
the analytic reduction and bound.

The rows in RESULTS.json use the ordered `columns` list.
Every row averages the two independent mixed coins and retains
its source condition explicitly. The two nonterminal window samples
are its registered flanks. No row resorts
to success-conditioned normalization.
