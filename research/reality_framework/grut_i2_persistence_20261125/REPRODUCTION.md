# Reproduce GRUT-I2

The archive is additive to [draft PR #21](https://github.com/ryangrvr/TestingGRUT/pull/21).
The inputs were frozen before production at
4dcc33c4720969b2898d7a0320a4a9ec3de9ea39.
[Freeze digests](FREEZE.json) · [Protocol](PROTOCOL.md) · [Parameters](PARAMETERS.json)

Run from a standalone copy of this directory using Python 3.12.14,
NumPy 2.3.5 and SciPy 1.17.0, the recorded environment.

~~~bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
timeout 240s python3 run_i2.py --verify RESULTS.json
timeout 240s python3 labelled_check.py --verify LABELLED_RESULTS.json
~~~

Both executables also accept --output and a destination path.
The labelled solver uses RESULTS.json next to its own script for
the numerical comparison. The main solver constructs the count
generator and uses sparse exponential action. The verifier builds
a labelled molecule chain with its own physical moves and sums a
Poisson uniformization with a bounded truncation tail. It also
compares the history statistic. These are separate algorithms in
the packet, without external peer replication.

The frozen outputs contain 64 main policy settings, four uniform
formation settings, a flat persistence control, nine one-patch size
rows, nine fuel-capacity controls, 21 canonical rows and six small
reference cases. Preserve the saved infeasible budgets and targets.
The post-run dominance summary is a consequence of the already
executed menu. It does not modify a frozen parameter or gate. The
full deterministic computations ran within the 240-second budget
with one BLAS thread. The recorded replay cross-check compares
89 scalars with maximum absolute difference $9.77\times10^{-15}$.
[Count results](RESULTS.json) · [Labelled check](LABELLED_RESULTS.json)

Temperature, rates and coupling supports are model inputs.
The lab assay in the novelty audit is prospective. The packet
is ready for independent external review. [Review packet](REVIEW_PACKET.md)
