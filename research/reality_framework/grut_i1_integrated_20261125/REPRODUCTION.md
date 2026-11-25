# Reproduce GRUT-I1

The design and stopping rules were published at
cee2d2197b700b1932cdba3db758bcbecd3b864e, before production.
[Freeze digests](FREEZE.json) · [Protocol](REALITY_TEST_PROTOCOL.md) ·
[Parameter grid](PARAMETERS.json)

## Run the packet

Use Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 for an exact
recorded replay. In a standalone copy of this directory run

~~~bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
timeout 300s python3 run_i1.py --verify RESULTS.json
timeout 300s python3 reviewer_check.py --verify REVIEW_CHECK.json
python3 cross_check.py
~~~

The first program constructs the finite phase transfer, fresh thermal
population recurrence, bounded continuous limit, reused three-qubit
unitary, independent product scaling, two-site shared decay and
four graph controls. The second, separately authored program
exponentiates two-site unitaries and builds the retained $B,M,Q,R$
density evolution and Liouvillian. The third validates the matched
row grid and compares 1,080 scalar quantities. It checks both
off-diagonal Bell blocks in addition to scored fidelities. The two
solvers also accept --output and a destination path for fresh
records. Replay uses no network or archived cosmology dependencies.

The reduced solver writes 108 bridge, six reuse, 15 organization,
six pair collective and 20 graph rows. The review writes 36
coarsest-grid density rows with all jitter scores and independently
evaluated continuous outputs. The graph results include five
registered times and a finite spectral Cesàro average. The reduced
solver preserves all row data before reporting any aggregate gate
failure. The same source input, time and materials are used for
assisted and blind scores in each case. [Results](RESULTS.json) ·
[Review data](REVIEW_CHECK.json) · [Cross-check](CROSS_CHECK.json)

Only the computational replay is complete. External RAI has a
standalone replication request, and the hardware protocol is
prospective. [External handoff](RAI_HANDOFF.md)

## Recover the archive survey

The survey records the complete parent tree and selection scans.
Use the parent commit from ARCHIVE_INDEX.json in an ancestor
checkout. The NUL-delimited output of
`git ls-tree -r --full-tree -z PARENT` reproduces its inventory digest.
Decision files are REPORT.md, VERDICT.md, SYNTHESIS.md,
CLAIM_LEDGER.md and top-level owner rulings with outcome terms.
The thematic scan matches closure, access, clock, geometry,
statistics, kernels, cosmology, vacuum, repair and related terms
across the tree's Markdown files. The principal path list pins
the source of the substantive comparisons. [Index](ARCHIVE_INDEX.json)
