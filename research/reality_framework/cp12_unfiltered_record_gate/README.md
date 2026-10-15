# CP-12 unfiltered relational-response gate

Read [REPORT.md](REPORT.md), [DERIVATION.md](DERIVATION.md), [THEORY_STATUS.md](THEORY_STATUS.md) and [CHAT_REVIEW.md](CHAT_REVIEW.md). This additive packet follows immutable parent `7d414961a56a49df1265b9501cd305331354ffc6`. No frozen law, coupling, main initial ray or canonical branch is changed.

The main response calculations apply **no initial loop filter**. Final circulation conditioning is a declared readout. Nonselective recording and selected branches are separate counterfactual controls with explicit work and information effects.

Install [requirements.txt](requirements.txt) with Python 3.11 or later. From this directory:

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python verify.py
```

The verifier uses a fresh temporary model cache and reconstructs its own N4/N5 wave snapshots. It reruns all three main unfiltered responses, both initial-record comparisons for each case, the independent raw-physical short-lag controls, both leading phase–circulation coefficients and the 24 finite unitary dilations. It checks hashes and saved consistency for the larger capacity scans and the two longer mixed-reference evolutions; it does not silently rerun them. `verify.py --record-manifest` intentionally replaces the source/output hash manifest after reviewed changes.

To independently reproduce every saved calculation, first reconstruct the already-published snapshots. A trusted private matrix cache is optional and only accelerates builds; do not load untrusted pickles.

```bash
python ../cp12_relational_record_gate/analyze.py --n 4 --q 4 --budget 256 --snapshot-dir /tmp/cp12-unfiltered-waves
python ../cp12_relational_record_gate/analyze.py --n 5 --q 5 --budget 128 --snapshot-dir /tmp/cp12-unfiltered-waves
python ../cp12_relational_record_gate/analyze.py --n 6 --q 6 --budget 128 --snapshot-dir /tmp/cp12-unfiltered-waves
python ../cp12_relational_record_gate/analyze.py --n 6 --q 12 --budget 10 --times 0 0.4 0.6 1 --snapshot-dir /tmp/cp12-unfiltered-waves
python audit.py --n 4 --budget 256 --wait 5 --snapshot-dir /tmp/cp12-unfiltered-waves
python audit.py --n 4 --budget 256 --wait 10 --snapshot-dir /tmp/cp12-unfiltered-waves
python audit.py --n 5 --budget 128 --wait 10 --snapshot-dir /tmp/cp12-unfiltered-waves
python references_audit.py
python short_lag.py --snapshot-dir /tmp/cp12-unfiltered-waves
python capacity.py --snapshot-dir /tmp/cp12-unfiltered-waves
python independent_checks.py
python plot.py
python verify.py
```

The `source=True` builder flag selects the exact subgroup of vertex permutations fixing source label zero. It is a computational basis choice, not a physical source term. Explicit source kicks appear only on the declared diagnostic paths.

The finite-lag intervals and bounded short-lag coefficient intervals control the full unbounded electric Hamiltonian through outgoing residuals. Their numerical evaluations are ordinary floating point, not formally rounded interval certificates. Extra short-lag display samples illustrate local asymptotics and are not long-time or large-size extrapolations.

Q=N spanning circulation is a ring. The new capacity theorem states what a richer size/density sequence and electric budget must admit before they can test an extended circulating network with matter. Neither this packet nor its small systems establishes spacetime, relativity, gravity or a complete theory of reality.
