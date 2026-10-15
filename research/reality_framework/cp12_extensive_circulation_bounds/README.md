# Extensive CP-12 circulation bounds

Read [REPORT.md](REPORT.md), [DERIVATION.md](DERIVATION.md) and [CHAT_REVIEW.md](CHAT_REVIEW.md). Parent: `36138c6140062b20c60864bce12a893c91dca8eb`. The canonical Hamiltonian, couplings and initial neutral ray remain unchanged.

This is an analytical probability/energy gate with independent exact-algebra and small complete-H-action controls. Supplied graph comparators establish energetic availability and the limits of moment-only selection. They are not formation runs from a supplied lattice.

With Python 3.11 or later, install [requirements.txt](requirements.txt). From this directory:

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python verify.py
```

The verifier freshly executes all four calculations in temporary output directories, checks their results against the saved scientific values, and validates source/output hashes. No prior wave snapshots, private pickles, source sessions, network download or large Hamiltonian are required.

To reproduce the individual files and figure:

```bash
python bounds.py
python checks.py
python prepared.py
python hosts.py
python plot.py
python verify.py
```

`verify.py --record-manifest` deliberately replaces the hash manifest after reviewed changes. Exact rational checks and analytic inequalities are separate from floating amplitude/norm calculations. There is no formal interval-arithmetic roundoff certificate.

The prepared-state controls retain every outgoing transition of one full unbounded-H action and use no amplitude pruning or electric cutoff. They compute energy, variance and local onset coefficients, not long-time dynamics. Large-host rows are exact graph/BFS calculations plus proved factorized quantum moments, never large-H evolutions.
