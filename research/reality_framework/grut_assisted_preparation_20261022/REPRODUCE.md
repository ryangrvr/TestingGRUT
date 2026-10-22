# Reproduction

```bash
cd research/reality_framework/grut_assisted_preparation_20261022
python3 code/verify.py
python3 code/verify.py --rerun
```

Requires Python 3.12,numpy 2.3.5,scipy 1.17.0. Exact environment in
ENVIRONMENT.json. No network or new dependency installation is needed.
Optional `code/plot.py` uses matplotlib 3.10.8. This is same-author internal
reproduction, not an independent scientific implementation.

The verifier checks an exact payload file set and SHA256 hashes, then reads
the five raw probability CSV files and recomputes normalization, positivity,
boundary mass and moments. `calibration.csv` records direct and eliminated
energies for the twelve declared c/N choices. No residual fitting occurs.

The optional fresh calculation runs in a temporary directory without mutating
this packet. It requires all 27 controls to pass, matches the complete result
leaf schema and compares a declared scientific projection at 2e-6 relative
tolerance with zero absolute allowance. Tiny physical densities cannot pass
by being replaced by zero. Runtime, solver counters, residual magnitudes and
tiny boundary probabilities are controlled by gates rather than bitwise agreement.

Regenerate in a separate working copy, since these overwrite generated files:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/calculate.py
python3 code/plot.py
python3 code/make_manifest.py
python3 code/verify.py
```

The calculation uses a probability-conserving sparse generator and matrix
exponential, not a time-stepping Monte Carlo fit. Its continuum/domain controls
do not complete reheating, validate the stochastic truncation at all masses,
or turn the dated parameter sensitivity exercise into a likelihood.
