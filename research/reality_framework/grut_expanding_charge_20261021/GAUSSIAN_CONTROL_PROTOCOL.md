# Additional Gaussian-spectrum control, version 1

2026-10-21, after exposure to the first four frozen-suite outcomes, before any dense
Gaussian run. This is an explicitly exploratory preparation comparison. The original
eleven-row suite, its source freeze, and every result are retained without replacement.

Motivation: finite 64-mode perturbations are not the same joint spatial ensemble as
a Gaussian random field with a populated power spectrum. One main seed selects one
vacuum while the other two form both. A finite mode set can have appreciable global
cross-correlations between S and X perturbations, even when each has zero mean. These
may source a coherent S zero mode; they are not an explicit Z2 bias in the action.
The result must not be attributed to a discrepancy between integrators before matching
the preparation sufficiently closely.

Keep the validated, frozen C++ evolution, diagnostics, expansion, kappa=2, amplitudes,
T=48 and numerical parameters. A separate wrapper reuses the exact frozen solver and
reads a supplied initial-noise field. Generate THREE independent Gaussian white-noise
arrays with NumPy PCG64, Fourier filter by exp(-k^2 ell^2/2), remove the zero mode and
normalize each realization to unit rms. This is still classical preparation data,
not a derived thermal spectrum. It is statistically like the smooth Gaussian class
used by Claude, without matching their particular seed or claiming bitwise IC identity.

Run, sequentially after the original suite:

| N | dx | dt | seed |
|---:|---:|---:|---:|
| 192 | 0.4 | 0.06 | 101 |
| 192 | 0.4 | 0.06 | 202 |
| 288 | 0.4 | 0.06 | 101 |

The Gaussian source is a new wrapper, not a replacement of the source certified by
SOURCE_FREEZE.json. Validate its initialization by charge, means, rms, covariance,
and a mirrored small-grid pair, and require a direct comparison against the frozen
implementation from identical explicit initial fields before production. Record a
source hash before launch. A cumulative 6-minute launch cap applies to these three
additional runs. No new observational target is adopted.

Use the original formation, topology and conservation definitions. Report differences
in finite-time behavior and the exact preparation covariance. These three runs cannot
establish a continuum or cosmological percolation probability, nor isolate every reason
the two perturbation families differ. No exponent near 10^-23 is tested.
