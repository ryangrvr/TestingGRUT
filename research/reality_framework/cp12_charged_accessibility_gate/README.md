# CP-12 charged-accessibility gate

Start with [REPORT.md](REPORT.md), then [DERIVATION.md](DERIVATION.md).
Independent review instructions are in [CHAT_REVIEW.md](CHAT_REVIEW.md).

The frozen source parent is a2f447dc37ea0a71f62c6e990a4630d727ebd85e.
This packet is additive research. It changes no physical-law file.

## Exact core

From this directory, using Node.js with BigInt support:

~~~sh
node verify.mjs
node exact_onset.mjs
node analytic_bounds.mjs
~~~

The exact onset program uses only integer arithmetic until it reports decimal
approximations. It enumerates all 960 minimum-order configurations via eight
permutation orbits, directly checks all 960 labelled forward amplitudes, and independently reverses every representative amplitude calculation.
The verifier also checks packet/dependency Git blob hashes and 680 unitary
remainder sanity cases. Those sanity cases are independent spin systems,
not CP-12 formation simulations.

## Numerical audit

The saved numerical commands completed in this session. The reusable wrapper
below was assembled after the container disconnected and has not itself been
rerun. A full packaged reproduction-suite PASS is not claimed.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit.py > NUMERICAL_RERUN.json
~~~

This builds fresh models unless --allow-cache is supplied. The N=5 closure
has 39,969 signed S5 states; it can take several minutes. To run only the
smaller initial-weighted spectral benchmark:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit.py --spectral-only
~~~

Compare against NUMERICAL_AUDIT.json with floating-point tolerances. Residuals
near machine precision can vary with the eigensolver, so compare the resulting
probability enclosures as well as raw residual numbers.

No wave snapshots or graph host files are needed. Python dependencies are
listed in requirements.txt and imported from the immutable parent research
tree. The numerical bounds control electric-space truncation analytically,
but do not certify machine rounding with formal interval arithmetic.
