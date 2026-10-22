# Verification actually performed

2026-10-22, one author. The owner's supplied review concerns 9d0362c and
independently hand-checks two earlier formulas; it is not a review of this
new packet. No sub-agents were used and no independent outside replication
of these new calculations is claimed.

## New calculation

`OPENBLAS_NUM_THREADS=1 python code/calculate.py` completes with **41 internal
conditions and zero failures**. All 34 scalar histories, two minimum-response
calculations, four refinements, four separately organized cosmic-time
controls and twelve pure-exponential phase-plane controls finish. The
scientific failure of generic field selection is retained as the result;
passing numerical checks is not evidence that this action describes nature.

| Check | Recorded outcome |
|---|---:|
| Integrated Raychaudhuri versus algebraic Friedmann H, all histories | max absolute log-H difference 3.94e-12 |
| Step/tolerance refinement of u, four controls | max sampled absolute difference 3.42e-12 |
| Step/tolerance refinement of H | max relative difference 3.47e-12 |
| Cosmic-time formulation versus refined N formulation, u endpoints | max absolute difference 5.55e-13 |
| Cosmic-time formulation versus refined N formulation, H endpoints | max relative difference 5.54e-13 |
| Scaled minimum response ODE versus separate double quadrature | absolute difference below 4e-15 |
| Pure-exponential endpoints | within the specified 2e-7 of analytic fixed points |
| Phase-plane stability | finite-difference Jacobian eigenvalues agree with analytic formulas within 2e-8; all real parts negative |

The cosmic-time formulation shares the same potential and calibration. It
tests numerical organization; it is not a blind independently authored
implementation. The minimum calculation resolves a scaled coefficient
R=-0.05015657 and multiplies by the tiny analytic mass ratio, rather than
subtracting rounded equal field values.

## Read-only verification and fresh replay

`python code/verify.py` checks the complete packet file set and hashes,
all 34 endpoint CSV/JSON records, the five raw histories' Friedmann energy
fractions and distance limits, independent quadratures of archived energy
and momentum identities, and analytic stability eigenvalues. On the
finite raw grid, the largest energy log-H quadrature discrepancy is
5.13e-10 and the largest momentum integral relative discrepancy is
5.01e-7; both satisfy the declared 2e-6 tolerance. These raw quadratures
have lower resolution than the adaptive evolution's auxiliary identity.

`OPENBLAS_NUM_THREADS=1 python code/verify.py --replay` copies the source
and input to a temporary directory and recomputes the new packet. It
compares **900 scientific JSON leaves** and **six CSV arrays**, requiring
identical key sets, columns and shapes. Maximum scientific relative
deviation is **zero** in this environment; all CSV comparisons pass.
Evaluation counts are computational diagnostics excluded from physical
leaf comparisons. Repeated leaves and checks are not independent physical
discoveries. This is same-code reproducibility.

The final manifest hashes 21 packet files and excludes itself. Source
syntax, JSON parsing, local Markdown targets and `git diff --check` were
checked. The figure was rendered, visually inspected, and its initially
overlong panel title shortened before finalization. No earlier scientific
packet was changed or replayed during this continuation. No new packages
were installed.

## Not performed

No inflation-to-reheating stochastic or quantum-state transport, spatial
gradient evolution, traceful thermal matter evolution, actual electroweak
transition, quantum corrections, data fit, likelihood calculation, held-out
observation, independent-agent review or all-actions theorem. Initial
temperature/field grids are explicit exploratory controls, not measured
initial conditions or a cosmological ensemble.
