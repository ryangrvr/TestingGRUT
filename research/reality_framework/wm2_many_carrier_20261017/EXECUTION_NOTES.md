# Execution and provenance

## Frozen boundaries

* Parent IG1: `650104cd9e2b3f0360b2f3eb3ad7cd96cb4422b5`.
* D1 charter/source freeze: `f4e85e7d58a93802692a441ac21234a47211eb42`,
  before finite D1 controls. Analytical discovery had begun; not a blind study.
* Predictivity addendum: `907272a7720a7fd8e99eb5b1fc5f67e77cf58e17`, before
  evaluating tau=32 holdouts. It responded to the user's broader predictivity
  directive without changing D1's physical law or prior bounded stop.

No formation worker, six-event target or large-receiver trajectory was run.
No external measurement or independent mathematical review occurred.

## D1 controls

The first D1 suite passed. Before final publication an independent two-site
tensor-matrix implementation was added to check the matter entries, and the
eigensolver starting vector was made deterministic. The suite was replayed;
the final archive is CHECKS.json. No time, state, error or resource gate was
relaxed. Final worker: 3.177845125 seconds, 121.9921875 MiB peak resident memory,
under the frozen 60-second / 1 GiB ceilings.

Fixed-path Fock and full-geometry builders share the matter_actions routine.
That comparison alone is not an independent implementation of the matter law.
The separate 12x12 local tensor construction does not use matter_actions;
all nine two-site number-sector matrix restrictions agree entrywise (maximum
difference zero). This is an independent code control, still by the same author.

Small complete physical models have N=3,4 and N_A=2,N_B=1. They include every
label ordering, storage mask, allowed two-mode occupation and internal state
in that number sector. Their electric kernel is eliminated by its exact
unique tree solution; all physical states are retained. Capacity enlargement
is exercised: their maximum electric magnitude is four, outside UC1's old
three-unit cap. Both reverse geometry and matter transitions are present.

The N=4 common-reference control deliberately has L=0 buffer width. Its
evaluated bound is loose but exercises full-state error, prefix loss and
nonzero storage feedback. It is not a numerical observation of the primary
large-N macroscopic-buffer regime. Samples at N=256,... evaluate analytical
constants only; no wavefunction at those sizes is propagated.

The very small distant commutator control is at floating-point roundoff scale;
it supplies no high-precision numerical proof of the all-size bound. Other
finite controls resolve nonzero hopping slopes, norms, flux changes and
storage errors. NumPy/SciPy complex128 is not outward-rounded arithmetic.

## Prediction calculation

Tau=32, the couplings and outcome channels were frozen before calculating
the holdouts. Twenty-five rational Bessel-series terms with alternating
remainders give exact rational mathematical bounds, then outward decimal
rounding gives PREDICTIONS.json. Independent quadrature and a 64-site matrix
evolution are ordinary-double cross-checks.

While completing the uncertainty derivation, the Hamiltonian-error coefficient
was conservatively changed to 66 to include the clock-calibration product.
The prediction values, time, couplings and finite-control gates did not change;
no experimental data existed. The final two-channel systematic budgets use
that coefficient. Additional exact-rational finite-boundary threshold assertions
were included before the final prediction archive was regenerated.

Exact rational certification concerns only these free-channel numerical target
values and their stated series/boundary evaluations. It does not certify the
S9, IG1 or D1 proof, a device, or physical reality. The proposed empirical error
limits are requirements for a future test, not observed device specifications.

Reproduction preserves the archives; hashes pin the executed scripts and all
upstream sources. Nothing in the canonical earlier GRUT ledgers was regraded.
