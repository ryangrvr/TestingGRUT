# CP09 bounded-below moving clock: full apparatus simulation

2026-10-11. Frozen CP09 unchanged; explicit apparatus extension, not CP10 or a Stage3 card. This replaces the previous vP ideal clock by P²/(2m), while retaining one prescribed source-pointer protocol. Mass and numerical resolution variations are disclosed apparatus-limit tests, not a hidden fundamental candidate search.

## Hamiltonian

H_total=P_c²/(2m)+H09-epsilon f_delta(q) B-zeta f_delta(q-.4) A Y_r.

Clock is a particle on a circumference-4 circle, simulated by a Fourier spectral grid. f_delta is the unit-area square window of width .15. Epsilon=zeta=.001. B,A are frozen disjoint relational-role operators, with norms <=1. H09>=0, so H_total >=-(|epsilon|+|zeta|)/delta. Unlike vP this full clock Hamiltonian is bounded below, including its finite-strength interactions. System-reference access is still global and supplied.

Prepare a product of CP09 ground state, pointer |0>, and clock Gaussian centered -.7 with position standard deviation .08 and mean momentum m (mean speed 1). This does not use CP09 energies to manufacture the clock. Total external duration 1.5; interactions are fixed position potentials, not externally timed pulses. The clock is periodic, not a fundamental monotonic time operator; external time remains in the simulation.

## Full simulation and results

Strang splitting alternates half-step exact free clock/system evolution in momentum and system-energy basis, with exact local source and pointer interactions in position basis. Clock wavepacket, matter/geometry and pointer are all retained. Separate source signs provide the actual pointer-X response estimator. No history state, injected ideal kicks or postselected geometry is used.

Ideal response .052588513539. Results:

| mass | grid | step | pointer response | absolute error |
|---|---:|---:|---:|---:|
|20|256|.005|.0476627993|.0049257142|
|80|512|.005|.0578050759|.0052165624|
|80|512|.0025|.0578049879|.0052164744|
|80|1024|.0025|.0597041453|.0071156318|
|320|1024|.0025|.0539540783|.0013655648|
|320|2048|.00125|.0539630583|.0013745448|

The mass80 case exposes grid sensitivity; time-step agreement alone was insufficient. It is retained, not treated as converged. At mass320 the combined refinement changes the response by 8.98e-6, much less than the remaining ideal-response difference. This is one numerical refinement comparison, not a rigorous discretization bound or arbitrary-mass convergence theorem. Sharp windows make grid treatment especially important.

All norm errors <1e-10. Maximum energy drift over these controls is 3.6e-7; refined mass320 drift 1.2e-10. The Hamiltonian analytically conserves total energy, while finite splitting introduces numerical drift. Initial total energies are 11.0713,40.3389,160.1558 for masses20,80,320. Increasing mass at fixed speed costs kinetic energy approximately m/2. No free ideal clock limit is claimed.

## What was overcome

A bounded-below, finite-energy, product-prepared autonomous apparatus produces an actual pointer response close to the ideal CP09 relational response. For mass320 it differs by about 2.61%. This improves on the unbounded ideal clock and on a supplied stationary correlated history.

What remains: exact response recovery or a controlled asymptotic error bound, robust numerical convergence across all controls, reliable clock reading/instrument calibration, local realization of global A,B, physical spacetime geometry, continuum scaling and gravity. The residual error cannot yet be separated into finite pulse width, dispersion, back-action and numerical error with a certified budget. No assertion that all limitations were overcome.

This is conventional quantum apparatus dynamics conditional on supplied CP09 and source/detector couplings; not a discovered law of reality. No target geometry or causal-link cap inserted. Code prints EXECUTION_COMPLETE rather than a scientific passing verdict: assertions check norm only, and physical differences are reported without suppressing failures.

Run python research/reality_framework/cp09_bounded_clock/check.py with NumPy/SciPy. Results include every mass/grid/step control. Author simulations, not independent validation or experiment.
