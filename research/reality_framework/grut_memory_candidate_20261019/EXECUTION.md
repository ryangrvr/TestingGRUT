# What was executed

All runs use local existing Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.
Matplotlib 3.10.8 generates the scientific figures. No package installation,
paid resource, empirical data or indefinite worker was used. The law is not
adjusted by an optimizer. An owner-directed independent code review prompted
the matched-control amendment below; this is not experimental validation.

## Core and attacks

1. The original `calculate.py`, now archived byte-for-byte as
   `calculate_original.py`, executed seven finite-system runs: L=3,4,5 at the core
   point, and L=4 controls raising the transverse field to F=1, removing feedback,
   freezing the memory, and halving softness. Exact spectral evolution
   propagates the supplied states. The complete run took 11.79 seconds;
   all 13 originally registered acceptance conditions passed.
2. `challenge.py` independently assembles the L=4 Hamiltonian through
   Kronecker products and propagates it with sparse Krylov exponentials.
   It then executes the predeclared independent-cell, bare writing,
   unconstrained-perturbation and centered-OTOC attacks. The run took
   7.36 seconds. The three numerical-method checks passed; the six stronger
   scientific attack conditions failed and remain false in the record.
3. `algebra.py` checks the exact zero-softness domain-wall commutator and
   evaluates ordered-sector weights in the separate formation preparation.
   The conservation and initial-weight identities pass. Finite formation
   weights are reported without claiming a persistent organized state.

## Matched-control amendment after review cc24967

The review's `uniform_field_control.py` was executed with one BLAS thread.
It reproduces all six rows of the archived `control_results.txt` to the
printed precision. No writer scattering replay was performed in this amendment.

The amended `calculate.py` then executed ten runs in 7.31415 seconds:
the original seven and F=epsilon=0.05 at L=3,4,5. It holds every other
coupling, preparation, observable and time grid fixed. The original source,
protocol and results remain unchanged. The new source and
`control_protocol.json` identities are recorded in `results_with_control.json`.
The control protocol explicitly records that the review outcome was known;
the added control is not a blind preregistered prediction.

At every size the bit-table and independent tensor Hamiltonians agree,
and the initial retention second derivative is the same as M1's. The
uniform-floor control is no worse at all 201 sampled times at each size,
within 1e-10 tolerance; its minima and t=40 values are better. Best readout
differences are 0.000179931, 0.000015581 and 0.000015652 at L=3,4,5.
The numerical result rejects the facilitation-benefit attribution in this
prepared regime. Passing the original 13 conditions did not test that
attribution, and cannot establish it after this comparison.

`assemble.py` now includes the matched control in the main figure and
summary. `verify.py` checks both source histories, all ten probability
traces, the independent review's rounded values, the original-run replay,
and the retained negative mechanism outcome. It does not overwrite any
original physical execution record.

## Writer and numerical refinement

The separate `writer_protocol.json` adds an internal quantum source and
massive particle to complete the failed preparation obligation. It keeps
the core fields, couplings, softness and readout operations fixed. The
physical Hamiltonian is time independent.

`autonomous_writer.py` first executed the full 2048-position scattering
calculation with step 0.0002 and a symmetric local-Pauli core split.
That completed in 414.59 seconds. It stores the completed result in
`writer_results_initial.json`. The subsequent slower step-0.0001 replay
was interrupted with Ctrl-C before completion. The completed coarse
record and original source were preserved. This is an execution
interruption, not an erased physical failure or a completed refinement.

`execution_amendment.json` records replacing that unfinished replay by
`writer_replay.py`. It uses exact core matrix exponentials, combines the
commuting core/kinetic half steps, and uses equivalent strided writer
rotations. It executes step 0.0001 at 2048 and 4096 clock positions.
This provides algorithm/time and spatial refinement without changing
physical parameters, gate shapes, preparations or acceptance thresholds.
Completed values and wall times are in `writer_results_refined.json`.

## Numerical domain and claim boundary

The scattering stage approximates a particle on R by an FFT box from -16
to 16 only until time 0.06. Both directions of travel and every clock
position enter the spin partial trace. Initial, final and interaction
energies, norm drift, edge weight, momentum tails and residual gate-region
weight are recorded. No transmitted-only normalization is used.

Storage uses the unconditional reduced spin state and its exact core
evolution after the particle has left the interaction region. Residual
gate interaction and possible later return are not given a rigorous
all-times scattering bound. Small reported overlap and numerical
refinement support this approximation but do not certify it. In
particular the FFT box is not claimed to be the physical particle's
finite periodic universe during arbitrarily long storage.

The original core dynamics have no such carrier-decoupling approximation.
The full writer's actual coupled scattering result, the subsequent
storage approximation and the original exact-core regime are separate
levels of evidence.
