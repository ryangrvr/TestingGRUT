# What was executed

All runs use local existing Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.
Matplotlib 3.10.8 generates the scientific figures. No package installation,
paid resource, empirical data, indefinite worker or external review was
used. The law is not adjusted by an optimizer.

## Core and attacks

1. `calculate.py` executed seven finite-system runs: L=3,4,5 at the core
   point, and L=4 controls removing saturation, removing feedback,
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
