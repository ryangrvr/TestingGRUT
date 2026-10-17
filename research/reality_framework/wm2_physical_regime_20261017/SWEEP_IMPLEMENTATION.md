# Executed production sweep certificate

The frozen SWEEP_CERTIFICATE.md specification is implemented in
`certificate.py` and `solver.py`. Sources are pinned in EXECUTED_SOURCES.json
after controls and before the two primary workers. No primary source changes
are permitted afterward.

`residual_enclosure` forms R=X-KZ and w_0=K^dagger R directly. Complex
matrix-product/subtraction bounds use absolute gamma reductions, including
underflow pads; positive stored-array norms use S8's checked S7 kernel.
The unknown physical prefix Gram correction is bounded from the sum of
checked local isometry defects. Every in-range defect is charged linearly.
Losses combine in square norm only in structurally nested one-way batches.
Direction/consecutive-index guards run on each eligible gate. Reversals,
boundary operations and physical steps begin new batches.
For a left sweep, adjoint-oriented matrices represent the complex-conjugate
suffix state; conjugation preserves the loss, in-range and Gram norms.

SWEEP_CONTROLS.json contains six direct dense nonorthogonal-prefix tests,
the explicit nested-projector/commutation diagnostics, and the coherent
return counterexample. It also reruns S8's independent operator/positive
arithmetic controls, six complete small physical wavefunctions, and the
full original N=3,L=576,T=36 wavefunction (333,505 amplitudes).

The physical reference error is .01359811, below the .19647972 account.
It is larger than S8's reference discrepancy because the new protocol also
uses its separately frozen square-sweep drop allowance. The in-range account
there is about 1.09e-9. This is not an equal-error speed comparison.
The complete control worker takes about 157.96 seconds; the physical
reference portion takes about 144.08 seconds, separately from primaries.

This is a production-arithmetic enclosure under the declared IEEE/BLAS
assumptions, not formal validation of the compiler/hardware stack. Direct
comparisons do not constitute external mathematical review.
