# B1-01 — derive, freeze, test, conditionally transfer

P1 is closed as physical matching unresolved. This new bounded test does not
extend a parent cutoff. The [derivation](DERIVATION.md), card, implementation and
analysis must be published before new corrected massless or massive histories.
The earlier M0-02 and P1 histories are known; this is not an untouched holdout.

Freeze one finite-volume dictionary: ghost endpoints, the standard neutral
boundary sector, the published bulk mass prescription, exact integrals of the
prescribed Gaussian over link cells, removal of the explicit alternating
c-number, and boundary weight at Phi=1/4. This choice is motivated by Gauss law
and quadrature, not uniquely forced by the finite Hamiltonian. No coefficients,
physical length, time scale, smearing width, mass or strength are fitted.

Use N=10,14,18,20 at m_R=0 and a refined N14 solve with tighter eigen tolerance
and two propagation substeps. At each N compute both the new complete instrument
and the unchanged old raw instrument for replay. Construct analytic massless
Gaussian histories for the actual cell smearing and for the ideal smooth target,
including the boundary mean and all ordered phases. Use 8192 modes and double to
16384. These mode sums are analytic reference accuracy checks, not extensions of
the interacting parent Hilbert space.

Separate first-readout and sequential TV, the measured disturbance of second
readouts, and the two dictionary comparisons (lattice-to-cell and cell-to-target).
Report means, covariances and the gap as diagnostics. Do not promote their
agreement to history matching. Also record the signed four-component diagnostic
from DERIVATION.md; its data-matched Gaussian auxiliaries are not the accepted
model and are not transferred.

Acceptance thresholds are fixed in the card. At N20 require first and sequential
target TV <=0.001, matched-cell TV <=0.001, cell-to-target TV <=0.00025,
disturbance difference <=0.001 and N18-to-N20 step <=0.00025, together with all
reference/solver/replay gates. No criterion is replaced by an auxiliary-summary
result. Preserve every miss, resource cap and exception.

Only if all massless gates pass, automatically execute the frozen rule at N20,
m_R=0.5, with no retuning. Its comparator is the already archived P1 E14/q128
smooth-target parent. Record first and sequential differences and keep the
interacting parent remainder and overflow unresolved. A finite comparison pass
is not a physical continuum certificate.

One BLAS thread, N<=20, 180 seconds per subprocess, 600 seconds total, including
the conditional transfer. No N22, E16, adaptive alternative dictionary or new
Gaussian massive model. Require a clean committed execution tree, preserve raw
probabilities, checkpoints, provenance and source/input hashes, and write outputs
outside the checkout. No merge, release or Zenodo upload.
