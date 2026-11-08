# P1-01 — physical dictionary and one further parent cutoff

This is a new bounded extension, not an amendment to M1, E1 or E2. Freeze this
document, card, derivation, point runner and analysis publicly before producing
new E14 Hamiltonians or histories. The old M1 E10-to-E12 result is already known;
these timings and the N20 comparator are not untouched holdouts. The extension
tests convergence of the already declared experiment, not a new prediction.

## Frozen experiment and construction

Use the unchanged M1 card: hbar=c=g=1, L=5, m_R/g=0.5, theta=0, Dirichlet
phi(0)=phi(L)=sqrt(pi)/4, neutral even-chain sector. Prepare the ground state of
each regulated interacting Hamiltonian as an approximation to the same continuum
ground-state preparation rule. Fix normalized Gaussian smearing sigma=0.8,
W=exp(2i F), F=integral f phi/sqrt(pi), phase zero, and times (0), (0,1), (0,2),
(0,1,2). No pointwise mean or variance matching.

Retain the unchanged normal-ordered Hamiltonian and exact oscillator Weyl matrix
elements from M1. The wrapper changes only its resource cap and labels. Higher
unoccupied modes supply the Weyl vacuum factor; no exponentiation of a compressed
field, survivor renormalization or fitted counterterm. Ground-state energies are
unrenormalized cutoff diagnostics; normalized histories and gaps are the targets.

Run exactly E12/q128 (replay), E14/q128, E14/q256. Preproduction basis counts are
1280 and 2894 respectively; E14 has 22 occupiable modes. The new dimension cap is
3000, explicitly larger than M1's 1300. No E16 or adaptive ladder. Reuse archived
E10/q128 and N20 without altering their preparations, dictionaries or numbers.

Check nested embeddings by padding occupation tuples with zero higher modes.
For E14/q128, rebuild E12/q128 using the unchanged code and compare the compressed
Hamiltonian and Weyl blocks. Check the parent-overflow decomposition in
DICTIONARY.md. Record embedded E12 ground-state overlap, E14 weight outside E12,
and the actual massive mean of F. These are descriptive diagnostics, not fitted
inputs or independent validation of the reused implementation.

## Frozen analysis

Always retain overflow as a distinct absorbing outcome. Compare complete
distributions, including overflow, with total variation (TV). Report every
protocol, not just the maximum.

1. Reproduce archived E12 histories within TV 1e-10. Check Hamiltonian Hermiticity,
   low eigenpair residuals, instrument contraction and total probability.
2. Compare E12-to-E14 differences against the already observed E10-to-E12
   differences, both per history and by their maximum. A decrease is a trend;
   do not infer an asymptotic power or extrapolated continuum value from it.
3. Apply the inherited necessary gates: E12-to-E14 max TV <=1e-4,
   E14 quadrature change <=1e-6, finest overflow <=1e-5, and N20-to-E14 max
   TV <=0.001. Take q128 as the primary point; q256 checks quadrature only.
4. Verify matching distance >= overflow (up to rounding). Report the gap
   between this necessary floor and the actual discrepancy. Do not condition
   the comparison on avoiding overflow.
5. Keep lattice discretization, boundary/mass matching, Hamiltonian truncation,
   quadrature, instrument truncation and solver diagnostics separate. A failure
   of a necessary prerequisite leaves physical matching unresolved. Even all
   numerical gates passing supplies no absolute parent-to-continuum bound.

## Resources and preservation

One BLAS thread; 450 seconds per point; 1050 seconds over all subprocesses.
Record complete points, exceptions, incomplete outputs and timeouts. Stop without
increasing a cap. Require a clean committed source tree; write new results outside
it; record source and input hashes, versions and elapsed times. Recomputing the
analysis of saved probabilities is not an independent physical implementation.

E2 remains the 72+1208-amplitude finite-parent reference, not compression. E1
misses and M0/M1 misses remain unchanged. No new boundary fit, lattice grid,
continuum certificate, distinctive GRUT law, merge, release or Zenodo upload.
