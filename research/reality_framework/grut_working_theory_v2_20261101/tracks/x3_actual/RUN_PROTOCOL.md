# X3 execution protocol, 5 November 2026

Continuation of `SCOPE.md` at `9f71d34`. The existing Codex claim is resumed on
`codex/grut-x3-closure-20261105`; all earlier research tracks remain read-only.

## Input issue found before new kinetic runs

The archived H55 diagnostic contains both per-species energy and number density.
Its N1 pair of moments fails the minimum-energy bound for a free, two-helicity
Fermi gas. This is a failure of representation at the proposed kinetic interface,
not an exclusion of the underlying field theory. We will compute the massive
bound at occupations 1 (physical Pauli bound) and 1/2 (the dephased saturation
prescription), cross-check its normalization, and preserve this failure.

## Conditional continuation, not a repair of the original point

1. Use the final *recorded* H55 row at t=2.9980937. Never extrapolate it to t=3.
2. Compare it with the saved H55 proxy handover that reused the H84 lattice.
3. The matched comparison uses the old kinetic-reference species energy shares
   and hard/thermal SM split. A separate actual-species variant uses H55's own
   energy shares. These are supplied spectra, not reconstructed quantum states.
4. On the destination momentum grid directly, fill a sphere with occupation 1/2
   in each of two helicities and one partially filled boundary bin to match each
   supplied energy exactly. Report the resulting number densities and their
   changes. A 1/4-occupation sphere is a spectrum sensitivity control. No clipping
   may silently discard energy. No intermediate number-only rebinning is used.
5. Evolve the unchanged dephased rates, couplings and natural-line prescription
   of `e2e_kinetic`. Compare 3200, 6400 and 12800 logarithmic bins; halve the three
   adaptive timestep limits in a separately compiled copy for a time control.
   This source modification is a numerical control only and must be hashed.
6. The engine has no separate spectator/gradient transport channel. Keep the
   inherited assignment of scalar gradient energy to the chi condensate explicit.
   Keep the unrepresented spectator energy and printed-column rounding residual
   in a separate splice ledger. Check conservation within the represented network
   and the full-source Hubble jump separately; do not call the full splice passed
   merely because the represented energy is exact.
7. Report K=a^2 H, its late drift, duration, SM reheating and raw run artifacts.
   No root, abundance or empirical point will be accepted from these runs.

## Fermion controls

Check the actual saved moments against the exact Fermi variational bound. Test
the inherited quantum line integrator with homogeneous pulse evolution against
an independently integrated two-level Hamiltonian, and separate vacuum,
unitarity, temporal and momentum errors. A full-3D controlled toy propagation
test, if performed, is only an algorithm test, not the H55 cosmological solution.
Quantify the resolution needed by the reported cosmological momenta. Missing
transverse dynamics, missing H55 calibration and the nonrepresentable moments
remain distinct from numerical errors in a line solver.

## Population accounting

Provide separate core, free-X, anti-X, coherent and neutral-reservoir variables,
with exact charge/energy balances for forward/reverse conversion and neutral
decay. Do not manufacture a population or supply missing collision rates by fit.

End with a self-contained working-theory update and a portable continuation
record, distinguishing executed results, inherited results and open obligations.
