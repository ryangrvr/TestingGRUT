# E1-01 — joint Hamiltonian and detector dressing

Freeze this card, derivation and implementation before computing the dressed scheme or the weak-detector histories. The old E=12 parent and strong-detector D1 outcomes are known. This is a method test within that finite parent, not independent continuum certification or a retrospective alteration of M1/D1.

## Fixed construction

- Parent: unchanged M1-01 bosonized Hamiltonian, m/g=0.5, Lg=5, theta=0, boundary phi/sqrt(pi)=1/4, E=12 free-energy cutoff, 128 quadrature nodes, dimension 1280. Reuse its frozen builder openly.
- Retained space: every free Fock state of that same parent with free energy <=6, dimension 72. No eigenvector chosen from the full parent is an input to the correction.
- Anchor: initialize with the bare retained ground energy. Apply equations (1)–(3) of DERIVATION.md, at most 12 iterations, stopping when successive anchor energies differ by <1e-12. Stop the construction if distance to spec(C) is <=1e-8. Record all iterations and norms. Diagonalizing C is allowed and its cost is counted; no computational speedup is claimed.
- Compare three 72-dimensional descriptions: bare A and PWP; the jointly dressed H_e and W_e with the dressed ground state; and H_e with that same state but the bare PWP detector. The last is an explicit detector-matching control, with exactly the same spectrum as the joint scheme.
- Target protocols: (0,0.7,1.9), (0,0.6,1.5,2.4), Gaussian width 0.8 and phase 0. Run c=2 as the known strong reference, and c=0.5 as a newly specified weak-detector regime. The weak experiment cannot replace a failed strong-detector verdict.
- Parent reference: diagonalize the full finite parent only after constructing J. Reproduce the archived D1 strong-parent histories before interpreting anything. All comparisons include absorbing overflow and never renormalize the ordinary probabilities.

## Frozen acceptance rules

1. Parent replay TV <1e-10; retained dimension exactly 72; isometry residual <1e-10; reduced Hermiticity <1e-10; defect eigenvalues >=-1e-10; distribution normalization error <1e-9. Ground residual <1e-8 and overlap with the parent's ground state must agree to 1e-8 in vector norm after phase alignment. Anchor convergence is required.
2. For each detector strength separately, the joint scheme's maximum history TV must be <= half the bare scheme's maximum and <= the Hamiltonian-only scheme's maximum. This is the method-improvement gate, not an absolute accuracy statement.
3. Separately require maximum joint TV <=0.001 for the finite-parent accuracy gate, and maximum bound from equation (8) <=0.001 for the certification gate. A measured pass with a failed bound is reported as such, not called certified.
4. Every computed bound must dominate its measured error up to 1e-10 numerical slack. Negative integral-Gram eigenvalues below -1e-9 or negative quadratic forms below -1e-10 are implementation failures. Small roundoff negatives may be clipped to zero and their extrema must be recorded.
5. Record preparation error, ground/gap errors, residual norm, norm of T, integrated dynamical-defect contributions, detector-defect contributions, all probabilities and per-branch bounds. A good ground energy does not substitute for any history gate.

## Resource and interpretation boundary

One 1280-state reconstruction and at most two full eigendecompositions of that order (parent and Q block), one BLAS thread, 600 seconds total. The parent cutoff is not enlarged. No fitted detector coefficient, changed preparation, new mass, altered boundary, increased retained cutoff or extra anchor is allowed in this registration. Archive exceptions and misses.

If the joint correction helps only in the weak regime, report that restricted result. If it fails, the eliminated sector cannot be replaced by this one-anchor autonomous approximation at the stated accuracy. Neither outcome addresses states above E=12. A successful finite-parent method would need a separate continuum test before updating M1-01's physical status. Zenodo remains paused.
