# E2-01: fixed-image memory diagnosis and continuation

Freeze this protocol, derivation and implementation before evaluating its
factorial controls, new bounds, memory scheme or new timings. E1's parent and
old autonomous outcomes are already known; they are calibration/replay, not
untouched holdouts. No result here can change E1's decisions or M1's continuum
status.

## Fixed inputs and comparisons

- Reuse the exact E1 E=12, 1280-state parent, m/g=0.5, Lg=5, theta=0,
  boundary phi/sqrt(pi)=1/4, 128 spatial quadrature nodes and its stored 72-column
  isometry. No new anchor, retained state, coupling fit, mass or boundary choice.
- Strong detector only: c=2, smearing width 0.8, phase zero. Its compressed Weyl
  contraction and absorbing overflow are unchanged.
- Diagnosis: old timings (0,0.7,1.9) and (0,0.6,1.5,2.4). Compute all four
  switches in DERIVATION section 3 and the absorbing-reset control. Archive all
  probabilities, signed factorial terms, per-prefix retained and complementary
  norms, outgoing detector leakage and detector return, and autonomous residuals.
- New timings, frozen here: (0,0.43,1.37,2.61) and
  (0,0.31,1.07,1.83,2.79). Test the unchanged autonomous scheme, exact parent,
  absorbing-reset control, and memory continuation. These timings have not been
  evaluated when this card is committed. They are a new numerical configuration,
  not independent observational data or a wholly unfamiliar physical model.
- Run the memory scheme at n=32,128,512 substeps per nonzero interval on all four
  histories. Keep every q coordinate and every singular value. Report all levels;
  n=512 is the predeclared accuracy/bound decision level. No adaptive enlargement.
- Report the autonomous bound in equation (6), including its separate initial,
  propagation and instrument budgets, without substituting it for E1's bound.
- Independent small-matrix verification: direct SciPy exponentials versus the
  low-rank exchange step; commuting limit; full instrument completion; local
  error bound; a deliberately expansive detector rejected by contraction check.

## Decisions fixed before execution

1. Implementation prerequisites: retained dimension 72; parent dimension 1280;
   isometry/unitarity, block reconstruction, Hermiticity and preparation residuals
   <1e-9; minimum instrument-defect eigenvalues >=-1e-10; probability normalization
   error <1e-9 and no probability below -1e-10. Never clip or renormalize ordinary
   outcomes. Abort interpretation if these prerequisites fail.
2. Replay: full exact old parent and both-off autonomous history TV <1e-10 from
   the archived E1 outputs. Factorial reconstruction error <1e-12, including
   overflow. Strong coupling and old misses remain unchanged.
3. Diagnosis is descriptive: no preferred sign or dominant component is chosen
   in advance. Rank all prefixes by unnormalized detector-leak norm and residual
   norm; also record conditional norms when branch weight is nonzero. Report
   cancellation, rather than allocating positive causal percentages.
4. Memory accuracy: maximum TV at n=512 <=0.001, separately for old and new
   timings. Method improvement: that maximum <= one half of the corresponding
   autonomous maximum. Measured accuracy and improvement are different gates.
5. Memory bound: equation (9) <=0.001, separately for old and new timings. Every
   bound must dominate measured TV within 1e-10 numerical slack. Use the full
   Frobenius-commutator upper bound, even if loose. Do not call a floating-point
   pass an interval-arithmetic certificate or a continuum result.
6. Autonomous diagnostic bound: equation (6) must dominate measured TV within
   1e-10; separately report whether it reaches 0.001. A miss is retained.
7. Report 72 retained plus 1208 memory amplitudes, matrices used, wall time and
   every reused component. Full auxiliary memory is accuracy without state-space
   compression. Exact parent diagonalization is a reference/diagnostic cost,
   not an input to the low-rank splitting propagation.

## Resources and failure handling

One parent reconstruction, one BLAS thread, 600-second external wall cap. Output
must be outside a clean committed checkout and must not already exist. Preserve
stage outputs and failure reports. A numerical exception/timeout is not a physics
verdict. No basis increase, new detector strengths, large grid, merge, release or
Zenodo deposit. Boundary/mass-dictionary review remains separate and open.
