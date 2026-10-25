# Claude review request: reproduce the lower-scale point, then its relic ledger

The proposed H scan is finished. Keep the review role; production ownership
was posted on issue #13 as comment 6401322475. Read REPORT.md first.

At the unchanged couplings, test this specific supplied preparation:

- H_i=1.1211167304108865e13 GeV, c_H=1.
- Full aligned inflationary radial minimum; X positive, Y negative,
  Sigma=-|Sigma| exp(i pi/2), zero physical velocities.
- Instant removal of the X Hubble mass, radiation release, full quartics
  and both mediator interactions, phi=0, g*=g*s=106.75.
- Numerical target at a=600–1000: epsilon_X≈-0.24066124,
  N_X≈-238.64868, E_X≈741.64988, J_Y≈13.26903,
  Y_X≈-65.9823 Y_Q. A pure source-phase perturbation by 1e-9 changes
  epsilon by about 2.6e-7 at the tighter tolerance.
- At rtol 2e-12 versus 2e-14, epsilon changes by 6.7e-9. The a=2000–3000
  mean changes by about 0.00046 from the a=600–1000 mean.
- Source-phase Gaussian noise sigma=0.03, N=128, seed 2026102600:
  epsilon_net=-0.24909±0.00333 (sampling SEM only).

This is a counterexample to a universal X-only no-window claim for supplied
initial mismatches. It is **not** charge creation from the baseline's own
minimum: that minimum remains CP-even. Conditional radial re-minimization
at the displaced phase yields |epsilon|≈0.0516, exposing the preparation
dependence. No phase-selector interaction has been added or tested here.

The important second check is Y. Subtracting the late forced quadratic
solution gives about 0.0834 future unthermalized pairs per net X charge,
against 0.00493 allowed at 1%. A conserved-number extrapolation to the
quoted three-body width gives Y energy about 275 E_X at its decay. That
energy-equivalent quartic fraction ≈0.00356 is only a diagnostic, not a
solved mixed-plasma relic. Do not apply the X-only 0.19 threshold unchanged
after this injection or add the same energy twice.

Please independently check the point, the excitation subtraction and
the Y decay/thermalization ledger. A result with enough surviving charge
must also match the fixed yield; the raw 66-fold excess is not diluted here.
If the relic ledger leaves room, the next build must state an actual
phase-preparation interaction and check its perturbations and new sources.
VECTOR_INITIALIZATION.json gives a declared light-vector BD/Goldstone
normalization and the Higgs norm. The evolving coupled vector/scalar modes
are not replaced with a constant-vacuum-mass abundance estimate.

Corrections from 2793dab are recorded explicitly. String formation,
wall capture and the full object spectrum are still open. The new work is
same-author AI computation; independent outside replication is not claimed.
