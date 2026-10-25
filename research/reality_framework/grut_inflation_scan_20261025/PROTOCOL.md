# Fixed-coupling inflation-scale test

Codex / AInstein, 25 October 2026. Written before the scan. Parent: b5174e5;
review: Claude 2793dab800794e78ebd74c43a512eb15622a7532, issue #13,
comment 6400868336. Codex takes the production scan; Claude retains review.

The question is whether the existing X/Y/Sigma sector has a regular,
productive release with a sufficiently large charge fraction under the
supplied instant radiation release. No new phase operator, coupling change,
target adjustment or erasure prescription is introduced.

## Inputs and preparation

M_Y=1e12 GeV, h=1e11 GeV, mu=5e11 GeV, kappa=1e-3,
lambda_X=1e-3, lambda_Y=lambda_Sigma=1, e=1e-7.
These are fixed classical tree coefficients. For this calculation their
stated reference scale is M_Y; no RG-improved release or UV matching is
claimed. The review's marginal high-field perturbativity remains a limit.
M_reduced=2.435e18 GeV, g_*=g_*s=106.75, and
Y_Q=4.8294009076988175e-9 are fixed.

Inflation adds -c_H H_i^2 |X|^2. At release this term is removed, the fields
have zero physical velocity, a=1, H=H_i/a^2, and phi is held at zero.
The physical potential, including both mediator interactions and all its
quartics, is retained. Initial amplitudes are obtained from the full
inflationary minimum. The exact aligned preparation is a control.

An optimistic supplied phase mismatch is also tested. In a gauge with
X positive and Y negative, Sigma=-|Sigma| exp(i delta). Thus the kappa
phase stays aligned while the mu source phase is displaced. This is an
additional initial condition, not a phase-selection mechanism in the
baseline action. Both CP signs must give opposite charge.

## Scan and adaptation

For each c_H in {1,10}, begin with 13 logarithmic H_i points from
3e11/sqrt(c_H) to 4.7e13 GeV. At each scale use delta=j*pi/12,
j=1,...,11, plus aligned controls. Hold the aligned radial amplitudes
fixed when supplying this displacement. Extend candidates with radial
minimization conditional on their supplied phase as a preparation control.

First integrate to a=1000. Record window averages, not an arbitrary
oscillation endpoint. Where results approach or cross yield, charge
fraction or phase-sensitivity thresholds, refine H_i and phase and run
ensembles with small phase seeds. A finite scan is not a mathematical
exclusion of every real H_i or initial condition.

## Quantities and gates

- Species charges N_i=a^3 n_i/H_i^3 and the compact Gauss residual
  N_X+3N_Y+6N_Sigma.
- N_eff=N_X+3N_Y=-6N_Sigma as well as N_X separately. Their time variation
  and late Y mixing must be exposed before calling either a final yield.
- The X-only quartic charge fraction epsilon_X=N_X/N_circ(E_X), with
  N_circ=2 sqrt(2 lambda_X) [E_X/(3 lambda_X)]^(3/4). This is a condensate
  diagnostic, not a complete particle census. Y/vector debris can worsen
  the leftover condition.
- Y_X=N_X H_i^3/s_i, Y_eff=N_eff H_i^3/s_i, with entropy from instant
  visible reheating. Compare both with the fixed Y_Q; do not refit it.
- Changes under delta -> delta+1e-9, compared against tighter numerical
  tolerances. Sensitivity of individual chaotic trajectories is not a
  converged ensemble mean.
- Y energy and its displacement from the late forced solution at a=1000;
  its approximate adiabatic number is an input to a decay calculation.
- Gauge Higgs norm sqrt(2 sum q_i^2 |z_i|^2), vector mass, and its varying
  background. A homogeneous gauge direction is not a longitudinal vector
  occupation number. No constant-mass abundance formula is silently
  promoted to a calculation in this Higgs history.

A candidate must be regular, reach the required net yield and have a
charge fraction above the review's optimistic homogeneous-equilibrium
threshold of about 0.19. The unprocessed cold remainder instead needs
about 0.99. A value near 0.01 counts only with independently demonstrated
erasure, which is absent here. Passing these necessary diagnostics alone
would not establish viability, phase preparation or acceptable relics.

## Numerical controls

Use conformal canonical fields with explicit eighth-order adaptive
integration. Check forces against finite differences, the gauge identity,
the conformal Hamiltonian's explicit-time-work ledger, CP reflection,
and an independent SciPy DOP853 implementation. For regular candidates
check tolerance, late windows, phase sampling and radial preparation.
For chaotic cases report ensembles and sampling error; do not certify
individual endpoints. Pinned review code and supplied reference results
are credited separately from original calculations.

No lattice fragmentation, finite-k vector production, Y decay,
thermalization, wall capture, string network, full stability spectrum,
observer dictionary or vacuum-density calculation is represented as done.
