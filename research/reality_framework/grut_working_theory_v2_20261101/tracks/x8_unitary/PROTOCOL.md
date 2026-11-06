# X8 · Unitary mass projection and driven H55 decay

6 November 2026 · Parent 8a5be8c741a6d2d27aec9a4052e72b2f87b61003.

## Spatial source

Project the Hermitian scalar phase into the longitudinal band before
exponentiating it. Use the resulting unitary mass phase and its exact
field derivative in the same action. Keep only the physical fermion.
Use one common zero-external-momentum kinetic subtraction and the
finite-grid homogeneous metric matching from X7.

Compare harmonic response one through three against the independent
unbounded finite-periodic integral, including the frequency-dependent
quadratic response. Check cancellation of the unpaired contact term
with the direct finite Hamiltonian. Quantify the remaining cutoff
order and repeat a finite-amplitude static comparison. Test massless
decoupling, static vacuum, exact homogeneous reduction, independent
small spatial ODE, fixed-scale reversibility, scalar force from energy
variation and pressure from volume variation.

Use the X7 scalar wave and second-order graph preparation on the
expanding synthetic slab. Run to time 3 at step 0.0025. Compare
cylinders (9,12), (12,16), (16,22), and increase the cutoff further
if necessary. Independently vary longitudinal cutoff, transverse cutoff,
transverse quadrature, timestep and scalar collocation and bandwidth.
For a reference cutoff, use Gaussian radial panels, with 8 and 12
nodes per panel as the independent resolution check. Compare net
and absolute transfer, cumulative work and peak-scaled instantaneous
work to 1%. Check the reciprocal pressure-work budget at 1e-5 and
total momentum at 1e-9. Scope is one resolved spatial direction
with a homogeneous metric.

## Transport

Compute the order-Yukawa-squared decay kernel of a Majorana mode in
a phase rotating at constant speed into a vacuum of massless lepton
and Higgs daughters. Use the inherited widths and a Yukawa vertex
without a scalar phase in the original mass-phase frame. Derive both
frequency branches, calibrate to the undriven rest width and test
helicity reversal, the kinematic threshold and positivity.

Construct a substep retaining the full Bloch vector. Assign the
two branches their physical energies and derive the scalar momentum
impulse from the same collision ledger. Test fixed-a energy
conservation, the rotating-frame charge ledger and the zero-drive
free result. Replay the archived H55 state for 0.12 in H_i units at
the same two timesteps as X7. Compare with the collisionless and
free-axis width replays. Assess the constant-rotation approximation
along the live trajectory using instantaneous eigenbasis and
bath-threshold variation indicators. Keep the X7 simultaneous gate
of 1e-6 over a 0.02 window. Treat a locally constant vacuum kernel
as a sensitivity wherever its assumptions fail.

The physical source, spatial H55 initial data and the transport
overlap must all qualify before updating a relic abundance.

## Follow-up after the first source run

The original unitary vertex cancels the growing contact edge, while
the residual paired bubble and the all-mode kinetic subtraction
have different supports. The original run misses instantaneous
work convergence between the two largest cylinders and has a
pressure-work budget miss at the original wide timestep. Keep
these results. Subtract the average of the two zero-external
kinetic integrands on each momentum pair, using the same Fourier
coefficient in the scalar kinetic and gradient action. The
additional cutoff term is determined by the finite support and
goes to zero at fixed scalar momentum as the cutoff increases.

Repeat the static ladder and independent variations of cutoff,
radial integration, scalar bandwidth and quadrature. Run the
paired-action comparison through time 3 at timestep 0.00125
with a further half-step control. Use (9,12), (12,16) and
(16,22) and independent cutoff-axis changes. Preserve both
subtraction schemes in the evidence and apply the same
convergence and budget criteria.

## Two-time vacuum follow-up

The locally constant transport run misses the bath-threshold
timescale diagnostic on most of its predicted power. For a direct
comparison, evaluate the order-Yukawa-squared two-time vacuum
emission functional on selected exact archived H55 projectors.
Use the fine collisionless X7 scalar trace, two smooth windows of
duration 0.04 starting at offsets 0.02 and 0.065, and the 25th and
75th percentiles of each species' free energy distribution with
both helicities. Evolve those spinors on the prescribed history.
Compare the vacuum functional with the same window's local
frequency-branch prediction and refine the conformal timestep,
frequency quadrature and spectral cutoff. Include a stationary
positive branch and a stationary negative branch control. The
switching envelope defines this perturbative diagnostic. Coupled
memory backreaction and inverse or medium processes are separate.

## Startup follow-up

The paired-action run passes its energy and momentum budgets and
reduces the static tail, while its instantaneous cutoff comparison
still misses at early times. Repeat (12,16) and (16,22) with the
next term of the adiabatic projector in the ultraviolet. Apply
the third-order term with a smooth gap window from 3 to 5 scalar
masses and vary it to 4 to 6 scalar masses. Keep the low-gap
initial state, scalar data, and paired action fixed. Repeat the
widest comparison at half step if needed. This is a separate
preparation check. Preserve the order-two failures.

The first startup comparisons pass at each fixed gap window
and show sensitivity to changing the window. Complete the
independent longitudinal cutoff, transverse cutoff, radial rule,
collocation and scalar bandwidth checks with the original
3-to-5-scalar-mass window. Test preparation sensitivity separately
from the discretization gate.
