# X4 · Source audit and kinetic tests

## Archived-source replay

The actual H55 scalar lattice was rerun with passive diagnostics at every
scalar timestep. Both fermion moments and the Standard Model energy are
captured before the fluid updates. The trace separately accumulates energy
and number assigned where the unblocked perturbative floor wins. The
published scalar run and all previous tracks remain the comparators.

For the free-streaming counterfactual, each step and channel contributes a
cohort at its native mean injection energy. The cohort's comoving momentum
is conserved and its decay is time dilated using its own massive energy.
Two controls use the stationary-parent two-body energy on the perturbative
branch or remove that branch entirely, respectively. These are held-source
diagnostics. The calculation intentionally leaves the Pauli limit unenforced
there so that it can test the source's required moments.

## Test of the full kinetic ledger

The dephased M2 engine evolves four explicitly bounded occupation spectra.
Number, energy and pressure use the same momentum grid and two helicities.
The transition probability tables are read from the archived X3 kinetic
input. No calibration against a relic yield enters the calculation.

Use the following two preparations.

- At the final recorded actual H55 time, carry the aligned scalar, spectator,
  Standard Model and native species energies. Initialize each species with
  the X3 energy-matched half-occupied sphere. The initialized counts are
  recorded against the native source counts. This tests the energy ledger
  while keeping the deliberate spectral replacement explicit.
- For a separate reduced closure, initialize at the first post-crossing
  maximum of the box-mean waterfall field. Return the existing energy of
  N1 to its sigma parent and of N2 through N4 to their chi parent. Start
  the mode occupations at zero and evolve the harmonic parent energies
  and spectra jointly. This is a supplied preparation on an actual H55
  slice. It replaces the later inhomogeneous scalar evolution by two
  cycle-averaged pressureless parent condensates. It tests a realizable
  reduced model and cannot establish actual H55 production.

Both use the X3 hard/thermal split prescription, the 0.05 eV seesaw rate
choice and the inherited finite-width perturbative kernel. The spectator
is treated as a pressureless oscillation average. Use that approximation
from each test's start time and carry its decay to radiation in the same
Friedmann equation. At the first waterfall maximum the approximation has
less scale separation than at the archived interface.

Apply visible-sector transfers at fixed scale factor before redshift.
Integrate the spectator/decay system exactly during each fixed-H step.
Track the pressure-work budget and the fixed-volume transfer residual.
Double the momentum resolution and halve all three adaptive timestep
ceilings separately. Carry an omitted-spectator control to isolate its
effect, plus a reduced-model run that stops at the archived splice time.

Acceptance of the numerical ledger requires physical occupations,
a full energy gap below 1e-6, a fixed-volume exchange error below 1e-10
and a pressure-work budget within 1%. The late total and visible K
must each change by less than 1% under the two refinements. The prior
actual scalar source budget error and the changed or reset particle
numbers are independent barriers to scientific acceptance.

For a conditional dark-radiation estimate after decay and thermalization,
use the late density ratio and standard visible entropy conservation.
This estimate assumes a decoupled, relativistic daughter sector and a
high-temperature visible bath with g* = g*s = 106.75.

## Direct mode probe on the H55 mean waveform

The source tape stores the box-mean inflaton fields and their velocities
at every scalar step. In flat FLRW, rescale the Majorana field by
$a^{3/2}$ and use the phase $\theta=\bar\phi/f$ of its parent. After a
chiral rotation the homogeneous helicity block is

$$
\mathsf H_h=M\sigma_x+(h q/a+\dot\theta/2)\sigma_z.
$$

The helicity sign labels follow the chosen chiral rotation. Use cubic
Hermite interpolation through the stored values and slopes. Propagate
all four species and both helicities with a fourth-order product of two
unitary rotations per substep. The evolved state starts at the supplied
first-order adiabatic vacuum. Project it at each recorded output onto
both the instantaneous and first-order dressed vacuum axes. In Bloch
variables the latter is the normalized direction
$-\hat b+(\hat b\times\dot{\hat b})/(2|b|)$ for
$b=(M,0,hq/a+\dot\theta/2)$.

These finite-time projections define particle diagnostics. The calculation
uses free-particle moments of the projected occupation. It omits collisions,
spatial gradients and backreaction, and fixes the scalar waveform even
when the new production differs from its original drain. The projector
comparison measures basis dependence on this waveform. The renormalized
stress tensor requires a separate calculation.

Use a static zero-production control and compare a time-varying control
with a direct integrator in the unrotated basis. Check momentum spacing,
cutoff and timestep separately on the H55 waveform. Numerical refinements
test the solver on this supplied source. The lattice and homogeneous
projection remain distinct physical approximations.
