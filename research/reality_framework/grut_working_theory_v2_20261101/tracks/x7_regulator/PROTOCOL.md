# X7 · Spatial cutoff completion and H55 overlap

5 November 2026. Parent 08b7e792c89dded6242a521c3cce581693a07e67.
[X6](../x6_slab/HANDOFF.md) establishes the sharp-cylinder edge
artifact and rejects the archived H55 transport boundary.

## Regulator test

Keep the X6 finite spatial Hamiltonian and its reciprocal scalar work.
Evaluate the physical response plus two signed auxiliary mass responses
with the same local phase. Let the coefficients satisfy
$\sum_j c_j M_j^2=\sum_j c_j M_j^4=0$, with physical coefficient one.
The auxiliary sectors are vacuum-polarization subtractions. Report
physical particle observables from the physical sector.

Test the static response at harmonics one, two and three against
an independent finite-periodic longitudinal calculation with
unbounded transverse momentum and frequency. Distinguish the
interior dispersion from the cutoff edge. Match the spatial and
temporal local kinetic coefficients together at zero derivative
order in the response coefficient. The finite-cylinder calculation
uses its discrete zero-external-momentum kinetic matching and its
second-order metric coefficient. Report residual differences against
the unbounded finite-volume result. Evaluate the limit in which
the physical mass vanishes and check regulator-mass sensitivity.

Use the X6 expanding slab control and its first scalar data at
$a=1$. Prepare second-order adiabatic graph states on a common
initial scalar and metric tangent. Auxiliary masses are 3.2 and
4.8, with a 4.8 and 7.2 sensitivity test. Use cylinders $(9,12)$,
$(12,16)$ and $(16,22)$, with the first entry the longitudinal
cutoff integer and the second the transverse cutoff. Transverse
spacing is 0.5 in the cutoff comparison. Run through $t=3$ at
step 0.0025. Halve the reference step and double its transverse
bins independently. Check both accumulated matched work and
peak-scaled instantaneous work at 1%, the pressure-work budget at
$10^{-5}$ and momentum conservation. Preserve the X6
benchmarks and any failures.

Use a smooth, compact pulse of duration 8 with phase
$0.22\sin^4(\pi t/8)\cos(x-2.6t)$ to check auxiliary excitations
after the drive has ended. The target is below 1% of the physical
excitation energy in the sum of the absolute weighted auxiliary
excitation energies. Distinguish this control from any claim
about heavy-field decoupling on actual H55.

## Coherent H55 overlap

Independently restart the exact archived X5 homogeneous H55
orbitals and fields. Evolve for $0.12 H_i^{-1}$ past its stored
boundary without discarding coherence. Use two timesteps on
the same saved grid. Compare collisionless evolution and
vacuum decay at the inherited widths, using the stated X5
instantaneous-basis decay prescription only as a sensitivity
control. Preserve the initial occupation state and energy.
This is a boundary diagnostic on the archived homogeneous
reduction.

Evaluate the energy, pressure and residual work tolerances at
each candidate time. Require simultaneous gaps below $10^{-6}$
in the corresponding normalized quantities across a
$0.02 H_i^{-1}$ window. Also compute the most favorable
phase-independent pressure envelope with the initial free
basis frozen and the drive removed, using the actual vacuum
decay rates. Treat that limit separately from the live overlap.
Keep the inherited spatial and terminal-refinement decisions
and carry no transport result into a new abundance calculation
unless the physical source and its collision overlap qualify.

[Derivative-coupling regularization context](https://arxiv.org/abs/2112.07645) ·
[Quantum coherence in transport](https://arxiv.org/abs/0811.0936)

## Radial follow-up after the first ten runs

The original uniform 0.5 to 0.25 transverse-bin comparison
missed the instantaneous-work target, and the heavier mass
pair changed the matched work history. Preserve both misses.
Hold the original fields, couplings, timestep and criteria
fixed. At cylinder (12,16), halve the transverse width again
to 0.125. Independently integrate the radial measure with
positive Gauss–Legendre weights in panels with edges
0, 1, 2, 4, 8, 12, 16 and the wider cutoff at 22. Compare
8 and 12 nodes per panel, test the wider cylinder at each
order, compare fine uniform and Gaussian rules, and repeat
the heavier-pair test at (16,22) with the 8-node rule.
Record the parameter sensitivity separately from discretization.
