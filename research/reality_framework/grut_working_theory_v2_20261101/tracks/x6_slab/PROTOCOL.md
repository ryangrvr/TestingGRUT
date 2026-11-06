# X6 · Spatial production and transport gates

5 November 2026. Parent 7a4c5aa9caf534d806564e89dced4cd66dda28c9.
Read the [X5 handoff](../x5_coupled/HANDOFF.md) and keep its H55
interface failures as the comparator.

## Gate 1 · Spatial source

Build a periodic Galerkin slab. Retain one resolved spatial
coordinate and integrate the full four-component fermion
Hamiltonian over a transverse radial momentum measure. Evolve
all occupied negative-frequency orbitals and the scalar field.
Use the same mass vertex for local scalar force and fermion
work. Evolve the isotropic spatially averaged scale factor from
the common energy. Report the pressure tensor and the longitudinal
momentum. This metric reduction leaves shear and metric
perturbations outside the calculation.

Use static vacuum subtraction, the X5 second-order temporal and
expansion subtraction on the slab quadrature, and a spatial
kinetic subtraction matched to the static vacuum response at
one comoving harmonic. Include their derivatives in the same
homogeneous metric variation. Use a finite cylinder regulator.
A covariant continuum calculation and actual H55 geometry require
separate controls.

Use a resolved synthetic source with one scalar, one Majorana
species, box length $2\pi$, $M=0.7$, phase scale $f=1.15$,
$m_\phi=2.6$, $M_P=16$, and vacuum density 1.5. Start at $a=1$
with $\phi=0.15\sin x$, $\dot\phi=0.9\cos x$ and the first-order
adiabatic graph state of the initial Hamiltonian. Run through
$t=3$ in these control units. Use longitudinal integers through
6 and a transverse cutoff 8 with 16 radial bins. Use 49 scalar
quadrature points, scalar modes through 3 and a step of 0.005.
Refine scalar quadrature, scalar bandwidth, transverse bins,
cutoff and timestep independently. Increase the cutoff to
longitudinal integer 9 and transverse cutoff 12 at fixed
transverse spacing. The independent timestep test halves the step.

Recover the X5 homogeneous helicity equations on the identical
momentum quadrature. Preserve static-vacuum, massless, independent
matrix ODE, reversibility, energy and momentum checks. Compare
integrated comoving scalar work and the full instantaneous work
time series, scaled by its peak absolute value, at the 1% level.
Use a pressure-work residual threshold $10^{-5}$ of initial
comoving energy. Test the finite spatial subtraction against a
static vacuum perturbation. Preserve unsuccessful criteria.

## Gate 2 · Transport adequacy

Keep the covariance and its coherences until the transition is
quantified. Test the inherited H55 interface from the X5 saved
state. An active source needs a continuous overlap. A pressure or
work discontinuity rejects a particle-only handoff.

In a separate controlled slab preparation at fixed scale factor,
drive the mass phase over $0<t<3$ with
$\theta=\sin^4(\pi t/3)\cos(x-2.7t)$ and set it to zero outside
that support. Measure work from the declared drive and evolve
all occupied states through it. Use the reference slab momentum
quadrature and halve the pulse timestep. The tail has constant
free mass and a stated decay law into two massless daughters.
Its rest rate is 0.015. Test an explicit elastic band-dephasing
channel of rate 0.6 and a zero-dephasing control. The channel
acts only in the free tail and preserves occupations, energy and
momentum. These are specified test interactions.

Use the pulse energy as the scale for errors. Require free
source, no work discontinuity, energy and momentum agreement
within $10^{-6}$, and a phase-independent bound below $10^{-6}$
for every volume-averaged pressure component. Require more than
half the original pulse energy to remain in the parent gas at
an accepted boundary. Inspect the tail on a fixed grid through
$t=80$ measured from pulse end. Test the earliest passing
boundary and two later boundaries, offset by 2.5 and 5.
Compare daughter energy, momentum and pressure to the full
coherent overlap over that downstream interval. Keep local
spatial coherence distinct from volume-averaged transport.

A successful control shows when the reduction is adequate under
its specified dynamics. Applying it to H55 requires a spatially
resolved source and a justified collision and matching model.
The prior relic scan and its physical interpretation stay in X5.

[Quantum transport with coherence](https://arxiv.org/abs/0811.0936) ·
[Adiabatic stress and source subtraction](https://arxiv.org/abs/1703.00908)
