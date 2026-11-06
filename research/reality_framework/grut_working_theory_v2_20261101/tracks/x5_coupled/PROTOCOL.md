# X5 · Coupled H55 production test

5 November 2026. Parent 83ce5341a3819585d250f6a6bd6f3da6390496f0.
The [X4 handoff](../x4_closure/HANDOFF.md) supplies the immediate scientific
state. The strict relic scan is 0/135.

Build a finite-comoving-cutoff, homogeneous semiclassical HFa model with
the same H55 couplings and initial mean fields as the archived scalar
run. Evolve the four Majorana species and both helicities. Use their
state to compute scalar work, stress and occupation moments. Evolve
the spectator as a field, apply its reciprocal force to the waterfall
field, and put the energy lost to spectator damping into dark radiation.
Use the X4 spectator width and initial spectator amplitude. The primary
fermion calculation is collisionless. A separate control tests vacuum
decay with the inherited widths and an instantaneous mass projector.
That control adopts an explicit approximate collision operator during
production.

Subtract the static vacuum and the second-order adiabatic polarization
in the homogeneous action. Derive the scalar and metric terms from the
same subtraction and specify the finite matching condition. Treat a
comoving regulator and this subtraction as a reduced model. A fully
covariant curved-space renormalization requires additional matching.
Use a first-order adiabatic initial state, following the X4 mode probe.
For comparison evolve the same scalar model with the fermion coupling
turned off.

Verify the homogeneous mode evolution against an independent integrator
in the unrotated frame on a time-dependent coupled test. Verify static
vacuum, the zero-mass limit, energy exchange, the pressure-work budget
and admissible fermion occupations. For the main H55 run halve the
timestep, double the radial resolution and increase the cutoff at fixed
spacing separately. Check each species' energy and count both at the
archived interface time and averaged over its final 0.10 inverse Hubble
units. A 1% relative change is the numerical target. A finite model
must have a maximum pressure-work residual below 0.01% of its initial
comoving energy. Preserve failures and any further resolution controls.

Compare the archived spatial phase variation to the scale used to
justify taking a homogeneous mean. Report stress and particle moments
separately and check the energy cost of discarding polarization at the
interface. A reduced history can support a conditional calculation
within its tested domain. Acceptance of the actual H55 source also
requires control of spatial production, collisions and renormalized
matching. No new abundance root is registered here.

The action and tests draw on
[fermion production with pseudoscalar backreaction](https://arxiv.org/abs/1803.04501)
and the consistency of
[adiabatic subtraction in coupled semiclassical equations](https://arxiv.org/abs/1703.00908).
The X5 regulator, Majorana normalization and discrete work law are
specified separately.
