# X6 · Spatial production and particle transport

5 November 2026 · Continuation of the H55 matter candidate

The four-component slab recovers the homogeneous helicity equations
and balances scalar work, fermion stress, and expansion. Its
instantaneous work still depends on the momentum cutoff at the
registered accuracy. The transport calculation finds a certified
volume-averaged boundary in a separate control with specified
interactions. H55's archived boundary fails the pressure and work
test even when its free-mass occupations preserve energy and count.

## Gate 1 · The spatial source

The periodic slab evolves the occupied fermion orbitals and one
spatial scalar field with full transverse momentum. The metric is
an isotropic spatial average. The local scalar force is the
variation of the same mass matrix that appears in the fermion
Hamiltonian. The expanding runs include the scalar gradient,
anisotropic fermion pressure and the scalar and metric
subtractions of the finite model. See [MODEL.md](MODEL.md).

The homogeneous limit matches an independent helicity ODE on the
same momentum grid to $4.14\times10^{-8}$ in occupation and
$3.28\times10^{-6}$ relative in fermion energy. The spatial
matrix ODE agrees to $7.63\times10^{-8}$ in state components.
Reversing the fixed-scale split restores the state to
$7.1\times10^{-14}$. The static-vacuum and massless controls
pass. Finite differences of the same energy recover the scalar
force to $3.5\times10^{-15}$ absolute and isotropic pressure to
$1.4\times10^{-10}$ relative in the action check.

The nine first-order source runs and four follow-ups satisfy
the combined pressure-work tolerance of $10^{-5}$ of their
initial comoving energy. The reference residual is
$1.71\times10^{-6}$. Its separate fermion work ledger closes
within $1.34\times10^{-6}$ after reallocating the known
kinetic subtractions. The longitudinal momentum ledger and
occupied-state orthonormality also pass.

| Independent change | Matched accumulated work at the end | Matched instantaneous work, peak scaled |
|---|---:|---:|
| Scalar quadrature 49 to 97 | Below $10^{-9}$% | Below $10^{-8}$% |
| Scalar bandwidth 3 to 5 | Below $10^{-9}$% | Below $10^{-8}$% |
| Transverse bins 16 to 32 | 0.945% | 1.23% |
| Transverse bins 32 to 64 | 0.233% | 0.294% |
| Halve reference timestep | 0.0069% | 0.0722% |
| Halve wide-cylinder timestep | 0.0311% | 0.376% |
| Increase longitudinal cutoff alone | 0.0218% | 2.62% |
| Increase transverse cutoff alone | 0.130% | 3.51% |
| Cutoff $(6,8)$ to $(9,12)$ | 0.148% | 6.41% |
| Cutoff $(9,12)$ to $(12,16)$, fine step | 0.0210% | 3.91% |
| Cutoff $(6,8)$ to $(9,12)$, second-order preparation | 0.122% | 4.87% |

Cutoff entries give the longitudinal integer and transverse cutoff.
The peak-scaled column is the maximum time-dependent difference
divided by the reference peak, including the oscillating phase.
The original bare wide-cutoff comparison is retained. Its
instantaneous difference is 7.62%. The matched comparison's
0.148% terminal integral difference accompanies a 0.412% maximum
difference over the accumulated-work history. The wider
cutoff and the second-order initial state leave an instantaneous
error above the 1% criterion. The spatial source gate is open
at that criterion.

### Why the static matching depends on the regulator

The static vacuum curvature has an interior contribution from
virtual neighboring momenta and an edge contribution from
neighbors removed by the Galerkin projector. The decomposition
in [MODEL.md](MODEL.md) is an exact identity for this regulator.
At the first harmonic the edge contributes 79.7%, 83.1% and
85.5% of the coefficient on the three evolved cylinders.
Independent static vacuum finite differences verify the
response at harmonics one through three on the two lower
cutoffs.

At fixed aspect ratio the edge contribution to the curvature
grows with the cutoff and scales as $|k|$. A single local
gradient term scales as $k^2$. At the largest analytic
cylinder the edge coefficients after division by $k^2$ have
ratios 1, 0.502 and 0.335. Matching one local coefficient at the
first harmonic leaves a regulator term at the others. This
identifies a problem with extending the present sharp-cylinder
matching to a spatial continuum. The measured work
differences retain contributions from both dynamics and
preparation.

## Gate 2 · Transport boundary

A separate traveling phase pulse produces the initial state for
a fixed-volume control. Its energy is 0.145756 and its
longitudinal momentum is 0.0348687 in control units. The
independent integral of the supplied drive matches that energy
and momentum to $1.17\times10^{-7}$ and $3.15\times10^{-8}$,
respectively, relative to the pulse energy. The drive and its
derivatives vanish at pulse end.

In the free tail, evolve the full covariance through a specified
vacuum decay and optional elastic band dephasing. Count,
energy and momentum stay continuous across the dephasing map.
The decay daughters receive the energy and momentum of the lost
parent occupations and their stress comes from the boosted
two-body decay. Their mean pressure is one-third of their
energy density. The elastic control uses a rate of 0.6 in these units
and the rest decay rate is 0.015. They are test interactions.

| Boundary time after pulse | Parent energy remaining | Pressure envelope / pulse energy | Maximum downstream pressure error / pulse energy |
|---|---:|---:|---:|
| 0 | 100% | $7.68\times10^{-2}$ | $6.09\times10^{-3}$ |
| 18.65 | 92.02% | $9.73\times10^{-7}$ | $3.00\times10^{-7}$ |
| 21.15 | 91.01% | $2.15\times10^{-7}$ | $3.34\times10^{-8}$ |
| 23.65 | 90.00% | $4.73\times10^{-8}$ | $1.16\times10^{-8}$ |

At the first row the exact particle projection preserves total
energy and momentum while the pressure criterion rejects the
cut. The other three boundaries pass the registered $10^{-6}$
criterion. The independently initialized continuations agree
on the downstream daughter and parent moments to
$3.9\times10^{-16}$ relative in this linear free-tail model.
The pressure comparison uses the full coherent reference,
and its envelope is independent of phase cancellation.

Without elastic dephasing, the pressure envelope remains 3.73%
of pulse energy at the parent half-energy time of 162.17.
The envelope decreases monotonically, so the control with
vacuum decay alone has no boundary that passes before half
the parent energy is lost. Same-band spatial covariance also
persists with the test dephasing. These results certify
volume-averaged transport in the stated control and leave
local spatial transport outside that certification.

### H55 interface

The H55 check uses the archived 5,760-bin X5 state. The dressed
energy gap is $9.76\times10^{-4}$ of total density. Free-axis
dephasing preserves the free energy and count and changes
pressure by $-7.57\times10^{-5}$ of total density. It also
removes the work carried by the coherent state. Even for an
energy-matched cut, the comoving fermion-energy derivative
would jump by

$$
\frac{\dot E_{F,\rm gas}-\dot E_{F,\rm coherent}}{a^3\rho}
=-\frac{\langle J\dot\phi\rangle}{\rho}
-3H\frac{P_{\rm gas}-P_{\rm coherent}}{\rho}
=-0.10831 H_i .
$$

This is the derivative mismatch on the saved slice. X5
already shows the instantaneous power is phase-sensitive
under refinement. The terminal N4 count miss is preserved.
The test interaction was confined to the controlled pulse.
H55 needs an evolved spatial source and an overlap using
its own physical collision and stress matching.

## Result for the working theory

X6 supplies a dynamical criterion for when particle occupations
can carry the relevant averaged observables. It also shows
why an occupation count and exact energy matching can be
insufficient at an interacting boundary. These results concern
the matter candidate and its standard-physics dynamics. The
strict relic decision remains zero of 135. The next source
and interface targets are in [HANDOFF.md](HANDOFF.md).

All 48 archived evidence checks and five fresh-replay comparisons
pass. These checks record the physical cutoff and H55 handoff
failures as well as the successful controls.

[Protocol](PROTOCOL.md) · [Refinements](REFINEMENT.md) · [Replay](REPLAY.md) · [Evidence](results/verification.json)
