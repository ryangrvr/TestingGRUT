# X5 · Coupled production in the GRUT matter candidate

5 November 2026 · HFa at $H_i=5.5\times10^{10}$ GeV

X5 implements the missing two-way response in a finite homogeneous
version of the H55 candidate. The evolving scalar phase produces
fermions. The same quantum state supplies the opposing force on
the scalar and its contribution to expansion. Particle counts and
energies follow from that state. In the tested homogeneous model the
fermion sector holds about 51% of the energy at the archived interface.
The archived 3D scalar field, however, has order-one spatial variation
in the fermion mass operator around the waterfall. The calculation
defines a working coupled reduction and identifies the physical step
needed to make the H55 history predictive.

## Coupled dynamics

The [model](MODEL.md) derives a finite-regulator evolution for the
inflaton, waterfall field, spectator and both helicities of the four
Majorana species. For each mode a unitary midpoint rotation changes
its axial component by $\Delta r_z$. The scalar receives the
canonical impulse $w\,\Delta r_z/(4f)$. The same rotation transfers
the opposite energy to the fermion. A discrete gradient implements
the scalar potential change. The spectator's mass coupling acts
reciprocally on the waterfall field, and its damping deposits its
energy into dark radiation.

The static vacuum is subtracted. A common second-order adiabatic
subtraction modifies the scalar inertia, pressure and gravitational
constraint. The finite scheme uses a comoving regulator and the full
second-order subtraction. A covariant stress calculation with
independent finite matching remains open. The main run is
collisionless. A vacuum amplitude-damping control tests sensitivity
to decay and carries its energy into a visible component.

Eleven H55 runs completed. The largest pressure-work residual is
$7.1\times10^{-10}$ of the initial comoving energy. The fully coupled
test agrees with a separate unrotated ODE to $1.14\times10^{-5}$
relative in the flat control and $6.27\times10^{-5}$ in the expanding
and decaying control. Both discrepancies decrease by about four per
step halving. The massless result coincides with the empty-species
control. The finite shift charge is conserved to $6.7\times10^{-16}$
in the flat symmetry control.

For a constant rotating vacuum, increasing the cutoff from 16 to
128 raises the bare energy from 0.0059484 to 0.0105054. The two
largest cutoffs give subtracted energies 0.000181758 and 0.000181877.
The independent quadrature tests the normalization and the leading
cutoff subtraction on this controlled state.

## Homogeneous H55 outcome

These runs start from the archived H55 mean fields and velocities
and the declared dressed initial vacuum. The table reports the
old diagnostic time $t=2.9980937H_i^{-1}$.

| History | Radial bins | Comoving cutoff | Fermion energy fraction | Scale factor |
|---|---:|---:|---:|---:|
| Scalar control | 720 | 45,000 | 0 | 3.70158 |
| Coupled reference | 720 | 45,000 | 0.49900 | 3.23430 |
| Coupled fourfold radial refinement | 2,880 | 45,000 | 0.51082 | 3.23068 |
| Coupled eightfold radial refinement | 5,760 | 45,000 | 0.50956 | 3.23108 |
| Wider cutoff at matched fine spacing | 4,320 | 67,500 | 0.51120 | 3.23060 |
| Optional vacuum decay control | 720 | 45,000 | 0.48178 | 3.23391 |

The fractions use the normal-ordered physical mass energy before
allocating the second-order subtraction to the fermions. Allocating
it to the fermions gives 0.50867 at 5,760 bins. The optional decay
control also puts 2.23% of total energy into visible radiation.

The finest N1 energy is 58.8% of its legacy value. All four species
have an admissible free-particle count and energy at the slice.
Their energies are 1.30–1.40 times their respective Pauli minima.
The model therefore supplies a physical occupation state in this
reduction without prescribing a separate fluid drain or replacing
the counts after evolution.

| Refinement | Largest terminal moment change | Largest final-window moment change |
|---|---:|---:|
| 720 to 1,440 bins | 5.64% | 4.92% |
| 1,440 to 2,880 bins | 3.12% | 3.26% |
| 2,880 to 5,760 bins | 1.28% | 0.980% |
| 45,000 to 67,500 cutoff at fine matched spacing | 0.146% | 0.133% |
| Half timestep at 2,880 bins | 0.0770% | 0.0938% |

Each row checks the energy and count of every species in both
projectors. The final window is the last $0.10 H_i^{-1}$. The
time-averaged moments meet the 1% criterion at the fine grids.
N4's terminal free count gives the remaining 1.28% miss. The
instantaneous net work is more sensitive to phase than these
moments. Its sign changes between the two finest grids.

## Why the spatial calculation matters

On the passive archived 64³ scalar replay, define the fractional RMS
difference of the fermion mass operator from its value on the mean
parent field as

$$
D_\phi =
\left\langle\left|e^{i\phi/f}-e^{i\langle\phi\rangle/f}\right|^2
\right\rangle^{1/2}.
$$

The replay matches the archived expansion history to $5.0\times10^{-10}$
at the compared times. The measured variation peaks around the
waterfall even when the field is highly coherent by the old scalar
energy diagnostic.

| Parent | Peak $D_\phi$ | Time in $H_i^{-1}$ | Scalar gradient energy fraction at that time |
|---|---:|---:|---:|
| Inflaton $\sigma$ | 1.4756 | 0.24775 | $1.75\times10^{-7}$ |
| Waterfall $\chi$ | 1.4956 | 0.25563 | $5.79\times10^{-7}$ |

The scalar coherent fraction at these times is 99.87% and 99.96%,
respectively. Small scalar gradient energy does not control the
dimensionless phase error of a mean-field fermion Hamiltonian.
This diagnoses the old scalar background. The response of the
fermions on a new spatially coupled H55 history remains to be
calculated.

X5 also evolves a flat periodic one-dimensional Dirac control with
all initially occupied negative-energy orbitals and reciprocal local
scalar work. A scalar wave with vanishing spatial mean yields
0.41449 units of fermion excitation energy, 18.32% of the initial
scalar energy in that example. The run initialized at the spatial
mean stays at vacuum to $1.4\times10^{-11}$ in energy. Doubling
spatial resolution changes the wave's terminal energy by 0.0244%.
Halving the timestep changes it by 0.00057%. A direct finite-matrix
ODE agrees to $3.9\times10^{-8}$ in state components, and reversing
the split evolution restores the initial state to $2.3\times10^{-13}$.

The control uses the full uniform vacuum Hamiltonian in the free
substep and integrates the residual mass force locally. This
preserves the static vacuum. The first kinetic-plus-local-mass
split produced a $2.73\times10^{-7}$ static excitation artifact.
Halving its timestep reduced that to $1.71\times10^{-8}$ before
the vacuum-preserving split replaced it. The separate rotating-vacuum
quadrature initially missed its tolerance at the widest cutoff,
with a $2.35\times10^{-5}$ discrepancy on 4,096 bins. Doubling
that control grid reduced the discrepancy to $5.88\times10^{-6}$.

## Coherent-to-kinetic interface

The direct free-mass projection at the interface gives Pauli-feasible
occupations. Dephasing along that same mass axis preserves count
and free energy but drops the off-diagonal terms. At 5,760 bins it
changes the pressure by $-7.57\times10^{-5}$ of the total density
and removes the instantaneous scalar work carried by the state.
The value of that work is phase-sensitive across grids.

Replacing the subtracted stress energy with the dressed particle
energy leaves a fractional gap $9.76\times10^{-4}$ at the finest
slice. This misses the inherited $10^{-6}$ energy criterion. The
phase drive $|\dot\theta_j|/(2M_j)$ ranges from 3.65 to 22.63 at
that time. The held interface is still inside active coherent
evolution.

The next physical calculation needs spatially resolved fermions
coupled back to the scalar fields, and a transition that evolves
the residual work and stress into the kinetic description. The
[handoff](HANDOFF.md) gives a slab construction that starts from
the verified one-dimensional work step and restores transverse
fermion momentum. The strict scan remains 0/135. No release,
accepted abundance point or distinctive observable was calculated
in X5.

All 42 X5 evidence checks and 39 prior X4 evidence checks pass.
These verify the calculations and retain the failed physical and
terminal-moment criteria. The comparisons use independent numerical
methods implemented within this investigation.

[Protocol](PROTOCOL.md) · [Refinements](REFINEMENT.md) · [Replay](REPLAY.md) · [Evidence](results/verification.json)
