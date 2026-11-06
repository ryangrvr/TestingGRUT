# Continue GRUT from X5

5 November 2026 · TestingGRUT branch codex/grut-h55-coupled-production-20261105.

Read [REPORT.md](REPORT.md), [MODEL.md](MODEL.md) and [REPLAY.md](REPLAY.md).
The [X4 handoff](../x4_closure/HANDOFF.md) links the earlier whole-theory
synthesis and records why the legacy fermion moments fail. The strict
relic scan is 0/135.

X5 builds a finite-regulator, self-consistent homogeneous H55 source.
It replaces the prescribed fermion drain with an evolved Majorana
state and reciprocal scalar force, includes the spectator and its dark
radiation, and derives the kinetic and constraint subtractions from
one matched homogeneous action. At the archived slice the finest
evolution has 50.96% of the energy in the normal-ordered fermion
sector and $a=3.23108$. The energy and count moments averaged over
the last $0.10H_i^{-1}$ agree within 0.981% on the tested fine controls.
The terminal N4 free count misses the 1% gate at 1.28%. The saved
spectral state has both helicities for each of the four species.
See [replay and state format](REPLAY.md).

The actual old H55 scalar lattice has RMS mass-phase errors of
about 1.48 and 1.50 for the two parent fields near the waterfall,
despite a gradient energy fraction below $10^{-6}$. A separate
fully occupied 1D spatial control proves that a scalar wave with
zero mean can produce fermions in a calculation with reciprocal
work. Neither the spatial control nor the homogeneous H55
evolution determines actual H55 production in three dimensions.

At the nominal interface the dressed particle energy differs
from the subtracted stress by $9.76\times10^{-4}$ of total density.
Free-axis dephasing preserves its energy and count while changing
the pressure and removing ongoing scalar work. The boundary
requires a dynamical transition with a checked common stress.

## Next implementation target

Build an intermediate slab with a spatially varying scalar phase
along one coordinate and transverse translation symmetry. For
each transverse momentum $\mathbf q_\perp$, evolve a four-component
block with kinetic operator
$\boldsymbol\sigma\cdot(-i\partial_x,\mathbf q_\perp)/a$ in the two
opposite chiral sectors and local off-diagonal mass
$M e^{i\theta(x,t)}I_2$. Integrate its occupied state over the
transverse measure. The homogeneous limit must recover the
helicity equations. The local scalar impulse must come from the
same mass-block variation as the fermion work. The spatial
subtractions require gradients of the mass phase as well as its
time derivative.

The 1D implementation demonstrates the discrete template. Split off
the full uniform-vacuum Hamiltonian and integrate the residual local
mass rotation and its scalar force together. The scalar gradient
evolution and the Fourier kinetic propagation can be controlled
independently. A later stochastic trace estimator can reduce cost
by estimating force and stress from the same propagated projector,
with separate noise and cutoff controls.

Only a resolved spatial source and a controlled conversion of
coherence into kinetic transport can justify a new release or
relic calculation. The fine homogeneous state is available as
a benchmark for those calculations. Scientific acceptance
remains false.
