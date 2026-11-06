# X8 · Spatial source and the limit of a local decay width

6 November 2026 · [Equations and derivations](MODEL.md)

X8 eliminates the growing sharp-cutoff term using the physical
Majorana field alone. An expanding one-axis slab converges at fixed
ultraviolet preparation. Its instantaneous work depends on that
preparation, which is a separate physical input. The archived
homogeneous H55 state continues to fail the tested transport
handoff. The vacuum Yukawa calculation identifies a bath threshold
that prevents a uniform local-rate approximation along that history.

## Spatial source

The projected unitary phase preserves $V_Q^2=M^2 I$ and removes
the unpaired contact term. The finite spectral calculation checks
three harmonics at two scale factors and at static, real and
imaginary frequencies. It agrees with the paired kernel within
$5.73\times10^{-17}$ absolute. The homogeneous covariance agrees
with the X6 solver within $1.95\times10^{-14}$. A separate matrix
exponential and Fréchet derivative in a spatial ODE give
$2.52\times10^{-7}$ maximum error at the finest step. Reversal,
static vacuum, massless decoupling, field variation and pressure
variation also pass. [Controls](results/controls.json)

The paired subtraction matches the kinetic coefficient on the
same momentum pairs as the vacuum response. The quadratic
remainder approaches its independent unbounded periodic value
with measured cutoff exponents of 1.993 to 2.002. Refining the
unbounded longitudinal sum removes its influence on the largest
cutoff comparison. The maximum error is $1.01\times10^{-6}$ relative.
At finite phase amplitude, the last two static paired values
differ by 0.068%. The original single-mass dynamic test retains
its energy-budget and instantaneous-work misses.

| Test | Instantaneous work difference | Cumulative work difference | Terminal transfer difference |
|---|---:|---:|---:|
| All-mode subtraction, (12,16) to (16,22) | 2.234% | 0.199% | 0.0546% |
| Paired subtraction, order-two state, same cutoff change | 2.053% | 0.0711% | 0.0189% |
| Paired subtraction, order-three UV state, same cutoff change | 0.273% | 0.0632% | 0.0206% |
| Same cutoff change with the later UV gap window | 0.273% | 0.0625% | 0.0201% |
| Order-three longitudinal cutoff alone | 0.0670% | 0.0137% | 0.00331% |
| Order-three transverse cutoff alone | 0.171% | 0.0370% | 0.0126% |
| Order-three radial quadrature | 0.0112% | 0.000799% | 0.000708% |
| Order-three wide timestep halved | 0.236% | 0.0539% | 0.0205% |
| Change the UV gap window at fixed (16,22) | 4.734% | 0.545% | 0.330% |

Instantaneous differences are scaled to the peak of the reference
work, and cumulative differences to its peak cumulative transfer.
The third-order term is applied above a smooth gap range of three
to five scalar masses, with a separate four-to-six-mass test.
Both choices retain the order-two low-gap state. The cutoff
result concerns a fixed choice. The difference between choices
quantifies state sensitivity even when the numerics converge.
The collocation and scalar-bandwidth controls also meet the 1%
work criterion. The largest order-three pressure-work and
comoving momentum residuals are $6.67\times10^{-6}$ and
$3.59\times10^{-11}$.
[All-mode source](results/source.json) ·
[Paired source](results/paired_source.json) ·
[Preparation tests](results/preparation.json)

This control isolates one spatial direction with an averaged
metric and a synthetic scalar wave. Applying it to H55 needs
the physical scalar and fermion spatial state and matching
for the corresponding geometry.

## Physical vacuum decay and the H55 overlap

At constant phase speed the two massless lepton–Higgs branches
have laboratory frequencies $\Omega\mp\dot\theta/2$. For
helicity opposite the drive one branch closes above
$p=M^2/(2|\dot\theta|)$. The independent angular phase-space
integration reproduces the rates to $3.30\times10^{-17}$
absolute in the benchmark units. At zero drive the total
width is the usual time-dilated width.

X8 evolves the exact archived four-species Bloch state for
another $0.12H_i^{-1}$ with this locally constant, flavor-diagonal
vacuum kernel. The energy ledger includes the scalar momentum
impulse and the laboratory energy delivered to the bath.
Both timestep runs conserve total energy with a pressure-work
residual below $5.14\times10^{-12}$ and obey the Pauli bounds.
The finer run transfers 0.1303% of its terminal total density
into visible radiation in this sensitivity calculation.
There are zero passing samples at either resolution for the
simultaneous energy, pressure and work criteria. The most
favorable $0.02H_i^{-1}$ window has a maximum normalized gap
of 0.882 on the fine run, against $10^{-6}$.

The rotating eigenbasis and secular diagnostics cover 93.3%
of the locally predicted integrated decay power. Once the
bath-threshold timescale check is included, only 21.45%
passes all the diagnostic bounds. The resulting radiation
fraction belongs to the locally constant sensitivity.
[Transport results](results/transport.json)

## Two-time vacuum control

The leading two-time emission functional was evaluated for
16 exact archived modes under the prescribed collisionless
H55 scalar trace. Each mode was tested in two smooth finite
exposures. The stationary control reproduces the integrated
width to $5.30\times10^{-11}$ relative. In the H55 windows the
two-time emitted energy is 2.48 to 27.76 times the local
estimate. Independent conformal timestep and spectral spacing
refinements change it by at most $4.42\times10^{-5}$ relative.
Halving the timestep of the prescribed background changes it
by up to 0.105%.

The same finite exposures on a frozen background give
2.24 to 32.17 times their own local estimate. The evolving
and frozen two-time energies differ by at most 10% for the
selected samples. Changing the exposure envelope also shifts
the ratios. This identifies the finite exposure as the main
source of the large discrepancy in this control. It leaves
the always-coupled H55 radiation yield open. That calculation
requires a renormalized two-time self-energy, the interaction
energy and the initial sterile–bath correlations.
[Two-time results](results/memory.json)

The archived 50.96% fermion energy is the interacting energy
fraction at the old slice. The inherited 1.28% terminal N4
refinement miss and zero of 135 strict acceptances remain part
of the evidence. The model supplies a conditional source and
an operational test for the emergence of particle transport.
A GRUT-specific empirical prediction requires a derived state
or coupling relation with a discriminating observable.
[Working theory](../../WORKING_THEORY.md) · [Verification](results/verification.json)
