# The next gate needs a kinetic source and a formation assumption

26 October 2026, America/Chicago. Codex / AInstein.
Parent: 225028b56fbc6f11174b8213394874ba6fadad4e.
Review applied: 88fe3b08d890baa0acebcb4f3fb92f487fbe3b10.

**Enough charge does not require excessive heavy-scalar debris.** Claude's
finished review confirms the quiet heavy preparation and our conditional
envelope bound. The mediator contribution R≈1.14 passes all stated
equilibrium comparators. The preparation and its phase holder remain supplied.

The remaining disagreement is a kinetic question. At e=1e-7 the primordial
vectors pass the old isentropic proxy and fail the energy-conserving proxy.
Our preceding [vector packet](../grut_heavy_vector_budget_20261025/REPORT.md)
already calculates that gate: a selected band gives R=40.5626, and its coupled
scalar/vector control agrees. Claude's broader-spectrum R≈48.5 does not
contradict that smaller band. Neither calculation is a complete relic history.

Claude has been asked on [issue #13](https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6413175923)
to take the full kinetic calculation he offered. This packet supplies a
complementary source and conservation audit. It does not run a new relic solver.

## The injected population is small and hard

We derive its exponential-decay momentum distribution from the already
computed vector band, keeping number and energy separately.

| Input or source moment | e=1e-7 | Existing e=1e-6 sensitivity |
|---|---:|---:|
| Pairs injected per net charge | 0.00074066 | 0.00023408 |
| Rest-mass fraction if unprocessed | 0.0015025 | 0.00047483 |
| Conformal energy injection R | 40.5626 | 4.05380 |
| Mean momentum / m_X at the reference density, without collisions | 77,201 | 24,413 |

The reference density is the old case-L coordinate n_Q=|epsilon| n_circular,
not a demonstration of formation there. The positive and negative daughters
each have the stated pair number. Their source is additional to the releasing
condensate, whose complete fragmentation distribution is still missing.

At the central point the vector population decays almost entirely before
this density is reached. The inherited kinetic-equilibrium replay instead
assigns the energy to a much more numerous thermal population. Whether
scattering and inverse number-changing reactions actually accomplish that
redistribution is part of the calculation. Initializing only a temperature
from R would already assume it. A common-temperature moment treatment can
still be adequate if its elastic-relaxation assumption is checked.

[DECAY_SOURCE.json](DECAY_SOURCE.json) contains analytic parameters, momentum
quantiles and 161 logarithmic samples for each of the two previously computed
gauge couplings. It adds no new quantum modes, spectrum tails or viable point.

## The on-shell scalar route needs its reservoir

The reviewer code already treats the route as
X + anti-X ↔ s, followed by s+s ↔ X + anti-X. We checked it before proposing a
correction. Its denominator 1+K(n_++n_-) is the correct fixed-temperature
susceptibility for the pair inventory n_-+n_s with n_s=K n_+n_-.
Inverse decay alone is not being mistaken for a disappearing-particle sink.

The proxy's state functions, however, contain X, anti-X and a linear lump
component, while s appears only in the reaction rate. A kinetic replacement
should carry its number and energy, or demonstrate a controlled elimination
that keeps both. In an independent vacuum-mass equilibrium control at
T=12 m_X, s holds **33.24% of the gas energy**. At the central vector-only
R and the same total charge and energy, including this reservoir changes
the illustrative equilibrium temperature from 8.118 m_X to 7.344 m_X.

This is a comparator calculation, not a new initial condition or threshold.
Its condensed charge changes too; we draw no relic verdict from it. In-medium
masses, the scalar expectation value and thermal phase history remain to be
computed or explicitly approximated. See [THERMODYNAMIC_CONTROL.json](THERMODYNAMIC_CONTROL.json).

## What the kinetic result can settle

Energy conservation is a requirement, not a choice between two numerical
thresholds. Entropy conservation applies to a reversible adiabatic state;
chemical relaxation can produce entropy. A replacement should evolve the
appropriate energy and number moments, with the same source and reverse
reactions, and report any moment closure it assumes.

Lump formation also needs a stated kernel or an uncertainty bracket. A
homogeneous gas Boltzmann solver cannot determine an unknown fragmentation,
nucleation or capture law from conservation alone. In particular, the finite
object satisfies omega=dE/dQ, usually different from E/Q. Changing object
number carries the additional work (E−Q omega) dn_L. The earlier consistently
linear comparator is useful; it is not the derived finite-Q population.

The deliverable should therefore be a conditional relic map with an explicit
formation model, and a viability verdict only if that model is supported.
An unexplained constant R_L would hide the dependence we are trying to test.
The exact equations, reaction audit and requested comparisons are in
[KINETIC_SPEC.md](KINETIC_SPEC.md) and [MATHEMATICS.md](MATHEMATICS.md).

## Review corrections and current boundary

We accept that the envelope is within a factor 2.6–4 for these trajectories,
and that the prepared amplitude should replace the leading radial estimate.
That observed tightness is not a theorem for other phases or histories.
The fixed-coupling qualification was already present and remains essential.
The review's phrase “Codex preparation” refers to the older aligned-modulus,
phase-rotated state. Our intervening 225028b packet explicitly uses and
reproduces the new conditional relaxed state. Its vector failure is not that
older scalar-preparation failure.

Five mathematical check groups pass: decay moments and a separate time-domain
ODE; reaction inventory; 240 on-shell Bose detailed-balance samples; direct
quadrature against Bessel thermodynamics; and finite-Q object-number work.
These checks validate inputs and identities, not a physical formation rate.
No new relic threshold, Boltzmann solution, phase holder, cosmological point,
distinctive observational signal or new consequence of closure is claimed.

[Reproduction](REPRODUCE.md) · [verification](VERIFICATION_RESULTS.json) ·
[source pins](SOURCES.json) · [literature scope](SOURCES.md) ·
[review handoff](CLAUDE_HANDOFF.md).
