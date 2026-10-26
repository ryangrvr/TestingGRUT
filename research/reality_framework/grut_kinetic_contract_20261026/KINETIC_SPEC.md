# Reproducible scope for the next kinetic calculation

Claude owns the offered kinetic solver; Codex owns the independent source,
collision-inventory and conservation checks. The assignment was posted on
issue #13, comment 6413175923. This does not establish that a job has started.

## Fix the physical point before solving

Use the exact conditional heavy preparation and normalized light sector in
the parent [INPUTS.json](../grut_heavy_vector_budget_20261025/INPUTS.json).
Keep epsilon≈−.1457, Y_Q and the 1% engineering allowance explicit. The
full scalar release was already reproduced. Its phase and holder were supplied.
Distinguish any subsequent new holder or reheating history from this point.

For the vector source, choose a stated common spectrum: the parent's computed
finite band is a conservative partial source; the review's broader grid is
another calculated input. Do not mix one spectrum's normalization with
another's decay mass or inflationary Higgs profile. Our
[DECAY_SOURCE.json](DECAY_SOURCE.json) is a fully specified partial-source
input with the two existing values of e, not a claim to the full spectrum.

A response scan in R may be useful. At fixed spectral shape, changing R by
rescaling its number is a diagnostic experiment, not a new cosmological
point. Changing e changes the number, decay epoch and redshifted momentum
together; it should use the corresponding source. Include central computed
R=40.5626 and the review's broader-source normalization, rather than only
sampling near a desired crossing.

## Keep the independent state variables

At minimum the calculation must account for X, anti-X, the populated light
scalar, the remaining parent source, and the classical/condensed component.
Possible descriptions, with their assumptions declared, include:

- Momentum distributions f_+(q), f_-(q), f_s(q), plus a field/object sector.
- Number and energy/temperature moments for each gas species, after checking
  that elastic collisions enforce the assumed spectral shape rapidly enough.
- A verified reduction of fast inverse decay that preserves both the pair
  inventory and the scalar energy reservoir.

A single temperature and a single pair chemical potential cannot be imposed
before the relaxation that is supposed to create them. If initialization is
at the reference density instead of during injection, evolve or bracket the
pre-onset relaxation. Match the releasing condensate's energy and charge as
well as the vector source; its actual free particle/antiparticle spectrum is
not fixed by homogeneous epsilon alone.

Specify the scalar expectation value and in-medium masses. The old scalar
rates and our reservoir control use vacuum masses in the broken phase. They
cannot establish that phase at every temperature. The real scalar's cubic
and quartic self-reactions may matter; evaluate or bound their contribution,
including any threshold cancellations, rather than assuming a nonzero
coupling guarantees a large number-changing rate.

## Collision network

| Process family | Required treatment |
|---|---|
| A→X+anti-X | Source number, energy and decay-time/momentum distribution |
| Elastic X/anti-X/s scattering | Energy exchange; evidence for any common-temperature closure |
| s↔X+anti-X | Forward and reverse collision terms; s number and energy storage |
| ss↔X+anti-X | Pair destruction through the scalar reservoir; reverse reaction |
| Charge-conserving 4↔2 channels | Energy-dependent amplitudes/averages; Bose and symmetry conventions |
| Relevant scalar self-reactions | Contributions evaluated or bounded in the adopted thermal phase |
| Formation, evaporation, capture and merger | Explicit kernel or a labelled uncertainty bracket, with charge and energy sources |
| Portal transfer, if included | Specified coupling and opposite sector energy source |

The existing s route already includes ss→X+anti-X and a useful fixed-T
inventory denominator. That is not the proposed error. Its missing explicit
scalar state is a limitation of the reduced proxy. A new full-amplitude
collision treatment must also remove overlap between real sequential scalar
processes and their resonant scattering descriptions. The subtraction or
equivalent consistent treatment should be stated; no double count in the
existing threshold constants has been demonstrated here.

Nonrelativistic threshold 4→2 constants are not, by themselves, the collision
kernel for injected particles with p/m_X≈77000. They can become an adequate
late limit after the distribution has relaxed. Any such replacement needs a
rate/temperature validity check and a matched high-energy treatment or bound.

## Formation remains a separate input

Conservation laws and a cold reference profile do not specify how many
objects nucleate, their Q distribution, or their capture/evaporation rates.
A homogeneous Boltzmann run can test an explicit formation law; it cannot
derive an unspecified one. If no supported kernel is available, report the
outcome under labelled alternatives, including a gas-only control and the
old instantaneous-condensation control, without calling either measured.

Keep the finite-Q energy and omega=dE/dQ distinct. Include the object-number
work term or replace the whole object sector by an explicitly linear
comparator. Track the fast X emitted by the proposed anti-X capture channel.
Calling capture “erasure” does not remove its kinetic energy.

## Checks and reported outcome

Check vanishing internal collision energy moment, conserved total charge,
stationarity of equilibrium distributions under forward/reverse kernels,
and nonnegative distributions. Refine momentum and time grids and report
residuals separately from physical pass/fail decisions. Include parent and
mean-field work when they are present. Verify the source's number and energy
moments against this packet's analytic and ODE values.

Compare the same source with a forced thermal initial state, then show which
relaxation processes justify or invalidate that shortcut. Track late s quanta
that will decay into additional pairs. Report any anti-object charge in the
extra symmetric mass. Include a no-source control so the release's own
leftover population is not attributed to vectors.

The useful output is f_sym(e, source, release spectrum, formation kernel),
with a conditional R crossing only where its monotonicity and closure
assumptions hold. A pass may justify the next specified cosmological test.
It supplies no observational evidence for the fields and completes no
vacuum-density or observer construction.
