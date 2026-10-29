# The source-mixing calculation that the mass test still needs

The supplied candidate has `n_ref/m_X^3=1.2141071423e-15`. The typical
number-density length is

    d = n_ref^(-1/3) = 93738 / m_X.

For a synchronous experiment starting from localized pulses with initial
support radius `R0 << 10000/m_X`, a causal signal at microscopic time
`t=10000/m_X` lies inside a sphere of radius at most `R0+t` around its source.
The union bound is

    V_covered/V_total <= n_ref (4 pi/3) (R0+t)^3
                       ~= 0.00509,  for R0 << t.

This does not presume a regular lattice or forbid overlap. If the injected
particles are all contained in these regions, their average density in the
covered region is at least about 197 times their global mean. This is **not**
a lower bound on their collision rate: directions, momenta, temporal overlap,
wave coherence and correlations all matter. A nonsynchronous formation
history or existing background changes the initial conditions and must be
specified. Cosmological expansion is negligible during this local stage under
the supplied 20 MeV normalization.

The mean energy of on-shell neutral flux in the finer two-window source is
`3.05056 m_X`. Inserting that energy into the vacuum decay-length formula gives

    ell_h = p_h / (m_h Gamma_h) ~= 3957 / m_X ~= 0.0422 d.

This is one representative length, not the mean of the distribution. At such
distances neutral decay can produce charged daughters in a still-localized
cloud, including inward-going particles that can encounter the parent core.

## Why the directions matter to the resonance

For `X Xbar -> h`, write `p_i/m_X=sinh(y_i)` and `E_i/m_X=cosh(y_i)`.
If both momenta are exactly outward and collinear, then

    s/m_X^2 = 2 + 2 cosh(y_1-y_2).

The `m_h^2=6 m_X^2` pole therefore requires
`|y_1-y_2|=arcosh(2)=1.31696`. For example, two particles each with
`p=m_X` instead require a relative angle of 90 degrees to hit the pole; they
cannot do so if their directions are both within 20 degrees of the same
outward direction. Neutral decay daughters or subsequent scattering can of
course produce broader directions.

`BEAM_KINEMATICS.json` applies this purely kinematic test to independent draws
from the time-integrated **direct charged** energy marginals of the fine
source. If both momenta are restricted to a hypothetical common cone of
half-angle 10 degrees, 0.142% of those energy pairs have *some* allowed
relative angle that reaches the exact pole. At 20 degrees this is 6.14%; in a
90-degree cone it is 82.05%. These are not collision probabilities or rate
ratios: real angles and time overlap have not been measured; h-decay daughters
are excluded; and finite-bin, wave/coherence and scattering weights matter.
The calculation demonstrates why the homogeneous isotropic rate cannot be
assumed to hold inside an unresolved plume.

## Concrete continuation

1. Reconstruct a time- and angle-resolved flux from the fine radial fields
   at a matching surface outside the core, with energy and charge conserved
   across the field/particle interface. Resolve outgoing and incoming fluxes;
   preserve finite-window, off-shell and signed-spectrum errors. Test more than
   one matching radius and source duration. The present compact spectra are
   insufficient to reconstruct their joint time-angle correlations.
2. Evolve `f_X(t,r,p,mu)`, `f_Xbar`, and `f_h` around an isolated reference
   core, including free streaming, vacuum h decay and inverse decay, relevant
   nonresonant collision kernels, and core boundary scattering. The spherical
   streaming term for velocity `v=p/E` is

       partial_t f + v mu partial_r f + v(1-mu^2)/r partial_mu f = C[f] + S.

   Expansion can be added later or retained as a controlled negligible term
   over this short stage. Evolve the core's energy and charge consistently,
   including interactions of inward-directed daughters. Treat the initially
   coherent classical field and any particle approximation explicitly; do
   not manufacture an isotropic gas at the matching surface.
3. In an explicit-h kinetic description, count `X Xbar <-> h` and h decay
   once. The parent reduced resonant `X Xbar -> h -> X Xbar` operator includes
   the same on-shell intermediate state and cannot be kept unchanged alongside
   explicit inverse decay. Off-shell and interference pieces need consistent
   real-intermediate-state handling rather than an arbitrary subtraction.
4. Compare an interacting cloud with a collisionless cloud, track all energy
   and charge to a matching radius comparable to the inter-core spacing, and
   propagate the outgoing joint distribution into the multi-core FRW problem.
   The fine-source modeled energy checkpoint has only about `135.6 m_X` of
   room below the 1% budget. The sensitivity criterion should therefore be
   the effect on the **full** redshifted energy and pressure ledger, not just
   the number of antiparticles or the cold mass at an intermediate time.

The 135.6 `m_X` figure is a margin of this particular conditional model; it
is not an observational uncertainty or a bound on the actual local correction.
An actual source, population distribution, formation epoch, subsequent core
excitation dynamics and independent replication remain necessary.
