# Local-plume calculation and conventions

All variables below are dimensionless in units of `m_X`; restore the stated
powers of `m_X` for physical units. The dimensionless neutral mass and
coupling are `M=sqrt(6)` and `g_h=0.2`. The quoted neutral width is
`Gamma=g_h^2 sqrt(1-4/M^2)/(16 pi M)`.

## Ballistic matching

Let `F_s(t,E)` be the positive-projected, windowed number flux of direct
species `s` per core per time and energy. It is an assumption that this
classical flux can be used as an incoherent particle source. If emissions
at `r0` travel radially and collisionlessly at speed `v(E)`, their local
number density per energy is

    f_s(r,t,E) = F_s(t-(r-r0)/v,E) / [4 pi r^2 v],    r >= r0.

For direct opposite-charge beams, set `s=2+2(E_+ E_- - p_+ p_-)` and

    sigma_s v_Moller = g_h^4 beta(s) /
      [32 pi E_+ E_- ((s-M^2)^2+M^2 Gamma^2)],
    beta(s) = sqrt(1-4/s).

Its narrow-width limit is `pi*g_h^2 delta(s-6)/(2 E_+ E_-)` for the
defined partial width. With `y_-=acosh(E_-)` and `a=acosh(2)`, allowed
positive rapidities are `y_+=y_-+a` and `y_-=a+y_+`. At each valid root,
energy integration gives the coefficient

    pi*g_h^2*p_+ / (4 sqrt(3) E_+ E_-).

Integrate that coefficient against `F_-` and the retarded `F_+` over the
antiparticle energy, emission time and radius `r0 <= r <= d/2`, with
`1/(4 pi r^2 v_- v_+)`. The root and retardation must be evaluated within
the energy bins; midpoint bins can alias a narrow time-overlap peak.
`continuous.py` uses Gauss quadrature across a piecewise-linear energy
model. `fullwidth_continuous.py` independently integrates the full
Breit–Wigner through an adapted plus-energy quadrature. Their endcap
conventions differ, and the full-width integral omits positive-energy
endcaps; see the report rather than interpreting their ratio as a bound.

## Explicit neutral daughters

In the neutral rest frame each daughter has `E*=M/2`, momentum
`p*=sqrt(M^2/4-1)` and speed `beta*=1/sqrt(3)`. For an outward parent of
boost `gamma, beta`, the daughter's lab energy and radial momentum are

    E = gamma E* (1 + beta beta* cos(theta*)),
    p_r = gamma E* (beta + beta* cos(theta*)).

Inward motion requires `cos(theta*) < -beta/beta*`; it is impossible when
`beta >= beta*`, equivalent to `E_h>=3`. Moreover, the same inequality
implies `E < E*/gamma <= E*`. Geometry gives a further impact-parameter
restriction for hitting an inner sphere. `neutral_return.py` integrates
that solid angle and the exponentially distributed parent path length.
`daughter_cloud.py` propagates both daughters and integrates their
opposite-charge primary exposure along their worldlines. Explicit daughters
and a full resonant kernel should not be counted a second time for the same
physical event.

## Poisson-center ensemble

For independent Poisson core centers of density `n_c` and identical,
synchronized radial source profiles `n_s^(1)(r,t)`, Campbell's identities
give

    <n_s> = n_c integral d^3r n_s^(1),
    Cov(n_+,n_-) = n_c integral d^3r n_+^(1) n_-^(1),
    Var(n_s) = n_c integral d^3r [n_s^(1)]^2.

The products are pointwise densities; neither covariance nor a first-passage
radius determines a collision rate. The **distinct-center** term of the
second factorial moment does factorize. Isotropy of independently placed
centers makes their relative direction uniform in this particular ensemble.
Its narrow-width angular kernel is

    K_NWA(E1,E2) = pi*g_h^2/(8 E1 E2 p1 p2)
                   if |(E1 E2 - 2)/(p1 p2)| < 1,
                   and zero otherwise.

Let `N_s(E,t)` be the cumulative spectrum emitted across `r0`. Then the
ensemble-mean distinct-center event rate per core is

    n_c integral dE1 dE2 N_+(E1,t) N_-(E2,t) K(E1,E2).

The full-width kernel is integrated in relative angle by changing variables
to `u=atan((s-6)/(M Gamma))`; it contains the pole once. The model includes
only direct primaries and takes their source to vanish after `m_X t=10000`.
The rate does not establish local mixing, an event-by-event bound or a
physical population model.
