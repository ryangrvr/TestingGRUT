# A defined conversion-product source and its streaming scales

Parent bb211d1680d8ee8540eb7a8996ddb9b2b8c1546e; C1-HP-03/v1 unchanged.
The initial source and the cold Q=200000 objects remain imposed. This is
not the missing primordial formation calculation.

## 1. Chemical energy, average energy and mass fractions

For the static background, omega_Q/m_X=Omega=0.8611951308217225.
Its average E_Q/(Qm_X)=0.9233161354 is different. The former is checked
against adjacent-charge energies in the parent daughter packet.

Conservation for object(Q)+anti-X -> object(Q-2)+X gives, neglecting
recoil, E_out=E_in+E_Q-E_(Q-2) approximately E_in+2 omega_Q. Thus

\[
K_{out}/m_X=K_{in}/m_X+2\Omega.
\]

Cold input yields gamma_out=1+2Omega=2.7223903 and p_out/m_X=2.532076.
For the parent's scalar-at-rest daughters, gamma_in=sqrt(1.5), giving
gamma_out=2.9471351 and p_out/m_X=2.7722925. Positive omega opens this
kinematic channel; it does not imply a nonzero or large amplitude for
every possible action or background. The parent computes our rate.

The parent's imposed source has 6260.9014 pairs per object, m_X=.0924090
GeV, object mass 17064.5375 GeV and total late rest mass 18221.7218 GeV
before the new conversion completion. Hence the initial anti-X rest-mass
fraction is 3.17513%, not half the binding-energy fraction.

Complete net conversion transfers 12521.8029 core charges and releases
996.5120 GeV of core chemical energy in the frozen-profile approximation.
This is 5.46881% of the original calibrated rest-matter budget. After
kinetic energy redshifts away, the rest budget is 94.53119% of that old
calibration; converted X comprise 3.35882% of the remaining rest matter,
and all free X comprise 6.71764%. Energy injection at birth and a late
warm rest-mass fraction are different quantities. The original charge
yield and v are not retuned. Changing the core profile remains necessary.

## 2. First conversions and inverse conversions

A particle which has never converted redshifts as u(a)=u_0 a_f/a.
Its first-event survival is exactly, within the fixed-object model,

\[
S_0(a)=\exp\left[-\int_{a_f}^{a}
\frac{s_0Y_Q}{Q_0 a'^3 H(a')}
(\sigma v)_{-+}(u_0a_f/a')\,d\log a'\right].
\]

This is a conversion fraction, not a permanent-accretion fraction.
Inverse reactions cannot affect the probability of a *first* event;
they do affect the subsequent population and the final energetic source.
The optical integral is performed piecewise at the rate-table spline
knots, with 32- and 64-point Gauss-Legendre controls.

For the reactive population define nu=E/m_X-qOmega and
delta=nu-(1+Omega). Then

\[
E_-/m_X=1+\delta,\quad p_-/m_X=\sqrt{\delta(\delta+2)},
\quad E_+/m_X=1+2\Omega+\delta.
\]

Forward and reverse conversion preserve delta. The exact positive
two-state reaction update is retained from the parent, with detailed
balance (sigma v)_+/(sigma v)_-=p_-E_-/(p_+E_+).
The grid includes delta=0, logarithmically resolved small delta, and a
fine high-energy section. Stable kinetic-energy arithmetic avoids
subtracting nearly equal rest energies for very cold anti-X.

## 3. A conserving receiver at the inverse threshold

Once a positive particle redshifts below delta=0, reverse conversion
into a free anti-X is closed. In this two-channel calculation it leaves
the reactive population at

\[
E_{esc}/m_X=1+2\Omega,\qquad
u_{esc}=\sqrt{(1+2\Omega)^2-1}=2.532076.
\]

For a half step beginning at scale factor a, a node with positive
momentum u crosses at a_esc=a*u/u_esc, if that point lies inside the
step. Deposit its number into an escape-time receiver with
q_comoving/m_X=a_esc*u_esc. Subsequent physical momentum is exactly
q_comoving/a. It is never rounded onto a finite energy-grid rest bin.

The receiver is an output of the imposed-source kinetic calculation,
not a receiver for the unresolved soft field of the formation simulation.
It neglects nonlinear capture, changing core profiles, additional free
particle reactions and core recoil. Elastic deflections may still occur
after reverse conversion closes, so this is the end of the modeled
inelastic stage, not a proof of all-interaction decoupling.

The receiver supplies a time/momentum source. Its compact numeric rows
have explicit column names. Original positive daughters begin below the
inverse threshold and are propagated separately from their source epoch.

## 4. Horizons, path bounds and the rejected diagnostic

For each escaped cohort,

\[
\lambda(a_{esc},q;a_o)=c\int_{a_{esc}}^{a_o}
\frac{q/a}{\sqrt{1+(q/a)^2}}\frac{da}{a^2H(a)}.
\]

These are collisionless post-escape horizons. Angular scattering can
reduce net spatial displacement; they are not a power-spectrum cutoff
or a Lyman-alpha likelihood. The reported means are number weighted
within the converted population, which has one common particle mass.

A direct radiation-era control is

\[
\lambda_{rad}=\frac{c q}{H_0\sqrt{\Omega_r}}
\left[\operatorname{asinh}(a_o/q)-
\operatorname{asinh}(a_{esc}/q)\right].
\]

An attempted active-grid path-moment diagnostic underestimates the
no-reaction mean path by 8.10096%. Its fields are retained as rejected
diagnostics and are not used for the physical conclusions. Conservation
alone did not validate that estimator.

Instead, bound the missing traveled path before escape using causality:

\[
0\le L_{pre}\le \chi_{light}(a_f,a_{esc})
=c\int_{a_f}^{a_{esc}}\frac{da}{a^2H(a)}.
\]

Thus post-escape path <= total traveled path <= post-escape path plus
the light-cone prefix, within the stated fixed-object kinetic model.
This is a bound on traveled distance, not a lower bound on net displacement.
The prefix contributes less than one percent of the mean post-escape
horizon in the three worked scenarios. No claim about the rejected
active-grid moments is needed for that bound.

Number, net charge, core chemical energy exchange and expansion work are
checked. The comoving momentum receiver and analytic horizon are tested
independently; three grid/time/source-bin refinements measure their stability.
Source-bin omission is disclosed rather than renormalized away.

## 5. Physical scope

The background, rate table, instantaneous source and object profile are
supplied. Radiation density is the same extrapolated background used in
the earlier packet, with entropy conversion used for the 20-MeV label;
it is not a new thermal-history or BBN solution. Narrow rate resonances,
evolving profiles and a complete perturbation calculation remain open.
The free-streaming scales establish a timing-dependent transport issue;
they alone neither establish nor exclude a mixed nonthermal cosmology.
