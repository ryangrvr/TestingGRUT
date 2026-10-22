# Same-action relation and a nonlinear preparation probe

M is reduced Planck mass; natural units. S1 is Camargo-Molina, Markkanen and
Scott (2019), literal Einstein action (5.2) in its PDF, (20) in HTML. This is
an analysis of a borrowed proposal, not GRUT-native physics or an adopted target.

## 1. Fix the normalization before calibrating

In the regime of negligible remnant energy and a stabilized Higgs, S1 has

\[
\mathcal L\supset-\frac34c^2(\partial\varphi)^2
-\alpha^2M^4(1-e^{-\varphi/M})^2,
\quad \psi_I=\sqrt{3/2}\,c\varphi.
\]

The metric coefficient is kappa=g c^2 with g=3/2, because canonical kinetic
energy is one half the metric. The paper's following text writes a different
canonical prefactor and eta; this calculation uses the **action as printed**.
In particular when the pole is suppressed, U_sigma=alpha^2 M^2 sigma^2 gives
eta_sigma=4/(3c^2), not 8/(3c^2), on its early plateau. This is a normalization
warning, not a claim to have independently reproduced the full model.

Let z=exp(varphi/M)>1. For single-adiabatic-field first-order potential slow roll,

\[
\epsilon=\frac{2}{\kappa(z-1)^2},\quad
\eta=\frac{2(2-z)}{\kappa(z-1)^2},\quad
r=16\epsilon,\quad 1-n_s=2\epsilon(z+1),
\]
\[
A_s=\frac{\alpha^2(1-1/z)^2}{24\pi^2\epsilon}.
\]

The end point epsilon=1 gives x_end=ln(1+sqrt(2/kappa)), x=varphi/M, and

\[
N=\frac\kappa2[e^{x_*}-x_*-e^{x_{end}}+x_{end}].
\]

The sensitivity table solves this equation, calibrates alpha from A_s, and
calculates the remnant with the **same c and kinetic normalization**. Alpha
cannot be held fixed while c is changed in this action's inflation calibration.
Finite canonical derivatives of the potential provide a separate numerical check.

## 2. Eliminate the apparent normalization freedom

For V0=lambda v^4/4 and sigma/M tiny, the pre-EW remnant potential is

\[
U_\sigma\simeq\alpha^2M^2\sigma^2+(1-2c\sigma/M)V_0,
\quad \sigma_0=\frac{cV_0}{\alpha^2M^3}.
\]

If this displacement is prepared with negligible canonical velocity and
retained after the Higgs reaches v, its initial positive energy is

\[
\rho_* =\frac{c^2V_0^2}{\alpha^2M^4}.
\]

Write delta=1-n_s. Since z=8delta/r-1, eliminating c and alpha gives

\[
\boxed{\rho_*=
\frac{64 V_0^2}{3g\pi^2 A_s M^4[8(1-n_s)-r]^2}}
=\frac{128 V_0^2}{9\pi^2 A_s M^4[8(1-n_s)-r]^2}\quad(g=3/2).
\]

This is a **conditional consistency relation**, not an unconditional density
prediction. It has the correct dimension four. The potential branch z>1
requires 0<r<4(1-n_s); r=0 in the table denotes a limit. The leading plateau
limit gives rho*≈V0^2 N^2/(18 pi^2 A_s M^4). No new free c is available to
adjust the amplitude once the same inflation data and action are used.
Changing relative remnant/inflaton potential amplitudes, kinetic coefficients,
initial population or perturbation transfer changes this relation and needs
a separately declared model/history. No such repair is adopted here.

Assumptions: single-adiabatic inflation, negligible entropy-to-curvature
transfer, standard vacuum/slow roll and tensor normalization, a prepared
sigma0, the specified common potential scale, zero absolute final vacuum
reference, and negligible loop/higher-derivative changes. These are not
consequences of closure. Reheating/history enters the N sensitivity, although
the observable relation eliminates N algebraically. If separately conserved
late scalar energy starts at rho* with negligible initial velocity, it can
only decrease in an expanding universe: rho_dot=-3H chi_dot^2. The relation
is then an upper bound on its later energy, not an allowance to increase it.

At rounded A_s=2.1e-9,lambda=.129,v=246 GeV,M=2.435e18 GeV and the dated
n_s=.9649 snapshot, r from zero to .036 gives 3.4530e-48 to 4.5433e-48 GeV^4.
Across n_s=.9649±2(.0042),r in [0,.036], the maximum is 8.6319e-48, .342536
of the known reference 2.52e-47. These are a **rectangular sensitivity exercise**,
not a joint confidence region or new likelihood. Reproducing the reference
would require n_s=.987007 at r=0, or .982507 at r=.036, with other inputs fixed.
No sigma-level exclusion is assigned, and no current-data global claim is made.

## 3. The pole changes what a fluctuation means

At fixed early inflaton value x=varphi/M and small positive s=sigma/M, S1's
remnant metric is approximately

\[
G_{\sigma\sigma}=\kappa+p^2/s^2,\qquad p=b e^{-x/2}.
\]

It can be nearly constant at sigma0 yet still has a pole at sigma=0. A global
Gaussian in sigma crosses that infinite-distance boundary and is not the
distribution specified by this action. Introduce a canonical coordinate q with

\[
q/M=F(s),\quad F(s)=\sqrt{\kappa s^2+p^2}
+p\ln\frac{\sqrt\kappa s}{p+\sqrt{\kappa s^2+p^2}}.
\]

An irrelevant additive constant fixes its origin. For large positive q relative
to the tiny transition scales, sigma≈q/sqrt(kappa). For negative q,
sigma approaches zero exponentially while remaining **positive**. This is not
an opposite-sign vacuum population. The full potential tends to a finite
constant as q→−infinity. Consequently a fixed-H equilibrium measure
exp[-8pi^2 U(q)/(3H^4)] is not normalizable on that half-line. A stationary
classical point alone does not supply a normalizable narrow quantum attractor.

## 4. Bounded stochastic limit, not a reheating result

When the minimum and pole width are negligible relative to S=H/(2pi), the
pre-EW canonical potential, after subtracting its constant, becomes

\[
U(q)\simeq\tfrac12m_{eff}^2(q_+)^2,\quad q_+=\max(q,0),
\quad m_{eff}^2=2\alpha^2M^2/\kappa,
\quad H^2\simeq\alpha^2M^2/3.
\]

For X=q/S the leading light-field stochastic approximation [S2] is

\[
dX=-a X_+dN+dW_N,\quad a=m_{eff}^2/(3H^2)=2/\kappa,
\]
\[
\partial_N P=a\partial_X(X_+P)+\tfrac12\partial_X^2P.
\]

The probe uses c=10,b=3 and alpha from that c's N=55 calibration, not an
observationally accepted point. It has m_eff/H=.2. Choose fixed x such that
p^2/(kappa s0^2)=e^-20, giving x=276.0004. The canonical minimum scale is
1.15e-50 S and the pole width 5.22e-55 S. Finite positive-noise excursions
remain sigma/M of order 1e-6, where the local potential approximation is good.
The fixed inflaton and H, light-field stochastic truncation and absent
reheating are explicit assumptions. This is not a solution through the end
of inflation or today's universe. It tests whether the assumed narrow state
is maintained in this limiting stage, even when initialized optimistically.

For a positivity-preserving finite-volume generator take cell coordinates X_i,
Phi_i=a(X_i+)^2 and B(z)=z/(e^z-1). Neighbor jump rates are

\[
Q_{i+1,i}=\frac{B(\Phi_{i+1}-\Phi_i)}{2\Delta X^2},\quad
Q_{i,i+1}=\frac{B(\Phi_i-\Phi_{i+1})}{2\Delta X^2},
\]

with diagonals chosen for zero column sum and reflecting distant endpoints.
The matrix exponential propagates the initial central delta. Refinement and
domain enlargement control discretization and artificial boundaries. Assign
half the central cell to either canonical half-line in the probability estimator.

At N=10,100 the refined positive canonical probabilities are .491670,.420515.
The asymptotic map gives rms sigma/M=1.1428e-6,2.1554e-6, respectively, or
1.81e50,3.42e50 times the selected sigma0. Its coarse-grained sigma energy is
only 1.31e-12,4.65e-12 of the background inflation energy, supporting negligible
backreaction for this probe. **These are not post-reheating relic densities.**
Their large ratio to sigma0 does not prove a complete-action empirical rejection:
later damping, changing geometry and energy transfer remain uncalculated.

## 5. Independent local mode controls

For a constant canonical mass on de Sitter, nu=sqrt(9/4-m_eff^2/H^2) and
the Bunch--Davies power is

\[
\mathcal P_q/H^2=\frac{z^3}{8\pi}|H_\nu^{(1)}(z)|^2,\quad z=k/(aH).
\]

The independent mode equation for v=sqrt(k)a q_k is
v_zz+[1-(nu^2-1/4)/z^2]v=0. Starting from the plane-wave subhorizon limit at
z=1000 and 2000 gives agreement with Hankel powers to better than 7e-7 here.
The canonical Wronskian is checked. A finite band z in [epsilon e^-N,epsilon]
has variance integral over d ln z of that positive power. Its local linear
scale vastly exceeds the tiny selected displacement; this diagnoses failure
of a narrow perturbative neighborhood, **not** permission to continue a Gaussian
through the pole. The Fokker--Planck probe handles the corresponding limiting
canonical geometry and remains an approximation to fixed-inflaton dynamics.
