# Position selection in the late pole

Natural units. M is the reduced Planck mass. Source S1 is the literal Einstein
action (5.2) of arXiv:1905.00045v2. This derivation keeps its kappa=3c^2/2
normalization, specializes to c=1, and treats a post-inflaton small-positive-
sigma interval. It does not derive the initial ensemble.

## 1. Controlled action and its minimum

Set x=varphi/M and s=sigma/M. Before electroweak breaking, h=0 has
V0=lambda*v^4/4. The parent remnant potential and log-coordinate metric at
x=0, c=1 are

\[
U_{\rm parent}=A(1-e^{-s})^2+V_0e^{-2s},\quad A=\alpha^2M^4,
\]
\[
\frac{d\ell^2}{M^2}=K(s)d(\ln s)^2,\qquad
K(s)=\kappa s^2+b^2\left[\frac{s}{2\sinh(s/2)}\right]^2
=b^2+(\kappa-b^2/12)s^2+O(s^4).
\]

The exact parent minimum is s_min=ln(1+V0/A). Use s0=V0/A and
u=ln(s/s0). The fractional minimum correction is -s0/2+O(s0^2), below
1.51e-56 here. The leading controlled action has K=b^2 and

\[
U=V_0+A s^2-2V_0s
=V_0-\rho_*+\rho_*(e^u-1)^2,
\qquad \rho_*=A s_0^2=V_0^2/A.
\]

It is everywhere positive since V0>rho*. The derivative, evaluated separately
from the large constant in the implementation, is

\[
U_u=2\rho_*e^u(e^u-1),\qquad U_{uu}(0)=2\rho_*.
\]

For s below the conservative tested bound 5.455e-15, the metric's relative
error is below 2(kappa/b^2+1/12)s_max^2=1.49e-29. Each nonconstant
potential component and its force has a Taylor correction bounded by
4s_max=2.19e-14. This is a component estimate; a relative-force error at a
stationary point would be undefined. The exact minimum formula handles that
comparison. It is unnecessary and inaccurate to subtract enormous floating
constants to resolve the tiny restoring force.

The inflaton is heavy, not exactly absent from its original equation when
sigma moves. At small x, its equation has the derivative source

\[
\ddot x+3H\dot x+m_I^2x
\simeq-\frac{b^2}{2\kappa}\dot u^2,
\quad m_I^2=2\alpha^2M^2/\kappa.
\]

An adiabatically sourced x has |x| at most 3H^2/(kappa*m_I^2), using the
positive Friedmann energy budget. The quoted conservative correction scale
6H_max^2/(kappa*m_I^2) is 8.95e-15 for these starts. No heavy inflaton
oscillation population is included; the interval begins after its decay.
This classical heavy-field estimate assumes slowly varying sources on
its mass timescale. It is not a loop bound.

## 2. Conserved scalar, radiation and gravity

For an expanding flat FRW spacetime the controlled system is

\[
\ddot u+3H\dot u+U_u/(b^2M^2)=0,
\quad \dot\rho_r=-4H\rho_r,
\]
\[
3M^2H^2=\rho_r+\tfrac12b^2M^2\dot u^2+U,
\quad \dot H=-\tfrac12b^2\dot u^2-\frac{2\rho_r}{3M^2}.
\]

Thus

\[
\dot\rho_{\rm total}
=-3H b^2M^2\dot u^2-4H\rho_r.
\]

There is no prescribed changing potential, omitted energy sink, or forced
radiation-only background. Let N=ln(a/a_i), v=u_N and
rho_r=rho_ri exp(-4N). Eliminating H by the constraint gives

\[
H^2=\frac{\rho_r+U}{M^2(3-b^2v^2/2)},\quad
\epsilon=-H_N/H=\tfrac12b^2v^2+2\Omega_r,
\]
\[
u_N=v,\qquad v_N=-(3-\epsilon)v-U_u/(b^2M^2H^2).
\]

An additional integrated log H, with derivative -epsilon, is compared against
the algebraic H; it is not used to impose that constraint twice. The maximum
log-H discrepancy over all 34 trajectories is 3.94e-12.

A second organization integrates cosmic-time Klein--Gordon and the scale
factor. Define H_ref^2=rho_ri/(3M^2), tau=H_ref*t, W=du/dtau and
z=ln(1+2tau). Then

\[
h^2=(H/H_{\rm ref})^2=e^{-4N}+U/\rho_{ri}+b^2W^2/6,
\]
\[
\frac{dN}{dz}=\tfrac12e^z h,\quad
\frac{du}{dz}=\tfrac12e^zW,\quad
\frac{dW}{dz}=\tfrac12e^z[-3hW-3U_u/(b^2\rho_{ri})].
\]

The independent clock stops at N=ln(T_i/T_stop). It shares the physical
potential and inputs; it is a numerical formulation control, not a blind
derivation. Four endpoints agree in u to 5.55e-13 or better and in H to
5.54e-13 fractionally. Separately, the read-only verifier quadratures the raw
Raychaudhuri and momentum identities,

\[
\ln(H_i/H_f)=\int\epsilon\,dN,
\quad [e^{3N}H u_N]_i^f
=-\int_i^f\frac{e^{3N}U_u}{b^2M^2H}\,dN.
\]

The bath uses rho_r=(pi^2*g/30)T^4 with constant g. Stabilizing h=0 is an
explicit assumption. For example, an equilibrium free energy
F=e^(-2s)V(h)+d*T^2*h^2/2-pi^2*g*T^4/90 has an exact h=0 solution;
entropy conservation there gives T proportional to a^-1 and precisely the
radiation and false-phase energies above. A positive d can keep that phase
stable down to the stop temperature. No value of d or realistic SM thermal
history is derived or fitted, and the transition is not executed.

## 3. Field-distance obstruction independent of slow roll

Positive nonkinetic energy implies

\[
\tfrac12 b^2M^2H^2u_N^2\leq 3M^2H^2,
\quad |\Delta u|\leq\sqrt6\Delta N/b.
\]

Consequently a necessary condition for reaching u=0 from u_i is
Delta N >= b|u_i|/sqrt(6). This remains necessary for any initial velocity
consistent with the positive-energy constraint in this regime. The full
metric and adiabatic-inflaton corrections above change its coefficient
only negligibly in the specified domain. It does not constrain an earlier
stage with an entirely different metric or extra negative-energy components.

At b=3, u_i=90 and 95 require at least 110.23 and 116.35 e-folds. The actual
N budgets are 15.65 and 22.56. This is an impossibility of sufficiently rapid
motion within the stated interval, not an assertion about all cosmic histories.

## 4. Resolve the local response below floating-point precision

Let chi=b*M*u. Near the false-phase minimum,
m_chi^2=2rho*/(b^2M^2). For an infinitesimal initially stationary displacement,
the linear equation is

\[
u_{NN}+(3-\epsilon_{\rm bg})u_N+
\frac{m_\chi^2}{H_{\rm bg}^2}u=0,\quad
H_{\rm bg}^2=\frac{\rho_{ri}e^{-4N}+V_0-\rho_*}{3M^2}.
\]

Replacing V0-rho* by V0 changes the quoted leading response by a relative
order rho*/V0=3.01e-56; the implementation uses that controlled approximation.
Write q=m_chi^2/H_f^2 and u/u_i=1+q R+O(q^2). Then

\[
R_{NN}+(3-\epsilon_{\rm bg})R_N=-H_f^2/H_{\rm bg}^2,
\quad R(0)=R_N(0)=0.
\]

The Green-function expression independently evaluated is

\[
R(N_f)=-\int_0^{N_f}\frac{dn}{H(n)/H_f}
\int_n^{N_f}\frac{e^{-3(\ell-n)}d\ell}{H(\ell)/H_f}.
\]

For pure radiation epsilon=2, and q_i=m_chi^2/H_i^2,

\[
\frac{u_f}{u_i}=1-q_i\left[
\frac{e^{4N_f}-1}{20}-\frac{1-e^{-N_f}}5
\right]+O(q_f^2).
\]

The leading contraction is q_f/20 when the interval is long. Retaining V0
gives R_f=-0.05015657297, q_f=1.024907234e-58 and
1-u_f/u_i=5.140583445e-60. The two organizations agree on R to better than
4e-15. These leading coefficients are resolved by scaling the response;
they are not inferred from subtracting two rounded equal endpoint positions.

## 5. Exponential tail, tracking and acceleration

For s much greater than s0 and A*s^2 much greater than V0, or after h=v,
the potential is W=A*s^2=rho* exp(2chi/(bM)). Define phi=-chi so that
W is proportional to exp(-ell*phi/M), ell=2/b. The standard S2 variables
X=phi_dot/(sqrt(6)MH), Y=sqrt(W)/(sqrt(3)MH), gamma=4/3 give

\[
X_N=-3X+\sqrt{3/2}\,\ell Y^2+\epsilon X,\quad
Y_N=-\sqrt{3/2}\,\ell XY+\epsilon Y,
\quad\epsilon=3X^2+2(1-X^2-Y^2).
\]

The stable radiation-scaling solution requires ell^2>4:

\[
X=\sqrt{3/2}\frac{4}{3\ell},\quad
Y=\frac{2}{\sqrt3\ell},\quad
\Omega_\phi=4/\ell^2=b^2,\quad w_{\rm total}=1/3.
\]

Its eigenvalues are -[1±sqrt(64/ell^2-15)]/2. Physical stable scaling needs
b<1; b=1 is a nonhyperbolic boundary. For ell^2<4 the stable solution instead
has Omega_phi=1,

\[
X=\ell/\sqrt6,\quad Y=\sqrt{1-\ell^2/6},\quad
w=-1+\ell^2/3=-1+4/(3b^2).
\]

Its eigenvalues are (ell^2-6)/2 and ell^2-4. Acceleration needs ell^2<2,
hence b>sqrt(2). At b=3, u_N=-2/9, epsilon=2/9 and w=-23/27.
These formulas reproduce S2; they are credited known results. The same
fixed slope cannot be both in the stable radiation-scaling range and in the
accelerating range. This is not a theorem for changing slopes, couplings,
noncanonical fluids or full modified actions. Phase-plane endpoints are
only controls on this tail limit, not calculations through the shifted minimum.

## 6. A bounded observation about missing matter traces

The original Jordan-to-Einstein transformation has a matter conformal factor
calA=exp(-s/2). For chi=b*M*u its dimensionless charge is
M*dln(calA)/dchi=-s/(2b). An additional bath trace Theta=rho-3p, after
isolating V0 in U, would add

\[
(u_{NN})_{\rm trace}=\frac{s\Theta}{2b^2M^2H^2}.
\]

If |Theta|<=4rho_total and Friedmann holds, its instantaneous magnitude is
at most 6s/b^2, or 3.64e-15 over the conservatively bounded tested domain.
This local suppression is not a complete thermal evolution or a bound on
the difference between self-consistent histories. Such a calculation must
also include the corresponding matter energy exchange; none is omitted
from an alleged traceful run here. All numerical histories in this packet
have exactly ideal, trace-free radiation.
