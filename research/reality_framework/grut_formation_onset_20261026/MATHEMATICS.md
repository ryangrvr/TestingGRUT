# Equations and derivations

Natural units; canonical complex field X=(q1+iq2)/sqrt2, rho=sqrt(q1²+q2²). The supplied vacuum potential is

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2+\frac{g^2}{2}\phi^2\rho^2+\frac{\lambda_X}{4}\rho^4.
\]

The background circular family below x_c has

\[
\phi_0^2=v^2-\frac{g^2}{\lambda_\phi}\rho_0^2,\quad
\omega^2=g^2\phi_0^2+\lambda_X\rho_0^2,\quad
n_Q=\omega\rho_0^2,\quad
x_c^2=\frac{\lambda_\phi v^2}{2g^2}.
\]

In the co-rotating frame the radial u, phase w and neutral b perturbations satisfy

\[
\ddot u-2\omega\dot w+(k^2+A)u+Bb=0,\quad
\ddot w+2\omega\dot u+k^2w=0,\quad
\ddot b+(k^2+D)b+Bu=0,
\]

where A=2lambda_X rho0², B=2g²rho0 phi0, D=2lambda_phi phi0². Thus y=alpha² obeys

\[
(y+k^2+D)\{(y+k^2)(y+k^2+A)+4\omega^2y\}-B^2(y+k^2)=0.
\]

The nonzero-k zero-growth edge satisfies

\[
k_{\rm edge}^2=\frac{-(A+D)+\sqrt{(A-D)^2+4B^2}}{2}
\]

when the right-hand side is positive. Maximizing the positive real growth exponent determines k_peak; neither this formula nor a plane-wave cell supplies a nonlinear object charge.

Eliminating b statically uses A_eff=A-B²/D. Its dynamical condition is |y+k²|<<D. At phi0=0 on the critical circle, B=D=0, so it cannot be used. The exact polynomial factorizes as

\[
(y+k^2)\{y^2+(2k^2+A+4\omega^2)y+k^2(k^2+A)\}=0.
\]

A>=0 and omega²>=0. The quadratic discriminant is (A+4omega²)²+16omega²k²>=0; both quadratic roots are nonpositive, as is y=-k². There is no exponential linear growth. At k=0 zero roots are marginal and do not prove nonlinear stability.

For a=|X|/x_c<=1, the dimensionless coefficients in this benchmark are

\[
\omega^2/m_X^2=1-0.7a^2,\quad A/m_X^2=0.6a^2,\quad
D/m_X^2=6(1-a^2),\quad B^2/m_X^4=12a^2(1-a^2).
\]

The finite-response matrix tends continuously to zero growth at a=1. The static elimination takes a different singular limit.

## Cold elliptic orbit and relative Floquet map

With phi=0 the charged background has a quartic central potential. In m_X units,

\[
\mathbf q''=-\lambda_X|\mathbf q|^2\mathbf q,\quad
\mathcal E_X=|\mathbf q'|^2/2+\lambda_X|\mathbf q|^4/4,\quad
\mathcal Q=q_1q_2'-q_2q_1'.
\]

At the critical radius rho_c/m_X=sqrt300, set q2=0, p2=Q/rho_c and choose inward p1 from E_X. This matches both supplied invariants. The radial period follows independently from the effective potential lambda_X rho^4/4+Q²/(2rho²).

The Cartesian Hessian for the two charged directions is lambda_X(rho² I+2q q^T). The neutral curvature is g²rho²-lambda_phi v². Cross derivatives with phi vanish at phi=0. A full6-by6 canonical transfer is evolved. After one radial period the background rotates by Delta theta; applying the inverse rotation to both coordinates and momenta gives the relative monodromy. Exponents are log|multiplier|/T_r. Its symplectic form is checked with J=[[0,I],[-I,0]]. The neutral mode has zero linear charge contrast; charge variation is

\[
\delta n_Q=\delta q_1p_2+q_1\delta p_2-\delta q_2p_1-q_2\delta p_1.
\]

The small numerical k=0 charged exponent from splitting unit Jordan multipliers is not a claimed physical instability; phase and charge identities directly test that sector. Only the declared finite momentum band is scanned. No interacting coherent mode is identified with an asymptotic free particle.

## Extended thermal-state inequality

Let p=P/n_Q, P=n_minus+n_s+n_anticore, eta>=1 bound every allowed object's rest energy per absolute charge in units of m_X. Charged rest energy is at most eta m_X times total charged inventory. Total charged inventory is n_Q+2(n_minus+n_anticore). Therefore

\[
\rho_{\rm rest}/m_X\le\eta n_Q+2\eta P+(m_s/m_X-2\eta)n_s
\le\eta n_Q+\max(m_s/m_X,2\eta)P.
\]

The gas number is at most n_Q+2P. For each positive Bose-series term its mean energy is r K1(r/T)/K2(r/T)+3T<=r+3T; replacing T by T/l for the lth term only tightens that inequality. Summing positive terms gives kinetic energy<=3T per gas quantum. Adding this to the rest bound proves the result in REPORT.md. It presupposes the thermal gas and the supplied object-energy cap; it bounds neither a general coherent field nor an unbooked false-vacuum reservoir.
