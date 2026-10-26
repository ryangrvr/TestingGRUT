# A source identity is not yet a lower bound on decay debris

This is a check of an inference in the proposed structural argument, not an
independent rerun of the review or a scan for a surviving cosmology.

## Exact currents

The fixed action contains canonical complex X,Y,Sigma with charges 1,3,6,
and V_int=K+K*+L+L*, with K=kappa Y*X^3 and L=mu Sigma*Y^2. In the published
number convention n_i=-2 Im(z_i* D_t z_i), the field equation implies

\[
(\partial_t+3H)n_i=2\operatorname{Im}
\left(z_i^*\frac{\partial V}{\partial z_i^*}\right).
\]

All radial terms drop out of these sources. The three source coefficients
are (-6 Im K, 2 Im K-4 Im L, 2 Im L). The gauge-weighted sum vanishes exactly.
These statements also hold for complex couplings. In inhomogeneous fields,
the corresponding current divergences include spatial flux; a local-cell
number is not generally conserved in the absence of a vertex.

Let Q_D=a^3(n_X+3n_Y) for a homogeneous comoving unit volume. Then

\[
\dot Q_D=-12a^3\operatorname{Im}(\mu\Sigma^*Y^2),\quad
|\Delta Q_D|\leq12\int dt\,a^3|\mu\Sigma|\,|Y|^2.
\]

This is necessary but not sufficient. Phase cancellation can reduce the
integrated charge. No finite-time duration, amplitude cap, yield requirement
or residual-occupation relation has been inserted into this inequality.

When m_-^2=M_Y^2-2|mu Sigma| is positive, the positive quadratic form

\[
\rho_Y^{(2)}=M_Y^2|Y|^2+2\operatorname{Re}(\mu\Sigma^*Y^2)
\geq m_-^2|Y|^2
\]

gives a conditional integral energy bound

\[
|\Delta Q_D|\leq12\int dt\,a^3
\frac{|\mu\Sigma|}{m_-^2}\rho_Y^{(2)}.
\]

This rho is not the late particle excitation energy: it includes the forced
Y configuration, omits its linear X coupling, and is evaluated during charge
transfer. If the gap closes, even this auxiliary bound needs replacing.

## Forced response in the already-published matching regime

Use A=M_Y^2, B=2 mu* Sigma, J=kappa X^3, D=A^2-|B|^2>0. Neglect the Y quartic
only while making the controlled quadratic matching. The Y equation at its
instantaneous minimum is AY+B Y*+J=0, so

\[
Y_* = \frac{-AJ+BJ^*}{D}.
\]

Substitution yields

\[
\Delta V=-\frac{A|\kappa|^2}{D}|X|^6+
\left[\frac{\mu\kappa^2\Sigma^*}{D}X^6+\mathrm{h.c.}\right].
\]

The effective phase torque is nonzero at generic X/Sigma relative phase.
The excitation quadratic form is instead

\[
\rho_{Y,\mathrm{exc}}^{(2)}=
|D_t(Y-Y_*)|^2+A|Y-Y_*|^2+
2\operatorname{Re}\{\mu\Sigma^*(Y-Y_*)^2\}.
\]

It vanishes when Y and its velocity coincide with their forced response.
Choosing that instantaneous state is not a proof that a complete solution
can track it for the required history. In a time-dependent release, the
deviation is sourced by derivatives of Y_* and by corrections to matching.
The acceleration/nonadiabaticity and any initially prepared free modes must
be computed. The old successful high-field releases did not lie in this
controlled quadratic matching regime.

The new code verifies this distinction at four small-amplitude points of
the same dimensionless coupling ratios. A full nonlinear instantaneous Y
minimum is also solved; its nonzero source is retained as the quartic
correction becomes small. These are algebraic states, not evolved releases.
They refute the inference from a nonzero source alone to a nonzero
instantaneous free occupation; they do not refute a proposed dynamical bound
restricted to adequate-yield histories.

## Higgs countercharge is not a radial-particle count

At leading derivative order after removing Y, write X=f exp(-i theta)/sqrt(2)
and Sigma=h exp(-i sigma)/sqrt(2). Eliminating homogeneous A0 gives
I=f^2 h^2/(f^2+36h^2), n_X=6I delta_dot=-6 n_Sigma, with delta=6theta-sigma.
The phase kinetic energy is

\[
\rho_\mathrm{phase}=\frac{n_X^2}{2f^2}+
\frac{n_X^2}{72h^2}.
\]

The second term is a genuine energetic cost of compensating charge. It is
not a number of radial quanta with mass sqrt(2 lambda_Sigma)h. Radial motion,
induced displacement and its final action must be evaluated separately.
The full published debris ledger already separates these quantities.

## The adiabatic alternative has a yield cost at the fixed point

The absence of an algebraic debris floor does not mean the simple adiabatic
limit solves genesis. There is a useful necessary bound within a stated
leading-EFT radiation history. Take constant broken Sigma and alpha6, zero
initial physical velocity, and the leading inflationary radial relation

\[
H_i=\sqrt{2\lambda_X/c_H}\,r_i,\qquad r_i=|X_i|.
\]

For the leading quartic comparator, chi=aX has conserved conformal energy
|chi'|^2+lambda_X|chi|^4. At a=1, chi'=H_i r_i, so its amplitude envelope is
R=r_i(1+2/c_H)^(1/4). **Assume the effective evolution obeys |X|<=R/a.**
This is not proved for the full mediator action outside perturbative
matching. With a'=H_i and dt=a da/H_i, triangle inequality gives

\[
|a^3 n_X|\leq\frac{12|\alpha_6|R^6}{H_i}
\int_1^\infty\frac{da}{a^2}
=\frac{12|\alpha_6|R^6}{H_i}.
\]

Use the same supplied prompt-radiation visible entropy as the old releases,
g_*=106.75 and reduced Planck mass 2.435e18 GeV. At fixed c_H the optimistic
yield bound scales as r_i^(7/2), including the H_i dependence of entropy.
The coefficient remains alpha6=3.55330040797e-32 GeV^-2.

| Initial r_i/M_Y | c_H | Conditional upper bound Y_X/Y_Q |
|---|---:|---:|
| 0.1 | 1 | 1.2201e-8 |
| 0.1 | 10 | 5.4890e-8 |
| 1 | 1 | 3.8583e-5 |
| 1 | 10 | 1.7358e-4 |

These are upper bounds within the envelope assumptions, not simulated
yields. They optimistically ignore phase cancellation. The r_i=.1 M_Y
examples have small matching, frequency and radial-response diagnostics;
at r_i=M_Y the Sigma frequency separation is already less strong.

Formally extrapolating the formula to Y_Q requires r_i/M_Y=18.24 or 11.87
for c_H=1 or 10. At those points, respectively, the estimated Y quartic force
ratio is about 477 or 9.16 and the Sigma fractional forced-shift estimate is
about 78344 or 1505. The extrapolations are outside the assumed regime and
are not predictions of the full theory. The formula therefore supplies no
controlled successful low-debris alternative at this point. It does not
extend that failure to changed couplings, other histories or all possible
adiabatic preparations.

The review's important task can now be stated sharply: show whether
sufficient charge necessarily pushes allowed histories out of the quiet
forced-response regime, and whether the resulting residual excitations
unavoidably exceed the processed leftover allowance.

## What a structural proof must add

Specify a domain of couplings, inflation/release histories and preparations.
Within that domain, establish that achieving the required lasting charge
forces nonadiabatic/free radial energy or some equivalent unavoidable
descendant population. Then propagate that population through the decays,
expansion, number-changing reactions and the actual-density leftover test.
The comparison concerns a late processed quantity, not simply the potential
energy present during the torque. A bound that holds only for the sudden
prepared release must retain that restriction.

No such complete bound is derived here. It remains the pending review's
scaling question. The closed-whole principle adds no missing step to it.
