# Mode, decay and relic dictionary

The full scalar equations and canonical convention are inherited explicitly
from 373d37e. No evolution equation is changed. Set chi=a z/H_i and
pi=dchi/da, a'=1 in radiation; physical z/H_i=chi/a and
dz/d(H_i t)/H_i=(pi-chi/a)/a^2. N_i=-2 Im(chi_i* pi_i).

## Two real Y modes

At fixed physical X and Sigma, let A=M_Y^2, B=2 mu Sigma, J=kappa X^3.
The forced response solves A y_p+B y_p*=-J. In real coordinates this is

\[
\begin{pmatrix}A+\Re B&\Im B\\\Im B&A-\Re B\end{pmatrix}
\binom{\Re y_p}{\Im y_p}=-\binom{\Re J}{\Im J}.
\]

Differentiating this matrix equation gives its velocity without a finite
difference. The independent C observable instead uses the complex inverse.
With delta y=Y-y_p and theta=arg Sigma, rotate by theta/2 and use the
co-rotating velocity delta dot y-i theta_dot delta y/2. The real eigenvalues
are omega_H,L^2=M_Y^2±2 mu |Sigma|. For each canonical real mode,
N_j=a^3 [qdot_j^2+omega_j^2 q_j^2]/(2 omega_j). This count differs from
the old single-mass E_Y,exc/M_Y convention.

Sigma's radial excitation is separated from its charge-bearing angular
motion: for r=|Sigma|,

\[
N_\sigma=\frac{a^3}{m_\sigma}
\left[\dot r^2+\lambda_\Sigma(r^2-h^2/2)^2\right],
\qquad m_\sigma=\sqrt{2\lambda_\Sigma}\,h.
\]

These are late adiabatic particle diagnostics. The new runs check their
stability across a=300–3000. They are not exact conserved quantum numbers
through the nonlinear early release.

## Vacuum widths and disjoint populations

At Sigma=h/sqrt(2),
m_H,L^2=M_Y^2±sqrt(2) mu h and m_A=6 e h. The tree formulas used are

\[
\Gamma_{L\to3X}=\frac{3\kappa^2m_L}{256\pi^3},\qquad
\Gamma_{A\to X\bar X}=\frac{e^2m_A}{48\pi}
(1-4m_X^2/m_A^2)^{3/2},
\]

\[
\Gamma_{\sigma\to AA}=\frac{m_\sigma^3}{32\pi h^2}
\sqrt{1-4r_A}(1-4r_A+12r_A^2),\quad r_A=m_A^2/m_\sigma^2,
\]

\[
\Gamma_{H\to LA}=\frac{9e^2}{16\pi m_H^3m_A^2}
\lambda(m_H^2,m_L^2,m_A^2)^{3/2}.
\]

The last formula follows from the charge-three derivative gauge vertex.
Its Goldstone limit agrees with the coupling sqrt(2) mu used in the review.
Every Y quantum is counted once in the eventual three-X population. A heavy
cascade additionally emits one vector; each radial Sigma emits two. Each
vector ultimately gives one X/anti-X pair. Energy carried by vectors is
removed from the parent sector, rather than added to an unmodified heavy
parent rest-mass budget.

The late N_H is only a central cascade estimate. A conservative bracket
allows zero to N_H+N_L cascade vectors, because early mixing/decay is not
solved on the undamped release. The lower injection in the decision is
Y daughters plus Sigma vectors, with no heavy-cascade vector at all.

## Radiation-era decay averages

For a nonrelativistic parent with width Gamma and conserved comoving
number, a^2≈2 H_i t and

\[
\langle a\rangle=\int_0^\infty da\,
\frac{\Gamma a}{H_i}e^{-\Gamma a^2/(2H_i)}a
=\sqrt{\frac{\pi H_i}{2\Gamma}}.
\]

Thus the radiation deposit in conformal units is N m <a>/H_i. This
corrects the old instantaneous H=Gamma estimate by sqrt(pi/2), before the
Y mass-eigenstate correction.

A vector born at a_b with comoving momentum q has conformal energy
epsilon_A(a)=sqrt(q^2+a^2 m_A^2) and optical depth

\[
\tau(a)=\frac{\Gamma_A}{H_i}
\int_{a_b}^a\frac{x^2\,dx}{\sqrt{x^2+(q/m_A)^2}}.
\]

The inverse optical depth and exponential survival are integrated directly.
An independent check uses
<epsilon_decay>=epsilon_birth+integral survival d epsilon_A/da da.
At the stated masses even vectors born as late as a=1000 become
nonrelativistic long before decay. The simple nonrelativistic deposit is
adequate to much better than a percent; explicit controls retain the actual
error. This does not certify the early condensate-to-particle conversion.

## Case-L proxy, reproduced independently

The review's light-scale onset is independently reconstructed from the
unchanged scalar potential:
X_c^2=lambda_phi v^2/(2g_X^2), rho_c=3 lambda_X X_c^4,
n_c=2 sqrt(2 lambda_X) X_c^3, V0=lambda_phi v^4/4.
At that onset use n_Q=epsilon n_c and total energy rho_c(1+R)+V0.
R contains disjoint injected energy, not another copy of the X remainder.

Bose integrals for X and anti-X use E=sqrt(p^2+m_X^2) and
f_±=1/[exp((E∓mu)/T)-1]. The net occupation difference is evaluated in an
algebraically odd-in-mu form, avoiding cancellation at tiny charge density.
Fixed Gauss integration is checked against adaptive reference quadratures.

The supplied lump proxy fixes mu_c=0.861195 m_X and E_c/Q=0.923316 m_X,
at Q=200000. Thermal gas plus condensed charge obey
n_Q=n_+-n_-+n_cond and rho=rho_gas+(E_c/Q)n_cond.
The subsequent proxy conserves entropy per net charge and uses the quoted
threshold 4->2 coefficients and cold-core conversion cross section.
This is the review's approximate finite-object prescription, not a new
derivation of a hot soliton equation of state or seed formation.

For p=n_anti/n_Q, a positivity-preserving exponential midpoint update solves
p_new=p_eq+(p_old-p_eq) exp[-Gamma(p_mid) Delta ln a/H]. The surviving
fraction is solved on [0,1], avoiding subtraction-induced bracket errors.
The threshold R_L is the first f=0.01 crossing, with
f=2 m_X Y_Q p/DM. Both Y_Q and the actual raw charge yield are tested.

The main calculation uses a monotone entropy-degree table. A separate
legacy replay retains the review's coarse table and endpoint update.
The old table jumps from g_s=3.5 to 3.91 and can imply negative time steps
under refinement; it is not treated as precision cosmology.

At the final temperature, an optimistic upper bound on all later
conversion optical depth is (n_Q sigma_v/(QH))_final in continued radiation
domination. The 4->2 contribution is bounded by its maximal final rate/H
divided by seven. Applying those bounds cannot rescue either case.

## Primordial vector correction

For a longitudinal mode B and Pi=z^2 B', the original Proca control evolves
(z^2 B')'+a^2m^2 B=0 with z^2=a^2m^2/(k^2+a^2m^2), a Bunch-Davies initial
state, stationary inflation and m(a)^2=m_vac^2+(m_inf^2-m_vac^2)/a^2 after
release. The frozen-angle limit gives the late suppression
(m_vac/m_inf)^2. The source normalization retains the review's approximate
GMR convention and is explicitly not a full oscillating-Higgs calculation.
Release-induced scalar/Goldstone production remains uncomputed here.
