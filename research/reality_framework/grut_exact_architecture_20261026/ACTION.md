# Covariant action, allowed operators and matching

This fixes the architecture's local mathematical content and identifies its
supplied inputs. It does not replace the canonical project action or claim
a numerical quantum matching calculation for the entire cosmological point.

## Scalar basis

Let z_i=(X,Y,Sigma) have U(1) charges (1,3,6), with real neutral phi. Before
field redefinitions, every nondifferentiated scalar monomial of mass dimension
at most four belongs to

\[
\begin{aligned}
V_{\rm ren}={}&c_0+t\phi+\frac12m_\phi^2\phi^2
 +\frac{a_\phi}{3!}\phi^3+\frac{\ell_\phi}{4!}\phi^4\\
 &+\sum_i(m_i^2+a_i\phi+b_i\phi^2)|z_i|^2
 +\sum_{i\le j}\lambda_{ij}|z_i|^2|z_j|^2\\
 &+[(\mu_0+\mu_1\phi)\Sigma^*Y^2+\kappa Y^*X^3+\mathrm{h.c.}].
\end{aligned}
\]

The norm-sector coefficients are real; mu_0,mu_1,kappa may be complex.
There are 26 invariant monomials including the constant and conjugate
monomials. With an exact phi→−phi symmetry the odd terms vanish and 19
remain. These counts precede field shifts, rephasings and parameter
redundancies; they are not counts of independent measured observables.
An explicit symmetry-breaking wall bias must be treated consistently with
that symmetry choice rather than silently keeping all even-only conclusions.

The enumeration in operators.py checks all seven-field exponent tuples
of total degree at most four against these generators. No extra charged
renormalizable monomial is omitted in this scalar content.

Standard Model gauge invariance additionally permits
(H_SM^dagger H_SM)[a_Hphi phi+b_Hphi phi^2+sum_i lambda_Hi |z_i|^2],
besides the terms already in the SM. Hidden/hypercharge kinetic mixing is
also an allowed datum unless an additional symmetry forbids it. Their
matching values must be stated. Merely supplying the field content does
not determine them.

The earlier tree potential is the reference choice
ell_phi=6 lambda_phi, m_phi^2=−lambda_phi v^2,
b_X=g_X^2, m_X,bare^2=0 at that tree reference,
m_Y^2=M_Y^2, m_Sigma^2=−lambda_Sigma h^2,
lambda_XX=lambda_X, lambda_YY=lambda_Y, lambda_SigmaSigma=lambda_Sigma,
mu_0=mu and kappa as listed. Other coefficients were set to zero at a
matching reference. The constants in the shifted quartics are absorbed
in c_0 together with the independent vacuum counterterm. The low-energy
X mass is a renormalized matching target, not an all-order identity enforced
by the tree condition m_X,bare^2=0.

## Gravitational and higher-dimensional terms

Write

\[
F=M_{\rm r}^2-2\xi_H H_{\rm SM}^\dagger H_{\rm SM}
 -\xi_\phi\phi^2-2\sum_i\xi_i|z_i|^2-2\zeta_\phi\phi.
\]

An odd phi curvature term can be forbidden by the same exact parity if it is
imposed. The scalar kinetic terms and F R/2, together with allowed curvature
counterterms and higher operators, define a gravitational EFT below a stated
cutoff. Coefficients are matched; a finite truncated list is not an exact
ultraviolet completion. Keeping metric and matter loops requires the
appropriate complete basis and orders for the observable being calculated.

A representative nonminimal choice xi_X=−c_H/12 reproduces the given
inflationary tachyonic mass. For ideal flat FLRW with equation of state w,

\[
R=6(2H^2+\dot H)=3(1-3w)H^2,\qquad m_{X,\rm curv}^2=\xi_XR.
\]

| Background | R/H^2 | m_curv^2/H^2 for c_H=10 |
|---|---:|---:|
| de Sitter | 12 | −10 |
| Matter | 3 | −2.5 |
| Radiation | 0 | 0 |

The sign conventions correspond to the action in THEORY.md. At the prepared
X amplitude of the normalized heavy point, this term shifts F/M_r^2 by
about 2.92e-11 and leaves F positive. This is a spectator-size check, not a
solved nonminimal-gravity background or perturbation spectrum. In a physical
transition the Ricci scalar evolves dynamically. The table does not supply
an inflaton potential, reheating width, thermal trace, clock, phase holder or
adiabatic switch-off history. It also shows why an arbitrary prescribed
−c_H H(t)^2 mass during matter domination is a different model.

## A checked scalar quantum correction

For X=(r+i s)/sqrt(2), Y=(u+i v)/sqrt(2), retain only kappa Y*X^3+h.c.
with real kappa for this coefficient check. Exact rational differentiation
of its canonical real Hessian M_s^2 gives

\[
\left.\operatorname{Tr}[(M_s^2)^2]\right|_{\kappa^2}
 =36\kappa^2|X|^4+72\kappa^2|X|^2|Y|^2.
\]

The one-loop scalar effective-potential logarithms consequently require an
X–Y mixed quartic counterterm. In the convention V includes
lambda_XY |X|^2|Y|^2, the isolated contribution is

\[
\left.\mu_R\frac{d\lambda_{XY}}{d\mu_R}\right|_{\kappa^2}
 =\frac{9\kappa^2}{4\pi^2}=2.72712\times10^{-6}
\]

at the supplied kappa. This is a full-theory running contribution where Y is
active. It must not be extrapolated below the heavy threshold without
integrating out and matching Y. No complete beta function, two-loop scalar
mass matching or physical pole correction is claimed.

Thus setting lambda_XY to zero at one scale is a legitimate matching input;
declaring it identically zero under quantum evolution is not. The original
gauge construction already identified additional portals and the light
hierarchy as open. This calculation makes one such obligation explicit.

## A gauge-potential sensitivity, not a mass exclusion

For fixed Sigma=h/sqrt(2), Y=0 and x=|X|^2,
m_A^2(x)=36 e^2 h^2+2 e^2 x. In Landau gauge and MS-bar the vector part is

\[
V_A^{(1)}=\frac{3m_A^4(x)}{64\pi^2}
 [\ln(m_A^2(x)/\mu_R^2)-5/6].
\]

At mu_R=m_A(0), its derivative is

\[
\left.\frac{\partial V_A^{(1)}}{\partial x}\right|_0
 =-\frac{e^2m_A^2}{16\pi^2}.
\]

At fixed h its magnitude scales as e^4. The ratios to the target m_X^2 are
.07596 at e=1e-7, 759.55 at e=1e-6 and 61523.8 at e=3e-6.
This is one gauge/scheme-dependent off-shell potential contribution,
**not** the pole mass, a full threshold, a predicted physical mass shift,
an imposed naturalness bound or a reason to discard those points.
Renormalized low-energy targets can be matched with counterterms. Increasing
e to reduce late vector injection therefore preserves the old light sector
only with stated matching, not by an automatic all-order identity.

No new e scan, corrected soliton, kinetic trajectory or recalibrated
cosmological point is produced in this packet. Martin's primary one-loop
formula is credited in [SOURCES.md](SOURCES.md); the derivatives and scalar
polynomial check are computed here.

## The vacuum coefficient is still an independent input

A constant shift of V changes the metric equation even when local scalar
forces are unchanged. At a constant-field vacuum and leading Einstein order,
M_r^2 G_mu_nu+rho_vac g_mu_nu=0, so R=4 rho_vac/M_r^2.
Spatially compact flat FLRW admits both a zero-density flat reference and
positive-density expanding solutions, with 3 M_r^2 H^2=rho_vac. These are
different closed-space solutions within the adopted equations. The word
closed does not choose the coefficient. No self-adjustment, one-state
encoding or de Sitter information-capacity rule is inserted into this action
to conceal the unresolved selection problem.
