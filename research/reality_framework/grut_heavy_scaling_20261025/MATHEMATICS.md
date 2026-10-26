# Dilation, source bound, quartic positivity and preparation

## 1. Early supplied-background homogeneity

The early potential used in the old releases is

\[
V=\lambda_X|X|^4+\lambda_Y|Y|^4+
\lambda_\Sigma(|\Sigma|^2-h^2/2)^2+M_Y^2|Y|^2
+[\kappa Y^*X^3+\mu\Sigma^*Y^2+\mathrm{h.c.}]
-c_H H_i^2 F(tH_i)|X|^2.
\]

The release profile F is a fixed function of dimensionless time. With
z_s(t)=s z(st), all masses and H_i multiplied by s, and
H_background,s(t)=s H_background(st), the equation
ddot(z)+3H_background dot(z)+partial V/partial z*=0 scales as s^3.
V scales as s^4, currents as s^3 and their sources as s^4. Equivalently
the conformal fields chi=a z/H_i have identical dimensionless equations.

For prompt visible radiation, T_i proportional to sqrt(M_P H_i) and
entropy density proportional to (M_P H_i)^(3/2). Thus Y_Q is proportional
to H_i^(3/2) at fixed normalized number. A conserved dimensionless hidden
particle number and nonrelativistic width Gamma proportional to mass give
the same a_decay because H_i/Gamma is invariant. Hence their conformal
energy deposit per initial X energy is invariant in that regime.

These are not full Einstein-matter scaling equations at fixed M_P.
The hidden-to-background energy fraction scales as s^2. Low-energy
v_phi, m_X and a fixed physical observation time also do not transform.
Finite daughter thresholds and changes in g_* can break the decay/entropy
approximations. No exact quantum renormalization-group symmetry is claimed.

At fixed core profile the included matched X6 rate is proportional to
|alpha6|^2, so it scales as s^(-4). A string solution with fixed
dimensionless couplings has tension proportional to h^2, so G tension
scales as s^2. The numerical rates and network spectrum still need the
actual physical point. These statements neither prove nor disprove survival
there.

## 2. Generalizing the earlier source bound

For real positive coefficients after phase choices and Sigma=h/sqrt(2),
let r=sqrt(2)mu h/M_Y^2 and u=M_Y/H_i. The quadratic Y determinant is
D=M_Y^4(1-r^2), positive for r<1. Then

\[
\alpha_6H_i^2=\frac{r\kappa^2}{2(1-r^2)u^2},\qquad
\beta_6H_i^2=-\frac{\kappa^2}{(1-r^2)u^2}.
\]

Complex phases are absorbed in the phase of alpha6; its magnitude appears
in the bound. The radial companion is retained. Near r=1 the smaller
vacuum gap M_Y sqrt(1-r) vanishes, invalidating any automatic inference
of adiabaticity from M_Y/H_i alone.

Write chi=aX/H_i and q=a^3 n_D/H_i^3, where n_D=n_X+3n_Y in the full
theory. At leading derivative order with Y eliminated and Sigma fixed,
the source equation in radiation (a'=H_i) is

\[
\frac{dq}{da}=-\frac{12}{a^2}
\operatorname{Im}(\alpha_6H_i^2\chi^6).
\]

Consequently, if |chi|<=C for the entire contributing epoch and alpha6
is constant,

\[
|\Delta q|\leq12|\alpha_6H_i^2|C^6.
\]

One can instead integrate the actual envelope and coefficient, including
the a^-2 weight. A full-model check must also quantify matching/derivative
corrections and deviations of Sigma. This effective expression is not the
exact full-theory source when real Y motion matters.

For the earlier leading quartic comparator with zero physical initial
velocity, r_i/H_i=sqrt(c_H/(2 lambda_X)) and
C=sqrt(c_H/(2 lambda_X))(1+2/c_H)^(1/4). The quartic-energy normalization is
N_star=2 sqrt(2 lambda_X) C^3/3^(3/4). A target epsilon therefore needs

\[
|\alpha_6H_i^2|\geq
\frac{|\epsilon|N_\star}{12C^6}.
\]

This is optimistic because it allows maximal same-sign torque throughout.
The actual E_X and C may differ from the comparator; these assumptions
must be tested before applying its numerical table to a full trajectory.

## 3. A separate exact quartic condition

Choose the X/Y phase that minimizes the kappa term and set |X|=x, |Y|=t x.
The X/Y quartic per x^4 is

\[
P(t)=\lambda_X+\lambda_Yt^4-2|\kappa|t.
\]

For positive lambda_Y its minimum occurs at
t=(|kappa|/(2 lambda_Y))^(1/3), giving
P_min=lambda_X-3|kappa|^(4/3)/(2^(4/3)lambda_Y^(1/3)). Thus

\[
P(t)\geq0\ \forall t\geq0
\quad\Longleftrightarrow\quad
27|\kappa|^4\leq16\lambda_X^3\lambda_Y.
\]

With lambda_Sigma>0 and strict positivity in the X/Y sector, the highest
quartic form is positive and the lower-degree trilinear cannot make the
full tree potential unbounded. The equality case has flat quartic
directions and requires its own lower-degree analysis. No claim of safe
full-potential boundedness at equality is used here.

Combining this necessary ceiling with the conditional envelope gives

\[
|\epsilon|\leq
3^{1/4}[c_H(c_H+2)]^{3/4}\sqrt{\frac{\lambda_Y}{\lambda_X}}
\frac{r}{(1-r^2)u^2}.
\]

Only the quartic positivity condition is exact for the stated quartic.
The combined epsilon inequality inherits every preparation/envelope and
matching assumption. It is not a structural no-go for the full field
content.

## 4. Exact linear preparation example

Consider a canonical real mode obeying
ddot(q)+m^2[q-q_*(t)]=0. During 0<t<tau let
q_*(t)=q0[1-3(t/tau)^2+2(t/tau)^3], constant q0 before and zero after.
Its first derivative vanishes at both endpoints. Put x=m tau, z=x/2.

For initially relaxed q(0)=q0, dot(q)(0)=0, the retarded solution gives

\[
q(\tau)+\frac{i\dot q(\tau)}m
=q_0e^{-ix}F(x),\qquad
F(x)=\int_0^1du\,6u(1-u)e^{ixu}
=e^{ix/2}\frac{3(\sin z-z\cos z)}{z^3}.
\]

With q(0)=0 instead, subtract q0 e^(-ix). Final free energies in units
m^2q0^2/2 are |F|^2 and |F-1|^2. For large x the relaxed energy envelope
falls as x^(-4), while the displaced case approaches one. Isolated
interference zeros are not robust preparation windows. The formula is
checked against independent quadrature and direct DOP853 integration.

This externally specified drive is a linear response illustration only.
It proves no complete closed-system history, charge yield, thermalization,
vector source or inflationary preparation for the GRUT fields. Its purpose
is to identify the initial-state assumption that an actual heavy-mediator
calculation must resolve.
