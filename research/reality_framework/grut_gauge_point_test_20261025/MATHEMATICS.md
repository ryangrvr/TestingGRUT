# Pair emission and the stated-point tests

## Independent pair phase-space reduction

Expand delta V=lambda_d X^d+h.c. about X0=x0(r) exp(-i omega t).
The quadratic term is (Lambda(r) exp(-i Wt) delta X^2/2+h.c.), with

\[
\Lambda=d(d-1)\lambda_d x_0^{d-2},\qquad W=(d-2)\omega.
\]

For vacuum outgoing modes, the identical-pair golden-rule rate is

\[
R_2={1\over2!}\int d\Pi_1d\Pi_2\,2\pi\delta(E_1+E_2-W)
       |\widetilde\Lambda(\mathbf p_1+\mathbf p_2)|^2.
\]

Use the invariant two-body phase-space integral at total momentum
P=(W,q): Phi2(P)=sqrt(1-4m_X^2/P^2)/(8pi). This reduces the six-dimensional
integral to

\[
R_2={1\over32\pi^3}\int_0^{\sqrt{W^2-4m_X^2}}
 dq\,q^2|\widetilde\Lambda(q)|^2
 \sqrt{1-{4m_X^2\over W^2-q^2}}.
\]

This is independently derived and implemented here. Claude's reviewed code
instead integrates the individual momentum magnitude and the relative
angle, with a tabulated/interpolated Fourier transform. Our source is the
original Codex profile, not a new independent profile solution. The two
implementations agree within 1.8 ppm (d=5) and 4.7 ppm (d=6).

The pair source can have small total q while the individual momenta are
large and opposite. Core charge lost per event is d-2, escaping charge is
-2, and the total approximate X-number violation is d. Energy W=(d-2)omega
agrees with the core chemical-energy ledger. In the full gauge model the
Higgs sector compensates the violation of the *approximate* X number.

Use x0=vF/sqrt(2), rho=m_X r and lambda_d=1/Mred^(d-4). The code factors
these physical units explicitly and divides the core charge rate by Q and
hbar. Its Gaussian Fourier control, radial refinement and smooth-endpoint
quadrature are included. The result is **vacuum Born**, not an exact pair
rate in a coupled Bogoliubov medium. The LDA and additional medium channels
in the review do not become independent confirmations merely by adding them.

For the engineering reference age t0, the included-channel constraint is

\[
|c_d|\le\sqrt{0.01\over t_0(\Gamma_{coh,d}+\Gamma_{pair,d})}.
\]

The old coefficient bounds covered the coherent channel only. A generic
instanton translation c=exp(-S) still assumes a unit prefactor. No
instanton action or quantum-gravity Wilson coefficient is selected here.

## Full high-field potential and failure of low-energy elimination

Keep the 7737487 point. Denote its mediator coupling y=0.001 (called kappa
in the code), distinct from the old Higgs portal. With amplitudes
R_X,R_Y,R_Sigma and phases at their common minimum, the potential includes

\[
V=\lambda_XR_X^4+\lambda_YR_Y^4
 +\lambda_\Sigma(R_\Sigma^2-h^2/2)^2+M_Y^2R_Y^2
 -2yR_YR_X^3-2\mu R_\Sigma R_Y^2.
\]

The parent Young-inequality proof reserves lambda_X/2 times R_X^4 and
lambda_Y/2 times R_Y^4. It proves the full tree potential is nonnegative.
If the *additional assumed* inflationary mass is -c_H H^2 R_X^2, the full
potential remains bounded below by -c_H^2 H^4/(2lambda_X). A runaway of the
truncated negative-sextic potential is therefore not a runaway theorem for
the full model.

At the review's H=4.7e13 GeV and c_H=1, minimize the full radial potential,
including Y and Sigma quartics. Its stationary amplitudes differ greatly
from the vacuum expansion. The vacuum mediator mass/H is 0.0213; the old
quadratic-Y expansion parameter |y|R_X^2/M_Y^2 is about 1100 at the
single-field starting amplitude. The vacuum Higgs and vector are also not
uniformly heavy. The review's effective X6 window cannot be inserted into
this point without a new full-field calculation.

An explicitly tagged curvature-coupling proposal -xi_X R|X|^2, with
xi_X=-c_H/12, realizes this mass on an approximately de Sitter background.
R is approximately zero in radiation domination. This supplies a possible
covariant origin for the assumed release; it is a **new input**, and the
inflaton, reheating and full scalar-tensor background are not solved here.
The effective Planck coefficient remains positive at the calculated point.

## Phase alignment cannot supply the missing mismatch

The constant real y and mu permit a gauge in which the potential minimum
has X real positive, Y real negative and Sigma real negative. The field
equations, a real metric, the proposed wall tilt and a real curvature
coupling preserve the real-field subspace. With real initial fields and
velocities,

\[
n_i=-2\operatorname{Im}(z_i^*\dot z_i)=0
\]

for every field and every subsequent classical time. Making the angle sit
in its *own unchanged* minimum therefore supplies zero charge. This does
not rule out Affleck-Dine generation with an inflaton-dependent phase
mismatch or a CP-odd initial state. Those are the additional inputs used
by the review and absent from the stated gauge action.

The short homogeneous diagnostic evolves all three complex fields in
temporal gauge A0=0. Initial Gauss charge is zero and remains zero:

\[
a^3(n_X+3n_Y+6n_\Sigma)=0,\quad
\dot\rho+6H\sum_i|\dot z_i|^2=0.
\]

Conjugate misaligned preparations generate opposite, compensated charges.
They are not evolved through heavy-particle decay, the phi transition or
late charge freeze-out, so their instantaneous X charge is not a final yield.
The background radiation metric and Hubble-mass switch are supplied.

The two physical phase masses at the assumed inflationary stationary point
are computed after removing the gauge direction, with canonical amplitude
kinetic weights 2R_i^2. A finite-difference phase Hessian independently
checks the matrix. One physical phase is light. These masses alone are
not a CMB isocurvature likelihood or a prediction of the later X perturbation.

## Screened binding bound and the discrete charge

At fixed approximate X density n(r), choose a Higgs countercharge.
The minimum positive extra energy for a fixed Higgs modulus is obtained from

\[
(-\nabla^2+m_A^2)A_0=e n,\quad m_A=6eh,\qquad
\Delta E_{min}={e^2\over2}\int{d^3k\over(2\pi)^3}
 {|\widetilde n(k)|^2\over k^2+m_A^2}.
\]

It obeys the useful variational upper bound

\[
0\le\Delta E_{min}\le\Delta E_{local}
 ={1\over72h^2}\int n^2\,d^3x.
\]

An admissible initial slice with E=0 and n_Sigma=-n/6 attains the local
trial cost and satisfies classical Gauss pointwise. Y can be zero in this
trial; subsequent Y relaxation is not assumed stationary. For the spectral
minimum, the relative local-approximation error is bounded by

\[
{\Delta E_{local}-\Delta E_{min}\over\Delta E_{local}}
 \le {\int|\nabla n|^2\over m_A^2\int n^2}.
\]

Using the supplied ground profile proves a large classical binding-energy
budget remains. It does not prove a complete coupled fluctuation/Floquet
spectrum or a quantum lifetime. The minimum new vacuum heavy-particle
mass also exceeds this object's entire energy, so production of an on-shell
new heavy quantum is closed for this cold isolated object.

Continuous screening does not remove residual discrete gauge charge.
Q=200000 is 2 modulo 6. A compact quantum whole needs additional sectors;
three such cores have Q_X=600000 and can be compensated by an integer
Sigma number charge -100000. The classical local countercharge need not
be a separate number eigenstate for each core. No complete quantum
population wavefunction is constructed here.

## A separately proposed bias and the leftover budget

Let x=phi/v and choose an **additional** Delta=1e-22 GeV^4:

\[
V_{bias}={\Delta\over4}(2-3x+x^3),\quad
V_0+V_{bias}={(x-1)^2\over4}
 [\lambda_\phi v^4(x+1)^2+\Delta(x+2)].
\]

Both +/-v stay stationary, with energies 0 and Delta. The bracket's
minimum is Delta-Delta^2/(4lambda_phi v^4)>0 at this point. The false
minimum has positive curvature. The true reference energy is zero by
construction, not by a derived solution of the vacuum-density problem.
Old exact virial/source identities change under this extension even though
the perturbation is tiny. Its extra trial energy on the old core is reported.

For a scaling wall network, t_ann≈C_ann A sigma/Delta is a literature
estimate. We vary A, C_ann and g* and compare with a stated 1-MeV target;
we do not simulate a biased network. Recent simulations find different
scaling in some models. The terminal wall energy is also reported; it is
not the entire phi-transition/fragmentation energy budget.

At unchanged calibrated rho_DM/s and net yield, an extra cold free pair
yield must satisfy 2m_X Y_pair <= delta (rho_DM/s), for an explicitly
allowed fractional excess delta. Anti-object binding, annihilation,
cannibalization and entropy transfer need their own dynamics. A sudden
quench of the phi potential into zero-momentum scalars gives a useful
debris benchmark, not a lower bound: adiabatic tracking could suppress
that source, and an existing total-energy source ledger may already contain
it. It must not be added a second time.
