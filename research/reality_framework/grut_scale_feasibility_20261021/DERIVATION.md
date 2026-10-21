# Scaling, identifiability and the actual population dependence

All derivations below are classical and at fixed dimensionless couplings, on
the declared m_D0=g_X v trajectory. They do not extend the action or produce
primordial charge. Gravity and cosmological history break the flat-space
similarity through independent dimensionless ratios.

## 1. Scalar similarity and the candidate velocity-ratio test

Use xbar=g_X v x, X=v F/sqrt(2), phi=v S and D=0. A fixed solution has

\[
Q=\bar Q/g_X^2,\qquad M=(v/g_X)\bar E=e g_X v Q,
\qquad R=\bar R/(g_Xv).
\]

At fixed Q and state, the same dimensionless collision problem gives

\[
\sigma_A(u;Q_i,\theta,\ldots)
=\frac{\bar\sigma_A(u;Q_i,\theta,\ldots)}{g_X^2v^2}.
\]

Here u is relative speed in units of c; A must name the relevant elastic or
inelastic observable. A suitable cross section divided by mass scales as v^-3.
For a fixed dimensionless population, phase ensemble and prescribed velocity
distributions, a ratio of the *same* averaged response in two velocity windows
cancels this factor. This is a conditional candidate test, not a prediction
of the population or halo response. If physical formation changes the population
with v, the cancellation does not remove that indirect dependence.

For distinguishable elastic scattering the usual transfer integral weights
d sigma/d Omega by (1-cos theta); a viscosity integral weights it by sin² theta.
Identical objects require exchange-consistent conventions. Merging or fragmentation
requires additional outcome kernels, for example

\[
K_{ij\to a}(u)=\int 2\pi b\,db\,u\,
\langle P_a(u,b,\Delta\theta)\rangle_{\Delta\theta},
\]

with separately conserved/accounted outgoing charge and energy, including
radiation. Neither integral is fixed by pi(R_i+R_j)^2. A halo analysis calibrated
with equal-mass elastic particles cannot be imported unchanged for arbitrary
composite inelastic outcomes.

For identical prototype masses the COM momentum is k=M u/2 in the nonrelativistic
limit, so

\[
kR=\tfrac12 e Q\bar R u.
\]

This has no v dependence. It is only an indicator of scattering kinematics,
not a bound on the quantum correction or a proof that the geometric limit applies.

## 2. Charge-density matching and expansion

Let r=Y_Q/Y_B, with separately conserved visible entropy after the specified
preparation. For nonrelativistic hidden matter with late energy-per-net-charge
e_eff m_X, the adopted all-DM calibration is

\[
\omega_h=\omega_b\frac{e_{\rm eff}g_Xvr}{m_b},\qquad
K\equiv\frac{\omega_c}{\omega_b}m_b=e_{\rm eff}g_Xvr.
\]

This uses a baryonic mass convention and a cosmologically inferred omega_c.
It does not derive r, v, e_eff or the initial state. If hidden radiation remains,
its background and perturbations need separate treatment.

In a restored homogeneous circular X orbit put y=|X|². Then omega²=2 lambda_X y,
n_Q=2 omega y. The phi spinodal occurs at 2g_X² y_c=lambda_phi v², giving

\[
n_c=\sqrt{\lambda_X}(\lambda_\phi/g_X^2)^{3/2}v^3.
\]

With n_Q=r Y_B (2 pi²/45)g_s T³ and radiation-dominated
H=sqrt(pi² g_*/90) T²/M_Pl, the crossing obeys

\[
T_c^3=\frac{n_c}{rY_B(2\pi^2/45)g_s}.
\]

At fixed e_eff, Y_B, g_s, g_*, all-DM calibration implies
r proportional to v^-1, T_c proportional to v^(4/3), and H/m_X proportional
to v^(5/3). The temperature here is the **visible** temperature at a cold-charge
spinodal; it is not the hidden thermal critical temperature near 1.56v.

The script's constant-g sensitivity bracket is deliberately not a precision SM
equation of state. It displays the size of the matching problem. A physical
population calculation requires entropy evolution and the actual supply/thermal
history, as well as a justified treatment of the much smaller H/m_X.

## 3. Two calibrations need not identify a third observable

In a monodisperse, fixed-density thin-wall limit, M proportional to vQ and
R proportional to Q^(1/3)/v. Density fixes vr when e_eff is fixed; the contact
proxy fixes v^-3 Q^-1/3. In logarithmic parameters (v,Q,r),

\[
J=\begin{pmatrix}1&0&1\\-3&-1/3&0\end{pmatrix},\qquad
J(1,-9,-1)^T=0.
\]

The matrix has rank two. Under v -> z v, Q -> z^-9 Q, r -> z^-1 r, both
combinations stay fixed, while M -> z^-8 M and R -> z^-4 R. The code checks
the rank with exact rational arithmetic and the finite transformation directly.
This is an explicit obstruction in a favorable limiting model, not proof that
all full-model observables are degenerate. An independently fixed population or
a further independent calibration can remove it.

## 4. Collision-count moments: a deliberately simplified counterexample

Let p(Q) be number-normalized. Approximate M=AQ and R=C Q^(1/3), with A and C
constant. Assume a common relative speed u and uncorrelated draws. The ordered-pair
contact area is pi C²(Q_i^(1/3)+Q_j^(1/3))². The mean per-object collision rate is
rho u Sigma_N, with

\[
\Sigma_N=\frac{2\pi C^2}{A}
\frac{\langle Q^{2/3}\rangle+\langle Q^{1/3}\rangle^2}{\langle Q\rangle}.
\]

Mass-weighting each object's collision **count** instead gives

\[
\Sigma_M=\frac{\pi C^2}{A}
\frac{\langle Q^{5/3}\rangle+
2\langle Q^{4/3}\rangle\langle Q^{1/3}\rangle+
\langle Q\rangle\langle Q^{2/3}\rangle}{\langle Q\rangle^2}.
\]

Neither is a drag, heat-transfer or viscosity coefficient. Those require the
actual collision kernel and the appropriate mass, velocity and angular weights.
The script checks these two moment expressions against explicit double sums.

The charge-weighted mean is Q_cw=<Q²>/<Q>. For a mixture at Q=a,b with desired
Q_cw=q0, choose number fraction

\[
p_b=\frac{a(q_0-a)}{(b-a)(b+a-q_0)}.
\]

For a=2e4, b=1e7, q0=1e6, this fixes exactly the same Q_cw as a monodisperse
sample at q0, yet changes Sigma_N by 3.35045 and Sigma_M by 9.52765. The
counterexample is purely analytical; the simplified thin-wall approximation is
not asserted accurate at the smaller charge. More realistic kernels introduce
additional dependence rather than making a single mean sufficient in general.

The old component CSV calculation uses the same proxy solely to illustrate
the loss of information under a mean substitution. A Q cut alone establishes
neither a soliton profile nor a cosmological population.
