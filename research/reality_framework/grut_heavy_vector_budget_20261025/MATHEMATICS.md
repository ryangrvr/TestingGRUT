# Primordial mode, the coupled transfer check and the energy budget

## Proca mode and normalization

Use a=1 and H_i=1 at the radiation release. The homogeneous scalar solution
defines the physical gauge-Higgs norm

\[
F^2(a)=2\left(|X|^2+9|Y|^2+36|\Sigma|^2\right)/H_i^2,
\qquad m_A(a)/H_i=eF(a).
\]

It starts at F=263.450421 and tends to 6h/H_i=240. The inflationary norm
is held at the prepared value as part of the supplied preparation. The
single-mode Proca action in dimensionless conformal time is

\[
L_k=\tfrac12\left[z^2|B'|^2-a^2\widehat m_A^2|B|^2\right],\quad
z^2=\frac{a^2\widehat m_A^2}{k^2+a^2\widehat m_A^2},\quad
\Pi=z^2B'.
\]

Thus B'=Pi/z^2 and Pi'=-a^2 mhat^2 B. During inflation a'=a^2, during
radiation a'=1. Deep inside the inflationary horizon initialize u=zB in
the Bunch-Davies vacuum and Pi=z[u'-(z'/z)u]. In the adiabatic late regime,

\[
n_k=\frac{|\Pi|^2/z^2+a^2\widehat m_A^2|B|^2}
{2\sqrt{k^2+a^2\widehat m_A^2}}-\frac12.
\]

There is exactly one longitudinal polarization. Its finite-band number is
integral dln(k) k^3 n_k/(2 pi^2); no extra factor of two is used. We integrate
k/k_star from .2 to 5, k_star=sqrt(m_A,vac/H_i). This is a subset of the
spectrum, not an extrapolation of an abundance fit. Signed numerical
occupations are retained; all computed band occupations are positive.

The actual background norm is used to a=40, with a matched constant conformal
excess F^2-F_vac^2 proportional to a^-2 afterward. A background to a=60,
smooth redshifting norm, constant-vacuum mode and mode endpoint/start/tolerance
controls are supplied. The primary calculation does not substitute a fixed
vacuum mass during inflation.

The Proca action and inflationary mechanism follow Graham, Mardon and
Rajendran, arXiv:1504.02102, section III A, especially equations (3)-(18).
Their constant-mass abundance formula is not used for our normalization.
The time-dependent Higgs and additional physical scalars motivate the next
check; a free Proca action alone is not the complete multiplet theory.

## Six-real-scalar Coulomb-gauge common-mode propagation

Let q be the six canonical real conformal scalar components, with each
complex chi=(q_1+i q_2)/sqrt(2). T is the real antisymmetric U(1) generator
with charges 1, 3, 6. For the background define g=T q and g'=T q'. Work in
Coulomb gauge, so the longitudinal fluctuation is in the scalar phases and
A0 is nondynamical. At linear order around the charge-neutral homogeneous
background,

\[
J=g\cdot\delta q'-g'\cdot\delta q,\quad
A_0=-\frac{eJ}{k^2+e^2g^2}.
\]

Eliminating A0 gives

\[
L_2=\frac{\delta q'^2}{2}
-\frac{\delta q^T(k^2+M_6)\delta q}{2}
-\frac{e^2J^2}{2(k^2+e^2g^2)},
\]

where M6 is the canonical Hessian of the full scalar conformal potential.
The equivalent Hamiltonian is

\[
H_2=\frac{p^2}{2}+\frac{\delta q^T(k^2+M_6)\delta q}{2}
+\frac{e^2}{2k^2}(g\cdot p-g'\cdot\delta q)^2.
\]

This retains the longitudinal constraint and scalar/phase mixing during
release; it does not set A0 to zero for perturbations. For numerical
stability subtract the exact common-angle zero-momentum solution:
delta q=g+w, p=g'+r, with one initial common angle normalized to 1. Set
lambda=e^2/k^2 and j=g.r-g'.w. Then

\[
w'=r+\lambda gj,\qquad
r'=-(k^2+M_6)w-k^2g+\lambda g'j.
\]

The Ward identity M6 g=T partial(V)/partial(q) follows from gauge covariance.
It is used to evaluate the common component without cancellation of large
Hessian entries. The analytic Hessian is checked against finite differences,
and the Schur reduction/Legendre transform against independent real-vector
algebra. Background and perturbations are integrated together with the
author's C DOP853 implementation; no interpolated background Hessian is used.

Inflation supplies the common-mode B and B'/B. At release w=0 and
r=g k^2/(k^2+e^2g^2) (B'/B). We evolve this one initially independent mode,
then project the common angle and its derivative at a=10,20,30,40 and
continue with the late Proca equation. The 49-node band uses a=40 matching.
The three momentum/time/tolerance controls change occupations by at most
0.05%; the coupled and Proca integrals on the same grid differ by 0.0011%.

These calculations propagate the initial common fluctuation through the
full scalar/longitudinal linear Hamiltonian during the release. They do not
evolve every initial physical mode or a full interacting quantum state.
In a stationary canonical, gauge-invariant held preparation the common
vacuum mode is initially independent of the orthogonal physical modes;
their omitted late covariance contributions are nonnegative. The late
projection and supplied preparation remain part of the calculation's
domain. Correlated nonstandard initial states or a different holding action
are not excluded.

## Decay deposit and its dependence on e

At the normalized point m_A=3.2003909e6 GeV. For the open scalar-pair channel,

\[
\Gamma_{A\to X\bar X}=\frac{e^2m_A}{48\pi}
\left(1-\frac{4m_X^2}{m_A^2}\right)^{3/2}.
\]

Nonrelativistic vectors in radiation have decay-averaged scale factor
<a>=sqrt(pi H_i/(2 Gamma_A)). A comoving particle number N_A therefore gives
R_A=N_A (m_A/H_i)<a>/E_X. At this point one vector in the dimensionless
number normalization gives R=19.641987. The momenta in our band are negligible
against the rest mass by decay. One vector produces one pair; its energy
is m_A, not another factor of two times m_A.

For fixed h/H and light vector during inflation, the number scales
approximately as e^-1/2. Each particle's deposit also scales as e^-1/2,
because m_A is proportional to e and Gamma_A to e^3. Thus R_A scales
approximately as e^-1. The energy abundance's sqrt(m_A) law must not be
mistaken for the particle-number law. A single e=1e-6 sensitivity run checks
this regime; it does not establish a viable replacement point.

## The corrected benchmark and finite objects

The old isentropic proxy continues a chemical-equilibrium thermal state
after number changing freezes and omits the heat from out-of-equilibrium
annihilation. The pinned review's ke comparator uses two chemical potentials,
keeps the actual antiparticle number and conserves energy in each reaction
step. We replay its implementation with mu_c=E_c=0.861195m_X and its
Maxwell-Boltzmann on-shell scalar route enabled. This is a credited code
replay, not our own independent kinetic solver.

For actual finite Q objects, E(Q)/Q need not equal dE/dQ=omega. For cold
object density n_L and charge density n_L Q,

\[
d[n_LE(Q)]=\omega\,d(n_LQ)+[E(Q)-Q\omega],dn_L.
\]

The missing term is object-number/surface work. At the reference Q it is
1148.11 GeV per object; E/Q exceeds omega by 7.21%. The review's 16% refers to
an evolution's energy-work defect, not this per-charge ratio. A fixed-Q
linear proxy cannot use both energies without that extra thermodynamics.
Setting E_c=mu_c is a consistent linear comparator, not a derived finite-Q
population equation of state.

The finite-band vector pair rest-mass fraction would be only about 0.00150
if its kinetic energy simply redshifted without reprocessing. The failure
here is therefore a **hidden-heating/relic-proxy** result. Thermalization,
actual formation and the subsequent kinetic history must not be silently
replaced by an observational exclusion or a raw pair-count bound.
