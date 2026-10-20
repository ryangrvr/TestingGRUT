# Cosmological sensitivity and wall scales

Use the unchanged benchmark g_X=.1, lambda_phi=.01, lambda_X=.001, h_D=.01,
lambda_D=.001, m_D0=100 GeV, v=1000 GeV and r=1. These are adopted inputs,
not empirical determinations of a hidden sector. The classical current,
charge normalization and sharp floor were derived in the parent audit.
With canonical fields and the stated positivity conditions,

\[
\rho_{\rm exc}\ge m_* |n_Q|,\qquad
m_*=v\min[g_X,(\lambda_\phi\lambda_X)^{1/4}]=56.234132519\ {\rm GeV}.
\]

This bounds proper **excitation** energy above the normalized vacuum. A
separate total constant vacuum energy contributes with its sign to gravity.
The positivity of gradients also holds in a local FRW orthonormal frame.
Conservation of comoving net charge and visible entropy is still an adopted
history; gravitational averaging/binding assumptions are unchanged.

Let H_f=113.251817518 km/s/Mpc be the conservative floor unit fixed by
omega_b=.02125 and m_b=.94738 GeV. The lower bound on the mean excitation
source, including ordinary baryons, scales as H_f^2(1+z)^3. For separate
constant-vacuum and curvature sensitivities write

\[
H^2(z)\ge H_f^2(1+z)^3+L+K(1+z)^2,
\quad L=100^2\omega_\Lambda,\quad K=100^2\omega_k.
\]

Here the squared rates use (km/s/Mpc)^2. K<0 denotes positive spatial
curvature, not negative matter energy. In the envelope-saturating construction
H_0^2=H_f^2+L+K, and Omega_k=K/H_0^2; this relation fixes how the reported
Omega_k threshold is normalized. It is not an independent measured H_0.

Set chi=c integral_0^z dz'/H_min(z'). Then

\[
D_{A,\max}=\frac{1}{1+z}
\begin{cases}
\chi,&K=0,\\
\dfrac{c}{\sqrt{-K}}\sin\!\left(\dfrac{\sqrt{-K}\chi}{c}\right),&K<0.
\end{cases}
\]

For the closed case this is an upper envelope on the portion of the light
cone before the sine maximum. The code checks that condition on every
evaluation. H_min squared remains positive on the tested interval. More
positive sources reduce chi; for the closed intervals used here this reduces
the transverse distance as well. Beyond that interval the argument cannot be
reused without revisiting the geometric monotonicity condition.

For L=0 the independent Mattig dust expression supplies a check:

\[
D_A=\frac{2c}{H_0\Omega_m^2(1+z)^2}
\left[\Omega_m z+(\Omega_m-2)(\sqrt{1+\Omega_m z}-1)\right],
\quad\Omega_m=H_f^2/H_0^2.
\]

It agrees with direct curved-distance quadrature within 2.8e-13 Mpc at the
thresholds. An independent 32-point Gauss–Legendre integral checks the
negative-vacuum quadrature. Root-finding uses the largest D_low-D_A across
all six galaxies. It retains the original exact multiplicative Doppler
correction and median-minus-three-lower-width construction; those widths do
not define a joint confidence region. No new data or likelihood is used.

The total-vacuum energy conversion is rho_crit,100=3 M_reduced^2 H_100^2,
H_100=100/(3.0856775814913673e19) s^-1 times hbar in GeV s,
M_reduced=2.435323e18 GeV. This gives the values in CORRECTIONS.md. Curvature
thresholds with L=0 are Omega_k=-.5630192054 (primary) and -.3397438305
(stress), at H0=90.5864 and 97.8440 respectively in the floor-saturated
construction. Neither is a complete observationally acceptable universe.

For X=D=0 the exact leading planar vacuum wall is

\[
\phi(z)=v\tanh(v\sqrt{\lambda_\phi/2}\,z),\qquad
\sigma=\int_{-v}^{v}\sqrt{2U(\phi,0,0)}\,d\phi
=\frac{2\sqrt{2\lambda_\phi}}3v^3
=9.4280904158\times10^7\ {\rm GeV}^3.
\]

Its tension cube root is 455.136 GeV; its tanh width in coordinates m_X z is
sqrt(2) at this benchmark. This is a planar classical solution, not a
network abundance, a wall decorated with charge, or a lifetime calculation.
In three spatial dimensions a finite closed bubble with degenerate vacua
has energy about 4 pi R^2 sigma and can lower its energy by contracting.
Planar topology does not make every finite bubble a persistent network.

For a common equilibrated hidden bath with zero chemical potential, the
leading high-temperature term is T_h^2 Tr(M_b^2)/24. The coefficient of
phi squared in the trace is 3 lambda_phi+2 g_X^2+h_D^2 (two real X modes).
Thus

\[
\Pi_\phi/T_h^2=\lambda_\phi/4+g_X^2/6+h_D^2/12=.004175,
\quad T_{\rm curvature\ zero}/v=\sqrt{.01/.004175}=1.547646465.
\]

This checks the review's leading estimate. It is not a resummed transition
temperature, a nonzero-charge phase diagram, or a relation between hidden
and visible temperatures. The original cold rotating condensate is not
that thermal ensemble. Its visible temperature alone does not restore a
decoupled hidden symmetry.

As a dimensional warning for a **persistent horizon-scale scaling network**,
G sigma/H0 of order DeltaT/T gives sigma of order 2.23e-9 GeV^3 for
H0=70 km/s/Mpc and DeltaT/T=1e-5. The benchmark planar tension is about
4.24e16 times that estimate. The unknown coefficient and network history
are not numerical error bars. No C_l spectrum, transient-bubble signal,
Planck likelihood, or new exclusion follows from this estimate alone.
