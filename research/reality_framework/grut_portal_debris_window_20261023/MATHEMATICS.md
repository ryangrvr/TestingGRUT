# Widths, charge accounting and the conditional coupling window

Natural units, except displayed lifetimes. These are leading classical/tree
calculations about the selected +v vacuum, with the proposed Higgs portal.
They are not an in-medium formation calculation or an all-orders matching.

## Physical masses and visible width

Set s=phi-v. With g=.1, lambda_X=.001, h_D=.01 and m_D0=.1v,

\[
m_X=gv,\qquad m_s=\sqrt{2\lambda_\phi}\,v,\qquad
m_D=\sqrt{g^2+h_D^2}\,v.
\]

The scalar mass cannot be scanned independently of lambda_phi and v.
The empirical extension is

\[
\Delta U={\kappa\over2}(\phi^2-v^2)(H^\dagger H-v_H^2/2).
\]

Its vacuum off-diagonal mass-squared is kappa v v_H. At small mixing,

\[
\theta={\kappa v v_H\over m_h^2-m_s^2},\qquad
\Gamma_{\rm vis}=\theta^2\Gamma_{\rm SM,unit}(m_s).
\]

Use m_h=125 GeV, v_H=246 GeV and hbar=6.582119569e-25 GeV s.
For each open charged-lepton channel,

\[
\Gamma_{\ell\ell}^{\rm unit}=
{m_\ell^2m_s\over8\pi v_H^2}
\left(1-{4m_\ell^2\over m_s^2}\right)^{3/2}.
\]

The pi-pi and KK inputs are the published hipsofcobra 0.0.8 calculation,
with both DGL and BTPZ matching and its tabulated bands. The package's
isospin-limit thresholds .268 and .994 GeV are retained. Above .5 GeV,
the old FORESEE total minus its listed pi/K/e/mu contributions supplies
a phenomenological additional-channel remainder. The added photon width is
C_gamma^2 alpha^2 m_s^3/(256 pi^3 v_H^2), with nominal C_gamma=10,
0/10 envelope and C_gamma=100 stress at the original mass point.
These photon choices and the combined method envelope are not statistical
confidence bands or a complete precision calculation of the total width.

## Existing hidden channels

The sXX* and sDD vertices are 2g^2v and 2h_D^2v. The constant-amplitude
two-body formula is |M|^2 beta/(16 pi m_s), with an additional factor
1/2 for identical daughters. Consequently,

\[
\Gamma_X={g^4v^2\over4\pi m_s}
\sqrt{1-{4m_X^2\over m_s^2}},\qquad
\Gamma_D={h_D^4v^2\over8\pi m_s}
\sqrt{1-{4m_D^2\over m_s^2}},
\]

where a closed channel has zero width. The threshold evaluations use
dimensionless coupling ratios to avoid opening an exact boundary through
floating-point roundoff.

The homogeneous binding condition is lambda_phi lambda_X<g^4. Therefore

\[
2g^2<\lambda_\phi<{g^4\over\lambda_X}
\]

is a possible simultaneous X-decay/bulk-binding interval. It exists iff
lambda_X<g^2/2. For the fixed couplings it is .02<lambda_phi<.1;
DD opens at .0202. Bulk binding is not finite-Q existence or stability.
Actual profiles are re-solved rather than assumed from this interval.

For a nominal prompt cut tau_cut,

\[
\kappa_{\rm prompt}=
\sqrt{\frac{\max(\hbar/\tau_{cut}-\Gamma_X-\Gamma_D,0)}
{(v v_H/(m_h^2-m_s^2))^2\Gamma_{\rm SM,unit}}}.
\]

With no hidden decay, tau_s kappa^2 is constant while recoil counts
scale as kappa^2. Their product is independent of kappa for a fixed
object/population/response. Changing kappa alone cannot repair a large
lifetime/recoil conflict. A one-second cut is not a substitute for BBN
or deposition calculations. Vacuum widths describe escaped particles,
not necessarily the excitation spectrum inside a dense object.

## Conserved-charge normalization with stable daughters

The conditional population is initially one sign of charge, all retained
net charge eventually resides in identical test objects, and each object
releases Delta=Qm_X-E_Q as initially slow free s particles. Then

\[
N_s={\Delta\over m_s},\quad N_X=N_{X^*}=B_XN_s,\quad
N_D=2B_DN_s,\quad
M_{late}=E_Q+2N_s(B_Xm_X+B_Dm_D).
\]

The pair daughters add zero net charge. Their rest mass is counted as
late matter. Their kinetic energy redshifts under the stated collisionless
assumption; it is not silently treated as rest mass or dark radiation.
For an initial scalar gamma_i, N_s=Delta/(gamma_i m_s); the reported
witnesses use gamma_i=1. No emission spectrum is predicted here.

With B_DM=rho_DM/s_vis=4.4e-10 GeV supplied,

\[
Y_Q^{required}={Q B_{DM}\over M_{late}},\qquad
f_{object}={E_Q\over M_{late}},\qquad
B_{visible}={Y_Q\over Q}\Delta B_{vis}.
\]

The main scan reports both this required Y_Q and the late matter obtained
at the original supplied Y_Q=4.8294009076988175e-9. These are distinct
questions. The two additional witnesses fix that original Y_Q and solve
for v from Y_Q M_late(v)/Q=B_DM. This is a density calibration of v,
not its prediction. Assuming a common halo distribution, f_object
multiplies the object's recoil spectrum. No inverse reaction, free-particle
accretion, annihilation, transport, population collapse or kinetic
decoupling calculation is contained in this rest-mass accounting.

## Finite-Q response and the statistical comparison

The stationary normalization is rho=gv r, F=f/v, S=phi/v and Omega=omega/(gv).
The solver uses

\[
F''+2F'/\rho=(S^2+(\lambda_X/g^2)F^2-\Omega^2)F,
\quad S''+2S'/\rho=[(\lambda_\phi/g^2)(S^2-1)+F^2]S,
\]
\[
g^2Q=4\pi\Omega\int\rho^2F^2d\rho.
\]

Regularity at the origin and F(infinity)=0, S(infinity)=1 are imposed.
The classical identity E_Q=lambda_phi v^2 integral(v^2-phi^2)d^3x
is checked directly from raw profiles. It gives the leading nucleon
coefficient

\[
c_n={\kappa f_Nm_NE_Q\over2\lambda_\phi v^2m_h^2},\qquad
\sigma_n={\mu_n^2c_n^2\over\pi},\qquad
F_Q(q)={\int\rho^2(1-S^2)\,\mathrm{sinc}(q\rho/(gv))d\rho
\over\int\rho^2(1-S^2)d\rho}.
\]

The nuclear recoil computation retains the parent's natural-xenon isotope
mixture, Helm factor, halo and official historical response arrays.
The Born row-sum diagnostics remain much smaller than one for the selected
portal-sized points. They do not constrain uncomputed inelastic processes.

Two LZ mappings are explicitly approximate: equal efficiency-weighted
counts relative to a same/heavy-mass point WIMP, and a normalized-shape
matched point WIMP. The latter minimizes squared Hellinger distance and
maps its normalization to the published point-particle limit. Neither is
an LZ profile likelihood. Reference masses outside the table are not
extrapolated. Masses above the heavy table use its last mass for the
count comparison and disclose that reference; soft fits below 9 GeV
need the separate low-mass analysis rather than that table.

The XENONnT calculation reruns all 1000 public background-only toys of
the historical 20 tonne-year approximate likelihood. A median projected
90% exclusion reach is not discovery significance. A provisional overlap
requires max(kappa_prompt,kappa_median_reach) below both approximate LZ
screens. This is an exploratory screen on already-known releases, not a
certified currently allowed interval or a held-out empirical test.

## New-parameter linear onset is still not a final charge

On the dilute positive-phi branch a=g^4/lambda_phi-lambda_X>0,
f_0^2 approximately equals n_Q/(gv), k_peak^2=a f_0^2 and
Gamma_peak approximately equals a f_0^2/(2gv). At the new lambda=.03
fixed-yield witness and 20 MeV visible bath, Gamma/H=1.3955e8.
The wavelength-cube charge is 1.4484e11. Both numbers are diagnostics of
a supplied coherent state. They do not produce Q=200000 or a distribution.
The controlled dilute equation must be matched to full-field core dynamics
before assigning a nonlinear outcome, collision kernel or detector population.
