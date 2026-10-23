# Equations, conventions and domains

Natural units unless a conversion is displayed. All limits below are
conditional on fixed profiles, halo, classical tree interactions and
the specified detector releases. Derivation methods have established
precedents; see [SOURCES.md](SOURCES.md).

## 1. Action and inputs

Use Einstein-Hilbert gravity, ordinary SM matter and canonical scalar
kinetics, with the parent's potential

\[
U_d=\frac{\lambda}{4}(\phi^2-v^2)^2+g^2\phi^2|X|^2
+\lambda_X|X|^4+\frac12(m_{D0}^2+h_D^2\phi^2)D^2
+\frac{\lambda_D}{4}D^4.
\]

Add the proposed subtracted portal
\(\Delta U=\kappa(\phi^2-v^2)(H^\dagger H-v_H^2/2)/2\).
This is an exploratory extension of the original zero-portal matching.
Inputs: g=.1, lambda=.01, lambda_X=.001, h_D=.01,
lambda_D=.001, m_D0=.1v. Profiles are the parent's 801-point solutions
for Q=19972,39944,79988 on rho=gvr in [0,60]. Their dimensionless
energies and shapes are held fixed while v varies. v=.5,.75,1,1.25 GeV
and kappa are supplied, not derived from electroweak physics.

SM reference values: m_h=125 GeV, v_H=246 GeV, m_N=.939 GeV,
f_N=.308; reduced Planck mass M=2.435e18 GeV. Adopted constants and
conversion factors are explicit in `model.py` and `history.py`.
The total renormalized constant vacuum term for the cosmological
history is supplied separately. No term in this calculation selects it.

At kappa=0, s=phi-v has

\[
m_s^2=2\lambda v^2,\quad m_X=gv,\quad
m_D=\sqrt{g^2+h_D^2}\,v.
\]

All physical s masses are below 2m_mu, 2m_X and 2m_D. They have no
tree hidden two-body decay in these examples. Finite-kappa profile
shifts, loop running, induced H-X/H-D portals and a full parameter
renormalization calculation are not performed.

## 2. Object recoil amplitude

The fixed-charge energy identity is

\[
E_Q=\lambda v^2\int(v^2-\phi_Q^2)d^3x.
\]

The h source from the portal is kappa v_H(phi_Q^2-v^2)/2. Thus

\[
g_{Qh}(q)=-\frac{\kappa v_HE_Q}{2\lambda v^2}F_Q(q),\qquad
F_Q(q)=\frac{\int_0^\infty r^2(v^2-\phi_Q^2)j_0(qr)dr}
{\int_0^\infty r^2(v^2-\phi_Q^2)dr}.
\]

The nucleon coupling is f_N m_N/v_H. For a Born, nonrelativistic,
spin-independent interaction define

\[
\sigma_{QN}(0)=\frac{\mu_{QN}^2}{\pi}
\left(\frac{\kappa f_Nm_NE_Q}{2\lambda v^2m_h^2}\right)^2.
\]

For isotope A, sigma_QA(0)=sigma_QN A^2(mu_QA/mu_QN)^2, and

\[
\frac{d\sigma_A}{dE_R}=
\frac{m_A\sigma_{QA}(0)}{2\mu_{QA}^2u^2}
F_A^2(q)F_Q^2(q)\left(\frac{m_h^2}{m_h^2+q^2}\right)^2,
\quad q^2=2m_AE_R.
\]

Integrate incident flux (rho_DM/E_Q)u against a truncated shifted
Maxwellian. Inputs: rho_DM=.3 GeV/cm^3, v0=220, v_E=232,
v_esc=544 km/s; all-DM object fraction initially 1. The code retains
unit conversions explicitly and returns events/(keV tonne year).
At fixed halo shape all signal strengths scale as f_DM kappa^2 and
detector coupling limits as f_DM^(-1/2). A different fraction also
changes the supplied cosmic population; it is not derived here.

Natural xenon isotope number fractions/atomic masses come from NIST.
Nuclear masses subtract 54 electron masses from atomic masses. Use
the Helm convention c=1.23 A^(1/3)-.6 fm, a=.52 fm, skin=.9 fm,
r1^2=c^2+7pi^2 a^2/3-5skin^2 and
F_A=3j1(qr1)/(qr1) exp[-(q skin)^2/2]. The parent used a simpler
r1^2=(1.2 A^(1/3))^2-5skin^2 and mono-131 xenon. The reference
1-50 keV ideal count changes by -.76% with the present Helm choice
and +.29% with natural isotopes, when each change is isolated.
These are declared convention changes, not a refit.

## 3. Detector likelihood and approximations

True bin i has rate integral R_i. For exposure T and public migration
probability M_ij, the reconstructed signal is

\[
S_j=T\sum_i R_i M_{ij}.
\]

M includes selection efficiency and energy resolution, **not**
exposure. T=1 tonne year for XENON1T and 20 tonne years for the named
XENONnT projection. Use 8-point Gauss-Legendre quadrature in every
true bin; the 16-point control is archived.

For a public likelihood toy k and signal strength alpha, construct

\[
t_k(\alpha)=\sum_j[-2\ln L_{kj}(\alpha S_j)]
-\min_{\alpha\geq0}\sum_j[-2\ln L_{kj}(\alpha S_j)].
\]

Linear interpolation of the supplied signal axes and a finite strength
grid implement the release's prescription. Observed limits use the
release's recommended signal-dependent threshold. The projection
uses chi-square_1 inverse-survival(.1)=2.70554, as in its supplied
method, and all 1000 background toys. Alpha=(kappa/1e-9)^2. The
observed control uses 1000 and 4000 strength samples; reported observed
coupling limits use 4000. Numerical agreement is not the statistical
accuracy of the approximate experimental release: its authors find
typical heavy-WIMP biases around 10%, with larger differences near
threshold. No Gaussian discovery significance is inferred here.

The first LZ estimate maps the published point-WIMP limit at mass E_Q
using the ratio

\[
r=\frac{\int\epsilon_{LZ}(E_R)\,dR_{object}/dE_R\,dE_R}
{\int\epsilon_{LZ}(E_R)\,dR_{point,same\ mass}/dE_R\,dE_R},
\quad
\kappa_{90}^{est}=10^{-9}
\sqrt{\frac{\sigma_{90}^{LZ}(E_Q)}{r\sigma_{QN}(10^{-9})}}.
\]

The alternative fits a free point-WIMP mass to the normalized,
efficiency-weighted true spectrum by minimizing squared Hellinger
distance, one-half the integral of (sqrt(p)-sqrt(q))^2. This metric
remains finite at different exact recoil endpoints; a true-spectrum
KL divergence would not. Then map its normalization through the
published limit at that fitted mass. Neither substitutes for the LZ
spatial/response/background likelihood. Their 22% coupling difference is material to the proposed
overlap. We do not assert these dated data are the latest full set of
constraints or that a globally allowed interval has been found.

Separately, fit a point WIMP directly in the historical XENONnT
reconstructed response. If p_j and q_j are normalized signal shapes,

\[
D_{KL}(p\Vert q)=\sum_{p_j>0}p_j\ln(p_j/q_j).
\]

With a free matching normalization and N pure signals the expected
twice-log likelihood ratio is 2N D_KL. The reported N=9/(2D_KL) is
only a signal-only heuristic information requirement; no real detector
background or nuisance test is implied.

## 4. Same-portal mixing, decays and low-energy production

Let delta= kappa v v_H. To first order the light eigenstate contains
-theta h_SM, with theta=delta/(m_h^2-m_s^2). The original terms
kappa v_H h_SM s^2/2 and lambda v s^3 give

\[
g_{hss}=\kappa v_H+6\lambda v\theta
=\kappa v_H\frac{m_h^2+2m_s^2}{m_h^2-m_s^2},
\quad
\Gamma_{h\to ss}=\frac{g_{hss}^2}{32\pi m_h}
\sqrt{1-4m_s^2/m_h^2}.
\]

For a virtual scalar of mass m=sqrt(s), define a unit-SM-Higgs
spectral response Gamma_SM(m). Lepton contributions are

\[
\Gamma_{\ell\ell}^{SM}(m)=\frac{m_\ell^2m}{8\pi v_H^2}
(1-4m_\ell^2/m^2)^{3/2}.
\]

Use Gamma_gammagamma=C_gamma^2 alpha_EM^2 m^3/(256pi^3 v_H^2).
C_gamma=10 nominal and 100 in the enlarged envelope are stated
uncertainty choices, not exact low-energy hadronic loop predictions.
The leading chiral two-pion response is

\[
\Gamma_{\pi\pi}^{SM}(m)=\frac{3}{32\pi v_H^2m}
\left(\frac{2m^2+11m_\pi^2}{9}\right)^2
\sqrt{1-4m_\pi^2/m^2}.
\]

Use m_pi=.1349768 GeV for this approximate three-pion-state
counting; exact charged/neutral and higher-order QCD corrections are
not calculated. The enlarged response multiplies this width by 1000
and caps the assumed response at Gamma_SM=m for m>=1 GeV. It is
a broad production stress test, not a rigorous QCD error bar.

The physical s width used is theta^2 Gamma_SM(m_s). Here only ee
and gamma gamma are relevant, since m_s<2m_mu and hidden pairs are
closed. Neutrino mass-suppressed decays are negligible at this accuracy.

Inclusive SM pair annihilation has a scalar production spectral factor

\[
P(s)=\frac{\Gamma_{h^*\to ss}^{eff}(s)}{(s-m_h^2)^2}
+\frac{\delta^2[\Gamma_{s^*\to XX^*}(s)+\Gamma_{s^*\to DD}(s)]}
{(s-m_h^2)^2(s-m_s^2)^2},
\]

where

\[
\Gamma_{h^*\to ss}^{eff}=\frac{\kappa^2v_H^2}{32\pi\sqrt{s}}
\left(\frac{s+2m_s^2}{s-m_s^2}\right)^2\sqrt{1-4m_s^2/s},
\]
\[
\Gamma_{s^*\to XX^*}=\frac{g^4v^2}{4\pi\sqrt{s}}\sqrt{1-4m_X^2/s},
\quad
\Gamma_{s^*\to DD}=\frac{h_D^4v^2}{8\pi\sqrt{s}}\sqrt{1-4m_D^2/s}.
\]

All expressions vanish below the relevant threshold. The first term
includes the off-shell s exchange, in addition to the portal contact
vertex. At the cold temperatures neither scalar pole is in the pair
production domain that contributes significantly.

For Maxwell-Boltzmann equilibrium SM initial particles the reaction
number density and injected power are

\[
\gamma=\frac{T}{\pi^3}\int dm\,m^4\Gamma_{SM}(m)P(m^2)K_1(m/T),
\quad
\mathcal Q_{pair}=\frac{T}{\pi^3}\int dm\,m^5\Gamma_{SM}(m)P(m^2)K_2(m/T).
\]

An independent spin-averaged ee->XX* calculation gives
sigma= kappa^2 m_e^2 g^4 v^4 beta_e beta_X /
[8pi(s-m_h^2)^2(s-m_s^2)^2], agreeing with
4pi Gamma_ee^SM P_XX/beta_e^2. This checks spectral/phase-space
normalization through a separate route.

Single s inverse decays contribute

\[
\mathcal Q_{inv}=\frac{\Gamma_s m_s^3T}{2\pi^2}K_2(m_s/T),
\quad
\dot n_{inv}=\frac{\Gamma_s m_s^2T}{2\pi^2}K_1(m_s/T).
\]

## 5. Cold-history energy and scattering estimates

Set T_i=.020 GeV, T=T_i exp(-N), g_*=g_*s=10.75,
rho=pi^2 g_* T^4/30, entropy density s_ent=2pi^2 g_*s T^3/45,
H=sqrt(rho/(3M^2)). Stop at 10 keV; production has become
exponentially irrelevant well before electron annihilation, so its
late entropy degrees of freedom do not change this production integral
appreciably. Fixed degrees of freedom and MB statistics are approximations.

For nonnegative particle pressure, no new source and no vacuum
instability, treating all injected energy as immediately pressureless
and stable overestimates its subsequent density:

\[
B(N)=\int_0^N\frac{\mathcal Q_{inv}+\mathcal Q_{pair}}
{s_{ent}H}dN',\qquad \rho_{new}\leq B(N)s_{ent}(N).
\]

Compare B_final to the supplied rho_DM/s_ent=4.4e-10 GeV.
Variance estimates for weak quasiparticles give a conservative-order
hidden tadpole diagnostic 1.5 rho_new/(lambda v^4), and a visible
trace-source bound 2kappa rho_vis/(m_h^2 lambda v^2). These are small
in the computed cold cases; nonlinear formation is not inferred.

For included coherent elastic collisions use |F_Q|<=1 and speed<=1.
For a lepton of momentum p the heavy-object contact estimate obeys

\[
\sigma_{Q\ell}\leq\frac{g_{Qh}^2(m_\ell/v_H)^2}{\pi m_h^4}
(m_\ell^2+p^2/2).
\]

Integrate its MB moments for e and mu, and the stated pion and
photon scalar-current envelopes, plus a generous baryon coherence
factor 64 with n_B/s_ent=8.6e-11. Tau=sum integral(n_i sigma_i/H)dN.
The JSON key `visible_any_collision_optical_depth_BOUND` means any
collision in these included **coherent elastic** channels. It does not
bound uncomputed inelastic matrix elements or all possible damage.

For nuclear injection set Y_phi=B_final/m_s, allowing *all* injected
energy, including non-scalar production, to become scalar rest mass.
ACROPOLIS starts at T0=10 keV with n_phi/n_gamma=7.04 Y_phi; this is
the post-electron-annihilation entropy/photon conversion. Use physical
ee/gamma branching ratios for the two main response choices. A separate
stress case pairs the enlarged population with a pure-photon decay at
the electron-only longest lifetime. Those independent extrema are an
additional test, not a physical branching-ratio point or proof of the
maximum effect at every intermediate lifetime.

Run the original cascade and 9-species nuclear matrix, including its
three supplied BBN abundance columns. Relative changes in the report
refer to the central column after the package's decay squeezing. Its
background temperature coverage is checked before running. Effective
resolution is 120 energy / 20 temperature points per decade; the final
stress case is rerun at 150 / 30. This tests incremental disintegration,
not a refit of primary BBN or CMB/structure likelihoods.

## 6. Hot energy transfer and the missing chemical dynamics

For an SM bath starting at 80 GeV with g_vis=106.75 and a hidden
sector initially empty, leading zero-temperature Higgs decay injects

\[
\mathcal Q=\Gamma_{hss}m_hn_h^{MB},\quad
n_h^{MB}=\frac{m_h^2T}{2\pi^2}K_2(m_h/T).
\]

Let V=rho_vis a^4/rho_i and R=rho_hid a^4/rho_i, with a_i=1.
Then dV/dN=-a^4 Q/(H rho_i), dR/dN=+a^4 Q/(H rho_i), and
H^2=(rho_vis+rho_hid)/(3M^2). Thus V+R=1 exactly for this
relativistic-daughter approximation. Thermal SM mass corrections and
changing g_vis are not included, and an earlier hotter epoch is not
being modeled.

Two different diagnostics:

* **Rapid chemical equilibrium assumed.** T_h=(30rho_hid/(pi^2*4))^(1/4),
  Pi=(lambda/4+g^2/6+h_D^2/12)T_h^2, with a leading local
  restoration criterion Pi>lambda v^2.
* **Collisionless, relativistic s daughters.** Their second moment obeys
  d< s^2 >/dt+2H< s^2 >=2Gamma_hss m_h T K0(m_h/T)/pi^2.
  The corresponding self-curvature diagnostic is 3< s^2 >/v^2.

The latter source follows by boosting isotropic two-body daughters:
the summed inverse-energy angular average is
2 ln[(E+p)/(E-p)]/p. With p=m_h sinh(y), its MB momentum integral
reduces to 2T K0(m_h/T). Direct quadrature independently confirms it.

The two diagnostics are not a complete hidden Boltzmann system. The
rapid-equilibrium temperature cannot be justified by an elastic
scattering rate alone: number-changing and species conversion rates
must also be solved. Conversely a collisionless diagnostic cannot
describe a sector that actually equilibrates. Their difference leaves
the actual hot history open. A phase criterion alone also does not
simulate wall formation, percolation or survival of pre-existing objects.
