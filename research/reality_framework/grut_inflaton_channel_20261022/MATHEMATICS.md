# Derivative decay and its conserved receiver

22 October 2026. Continuation of d601fd1. Natural units, reduced Planck
mass M, signature (-,+,+,+). This computes a conditional limit of the
Camargo-Molina–Markkanen–Scott action (1905.00045v2), not a new GRUT law.

## 1. Start from the displayed metric

The literal Einstein action (5.2) gives

\[
G_{\varphi\varphi}=\kappa={3c^2\over2},\quad
G_{\sigma\sigma}=\kappa+e^{-\varphi/M}P(\sigma),\quad
P={b^2e^{-c\sigma/M}\over(e^{-\sigma/M}-1)^2}.
\]

The potential separates the inflaton term
alpha² M⁴(1-exp(-varphi/M))² from the sigma/Higgs terms.
Work near varphi=0 and a positive, small, slowly changing sigma_bar,
after the pole term dominates. No crossing of sigma=0 is used.
Let I=sqrt(kappa) varphi and locally normalize the remnant fluctuation
chi=sqrt(G_sigma_sigma(0,sigma_bar)) delta_sigma. Define

\[
F={P(\bar\sigma)\over\kappa+P(\bar\sigma)}.
\]

At first order in I,

\[
\mathcal L_{\rm int}={F\over2\sqrt\kappa M}
 I(\partial\chi)^2,\qquad m_I^2={2\alpha^2M^2\over\kappa}.
\]

The kinetic-pole limit is F→1. Using the log coordinate
chi=bM ln(sigma/sigma_ref) in that limit gives the whole kinetic
factor B(I)=exp[-I/(sqrt(kappa)M)]. The regular metric contribution,
finite sigma/M and the remnant potential give small corrections only
within the stated background regime. The leading width does not depend
on b or on choosing sigma exactly at its tiny selected minimum.

For two identical daughters with mass m_chi, the derivative amplitude is

\[
|\mathcal A|={F\over2\sqrt\kappa M}(m_I^2-2m_\chi^2),
\]
\[
\Gamma_{I\to\chi\chi}={F^2m_I^3\over128\pi\kappa M^2}
 (1-2r_\chi)^2\sqrt{1-4r_\chi},\qquad
r_\chi=m_\chi^2/m_I^2.
\]

The calculations use light daughters, r_chi≈0. At c=1 this is
m_I³/(192 pi M²), the familiar single-real-scalar scalaron coefficient
in Gorbunov–Panin (1009.2448), independently derived here for this metric.
It does not imply their universal visible-sector couplings are present
in the action under study.

**Coordinate check.** In the pole limit let g=1/(sqrt(kappa)M).
The nonzero connection at the expansion point has Gamma^I_chi_chi=g/2.
Normal coordinates give I_old=I_RNC-g chi_RNC²/4 and remove the cubic
kinetic vertex. The inflaton potential then supplies a cubic amplitude
g m_I²/2 in the massless limit, identical to the derivative amplitude.
A field redefinition does not remove the physical decay.

The supplied N=55 single-adiabatic inflation normalization follows
d601fd1. Recomputed here, c=1 gives alpha=1.05623e-5,
m_I=2.96981e13 GeV and Gamma_chi=7.32378 GeV. Other c values are recorded,
with alpha recalibrated consistently for each. No dark-energy density
was fitted to obtain these widths. The common potential scale, assumed
inflationary state and absolute final vacuum convention are inherited.

There are additional common-action consistency relations. If the tiny
prepared displacement of d601fd1 is used, rho_star=c² V0²/(alpha² M⁴),
and therefore

\[
\boxed{\rho_* m_I^2={4V_0^2\over3M^2}}.
\]

This links the putative prepared remnant to the inflaton mass without
adjusting c and alpha separately. It does not select that preparation
or identify rho_star with today's density. Under the parent's
single-adiabatic first-order inflation assumptions, define
D=8(1-n_s)-r and E=4(1-n_s)-r. The same model gives

\[
\kappa={8r\over E^2},\quad {m_I^2\over M^2}={3\pi^2 A_s\over32}D^2,
\]
\[
\Gamma_\chi=F^2\left({3\pi^2 A_s\over32}\right)^{3/2}
 {M D^3 E^2\over1024\pi r}.
\]

Thus the partial rate can be expressed through the same inflation
observables instead of a separately fitted c. The stated branch has
0<r<4(1-n_s), the effective theory must remain perturbative, and
entropy-to-curvature transfer or changed potential normalizations would
require a new relation. These algebraic consequences of borrowed
physics are not new independent measurements or a GRUT uniqueness test.

## 2. An independent quantum-wave normalization check

In flat spacetime, prepare a smooth external inflaton pulse
z(t)=g I(t)=q exp[-t²/(2T²)] cos(t), using m_I=1 units. The asymptotic
metric tends to one at both endpoints, so in/out particle number is
unambiguous in this test. For v_k=sqrt(B) chi_k,

\[
\ddot v_k+[p^2+\ddot z/2-\dot z^2/4]v_k=0,
\quad p=k/m_I.
\]

The code integrates positive-frequency incoming modes and extracts
the late Bogoliubov beta. A separate selected-mode integration uses
the original equation chi_ddot-z_dot chi_dot+p² chi=0. The canonical
Wronskian is checked. The produced radiation energy is

\[
\rho_\chi/m_I^4={1\over2\pi^2}\int_0^\infty dp\,p^3|\beta_p|^2.
\]

At first Born order,

\[
|\beta_p|^2=q^2p^2{\pi T^2\over2}
 [e^{-T^2(2p-1)^2/2}+e^{-T^2(2p+1)^2/2}]^2.
\]

For T>>1 the integrated answer, up to exponentially small carrier-tail
terms, is

\[
\rho_{\chi,\rm Born}/m_I^4=
{q^2T\over256\sqrt\pi}
 [1+5/T^2+15/(4T^4)].
\]

The leading coefficient equals Gamma_chi times the cycle-averaged
inflaton pulse energy integrated over time. The explicit finite-envelope
terms prevent mistaking a finite-duration effect for a decay-rate error.
For q=.001/.002 and T=40, direct wave energies differ from Born by
.02225%/.08901%; smaller q approaches the perturbative limit.
The spectra peak near p=1/2, as expected for two-body decay. Momentum
grid, time-step/tolerance refinement and the original-coordinate
calculation are recorded.

This pulse is an **external-pump scattering check**. It is not an
energy-conserving backreacting reheating simulation, and its finite
duration is not assigned to a cosmological history. The conserved
coarse-grained calculation below is a separate perturbative limit.
Combining their evidence checks the source coefficient; it does not
solve early nonperturbative preheating or the full quantum state.

## 3. Energy and momentum accounting

Before averaging, a homogeneous I and remnant waves have, neglecting
their tiny potential during this stage,

\[
\rho_\chi={B\over2}(\dot\chi^2+a^{-2}|\nabla\chi|^2),\quad
p_\chi={B\over2}(\dot\chi^2-a^{-2}|\nabla\chi|^2/3),
\]
\[
\dot\rho_\chi+3H(\rho_\chi+p_\chi)
 ={\dot B\over2}(a^{-2}|\nabla\chi|^2-\dot\chi^2).
\]

The inflaton equation has the opposite energy transfer. This follows
from the two-field action; energy is not removed from an evolving
constant without a receiver.

In the weak-amplitude, incoherent perturbative limit m_I>>H,Gamma,
average the quadratic inflaton oscillations to dust. For a total rate
Gamma_tot=Gamma_chi+Gamma_vis the receiver equations are

\[
\dot\rho_I+3H\rho_I=-\Gamma_{\rm tot}\rho_I,
\]
\[
\dot\rho_d+4H\rho_d=\Gamma_\chi\rho_I,\qquad
\dot\rho_v+4H\rho_v=\Gamma_{\rm vis}\rho_I,
\]
\[
3M^2H^2=\rho_I+\rho_d+\rho_v.
\]

Here subscripts d and v mean dark and visible **radiation**, not dust
and vacuum. Only Gamma_chi is derived from the displayed terms in this
packet. The Gamma_vis comparisons are explicit diagnostic inputs;
they do not constitute an action for visible-sector reheating.

The main solver uses tau=Gamma_tot t and densities scaled by
3M² Gamma_tot². It evolves H independently with the acceleration
equation H_tau=-3H²/2-(r_d+r_v)/2. Its Friedmann residual obeys
C_tau=-3H C, and vanishes on physical solutions. Both sources sum to
the inflaton loss. An independent formulation evolves Omega_I,
log(H/Gamma_tot) and tau in N:

\[
\Omega_I'=\Omega_I[1-\Omega_I-\Gamma_{\rm tot}/H],\quad
(\ln H)'=-2+\Omega_I/2,\quad\tau'=\Gamma_{\rm tot}/H.
\]

The declared start H/Gamma_tot=1e6 is already well within the weak
oscillation regime for the recorded rates. Earlier inflation and
nonperturbative preheating are omitted. Initially only the inflaton
condensate is populated. The selected spectator potential is negligible;
a different large remnant population would require solving it too.

The two implementations agree on the final H to about 2.2e-12 relative
in the refined comparison. Inflaton/radiation equality occurs near
Gamma_tot t=1.05380. The calculation produces radiation domination
for the computed-channel-only case. It does not produce a visible
plasma, a thermal spectrum, baryogenesis or the selected homogeneous
dark-energy value by that fact alone.

## 4. A linked radiation relation

If a visible bath is supplied and thermalizes, while the chi quanta
remain decoupled and relativistic, the final high-temperature energy
ratio is rho_d/rho_v=Gamma_chi/Gamma_vis for these constant-width
zero-initial-radiation histories. Standard visible entropy conservation
then gives

\[
\Delta N_{\rm eff}={4\over7}g_{*\rho,\rm rh}
 \left({10.75\over g_{*s,\rm rh}}\right)^{4/3}
 {\Gamma_\chi\over\Gamma_{\rm vis}}.
\]

For gstar_rho=gstar_s=106.75 the coefficient is 2.85792. This is the
standard decoupled-radiation result, e.g. Cicoli–Conlon–Quevedo
(1208.3562), not a new GRUT relation. The examples Gamma_vis/Gamma_chi
=1,4,20 give Delta_Neff=2.858,.7145,.1429. If Gamma_vis=0 there is
no visible reference bath and Delta_Neff is undefined, not a finite
cosmological prediction. Late visible entropy injection, appreciable
initial radiation, nonrelativistic daughters or additional channels
change this formula.

For an **illustrative tolerance**, Delta_Neff<.2 would require
Gamma_vis/Gamma_chi>14.29. The .2 is a declared diagnostic input,
not an observational limit extracted here. No current-world bound or
global likelihood is claimed. A specified visible interaction is
needed before this ratio becomes a prediction of a complete model.

## 5. Why a radiation receiver does not select the homogeneous value

In the post-pole kinetic limit there is a shift symmetry in chi. Its
homogeneous equation is

\[
{d\over dt}(a^3 B\dot\chi_0)=-a^3W_{,\chi}.
\]

At negligible W_chi, a patch with chi_dot=0 keeps **any** initial chi.
Changing the inflaton metric, damping momentum and producing
relativistic waves do not provide an attracting field position.
The numerical control checks constant stationary positions and the
conserved canonical momentum in a varying B and expanding background.
This statement applies to the stated shift-symmetric stage. It does
not exclude potential-driven relaxation earlier, large fluctuations,
nonperturbative dynamics, or new interactions that break that symmetry.

As a useful conditional input check, if a homogeneous small sigma
remains frozen with energy alpha² M² sigma² until H≈Gamma_chi, avoiding
its domination there requires approximately

\[
\sigma/M\ll\sqrt3\,\Gamma_\chi/(\alpha M).
\]

The c=1 scale is 4.93e-13. This is a frozen-state diagnostic, not a
general fluctuation bound: if the field rolls earlier, its energy
must be evolved. The previous c=10 fixed-inflaton stochastic probe
is not a post-inflation state at this c=1 point and is not inserted
into these receiver equations. A genuine joint-distribution transport
through the earlier changing metric is still required.
