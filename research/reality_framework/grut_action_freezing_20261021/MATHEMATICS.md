# Declared actions, equations and bounds

Natural units, signature (-,+,+,+), reduced M. These proposed extensions are
borrowed kinetic-pole constructions, not derived from GRUT closure. S1--S3
are identified in SOURCES.md. The original q action at 10e66b7 is preserved.

## Action and normalization

Let a=sqrt(2/3), V_H=lambda(h^2-v^2)^2/4, m^2=4 alpha^2 M^2/3 and

\[
S=\int\sqrt{-g}\left[\frac{M^2R}{2}
-\frac12G_{AB}\partial q^A\partial q^B-U\right]+S_r,
\]
\[
U=\frac34m^2M^2(1-e^{-a\phi/M})^2+e^{-2ac\phi/M}V_H(h).
\]

At c=1 the potential is S1 (2.3). Its combination with the kinetic
illustrations below is our positive hybrid completion, **not the literal full
action (3.6)**. The c-dependent positive exponential continuation is an explicit
choice. The canonical Higgs kinetic term is included. In the phi>0 chart,

\[
ds^2_{KP0}=d\phi^2+dh^2+b^2M^2(h/v)^2d\ln(\phi/M)^2,
\]
\[
ds^2_{KP1}=d\phi^2+dh^2+b^2M^2d[(h/v)\ln(\phi/M)]^2.
\]

Both metrics are positive definite at finite positive phi. The absolute
reference U(0,v)=0 is supplied. Adding these kinetic terms changes the model;
they do not follow from the old linear K/quadratic epsilon q action. These are
QV3-KP0 and QV3-KP1 proposals, not approved new empirical targets.

At h=0, the minimum for c=1 obeys exactly

\[
\frac{a\phi_0}{M}=\ln\left(1+\frac{V_0}{\alpha^2M^4}\right),
\qquad V_0=\lambda v^4/4.
\]

For tiny V0/(alpha^2 M^4) and general c,

\[
\phi_0\simeq\sqrt{3/2}\frac{cV_0}{\alpha^2M^3},\qquad
\rho_*=U(\phi_0,v)\simeq\frac{c^2V_0^2}{\alpha^2M^4}
=\frac{c^2\lambda^2v^8}{16\alpha^2M^4}.
\]

phi0 has dimension one and rho* four. With M=2.435e18 GeV,v=246 GeV,
lambda=.129,alpha=9.97e-6, m=2.80326e13 GeV, phi0/M=4.13933e-56 and
rho*=3.99167e-48 GeV^4 at c=1. c=2.5 multiplies the density by 6.25;
S1 chose that approximate value using the observed density. At c=1 the action
is a benchmark, not a closure-derived or blind parameter choice.

Pure E^8/M^4 has fourth root E^2/M=.01051,.02485,.41068 meV for
E=160,246,1000 GeV respectively. The alpha^-2 enhancement above is essential;
choosing a TeV epoch alone neither fixes the SM transition nor its remnant.
Alpha is inflation-model calibrated. Indicative sensitivity is

\[
d\ln\rho_*=2d\ln|c|+2d\ln\lambda+8d\ln v-2d\ln\alpha-4d\ln M.
\]

No covariance or confidence interval for these theoretical assumptions is
supplied by solver refinement. The assumed absolute zero is not a quantum
solution of the cosmological-constant problem.

## Coupled homogeneous equations

For a spatially flat expanding FRW universe, with prime d/dN and N=ln a,

\[
D_t\dot q^A+3H\dot q^A+G^{AB}U_{,B}=0,\quad \rho_r'=-4\rho_r,
\]
\[
3M^2H^2=\tfrac12G_{AB}\dot q^A\dot q^B+U+\rho_r,
\quad\dot H=-\frac{G_{AB}\dot q^A\dot q^B+4\rho_r/3}{2M^2},
\]
\[
\dot\rho_s=-3H G_{AB}\dot q^A\dot q^B.
\]

The full solver evolves log H independently by Raychaudhuri and accumulates
W'=3H^2 G_AB q'^A q'^B. It checks Friedmann and rho_s+W=rho_s(initial).
Radiation is separately conserved. No potential is deleted by an imposed
switch or sent into an untracked receiver. Higgs potential energy is retained
as scalar motion and expansion work. The zero-temperature compressed runs
are not a physical thermal SM history; h is a homogeneous mean field.

For u=ln(phi/phi0),L=ln(phi0/M),B=bM/v the original metric has
G_uu=phi^2+B^2 h^2,G_uh=0,G_hh=1. The completed metric has

\[
G_{uu}=\phi^2+B^2h^2,\quad G_{uh}=B^2h(u+L),\quad
G_{hh}=1+B^2(u+L)^2,
\]
\[
\det G=\phi^2[1+B^2(u+L)^2]+B^2h^2>0.
\]

Connection covectors C_A=Gamma_ABC q'^B q'^C are

\[
C_u=\phi^2u'^2+2B^2hu'h',\quad
C_h=-B^2hu'^2\ (KP0),\quad C_h=2B^2(u+L)u'h'\ (KP1).
\]

Thus q''=-(3+H'/H)q'-G^-1(C+U_gradient/H^2). The metric finite-difference
control checks the connection separately. Dropping mixed terms in KP1 changes
the action; its fixed-h canonical trajectory cannot be used while h is moving.

## Exact KP1 preparation bound

Define z=bM(h/v)ln(phi/M). KP1 has ds^2=dphi^2+dh^2+dz^2, hence ell>=|Delta z|.
Starting at h=0,phi=phi0 and ending at h=v,phi=phi0 requires

\[
\ell/M\ge b|\ln(\phi_0/M)|.
\]

Positive-energy flat Einstein gravity implies

\[
d\ell/dN=\sqrt{G_{AB}q'^Aq'^B}\le\sqrt6M,\qquad
\Delta N\ge b|\ln(\phi_0/M)|/\sqrt6.
\]

For scalar energy fraction everywhere at most f, replace sqrt(6) by sqrt(6f).
For other endpoints use their actual Delta z. At c=1,b=3 this is 156.185
e-folds, or 1561.85 for f=.01. It is independent of the compressed simulation
and slow-roll approximation. It obstructs a rapid transition with these
endpoints, positive energy and flat GR. It is not a bound on KP0, all actions,
negative-energy histories, curvature, prolonged extra inflation, or a different
preparation. The phi>0 chart and endpoints are part of the domain.

For scale only: requiring initial m_chi<=H_ref, with m_chi^2=4rho*/(b^2 M^2),
gives b>=1.141 at c=1 (2.853 at c=2.5). A KP1 path within one e-fold demands
b<=.01921 (.01935). The first criterion is approximate freezing, not a general
acceleration theorem; the second follows exactly for its declared endpoints.

## KP0 Higgs invisible width

At h=v+eta,phi near phi0, chi=bM ln(phi/phi0) is the canonical light mode to
negligible O(phi0^2/M^2) corrections. KP0 contains

\[
\mathcal L\supset-\tfrac12(1+\eta/v)^2(\partial\chi)^2
\supset-\eta(\partial\chi)^2/v.
\]

For a nearly massless chi the matrix element is m_h^2/v; including identical
particles gives

\[
\Gamma(h\to\chi\chi)=m_h^3/(32\pi v^2)=.322584\ {\rm GeV}.
\]

The coupling does not decrease with b. A rounded reference SM visible width
.0041 GeV yields BR_inv=.987450. The named ATLAS 2023 observed 95% limit
.107 assumes SM Higgs production [S2]. With an escaping nearly massless mode
these assumptions reject KP0. We have not refitted the detector likelihood or
claimed the latest worldwide limit. KP1 was proposed to avoid this leading
portal; we do not certify its loop phenomenology, and its rapid-preparation
obstruction already matters here.

## Conditional late solution

Assuming h=v,phi=phi0 was prepared, both metrics along fixed h give

\[
\chi=bM\ln(\phi/\phi_0),\quad U\simeq\rho_*e^{2\chi/(bM)},\quad
m_\chi^2=4U/(b^2M^2).
\]

This is thawing quintessence, not permanent Lambda. With y=chi/M,u=y',
C=3M^2 H_ref^2,V=U/C,D=Omega_m e^-3N,R=Omega_r e^-4N,

\[
E^2=(D+R+V)/(1-u^2/6),\quad
\epsilon_H=-u^2/2-(3D+4R)/(2E^2),
\]
\[
y'=u,\quad u'=-(3+\epsilon_H)u-6V/(bE^2),\quad W'=E^2u^2.
\]

Test E^2 u^2/6+V+W=V(initial). Friedmann is algebraic in this solver and
is not counted as an independent check. Start at a=1e-8 with y=u=W=0.
Known H_ref=67.4 km/s/Mpc,Omega_m=.315,Omega_r=9.2e-5 fix the physical matter
densities; H(a=1) is an output, not reset to H_ref. The reference includes
existing cosmological inference assumptions and is not held-out data.

## Specific different action to test next

S3's simple assisted-relaxation example has a separate inflaton varphi and
remnant sigma, with

\[
G_{\varphi\varphi}=1,\quad G_{\sigma\sigma}=b^2M^2/(\varphi^2+\sigma^2),
\quad U=\tfrac12m_I^2\varphi^2+\tfrac12m_\sigma^2\sigma^2
+(1+c\sigma/M)V_H(h).
\]

S3 takes equal masses; its linear local potential needs a global completion,
for which a different action is proposed in its (19). The conditional minimum
is sigma0=-cV0/(M m_sigma^2), giving retained energy c^2 lambda^2 v^8/(32 M^2
m_sigma^2). When varphi>>|sigma| its canonical mass is approximately
m_sigma varphi/(bM); when varphi=0 the pole is b^2M^2/sigma^2. Thus our two
objections to Higgs-controlled realizations do not directly reject this action.

During inflation its local canonical coordinate is chi_inflation approximately
bM sigma/varphi, not the late logarithmic coordinate. Fluctuations must be
normalized in the appropriate stage and their distribution transported through
reheating; a tiny mean does not establish a tiny remnant population. No such
numerical, stochastic or global-viability calculation is claimed here. S3
still assumes zero absolute broken-phase vacuum energy, a selected mass scale,
sufficient inflation and Higgs stabilization. This is a bounded next question,
not a rescue already earned by the present runs.
