# Equations and limits

22 October 2026 (America/Chicago). Same author as the code.
Natural units, reduced Planck mass M=2.435e18 GeV.

## 1. Adopted model and state

This is a conventional asymmetric scalar-composite / Q-ball-like dark-matter
model with a Higgs portal, of the Friedberg--Lee--Sirlin type. Its potential is

\[
 U={\lambda\over4}(\phi^2-v^2)^2+g^2\phi^2|X|^2+\lambda_X|X|^4
 +{m_{D0}^2+h_D^2\phi^2\over2}D^2+{\lambda_D\over4}D^4,
\]

with the proposed separate coupling
\(\Delta U=(\kappa/2)(\phi^2-v^2)(H^\dagger H-v_H^2/2)\).
Here v=1 GeV, g=.1, lambda=.01, lambda_X=.001, h_D=.01,
lambda_D=.001, m_D0=.1 GeV; kappa=6.209329563372114e-10 is inherited
from an illustrative detector anchor, not a demonstrated allowed point.

The supplied single-sign coherent state has charge density
n0=1.8218345666575245e-13 GeV^3 at visible T=20 MeV. This is a stated
initial condition, not production from an empty sector. Its primordial charge
yield YQ=4.8294009076988175e-9 was originally calibrated to all dark matter
in the chosen Q=39944 profile. We keep YQ fixed; final abundance need not match.
We do not use the reference Q as a dynamical constraint.

The prior packet derived the rotating background and full linear dispersion.
We consume that result, not rerun the original soliton/profile solver.

## 2. Controlled dilute evolution

Minimizing phi adiabatically at small X gives
phi^2=v^2-2g^2|X|^2/lambda and

\[
 U_{\rm eff}=m_X^2|X|^2-a|X|^4,\quad
 a=g^4/\lambda-\lambda_X=.009.
\]

For X=e^{-im_X t}Psi/sqrt(2m_X), neglecting fast antiparticles and the
mediator's gradients/time derivatives gives

\[
 i\partial_t\Psi=-{\nabla^2\Psi\over2m_X}
       -{a\over2m_X^2}|\Psi|^2\Psi,\qquad n=|\Psi|^2.
\]

Put k0^2=a n0/m_X, Gamma0=a n0/(2m_X^2),
x_tilde=k0 x, tau=Gamma0 t, Psi=sqrt(n0) psi. Then

\[
 i\partial_\tau\psi=-\widetilde\nabla^2\psi-|\psi|^2\psi,
 \quad \psi_{\rm uniform}=e^{i\tau},\quad
 \Gamma_k/\Gamma_0=\widetilde k\sqrt{2-\widetilde k^2}.
\]

k0=1.28049e-7 GeV and Gamma0=8.19826e-14 GeV. One unit of tau is
8.02868e-12 seconds. Gamma0/H=4.5965e8 in the original 20 MeV bath;
expansion changes by less than 1.5e-8 over these runs. Omitting expansion
here does not license omitting it during later accretion or relaxation.

Conserved integrals of the reduced equation are
Q=(n0/k0^3) integral |psi|^2 d^3x_tilde and
H_tilde=integral (|grad psi|^2-|psi|^4/2) d^3x_tilde.
The homogeneous energy's dimensionless sign is negative; energy drift is
reported relative to its nonzero initial value.

The box side is 4pi, with total Qbox=1.72191361e11, not one object.
The fastest-wavelength cube contains about 2.15239e10 charge, a length-scale
diagnostic. This is not a final object charge. The imposed RMS 1e-3 seed
and its .5<=k_tilde<=1.5 support are numerical initial conditions, not a
calculated primordial perturbation spectrum. Two seeds are not an ensemble.

Norm-conserving Strang splitting alternates exact local nonlinear phases
and exact spectral free evolution. dt<=min(.02,.04/max|psi|^2); the control
has maximum .01. Spectral stopping uses radial |k|>(2/3)pi/dx, not the
diagonal corner Nyquist. The last safe sample has both high-shell norm
fraction<=1e-6 and gradient fraction<=1e-3. Once peak density exceeds 8,
the code inspects these every five steps, as well as at every .1 output.
The first exceeding sample is archived as a stop, excluded from resolved
physical claims. There is no artificial stabilizing higher-power term.

The quartic NR truncation cannot describe the eventual dense core or massive
phi emission. Its continued collapse would not establish a singularity in
the bounded full potential. Even the stopped 128^3 run has only
f^2/v^2~2.8e-10, well below the full-model nonlinear scale.

For an illustrative bulk object with Q equal to a fastest-wavelength cube,
omega_floor=(lambda lambda_X)^(1/4)v=.0562341 GeV,
n_bulk=omega_floor sqrt(lambda/lambda_X)v^2=.177828 GeV^3,
R_bulk~3069 GeV^-1. The phi wall's representative width 2/ms is 14.14
GeV^-1. A direct uniform box of length 9.8137e7 GeV^-1, four cells
across that wall, would require 2.78e7 cells per side, 2.14e22 cells total,
and at least 1.03e24 bytes for six real field/momentum arrays alone.
Resolving full carrier oscillations over the onset is also prohibitive.
This is a direct-grid estimate, not a lower bound on adaptive or matched methods.
The next calculation needs a conservation-preserving NR-to-full-field
handoff and core refinement, plus box/seed convergence.

## 3. Energy budget changes with the actual population

For a supplied charge yield YQ and relaxed final population with charge-weighted
mean energy per charge ebar, B_object=YQ ebar and
B_release=YQ(m_X-ebar), ignoring initial kinetic/interaction corrections of
order 1e-12 here. Charge lost to free carriers must also be counted.
Energy can remain in excited objects; B_release is not a derived instant
emission history.

The old reference gives B_release=4.29400908e-11 GeV per visible entropy.
Entropy per photon today is 7.03943, so this is 3.02274e-10 GeV per photon
in that convention. If all charge instead relaxes into asymptotically large
ground-state objects, fixed YQ gives
B_object=2.71577e-10 and B_release=2.11363e-10 GeV per entropy.
Objects then supply .617221 of the previous all-DM normalization. We do not
silently recalibrate YQ; this bulk limit is not an evolved final population.

For prompt electromagnetic heat at post-annihilation redshift z,

\[
 {\Delta\rho_\gamma\over\rho_\gamma}
 = {B_{\rm EM}\over(\rho_{\gamma0}/s_0)(1+z)},\quad
 {\rho_{\gamma0}\over s_0}={3T_0\over2g_{*s0}}.
\]

With T0=2.7255 K and g*s0=43/11, the reference budget at z=1e5 gives
.00476457, before a deposition fraction. This refines the relayed rough
.003 estimate. Comparing with the named historical 6e-5 heat screen needs
an EM fraction below .0126 for prompt heating there. This is not a BBN bound.

## 4. Delayed scalar response

The free spectrum has ms=.141421, mX=.100000 and mD=.100499 GeV.
There is no massless hidden free particle. The s->XX* and s->DD decays
are closed. Higgs mixing permits visible decays. We inherit the prior
approximate SM-response lifetimes, tau_s=1.23796498e11 seconds nominal
and 2.17961373e9 seconds in the enlarged photon-response envelope.
These are declared response choices, not a probability distribution for
the exact lifetime.

For instant emission at a_f with initial gamma_i and scalar energy B_s,
Ys=B_s/(ms gamma_i), u=p/ms=sqrt(gamma_i^2-1)a_f/a, gamma=sqrt(1+u^2).
The surviving number fraction is

\[
 S(a)=\exp[-\int_{a_f}^a{d\ln a'\over\tau_s\gamma(a')H(a')}].
\]

The lab decay probability carries 1/gamma, cancelling the gamma in
the emitted energy. Assuming all visible energy immediately becomes heat,

\[
 {dR_\gamma\over d\ln a}
 = {m_sY_s a\over\rho_{\gamma0}/s_0}\,{S(a)\over\tau_s H(a)}.
\]

We use the adopted flat background H0=67.4 km/s/Mpc, Omega_m=.315,
Omega_r from T0 and N_eff=3.046, Omega_L fixed by flatness. This is a
calibrated cosmological input, not a derived vacuum density. The late
radiation density is extrapolated at high z; a_f includes the entropy
ratio (g*s0/10.75)^(1/3) at 20 MeV. The omitted e+- transition in H is
irrelevant to the late scalar decays but precludes using this as a BBN solver.

The heat proxy integrates z>=1000 with
Jbb=exp[-(z/1.98e6)^(5/2)], compared with historical 6e-5. It assumes
unit prompt heat deposition from the visible daughters; no electromagnetic
cascade, photon injection spectrum, mu/y likelihood, recombination fit or
updated experimental limit is computed. At grossly failing amplitudes a
linear small-distortion kernel is diagnostic, not an accurate distortion
prediction. A passing proxy is not a viability certificate.

For a 20 MeV source, reference B_s, gamma_i=1, nominal response:
median decay z=16135, heat proxy=.0317689, source fraction below .00188864.
The enlarged response gives .00405335 and fraction .0148026. Thus the
conditional slow freely surviving scalar completion fails this heat screen
even for early formation. Subsequent reactions or actual spectra can change it.
For gamma_i=1000 the nominal heat is 3.17689e-5, because early expansion
redshifts most kinetic energy before decay. This is a response kernel, not
a mechanism producing that hard spectrum. The gamma_i=1e8 row is formal:
one reference object's 355 GeV binding energy can give a single scalar at
most gamma~2511, so this row requires a different source or much larger objects.

Energy accounting explicitly verifies the fractions decayed, redshifted
and remaining add to one. In pure radiation domination, cold early
emission gives the alternative closed form
R_gamma=B_s/(rho_gamma0/s0) sqrt(2 H0 sqrt(Omega_r) tau_s) Gamma(3/2).

## 5. Existing massive conversion channel

The vertices in the stated potential are sss=6lambda v,
sDD=2h_D^2v, ssDD=2h_D^2 (and g replaces h_D for X X*).
Up to a common irrelevant sign, at initial s particles at rest,

\[
 {\cal M}_{ss\to DD}=2h_D^2+{(6\lambda v)(2h_D^2v)\over3m_s^2}
 -{2(2h_D^2v)^2\over m_s^2}
 =4h_D^2-{4h_D^4\over\lambda}=.000396.
\]

The identical final-state factor yields
sigma v_rel=beta_D |M|^2/(64pi ms^2)=2.74368e-8 GeV^-2.
ss->XX* has distinct daughters and amplitude 4g^2-4g^4/lambda=0
at this parameter point ONLY at threshold. The general COM angular
integral confirms that p/ms=.1 gives sigma v(XX*)=4.30060e-8 GeV^-2,
already exceeding DD=2.69475e-8. No all-temperature cancellation is claimed.

A cold, homogeneous, irreversible, DD-only radiation-era diagnostic obeys
dY_s/dln a=-(s sigma v/H)Y_s^2. For constant g*=g*s,

\[
 {Y_{s,\infty}\over Y_{s,f}}=
 [1+Y_{s,f}(s_f\sigma v/H_f)]^{-1}.
\]

The reference 20 MeV source has initial rate/H=1.76199 and residual
scalar fraction .362058. Multiplying its collisionless heat by this
residual is about .01150, still above the screen. This composition is
only a diagnostic: it neglects inverse DD->ss, finite-momentum ss->XX*,
kinetic evolution, accretion into objects and correlated/coherent debris.
It is not a rigorous extremum or a complete reaction network. At 35 GeV
with an assumed constant g*=86.25 the residual would be .0001145;
that is a different history from the fixed 20 MeV maximum and is not a
validated solution. Stable D adds matter and changes abundance accounting.

## 6. Dark radiation qualification

If an actual massless decoupled receiver carried B_DR emitted at a_f,
Delta_Neff=B_DR a_f/[(rho_gamma0/s0) c_nu],
c_nu=(7/8)(4/11)^(4/3). Reference values are 1.7585e-8 for 20 MeV,
.0209794 at z=1e5 and .209775 at z=1e4. These are counterfactual
bookkeeping coefficients, not a branch present in the action. Gravitational
waves are in principle massless; no order-one graviton branching is derived.

For the actual stable D field, w=u^2/[3(1+u^2)]. Emission at 20 MeV
requires gamma_i~1.08e8 just to reach p/m=1 at z=1100; even gamma_i=1e8
gives w=.1533, not radiation's 1/3. Modest-boost hidden waves become matter.
The formation and reaction calculation must supply their spectrum and
abundance before invoking N_eff or a dark-radiation escape.
