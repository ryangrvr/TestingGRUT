# Point-specific rates, capture and the core prototype

These calculations use C1-HP-03/v1, the proposed scalar/Higgs-portal EFT
at the explicit parameters in PARAMETER_POINT.json. The original slow-s
point, its 62% bulk figure and its detector/formation numbers are not used
as predictions of this point. Natural units unless stated.

## 1. Mass condition and free-pair annihilation

The vacuum masses are

\[
m_X=g_Xv=.0924090\ {\rm GeV},\quad
m_s=\sqrt{2\lambda_\phi}v=.226355\ {\rm GeV},\quad
m_D=\sqrt{g_X^2+h_D^2}v.
\]

The existing interaction gives s -> XX* only if
lambda_phi>2g_X^2. Necessary homogeneous binding requires
lambda_phi lambda_X<g_X^4. Both hold here. This is a restriction on
parameters; it is not a property of the failed lambda_phi=.01 example.

Let W be the invariant mass of a free X/X* pair and a=2g_X^2v its
trilinear coupling to s. At leading small mixing,

\[
\theta={\kappa v v_H\over m_h^2-m_s^2},\qquad
\mathcal P(W)=\theta\left[
{1\over W^2-m_s^2+i m_s\Gamma_s}
-{1\over W^2-m_h^2+i m_h\Gamma_h}\right].
\]

Both poles are retained with their relative sign. With Gamma_h=.0041 GeV,
the parent's stated visible-width model and a fixed narrow light width,

\[
\boxed{\sigma v_{\rm M}(XX^*\to{\rm SM})
={2a^2\over W}|\mathcal P(W)|^2\Gamma_{\rm SM}^{unit}(W).}
\]

This follows from flux 1/(4E_1E_2)=1/W^2 and the decay identity
integral dPhi |M_decay|^2=2W Gamma_unit. A separate spin trace gives
sigma v(ee)=a^2 |P|^2 (m_e/v_H)^2 beta_e^3/(4pi), reproducing it.
The on-shell hidden width regularizes the light pole. This leading
Breit-Wigner screen is not a full finite-density/early reaction network.

If a total late matter fraction f_X is a symmetric free population,
n_X=n_X*=f_X rho_DM/(2m_X). Consequently

\[
\mathcal Q_{inj}=2m_X\langle\sigma v\rangle n_Xn_{X^*}
={f_X^2\langle\sigma v\rangle\over2m_X}\rho_{DM}^2.
\]

The factor 1/2 is important for a complex species when using a constraint
written in the usual self-conjugate p_ann convention. The nominal cold
source has f_X=.0635026062, not the full binding fraction .0766838646:
some released energy is daughter kinetic energy rather than late rest mass.

The named Planck 2018 TT,TE,EE+lowE+lensing+BAO screen is
p_ann<3.2e-28 cm^3 s^-1 GeV^-1. We count all injected SM energy in
the comparison, without claiming a new deposition/recombination fit.
The nonrelativistic and an explicit light-resonance rate envelope are
recorded separately. The latter is restricted to 2m_X<=W<=2 GeV;
it is not an all-energy bound or a bound on trapped-core emission.

## 2. The proposed cold recycling cycle is not an open free-s channel

Let E(Q) be the mass of a stable ground-state object and
mu_-=E(Q)-E(Q-1). Stability against one free positive carrier implies
mu_-<=m_X. For an incident antiparticle at rest,

\[
E_{available}\leq m_X+\mu_-\leq 2m_X<m_s.
\]

Thus object(Q)+X* -> object(Q-1)+ONE FREE s is kinematically closed
throughout this stable-bound-object/hidden-decay class. Neglecting recoil
only enlarges the available energy, so it cannot reverse the conclusion.
An antiparticle born from s at rest has energy m_s/2; the same argument
gives available energy m_s/2+mu_-<m_s, so it is also insufficient.

At the test object Q=200000, neighboring integer-charge solutions give
mu_-=.0795821601 GeV, compared with omega=.0795821496 GeV. The cold
energy budget is .171991124 GeV, leaving a .054363685 GeV deficit for
a free s. At least that much incident kinetic energy is necessary before
this particular channel can open, up to a tiny extra recoil cost.

This does NOT prove absence of capture. Bound excitations, changes in
occupation, multi-particle channels, off-shell emission and hot incident
particles need their own dynamics. A vacuum branching fraction alone is
not a lifetime bound for a trapped excitation repeatedly interacting
inside an object. Those rates are not supplied by the argument above.

## 3. Daughter momenta and collisionless travel

For an s at rest,

\[
E_*={m_s\over2},\quad p_*={1\over2}\sqrt{m_s^2-4m_X^2},
\quad\gamma_*={m_s\over2m_X}=\sqrt{1.5},\quad\beta_*={1\over\sqrt3}.
\]

Each X daughter receives 20.7684 MeV of kinetic energy. For a parent
boost gamma_s and isotropic c=cos(theta*) in its rest frame,

\[
E_X=\gamma_s(E_*+\beta_s p_*c),\qquad
p_X=\sqrt{E_X^2-m_X^2}.
\]

Gauss-Legendre integration checks normalization, mass shell and energy.
At fixed emitted energy, increasing gamma_s reduces the scalar/daughter
number by 1/gamma_s. The hard-source sensitivity table explicitly reduces
the free rest-mass contribution at fixed original charge; it does not
silently keep the old abundance or refit the charge yield.

For collisionless p/m=u_i a_i/a, the comoving travel distance is

\[
\lambda_{FS}(a_f)=c\int_{a_i}^{a_f}
{u_i a_i/a\over\sqrt{1+(u_i a_i/a)^2}}{da\over a^2H(a)}.
\]

The adopted flat background uses the parent's H0, Omega_m, radiation and
calibrated vacuum. Early temperature labels include visible entropy
conversion. The nominal integral extrapolates the late radiation density;
an instantaneous e+- entropy-change sensitivity is supplied separately.
It is not a BBN background solution. In pure radiation domination,

\[
\lambda_{FS}={c u_i a_i\over H_0\sqrt{\Omega_r}}
\left[\operatorname{asinh}{a_f\over u_i a_i}
-\operatorname{asinh}{1\over u_i}\right],
\]

providing an independent analytic quadrature control. A free-streaming
length is not a Lyman-alpha likelihood or a complete matter transfer
function. Nonthermal mixed populations must be tested with their actual
production history and spectrum, not a thermal-WDM mass relabeling.

## 4. Full-field radial patch and an explicit initial-data map

Use tau=g_X v t, rho=g_X v r, X=v psi/sqrt(2), phi=v S. At D=0 and
leading flat-space/tree order the exact bounded two-field equations are

\[
\ddot\psi=\nabla_\rho^2\psi-(S^2+.1|\psi|^2)\psi,
\qquad
\ddot S=\nabla_\rho^2S-[3(S^2-1)+|\psi|^2]S.
\]

\[
\bar Q=-\int\!\operatorname{Im}(\psi^*\dot\psi)d^3\rho,\qquad
Q=\bar Q/g_X^2,\qquad E=(v/g_X)\bar E,
\]
\[
\bar E=\int\left[{|\dot\psi|^2+\dot S^2+|\nabla\psi|^2+|\nabla S|^2\over2}
+{3\over4}(S^2-1)^2+{S^2|\psi|^2\over2}+{.1\over4}|\psi|^4\right]d^3\rho.
\]

For the supplied Gaussian u=A exp(-rho^2/(2R^2)), solve the neutral static
elliptic equation with S'(0)=0 and S(infinity)=1. Set psi=u,
dot(psi)=-i u, dot(S)=0. This preserves the supplied canonical phase,
gradients and momenta rather than setting a target soliton charge.
The exact full-state charge is the NR norm N=integral |u|^2. Define

\[
H_{ad}=\int\left[{|\nabla u|^2+|\nabla S|^2\over2}
+{3\over4}(S^2-1)^2+{(S^2-1)|u|^2\over2}+{.1|u|^4\over4}\right]d^3\rho.
\]

Then Ebar=N+H_ad at this initial-data map, including the neutral gradient
energy. The same identity is exact on the shared discrete grid. The
reported continuum Gaussian norm is a separate quadrature reference.
This construction is not yet a match to the earlier cosmological 3D NR
trajectory: the supplied patch's amplitude/density hierarchy is reduced,
and the omitted higher-order slow momentum correction is not a completed
relativistic NR-EFT matching calculation.

Case A uses A=.18, R=60; case B halves A and doubles R. Charges follow
from those initial data, rather than being set to Q=200000. They are
about 3.90 and 7.79 million. Both are far denser and smaller than the
actual primordial-density onset patch. The same scalar parameters are
retained; cosmic abundance is not inferred from the scaled prototypes.
Local expansion and gravity are negligible at their microscopic scales.
The portal's background correction is higher order here. Classical D=0
stays zero; quantum D production inside cores is not computed.

The finite-volume radial gradient uses face area/distance and shell
volumes, so its force is the exact derivative of the discrete energy.
Complex-field drift and real-potential kick separately preserve total
Noether charge. A Strang velocity absorber lies outside the probe region.
Its lost charge and kinetic energy are integrated into a numerical
receiver, with every reported conservation residual including it. This
is an open numerical boundary treatment, not a new physical dissipator.

## 5. Field flux versus a complete particle source

At a probe radius, energy and charge fluxes are

\[
\mathcal F_E=-4\pi\rho^2\,[\operatorname{Re}(\dot\psi^*\psi')+\dot S S'],
\qquad\mathcal F_Q=4\pi\rho^2\operatorname{Im}(\psi^*\psi').
\]

Their windowed Fourier sums obey Parseval, checked directly against the
time samples. For a linear outgoing mode of signed frequency omega,
k=sign(omega)sqrt(omega^2-m^2) and
u_out=(F+(F'+F/rho)/(i k))/2. Incoming modes are kept separately.
Complex positive/negative frequency branches encode opposite charges;
the real scalar has both Fourier signs but only one physical species.
Energy divided by |omega| gives the corresponding number in g_X^2 units.
Free plane-wave controls verify these normalization factors.

Near threshold the local mass perturbation can be comparable to k^2.
The declared restricted band therefore requires k^2>10 max|delta m^2|
and k^2>1e-4. This controls one local coefficient, not every approximation
or the total particle-source uncertainty. The remaining soft/coherent
flux is explicitly unresolved. Increasing the detector radius improves
coverage but does not complete it. Incoming, window and finite-time
effects are recorded rather than absorbed into an artificial normalization.

Only that restricted outgoing band is folded through isotropic vacuum
s -> XX*/DD decays. Stable particle abundance and momentum bins include
direct charged waves and both decay daughters with exact decay-energy
and net-charge accounting. The calculation assumes the emitted modes
linearize before decay and neglects subsequent scattering/recapture.
Classical s is propagated without decay to diagnose its source; that
counterfactual propagation is not counted as an additional surviving
cosmological species after the decay fold. A full quantum/kinetic outer
receiver and a complete asymptotic particle population remain missing.
