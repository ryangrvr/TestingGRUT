# Fixed-mu q-theory: one response in three observational channels

21 October 2026. Codex contribution, based on Claude checkpoint
`8ec108c9807892a1c6e7afb81ce2b8a65864e4f1`. See [REPORT.md](REPORT.md),
[INPUTS.json](INPUTS.json), [calculate.py](calculate.py), and [SOURCES.md](SOURCES.md).
This is an author derivation and implementation of standard scalar-tensor results.
No mathematical or physical novelty, state-selection rule, or empirical confirmation is claimed.

## 1. The action and the local scalar

Use signature (-,+,+,+), the curvature convention R=12H² in de Sitter, and
the same action as the pinned brief:

\[
S=\int d^4x\sqrt{-g}\,[K(q)R/2-\epsilon(q)]+S_m[g,\psi],\qquad
q=\widetilde F/\sqrt{-g},\quad F=dA.
\]

The three-form equation fixes a spacetime constant
\(\mu=\epsilon'-K'R/2\). A prime on K or epsilon here denotes d/dq.
For a fixed value of that constant, the local equations are equivalently
obtained from a scalar action with potential \(U=\epsilon-\mu q\).
This equivalence and the f(R) connection are explicitly present in KV
0806.2805, section II; they are not GRUT additions.

Where K>0 and K' is nonzero, choose a fixed reference mass M and define

\[
F_s=K/M^2,\qquad g^E_{\mu\nu}=F_sg_{\mu\nu},\qquad
\chi=\sqrt{3/2}\,M\ln F_s.
\]

Up to the usual boundary term the Einstein-frame action is

\[
S=\int\sqrt{-g_E}\,[M^2R_E/2-(\partial\chi)^2/2-U/F_s^2]
  +S_m[F_s^{-1}g_E,\psi].
\]

Thus the matter conformal factor is
\(A(\chi)=\exp[-\chi/(\sqrt6 M)]\). Its dimensionless coupling is
\(\alpha=d\ln A/d(\chi/M)=-1/\sqrt6\); the unscreened scalar exchange
has relative Yukawa strength \(2\alpha^2=1/3\). There are two tensor
polarizations and one scalar on this branch. K'=0 is a separate
nonpropagating limit, not a case to which the scalar formulas can be
extended by dividing by K'. Positivity of K supplies the local kinetic
sign; a positive scalar mass is also required for the stable linear branch.

The metric field equations obey the identity
\(\nabla^\mu E_{\mu\nu}=-q\partial_\nu\mu\). A proposed time-dependent
mu requires a receiver/source stress tensor: ordinary matter conservation
cannot simultaneously be kept as an independent equation with mu changing.
For homogeneous matter the identity gives
\(\dot\rho_m+3H(\rho_m+p_m)=q\dot\mu\).
All calculations here keep mu constant.

## 2. Curvature response, mass conventions, and approximations

Where q can be eliminated, set

\[
f(R)=K(q_R)R-2\epsilon(q_R)+2\mu q_R,\qquad
\epsilon'(q_R)-\mu=K'(q_R)R/2.
\]

This f has mass dimension four; the gravitational action is f/2. Then

\[
f_R=K,\qquad K_R=f_{RR}={K'^2\over2(\epsilon''-K''R/2)}.
\]

The trace perturbation about a constant-curvature background has mass

\[
m_{\rm tr}^2={K-RK_R\over3K_R}.
\]

At Minkowski equilibrium this is \(2K\epsilon''/(3K'^2)\), in agreement
with Claude's stage-1 result. Distinguish it from the leading
quasi-static coefficient \(\nu^2=K/(3K_R)=m_{\rm tr}^2+R/3\).
Confusing these two masses silently drops a curvature term. Both agree
to the stated order when that curvature is negligible compared with
\(p^2+m_{\rm tr}^2\), with physical wavenumber p=k/a.

Below assume linear perturbations, pressureless minimally coupled matter,
negligible matter anisotropic stress, wavelengths well inside the horizon,
and a matter-induced scalar mode dominating the free oscillatory mode.
These conditions are additional state/regime assumptions. The full
relativistic scalar perturbations are not integrated in this packet.

## 3. Matter and lensing

Write the Jordan metric as
\(ds^2=-(1+2\Psi)dt^2+a^2(1-2\Phi)d\mathbf x^2\), so GR has Phi=Psi.
Define \(G_{\rm ref}=1/(8\pi M^2)\), and define the response coefficients by

\[
-p^2\Psi=4\pi G_{\rm ref}\,g_{\rm dyn}\rho_d\delta_d,\qquad
-p^2(\Phi+\Psi)/2=4\pi G_{\rm ref}\,\Sigma\rho_d\delta_d,\qquad
\eta=\Phi/\Psi.
\]

G_ref is a reference normalization; identifying it with a Cavendish
measurement needs the local scalar response and environment.
The leading quasi-static constraint equations can be expressed as

\[
2Kp^2\Psi\simeq-\delta\rho_d-p^2\delta K,\qquad
K(\Phi-\Psi)\simeq\delta K,\qquad
(p^2+\nu^2)\delta K\simeq\delta\rho_d/3.
\]

Define \(s=p^2/(p^2+\nu^2)\). Their solution is

\[
g_{\rm dyn}={M^2\over K}(1+s/3),\qquad
\eta={1-s/3\over1+s/3},\qquad
\Sigma={M^2\over K}.
\]

In particular,

\[
\boxed{g_{\rm dyn}(1+\eta)/2=\Sigma},\qquad
{g_{\rm dyn}\over\Sigma}-1={s\over3}.
\]

The first equality is also inherent in the definitions of the two
potentials. The theory content is the scale-independent leading
\(\Sigma=M^2/K\), the fixed scalar strength, and the particular mass
dependence of s. A general model can satisfy the definitional equality
and still violate these additional restrictions. These expressions agree
with De Felice--Kobayashi--Tsujikawa 1108.4242, sections IV and V.A, after
translating their opposite sign for the spatial potential.

Dust growth consequently obeys

\[
\ddot D+2H\dot D-{\rho_d\over2K}(1+s/3)D\simeq0.
\]

It needs the background H, initial D and its derivative, and the scalar
mass/history. Neither a primordial amplitude nor the background cosmology
is selected by this equation. The numerical values below use matched
supplied initial growth data, not a CMB-normalized power spectrum.

## 4. Tensor propagation and a conditional cross-channel relation

Expansion of the same action in transverse traceless h_ij gives

\[
S_T={1\over8}\int dt\,d^3x\,a^3K
[\dot h_{ij}\dot h_{ij}-(\partial_kh_{ij})^2/a^2].
\]

An overall polarization convention does not affect the wave equation:

\[
\ddot h+(3H+\dot K/K)\dot h+(k^2/a^2)h=0,\qquad c_T=1.
\]

In conformal time, \(v=a\sqrt K\,h\) obeys
\(v''+[k^2-(a\sqrt K)''/(a\sqrt K)]v=0\).
In geometric optics its amplitude is proportional to \(1/(a\sqrt K)\).
After the emission waveform is independently accounted for, the
propagation distance ratio on a homogeneous branch is

\[
{d_L^{\rm GW}(z)\over d_L^{\rm EM}(z)}=\sqrt{K(0)/K(z)}.
\]

Combining the tensor result with the preceding quasi-static limit gives

\[
\boxed{\left({d_L^{\rm GW}(z)\over d_L^{\rm EM}(z)}\right)^2
\simeq{\Sigma(z)\over\Sigma(0)}}.
\]

This is a conditional consistency relation of the familiar conformally
coupled scalar sector. It is shared with the corresponding f(R)
comparator; satisfying it would not identify GRUT uniquely. Screening
can make source and detector K differ from the cosmological K used in
lensing. Source dynamics, scalar radiation, and propagation through such
environments require separate modeling before substituting observed
standard-siren distances. The unscreened homogeneous relation is not a
universal observable equality for every screened q-theory realization.

## 5. An exactly specified local branch for a worked calculation

Supply K=K_*+kappa(q-q_*), with K_*>0, and
\(U_\mu=c+b(q-q_*)^2/2\), b>0. Here epsilon'(q_*)=mu.
This is an exact quadratic/linear representative, or a local truncation
of generic functions if their higher derivatives are negligible; it is
not derived from closure. Eliminating q gives

\[
f(R)=K_*[R+R^2/(6m^2)-2\Lambda],\quad
m^2={2K_*b\over3\kappa^2},\quad\Lambda=c/K_*.
\]

For a fixed epsilon, q_* and c depend on the chosen mu. Taking c positive
by that choice calibrates a residual; the action has not predicted it.
For a q-independent addition C at fixed mu, c shifts by C. These
response formulas supply no protection against that shift.

Set K_*=H_ref=1 in the integration, F_s=1+x, and v=dot F_s. The exact
homogeneous equations, with no imposed GR background, are

\[
\dot x=v,\quad
\dot v=-3Hv-m^2x+(\rho_d+4\Lambda)/3,\quad
\dot H=m^2x/2-2H^2,\quad
\rho_d=3\Omega_{m,\rm ref}a^{-3},\quad
\rho_r=3\Omega_{r,\rm ref}a^{-4}.
\]

Their Friedmann constraint is

\[
C=3(1+x)H^2+3Hv-\rho_d-\rho_r-\Lambda-3m^2x^2/4=0,
\qquad\dot C=-4HC.
\]

We solve it for the initial H and evolve these equations directly,
alongside the approximate dust growth equation. The numerical label
Omega_Lambda,ref selects Lambda in units; the finite-mass endpoint H
is an output and need not equal H_ref. All inputs and initial data are
in INPUTS.json. The positive residual is supplied in those runs.
An additional run sets Lambda=0 with the same matter/radiation reference
densities and analogous initial state. It tests this late-time branch,
not every possible high-curvature or nonadiabatic q-theory state.

## 6. Limits of the calculation

Exact background integration is not exact relativistic perturbation
evolution. The recorded H²/(p²+nu²) values are dimensional indicators of
one quasi-static correction, not certified total error bars. Also record
p/H; discarded horizon terms can scale as H²/p². Rapid free scalar
oscillations, initial isocurvature, radiation perturbations and screening
can require going outside this approximation. Numerical refinement only
controls integration error within the specified equations.

The wave control uses a prescribed smooth scale factor and K to test
the propagation equation independently of the cosmological solver. It
does not supply another solution of the full background theory.
No phase-transition mechanism, microscopic dissipation, loop stability,
hidden-sector abundance or consciousness dynamics is inferred here.
