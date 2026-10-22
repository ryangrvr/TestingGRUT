# One action: mass, organization and readout

All new results here concern the frozen Candidate 1 classical flat-space
saddle and linear elastic D scattering. Natural units; v below means the
**new scalar scale v_phi**, not the electroweak Higgs value v_H. The weak
Einstein response is calculated to leading order in compactness. No quantum
vacuum level, initial cosmic state or measurement postulate is derived.

## 1. Stationary objects at fixed charge

Use X=exp(-i omega t)f(r)/sqrt(2), phi=phi(r), D=0 and

\[
U_0=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+\frac{g_X^2}{2}\phi^2f^2+\frac{\lambda_X}{4}f^4.
\]

The charge and vacuum-relative energy are

\[
Q=\omega I_f,\quad I_f=\int f^2d^3x,\qquad
E_Q=\frac{Q^2}{2I_f}+\int\left[
\frac{(\nabla f)^2+(\nabla\phi)^2}{2}+U_0\right]d^3x.
\]

Stationarity at fixed Q gives

\[
f''+\frac2r f'=(g_X^2\phi^2+\lambda_Xf^2-\omega^2)f,
\quad
\phi''+\frac2r\phi'=[\lambda_\phi(\phi^2-v^2)+g_X^2f^2]\phi.
\]

Set rho=g_X*v*r, F=f/v, S=phi/v, Omega=omega/(g_X*v),
a=lambda_phi/g_X^2 and b=lambda_X/g_X^2. Their equations are

\[
F''+2F'/\rho=(S^2+bF^2-\Omega^2)F,
\quad S''+2S'/\rho=[a(S^2-1)+F^2]S.
\]

An additional accumulator J'=rho^2 F^2 and the condition
4*pi*Omega*J(L)=g_X^2 Q fix the integer charge. Regularity and vacuum
conditions complete the boundary problem. The first solve uses the pinned
archived profile only as an initial guess, not as a rescaled answer.
Continuation reaches every prescribed charge; larger-domain/tighter solves
check the three main backgrounds. Dimensionless integrals obey

\[
\bar E=4\pi\int\rho^2\left[
\tfrac12(F'^2+S'^2)+\tfrac12\Omega^2F^2
+\tfrac a4(S^2-1)^2+\tfrac12S^2F^2+\tfrac b4F^4\right]d\rho,
\quad E=\frac v{g_X}\bar E.
\]

The envelope identity dE/dQ=omega becomes
d[E/(g_X v)]/dQ=Omega. The central three-charge difference agrees with
Omega to 1.12e-10 in this run. Positive finite-grid fixed-Q radial amplitude
Hessians support the local stationary branch; they do not prove all-channel
quantum stability, fragmentation resistance or a lifetime.

## 2. An exact classical mass–depletion sum rule

There is only one dimensional scale v in this isolated classical scalar
energy. At fixed dimensionless couplings and fixed Q, scaling the solution
shows E(v,Q)=v*epsilon(Q). The fixed-Q envelope theorem gives

\[
\frac Ev=\left.\frac{\partial E_Q}{\partial v}\right|_{f,\phi}
=\lambda_\phi v\int(v^2-\phi^2)d^3x.
\]

Therefore every localized stationary solution in this class satisfies

\[
\boxed{E=\lambda_\phi v^2\int(v^2-\phi^2)d^3x.}
\]

This is not limited to a thin-wall or large-Q approximation. The numerical
pattern motivated this strengthening after the initial run; see
[AMENDMENT_01.md](AMENDMENT_01.md). Standard scaling identities are the
method, and literature novelty is not asserted.

An independent derivation uses the integrated equations and virial theorem.
Define G_f=integral (grad f)^2, G_phi=integral (grad phi)^2,
C=g_X^2 integral phi^2 f^2, L=lambda_X integral f^4,
B=lambda_phi/4 integral (phi^2-v^2)^2 and T=omega^2 I_f/2.
Then

\[
G_f+C+L=2T,\qquad
G_\phi+C+\lambda_\phi\int\phi^2(\phi^2-v^2)=0,
\]
\[
\tfrac12(G_f+G_\phi)+3(B+C/2+L/4-T)=0.
\]

Eliminating T and G_f gives E=G_phi+C+4B, and the phi equation gives
the boxed relation. The corresponding off-shell integral residual is the
sum of the two integrated equation residuals minus the virial residual;
all three are separately evaluated. In the refined backgrounds the mass
identity holds numerically within 2.2e-12 relative.

**Scope matters:** an independent bare X mass, another dimensionful coupling,
extra fields with nonzero backgrounds, higher operators, substantial gravity
or quantum/running corrections can change the identity. All symmetry-allowed
counterterms are still required by the underlying EFT. The calculation does
not promote a tree-level scaling identity to an exact quantum Ward identity
or protection of vacuum energy.

## 3. The same depletion controls the probe

The linearized D equation is

\[
[\Box-m_{D0}^2-h_D^2\phi(r)^2]D=0.
\]

After subtracting its asymptotic mass, the radial scattering problem is

\[
u_l''+\left[k^2-\frac{l(l+1)}{r^2}-V_D(r)\right]u_l=0,
\quad V_D=h_D^2(\phi^2-v^2).
\]

Here k is physical momentum; the code uses k/(g_X v), and dimensionless
V=(h_D/g_X)^2(S^2-1). The zero-energy scattering length convention is
u_0 proportional to r-a_D at large r, so delta_0(k)=-k*a_D+O(k^3).
Its first-Born value is

\[
a_D^{\rm B}=\frac1{4\pi}\int V_Dd^3x
=-\frac{h_D^2E}{4\pi\lambda_\phi v^2},\qquad
\sigma_D^{\rm B}(k\to0)=\frac{h_D^4E^2}{4\pi\lambda_\phi^2v^4}.
\]

Thus the leading low-energy readout amplitude is fixed by the object's mass
and the common couplings, with no independently adjustable shape parameter
at this order. This is a conditional relation for the proposed sector; no
such particles or experiment have been observed here.

The finite-coupling correction is controlled when the potential is attractive
and I=integral_0^infinity r|V_D|dr<1. The positive solution S cannot exceed
one by the scalar maximum principle, hence this branch is attractive.
For the zero-energy normalized wave psi, the integral equation is

\[
\psi(\mathbf x)=1+\int
\frac{|V_D(\mathbf y)|\psi(\mathbf y)}{4\pi|\mathbf x-\mathbf y|}\,d^3y.
\]

Its positive kernel has sup norm I for a radial potential: the response to
one is maximal at the origin. The Neumann series gives
1<=psi<=1/(1-I). Since a_D=-(1/4pi) integral |V_D|psi,

\[
\boxed{a_D^{\rm B}/(1-I)\le a_D\le a_D^{\rm B}<0.}
\]

All three exact zero-energy lengths lie inside this interval. I<1 also
suffices to exclude a negative eigenvalue of the radial operator in every
angular sector, using |u(r)|^2<=r integral |u'|^2. This does not exclude
gravitationally bound states, nonlinear D structures or inelastic processes.
The bound applies to the specified continuum potential; the numerical
background is convergence-tested rather than enclosed with interval arithmetic.

The scattering length here belongs to the short-range scalar potential in
the flat-space matching problem. In the full metric, gravity has a long-range
Coulomb tail and a plain short-range scattering length for the total force
is not defined. The G→0 limit precedes k→0 in the displayed sum-rule matching.
At the *finite* momenta used for readout, the exterior monopole gives the
small expansion parameter |eta_G|=G*E*(2omega_D^2-m_D^2)/k, at most 1.35e-7
in these backgrounds. This estimates a hierarchy, not a computed gravitational
phase correction or an apparatus-distance-independent total wavepacket phase.
An actual full-metric measurement must model or subtract that contribution.

The phase equation is independently organized from the regular radial equation:

\[
\delta_l'=-\frac V k
[\cos\delta_l\,\hat j_l(kr)-\sin\delta_l\,\hat y_l(kr)]^2,
\quad \delta_l(0)=0.
\]

Both equations are integrated as different variables in one adaptive solver;
they share its time steps and are an internal formulation check. The
radial endpoint is matched to Riccati–Bessel functions. A separate
zero-energy integration checks delta_0/k, and a zero-potential control
checks the phase convention. The source profile and scattering tolerances
are refined together on a larger domain.

## 4. Gravitational mass, cosmic matter and vacuum

At leading weak Einstein order a stationary localized conserved stress obeys
the Laue identity integral T_ij=0. Its active monopole mass is therefore
integral (T00+sum_i Tii)=E. This is checked against the same profiles;
the central clock shift is a separate profile-dependent output, not fixed
by E alone. Einstein gravity is adopted, not derived from this sum rule.

For a dilute cold population of well-separated objects with no charge exchange,
the long-distance matter action begins as -sum_A E(Q_A) integral d tau_A.
Number conservation then gives an approximately dust-like density scaling
as a^-3, with velocities, interactions and radiation requiring corrections.
The objects' interior potential energy contributes to their mass; it does
not create an extra volume-filling cosmological constant. Nor does clumping
remove the classical charge-energy bound rho_exc>=m_star|n_Q|.

These statements are vacuum-relative. The cosmological coefficient is the
total renormalized vacuum density, including matter contributions and its
counterterm. The frozen local benchmark matching does not derive it. No
new cosmological matching or full population history is adopted in this packet.

## 5. A prepared charge label is readable in the elastic approximation

Two preparations Q0=19972 and Q1=39944 encode a classical binary label.
Stationarity and spatial isolation are supplied; global charge conservation
does not prove protection against charge escaping a finite region. The
travelling D field is a readout carrier, rather than the storage mode.

For a monochromatic input (|l=0>+|l=1>)/sqrt(2), the outgoing state is

\[
|\psi_Q\rangle=\frac{e^{2i\delta_0(Q)}|0\rangle
+e^{2i\delta_1(Q)}|1\rangle}{\sqrt2}.
\]

One partial wave alone has only an overall phase and carries no such
monochromatic label information. With the declared finite normalized
momentum density w(k), the overlap is

\[
C=\frac12\int w(k)\sum_{l=0,1}
e^{2i[\delta_l(Q_1,k)-\delta_l(Q_0,k)]}\,dk.
\]

For equal priors and N ideal independent probes, the standard Helstrom result is

\[
P_{\rm opt}(N)=\frac{1+\sqrt{1-|C|^{2N}}}{2}.
\]

No apparatus implementing this optimum is claimed. An explicit fixed
two-mode projection has
p_plus(Q)=integral w(k)[1+cos(2theta_Q(k)-beta)]/2,
theta_Q=delta_0-delta_1. Beta is fixed at the two central phases' sum plus
pi/2. The exact binomial likelihood decision then gives its actual
N-probe error. This supplies a measurement prescription, not a constructed
autonomous physical detector.

The run obtains |C|^2=.999790725574, P_opt(1)=.50723316, N_opt(.95)=7935,
and N_fixed(.95)=13270. The repeated-probe counts presume the background
stays unchanged. The total incident kinetic budgets are .246%/.411% of
the smaller object's mass; these are resource diagnostics, not a computed
backreaction or lifetime certificate. Rest energy of transmitted probes
also requires a source. No dynamical writing or cumulative disturbance was
simulated, and the original SM portals do not supply a practical D detector.

For a coherent joint register/reservoir state with fixed *total* charge,
the conditional outgoing probes multiply its joint off-diagonal coherence
by C^N. A compensating reservoir is essential if its two object charges
differ; this is not a locally observable superposition of different total
charge sectors. Defining visibility V_N=|C|^N gives

\[
(2P_{\rm opt}-1)^2+V_N^2=1.
\]

This is standard unitary information transfer with adopted Born probabilities,
not a new collapse law. It links the readable information and dephasing to
the same computed scattering matrix. It neither constructs a conscious
observer nor explains subjective experience.
