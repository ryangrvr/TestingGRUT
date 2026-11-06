# X7 · Spatial source and coherent transport

The X6 Galerkin reduction resolves one periodic spatial coordinate
of length $L$. The transverse momentum has two components and is
integrated radially. A homogeneous metric responds to the spatial
average through one scale factor. Directional fermion stress is
available as a diagnostic.

## Signed vacuum response

For each transverse node, the one-particle Hamiltonian is the finite
projection of

$$
 h_j=a^{-1}\boldsymbol\alpha\cdot\mathbf q+
       \beta M_j e^{i\gamma^5\phi/f}.
$$

Each occupied projector follows its own Hermitian Hamiltonian. The
Majorana energy carries the factor $1/2$ in the trace. For physical
mass $M=M_0$ and auxiliary masses $M_1<M_2$, define $w_j=c_j M_j^2$
with $c_0=1$. The other weights are

$$
 w_1=-M^2\frac{M_2^2-M^2}{M_2^2-M_1^2},\qquad
 w_2=M^2\frac{M_1^2-M^2}{M_2^2-M_1^2}.
$$

They obey $\sum_j w_j=\sum_j w_j M_j^2=0$. Energy, pressure, force
and momentum use the same signed weights. The physical projector
alone defines physical particle occupations. Each sector subtracts
its zero-phase, fixed-$a$ vacuum in stress.

Write $E_n=\sqrt{M_j^2+(q_n^2+q_\perp^2)/a^2}$ and
$k=2\pi r/L$ for external harmonic $r$. At fixed $a$ the quadratic
response on a sharp cylinder is

$$
 \Pi_Q(\omega,k)=\frac{1}{f^2a^3}\sum_{j,\ell} w_j W_\ell
 \left[
   \frac12\sum_{b=0}^{r-1}(E_{-N+b}^{-1}+E_{N-b}^{-1})
   +s\sum_{n=-N}^{N-r}
       \frac{E_n+E_{n+r}}{2E_n E_{n+r}[(E_n+E_{n+r})^2-(\omega+i0)^2]}
 \right],
 \qquad s=(k/a)^2-(\omega+i0)^2 .
$$

$W_\ell$ is the comoving radial integration weight including $1/L$.
The first term is the unpaired contribution at the longitudinal
boundary. It is independent of frequency. The particle-hole numerator
$E_n E_{n+r}+M_j^2+(q_nq_{n+r}+q_\perp^2)/a^2$ equals
$[(E_n+E_{n+r})^2-(k/a)^2]/2$. This identity leaves the term
proportional to $s$ after the virtual transition and the quadratic
contact vertex are combined on paired modes. The spectral check
computes the same kernel from the contact vertex and the finite
Hamiltonian's occupied-to-empty transitions.

At large physical loop momentum $p$, the unpaired summand has
expansion $1/E=p^{-1}-M_j^2/(2p^3)+3M_j^4/(8p^5)+O(p^{-7})$.
The signed constraints cancel the first two terms at every fixed
harmonic. At fixed cylinder aspect ratio and auxiliary masses the
remaining edge scales as $Q^{-3}$ because the transverse measure
supplies two powers of $Q$. The physical-only edge grows as $Q$.

An independent reference removes the cylinder while retaining the
periodic direction. Combining the propagators with one Feynman
parameter and integrating frequency and transverse momenta gives

$$
 \Pi_L(\omega,k)-s Z_L(0,0)
 =\frac{s}{8\pi a L f^2}\sum_{j,n}w_j\int_0^1\!dx\,
  \left[
    \frac{1}{\sqrt{M_j^2+(1-x)(q_n/a)^2+x((q_n+k)/a)^2-x(1-x)(\omega+i0)^2}}
    -\frac{1}{\sqrt{M_j^2+(q_n/a)^2}}
  \right].
$$

The real-frequency auxiliary checks stay below the pair threshold.
This unbounded periodic reference retains finite dependence on the
auxiliary masses after local kinetic matching. In the infinite-volume
limit the matched kernel is

$$
 -\frac{s}{16\pi^2f^2}\sum_j w_j\int_0^1\!dx\,
       \log\left(1+\frac{x(1-x)s}{M_j^2}\right).
$$

For an auxiliary mass large compared with the external frequency and
wavenumber, its leading residual is $-c_j s^2/(96\pi^2 f^2)$.
This describes reversible response below threshold. The measured
time-domain sensitivity can also include regulator-dependent
preparation and finite-cutoff stress. The fixed-metric quadratic
calculation isolates one contribution to the expanding result.

## Force and stress ledger

At finite cutoff let

$$
 Z=\frac{1}{4 f^2 a^3}\sum_{j,\ell,n} w_j W_\ell E_n^{-3},
 \qquad A=1-Z,\qquad
 Z_g=\frac{1}{4 a^3}\sum_{j,\ell,n} w_j W_\ell p^2 E_n^{-5}.
$$

The scalar kinetic and gradient terms use coefficient $A$.
The canonical momentum is $a^3A\dot\phi$.
Fermion force density $J$ is the scalar derivative of the signed
mass interaction energy. Its reciprocal work is
$\langle J\dot\phi\rangle$. The split propagator integrates
the local mass rotation and its opposite scalar impulse together.

For $A'=dA/d\log a$ and $g=a^{-1}\partial_x\phi$, the scalar
energy and pressure are

$$
 \rho_\phi=V+\tfrac12 A\langle\dot\phi^2+g^2\rangle,\qquad
 P_\phi=\tfrac12 A\langle\dot\phi^2\rangle-V
       +\tfrac16 A'\langle\dot\phi^2\rangle
       -\tfrac16(A+A')\langle g^2\rangle.
$$

$V$ denotes the spatially averaged potential here. Expansion uses
$H^2=\rho/(3M_P^2+Z_g/2)$ and $P=P_\phi+(P_{F,x}+2P_{F,\perp})/3$.
The volume variation at fixed scalar canonical momenta and fermion
state checks the pressure. The total longitudinal momentum checks
the translational Noether quantity $a^4 P_x^{\rm mom}$.

For matched transfer, move the known scalar kinetic and gradient
subtractions into the fermion ledger. The changes give

$$
 \rho_{F,R}=\rho_F-\tfrac12 Z\langle\dot\phi^2+g^2\rangle,\qquad
 P_{F,R}=P_F-(Z/2-A'/6)\langle\dot\phi^2\rangle
                    +(Z-A')\langle g^2\rangle/6 ,
$$

$$
 J_R=J-Z\ddot\phi-H(3Z-A')\dot\phi+Za^{-2}\partial_x^2\phi .
$$

The code compares the integral of $a^3\langle J_R\dot\phi\rangle$
with the change in $a^3\rho_{F,R}$ plus $\int P_{F,R}\,d(a^3)$.
The metric matching coefficient stays in the Friedmann equation.
This is a finite-slab prescription. General curved spacetime also
requires matching gravitational and phase-derivative operators,
including anomaly terms where applicable.
[Axion anomalies and derivative-coupling regularization](https://arxiv.org/abs/2112.07645)

## Exact homogeneous transport variables

In the X5 chiral frame let $k=hq/a$ for helicity $h$,
$E=\sqrt{M^2+k^2}$, $\mathbf n=(M,0,k)/E$, and
$\mathbf e=(k,0,-M)/E$. For the Bloch vector $\mathbf r$,
the free occupation is $f_0=(1+b)/2$ with
$b=\mathbf n\cdot\mathbf r$. The two coherence components
are $u=\mathbf e\cdot\mathbf r$ and $y=r_y$. Positivity
gives $u^2+y^2\le 4 f_0(1-f_0)$.

Set $\delta=\dot\phi/(2f)$, $h_\parallel=E+k\delta/E$,
$h_\perp=-M\delta/E$, and $s_b=HMk/E^2$. The coherent
equations expressed in these variables are

$$
 \dot b=2h_\perp y+s_bu,\qquad
 \dot u=-2h_\parallel y-s_b b,\qquad
 \dot y=2h_\parallel u-2h_\perp b .
$$

The vacuum decay sensitivity adds $-\gamma(b+1)$,
$-\gamma u/2$, and $-\gamma y/2$ to these equations with
$\gamma=\Gamma_{\rm rest}M/E$. The lost occupation energy
enters the visible radiation ledger. Scalar and metric driving
remain live.

The mode contribution to $P_F-P_{\rm gas}$ is
$-W kM u/(6a^3E)$, and its contribution to $J\dot\phi$
is $-WM\delta y/a^3$. For collisionless evolution,
$d(a^3\rho_F)/dt=a^3[J\dot\phi-3HP_F]$ follows from
$\dot E=-H k^2/E$ and $\dot f_0=\dot b/2$.
Reflecting $y$ gives a valid instantaneous covariance with
the same free occupations, energy and pressure and opposite work.
An instantaneous closure supplied only occupations and stress
cannot predict the work for both states. The replay evolves the
archived state with both coherences present.

The frozen-tail calculation fixes $a$, removes the scalar drive
and retains the vacuum decay law. Its exact worst-phase coherent
pressure envelope is

$$
 \sum_{j,h,q}\frac{W |k|M}{6a^3E}
       \sqrt{u^2+y^2}\,e^{-\gamma t/2}.
$$

This supplies a sufficient pressure criterion for arbitrary phases
in the frozen control. The live H55 run measures actual pressure
and work in the evolving background.
[Coherent quantum transport](https://arxiv.org/abs/0811.0936)
