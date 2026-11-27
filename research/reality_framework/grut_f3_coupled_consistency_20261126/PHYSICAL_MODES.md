# Physical modes and the controlled matter limit

## Linear dynamics on the common FLRW solution

On the rolling branch, choose comoving gauge for the active
scalar. The physical curvature perturbation is $\zeta$,
the spectator perturbation is $\delta\sigma$, and
$\gamma_{ij}$ is transverse and traceless.
Foliation time is proper time on the homogeneous normal
congruence. Put $\epsilon_H=-\dot H/H^2=3(1+2x)/(1+3x)>0$.
The constant spectator contributes stress first at
second order and decouples in the linear mode equations.
Variation of the scalar lapse and shift constraints at
first order gives $\delta N=\dot\zeta/H$ and, with
$N_i=\partial_i\psi$ in the scalar sector,
$$
 \psi=-\zeta/H+\psi_1,\qquad
 \partial^2\psi_1=a^2\epsilon_H\,\dot\zeta/c_{\rm a}^2 .
$$
The inverse spatial Laplacian acts on nonzero periodic
modes. The homogeneous part belongs to the background data.

Substitution into the expanded action, integration by parts
on the periodic slices, and the background equations give
the curvature action
$$
 S^{(2)}_\zeta=M^2\int a^3 d^3x\,dt\,\frac{\epsilon_H}{c_{\rm a}^2}
    \left[\dot\zeta^2-c_{\rm a}^2 a^{-2}(\partial_i\zeta)^2\right].
$$
In particular the kinetic factor is $X(P_{1,X}+2XP_{1,XX})/H^2$
and the gradient factor is $XP_{1,X}/H^2$.
The canonical control replaces $c_{\rm a}$ by 1 and has the
same $\epsilon_H$. The spectator action has kinetic coefficient
$1/2$ and metric speed. The tensor action has prefactor
$M^2/8$ multiplying
$a^3[\dot\gamma_{ij}^2-a^{-2}(\partial_l\gamma_{ij})^2]$.
These reduced actions display all four propagating degrees.
The FLRW lapse, shift and their scalar constraints add no
kinetic modes. Positive coefficients and the nested
cones give the perturbative hyperbolicity needed in
the experiment. Global nonlinear stability and long-mode
cosmological growth are separate questions.
[Garriga and Mukhanov](https://arxiv.org/abs/hep-th/9904176)

Use the future-directed fluid normal
$u^\mu=-\nabla^\mu\phi/\sqrt{2X}$ on the active solution.
Local frequencies use $-u^\mu k_\mu$ and the physical
wavenumber uses the projection orthogonal to $u^\mu$.
The tensor and spectator speeds in this frame both equal
one. The same normal results using the monotone matched
canonical scalar. Changing from proper time to either
field clock on the common homogeneous branch rescales
all rates in the same way and keeps their speed ratios.
The acoustic effect is a response of the matter state
and kinetic law in its rest frame.

At finite perturbative spatial cutoff, each independent
real canonical mode has a nondegenerate positive kinetic
term and admits canonical quantization on a positive
oscillator Hilbert space. For example the curvature
coordinate is $v=z_\zeta\zeta$ with
$z_\zeta^2=2a^2 M^2\epsilon_H/c_{\rm a}^2$.
On the registered high-wavenumber band its quadratic
frequency squared is positive. The real time-dependent
quadratic Hamiltonians generate unitary linear evolution,
and their retarded commutators solve the same linear
mode equations. The continuum linear equations have the
characteristics of the reduced classical symbols.
Choice of a Gaussian state fixes covariances without
moving these characteristics. The untruncated hyperbolic
equations give finite propagation for this classical
linear EFT. A finite Fourier cutoff is used for the
oscillator quantization and smears local observables.
These statements concern reduced perturbations. They do
not establish the full quantum Hamiltonian constraint
algebra, a unique vacuum, or nonperturbative quantum gravity.

In the vanishing active-gradient limit the two kinetic
coefficients approach one. The direct scalar fluctuations
then obey the usual metric wave law. The $\zeta$ clock
description at exactly zero gradient is singular because
the clock has stopped. The fixture at $x=0$ checks the
kinetic limit only.

## Canonical finite-density matter comparison

A standard local relativistic matter model has a complex
scalar with radial coordinate $r$, angle $\theta$, and
potential $\lambda(r^2-v_0^2)^2/4$, with $\lambda,v_0>0$.
Its kinetic terms are
$-(\nabla r)^2/2-r^2(\nabla\theta)^2/2$. It too can be
coupled to the same Einstein–Hilbert action and spectator.
This supplies a matter completion check below the gravity
cutoff. In the slow-radial limit put $\phi=v_0\theta$ and
$\Lambda^4=\lambda v_0^4$. At fixed $X$ the stable radial
stationary value is
$$
 r_b^2=v_0^2+\frac{2X}{\lambda v_0^2}.
$$
Substitution gives $P_1=X+X^2/\Lambda^4$ at tree level
when radial derivatives are dropped. On the local flat
rotating background $\dot\theta=\mu$, the relation is
$r_b^2=v_0^2+\mu^2/\lambda$. This is a solution of the
constant-density flat matter equations.
[Relativistic scalar effective theory](https://doi.org/10.1007/JHEP09(2022)066)

For local radial perturbation $h=\delta r$ and
$s=r_b\delta\theta$, the quadratic Lagrangian density in
flat space is
$$
 \tfrac12(\dot s^2-|\nabla s|^2)+
 \tfrac12(\dot h^2-|\nabla h|^2-m_r^2 h^2)+2\mu h\dot s ,
 \qquad m_r^2=2\lambda r_b^2 .
$$
It has positive Hamiltonian from the squares
$(\pi_s-2\mu h)^2/2,\pi_h^2/2$ and the positive spatial
gradient and radial mass terms. The Fourier characteristic
polynomial for its frequencies is
$$
(k^2-\omega^2)(k^2+m_r^2-\omega^2)-4\mu^2\omega^2=0.
$$
Set $m_{\rm gap}^2=m_r^2+4\mu^2$, $\beta=4\mu^2/m_{\rm gap}^2$
and $q=k^2/m_{\rm gap}^2$. Both branches have nonnegative
squared frequency and the exact light branch is
$$
\omega_-^2 = k^2+\frac{m_{\rm gap}^2}{2}
 -\frac{m_{\rm gap}^2}{2}\sqrt{1+4\beta q}.
$$
Its leading speed squared is $1-\beta=m_r^2/m_{\rm gap}^2$,
the same as $(1+2x)/(1+6x)$ using
$X=v_0^2\mu^2/2$. Rationalizing the square root yields
the exact bound
$$
0\le\frac{\omega_-^2}{k^2}-(1-\beta)
=\frac{4\beta^2 q}{(1+\sqrt{1+4\beta q})^2}
\le\beta^2 q .
$$
This controls the local flat tree-level dispersive
approximation at nonzero $k/m_{\rm gap}$. Both microscopic
matter fields have the metric principal symbol because
the mixing has one time derivative. Their high-frequency
front is metric-null. Coupling them to the metric leaves
that principal statement unchanged. The observed
low-energy phonon cone and the microscopic front are
different levels of this effective description.

The gravity calculation uses the $P(X)$ truncation on an
exact FLRW solution. The separate rotating-state result
controls flat-space matter matching at finite momentum.
Slow evolution relative to the gap is needed to combine
these descriptions. The packet has no uniform error
certificate for the full time-dependent coupled radial
model or its quantum-loop corrections. It makes no
claims about physics at the Planck scale.
