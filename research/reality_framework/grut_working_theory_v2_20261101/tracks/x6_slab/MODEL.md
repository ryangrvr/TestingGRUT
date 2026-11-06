# Slab action and transport overlap

## Spatial state and reciprocal work

The scalar has a finite Fourier band along a periodic comoving $x$
coordinate. Its mass phase $\theta=\phi/f$ multiplies the fermion
mass. For a transverse momentum of magnitude $q_\perp$, align it
with one transverse spin axis and restore transverse isotropy by
averaging azimuth. For the finite longitudinal plane-wave projector
$G$, the Hamiltonian is

$$
H=
\begin{pmatrix}
(\sigma_z q_x+\sigma_x q_\perp)/a & M G e^{i\theta}G\\
M G e^{-i\theta}G & -(\sigma_z q_x+\sigma_x q_\perp)/a
\end{pmatrix}.
$$

The matrices on the off-diagonal act as the identity on spin.
Longitudinal Fourier truncation, scalar bandwidth, real-space
quadrature and transverse bin spacing are independent settings.
The radial annular weight per length is
$W_\ell=\int_\ell q_\perp\,dq_\perp/(2\pi L)$.
Evolve every occupied column of $\Psi$ by the same unitary operator.
Then $C=\Psi\Psi^\dagger$ stays a projector. The Majorana trace is
half the four-component trace. For a constant phase the spin
helicities split into exactly the two X5 Hamiltonians on this
quadrature.

The uniform reference $C_0(a)$ is the negative projector of
$H_0(a)$ at $\theta=0$. The normal-ordered comoving energy and
physical pressure components use

$$
E_F=\frac12\sum_\ell W_\ell\,
{\rm Tr}(HC-H_0C_0),\qquad
P_{F,i}=\frac1{2a^3}\sum_\ell W_\ell\,{\rm Tr}(K_i(C-C_0)).
$$

Here $K_i$ is the kinetic contribution in direction $i$. The
azimuth-averaged pressure in each transverse direction is half
the transverse kinetic trace. The longitudinal momentum density
has the same measure and uses $q_x/a$ in place of $K_i$.
The force density is the field variation of $E_F/a^3$ at fixed
state. Differentiating the same energy gives

$$
\frac{dE_F}{dt}+\frac{P_{F,x}+2P_{F,\perp}}3
\frac{d a^3}{dt}=a^3\langle J\dot\phi\rangle.
$$

The scalar equation carries the opposite force. The force and
stress both use the same orbitals. The pressure of a freely
projected particle distribution can differ from this pressure.

## Finite matching and expansion

The homogeneous temporal and metric coefficients are those of
the X5 second-order expansion on the slab quadrature,

$$
Z_t=\frac{M^2}{4f^2a^3}\sum_{\ell,n}\frac{W_\ell}{E_{\ell n}^3},
\qquad
Z_g=\frac{M^2}{4a^3}\sum_{\ell,n}\frac{W_\ell p^2}{E_{\ell n}^5}.
$$

Match a spatial coefficient $Z_x$ to the curvature of the
instantaneous static vacuum energy for
$\phi=\epsilon\cos(2\pi x/L)$. The scalar action uses
$A=1-Z_t$ in the kinetic term and $B=1-Z_x$ in the gradient term.
Both coefficients are kept in the same force and metric variation.
With primes for derivatives in $\ln a$, the corresponding scalar
pressure is

$$
P_\phi=\frac A2\langle\dot\phi^2\rangle
       -\langle V\rangle
       +\frac{A'}6\langle\dot\phi^2\rangle
       -\frac{B+B'}{6a^2}\langle(\partial_x\phi)^2\rangle .
$$

The metric is the isotropic minisuperspace reduction,
$H^2=(\rho_\phi+\rho_F)/(3M_P^2+Z_g/2)$.
The pressure average is varied in this reduction. The directional
scalar diagnostic treats the kinetic coefficients as functions of
the volume scale. Shear and spatial metric perturbations require
additional equations.

The matched fermion diagnostic reallocates
$\frac12 Z_t\langle\dot\phi^2\rangle+
\frac{Z_x}{2a^2}\langle(\partial_x\phi)^2\rangle$ from the scalar
subtraction into the fermion ledger. This leaves the total
unchanged. Its scalar-work term is obtained using

$$
J_{\rm matched}
=J-Z_t\ddot\phi-H(3Z_t+Z_t')\dot\phi+Z_x a^{-2}\partial_x^2\phi.
$$

The extra fermion pressure terms are the corresponding variation
of that subtracted energy. The gravity coefficient stays in the
constraint. The two separate pressure-work identities test this
allocation. A finite-cylinder prescription has a smaller domain
than a covariantly matched continuum action.

## Integration and preparation

A symmetric fixed-scale step splits off the uniform fermion
Hamiltonian and the free scalar oscillator, including the
matched gradient. The residual mass Hamiltonian $B_F$ is held
fixed for a middle step. Its eigenvalues are $\epsilon_i$.
In that basis the exact integrated covariance for the scalar
impulse is

$$
\int_0^{\Delta t} C_{ij}(t)\,dt
=C_{ij}(0)\Delta t\,
e^{-i(\epsilon_i-\epsilon_j)\Delta t/2}
\,{\rm sinc}\!\left((\epsilon_i-\epsilon_j)\Delta t/2\right).
$$

Vary the same residual mass matrix with the resulting covariance
to obtain the opposite scalar momentum impulse. The uniform
vacuum propagates exactly in this split when its source and
metric are static. Outer expansion substeps hold the scalar
canonical momenta and spin state fixed. They integrate the
pressure work independently of the matter integrator.

At the initial slice the first adiabatic graph from occupied to
positive states is $D_{1,+-}=-i\dot H_{+-}/\Delta_{+-}^2$.
The separate second-order preparation uses
$D_{2,+-}=i(U_+^\dagger\dot X U_-)_{+-}/\Delta_{+-}$ with
$X=U_+D_1U_-^\dagger+{\rm h.c.}$.
The time derivative uses a short centered difference on the
initial tangent and acceleration of the scalar and metric.
Normalize the graph so the initialized state is an exact
projector. The rotating-vacuum control tests the derivative
order. This choice of preparation is an input to the model.

## The sharp longitudinal edge

For a spatial harmonic of index $r$, negative and positive
vacuum states mix between $n$ and $n+r$. Let
$E_n=\sqrt{M^2+\mathbf p_n^2}$ and
$S_{nn'}=(1+(M^2+\mathbf p_n\cdot\mathbf p_{n'})/(E_nE_{n'}))/(E_n+E_{n'})$.
The response has the following exact decomposition.

$$
\begin{aligned}
\sum_{n=-N}^N\frac1{E_n}-\sum_{n=-N}^{N-r}S_{n,n+r}
={}&\frac12\sum_{j=0}^{r-1}\left(\frac1{E_{-N+j}}+\frac1{E_{N-j}}\right)\\
 &+\sum_{n=-N}^{N-r}\frac{(k_r/a)^2}{2E_nE_{n+r}(E_n+E_{n+r})}.
\end{aligned}
$$

The first line counts the missing virtual neighbors. The second
is the interior response. Multiplication by
$M^2 W_\ell/(f^2 a k_r^2)$ and summation over transverse bins
gives the spatial coefficient for that harmonic. At fixed aspect
ratio and large cutoff, the edge contribution to the static
curvature scales as $|k_r|$ and grows with the cutoff. A single
local gradient coefficient gives curvature proportional to
$k_r^2$. Matching it at one harmonic therefore leaves an edge
term at other harmonics. This diagnoses the static regulator
artifact. The time evolution measures the residual cutoff
sensitivity directly.

## Exact overlap and a tested transport boundary

Let $C^b=U^\dagger C U$ in a moving energy basis, and
$\mathcal A=iU^\dagger\dot U$. The exact collisionless overlap is

$$
i\dot C^b=[E-\mathcal A,C^b].
$$

The diagonal and off-diagonal blocks are parts of this same
state. While the basis moves, the off-diagonal blocks can feed
occupation and work. The fixed-volume transport control uses
a compact driven phase, records its injected energy and momentum,
and lets the entire covariance enter a stationary free tail.

The tail specifies a parent rest decay width $\Gamma_0$ and an
optional elastic dephasing channel. The latter is a convex
random-unitary channel generated by the free band sign. It damps
opposite-band coherence without changing occupations or the
conserved energy and momentum. In the free eigenbasis, if
$\Delta C=C-C_0$,

$$
\Delta C_{ij}(t)=\Delta C_{ij}(0)
e^{-i(E_i-E_j)t-(\Gamma_i+\Gamma_j)t/2}
e^{-\gamma_\phi t\,{\bf1}_{\text{opposite bands}}},\qquad
\Gamma_i=\Gamma_0 M/|E_i|.
$$

The affine amplitude-damping map has
$C(t)=D C(0)D^\dagger+(I-DD^\dagger)C_0$ before band dephasing.
It preserves $0\leq C\leq I$. The decay daughters receive the
energy and momentum of the lost parent occupations. Their
directional stress is calculated by boosting isotropic rest-frame
massless daughter pairs. Opposite-band coherence is retained in
the parent pressure until its omission is certified.

Volume averaging removes different-momentum matrix elements of
free stress. At each momentum the bound on the remaining
off-diagonal pressure is a sum of absolute coefficients times
$e^{-(\Gamma+\gamma_\phi)t}$. It does not rely on catching a
zero of an oscillation. A passing cut preserves the actual
occupations and the previously accumulated daughter state.
Local same-band spatial coherence can persist. The controlled
comparison certifies averaged stress for its stated channel.
Deriving a cosmological collision operator requires separate
physics. [Kinetic transport with coherence](https://arxiv.org/abs/0811.0936)
