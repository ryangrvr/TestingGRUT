# X10 · A spatial source, joint preparation, and a transport boundary

## Hamiltonian and preparation

The finite control uses X8's unitary mass operator in one
longitudinal Galerkin band. The transverse momentum is integrated
in radial shells. At each retained momentum and helicity, X9's
two chiral bath channels carry signed frequencies. Write
$S$ for the sterile subspace and $B$ for the bath. In physical
time the finite Hamiltonian is

$$
 h(a,\phi)=a^{-1}K_{\rm com}+m_0(a) V(\phi),\qquad
 K_{\rm com}=\begin{pmatrix}K_S&G^\dagger\\G&D_B\end{pmatrix}.
$$

$D_B$ contains the conformal bath frequencies. $G$ is the
square root of the quadrature weight times the signed
vacuum spectrum. The real sterile mass is matched on the
zero-drive shell as in X9. $V$ acts only on $S$ and has the
projected unitary off-diagonal phase. Spatially uniform
$\phi$ reduces to the chiral rotation of the homogeneous
sterile and bath system. The physical bath closure here is
a quadratic representation of the vacuum Born spectrum.
Later Standard Model scattering and a finite-density bath
are separate dynamics.
[X8 mass operator](../x8_unitary/MODEL.md) · [X9 bath spectrum](../x9_memory/MODEL.md)

The initial covariance is the negative-energy projector of
the entire coupled Hamiltonian at $a=1,\phi=0$. It has
nonzero sterile–bath correlations even before any drive.
Keep $G$ present throughout the preparation. With $u=(t+1.5)/1.5$,
the prescribed profile for $-1.5\le t\le0$ is

$$
 \phi(x,t)=(10u^3-15u^4+6u^5)
           [0.15\sin x+t(0.9\cos x+0.30)] .
$$

The final homogeneous velocity breaks the half-period
symmetry and gives a channel difference to test. The
covariance is carried intact into reciprocal evolution
at release. The scalar profile and prior drive are choices
for the finite control. The physical H55 probe uses the
archived seed and an independently specified initial
adiabatic state below.
[Initially coupled quantum systems](https://doi.org/10.1088/2058-9565/ad98be)

The propagation evolves weighted frames $C=OO^\dagger$.
The fixed-field residual mass kick uses the Fréchet
derivative of the matrix exponential. It integrates the
same covariance that receives the quantum kick. The
opposite force updates the scalar canonical momentum.
The uniform Hamiltonian advances by a Schur-complement
Cayley step, which is unitary and preserves its static
vacuum. An independent dense solver checks the split.

## Energy, pressure and momentum

Let $F_0(a)$ be the weighted filled-vacuum energy of
$h(a,0)$, with the common Majorana factor accounted for
below. For a specified covariance the weighted raw
expectations of $K_{\rm com}$ and $V$ are $K$ and $M_\phi$.
Then

$$
 \rho_Q = \frac{K/a+m_0 M_\phi-F_0}{2a^3},\qquad
 P_Q = \frac{K/a-m_0' M_\phi+F_0'}{6a^3},
 \quad {}'=\frac{d}{d\log a}.
$$

The volume derivative is at fixed scalar canonical data and
covariance. Bath, sterile–bath interaction, and vacuum terms
enter the same derivative. Tabulation of $F_0$ and its
Hellmann–Feynman derivative keeps the subtraction consistent.
The transverse angular average here defines a scalar
pressure. A full directional stress tensor for the
massless daughter pair needs its internal momenta.

The finite zero-frequency temporal matching is
$Z_t=a^{-3}\sum W |\langle u|\partial_\phi h|o\rangle|^2/\Delta^3$
over occupied and empty states of $h(a,0)$. The same
coefficient is used in the scalar gradient term. Replacing
the vertex by $\partial_{\log a}h$ defines $Z_g$. The
classical metric constraint uses
$H^2=\rho/(3M_P^2+Z_g/2)$.
The scalar pressure includes the derivative of its
coefficient $A=1-Z_t$. The longitudinal Noether charge is
the covariance momentum minus the scalar momentum-gradient
contraction. The replay checks it alongside $a^3\rho+\int
P\,d(a^3)$ and the separate reciprocal quantum work.

## Compression and continuation

The transport tests use pinching in orthogonal subspaces.
For projectors $P_j$, $C'=\sum_j P_j C P_j$ preserves
$0\le C'\le I$. The bare occupation map is diagonal in the
free sterile and bath bases. The dressed occupation map
is diagonal in the instantaneous interacting Hamiltonian.
It preserves instantaneous energy. The covariance of a
separated sterile block can retain scalar work at the
boundary, while removal of the sterile–bath block changes
the interaction energy and channel currents.

For each signed bath leg with spectral matrix $D_\ell$
and coupling vector $g_\ell$, retain the span of

$$
 g_\ell,\ D_\ell g_\ell,\ \ldots,\ D_\ell^{r_\ell-1}g_\ell .
$$

The coherent block is the whole sterile space and all
retained directions. Diagonalize the free bath
Hamiltonian in each remaining leg subspace and keep its
populations. All terms in the sterile mass force and
sterile–bath interaction are preserved when each
$r_\ell\ge1$. So are the signed channel counts and the
longitudinal momentum. The only jump in the quantum
energy is the omitted covariance across the remaining
bath-energy boundary. Its signed value per leg is
$-W_\ell\Re\operatorname{tr}(D_{A Q}C_{Q A})/a^4$
in energy density units. The pressure jump is one
third of that value. The adaptive rule carries the
smallest number of actual directions that gives each
leg an equal share of half the permitted energy jump.
The plus-one rule keeps a further direction per leg.
Both use the measured covariance, without adjusting
any particle population.

There is a local derivative guarantee. If each leg
retains at least $r$ directions, $h^jS$ lies inside
the retained block for $j\le r$. Repeated commutators
of a sterile observable therefore keep their
expectations through derivative order $r$ at frozen
$h$. The first commutator of a channel population
has support on the sterile space and the first
coupling direction. The same induction applies to
its derivatives. This guarantee is local in time.
For downstream predictions the code replays the
pinched state with both the full dynamics and a
hybrid closure.

In the hybrid closure the active block stays coherent.
The passive sector contains fixed comoving occupations
in its projected bath eigenmodes. The Hamiltonian
terms between the two blocks are dropped. Their
expectation vanishes in the pinched state, so the
switch of generator adds zero energy or pressure.
The same matching table is kept. At later times the
missing term is the response of the passive bath to
the retained block, of the form
$D_{A Q}e^{-iD_{Q Q}(\eta-\eta')}D_{Q A}$ in conformal
time. Moving the boundary tests whether the chosen
active space captures enough of this response for
the recorded observables.
[Reaction-coordinate transformations](https://doi.org/10.1103/PhysRevB.97.205405)

## Physical H55 axial preparation and feedback

The physical input check regenerates the original X2
scalar vacuum seed. An axial projection takes the
plane average perpendicular to one axis and retains
its first Fourier harmonic. The homogeneous scale and
mean fields use X9's early H55 history. A single
sterile shell has the physical N1 mass, phase scale
and width. Its initial covariance uses the coupled
negative spectral subspace at that scalar profile
and its first adiabatic graph correction. For
positive and negative energies separated by $\Delta$,
the off-diagonal graph is
$-i\langle u|\dot h|o\rangle/\Delta^2$.
The derivative uses the prescribed H55 velocities
and scale at the reconstructed slice. The bath
remains coupled during evolution.

On the fixed metric the three axial scalar
coefficients $z$ use the Hessian
$A_{\rm ax}=k^2/a^2+V''(\overline{\boldsymbol\phi})$.
Their canonical momenta are $p=a^3\dot z$. The
feedback test applies the exact first Fourier
coefficient of the N1 mass kick to $p_\sigma$.
The quadratic axial energy is
$a^{-3}p^\dagger p+a^3 z^\dagger A_{\rm ax}z$.
The external ledger includes the time variation
of the metric and Hessian and the work of the
prescribed homogeneous fields. A simultaneously
evolved homogeneous quantum reference removes
the common vacuum energy. The axial work then
cancels between the scalar coefficient and the
quantum difference. This probe covers one early
spatial Fourier mode and one momentum shell,
with the H55 mean and metric prescribed.
