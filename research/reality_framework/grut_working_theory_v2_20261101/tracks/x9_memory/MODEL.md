# X9 · State memory in H55 production and a continuous vacuum bath

## The early phase family

Use the X5 homogeneous chiral mode with Hamiltonian
$h_j=M_j\sigma_x+(k+\dot\phi_j/(2f_j))\sigma_z$ and
$k=hq/a$ for helicity $h$. Its Bloch vector is
$\mathbf r=(r_x,r_y,r_z)$. The physical free occupation is
$f_0=(1+(M_j r_x+k r_z)/E)/2$ with $E=\sqrt{M_j^2+k^2}$.
The free-mass stress uses $r_x,r_z$ and the scalar work
uses $r_y$. At fixed scalar and metric data, a phase
change of $r_y$ can leave the energy, pressure and
occupations identical. [X7 moment identities](../x7_regulator/MODEL.md)

X9 evolves from the inherited slow-roll state and forks at
$t=0.16H_i^{-1}$. Write $F(p)$ for the compact momentum
weight. It is one up to $p=3000$ and smoothly declines to zero
at $p=6000$, in $H_i$ units. The evolved reflection sets
$r_y^{\rm new}=(1-2F)r_y$. The new covariance lies in the
Bloch ball. Inside the window it is an opposite pure phase
up to numerical roundoff. The interpolating region can
be mixed. Each fork evolves its own scalar force and expansion.

For a fixed scalar trajectory the mode response is linear.
Let $R_i(t,t_*)$ be the orthogonal Bloch rotation supplied
by the same midpoint Hamiltonian on that trajectory. For
each mode $i$ the allowed component before tapering satisfies
$|y_i|\le y_{\max,i}=\sqrt{1-r_{x,i}^2-r_{z,i}^2}$.
The preparation interpolates between $y_i$ and the
reference component with weight $F_i$. For any linear
observable with time-dependent coefficient $B_i(t)$,
define the final-window response coefficient

$$
 C_i=\int dt\,g(t) B_i(t)\cdot R_i(t,t_*)\hat y ,
$$

where $g(t)$ is the normalized average window or a delta
at the endpoint. The largest possible pairwise difference is

$$
 2\sum_i F_i y_{\max,i}|C_i| .
$$

Independent choices of the endpoint signs attain this
maximum. The exact discrete response uses the velocities
and scale at the midpoint of each actual X5 matter step.
The bound applies with that history held fixed. The
coevolved pair is the separate test of nonlinear feedback.

## A quadratic representation of the vacuum self-energy

At leading order in the Yukawa coupling the massless
lepton–Higgs continuum supplies a chiral spectral
density. The positive-frequency density follows the X8
two-body phase-space calculation, and the negative leg
is fixed by vacuum charge conjugation. In conformal time
use $k=hq$, chiral sign $b=\pm1$, and
$\gamma=\Gamma_{\rm rest}/(2\pi M)$. The regulated density is

$$
 J_b(\omega;q)=\gamma\,|\omega-bk|\,\mathbf1_{|\omega|>q}
              R\big((\omega^2-q^2)/\Lambda^2\big).
$$

The taper is one up to invariant argument $1/2$, falls
with a cubic smoothstep, and reaches zero at one.
The bath frequencies are signed. In the Gaussian
representation the occupied negative spectral subspace
is the interacting vacuum at the initial laboratory
Hamiltonian. A sterile excitation is a positive spectral
addition, and its Pauli bound follows from orthogonality
to that vacuum. This is a quadratic representation of
the vacuum Born spectral function. Full lepton–Higgs
scattering, medium effects, flavor correlations and
higher bath correlations have separate dynamics.
[X8 vacuum kernel](../x8_unitary/MODEL.md) ·
[Fermionic environments and two-point correlations](https://doi.org/10.1103/PhysRevResearch.5.033011)

The finite quadrature Hamiltonian at scale $a$ has a
sterile block $m_0(a)\sigma_x+k\sigma_z$, bath diagonal
$\omega_\ell$, and a real coupling
$g_{b\ell}=\sqrt{J_b(\omega_\ell)\Delta\omega_\ell}$
between sterile chirality $b$ and its bath channel.
Let $Z$ have eigenvalues $b$ on both parts of each
channel. The rotating-frame Hamiltonian is

$$
 h_{\rm rot}=h_{\rm lab}+\delta Z,\qquad
 \delta=\theta'/2,\qquad \theta=\phi/f .
$$

Here $h_{\rm lab}$ denotes the physical energy matrix
expressed in the rotating basis. The Yukawa couplings
are present for the entire evolution. Eliminating a bath
amplitude gives

$$
 b_{b\ell}(\eta)=
 e^{-i\omega_\ell(\eta-\eta_0)-ib[\theta(\eta)-\theta(\eta_0)]/2}
 b_{b\ell}(\eta_0)
 -i g_{b\ell}\int_{\eta_0}^{\eta}\!d\eta'\,
 e^{-i\omega_\ell(\eta-\eta')-ib[\theta(\eta)-\theta(\eta')]/2}
 s_b(\eta').
$$

The second term produces the retarded convolution in
the sterile equation. The first carries the correlations
in the initial bath amplitudes. A reduced sterile
covariance alone cannot determine their interference
with the sterile source. The finite Hermitian model
keeps both. [Continuous-environment unitary mappings](https://doi.org/10.1088/1367-2630/ab5ec5)

The two signed legs combine at zero drive into
$\Sigma_b(z)=\gamma(z-bk)I((z^2-q^2)/\Lambda^2)$ with
$I(y)=\int_0^1 dx\,R(x)/(y-x)$. The code matches the
dispersive mass on the physical zero-drive shell by
$m_0=aM[1-\gamma\,\Re I((aM/\Lambda)^2+i0)]$.
The retarded imaginary part recovers the Born
time-dilated width. The exact pole includes higher-order
residue corrections of this quadratic model. This
prescription matches the real shell mass. The tested
cutoff variation records the residual finite-band and
field-strength sensitivity.

At a fixed drive with $|\delta|<q$, the common central
gap has $|z|<q-|\delta|$. Both shifted spectral arguments
are spacelike in the gap and both real factors
$1-\gamma I$ are positive. The inverse sterile
propagator determinant there is

$$
 [z^2-(k+\delta)^2]\prod_{b=\pm1}
 \left[1-\gamma I\!\left(\frac{(z-b\delta)^2-q^2}{\Lambda^2}\right)\right]
 -m_0^2<0 .
$$

Thus there is no sterile pole in that common central
gap. Band edges and recurrence of a finite quadrature
still matter for a late-time test.

## Matched joint states and the energy ledger

Construct real orthonormal positive spectral packets
$s$ and $w$. The first is the positive projection of
a free sterile spinor. The second starts from a bath
packet centered on the free mass shell and is
orthogonalized against both projected sterile
coordinates, so $w$ has zero sterile component. The two
excitation states are

$$
 |\psi_\pm\rangle=(|s\rangle\pm i|w\rangle)/\sqrt2 .
$$

They have the same sterile reduced covariance, bath
energy-bin occupations, total number and total energy.
Their bath coherences and sterile–bath correlations can
differ. For any Hermitian quadratic observable $O$ their
difference at time $t$ is exactly
$-2\Im\langle s_t|O|w_t\rangle$ on a common trajectory.
The common filled vacuum cancels in this difference.
The driven test therefore propagates the two added
excitations with the same Hamiltonian and records
changes per matched excitation. It is separate from
integrating the physical H55 distribution.

The finite system tracks sterile, bath and interaction
energy. Since the interaction commutes with $Z$,
the comoving laboratory energy in a prescribed
background obeys

$$
 E'_{\rm lab}=m_0'\langle\sigma_x\rangle
               -2\delta m_0\langle\sigma_y\rangle .
$$

The midpoint rotating step can also track exactly
$\Delta E_{\rm lab}=
\Delta m_0(\langle\sigma_x\rangle_i+\langle\sigma_x\rangle_f)/2
-\delta_{\rm mid}\Delta\langle Z\rangle$.
The Padé evolution preserves the midpoint Hamiltonian.
The continuous torque integral and a dense matrix
exponential are independent checks.

In the finite live scalar control, $\delta=\bar v/(2f)$
uses the self-consistent midpoint scalar velocity.
The scalar impulse is
$\Delta P=\Delta\langle Z\rangle/(2f)-dt\,\overline{V'}$.
The first term exchanges energy with the Hermitian
quantum system. The second uses the discrete potential
gradient. This conserves total energy and tracks the
shift charge including the potential term. The filled
vacuum and the added excitation are all evolved in
this control.

## What dephasing can erase

Under a stationary real Hamiltonian, conjugate
preparations have the same weight in every real
energy spectral projector. The infinite-time average
of a real quadratic observable is its expectation
after projection onto energy eigenspaces, so those
averages agree exactly in a finite stationary model.
Before recurrence, the finite bath also approximates
escape into outgoing continuum channels. A sterile
count can lose its phase dependence while a difference
survives in the channel-resolved outgoing occupations.
The two orders of limits concern different observables.
Neither determines the result with a scalar still
doing work, or with later interactions between
daughter channels. The numerical controls test those
distinctions without interpreting an auxiliary
chirality count as an observed baryon abundance.
