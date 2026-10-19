# Fixed-charge rotating one-loop energy: derivation

This is the original Candidate 1 scalar action and MSbar matching, not a new
potential, parameter fit or change of the gravitational vacuum parameter. All
coordinates below are dimensionless, `rho=gv r`, `tau=gv t`. The classical
action has an overall `1/g^2`; one scalar-loop energy is measured in `gv`, while
the classical energy is `(v/g) Ebar`. The ratio is essential to loop counting.

## 1. Background and full operator

Write `sqrt(2) X/v=exp(-i Omega tau)(F+u+iw)`, `phi/v=S+eta`.
The background is the hashed original solution at
`Omega=0.8000000019763136`, with reporter background zero. After replacing each
radial fluctuation by `rho` times its value, the measure is `d rho`.

For `q=(u,w,eta)` in that transformed measure,

\[
\ddot q+G\dot q+K_\ell q=0,\quad
G=\begin{pmatrix}0&2\Omega&0\\-2\Omega&0&0\\0&0&0\end{pmatrix},
\]

\[
K_\ell=\begin{pmatrix}
-\partial_\rho^2+\ell(\ell+1)/\rho^2+S^2+.3F^2-\Omega^2&0&2FS\\
0&-\partial_\rho^2+\ell(\ell+1)/\rho^2+S^2+.1F^2-\Omega^2&0\\
2FS&0&-\partial_\rho^2+\ell(\ell+1)/\rho^2+3S^2-1+F^2
\end{pmatrix}.
\]

The reporter block is `-partial_rho^2+ell(ell+1)/rho^2+1+.01 S^2`.
For a Euclidean Fourier convention the pencil is
`P_l(zeta)=K_l+zeta^2+zeta G`. Reversing the Fourier sign exchanges the two
conjugate charged branches and leaves its determinant unchanged. Their free
inverse lengths are

\[
\kappa_\pm=\sqrt{1+(\zeta\pm i\Omega)^2},\quad
\kappa_\eta=\sqrt{2+\zeta^2},\quad
\kappa_D=\sqrt{1.01+\zeta^2},\quad \Re\kappa>0.
\]

The full gyroscopic operator has not acquired an additional rank-one term.

## 2. Fixed charge, not a naive grand-canonical Gaussian

The unfixed-charge static curvature contains a negative amplitude direction.
One must not define this calculation by integrating a real Gaussian over that
unconstrained direction. Let `e=F/sqrt(I)` be the normalized phase zero vector
in the original three-dimensional inner product, and split
`w=w_perp+theta e`. Equivalently use `rho F` in radial measure. Up to a total
time derivative, the part involving the parallel phase is

\[
L_\theta=\tfrac12\dot\theta^2-2\Omega\langle e,u\rangle\dot\theta,
\quad p_\theta=\dot\theta-2\Omega\langle e,u\rangle=-Q_1/\sqrt I.
\]

At `Q1=0`, Routh reduction adds `4 Omega^2 |e><e|` to the **reduced amplitude
stiffness**, and projects the remaining gyro coupling onto `w_perp`.
Eliminating the parallel phase in the full Euclidean pencil gives, in any
orthonormal finite kinetic basis,

\[
\det P_{\rm full}(\zeta)=\zeta^2\det P_{Q_1=0}(\zeta).
\]

The coefficient is one. There is no extra `I` or `dQ/dOmega` factor in this
identity. Retain the `zeta^2` when normalizing by the **original complete
vacuum** determinant. Otherwise one deletes a vacuum oscillator and changes
the ultraviolet degree. This is a charge-projected/rest-frame semiclassical
construction, conditional on stability of its reduced continuum Hamiltonian;
the previous finite-box stability tests are numerical evidence, not a proof of
that continuum premise.

`theory_toy.py` verifies the determinant polynomial identity using exact rational
arithmetic, reconstructs the full nonzero charge-compatible modes, and compares
the complete-vacuum-normalized logdet integral with half the frequency-sum
difference. They agree to `4.44e-16`. A second toy retains negative unfixed
curvature and positive reduced curvature, explicitly testing why the naive
unfixed Gaussian is not the prescription. This is an algebra/normalization
test, not a continuum theorem or a second physical model used for the result.

Translations are collective coordinates with rest momentum zero. Their
normalization and the phase/charge collective normalization can change
subextensive measure factors, but not the coefficient of the long-time energy.
Their zero-frequency logarithms must still be retained in the determinant.

## 3. Definition of the relative one-loop Routhian

For the stable reduced saddle and the same scalar vacuum,

\[
\frac{K_1^{\rm ren}}{gv}=
\frac1{2\pi}\int_0^\infty d\zeta\sum_{\ell=0}^\infty(2\ell+1)
\left[\log\frac{\det P_\ell}{\det P_{\ell,0}}
-B_{1\ell}-B_{2\ell}\right]
+\frac{K_{1,\overline{\rm MS}}^{(1)}+K_{1,\overline{\rm MS}}^{(2)}}{gv}.
\]

The reporter determinant is included. `B1=Tr(G0 V)` and
`B2=-Tr(G0 V G0 V)/2` are the first two coefficients of **the same** continuum
determinant. Subtracting both makes the remainder ultraviolet integrable in
three spatial dimensions. Restoring the two graphs with MSbar counterterms at
`mu0=gv` implements the original matching. It is not subtraction at an
adjustable on-shell normalization point. The counterterm constants are not
fitted to the soliton energy.

The Born subtraction/restoration principle follows the continuum
renormalization method described in [Graham, hep-th/9912286](https://arxiv.org/pdf/hep-th/9912286).
The single-potential approximation and finite on-shell counterterms in
[hep-th/0105009](https://arxiv.org/pdf/hep-th/0105009) are not imported into this
coupled rotating model.

## 4. Restored graphs, including the laboratory harmonics

Relative to vacuum squared masses `(2,1,1,1.01)`, define

\[
A=3(S^2-1)+F^2,\quad B=S^2-1+.2F^2,\quad C=.1F^2,
\quad D=2FS,\quad E=.01(S^2-1).
\]

The laboratory charged block is `B I+C M(2 Omega tau)` with a traceless rotating
matrix `M`, while the eta/charged vector is `D(cos Omega tau,-sin Omega tau)`.
Thus the two nonstatic invariants are `4 Omega^2-k^2` and `Omega^2-k^2`, not
their static limits or negative frequency squares. Define

\[
\mathcal A(t)=\frac{t(\log(t/\mu^2)-1)}{16\pi^2},\quad
B_0(s;a,b)=-\int_0^1dx\log\frac{xa+(1-x)b-x(1-x)s}{\mu^2},
\quad \widetilde V(k)=4\pi\int d\rho\,\rho^2\frac{\sin k\rho}{k\rho}V(\rho).
\]

In the energy units used here, `mu=1`, and

\[
K^{(1)}=\tfrac12\left[\mathcal A(2)\int A+2\mathcal A(1)\int B+\mathcal A(1.01)\int E\right],
\]

\[
\begin{split}
K^{(2)}=-\frac1{64\pi^2}\int_0^\infty\frac{k^2dk}{2\pi^2}\big[&
\widetilde A^2 B_0(-k^2;2,2)+2\widetilde B^2 B_0(-k^2;1,1)\\
&+2\widetilde C^2 B_0(4\Omega^2-k^2;1,1)
+2\widetilde D^2 B_0(\Omega^2-k^2;2,1)
+\widetilde E^2 B_0(-k^2;1.01,1.01)\big].
\end{split}
\]

Both factors of two arise from actual species/harmonic traces, not an added
polarization assumption. The kernels stay below their pair cuts. The local
second-order contraction is
`tr(V^2)=A^2+2B^2+2C^2+2D^2+E^2`, as required for the scalar potential
counterterms. `graphs.py` checks that identity, all five Parseval transforms,
refinement, pair cuts and the explicit MSbar scale derivative. A scale
derivative of these restored graphs alone is not a full RG-invariance test.

## 5. Infrared and fixed-charge Legendre conversion

The phase sector in `ell=0` and each translation sector in `ell=1` contributes
`2 log zeta` near zero. For numerical integration subtract
`log[zeta^2/(zeta^2+1)]` and restore its exact integral up to the same cutoff Z,

\[
Z\log\frac{Z^2}{Z^2+1}-2\arctan Z.
\]

The multiplicities are one and three. Over the full frequency axis each
restored factor contributes `-1/2` energy units. Omitting it would shift the
answer by `+1/2` per zero mode. No arbitrary primed determinant or fitted
zero-mode normalization replaces this step.

Write `K(Omega)=E-Omega Q=K0+K1+...`. At fixed reference Q, stationarity gives
`K0'(Omega0)+Q=0`, so the O(hbar) frequency shift cancels out of the energy:

\[
E_1(Q)=K_1(\Omega_0).
\]

The reference `Qref=39943.84624223011` is a continuous semiclassical branch
point, not an integer charge eigenvalue. The predeclared exact-charge target is
the nearest integer `N=39944`. To the calculated order,

\[
\frac{E_{\rm sol}^{0+1}(N)}{gv}=
\frac{\bar E_0}{g^2}+\Omega_0(N-Q_{\rm ref})+
\frac{K_1^{\rm ren}(\Omega_0)}{gv}.
\]

The next classical charge-displacement term is
`g^2 (N-Qref)^2/(2 dQbar/dOmega)`; it is checked separately. The shift of K1
under this O(g^2) displacement is beyond strict one-loop order.

Use the already calculated **strict one-loop pole mass** for a like-order
free-carrier threshold, not a curvature mass or a square root that resums
selected higher-order terms:

\[
\Delta E_{\rm freeX}=E_{\rm sol}^{0+1}(N)-N m_X^{0+1}.
\]

This is one specified breakup channel. Its sign is not a calculation of
single-particle emission leaving a residual soliton, fission, tunneling,
annihilation rates, or an all-order lifetime. All energies subtract the same
scalar vacuum; this does not retune the gravitational vacuum density.
