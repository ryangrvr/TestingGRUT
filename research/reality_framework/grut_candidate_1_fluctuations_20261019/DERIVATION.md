# Rotating finite-charge scalar fluctuations

19 October 2026. Same Candidate 1 action, matching and prepared branch; no new field or evolution law. This file derives the **tree-level flat-space scalar** fluctuation problem. Gravitational/Standard Model fluctuations and the renormalized inhomogeneous one-loop determinant are outside its scope. The existing real reporter \(D\) is included as a decoupled quadratic block. Numerical mode searches and the vacuum pole calculation are performed separately; this derivation does not promote a finite eigenvalue search into a complete continuum spectrum.

## 1. Background, units and sign conventions

Use \(\tau=g_Xvt\), \(\mathbf y=g_Xv\mathbf x\), \(A=\sqrt2X/v\), \(S=\phi/v\), and \(\chi=D/v\). With \(a=1,b=0.1\), the background is

\[
A_*=e^{-i\Omega\tau}F(r),\qquad S_*=S_0(r),\qquad\chi_*=0,
\quad \Omega\simeq0.8000000019763.
\]

The dimensionless potential is

\[
V=\frac a4(S^2-1)^2+\frac12S^2|A|^2+\frac b4|A|^4
+\frac12(\mu_D^2+\kappa_D S^2)\chi^2+\frac{b_D}{4}\chi^4,
\]

where \(\mu_D^2=m_{D0}^2/(g_X^2v^2)=1\), \(\kappa_D=h_D^2/g_X^2=0.01\), and \(b_D=\lambda_D/g_X^2\). The matching-scale \(|A|^2\chi^2\) coefficient is zero. The quartic \(b_D\) does not enter the Hessian at \(\chi=0\). This pass uses these original classical coefficients, rather than a partially RG-improved fluctuation operator.

Energy and charge conventions are

\[
E=(v/g_X)\bar H,\quad Q=\bar Q/g_X^2,
\qquad \bar Q=-\int\operatorname{Im}(A^*\dot A)\,d^3y.
\]

Thus \(\bar Q_*=\Omega I\), \(I=\int F^2d^3y\), is positive. Write

\[
A=e^{-i\Omega\tau}(F+u+iw),\qquad S=S_0+z,
\qquad q=(u,w,z,\chi)^T.
\]

All inner products below use the full spatial measure unless a radial measure is explicitly specified.

## 2. Rotating quadratic Hamiltonian and gyroscopic equations

Let

\[
J=\begin{pmatrix}0&1&0&0\\-1&0&0&0\\0&0&0&0\\0&0&0&0\end{pmatrix},
\quad J^T=-J,\quad J^2=-P_X,
\]

where \(P_X\) projects onto the \((u,w)\) components. Multiplication of the complex field by \(-i\) corresponds to this \(J\). If \(B=(F+u,w,S_0+z,\chi)\) denotes rotating coordinates, its scaled canonical momentum is \(p=\dot B+\Omega JB\). The background has \(p_*=(0,-\Omega F,0,0)\). Define the canonical momentum perturbation \(r=p-p_*\). Then

\[
r=\dot q+\Omega Jq.
\]

The time-dependent canonical rotation changes the Hamiltonian to \(\bar K=\bar H-\Omega\bar Q\). The background is a critical point of \(\bar K\), not of \(\bar H\) alone. Let \(H_V\) be \(-\nabla^2\) plus the Hessian of \(V\) at the background. Its quadratic Hamiltonian is

\[
\boxed{K_2=\tfrac12\|r\|^2-\Omega\langle r,Jq\rangle
+\tfrac12\langle q,H_Vq\rangle.}
\]

Set \(L=H_V-\Omega^2P_X\). Completing the square gives the equivalent velocity form

\[
\boxed{K_2=\tfrac12\|\dot q\|^2+\tfrac12\langle q,Lq\rangle.}
\]

The canonical and velocity generators are, respectively,

\[
\partial_\tau\binom q r=
\begin{pmatrix}-\Omega J&I\\-H_V&-\Omega J\end{pmatrix}\binom q r,
\]

\[
\boxed{\ddot q+2\Omega J\dot q+Lq=0,\qquad
\partial_\tau\binom q{\dot q}=
\mathcal A\binom q{\dot q},\quad
\mathcal A=\begin{pmatrix}0&I\\-L&-2\Omega J\end{pmatrix}.}
\]

Since \(J\) is skew, \(dK_2/d\tau=-2\Omega\langle\dot q,J\dot q\rangle=0\). In components,

\[
\ddot u+2\Omega\dot w+L_+u+2FS_0z=0,
\quad \ddot w-2\Omega\dot u+L_-w=0,
\]

\[
\ddot z+L_Sz+2FS_0u=0,\qquad \ddot\chi+L_D\chi=0,
\]

\[
L_+=-\nabla^2+S_0^2+3bF^2-\Omega^2,
\quad L_-=-\nabla^2+S_0^2+bF^2-\Omega^2,
\]

\[
L_S=-\nabla^2+a(3S_0^2-1)+F^2,
\qquad L_D=-\nabla^2+1+0.01S_0^2\ge1.
\]

For \(e^{s\tau}\), the quadratic pencil is \(s^2I+2\Omega sJ+L\). With the convention \(e^{-i\nu\tau}\), it is \(-\nu^2I-2i\Omega\nu J+L\). A positive imaginary generator eigenvalue \(s=i\nu\) is the complex-conjugate frequency convention; stable frequencies occur in the expected conjugate/sign pairs. Static Hessian eigenvalues are not generally squared dynamical frequencies because of the gyroscopic term.

The actual dimensionless action is \(g_X^{-2}\int\mathcal L_*\), so the Hamiltonian generating \(\tau\)-translations in canonical quantum units is \(K_2/g_X^2\). Our formulas use scaled classical momenta/brackets. Quantum commutators of these scaled variables carry \(g_X^2\), or one can rescale perturbations by \(g_X\) to use standard canonical normalization. This bookkeeping does not change the classical frequencies; it matters before a quantum mode sum is attempted.

## 3. Charge constraint and relation to the fixed-charge energy Hessian

The linear charge variation is

\[
\boxed{Q_1=\int(\Omega Fu-Fr_w)\,d^3y
=2\Omega\langle F,u\rangle-\langle F,\dot w\rangle.}
\]

It is conserved by the linear equations: \(\dot Q_1=\langle F,L_-w\rangle=\langle L_-F,w\rangle=0\). Thus every eigenmode with \(s\ne0\) satisfies \(Q_1=0\), while zero and generalized-zero modes require separate treatment. Fixing a conserved global charge is not a new local gauge constraint.

On \(Q_1=0\), decompose

\[
\dot w=\frac{2\Omega\langle F,u\rangle}{I}F+\dot w_\perp,
\qquad\langle F,\dot w_\perp\rangle=0.
\]

Let \(\eta=(u,z)\) and

\[
L_{\rm amp}=\begin{pmatrix}L_+&2FS_0\\2FS_0&L_S\end{pmatrix},
\quad f=(F,0),\quad
L_Q=L_{\rm amp}+\frac{4\Omega^2}{I}|f\rangle\langle f|.
\]

Then the constrained quadratic energy is exactly

\[
\boxed{K_2\big|_{Q_1=0}=\frac12\left[
\|\dot u\|^2+\|\dot z\|^2+\|\dot\chi\|^2+\|\dot w_\perp\|^2
+\langle\eta,L_Q\eta\rangle+\langle w,L_-w\rangle
+\langle\chi,L_D\chi\rangle\right].}
\]

This is the precise connection to the earlier fixed-\(Q\) amplitude Hessian. Omitting the positive rank-one term, or arbitrarily setting \(\dot w=0\) for a spherical amplitude perturbation, solves a different constrained problem. An unconstrained radial negative eigenvalue of \(L_{\rm amp}\) alone is consequently not a demonstrated growing mode of the fixed-charge dynamics.

At second order an exact fixed-charge perturbation includes background corrections. If the full configuration is \(z_*+\epsilon z_1+\epsilon^2z_2+\cdots\), then \(Q_1[z_2]+Q_2[z_1]=0\) and \(H_1[z_2]=\Omega Q_1[z_2]\). The physical order-\(\epsilon^2\) fixed-charge energy is therefore \(H_2[z_1]-\Omega Q_2[z_1]=K_2[z_1]\). Mode-only \(H_2\) and \(Q_2\) need not separately be conserved or equal the complete physical perturbation energy/charge. This distinction is relevant before assigning quantum energy to individual modes.

## 4. Angular reduction and symmetry modes

Use \(\int|Y_{\ell m}|^2d\Omega=1\) and transform a radial perturbation component as \(u=\psi_u(r)Y_{\ell m}/r\), and similarly for the other components. The radial Hilbert measure is \(dr\), with the regular self-adjoint boundary condition at the origin; regular transformed fields behave as \(r^{\ell+1}\). Infinite-volume bound modes are square integrable; scattering states require their separate generalized normalization.

The amplitude matrix is

\[
H_{{\rm amp},\ell}=-\frac{d^2}{dr^2}I_2+
\frac{\ell(\ell+1)}{r^2}I_2+
\begin{pmatrix}
S_0^2+3bF^2-\Omega^2&2FS_0\\
2FS_0&a(3S_0^2-1)+F^2
\end{pmatrix}.
\]

Only \(\ell=0\) receives the charge rank-one term:

\[
H_{Q,0}=H_{{\rm amp},0}+\frac{4\Omega^2}{I}|f_0\rangle\langle f_0|,
\quad f_0=(\sqrt{4\pi}\,rF,0),\quad \|f_0\|^2=I.
\]

For \(\ell>0\), \(\int Y_{\ell m}d\Omega=0\), so \(H_{Q,\ell}=H_{{\rm amp},\ell}\). Every \(\ell\) has multiplicity \(2\ell+1\).

The phase zero mode is \(w\propto F\), since \(L_-F=0\). In generator variables it has \(q_{\rm ph}=(0,F,0,0)\), \(\dot q_{\rm ph}=0\). Translations have \(\eta_i=(\partial_iF,\partial_iS_0)\), \(w=\chi=0\), \(\dot q=0\), and lie in \(\ell=1\). They satisfy the charge constraint because \(\int F\partial_iF=0\).

Along a differentiable background branch,

\[
H_{{\rm amp},0}\eta_\Omega=2\Omega f_0.
\]

In untransformed variables a phase generalized vector is \(q_g=(-\partial_\Omega F,0,-\partial_\Omega S_0,0)\), \(\dot q_g=q_{\rm ph}\), for which \(\mathcal A(q_g,\dot q_g)=(q_{\rm ph},0)\). Its charge is \(Q_1=-d\bar Q/d\Omega\). It is excluded from the fixed-charge tangent space when that derivative is nonzero. A turning point can introduce an additional constrained zero direction. Boosts likewise generate center-motion Jordan chains: a moving soliton can drift linearly along translation modes without internal instability. Global-phase/translation drift must not be counted as exponentially growing internal modes.

## 5. Precise sufficient criterion for all angular momenta

Assume an exact smooth finite-energy radial background with \(F>0\), \(S_0>0\), \(F\to0\), \(S_0\to1\), \(0<\Omega<1\), and bounded sufficiently decaying potential deviations. Consider the regular infinite-volume self-adjoint radial operators. The following are sufficient conditions, rather than conclusions of a finite eigenvalue search.

### Phase and reporter sectors

Let \(g_0=rF>0\). The phase operator satisfies \(H_{-,0}g_0=0\). For compactly supported regular test functions,

\[
\boxed{\langle\psi,H_{-,0}\psi\rangle
=\int_0^\infty g_0^2\left|\left(\frac\psi{g_0}\right)'\right|^2dr\ge0.}
\]

Form closure gives nonnegativity and the single phase kernel \(\operatorname{span}g_0\). Higher \(\ell\) add \(\ell(\ell+1)/r^2\) and cannot introduce a negative direction. The reporter satisfies \(H_{D,\ell}\ge1\), so it introduces no unstable quadratic direction.

### Amplitudes for \(\ell\ge1\)

It suffices to establish \(H_{{\rm amp},1}\ge0\) with kernel exactly the translation vector \((rF',rS_0')\). A useful **analytic sufficient hypothesis** is strict monotonicity \(F'<0,S_0'>0\) for \(r>0\). After the fixed component sign flip \(\psi_u\mapsto-\psi_u\), the off-diagonal potential is \(W=-2FS_0<0\), and

\[
h=(h_1,h_2)=(-rF',rS_0')>0
\]

is a zero solution of the sign-flipped \(\ell=1\) operator. Its coupled ground-state identity is

\[
\boxed{\langle\psi,\widetilde H_{{\rm amp},1}\psi\rangle
=\int_0^\infty\left[
h_1^2\left|\left(\frac{\psi_1}{h_1}\right)'\right|^2
+h_2^2\left|\left(\frac{\psi_2}{h_2}\right)'\right|^2
+(-W)h_1h_2\left|\frac{\psi_1}{h_1}-\frac{\psi_2}{h_2}\right|^2
\right]dr.}
\]

To derive it, substitute \(\psi_i=h_iv_i\), integrate each diagonal kinetic term once, and use \(-h_i''+V_{ii}h_i+Wh_j=0\). Boundary terms vanish first on compactly supported tests and the identity extends by the quadratic-form closure. Every term is nonnegative, and equality forces both ratios to be the same constant. Thus the only kernel is the translation mode in each \(m\) sector.

For every \(\ell\ge2\),

\[
\boxed{H_{{\rm amp},\ell}=H_{{\rm amp},1}
+\frac{\ell(\ell+1)-2}{r^2}I_2.}
\]

The added term is strictly positive on every nonzero radial test vector. Hence no higher angular momentum can introduce a negative or zero amplitude direction under the \(\ell=1\) hypothesis. This is an algebraic all-\(\ell\) implication; it is not an extrapolation from five computed \(\ell\) values. Numerical monotonicity samples and approximate translation overlaps support, but do not by themselves certify, the exact monotonicity/operator hypothesis. If monotonicity is not established, one may instead prove the stated \(\ell=1\) operator condition directly.

### Remaining radial amplitude condition

Require \(H_{Q,0}>0\), with a positive lower bound on its form domain. An alternative conditional criterion uses its local operator \(A_0=H_{{\rm amp},0}\): assume \(A_0\) is invertible and has **exactly one** negative eigenvalue. Along a differentiable branch,

\[
1+\frac{4\Omega^2}{I}\langle f_0,A_0^{-1}f_0\rangle
=\frac{d\bar Q/d\Omega}{I}.
\]

If \(d\bar Q/d\Omega<0\), the positive rank-one update removes that unique negative direction. For clarity, the inertia proof uses the block operator

\[
\mathcal B=\begin{pmatrix}A_0&f_0\\\langle f_0,\cdot\rangle&-1/\beta\end{pmatrix},
\qquad\beta=4\Omega^2/I>0.
\]

Taking the two Schur complements gives

\[
n_-(\mathcal B)=1+n_-(H_{Q,0})
=n_-(A_0)+n_-[-(1+\beta\langle f_0,A_0^{-1}f_0\rangle)/\beta].
\]

When the bracket is negative and \(n_-(A_0)=1\), this yields \(n_-(H_{Q,0})=0\); invertibility excludes a zero direction. The asymptotic positive essential threshold then supplies a gap. At \(d\bar Q/d\Omega=0\), \(H_{Q,0}\eta_\Omega=0\). Neither a negative charge slope alone nor a short list of positive finite-box eigenvalues establishes the needed index and domain hypotheses.

### Coercivity and the no-growth conclusion

Under these conditions, the only constrained energy kernel consists of phase and translations. The asymptotic static thresholds are \(1-\Omega^2>0\), \(2a>0\), and \(1.01>0\). Decaying matrix potentials and the finite-rank charge term leave those essential thresholds unchanged. With no additional zero mode, each relevant lowest sector has a positive gap after removing its symmetry kernel; higher sectors are bounded below by the fixed lower-sector centrifugal comparisons. Ellipticity and bounded potentials then give an \(H^1\)-configuration/\(L^2\)-velocity coercive bound modulo those symmetries.

The constrained energy decomposition also controls the omitted parallel velocity: \(\dot w_\parallel=2\Omega\langle F,u\rangle F/I\), and translations have \(\langle F,u_{\rm trans}\rangle=0\). Thus, for the well-posed linear flow, there is a constant \(c>0\) such that

\[
K_2\ge c\,\operatorname{dist}_{H^1\times L^2}
((q,\dot q),\mathcal Z)^2,
\]

where \(\mathcal Z\) is the phase/translation zero subspace with zero velocities. Since \(K_2\) is conserved, the distance to this subspace stays bounded for all linear times. This excludes linear growth in that quotient norm. It does not require the center or phase themselves to stop drifting.

Equivalently, conservation gives \(\mathcal A^*\mathbb H+\mathbb H\mathcal A=0\), where \(K_2=\langle Z,\mathbb H Z\rangle/2\). An eigenvector with \(\operatorname{Re}s\ne0\) must have \(Q_1=0\), and

\[
2\operatorname{Re}s\,\langle Z,\mathbb H Z\rangle=0.
\]

Positivity on the quotient forces \(Z\) into the symmetry kernel, on which \(\mathcal A Z=0\), contradicting \(s\ne0\). Hence no finite-energy exponentially growing eigenmode exists under the stated hypotheses. Mere nonnegativity without coercivity does not automatically control the usual \(H^1\times L^2\) distance, and this linear argument is not a nonlinear orbital-stability or lifetime theorem.

## 6. Continuum thresholds, finite grids and quantum limits

At spatial infinity, the rotating complex sector obeys \((\nu\pm\Omega)^2=k^2+1\). The generator frequencies are therefore

\[
\boxed{\nu=\pm(\sqrt{k^2+1}\pm\Omega),\qquad
\nu_S=\pm\sqrt{k^2+2a},\qquad
\nu_D=\pm\sqrt{k^2+1.01}.}
\]

At \(\Omega\simeq0.8\), the positive charged thresholds are \(0.2\) and \(1.8\), the real-\(S\) threshold is \(\sqrt2\), and the reporter threshold is \(\sqrt{1.01}\). The number \(\sqrt{1-\Omega^2}\simeq0.6\) is a static-Hessian/tail scale, **not** the smallest dynamical frequency threshold. To express these rotating-frame frequencies in physical units multiply by \(g_Xv\). They are not vacuum particle pole masses: the laboratory charged frequencies include the additional \(\pm\Omega g_Xv\) shift, and the tree vacuum \(X\) particle mass remains \(g_Xv\). Loop-corrected vacuum pole thresholds are a separate calculation.

### A conditional angular cutoff for real subgap point modes

For the \((u,w,z)\) block, separate the constant vacuum potential \(V_{\rm vac}=\operatorname{diag}(1-\Omega^2,1-\Omega^2,2a)\) from the local profile perturbation

\[
\Delta V(r)=\begin{pmatrix}
S_0^2+3bF^2-1&0&2FS_0\\
0&S_0^2+bF^2-1&0\\
2FS_0&0&3a(S_0^2-1)+F^2
\end{pmatrix}.
\]

Suppose a true half-line bound is available:

\[
B_{\rm true}=\sup_{r>0}r^2[-\lambda_{\min}(\Delta V(r))]_+\le B.
\]

If \(\ell(\ell+1)\ge B\), then the exact radial stiffness obeys

\[
L_\ell=-d^2/dr^2+\ell(\ell+1)/r^2+V_{\rm vac}+\Delta V
\succeq-d^2/dr^2+V_{\rm vac}.
\]

For real \(\nu\), use the \(e^{i\nu\tau}\) pencil \(P_\ell(\nu)=L_\ell-\nu^2I+i\nu G\), \(G=2\Omega J\). The charged vacuum \(2\times2\) block has minimum eigenvalue \(1-(\Omega+|\nu|)^2\), while the real-field block has \(2a-\nu^2\). Therefore

\[
\boxed{P_\ell(\nu)\succeq-d^2/dr^2+
\min\{1-(\Omega+|\nu|)^2,2a-\nu^2\}\,I\succ0}
\]

whenever \(|\nu|<\min(1-\Omega,\sqrt{2a})\). There is then no square-integrable eigenfunction at a real frequency in that common vacuum gap for such \(\ell\). This is a direct quadratic-form proof; no finite eigenvalue enumeration is needed once the global pointwise bound is established. The reporter is excluded from this gap independently by \(L_D\ge1\). The same stiffness comparison is positive at \(\nu=0\), so it also supplies the high-\(\ell\) energetic no-growth condition. It does not enumerate above-threshold scattering modes or resonances.

The integrating worker's `angular_bound_results.json` records a numerical maximum \(28.0100296087\) at \(r\simeq4.2459486063\) on the saved profile domain \([0,50]\). The maxima on 4,001, 8,001 and 16,001 sampled points approach that value. The declared margin comparison value is \(B=31\); **if** \(B_{\rm true}\le31\), the theorem excludes all \(\ell\ge6\), since \(6\cdot7=42\ge31\) while \(5\cdot6=30<31\). The extra \(\ell=5\) numerical calculation covers the remaining sector under that conditional cutoff.

This is not a certified half-line bound. The recorded flag is `half_line_bound_certified: false`. Sampled tail integrands are approximately \(0.00953571\) for \(r\ge20\) and \(0.000593414\) for \(r\ge25\) within the saved domain; neither those values nor the finite-domain endpoint establish an enclosure on \((50,\infty)\), interpolation error, or the exact continuum profile. The explicit margin is a numerical robustness choice, not a proof of the missing bound. Thus the high-\(\ell\) spectral exclusion remains conditional even after the finite low-\(\ell\) searches are completed.

### Finite-grid and quantum interpretation

The numerical `modes.py` convention \(q=(ru,rw,rz)\), \(G_{uw}=+2\Omega\), \(G_{wu}=-2\Omega\), its charge left-null vector, rank-one normalization and decoupled reporter block agree with this derivation. A consistent discrete background makes the discrete phase identity especially useful. A radial mesh/finite boundary does not preserve exact spatial translations; near-translation real or imaginary splittings require residual/overlap/refinement analysis. Identifying a symmetry artifact cannot justify discarding an unrelated unstable root.

Finite-box modes discretize continua; sparse shifts do not enumerate every eigenvalue. The exact centrifugal comparison extends positivity from an established lower sector for that operator, while extending it to the continuum requires the corresponding continuum hypotheses. Missing a subthreshold mode in bounded searches does not prove that none exists.

Quantizing this quadratic system requires symplectic mode normalization, treatment of charge/phase and center collective coordinates, continuum completeness and scattering data, and the correct fixed-charge reduction. A renormalized one-loop energy additionally needs the vacuum-subtracted inhomogeneous determinant, counterterms and the appropriate Legendre/fixed-charge matching. A raw sum of finite-box frequencies, or the constant-field Coleman–Weinberg potential evaluated on the profile, is not that calculation. The conditional no-growth argument concerns the leading classical linearized flow; it establishes no full interacting infinite-volume quantum spectrum, fission exclusion or quantum lifetime.

## Primary antecedents and derivation status

The equations, constrained energy identity, coupled ground-state factorization and rank-one inertia argument above are direct derivations for the stated Candidate 1 action. Established conserved-charge soliton antecedent: S. Coleman, *Q-balls*, https://doi.org/10.1016/0550-3213(85)90286-X . For the need to account for second-order corrections when assigning perturbation energy and charge, see M. N. Smolyakov, *Perturbations against a Q-ball: Charge, energy, and additivity property*, https://arxiv.org/abs/1711.05730 , and the continuation on nonoscillation modes, https://arxiv.org/abs/1906.02117 . No numerical bounds or quantum-stability conclusions from those different models are imported here.
