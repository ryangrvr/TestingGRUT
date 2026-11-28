# C2 continuum certificate method

The object is the leading flat scalar saddle used in C1. The action and matching are those in the [C1 report](../grut_c1_quantum_channels_20261127/REPORT.md). A certificate of this saddle and its linear fluctuations supplies a premise of C1. Its one-loop numerical allocations and its omitted corrections remain separate.

## Operator and normalization

Set $X=vF(\rho)e^{-i\Omega\tau}/\sqrt2$ and $\phi=vS(\rho)$, with $\rho=g_Xvr$ and $\tau=g_Xvt$. Put $y=(u,w)=(\rho F,\rho(S-1))$, $I=\int u^2\,d\rho$ and $\bar Q=\Omega I=g_X^2 Q/(4\pi)$. The stationary system on the half-line is
\[
G(y,\Omega)=
\binom{-u''+\rho(S^2+F^2/10-\Omega^2)F}
{-w''+\rho(S^2-1+F^2)S}=0.
\]
The derivative in $y$ is $L_+=-\partial_\rho^2+V_+$, with entries $S^2+0.3F^2-\Omega^2$, $3S^2-1+F^2$ on the diagonal and $2FS$ off diagonal. The centrifugal term is $\ell(\ell+1)/\rho^2$ for all components. The phase block has potential $S^2+0.1F^2-\Omega^2$. In amplitude, phase and real-field coordinates the velocity term couples the first two with antisymmetric coefficients $\pm 2\Omega$. The reporter form has potential $1+0.01S^2$ and is bounded below by the identity.

In the radial sector let $e=(u,0)/\sqrt I$. The momentum conjugate to the normalized phase coordinate is its velocity minus $2\Omega\langle e,\eta\rangle$, up to the choice of rotating-phase sign. Fixing this momentum and quotienting phase gives $L_Q=L_++4\Omega^2|e\rangle\langle e|$. Differentiation on a fixed-frequency branch gives $v=\partial_\Omega y=L_+^{-1}(2\Omega u,0)$ and
\[
\bar Q'=I+2\Omega\langle u,v_1\rangle,\qquad
1+4\Omega^2\langle e,L_+^{-1}e\rangle=\bar Q'/I.
\]
The block form with diagonal $(L_+,-1)$ and off-diagonal $2\Omega e$ has two Schur reductions. Its first removes a negative scalar and leaves $L_Q$. Its second leaves $L_+$ and the scalar $-\bar Q'/I$. If $L_+$ has index one, is invertible with positive essential spectral floor, and $\bar Q'<0$, these reductions give zero negative index and zero kernel for $L_Q$. The positive essential floor persists under finite rank perturbation. This is the constrained quadratic energy of the classical background. [Conventional two-field comparison](https://arxiv.org/abs/1612.00737)

## From local polynomial bounds to a half-line mode count

The predictor uses the fine C1 knot array, explicit Hermite interpolation data and a short final polynomial joining the vacuum through second derivative. Every interval operation encloses the polynomial defined by those binary nodal data. The stationary and tangent defects are polynomial in the normalized cell coordinate. Power-to-Bernstein conversion bounds their suprema and provides the interval integral of the charge and tangent pairing. At the origin $F',S'$ vanish. Outside the tapered support the predictor is the vacuum.

For an eigenvalue test at $\lambda$, split the exact predictor operator at the taper radius with a Neumann form on each side. The exterior is a positive vacuum form for the frozen shifts. On a cell of length $h$, choose a constant matrix $B$ below $V-\lambda I$, including a lower centrifugal value at the outer edge. Component intervals and maximum row-sum errors supply this lower bound. For the nodal affine interpolant $p$ and the endpoint-zero remainder $z$, derivative energies split exactly. With $M\ge\|B\|$ the remainder is bounded below by $(9/h^2-M)\|z\|^2$. Its cross term with $p$ is bounded by completing the square. Hence the lower nodal form uses
\[
B-\frac{M^2}{9/h^2-M}I.
\]
The Neumann endpoint remains a free nodal variable. The ghost exterior creates no negative index. A nonnegative subspace of the lower nodal form pulls back to a subspace of the continuum form of at most the same codimension by nodal evaluation. Thus its interval LDL negative index bounds the number of continuum eigenvalues below the shift from above. An interval pivot whose sign is unresolved leaves the count unresolved. A continuous compactly supported Dirichlet trial with a verified negative shifted form is a separate upper bound on the first eigenvalue.

## Validation of the exact stationary saddle

On $(H^2\cap H_0^1)(0,\infty)^2\times\mathbb R$ use the maximum of $\|\eta\|_2$, $\|\eta''\|_2/16$ and $|\delta\Omega|$. The value space uses the maximum of the $L^2$ defect and scalar charge defect. For a perturbation of norm at most $r$, integration by parts gives $\|\eta'\|_2\le4r$. The one-dimensional Sobolev estimate gives $\|\eta/\rho\|_\infty\le\|\eta'\|_\infty\le\sqrt{128}r<12r$ and $\|\eta\|_\infty\le3r$.

Let $L_0$ and $\Omega_0$ be those of the predictor. A negative trial below $-g$ and an index at most one below $g$ give $\|L_0^{-1}\|_{2\to2}\le1/g$. A separately solved tangent predictor $\tilde v$ is enclosed by its residual, $\|v_0-\tilde v\|_2\le\|L_0\tilde v-(2\Omega_0u_0,0)\|_2/g$. Invert the bordered stationary and charge derivative by the Schur denominator $S_0=I_0+2\Omega_0\langle u_0,v_{0,1}\rangle$. The solution of a unit-bounded right hand side $(f,b)$ has
\[
|\delta\Omega|\le d=(1+\|v_0\|_2)/|S_0|,
\quad
\|\eta\|_2\le K=1/g+\|v_0\|_2d,\quad
\|\eta''\|_2\le \|V_0\|K+2|\Omega_0|\|u_0\|_2d+1.
\]
Write $\gamma$ for the resulting norm bound on the bordered inverse. The stationary map uses the normalized charge $\Omega I-\bar Q$. On the prescribed tube, the maximum row-sum variation of $V_+$ is at most $336r$. The same bound also covers the phase potential. Coupling and charge rows yield a full derivative bound $360r$ using the displayed tube and norm caps. With full predictor residual $\epsilon$, the frozen radius is $r=2\gamma\epsilon$. Verify $360\gamma r<1/2$ and the tube caps. The exact frozen derivative inverse then makes the Newton map a contraction of the closed ball. This certifies the charge, regularity, stationary profile, and frequency on the half-line. The charge-normalized branch derivative is obtained by inverse perturbation and the tangent defect on this exact ball.

The weak solution has $u,w=O(\rho)$ and the stationary equations then give $u'',w''=O(\rho)$ at the origin. At the exterior, $y''=Ay$, with diagonal $S^2+0.1F^2-\Omega^2$ and $2+3w/\rho+(w/\rho)^2$, lower off-diagonal $FS$ and zero upper off-diagonal. If the symmetric part exceeds $\kappa^2 I$ and $\|A\|\le4$, the distributional inequality $|y|''\ge\kappa^2|y|$ gives $|y(\rho)|\le3r e^{-\kappa(\rho-R)}$. The end conditions follow from $H^2$, and integration of $y''$ gives $|y'|\le (12r/\kappa)e^{-\kappa(\rho-R)}$. The spectral comparison is global using the certified profile tube.

## Symmetries and gyroscopic evolution

An upper count of one below a positive shift in the phase $\ell=0$ form and the amplitude $\ell=1$ form, with perturbation smaller than the shift, proves a positive complement once the exact nonzero modes are identified. They are $u$ and $\rho(F',S')$, respectively. The second vector corresponds to translations with three angular components. Higher sectors add a strictly positive centrifugal difference to an already nonnegative operator. Their essential floor is positive, and angular ordering reduces the gap question to the first higher sector.

The energy of the linearized Hamiltonian equations is conserved because the gyro term is skew. On the fixed-charge, zero-rest-momentum phase space with the symmetry coordinates quotiented, the positive stiffness and kinetic forms control the energy norm. In the translation sector, imposing zero momentum adds a nonnegative square by the same cyclic-coordinate reduction. This proves classical continuum linear stability in the specified scalar dynamics with the identified collective modes removed. Nonlinear orbital dynamics, quantum lifetimes, gravitational perturbations and higher orders are separate.
