# Linear scalar/longitudinal equations

The homogeneous equations and potential are those pinned from e160356.
Here the metric is unperturbed FLRW; its scale factor and Ricci curvature
are evolved with the background. Define canonical real background
`r=sqrt2(Re X,Im X,Re Y,Im Y,Re Sigma,Im Sigma)`, rescaled field
`chi=a r`, and conformal time `du=d tau/a`, tau=H_i t. All field/momentum
units follow the dimensionless H_i normalization in the parent.

Let `T=diag(J,3J,6J)`, `J=[[0,-1],[1,0]]`, `g=T chi`, `g'=T chi'`.
For fixed nonzero comoving k, take Coulomb gauge and define perturbations
q=delta chi. The full canonical potential matrix is

\[
M=a^2[V_{,rr}+\xi_X R P_X]-\frac{a''}{a}1,
\qquad\frac{a''}{a}=\frac{a^2R}{6}.
\]

R in this equation is in H_i units. P_X projects on the two real X
components. V_{,rr} is the entire six-real-field Hessian, including the
kappa and mu off-diagonal couplings. Holding Y fixed and checking only
the X block would not give these equations.

Before eliminating A0 the quadratic scalar Lagrangian contains
`e A0(g.q'-g'.q)+(k^2+e^2g^2) A0^2/2`. Its elimination gives

\[
L_2=\tfrac12q'^2-\tfrac12q^T(k^21+M)q
-\frac{e^2}{2(k^2+e^2g^2)}(g\cdot q'-g'\cdot q)^2.
\]

Legendre transforming, with p the canonical momentum, yields

\[
H_2=\tfrac12p^2+\tfrac12q^T(k^21+M)q
+\frac{e^2}{2k^2}(g\cdot p-g'\cdot q)^2.
\]

Thus, with c=e^2/k^2 and j=g.p-g'.q,

\[
q'=p+cgj,\qquad p'=-(k^21+M)q+cg'j.
\]

The six scalar coordinates include the longitudinal gauge degree of
freedom; the two transverse vector modes are separate and not calculated.
The generator divided by a evolves tau derivatives. The C engine uses
the gauge Ward identity `M g=-g''` to evaluate the common component
without cancellation of large Higgs Hessian terms. The Python Radau
control instead directly constructs the block generator, providing a
distinct numerical implementation. No quantum-gravity constraint
solution is presumed by eliminating this ordinary gauge constraint.

The canonical12-by12 transfer matrix satisfies `Tcan(0)=1` and the
time-dependent Hamiltonian equation above. Its symplectic condition is
`Tcan^T J12 Tcan=J12`. We monitor absolute and scaled residuals. The
k=0,e=0 global symmetry tangent is an algebraic control, not a Coulomb
gauge definition at k=0.

## Gauge-invariant response rows

The fractional amplitude row is
`delta|X|/|X|=chi_X.delta q_X/|chi_X|^2`. The local comoving charge is
`N_X=-(J chi_X).chi'_X`, with perturbation

\[
\delta N_X=-(Jq_X)\cdot\chi'_X-(J\chi_X)\cdot p_X.
\]

The canonical momentum already contains the eliminated A0 term. A pure
common gauge tangent `(q,p)=(g,g')` gives zero amplitude and charge
response, which also checks these observables. Raw common charged-field
coordinate growth would not be sufficient evidence of physical breakup.

For the generic input norm, each species' two real fields are scaled by
its initial norm D_i, and velocities by H_i D_i. At a=H/H_i=1 and zero
background velocities, the map from such physical inputs to(q,p) is

\[
S_0=\begin{pmatrix}D&0\\D&A^{-1}D\end{pmatrix},
\quad A^{-1}=1-\frac{e^2gg^T}{k^2+e^2g^2}.
\]

The reported response coefficient is the Euclidean norm of the observable
row times Tcan S_0. This norm fixes units; it does not select a state.

The supplied phase-cosine test differentiates the same conditional
preparation at Sigma phase `-pi/2+delta`, keeping X positive real and
readjusting Y/radii. Its initial q derivative is obtained at
delta=+/-1e-5. Initial p=q because velocities vanish and g'=g initially.
The zero-k tangent is checked against separately evolved +/-1e-4 full
homogeneous histories. Cosine amplitude .03 is a chosen seed, not an
inflationary fluctuation prediction or an integrated stochastic spectrum.

For the auxiliary six-field Gaussian control, diagonalize M(0). Set
`W=sqrt(k^2+max(M(0),0)+1)` by spectral calculus and
`C_0=diag(W^-1,W)/2`. This is a positive canonical covariance with hbar1.
It is not derived from a phase holder or a de Sitter vacuum. The displayed
variances use k^3/(2pi^2) per log k, include zero-point fluctuations and
are not a renormalized abundance or interacting backreaction estimate.

## Isolated quartic-orbit diagnostic

For the matched two-real X orbit,

\[
q''=-\lambda_X(q^2)q,\qquad
\delta q''+[k^21+\lambda_X(q^21+2qq^T)]\delta q=0.
\]

Its radial motion is periodic while its angle can advance by Delta theta
per radial period P. If Phi(P) is the perturbation transfer over that
fixed period, close by the same O(2) symmetry:
`M_rel=diag(R(-Delta theta),R(-Delta theta)) Phi(P)`.
Its eigenvalue moduli give the relative-period growth rates
`mu=log|lambda|/P`. This is a true periodic-orbit diagnostic for this
auxiliary central potential, unlike a constant Floquet exponent assigned
to the expanding six-field background. The matched period is.81404 in
H_i eta and the angular advance2.99247rad. The stable circular control
and symplectic residual are reported with the computed band.

## Binding-field perturbations

Since phi background0, it decouples linearly from the other perturbations:

\[
v_\phi''+[k^2+a^2(2g_X^2|X|^2-\lambda_\phi v^2)-a''/a]v_\phi=0.
\]

xi_phi=0 is the reference matching choice, not an all-scale symmetry.
We compare initial canonical instantaneous and first-order WKB Gaussian
seeds. For initial positive omega_0, the second has
`u_0=1/sqrt(2omega_0)`, `u'_0=(-iomega_0-d_0)u_0`,
`d_0=omega'_0/(2omega_0)=M_phi,conf(0)/(2omega_0^2)` for the supplied
initial zero-velocity background. Its covariance includes the real
q-p correlation. At later positive omega, we report the instantaneous
projection `n=(C_pp+omega^2 C_qq)/(2omega)-1/2` for both seeds and
integrate `k^3 omega n/(2pi^2)` over the stated band. The projection
convention and zero-point subtraction are explicit; this is not a
renormalized interacting stress tensor. A ratio reaching unity would
invalidate the prescribed-background approximation rather than create
physical energy from nothing. No such crossing is seen here through80.
