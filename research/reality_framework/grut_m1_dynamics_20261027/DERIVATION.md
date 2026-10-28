# Equations and what the experiment measures

Use the exact fixed GRUT-M1 light parameters in the preceding candidate packet.
Put m_X=1 and canonical q=(sqrt(2) Re X,sqrt(2) Im X,phi), so v=10,
g_X²=.01, lambda_X=.001, lambda_phi=.03 and epsilon=.002. The local experiment is

\[
 H=\int d^3x\left\{\tfrac12\sum_i(p_i^2+|\nabla q_i|^2)
 +\frac{\lambda_\phi}{4}(\phi^2-v^2)^2-\epsilon v(\phi-v)^3
 +\tfrac12g_X^2\phi^2(q_1^2+q_2^2)
 +\tfrac14\lambda_X(q_1^2+q_2^2)^2\right\}.
\]

The phi force is

\[
 \ddot\phi=\nabla^2\phi-
 [\lambda_\phi(\phi^2-v^2)+g_X^2(q_1^2+q_2^2)]\phi
 +3\epsilon v(\phi-v)^2.
\]

The charged equations have effective mass squared
g_X² phi²+lambda_X(q1²+q2²). Their exact conserved scalar current is
j0=q1 p2-q2 p1. The tiny charge-six term and gauge dynamics are omitted here,
under the previously stated cold/local approximation. This experiment does not
calculate the actual early gauge plasma or the full hot-core lifetime.

For a circular uniform field q1+i q2=v F exp(i Omega t), phi=v S, put
a=lambda_X/g_X²=.1, b=lambda_phi/g_X²=3 and c=epsilon/g_X²=.2. Then

\[
 b(S^2-1)S-3c(S-1)^2+F^2S=0,\qquad
 \Omega^2=S^2+aF^2.
\]

For every F>=0 this has just one equilibrium phi value, 0<S<=1. On that interval
the equation is equivalent to

\[
 F^2=b(1-S^2)+3c(S-2+S^{-1}),
\]

whose derivative is -2bS+3c(1-S^-2)<0. For S<0 both the vacuum-potential derivative
and the added F²S are negative; for S>1 they are positive. Consequently adding a
constant charged density does not resurrect two sign vacua. This does not prove
that a driven field never crosses zero or that all transient interfaces disappear.

Define U_phi(S)=b(S²-1)²/4-c(S-1)³. The infinite, uniform zero-pressure core obeys

\[
 \Omega_0^2=\min_{S\in[0,1]}
 \left[S^2+2\sqrt{aU_\phi(S)}\right],\qquad
 F_0^2=2\sqrt{U_\phi(S_0)/a}.
\]

This follows by minimizing 2U/F² first over F². Negative S has larger U_phi than
the corresponding positive S, while S>1 already gives Omega²>=1. Completing the
charged kinetic squares gives the classical global bound E>=Omega0*|Q| in m_X
units. Positive gradients only strengthen it. The minimum is a bulk limit, not
the mass or selected charge of a finite object. Surface energy raises E(Q)/Q;
the M1 reference object at Q=200000 has E/(Q m_X)=.93739999 rather than the bulk
value near .776. Coarsening therefore makes the assumed final size distribution
an essential input to the mass normalization.

For these exact a,b,c, the stationary condition reduces to
2.1 S^4+3.76 S^3+2.864 S²-.288 S-.036=0. It has exactly one positive root by its
single coefficient-sign change, and the unsquared equation has positive sides
on (0,1). Its root therefore selects the unique interior minimum; both endpoints
are also checked. The reported bulk value does not rely on assuming a numerical
minimizer found the global rather than a local minimum.

The fourth-order spatial operator is the five-point centered derivative in each
direction. The engine's signed-link energy has the equivalent positive local
gradient assignment

\[
 E_{\nabla,i}=\frac1{2h^2}\sum_a
 \left\{\tfrac12[(q_{i+e_a}-q_i)^2+(q_i-q_{i-e_a})^2]
 +\frac{(q_{i+e_a}-2q_i+q_{i-e_a})^2}{12}\right\}.
\]

The two assignments sum to the same torus Hamiltonian; regional energies differ
by boundary terms. A low regional E/Q embedded in a dense background is not by
itself an isolated-object certification. We report shell field values, radius
dependence and persistence. There is no subtraction of an arbitrarily chosen
ambient energy to force a bound verdict.

The integrator is a fourth-order symmetric composition of conservative kicks
and drifts. The negative middle substep is part of that standard composition.
No friction or cooling is inserted. Periodic radiation can return to a core;
the finite-time survival test is not an asymptotic scattering calculation.

The elliptical control is a canonical stretch q1,p1 by sqrt(r) and q2,p2 by
1/sqrt(r), followed by the same exact total-Q normalization as the circular
control. r=(9+sqrt(77))/2 gives (r+1/r-2)/4=1.75 in a harmonic reference. In the
interacting model it changes the energy and excites the binding field. It is
not the actual quantum release spectrum, its chemical potential, or a physical
antiparticle preparation.

For exterior object waves only, taper the charged field and momentum by a
time-independent W(r). In the vacuum reference with the actual lattice symbol,
omega_k=sqrt(1+k_FD4²),

\[
 N_\pm=\tfrac12\sum_k\left[
 \frac{|p_{1k}|^2+|p_{2k}|^2+\omega_k^2(|q_{1k}|^2+|q_{2k}|^2)}{2\omega_k}
 \ \pm\operatorname{Re}(q_{1k}^*p_{2k}-q_{2k}^*p_{1k})\right]V.
\]

The rFFT multiplicities are included once, and Nplus-Nminus agrees with the
real-space tapered charge. The interaction residual and taper-gradient energy
change are reported. A weak exterior with a mass gap supports a free-wave
occupation diagnostic. It does not turn a coherent dense orbit into particles
or yield the cosmological p. Thresholds and several cut/taper choices stay visible.

The frozen M1 mass criterion remains
.0667805(1-eta)+2.13356 p<=.01 for its stated cold, positive-Q200000-object plus
free-X/pair completion. The simulated mask fractions are not eta. A different
object spectrum needs its actual charge-weighted E(Q)/Q, and the unredshifted
heat in this conservative box is not silently discarded. Actual expansion,
timed sources and quantum number-changing reactions are needed for the final
mass and antiparticle comparison.
