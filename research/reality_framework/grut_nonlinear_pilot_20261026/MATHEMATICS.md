# Discrete equations and interpretation

The canonical fields obey qddot_a=Lap q_a-dU/dq_a, with

\[
U=\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
+\frac{g^2}{2}\phi^2(q_1^2+q_2^2)
+\frac{\lambda_X}{4}(q_1^2+q_2^2)^2.
\]

The local current is j0=q1 p2-q2 p1. Its volume integral is conserved by the global U(1) symmetry of this truncated classical action.

For periodic spacing d, the spatial operator in each direction is

\[
L_4 f_i=\frac{-f_{i+2}+16f_{i+1}-30f_i+16f_{i-1}-f_{i-2}}{12d^2}.
\]

Its positive symbol is k_eff²=4 sin²(kd/2)[1+sin²(kd/2)/3]/d². Summing over three directions gives the implemented dispersion. Its matching gradient energy, summed over all fields and directions, is

\[
H_{\rm grad}=\frac{d^3}{2d^2}\sum_i
\left[\frac43(f_{i+1}-f_i)^2-\frac1{12}(f_{i+2}-f_i)^2\right].
\]

This equals -d³ sum f L4 f/2 globally. The difference form avoids cancellation of a large constant background in the energy diagnostic. Local allocation of this discrete edge energy across a candidate boundary is not a unique physical subsystem mass.

Kick and drift maps are p+=h(L4 q-grad U) and q+=h p. Each conserves the global U(1) charge up to floating-point arithmetic. The fourth-order composition uses w1=1/(2-cuberoot2), w0=1-2w1 and the symmetric three-step Verlet composition. Adjacent kicks are combined. State evolution is double precision; registered snapshots are double precision for production.

## Matched preparation

For chosen radius r(x) and phase theta(x), set q=(r cos theta,r sin theta) and tangential momentum Q0/r. This fixes j0(x)=Q0 at the initial time. A uniform inward radial momentum is selected from the full gradient, potential and tangential kinetic energies to give the supplied mean total energy. The scalar momentum initially vanishes. Fourier coefficients, phases and the conserved targets are retained in RUNS.json.

## Candidate geometry

The diagnostic smooths j0 and phi² periodically with a Gaussian of width1/m_X. It selects smoothed j0 above2,4 or8 times its mean and smoothed phi² below .25v². Six-neighbour connected components are joined across periodic faces. Charge-weighted circular coordinates define centers; a minimal-image radius and coherence flag identify noncompact components. None of this imposes Q or establishes E(Q).

Persistence means an overlapping pattern, not a conserved particle identity. Successive masks need Dice overlap>0.2 and charge ratio between.5 and2. At least four stored snapshots spanning about60/m_X are required. Snapshot cadence does not resolve every short collision or threshold excursion. Expanded masks use periodic six-neighbour dilation and can overlap other structures, so their charges are not added into a population ledger.

## Spectra and reference occupations

Real-field Fourier transforms use forward normalization. The half-space weights are2 except at the real-axis zero/Nyquist planes. Parseval identities check energy and charge normalization. The positive vacuum-reference frequency is omega_k=sqrt(1+k_eff²). For the charged pair,

\[
n_{\rm ref}(k)=\frac{|p_1(k)|^2+|p_2(k)|^2+
\omega_k^2(|q_1(k)|^2+|q_2(k)|^2)}{2\omega_k},
\quad Q(k)=\operatorname{Re}[q_1^*p_2-q_2^*p_1].
\]

Reference n_plus/minus are half the weighted sum of n_ref±Q. Their difference reproduces the actual conserved charge even in strongly interacting fields. Their sum is not thereby a physical particle count. The diagnostic interaction norm is

\[
I_X=\frac{\left\langle\tfrac12|g^2\phi^2-1|(q_1^2+q_2^2)
+\tfrac14\lambda_X(q_1^2+q_2^2)^2\right\rangle}
{\rho_{X,\rm vacuum\ reference}}.
\]

It is reported without using cancellations to hide coupling energy. A large global value leaves the vacuum-reference particle interpretation uncertified; it does not say that every exterior region lacks weak waves. No scalar particle projection is assigned across both signs of the binding vacuum, and no local negative current is identified with an antiparticle density.

Boundary-shell spectra use either continuum radial |k| or max component |k_i|, as labelled, and include kinetic-plus-gradient energy only. A classical lattice cutoff, statistical seed and periodic box are supplied approximations. Quantum collision rates and a regulator-independent plasma state require a separate match.
