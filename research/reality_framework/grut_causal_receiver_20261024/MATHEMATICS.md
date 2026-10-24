# A causal massive-wave receiver

All fields and times use the same dimensionless action and units as the
parent prepared radial runs. This is numerical boundary construction, not
an additional physical interaction or a selected cosmological source.

## 1. Direct causal realization

In a homogeneous exterior let U=rX and Z=r(S-1). Each free real component
has mass m=1 for X or m=sqrt(6) for Z. For Laplace variable s choose

\[
\Gamma(s,x)=1+\frac{\sigma(x)}{s+\alpha},\qquad
C=\Gamma^{-1}=\frac{s+\alpha}{s+\alpha+\sigma},
\quad\alpha=.1,\quad\sigma=20[(r-R)/D]^3.
\]

The strength is zero in the physical region and clipped at the layer end.
The vacuum matched equation is U_tt=C_node div(C_face grad U)-m^2 U.
Introduce face memory psi and node memory q. Its direct time-domain form is

\[
\begin{aligned}
\dot U&=V,\\
\dot V&=L U+\operatorname{div}\psi+q-F(U),\\
\dot\psi&=-(\alpha+\sigma_f)\psi-\sigma_f\operatorname{grad}U,\\
\dot q&=-(\alpha+\sigma_n)q-
\sigma_n(LU+\operatorname{div}\psi).
\end{aligned}
\]

Here L=div grad and F=m^2 U in a free component. In the nonlinear runs F
retains the actual quartic force everywhere, including the layer. The
homogeneous free matching theorem consequently is not a nonlinear matching
theorem. Radius, strength and finite-time controls remain necessary.

The initial memory is psi=q=0, specifying the direct equation's causal zero
prehistory. All initial physical Gaussian tails are retained, and their
charge fraction inside the layer is measured. At the formation layer start
R=360 that fraction is about 1.63e-15. Zero initial memory does not silently
replace nonzero physical fields by zero. An exploratory first-order
realization with extra memory states was replaced before the released runs;
those pilot traces are excluded from all production counts and conclusions.

## 2. Spatial Hamiltonian and exact lattice leads

U and V live at nodes with zero Dirichlet endpoints and odd continuation.
Face gradients use

\[
(GU)_{i+1/2}=\frac{9}{8h}(U_{i+1}-U_i)
-\frac{1}{24h}(U_{i+2}-U_{i-1}).
\]

The divergence is the negative transpose, giving a self-adjoint L=-G^T G.
Its seven-point coefficients are -365/144 at the center, 87/64 at each
first neighbor, -3/32 at each second neighbor, and 1/576 at each third
neighbor, all divided by h^2. Face ghosts are even at the reflecting wall.
This differs from the parent's five-point nodal Hamiltonian; the action and
smooth preparation are unchanged and the resulting core is directly refined.

For lattice phase k h, the physical-region frequency is

\[
\omega_h^2=m^2+\left[\frac{2}{h}\sin(kh/2)
\left(1+\frac16\sin^2(kh/2)\right)\right]^2.
\]

The boundary-value test uses this frequency, the actual node/face C factors,
the outer wall and all six spatial roots. The incident/reflected oscillatory
roots are exp(plus/minus i k h); the two left-decaying evanescent roots are
also matched, without an approximate continuum Robin condition. The
zero-layer closed-wall amplitude is independently -exp(2 i k h N).
The ratio |R|^2 is an energy reflection ratio because the incident and
reflected leads share their frequency and opposite group velocities.

## 3. Causal stability and time integration

For constant coefficients, a=alpha+sigma, the dispersion polynomial is

\[
(s^2+m^2)(s+a)^2+k^2(s+\alpha)^2=0.
\]

Its coefficients are 1, 2a, a^2+m^2+k^2,
2(am^2+alpha k^2), and a^2m^2+alpha^2k^2. For sigma>0 and k>0, the two
nontrivial Routh-Hurwitz determinants are

\[
\Delta_2=2(a^3+\sigma k^2)>0,\qquad
\Delta_3=4\sigma k^2(am^2+\alpha a^2+\alpha k^2)>0.
\]

This establishes constant-coefficient continuum decay. It is not a proof
for every varying profile, nonlinear state or grid. The actual C++ RHS is
therefore exported as a finite-grid linear matrix for two representative
resolutions and both masses; its complete spectra and RK4 amplification
factors are recorded. Propagating free packets and the nonlinear production
runs provide additional time-domain controls.

The executable uses classical fourth-order Runge-Kutta. Its small bulk
energy/charge drift is measured, not booked as physical emission. A
large-domain sine-mode reference is evolved both with exact time phases
and the same RK4 bulk multiplier. The latter comparison isolates the layer
from bulk time error; the former tests their combination. The long soft
packet uses a separate coarser free grid which is also reflection-tested.
The free_s test uses complex bookkeeping for two independent real free
components of mass sqrt(6). This tests the neutral-mass operator; it does
not introduce a charged neutral-sector particle into the physical action.

## 4. A signed field-to-boundary ledger

Use physical Hamiltonian H with positive discrete gradient energy
4 pi h sum |GU|^2/2 and the same quartic potential. Physical canonical
charge is Q=4 pi h sum(y V_x-x V_y). Let

\[
B_f=\operatorname{div}\psi_f+q_f.
\]

Only that term changes the physical field's boundary ledger:

\[
\dot E_{rec}=-4\pi h\sum_f\sum_i V_{f,i}B_{f,i},\qquad
\dot Q_{rec}=4\pi h\sum_i(x_iB_{y,i}-y_iB_{x,i}).
\]

The same RK stages integrate these rates. Thus H+E_rec and Q+Q_rec
have measured integration errors. The receiver work is signed and can
include reversible auxiliary storage; it is not a positive spectral count
or a quantum particle source. Independent finite variations of the actual
executable's H and Q check these identities, and a separate direct-gradient
sum checks the seven-point energy against its face-gradient form.

An instrumented copy adds an exact discrete extraction cut between nodes
c and c+1 inside the physical region. It changes no field equation. With
L_ij the symmetric hopping matrix, the outward signed currents are

\[
J_Q=4\pi h\sum_{i\le c<j}L_{ij}
\operatorname{Im}(U_i^*U_j),\qquad
J_E=2\pi h\sum_f\sum_{i\le c<j}L_{ij}
(U_{f,i}V_{f,j}-V_{f,i}U_{f,j}).
\]

Only the three hopping distances enter. E_inside uses the corresponding
node energy assignment 4 pi h [V^2/2-U L U/2+potential]. Then
E_inside+integral J_E and Q_inside+integral J_Q close with measured RK4
errors, without a frequency cut or free-particle assumption. The assignment
is a discrete Hamiltonian subregion, not a claim that every local gradient
term is positive. For a probe node at 240 and h=.25 the cut is at 240.125.
Field, trace and memory hashes are directly compared with the otherwise
identical uninstrumented run, so instrumentation cannot silently change
the measured dynamics. Signed currents retain any real returning wave.

## 5. What is and is not obtained

The reflection benchmark tests the .01-to-1 band and extends samples to
k=3 for each physical mass. These are sampled linear vacuum controls, not
an unsampled-frequency or nonlinear theorem. Nonlinear mass shifts in the
layer are measured. Their size relative to k^2 still matters for treating
the field as asymptotic vacuum particles, even when the numerical layer is
accurate in the vacuum test.

The prepared spherical runs can measure finite-time core charge, a boundary
ledger and extracted classical fields. Formation time in these units does
not select a cosmic epoch; spherical symmetry does not produce a 3D charge
distribution. Quantum scalar decays, general hot-core leakage, a complete
bound-core subtraction, the cosmological source and vacuum density are not
derived by adding this numerical receiver.
