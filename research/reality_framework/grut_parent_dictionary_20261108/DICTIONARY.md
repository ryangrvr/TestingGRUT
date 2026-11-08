# The physical dictionary behind the finite parent

This derivation reviews the adopted Schwinger target and its implementation.
It is by the same Codex worker; it is not external or independent review. The
bosonized implementation was constructed independently of the lattice solver,
but P1 reuses that implementation. E1/E2 only compare reductions with its E12
truncation. Nothing below promotes that truncation to continuum truth.

## 1. Declared continuum theory and units

The target is the Schwinger model with one charged Dirac field on a fixed interval,
not 3+1 dimensional electrodynamics or an emergent spacetime. Adopt its bosonized
Hamiltonian at theta=0,

\[
H=\int_0^L dz\left[\tfrac12\Pi^2+\tfrac12(\partial_z\phi)^2
+\tfrac12\mu^2\phi^2
-\kappa:\cos(2\sqrt\pi\phi):_{\infty,\mu}\right],
\quad \mu=\frac{g}{\sqrt\pi},\quad
\kappa=\frac{m_R g e^\gamma}{2\pi^{3/2}}.
\]

Here phi is dimensionless in 1+1 dimensions, g and m_R have mass dimension one,
z and time have dimension minus one, and Pi has dimension one. The colons denote
infinite-line vacuum normal ordering at mass mu. This convention is part of the
mass dictionary, not an adjustable detector normalization. It follows the adopted
bosonization convention in [Okuda, equations (9)–(13)](https://arxiv.org/html/2210.00297v2).
The electric-field dictionary in this convention is E/g=phi/sqrt(pi) at theta=0.
At nonzero theta the explicit background term must be restored; this packet does
not test that case. The numerical builder specializes to g=1, theta=0 and equal
one-quarter Dirichlet boundaries. A generic-looking card does not establish
support for arbitrary g, theta or boundary values.

## 2. Boundary sector: zero incoming lattice flux does not mean zero smooth field

For the even neutral staggered chain, write Z_j=2 n_j-1 and

\[
q_j=\tfrac12[Z_j+(-1)^j],\qquad L_j-L_{j-1}=q_j,
\qquad L_{-1}=0.
\]

Direct summation gives the operator identity

\[
L_j=\tfrac12\sum_{i=0}^j Z_i+\tfrac14[1+(-1)^j].
\]

Thus the smooth link field contains a one-quarter offset as well as a staggered
c-number. Removing all of the mean would remove physical boundary information.
The continuum boundary map for this Gauss law is
phi(0)/sqrt(pi)=1/4 and phi(L)/sqrt(pi)=Q/2+1/4, where Q=sum Z_j. In our Q=0
sector both endpoints equal one quarter. This is the standard boundary map
adopted from [Okuda, equations (50), (62), (66)–(69)](https://arxiv.org/html/2210.00297v2),
including fractional boundary charges +1/2 and -1/2 in his convention. Neutral
bulk charge and a nonzero boundary-induced electric profile are compatible.

This is a continuum boundary-sector identification. It does not prove the exact
location of the continuum endpoints at finite lattice spacing. Okuda uses
L=(N-1)a; P1 inherits the M0-02/M1 ghost-endpoint convention L=(N+1)a, with sites
(j+1)a and link midpoints (j+3/2)a. The open free-hopping spectrum has lowest
particle-hole gap 2 sin[pi/(2(N+1))]/a, motivating that improved length. This
free-spectrum argument does not prove an exact interacting boundary matching.
The conventions approach the same continuum interval but differ by O(a) at
fixed physical L. No new endpoint displacement is fitted here.

## 3. Derive the massless mean; do not impose it in the massive theory

With fixed Dirichlet values, minimize the quadratic energy. Its first variation
requires (-partial_z^2+mu^2) bar(phi)=0, with both endpoints sqrt(pi)/4. The unique
solution is

\[
\bar\phi(z)=\frac{\sqrt\pi}{4}
\frac{\cosh[\mu(z-L/2)]}{\cosh(\mu L/2)}.
\]

Set phi=bar(phi)+psi with psi zero at the endpoints. Integration by parts makes
the cross term vanish: integral (bar(phi)' psi'+mu^2 bar(phi) psi)=0.
The remaining quadratic theory has modes

\[
u_n=\sqrt{2/L}\sin(n\pi z/L),\quad
\Omega_n=\sqrt{\mu^2+(n\pi/L)^2},\quad
\psi=\sum_{n\ge1}\frac{u_n}{\sqrt{2\Omega_n}}(a_n+a_n^\dagger).
\]

The massless Gaussian vacuum therefore has mean bar(phi). Once m_R is nonzero,
the cosine couples these modes. The same bar(phi) is only a coordinate shift;
the interacting ground state determines additional oscillator expectations.
Its mean is not supplied by the massless profile. In a stationary state, the
formal mean equation is

\[
(-\partial_z^2+\mu^2)\langle\phi\rangle
+2\sqrt\pi\kappa\langle:\sin(2\sqrt\pi\phi):_{\infty,\mu}\rangle=0,
\]

with the boundary values fixed. It is not a closed classical equation because
the sine expectation is not generally the sine of the expectation. The finite
parent diagonalizes the cosine; it does not make that mean-field replacement.
Cutoff versions can have residual violations of continuum identities.

## 4. Normal ordering and mass conventions

Let C_infty(r)=K_0(mu|r|)/(2pi) and C_D be the Dirichlet vacuum covariance.
Changing Wick conventions for beta=2 sqrt(pi) gives

\[
:e^{i\beta\phi}:_{\infty,\mu}
=e^{i\beta\bar\phi}\lambda(z):e^{i\beta\psi}:_D,\quad
\log\lambda=2\pi[C_\infty(0)-C_D(z,z)].
\]

The coincident difference is finite in the interior. Applying the image formula
for C_D yields the Bessel sum implemented in M1; its independent integral check
tests this subtraction. Near z=0, lambda~exp(-gamma)/(mu z), while
beta bar(phi) approaches pi/2. For each fixed finite occupation matrix element,
the vacuum cosine and nontrivial sine-mode terms remove that endpoint singularity.
This justifies ordinary interior quadrature for each such finite matrix. It is
not a uniform bound as the occupation or energy cutoff increases.

The lattice input is

\[
\frac{m_{\rm lat}}g=\frac{m_R}g-\frac{ag}{8},\qquad
H_{\rm spin}=\frac{g^2a}{2}\sum L_j^2
+\frac1{4a}\sum(X_jX_{j+1}+Y_jY_{j+1})
+\frac{m_{\rm lat}}2\sum(-1)^jZ_j+\text{constant}.
\]

The solver's dimensionless matrix is W=2H/(a g^2), x=(ag)^(-2), so H/g=W/(2sqrt(x)).
The mass shift is supported by the discrete chiral-symmetry argument and numerical
open-boundary improvement in [Dempsey et al., equation (3) and section IV.5](https://arxiv.org/html/2206.05308v3).
Open boundaries break the translation symmetry used in that argument. It does
not establish our interacting finite-interval mass calibration exactly; no
boundary counterterm or finite-cutoff mass fit is introduced. m_R=0.5g is a
Lagrangian mass parameter, not the measured bosonic gap (~1.61g at old cutoffs).

Absolute interacting vacuum energies here are unrenormalized cutoff diagnostics.
A state-independent constant cancels from normalized histories and energy gaps.
This observation does not justify removing state-dependent truncation errors,
nor prove ultraviolet convergence. No fitted vacuum or interaction counterterm
is used to obtain P1's verdict.

## 5. The detector is the same physical observable, not the same matrix

Fix

\[
f(z)=\frac{e^{-(z-L/2)^2/(2\sigma^2)}}
{\sigma\sqrt{2\pi}\,\mathrm{erf}[L/(2\sqrt2\sigma)]},\quad
F=\int_0^L f(z)\frac{\phi(z)}{\sqrt\pi}dz,\quad
W=e^{icF},\quad K_\pm=\frac{I\pm W}{2}.
\]

F and c are dimensionless; integral f=1. In P1 sigma=0.8/g, c=2, phase zero.
The physical instrument is complete because W is unitary. It can be defined by
a qubit-controlled Weyl interaction and a qubit readout. This extended,
instantaneous ideal instrument has no demonstrated local causal implementation
error budget in the present field-theory calculation.

The lattice observable is the normalized Gaussian-weighted raw electric link
field. Its exactly identifiable staggered c-number is (1/4) sum_j w_j(-1)^j;
the physical one-quarter boundary term must not be confused with it. P1 keeps
the existing raw-link choice rather than subtracting a new mean. Finite spacing,
edge sampling and boundary placement leave detector-dictionary errors; none is
bounded merely by matching Gaussian variance. The smearing tends toward the
declared continuum f, but finite equality has not been established.

Let h_n=integral f u_n/sqrt(2pi Omega_n), and bar(F)=integral f bar(phi)/sqrt(pi).
Since W is an infinite product of oscillator displacement operators, its matrix
between finite occupation states contains the factor
exp(i c bar(F)-c^2 sum_n h_n^2/2) multiplying normal-ordered oscillator matrices.
This accounts for unoccupied higher-mode vacuum fluctuations. The implementation
uses 8192 reference modes, numerical smearing integrals, and an analytic remaining
variance-tail estimate. Those quadrature and tail controls are distinct from
errors of the interacting state and Hamiltonian.

## 6. Two different truncations and their overflow

Let P_E retain all free-energy occupation tuples with sum Omega_n n_n<=E. Fixed
physical cutoffs are nested by padding higher-mode occupations with zero. The
implemented H_E=P_E H P_E and W_E=P_E W P_E use matrix elements of the same H,W.
In particular W_E is not exp(i c P_E F P_E). Define

\[
D_E=\frac{I_E-W_E^\dagger W_E}{2}\succeq0,
\quad K_{\pm,E}=\frac{I_E\pm W_E}{2}.
\]

K_+ and K_- together have total effect I_E-D_E. Complete them with an overflow
map of effect D_E into an absorbing orthogonal record. This is an honest
approximation instrument, not a new outcome of the intended continuum detector.
For an isometric inclusion I:E->F into a larger parent, put
R=W_F I-I W_E. Then I^dagger R=0 and

\[
D_E=I^\dagger D_F I+\tfrac12R^\dagger R.
\]

Proof: W_F I=I W_E+R and the two terms have orthogonal ranges. Expand its norm
square and subtract from I_E. The second term is detector transfer into states
newly represented between E and F; the first is still-unrepresented parent
overflow. This matrix identity does not attribute history probabilities
additively, because histories also depend on preparation, propagation and
interference. Parent truncation and the later E1 reduction are separate:
for E1's J:C^72->C^1280 replace I with J and W_E with J^dagger W_12 J to obtain
an analogous extra reduction defect. E2 restores that latter complement, not
states above E12 and not the interacting continuum.

For any complete physical distribution q with q(overflow)=0 and any normalized
truncated p with p(overflow)=o, the triangle inequality gives

\[
\mathrm{TV}(p,q)=\tfrac12\left(o+\sum_{s\ne o}|p_s-q_s|\right)\ge o,
\]

because the ordinary probabilities differ in total by -o. This necessary floor
needs no assumption about the dynamics. It is not an upper error bound.

## 7. What parent refinement can establish

The preparation rule remains the same (interacting ground state), while each
cutoff's approximation to that state changes. A fixed vector at all cutoffs would
instead change the physical preparation question. P1 measures nested-state
overlap as a diagnostic, without substituting a matched or refitted state.

Successive complete-history differences can show whether the present sequence
is stabilizing. They do not bound its infinite tail without additional control.
In particular, small overflow bounds only this instrument's instantaneous loss,
not omitted virtual dynamics or coherent returns. A massless Gaussian check
does not supply an interacting truncation theorem. Lattice-to-parent differences
cannot be identified uniquely with a boundary, mass, detector or continuum
failure while both approximations have unresolved error.

The two independent questions remain reduced-to-parent accuracy and
parent-to-continuum accuracy. P1 finishes an explicit dictionary and tests one
further regulated comparison. It neither supplies a new GRUT restriction nor
derives the adopted quantum field theory from the operational framework.
