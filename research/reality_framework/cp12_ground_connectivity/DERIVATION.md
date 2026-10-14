# CP-12 ground connectivity: full-sector reduction and electric-tail control

14 October 2026 UTC. Additive follow-up to connected-bridge checkpoint
3b914e5a730a8d246175159cf5bd4d3f6870404a. All formulas below use the frozen
CP-12 Hamiltonian and the connection sector Q=N. No thermal state is assumed.

## 1. Declared physical problem

Let M=N(N-1)/2. A basis state specifies an N-vertex simple graph G with Q
occupied connections, fermion occupations, and integer electric fields on its
two incidence legs per present edge. Absent edges have their center fermion
filled and both fields zero. The physical subspace satisfies every vertex
and edge-center Gauss constraint and total fermion number M.

At the frozen parameters, in dimensionless model units,

\[
H=\tfrac12 Q+2N_v+\tfrac18\sum_{e,i\in e}E_{ie}^2+T+R,
\qquad t_{\rm hop}=1,\quad\Gamma=1/10.
\]

T is the signed gauge-covariant fermion hopping operator. R is the positive
neutral-rewiring Laplacian with off-diagonal amplitude -Gamma/M, including its
full diagonal degree. Its positivity holds in each fixed fermion/electric
fiber, also for unbounded electric fields.

We compute the zero-temperature ground state of H, not the ground state of a
degree-only energy or the disconnected matching trial from earlier work.
The rotor regulator P_L retains all physical states with |E_ie|<=L and uses
the compression A_L=P_L H P_L. The regulator removes hopping matrix elements
through its boundary; it does not refit the Hamiltonian.

## 2. All permutation sectors, without imposing a singlet ground

Simultaneous event permutations, oriented-edge signs, and fermionic Fock signs
give a signed unitary representation commuting with H. We use the same
neutral-background character convention as the preceding invariant solver.
Multiplying by this one-dimensional character relabels irreducible
representations but does not change energies or invariant graph observables.
Here this symmetry is a computational decomposition; no decision that event
permutations are a gauge redundancy is needed.

For an orbit representative s, let its orbit size be o_s and its signed
stabilizer be K_s. For a real orthogonal irreducible representation D_lambda,

\[
P_s^\lambda=\frac1{|K_s|}\sum_{h\in K_s}\sigma_s(h)D_\lambda(h)
\]

is an orthogonal projector. Retain an orthonormal basis V_s of its image,
with rank r_s. Its null directions must be removed; retaining them would
introduce artificial zero eigenvalues and misidentify sectors with positive
lowest energy.

Suppose a physical transition of coefficient h_ts takes the representative s
to a basis state whose canonical representative is t, and g maps that state
to t with Fock sign sigma. The reduced contribution is

\[
\sqrt{\frac{o_s}{o_t}}\,h_{ts}\sigma\,
V_t^T D_\lambda(g)V_s.
\]

The total multiplicity-space dimension is m_lambda=sum_s r_s. Each eigenvalue
of that reduced block occurs d_lambda times in the original Hilbert space,
and

\[
\sum_\lambda d_\lambda m_\lambda=\dim\mathcal H_L
\]

must hold. The implementation checks this complete dimension identity.
Young orthogonal matrices are generated from adjacent transpositions and
checked against independently chosen group products.

Independent full matrices at (N,L)=(3,2) and (4,1) verify every eigenvalue,
including multiplicities. The latter has 1,875 states and includes a sector
whose minimum is positive. A full 427,986-state N=5,L=3 sparse matrix
independently checks the lowest eigenvalues and ground observables.
For N=6,L=3 all eleven irreducible sectors are solved; together they account
for 36,801,265 states. Floating Lanczos residuals are numerical checks,
not interval-certified eigenvalue lower bounds.

## 3. Connectedness bounds that retain negative hopping and graph mixing

Let C project onto connected graphs and D=I-C. Work either in the complete
physical Hilbert space or in a declared invariant subspace preserved by C.
Write

\[
H=\begin{pmatrix}H_C&B\\B^\dagger&H_D\end{pmatrix},
\quad E_C=\inf\sigma(H_C),\quad E_D=\inf\sigma(H_D),
\quad b=\|B\|.
\]

Each compression retains the full rewiring diagonal. Rebuilding a Laplacian
using only internal transitions would incorrectly remove the boundary cost.
The supplied code therefore slices the complete matrix.

For a normalized ground eigenvector psi=(u,v) of energy E below both block
minima, set p=||u||^2, delta_C=E_C-E, and delta_D=E_D-E. Its two projected
eigen-equations give

\[
\delta_C\sqrt p\le b\sqrt{1-p},\qquad
\delta_D\sqrt{1-p}\le b\sqrt p.
\]

Consequently,

\[
\boxed{\frac{\delta_D^2}{\delta_D^2+b^2}\le p\le
\frac{b^2}{\delta_C^2+b^2}.}
\]

Also, a direct quadratic-form comparison yields

\[
\boxed{\frac{E_C+E_D-\sqrt{(E_C-E_D)^2+4b^2}}2
\le E\le\min(E_C,E_D).}
\]

These statements are exact for the specified operators. Their numerical
evaluation here uses floating estimates of the block minima and cross norm.
For N=6 the reported block values are those of the invariant sector.
This suffices for its ground vector, and the all-irrep calculation separately
identifies that vector as the full regulated ground. We do not silently
identify every invariant-sector block minimum with its full-space analogue.

The graph-mixing contribution is

\[
E_{\rm mix}=2\operatorname{Re}\langle u,Bv\rangle .
\]

Thus E=p E_{C|psi}+(1-p)E_{D|psi}+E_mix. Dephasing only the binary
connected/disconnected information changes the energy by -E_mix.
This is a mathematical record/dephasing operation, not a supplied physical
apparatus. It makes explicit why connected and disconnected variational
energies cannot be compared while omitting their coherent coupling.

At Q=N a connected graph has one cycle. Leaf pruning gives
ell(G)=N-b(G), where b(G) is its number of bridges. We compute the full
conditional distribution of this cycle length in the same ground state.
Connectedness by itself does not diagnose an extended geometry.

## 4. A lower bound for the omitted, unbounded electric sector

For a fixed graph, group the hopping terms by their vertex endpoint. At
vertex i the dressed edge annihilator

\[
B_i=\sum_{e\ni i}s_{ie}U_{ie}b_e
\]

obeys {B_i,B_i^dagger}=d_i. Together with the vertex fermion it gives the
two-mode operator 2 n_i+a_i^dagger B_i+B_i^dagger a_i with lower bound

\[
1-\sqrt{1+d_i}.
\]

Different stars need not commute: summing their operator lower bounds is
still valid. Since sum_i d_i=2N, concavity of the square root implies

\[
\tfrac12Q+2N_v+T
\succeq N\bigl(\tfrac32-\sqrt3\bigr).
\]

This bound holds graph by graph and hence on graph superpositions. It uses
the untruncated unitary U(1) shifts. Positive R may be dropped when making
a lower bound.

Outside P_L at least one leg has magnitude >=L+1. Edge Gauss law fixes the
sum of the paired legs to 0 or 1, so their squared fields sum to at least
(L+1)^2+L^2. The complete outside compression C_L therefore obeys

\[
\boxed{C_L\succeq c_L I,\qquad
c_L=N(\tfrac32-\sqrt3)+\frac{(L+1)^2+L^2}{8}.}
\]

For rational certificates replace sqrt(3) by 1732051/1000000, whose square
is greater than 3. This gives a slightly weaker, exactly rational lower
bound, not a floating approximation passed off as a proof.

## 5. The cutoff coupling has an exact diagonal square

Let B_L=(I-P_L)H P_L. Only fermion hopping can cross the electric cutoff:
rewiring leaves the fields unchanged. A reachable outside state has exactly
one out-of-cutoff leg. Its value is L+1; a negative value -(L+1) would,
by edge Gauss law, require its partner to be outside too.

The unique bad leg specifies the hopping term and its unique predecessor.
Different inside states therefore cannot reach the same outside state in
one hop. It follows, at frozen unit hopping, that

\[
\boxed{B_L^\dagger B_L=D_L,}
\]

where D_L is diagonal and counts allowed outward hopping moves from each
inside physical basis state. This is stronger than estimating B_L by the
total number of Hamiltonian terms. The coefficient is invariant under
event permutations, so D_L is equally simple in each multiplicity block.

For any z<c_L, completing the square in the outside component gives the
sufficient lower-energy certificate

\[
\boxed{A_L-zI-\frac{D_L}{c_L-z}\succeq0
\quad\Longrightarrow\quad H\succeq zI.}
\]

Equivalently one may use the Schur complement and
(C_L-z)^(-1)<=1/(c_L-z). This is an application of established
Feshbach-Schur comparison, not a new spectral theorem.

At N=3,L=2 and L=3 the finite comparison matrix has rational entries.
Exact Fraction Schur elimination gives strictly positive pivots, proving
positivity without relying on a floating eigenvalue sign. An integer-valued
trial vector, suggested by a floating eigenvector but evaluated with an exact
rational Rayleigh quotient, supplies a rigorous upper bound for H.

The tighter result is

\[
-0.46729956\le E_0(N=3,Q=3)
\le-\frac{18691982041980411}{39999999816851996}.
\]

The upper endpoint is approximately -0.4672995531891348. This is a
certificate for the unbounded-rotor Hamiltonian at N=3, not an assertion
about finite-cutoff convergence at every N.

For N=5 and N=6 we also evaluate the same comparison matrix in every
permutation sector, at z=-0.728494 and z=-0.85832 respectively. These larger
computations use floating smallest-eigenvalue searches. They are explicitly
recorded as numerical positivity evidence; exact rational positivity has
not been certified for those larger blocks.

## 6. Limits and provenance

The finite connectedness probabilities do not determine their N->infinity
limit. Large component fraction, complete connectedness, local regularity,
and a physical propagation law are different observables. No spacetime
continuum, massless graviton, or universal matter coupling follows here.
An N=7,L=2 invariant result is retained as provisional because neither its
other permutation sectors nor its cutoff convergence were checked.

The Schur method has established antecedents, for example G. Dusson,
I. M. Sigal, and B. Stamm, “The Feshbach-Schur map and perturbation theory,”
[arXiv:2105.02058](https://arxiv.org/abs/2105.02058).
The CP-12-specific content is its physical-state decomposition, connected
sector evaluation, and the explicit electric-tail/outward-hop certificate.
