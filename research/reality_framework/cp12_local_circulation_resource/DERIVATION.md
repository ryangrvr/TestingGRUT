# CP-12: local circulations, graph entropy, and preparation resources

2026-10-14. Additive derivation under the frozen CP-12 Hamiltonian. Based on the
holonomy audit `ae91c9beb0d8dda0f7692ea8bab0c5f808355530`, which contains the ground
connectivity checkpoint `c3c3de35729b175aab01fe45624ae65a10ad0e9d`.
This proof and its executable controls have not received an independent partner audit.

## 1. Ensemble and law

Use the complete **labelled**, physical U(1)-Gauss-law Hilbert space at fixed
N vertices and Q occupied connections, with M=N(N-1)/2. All integer electric
fields are allowed. No projection onto an event-permutation singlet sector is
made in the thermal or preparation theorems. A permutation-invariant observable
is not the same as a quotient equilibrium ensemble.

Write

\[
H=A+V,\qquad A=\kappa Q+2mN_v+\frac{\epsilon_E}{2}\sum_{e,i\in e}E_{ie}^2,
\qquad V=T_{\rm hop}+R_{\rm rewire}.
\]

Frozen parameters: m=t_hop=1, kappa=1/2, epsilon_E=1/4, Gamma=1/10.
R is the positive neutral-connection graph Laplacian. Its off-diagonal swap
amplitude is -Gamma/M, and its diagonal counts **all** permitted swaps, including
those leaving a chosen graph family. Removing that diagonal would change the law.

For every physical basis state, at most Q(M-Q) swaps are permitted. Consequently

\[
0\le R\le\frac{2\Gamma Q(M-Q)}M I,
\quad \|T_{\rm hop}\|\le2|t_{\rm hop}|Q.
\]

The hopping bound follows by summing the norms of the 2Q Hermitian incidence
hops. At each incidence the fermion transfer, multiplied by its unitary electric
shift, is a partial isometry plus adjoint of norm at most one. The rewiring bound
follows from the maximum weighted degree of a positive graph Laplacian, also
valid for the countably infinite physical basis. Thus

\[
v_-I\le V\le v_+I,\qquad
v_-=-2|t_{\rm hop}|Q,\quad
v_+=2|t_{\rm hop}|Q+2\Gamma Q(M-Q)/M.
\]

Define the usable spectral-width bound

\[
\boxed{\Omega_N=v_+-v_-=4|t_{\rm hop}|Q+2\Gamma Q(M-Q)/M.}
\tag{1}
\]

For Q=nu N+O(1), fixed nu>0, Omega_N=O(N). This is not a small-coupling expansion.

## 2. Global cycles versus local gauge generators

The audit's Gauss-incidence result is correct: at fixed graph and fermion
occupations the affine integer circulation lattice has rank Q-N+C(G), provided
each component has zero net Gauss charge. For a connected Q=N graph this is one.

For a d-dimensional periodic nearest-neighbor grid with all circumferences
at least five, Q=dN. Its elementary square boundaries span a rational cycle
subspace of dimension

\[
r_{\rm faces}=Q-N+1-d=(d-1)N+1-d.
\tag{2}
\]

The remaining d independent classes are torus winding circulations. This follows
from the first homology of the cellulated torus; here it is also verified by
**exact Fraction-arithmetic elimination of the actual boundary vectors**, with
no floating-point rank tolerance:

| Supplied host | N | Q | Cycle rank | Elementary-face span | Remaining winding rank |
|---|---:|---:|---:|---:|---:|
| 6 x 6 torus | 36 | 72 | 37 | 35 | 2 |
| 5 x 5 x 5 torus | 125 | 375 | 251 | 248 | 3 |

These are capacities of prepared hosts, not a ground-state distribution or a
claim that the Hamiltonian selects dimension d. Algebraic independence of flux
generators must not be confused with the number of independently excitable
Wilson-loop states: products of loop shifts can satisfy relations even though
individual one-loop flux states are orthogonal.

### A full-law physical response check

On a supplied simple cycle c of length ell, let W_c shift the oriented integer
electric flux by one circulation, with no change of fermion occupations. Set it
to zero on graphs missing any cycle edge. W_c is invariant under the internal
U(1) gauge transformations. A choice of a particular cycle still supplies a
reference; we do not derive a localized observer or an event-label gauge fixing.

In the unregulated theory define

\[
X_c=(W_c+W_c^\dagger)/\sqrt2.
\]

Let |G,0> be the neutral sea on a definite graph containing c. X_c|G,0> has unit
norm and consists of the two opposite flux orientations. The diagonal energy
of each flux state exceeds that of the neutral state by

\[
\Delta_c=\ell\left[\epsilon_E-\frac{\Gamma(M-Q)}M\right].
\tag{3}
\]

The first term is the electric cost: each occupied cycle edge has two unit fields
and therefore costs epsilon_E. The second term is the decrease in the rewiring
diagonal when those ell edges are locked. Mass and connection costs cancel.
Hopping changes occupations and has no matrix element between the three states
used in this expectation. Rewiring off-diagonal terms change the graph and have
zero overlap with them. They have not been dropped from H.

The exact double-commutator identity gives

\[
\boxed{\langle G,0|[X_c,[H,X_c]]|G,0\rangle=2\Delta_c.}
\tag{4}
\]

With the response convention i<[X_c(T),X_c]>, its expansion is

\[
R_c(T)=2\Delta_c T+O_N(T^3).
\tag{5}
\]

The absence of the even term follows from the real Hamiltonian, real observable,
and real neutral preparation. If the supplied graph lacks c, the first moment
is zero because X_c annihilates that state. These are prepared-state moments,
not equilibrium spectral densities. A nonzero first moment neither proves a
mass gap nor excludes a low-frequency collective pole.

The full transition generator was applied once to the neutral state and both
flux states on the two hosts above, retaining every destination and fermion
sign. The exact slopes are 226/175 and 192/155. On the cubic host the neutral
state has 2,765,625 rewiring outputs and 750 hopping outputs; these were streamed,
not stored as a huge diagonalization. All source/flux states obey Gauss's law.
A separate complete 1,875-state (N,Q,L)=(4,4,1) matrix and finite source impulses
recover slope 26/15 with quadratic time error. Its first moment is unregulated
exact because the contributing matrix elements fit within L=1; the later-time
finite-impulse data remain cutoff calculations.

## 3. A graph family large enough to avoid the template loophole

Fix a microscopic cycle-length bound R>=3. Let K_R(G) be the maximum number of
**vertex-disjoint** simple cycles of lengths 3 through R. Define

\[
\mathcal S_{N,Q}(R,k)=\{G: |E(G)|=Q,\ K_R(G)\ge k\}.
\]

This family has no embedding, lattice template, homogeneity requirement, or
preferred dimension. It includes all graphs with the stated loop packing,
including many nonspatial graphs. Its projector retains **all** physical matter
and electric states on each included graph. It is invariant under relabelling.

Let U_{N,Q} be uniform graph counting on all C(M,Q) labelled graphs, solely as an
intermediate counting device. Define

\[
A_{N,Q,R}=\sum_{r=3}^{R}\frac{(NQ/M)^r}{2r}.
\]

Then

\[
\boxed{U_{N,Q}(\mathcal S_{N,Q}(R,k))\le
\min\left(1,\frac{A_{N,Q,R}^{\,k}}{k!}\right).}
\tag{6}
\]

**Proof.** For an ordered cycle-length list r_1,...,r_k, let ell=sum r_i. The
number of collections of disjoint unoriented cycles, summed over ordered length
lists and divided by k!, is (N)_ell/[k! product_i(2r_i)]. A specified collection
requires ell different edges. Its fixed-Q inclusion probability is
(Q)_ell/(M)_ell <= (Q/M)^ell. Replace (N)_ell by N^ell and sum over length lists.
The expected number of such packings is at most A^k/k!. If any packing exists,
its count is at least one. This proves (6). Lists with ell>N or ell>Q contribute
zero to the exact expression; retaining them in the upper sum is harmless.

The code independently enumerates complete graph ensembles and all allowed
packings, and checks their sum against

\[
\sum_{r_1,...,r_k=3}^{R}\frac{(N)_\ell}{k!\prod_i2r_i}
\frac{\binom{M-\ell}{Q-\ell}}{\binom MQ}.
\tag{7}
\]

For Q=nu N+O(1), fixed nu and R, A=O(1). If k=ceil(alpha N), alpha>0,

\[
\log U_{N,Q}(\mathcal S)\le-\alpha N\log N+O(N).
\tag{8}
\]

This is stronger in scope than assuming |S_N|<=N! exp(O(N)): it obtains the
rarity directly from the presence of many bounded-length local cycles. The
counting is elementary random-graph mathematics, not a new physical principle.

### Why the condition covers a conventional local gauge lattice

Square and cubic tori with even first two circumferences can be tiled by N/4
vertex-disjoint elementary squares. The code explicitly constructs these tiles.
Their loop density therefore satisfies alpha=1/4, including Q=2N and Q=3N.

More generally, suppose maximum degree <=D and at least fN vertices lie on
cycles of length <=R. Select a maximal vertex-disjoint family. Every other
short cycle intersects a selected one. Every covered vertex is within graph
distance floor(R/2) of an intersection vertex. With

\[
B_D(r)=1+D\sum_{j=0}^{r-1}(D-1)^j,
\]

a selected cycle accounts for at most R B_D(floor(R/2)) covered vertices. Thus

\[
K_R(G)\ge\frac{fN}{R B_D(\lfloor R/2\rfloor)}.
\tag{9}
\]

So any fixed positive density of bounded-length circulation paths in a uniformly
bounded-degree graph is covered, with loose but N-independent constants. The
condition is sufficient for this obstruction, not a universal definition of
physical space. It does not impose finite valency on Lorentzian causal links.

## 4. Direct unregulated quantum Gibbs bound

For fixed beta>0 let sigma_beta=exp(-beta A)/Z_A and
rho_beta=exp(-beta H)/Z_H. A has trace-class heat kernel at finite N because
its electric energy is positive quadratic and its fermion and graph factors
are finite. Bounded V preserves self-adjointness and trace class. No electric
cutoff or limit exchange is needed here.

For one graph the active N+Q fermion modes have exactly Q occupied states;
absent edge modes are filled and have zero fields. Ignoring Gauss constraints
only overcounts. With Theta(s)=sum_{z in Z} exp(-s z^2),

\[
Z_A(G)\le e^{-\beta\kappa Q}\binom{N+Q}{Q}
\Theta(\beta\epsilon_E/2)^{2Q}.
\]

Every graph also supplies one neutral sea with A-energy kappa Q. Hence
Z_A>=C(M,Q)exp(-beta kappa Q). Define

\[
F_N(\beta)=\log\binom{N+Q}{Q}+2Q\log\Theta(\beta\epsilon_E/2),
\quad
\mathcal L_N=\log k!-k\log A_{N,Q,R}-F_N(\beta).
\tag{10}
\]

Equations (6) and (10) give
sigma_beta(S)<=exp(-L_N) whenever L_N>0.
For fixed nu,beta,R and k=ceil(alpha N), L_N=alpha N log N+O(N).

The bounded-perturbation Gibbs identity is

\[
D(\rho_\beta\|\sigma_\beta)
=-\beta\langle V\rangle_{\rho_\beta}+\log Z_A-\log Z_H
\le\beta\Omega_N.
\tag{11}
\]

For the last inequality, v_-<=V<=v_+ gives Z_H>=exp(-beta v_+)Z_A by eigenvalue
min-max comparison, while <V>>=v_-. This is a trace statement; it does **not**
use the false general inference A<=H implies exp(-beta H)<=exp(-beta A).

Measurement of {Pi_S,1-Pi_S} cannot increase quantum relative entropy. If
p=Tr(rho Pi_S), q=Tr(sigma Pi_S), its binary value satisfies

\[
D(\rho\|\sigma)\ge d(p\|q)\ge p\log(1/q)-\log2.
\tag{12}
\]

No commutation of Pi_S with H is assumed. Combining (10)-(12),

\[
\boxed{P_\beta(\mathcal S_{N,Q}(R,k))\le
\min\left(1,\frac{\beta\Omega_N+\log2}{\mathcal L_N}\right)}
\quad(\mathcal L_N>0).
\tag{13}
\]

For extensive k this is O(1/log N). The theorem retains full fermions, all
unbounded rotors and coherent rewiring. Unlike the older unregulated extension,
there is no electric low/high split and no N log log N rank estimate. The same
direct argument also simplifies that older narrow-family theorem.

A cooling consequence: if beta_N>=beta_min>0 but beta_N=o(log N), F_N=O(N)
uniformly and (13) still tends to zero. To escape this sufficient obstruction
by cooling alone requires beta_N at least of logarithmic order. That is **not**
a claim that logarithmic cooling suffices to select a spatial ground state.

The finite-N constants are very loose. Some small-system bounds exceed one
before clipping. This theorem must not be used as a useful numerical prediction
at N=5 or N=6 or as a contradiction of their connected ground-state weights.

## 5. Nonequilibrium does not automatically evade the obstruction

Use any one fixed beta>0 to define sigma_beta. For a density matrix rho_0 of
finite reference relative entropy and energy, evolve with the **full closed**
CP-12 Hamiltonian: rho_T=exp(-iHT)rho_0 exp(iHT). Entropy and <H> are conserved,
so exactly

\[
\boxed{D(\rho_T\|\sigma_\beta)
=D(\rho_0\|\sigma_\beta)
+\beta(\langle V\rangle_0-\langle V\rangle_T)
\le D(\rho_0\|\sigma_\beta)+\beta\Omega_N.}
\tag{14}
\]

Apply the same binary measurement inequality:

\[
\boxed{\sup_T P_T(\mathcal S)\le
\min\left(1,\frac{D(\rho_0\|\sigma_\beta)+\beta\Omega_N+\log2}
{\mathcal L_N}\right).}
\tag{15}
\]

If the initial D is O(N), this is O(1/log N), **uniformly over all evolution
times**, including times growing with N. No approximation to rewiring or matter
evolution is involved. This is stronger than a short-time graph departure bound,
but its preparation assumption is essential.

One explicit preparation is the incoherent neutral mixture

\[
\rho_{\rm mix}=\binom MQ^{-1}\sum_G|G,0\rangle\langle G,0|.
\]

Its entropy is log C(M,Q), its A-energy is kappa Q, and

\[
0\le D(\rho_{\rm mix}\|\sigma_\beta)
=\beta\kappa Q+\log Z_A-\log\binom MQ\le F_N(\beta)=O(N).
\tag{16}
\]

It cannot assemble an extensive collection of bounded-length, disjoint loops
with nonvanishing probability at any time in the thermodynamic limit.

The coherent neutral state

\[
|\Omega\rangle=\binom MQ^{-1/2}\sum_G|G,0\rangle
\]

has the same graph probabilities and A-energy, but zero von Neumann entropy.
At preparation, exactly

\[
D(|\Omega\rangle\langle\Omega|\|\sigma_\beta)
-D(\rho_{\rm mix}\|\sigma_\beta)=\log\binom MQ
=\nu N\log N+O(N).
\tag{17}
\]

Both preparations have extensive full energy, but these energies are not identical:
<E>_mix=kappa Q+Gamma Q(M-Q)/M, whereas <E>_coherent=kappa Q. The coherent
neutral vector is annihilated by the rewiring Laplacian. The O(N log N) entropy
difference is therefore not an unnoticed superextensive energy investment.
These energy identities are independently checked in the complete small matrices.

Thus it escapes (15). This is an **available information resource**, not proof of
assembly, phase locality, or a physical preparation mechanism. It supplies
coherence across the enormous graph ensemble but no specified spatial ordering.
Escaping this entropy bound also does not cool the system: the full energy remains
conserved, and a neutral coherent sea is not a prepared ground state.
A low-entropy incoherent, specially prepared graph can also have large D; the
resource is entropy deficit relative to the reference, not uniquely coherence.

Conversely, if P_T(S)>=p_*>0 for extensive k, (15) requires

\[
D(\rho_0\|\sigma_\beta)\ge
p_*\alpha N\log N-O(N).
\tag{18}
\]

This converts the preparation question into a quantitative obligation. A thermal
reservoir, external work or entropy export would change the closed-system
premises and must be modeled if used. Retaining graph-record coherence is not
an automatic entropy sink.

## 6. What is and is not resolved

- Higher density gives the required *capacity*: explicit extensive independent
  local circulations and gauge-loop response exist on supplied Q=2N and Q=3N hosts.
- Full equilibrium selection of a broad bounded-loop class fails at fixed
  positive temperature in the stated labelled CP-12 ensemble.
- Closed evolution from the specified high-entropy preparations fails even at
  arbitrarily late times. A low-entropy or coherently purified preparation remains open.
- Zero-temperature ground states, an event-permutation quotient, R increasing
  with N, unbounded-degree operational geometries, relational coarse-graining,
  external cooling and different physical ensembles are not eliminated here.
- No two-polarization photon, relativistic light cone, three-dimensional spatial
  phase, graviton, Born-rule derivation or new empirical law is established.

The useful constructive target is now a **declared low-entropy preparation under
one Q-density sector**, tested for local cycle packing, intrinsic volume growth
and gauge-invariant propagation in the same evolved state. Supplying such a
resource only makes the test admissible; it does not make its outcome favorable.

## 7. Literature and priority

The graph-entropy difficulty is established in quantum graphity: T. Konopka,
[Statistical Mechanics of Graphity Models](https://arxiv.org/abs/0805.2283),
Phys. Rev. D **78**, 044032 (2008), DOI
[10.1103/PhysRevD.78.044032](https://doi.org/10.1103/PhysRevD.78.044032).
We do not claim the graph-entropy mechanism or elementary cycle counting as novel.

Quantum relative-entropy data processing is established mathematics. A primary
source covering separable Hilbert spaces is A. Mueller-Hermes and D. Reeb,
[Monotonicity of the Quantum Relative Entropy Under Positive Maps](https://arxiv.org/abs/1512.06117),
Ann. Henri Poincare **18**, 1777-1788 (2017), DOI
[10.1007/s00023-017-0550-9](https://doi.org/10.1007/s00023-017-0550-9).
The calculations here specialize these tools to CP-12's physical constrained
Hilbert space; they are not a new selection principle.
