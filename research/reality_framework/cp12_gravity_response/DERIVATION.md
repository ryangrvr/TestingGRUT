# CP-12: local back-reaction, graph harmonics, and the tensor target

This is an additive analysis of frozen CP-12 at the parent research commit
`553aeffc10faef6928a7646d4a9173019273e2d7`. No coupling is changed.
The Hamiltonian and physical constraints are those of
`research/constructive_physics/attempt12/LAW.md` and `model.py`.

## 1. Definitions and what is being compared

There are N vertices, M=N(N-1)/2 potential connections, Q occupied
connections, and compact electric fields on the two incidence legs of each
present edge. The graph-changing off-diagonal term is V. It has matrix
element -Gamma/M between allowed neutral rewiring states. H_B=H-V retains
every fixed-graph block's **actual** rewiring diagonal.

The frozen constants are m=t=1, kappa=1/2, epsilon_E=1/4, Gamma=1/10.
Electric cutoffs regulate numerical calculations; they are not changes to
the law. Exact finite-step identities and analytical inequalities below
do not use a weak-hopping expansion.

Three different constructions must remain distinct:

* A definite or coherent prepared host supports an operational response
  comparison, but its geometry is supplied.
* A graph-harmonic projection characterizes connection coherence. It is
  not an added physical constraint and is not an invariant sector of H.
* A signed event-permutation quotient gives exact symmetry-sector spectra.
  Treating event relabelling as gauge is an additional physical interpretation.

Time T remains the microscopic external Hamiltonian parameter. No clock,
Lorentz symmetry, stress tensor, or Einstein constraint is derived here.

## 2. Exact leading connection response to a matter intervention

Prepare the neutral sea on one definite graph F: all vertex fermions are
empty, every edge fermion is filled, and every electric field is zero.
Choose vertex i with occupied degree d_i>0. Let B_i be the sum of **all
frozen hopping terms at i**, divided by sqrt(d_i). Presence factors in H
activate occupied edges. The source does not inspect an independent
neighborhood pool. Its norm on the initial graph sector is |t|=1.

For the connection readout n_e define

    R_ei(T) = < i[n_e(T), B_i] >.

Then the exact cubic coefficients are

    e present, incident on i:
      R_ei(T) = -Gamma^2 t^2 (M-Q)/(M^2 sqrt(d_i)) T^3 + O_N(T^5);

    e present, not incident on i:
      R_ei(T) = O_N(T^5);

    e absent:
      R_ei(T) = +Gamma^2 t^2 sqrt(d_i)/M^2 T^3 + O_N(T^5).

The final line holds even when the absent edge is not incident on i.
No physical distance enters that coefficient. Summing over all edges
gives zero, as required by conserved Q.

Here is the path-counting reason. A graph-diagonal readout must encounter
two rewiring amplitudes to contribute to this expectation in a definite
graph. The source changes the matter state by one hopping operation, so
one further Hamiltonian hopping operation is needed. At cubic order,
terms that remove an edge unrelated to the source cancel. For an incident
edge, the source-created charge disables removal in one of the two
operator orderings. Every M-Q possible replacement then contributes the
same factor (Gamma/M)^2 t^2/sqrt(d_i). For one specified absent edge,
the d_i possible incident removals supply the positive coefficient.
Diagonal mass, electric, and rewiring energies require an extra operation
and cannot affect this coefficient. All even response coefficients vanish
because H, B_i, n_e and the initial vector are real in this basis.

`exact_short_time.py` assembles Hamiltonian steps independently with exact
integers, including the fermionic and incidence signs. It enumerates every
reachable unregulated state through the needed order, checks Gauss's law,
and obtains these rational coefficients at N=4,...,7. Independently,
complete frozen-builder evolutions and source-pulse derivatives check the
response on a path, a cycle, and a star.

For Q=nu N and bounded d_i, a specified existing connection has cubic
coefficient O(N^-2); a specified absent connection has coefficient
O(N^-4). This is a short-time statement, not a uniform Taylor expansion at
large N or a measured Newton constant.

## 3. A nonperturbative bound for one local connection

Now allow arbitrary physical matter/electric states, initially on the
single graph F, with |F|=Q. Write L=M-Q, and k for occupied edges outside F.
The previous full-law pool theorem gives

    <k(T)> <= g^2 T^2,      g = Gamma Q sqrt(L)/M.

For a specified initially present edge e, split the transitions that remove
it into replacements inside F and replacements outside F. Their k-sector
operator norms are bounded by

    Gamma/M sqrt((Q-k) k),
    Gamma/M sqrt((L-k)(k+1)),

respectively. These are rectangular inclusion-matrix bounds. Fixed
fermion/electric fibers and the neutral-removability restriction produce
submatrices of these matrices and cannot increase their norm.

Integrating the off-diagonal forcing in the n_e=0/1 decomposition, using
Cauchy-Schwarz and the <k> theorem, bounds the amplitude to lose e by

    a_in(T) = Gamma/M [ sqrt(Q) g T^2/2 + sqrt(L) I_g(T) ],

where

    I_g(T) = integral_0^T sqrt(1+g^2 s^2) ds
           = [T sqrt(1+g^2 T^2) + asinh(g T)/g]/2.

For a specified initially absent edge f, removing an inside edge has norm
Gamma/M sqrt((Q-k)(k+1)), and removing another outside edge has norm
Gamma/M sqrt(k(L-k)). Thus

    a_out(T) = Gamma/M [ sqrt(Q) I_g(T) + sqrt(L) g T^2/2 ].

The corresponding occupation-change probabilities are at most a_in^2 or
a_out^2, capped at one. These are full-H inequalities; matter hopping is
not switched off.

For any graph-diagonal source B with initial-sector norm <=1, the linear
response of the specified connection is at most 2 a_in^2 or 2 a_out^2.
To see this, evolve both psi and B psi. The initially constant part of n_e
contributes a real conserved inner product. The remaining Kubo matrix
element involves two leaked components, each bounded by the respective
amplitude. B's possibly larger norm on other graph sectors is irrelevant.

At fixed Q/N, a specified existing edge therefore has response upper bound
O((T+T^2)^2/N^2). On the supplied three-dimensional ballistic scale
T~N^(1/3), it tends to zero as O(N^-2/3). A bounded, normalized average of
such connection readouts inherits this suppression. This sharpens the
earlier whole-state estimate.

It does **not** prove that every possible emergent metric has zero
gravitational coupling. An emergent metric might involve gauge or matter
correlations, an extensive source, a different state/limit, or a nonlocal
relational encoding. It shows that a finite local metric back-reaction
cannot simply be identified with these bounded connection readouts in the
certified prepared window.

## 4. All-time neutral geometry retention and its missing dispersion

Isolate the unchanged neutral rewiring term R. This is a diagnostic sector,
not full CP-12 evolution: t!=0 creates matter excitations from the neutral
sea. R is Gamma/M times the Laplacian of the Johnson graph J(M,Q).

For Q<=M/2 its eigenvalues and multiplicities are

    lambda_j = Gamma/M j(M+1-j),
    multiplicity_j = binom(M,j)-binom(M,j-1),   j=0,...,Q.

One direct derivation uses f_A(X)=1_{A subset X}, |A|=j, on Q-subsets X.
The unscaled Laplacian acts as

    L f_A = j(M+1-j) f_A - (Q-j+1) sum_{a in A} f_(A minus a).

Successive inclusion spaces therefore supply the stated eigenvalues; the
dimensions give the multiplicities. Vertex transitivity makes the spectral
weight of any one definite graph multiplicity_j/binom(M,Q).

Its highest-energy sector has weight

    w_Q = 1-alpha,       alpha = Q/(M-Q+1).

For alpha<=1/2, the reverse triangle inequality gives, at **every** time,

    P(depart the initial graph) <= 4 alpha(1-alpha).

For Q=nu N this is O(1/N), uniformly in T. This is a stronger result than
the fixed-time perturbative survival estimate. The definite preparation
is concentrated near the top of this neutral spectrum, not its ground
state. Coherent preparations need not have this spectral concentration.
The bound must not be extended to full matter/gauge dynamics.

Within a supplied pool of S=Q+1 short edges, there is one vacancy. The
conditional neutral kinetic operator, after its scalar diagonal is
removed, is

    H_vac = Gamma/M ( S I - 1 1^T ).

Every nonuniform vacancy mode has the same frequency Gamma S/M. If the
pool has a translation structure, all nonzero momenta still share this
frequency. A localized vacancy has departure probability at most
4(S-1)/S^2 at any time. A time rescaling does not turn this flat spectrum
into a linear or quadratic spatial dispersion. Matter/gauge correlations
would have to generate that dispersion.

## 5. Exact charged-bond transitions between graph harmonics

Fix the fermion bits and electric fields. If S_l edges are locked by a hole
or nonzero electric field, K=M-S_l positions remain eligible and
r=Q-S_l neutral edges are occupied. The highest Johnson harmonic is

    W_r^K = ker D_r,

where D_r sums r-subset amplitudes over extensions of each (r-1)-subset.
The following statements concern this graph harmonic, not spatial angular
momentum or the event-permutation representations in section 7.

Charge-locking an already occupied neutral bond restricts an r-subset to configurations
containing a specified eligible edge. For f in ker D_r, its restricted
amplitude f_e(Z)=f(Z union {e}) satisfies D_(r-1) f_e=0. Thus locking maps
the highest old harmonic exactly into the highest new harmonic.

Unlocking is the adjoint inclusion: add one eligible position e and
support the lifted amplitude only on subsets containing e. For an old
normalized highest harmonic f in W_r^K,

    ||(I-P_top,new) lift(f)||^2 = 1/(K+1-2r).

The surviving top-to-top transport has squared singular value
1-1/(K+1-2r). One proof is to apply the new down map. It gives f extended
by zero on r-subsets containing e, a highest harmonic of degree r. The
non-top component of the lift therefore lies only in the new degree-r
harmonic. On this component, D_(r+1)^dagger D_(r+1) has eigenvalue
K+1-2r, which fixes its norm.

At extensive sparse Q, K is O(N^2) and r is O(N). Leakage per unlock
amplitude is O(1/N), and the top-to-top hopping amplitude differs from t
by O(1/N^2). This is a matter/connection coherence relation generated by
the declared microscopic rule, not a fitted transport kernel.

The sum of many channels is important. It is invalid to call the whole
sector invariant from this small per-channel factor. There are at most
Q unlock channels per source physical fiber and at most 2M incoming
channels per target fiber. Block Schur counting yields the safe full-H bound

    ||(I-P_top) H P_top||^2 <= 2 M Q t^2/(M-2Q+1),  Q<=M/2.

This grows as O(N); it is not a vanishing global leakage bound.

Nevertheless it gives a useful ground/preparation separation. On each
fiber,

    P_top V P_top = Gamma r/M P_top >=0.

The pinned parent proved H_B >= N ell_N at Q=2N, and supplied a full-H
matching-state trial upper bound U_N, with

    ell_N -> -0.03262004145,
    U_N/N -> -0.08893517460.

For delta_N=N ell_N-U_N>0, project the full ground eigenvalue equation onto
P_top. With b_N^2 the preceding crossing bound,

    ground probability in P_top <= b_N^2/(delta_N^2+b_N^2) = O(1/N).

Conversely, every physical state on any one definite graph has

    probability in P_top >= 1-Q/(M-Q+1) = 1-O(1/N).

These full-law analytical statements use the parent's unregulated bounds,
not a small-cutoff ground extrapolation. The accompanying script evaluates
them with outward decimal arithmetic. They identify distinct connection
coherence structures for the absolute ground state and a definite prepared
geometry. They do not select a dimension or prove a prepared-sector lifetime.

There is a record consequence of this exact V identity. A retained ideal
full-graph record dephases the graph basis, removing V while leaving H_B.
For a state supported in P_top its system-energy change is

    Delta E_record = -Gamma <r>/M,   |Delta E_record| <= Gamma Q/M.

At sparse extensive Q this is O(1/N), and its sign is nonpositive. This
contrasts with the extensive positive record cost established for the
absolute ground state. It is not a free-energy or total-apparatus claim.
The complete N=4 control confirms a -0.05 system-energy shift for a
normalized projected neutral preparation.

Projecting observables also makes some originally commuting geometry and
source operators fail to commute. For commuting Hermitian A and B,

    [P A P, P B P] = -P A (I-P) B P + P B (I-P) A P.

For one connection A=n_e in P_top, its non-top occupation amplitude has
norm squared at most 1/(M-2Q+2). On a fiber with eligible K,r its exact
nonzero eigenvalue is s(1-s), s=1-1/(K-2r+2). A full vertex hopping source
normalized by sqrt(N-1) has top-crossing norm at most
|t| sqrt(N-1)/sqrt(M-2Q+1), by the local unlock-channel Schur count.
Consequently the projected commutator norm is at most

    2|t| sqrt(N-1)/sqrt[(M-2Q+2)(M-2Q+1)] = O(N^-3/2).

The complete control exhibits nonzero finite-size commutators and checks
the complementary-projector product bound. This particular nonlocal
projection therefore does not automatically supply a finite gravitational
connection algebra in the sparse limit. It is graph-harmonic dressing,
not the energy dressing required by Einstein gravity.

## 6. Mean bonds miss tensor components; physical correlations restore them

Use four supplied unit tetrahedral directions d_a. A weighted diamond
Laplacian has lower Bloch band

    lambda_-(k)=sum_a w_a-|sum_a w_a exp(i k.d_a)|.

Its quadratic stiffness is

    C(w) = 1/2 [sum_a w_a d_a d_a^T
                - (sum_a w_a d_a)(sum_b w_b d_b)^T/sum_a w_a].

At w_a=1, the linear stiffness variations are d_a d_a^T/2. They span only
four of the six symmetric-tensor components: a trace and three off-diagonal
shears. The diagonal traceless two-dimensional component is absent. For a
wave vector along a coordinate axis, this construction has only one
transverse-traceless response direction. This is a kinematic limitation of
the **supplied weighted-Laplacian identification**, not CP-12's derived
continuum dispersion.

Six pair tensors (d_a-d_b)(d_a-d_b)^T span all six components. CP-12 already
contains physical operators with the corresponding tetrahedral symmetry.
Let u_a be the bounded indicator of nonzero electric field on a central
incidence leg. The four u_a contain a singlet and a triplet. The six u_a u_b
contain a singlet, a doublet E, and a triplet T2. For example, normalized
doublet components are

    O_E1=(u_0u_1-u_0u_2-u_1u_3+u_2u_3)/2,
    O_E2=(u_0u_1+u_0u_2-2u_0u_3-2u_1u_2+u_1u_3+u_2u_3)/sqrt(12).

All are bounded, even, gauge invariant operators of the existing theory.
On a five-vertex four-arm star the complete unregulated physical block has
126 states: leaf Gauss law fixes E_leaf=n_leaf in {0,1}, and central fields
are hole-n_leaf in {-1,0,1}. Cutoff one therefore already contains every
physical electric state of this particular graph. The **actual** rewiring
diagonal is 0.06 per removable edge.

Exact diagonalization finds nonzero doublet and triplet response. The first
E excitation has gap 3.0522905057 and degeneracy two; the first T2 excitation
has gap 2.5505970919 and degeneracy three. All 24 simultaneous leaf, rotor,
and fermion permutations, including Fock signs, verify these representations.
Kubo evolution, spectral sums, and independent finite pulses agree.

This establishes a response channel missed by mean bond variables. It does
not establish that these operators are an operational metric, that their
excitations are massless, or that they couple to the full matter stress.

## 7. Full geometry mixing and the event-pair representation

The complete N=5,Q=4 calculation has 37,965 states at electric cutoff two
and 49,665 at cutoff three. The ground gap is 0.0684711204. Both E and T2
probes couple to its fivefold shell, with very different overlaps. Complete
static susceptibilities are obtained from a projected linear solve, not
from truncating a few low eigenvectors.

The fivefold shell has event-permutation character

    chi_pair(p) = binom(fixed_vertices(p),2)
                  + number_of_two_cycles(p) - fixed_vertices(p).

This is the unordered-pair permutation space with its singlet and standard
vertex representation removed: the [N-2,2] representation. Its dimension is
N(N-3)/2, not five at every N. Restriction to the four-neighbor tetrahedral
group gives the doublet plus triplet at N=5. This explains the apparent
tensor degeneracy without establishing a graviton.

For larger N the exact signed quotient under S_2 x S_(N-2) keeps an unordered
marked pair. Its basis is built from full physical-state orbits. Negative
stabilizer signs give zero invariant orbit; allowed orbits carry their exact
Fock/incidence phases. The quotient does not retain only trees, a charge
pair, or a fixed geometry. It retains all physical states in this declared
symmetry subspace at the numerical electric regulator.

The central transposition sum commutes with H and has only three allowed
eigenvalues in this subspace:

    trivial: N(N-1)/2;
    standard: N(N-3)/2;
    pair: (N-1)(N-4)/2.

These follow by viewing a transposition on unordered vertex pairs as a
Johnson J(N,2) adjacency plus its diagonal. Frobenius reciprocity gives
one subgroup-invariant vector per copy of these three representations.
Integer similarity matrices verify H/class-sum commutation and the cubic
minimal polynomial without floating roundoff at moderate sizes. Large
calculations use recorded deterministic-column checks of those identities;
the exact group construction is the same.

The independent frozen builder verifies the intertwiner H T=T H_quotient
for complete N=4 and N=5 matrices, including the actual rewiring diagonal.
Larger eigenvectors are therefore genuine eigenstates in the declared
symmetry sector. The lowest singlet is the complete permutation-invariant
ground state. A global ground in another irrep is not ruled out at sizes
where the full spectrum is not searched.

Numerical tables and electric-cutoff comparisons are in REPORT.md and the
JSON results. If relabelling is gauge, these nontrivial pair modes alone
are not physical states; a reference system must carry the compensating
representation. If labels are physical, they are legitimate model modes,
but still have no demonstrated wave vector or Lorentz helicity.

## 8. The existing Gauss law is not the gravitational energy constraint

On the complete star, compare the neutral sea with a central fermion and
an edge hole linked by E_center=1. All exterior leaf occupations and leaf
electric fields are identical (zero) in the two states. Both satisfy every
Gauss constraint. Their energy expectations differ by

    2m+epsilon_E/2-a = 2.065,     a=0.06.

An arbitrary operator supported only on that unchanged exterior has the
same expectation in both states. The star Hamiltonian cannot therefore be
identified with such a boundary energy operator on its whole physical
Hilbert space. This is a direct check that gauge-charge Gauss law does not
already implement a Gauss law for energy.

A different low-energy code, gravitational dressing, or nonlocal metric
identification could change that conclusion. Marolf's boundary-Hamiltonian
argument [1] assumes a compatible locality map and the same time evolution;
it is not a universal exclusion of every CP-12 encoding. Donnelly and
Giddings [2] establish why energy-carrying observables in perturbative
gravity require long-range metric dressing. Merely having nonlocal graph
coherence is not a derivation of that dressing.

## 9. A sufficient target for the next genuine gravitational test

An independently specified tensor operator family should have, in a
controlled physical state and continuum regime, a retarded TT pole

    chi_TT(omega,k) ~ Z P_TT(k)/[c^2 k^2-(omega+i0)^2],

with two physical polarizations, not just a finite point-group degeneracy.
The six h_ij and six canonical momenta require four first-class constraints
to leave a four-dimensional physical phase space. The conventional linear
constraints are k_j pi_ij=0 and
(k_i k_j-k^2 delta_ij)h_ij=0. Their Poisson brackets vanish; the TT projector
has rank two. `spin2_consistency.py` checks this benchmark explicitly.

With Lorentz-invariant soft-pole factorization, longitudinal graviton
decoupling gives sum_a eta_a g_a p_a^mu=0. Generic scattering together
with momentum conservation requires universal g_a. The script independently
checks this against unequal-coupling counterexamples. This is Weinberg's
established result [3], not a CP-12 Ward identity. Consistent self-coupling
of a massless spin-two field supplies the Einstein nonlinear structure
under the usual locality/derivative assumptions [4,5]. These conditions
cannot be invoked before the microscopic model earns the field and its
constraints. Higher derivative EFT terms and microscopic selection remain.

Existing lattice tensor models also distinguish helicity-two structure from
Einstein gravity: Gu and Wen's controlled L-type model has cubic dispersion,
while their linear-dispersion N-type claim lacks the same controlled
approximation [6]. This is a relevant precedent, not a completion to import
into CP-12.

The next test should therefore hold the frozen law fixed, specify a
relational reference/state independently, measure momentum-resolved
mixed tensor–matter response, and verify a finite pole with the necessary
constraints and energy dressing. A connectivity statistic, a flat neutral
band, or an event-pair multiplet cannot substitute for it.

## Primary literature

1. D. Marolf, *Emergent Gravity requires (kinematic) non-locality*,
   [arXiv:1409.2509](https://arxiv.org/abs/1409.2509), especially pp. 3–6.
2. W. Donnelly and S. B. Giddings, *Observables, gravitational dressing, and
   obstructions to locality and subsystems*,
   [arXiv:1607.01025](https://arxiv.org/abs/1607.01025), section III.
3. S. Weinberg, *Photons and Gravitons in S-Matrix Theory*,
   [Phys. Rev. 135, B1049 (1964)](https://doi.org/10.1103/PhysRev.135.B1049).
4. S. Deser, *Self-Interaction and Gauge Invariance*,
   [arXiv:gr-qc/0411023](https://arxiv.org/abs/gr-qc/0411023), original 1970 paper.
5. N. Boulanger, T. Damour, L. Gualtieri and M. Henneaux,
   *Inconsistency of interacting, multi-graviton theories*,
   [arXiv:hep-th/0007220](https://arxiv.org/abs/hep-th/0007220), theorem 1.1.
6. Z.-C. Gu and X.-G. Wen, *Emergence of helicity +/-2 modes from qubit models*,
   [arXiv:0907.1203](https://arxiv.org/abs/0907.1203), sections II, III and conclusions.

The Johnson/inclusion identities and permutation representations are
established finite combinatorial structures. Their frozen CP-12 applications
are reproduced here; no priority claim for the underlying mathematics is made.
