# CP-12 independent tensor-symmetry audit: event-pair degeneracy is not helicity or momentum

**2026-10-13.** Separate independent additive work rooted at [11fbaec7](https://github.com/ryangrvr/TestingGRUT/commit/11fbaec7b28a1bbe88e6ec54808829ebbc804c48), which is two commits after [d2fbd05](https://github.com/ryangrvr/TestingGRUT/commit/d2fbd05c420a209faae618c88fbe4d1902bc999c). The frozen CP-12 law, previous results, canonical GRUT-RAI, Stage-3 Card 3 and other branches are **not changed**. This is a symmetry/theory audit and a *separately declared uniform-graph comparator*, NOT an independent re-eigensolve of all CP-12 states or a discovery of a new physical law. Permutation representation theory and spectral graph facts used here are established mathematics.

## 1. Two successful owner results concern DIFFERENT physical sectors

The full-Hamiltonian variational theorem of d2fbd05, at **Q=2N**, excludes a ground-state phase in which absolute defects from degree four are subextensive:
\[
 \Pr_{\rm ground}[\sum_i|d_i-4|\le K_N]\to0\quad(K_N=o(N)).
\]
It includes graph mixing, all integer electric fields and the actual rewiring diagonal. The trial-state energy advantage over a fixed-graph block lower bound is \(\simeq0.056315N\). **A spatial phase with a nonzero density of defects remains entirely open.** This reviewer did not independently reproduce the 80-digit certificates and does not upgrade them to empirical validation.

The gravity-response report 11fbaec, at **Q=N-1** for its \(N=5,6,7\) event-pair shells, finds lowest measured pair multiplet gaps 0.06847112, 0.06938808, 0.07012838 and dimensions **5, 9, 14**. It explicitly warns that the dimensions follow the event-permutation representation \([N-2,2]\), not physical Lorentz spin. The seven-event result uses a signed quotient at one electric cutoff, not a complete full-space eigensolve or checked cutoff convergence at N=7.

These **different Q sectors** cannot be sewn into a single demonstrated gravitational phase. In particular the pair gaps **increase** slightly over N=5–7, not providing evidence of a massless mode; that tiny sample does not prove a positive thermodynamic gap either.

## 2. Explicit projector onto the measured pair representation

Let \(\mathcal E=\{\{i,j\}:i<j\}\) contain \(M=N(N-1)/2\) unordered event-pair labels, N>=4. Define the **unoriented** incidence matrix \(B\in\mathbb R^{N\times M}\), \(B_{i,e}=1\) when event i belongs to pair e. It obeys
\[
 BB^T=(N-2)I_N+\mathbf1_N\mathbf1_N^T
\]
and has rank N for N>=3. Hence
\[
 \boxed{P_{[N-2,2]}=I_M-B^T(BB^T)^{-1}B}
\]
is an **exact orthogonal projector** of rank
\[
 M-N=N(N-3)/2.
\]
Every \(h\) in its image obeys \(\sum_{j\ne i}h_{ij}=0\), one constraint per event. The full pair-label space decomposes as
\[
 \mathbb R^M=\mathbf1\oplus\mathbf{(N-1)}\oplus\mathbf{[N-2,2]}.
\]
The row-sum constraint is a **permutation-incidence constraint** and must not be identified with a continuum constraint \(k_i h_{ij}(k)=0\) for a **physical** spatial vector k.

As an independent character check, for any permutation p with f fixed events and c2 two-cycles,
\[
 {\rm tr}(P_{[N-2,2]}R(p))=\binom f2+c_2-f.
\]
This exactly matches the character formula reported by the owner. [check.py](check.py) verifies projection rank, idempotence, incidence transversality and this character for random permutations N=4..7, plus N=10,16 projection checks.

## 3. Three-orbit response theorem for an unbroken event-permutation state

**THEOREM.** Let \(U(p)\) represent the full event permutation group, assume \([H,U(p)]=0\) and \(U(p)\rho U(p)^\dagger=\rho\). Let \(A_e,B_e\) be any uniformly covariant families of bounded operators labelled by unordered event pairs e. Then the retarded response
\[
 K_{ef}(T)=i\,{\rm Tr}\bigl(\rho[A_e(T),B_f]\bigr)
\]
can depend on **only three orbits** of the ordered pair (e,f):
(i) e=f, (ii) e,f share exactly one event, or (iii) they are disjoint. This holds at **every T**, irrespective of the many-body interactions, degenerate energies, electric rotors or graph coherence. For vertex-indexed covariant operators, there are only **two** orbits: same vertex versus different vertices.

**Proof:** simultaneous conjugation by U(p) maps \(K_{ef}\) to \(K_{p(e),p(f)}\). The event-permutation group acts transitively on each intersection class and has no other orbits on two unordered pairs. Therefore the coefficients are identical within each class. QED.

For an edge-pair-indexed kernel with class coefficients a(T),b(T),c(T), write
\[
 K(T)=aI+bA_{\rm share}+cA_{\rm disjoint},
\]
where \(A_{\rm share}\) is the adjacency of the Johnson graph J(N,2) on pair labels, and \(A_{\rm disjoint}=\mathbf1\mathbf1^T-I-A_{\rm share}\). In the three irreducible components the kernel has eigenvalues
\[
 \begin{aligned}
 \lambda_{\bf1}&=a+2(N-2)b+\tbinom{N-2}{2}c,\\
 \lambda_{\rm standard}&=a+(N-4)b-(N-3)c,\\
 \boxed{\lambda_{[N-2,2]}&=a-2b+c.}
 \end{aligned}
\]
Thus **one copy of the pair representation carries one identical kernel eigenvalue**, with degeneracy N(N-3)/2. Schur's lemma gives the same result, but the Johnson-adjacency construction is an explicit finite-N proof. The spectral functions a,b,c can have arbitrarily rich frequency dependence, so this is **not** a theorem that the full many-body spectrum is flat or gapped.

**Crucial scope:** the full CP-12 Hilbert space can contain **multiple copies** of \([N-2,2]\). An S_N-symmetric full Hamiltonian can act nontrivially among those copies while remaining scalar in each *irrep index*. The theorem only forbids extracting **spatial momentum labels or gravitational helicity** directly from the individual irrep coordinate of a single pair-label space. A physically selected, relational graph-sector or reference-carrying probe can provide additional structure and escape this labelled-orbit limitation.

## 4. Independent uniform fixed-Q annealed graph control

To prevent a common interpretive mistake—averaging graph geometry away before seeking its wave equation—consider a **separate** uniform classical ensemble of simple graphs with exactly Q edges. Each unordered edge has occupation mean \(p=Q/M\). Its **averaged vertex Laplacian** is exactly
\[
 \boxed{\langle L_G\rangle_{\rm uniform}
  =p(NI-\mathbf1\mathbf1^T)
  =\frac{2Q}{N-1}\left(I-\frac{\mathbf1\mathbf1^T}{N}\right).}
\]
It has one zero mode and **N−1 identical positive modes**, each of gap \(2Q/(N-1)\). For Q=N−1 this gap is identically **2** at *every* N. There is no momentum-dispersive acoustic band in this **average**.

Contrast a particular selected, unit-weight N-cycle, whose vertex Laplacian has smallest positive eigenvalue
\[
 \boxed{\lambda_{\rm ring}=4\sin^2(\pi/N)\sim4\pi^2/N^2.}
\]
These different operators give different long-wavelength physics. The uniform graph ensemble itself typically has a giant random-graph component, as the owner's prior preparation report correctly explained. **Neither average is the actual dynamical CP-12 quantum response**, and in general \(\exp(-it\langle H\rangle)\ne\langle \exp(-itH)\rangle\). The comparator is an exact warning against substituting a permutation-invariant **annealed kernel** for a physically conditioned geometry.

| N | Rank of pair irrep | Uniform mean-L gap, Q=N−1 | Prepared N-cycle gap |
|---:|---:|---:|---:|
| 4 | 2 | 2.000000 | 2.000000 |
| 5 | 5 | 2.000000 | 1.381966 |
| 6 | 9 | 2.000000 | 1.000000 |
| 7 | 14 | 2.000000 | 0.753020 |
| 10 | 35 | 2.000000 | 0.381966 |
| 16 | 104 | 2.000000 | 0.152241 |

The coincidence at N=4 is not a continuum convergence. The uniform mean-L has no linear or quadratic spatial dispersion for any N. The ring's N^-2 gap results from the **specific cyclic neighborhood**, which was not selected by the full CP-12 ground state.

## 5. Explicit contrast with a *physical* transverse-traceless projector

Given a **supplied three-dimensional Euclidean metric and nonzero physical wavevector** \(k\), define \(P_{ij}=\delta_{ij}-\hat k_i\hat k_j\). The conventional symmetric spatial rank-2 TT projector is
\[
 \Pi^{TT}_{ij,kl}=\tfrac12(P_{ik}P_{jl}+P_{il}P_{jk})
                   -\tfrac12P_{ij}P_{kl},
\]
and \({\rm rank}\,\Pi^{TT}=2\). [check.py](check.py) independently checks its rank, Hermiticity and idempotence for \(k\parallel z\).

**The two constraints differ physically:** the pair graph projector sums edge labels incident at a vertex, while \(\Pi^{TT}\) depends on *physical momentum* and a spatial metric. The dimensions 2,5,9,14 do not count helicities: even if many discrete spatial momenta existed, it would take a separately earned geometry and dispersion relation to associate them with TT modes. No such association is established by a pair-shell degeneracy.

## 6. Next test that could actually reveal new gravitational physics

Do not declare success from a low multiplet gap, a tetrahedral doublet, or an S_N tensor-looking label. Use frozen CP-12 **at the same Q/N density and in the same physical state**:

1. Select a state, not just an energy sector: a prepared geometry, a genuine ground state, or a mixture, with the state preparation and reference energy tracked.
2. Define physical relational source and detectors by actual **charge/flux records** and bond traversal, not by globally supplied event labels. At minimum record graph-conditional shells \(d_G(i,j)=r\), normalized shell occupancy and the reference back-reaction. Ideal graph-controlled shells used in calculation should be distinguished from a constructed autonomous local detector.
3. Characterize volume growth \(V(r)\), heat-kernel/spectral dimension and graph-conditioned matter response on the **same ensemble**. Avoid calling logarithmic-distance random graph growth an emergent 3D manifold.
4. Only if a spatial continuum regime emerges, build geometry-conditioned long-wavelength modes \(\phi_k(G)\) from its own transport operator and test an independently normalized retarded response with
   \[
   \chi_{TT}(\omega,k)\stackrel{?}{\sim}
   \frac{Z(k)P_{TT}(k)}{c^2 k^2-(\omega+i0)^2},
   \]
   including finite nonzero residue, exactly two independent physical tensor channels, and the appropriate energy/stress constraints.
5. Compare mixed tensor response to **at least two inequivalent matter species** to probe universal coupling; a single fermion–geometry response cannot establish equivalence. Keep scalar/vector counterchannels and a label-shuffled null control.

A nontrivial relational momentum profile would be a new physical achievement for CP-12. If it remains only a permutation-irrep eigenvalue, the tensor-like modes are not identified as gravitational. If an irregular but manifoldlike ground phase develops, it is allowed by the Q=2N no-regularity theorem. Both directions must be evaluated with one frozen Hamiltonian, not spliced from \(Q=N-1\) and \(Q=2N\) calculations.

## Decision

The new Q=2N variational theorem **rules out an overly clean regular ground phase** at those frozen couplings. The new Q=N−1 pair-shell analysis **isolates genuine many-body symmetry responses** but does not show a graviton.

Independent exact permutation representation tests establish why the shell multiplicities alone do not certify physical tensor helicity or spatial momentum. The next calculation must produce a **geometry-resolved, graph-conditioned quantum response**, not another unsituated irrep degeneracy.

No changes to the frozen law, no CP-13, no Card 3, and no publication-level novelty claim.
