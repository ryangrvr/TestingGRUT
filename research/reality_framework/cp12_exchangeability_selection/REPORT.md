# CP-12 symmetry gate: internal records do not imply labelled spatial selection

**14 October 2026.** Analytic continuation of 8ece3e2; frozen Hamiltonian unchanged. No new numerical time evolution claimed.

## Exact finite-N symmetry theorem

Let U_pi be the *physical* unitary representation of vertex permutations, including all required fermionic signs and gauge/rotor relabellings. The frozen labelled CP-12 Hamiltonian is permutation covariant, [H,U_pi]=0, and the equal-amplitude neutral preparation used in 8ece3e2 is permutation invariant at the level of its density operator (its ray may transform by a one-dimensional phase). Therefore
\[
 U_\pi \rho_T U_\pi^\dagger=\rho_T
\]
for every permutation pi and every time T.

For any fixed potential edge e, let n_e^G denote its graph-presence projector. The action of S_N on edges is transitive and \(\sum_e n_e^G=Q I\) on the fixed-Q sector. Therefore
\[
\boxed{\langle n_e^G\rangle_T=\frac{Q}{M}\quad\text{for every e, every T}.}
\]
Likewise, let p_e denote the projector that edge e is present and nonremovable (occupied-edge flux or hole). Then \(\sum_e p_e=Q-r\), so
\[
\boxed{\langle p_e\rangle_T=\frac{Q-\langle r\rangle_T}{M}\quad\text{for every e, every T}.}
\]
For \(Q=\nu N+O(1)\), these are \(2\nu/N+O(N^{-2})\) and \(2(\nu-\langle r\rangle/N)/N+O(N^{-2})\), respectively. Even substantial total protection cannot assign a preferred *labelled* edge in the unconditional state. The one-vertex matter density is similarly uniform: \(\langle n_{a_i}\rangle=\langle n_v\rangle/N\).

The theorem does **not** forbid nontrivial pair, triangle, cycle, graph-distance, or relational correlations; the full graph distribution may change substantially while all edge marginals remain equal. It also does **not** forbid a symmetric quantum superposition of individually spatial graphs, spontaneous symmetry breaking in an infinite-volume limit, or spatial order visible after conditioning on internal records. An emergent relational geometry need not choose external vertex labels.

## Exact relational diagnostic / next scientific gate

The correct comparison is **not** whether an unconditional edge marginal singles out an edge; symmetry forbids it by construction. Instead measure the *joint graph–environment correlations* in the same evolved state, without applying a global measurement during the dynamics:

1. For complete matter/flux configuration mu, its protected edge set P(mu) is a kinematic internal record. Characterize its induced graph **up to relabelling**, e.g. short-cycle packing, degree distribution, intrinsic ball growth, and connected-component structure. Because the environment may be in superposition, treat this as a diagnostic classical joint distribution over final complete configurations, not an autonomous observer selecting an outcome.
2. For full graph g in the same joint distribution, measure the fraction of protected edges, conditional graph distances and the number of removable destinations whose addition would shorten long relational distances. Do not conflate a *candidate* edge pool chosen after observing a graph with a fixed pool in the earlier norm theorem.
3. Compare size-scaled distribution of protected-core and full-graph observables with the uniform fixed-Q labelled-graph null and an energy-matched graph/flux null, not just with an imposed lattice.
4. Evaluate physical source–receiver susceptibility **conditioned on a locally accessible record** (if one can be constructed); account for recording work and backaction. Do not treat postselected geometric subensembles as freely prepared spatial phases.
5. Seek a concentration statement for *permutation-invariant spatial observables* (e.g. probability of extensive disjoint bounded loops and polynomial intrinsic ball growth), rather than expecting one labelled graph to dominate.

These requirements expose a missing bridge between demonstrated flux/entanglement and spatial order. No spontaneous spatial selection, 3D growth, causal cone, or state-selection mechanism is derived by the symmetry theorem.

## Important control on claims

The entropy 6.344 bits reported in 8ece3e2 is entanglement of the graph factor with the matter/flux factor in the declared kinematic embedding. It is not an observer record, measurement outcome, a classical selection of a particular graph, or a loss of global purity. The exact symmetry theorem applies equally when that entropy is large.

A source n_0 explicitly distinguishes vertex 0 for a diagnostic response and is not an unforced symmetry-breaking mechanism in the source-free state. Conditioning on a source, a chosen receiver, or a measurement record also breaks the exchangeability of that *conditional* description. This must not be mistaken for spontaneous selection.
