# Coherent destinations, records, and higher-order correlations

Frozen law: [CP-12 at 11fbaec](https://github.com/ryangrvr/TestingGRUT/blob/11fbaec7b28a1bbe88e6ec54808829ebbc804c48/research/constructive_physics/attempt12/LAW.md). The host-only first-rewire theorem at [8cd95db](https://github.com/ryangrvr/TestingGRUT/blob/8cd95db20f18764eed12f1868f2da54b13131a3f/research/reality_framework/cp12_destination_uniformity_audit/REPORT.md) is correct. Its scope can be extended without changing the Hamiltonian.

## 1. Common-vacant destinations remain uniform under graph coherence

Let \(\rho\) be any physical initial state in the fixed-\(Q\) sector, with \(n_f\rho=0\) for every edge in a candidate set \(\mathcal V\). This allows arbitrary graph, fermion, electric, and auxiliary-reference coherence. Every candidate edge is absent in every initial branch; CP-12 then requires its edge fermion filled and both electric legs zero.

Define a mathematical deletion map \(D=\sum_e D_e\): \(D_e\) removes edge \(e\) if it is neutral/removable, otherwise giving zero. This maps to a formal \(Q-1\) space for bookkeeping; the Hamiltonian itself still conserves \(Q\). Define \(C_f\) to restore a neutral candidate \(f\). On the range of \(D\rho\), \(C_f\) is an isometry, independently of \(f\in\mathcal V\).

Let \(V\) be the bounded graph-changing part of the Hamiltonian. Only it can create \(f\), so

\[
P_f V\rho^{1/2}=-\frac\Gamma M C_fD\rho^{1/2}.
\]

The graph-preserving block evolution is strongly continuous and \(V\) is bounded. Its Duhamel integral gives the following limit at fixed finite \(N\), even for normal states that lack moments of the unbounded electric energy:

\[
\boxed{\lim_{T\to0}\frac{\Pr(n_f=1;T)}{T^2}=\frac{\Gamma^2}{M^2}
\operatorname{Tr}(D\rho D^\dagger),\qquad f\in\mathcal V.}
\]

For the finite-support states used here, or suitably regular states with the required Hamiltonian moments, the expansion has an \(O_N(T^3)\) remainder. For a real such state, Hamiltonian, and graph readout, evenness in time improves it to \(O_N(T^4)\). A quantitative remainder must not be asserted for arbitrary infinite-energy physical states. For a definite graph, the deletion images of distinct removed edges are orthogonal and the coefficient reduces to the expected removable-edge number. For different initial graphs, common deletion images interfere. They change the coefficient equally for **all** common-vacant destinations.

The argument also covers a joint reference state: include its register in \(\rho\) and act with \(D\otimes I\). Partial or complete graph records change the deletion norm but do not change this destination equality. If a complete initial-graph record is retained, each individual record branch recovers the original definite-graph theorem for all of that branch's absent edges.

This is not a theorem about arbitrary destinations already present in some branches. A superposition over all graphs has no universally absent edge. Nor is it an all-time theorem for full CP-12.

## 2. Coherence and record visibility in a sparse pool

For a supplied pool \(F\) of \(S\) possible edges, prepare the neutral equal superposition of all its \(Q\)-edge subsets. The unsigned deletion incidence matrix satisfies \(D^\dagger D=QI+A_{J(S,Q)}\), giving

\[
\|D\psi\|^2=Q(S-Q+1).
\]

If graph records have a common real off-diagonal overlap \(\eta\in[0,1]\), the reduced state is the mixture of the coherent state and its graph-dephased version. Then

\[
\|D\rho_\eta^{1/2}\|^2=Q[1+\eta(S-Q)],
\qquad E_\eta=\kappa Q+\frac\Gamma M[Q(M-Q)-\eta Q(S-Q)].
\]

Thus coherence can accelerate common-vacant creation while lowering the initial rewiring energy. It does not make the creation destinations preferentially local. The record energy shift concerns the CP-12 system; a physical apparatus and total energy exchange are not modeled here.

For \(F=C_N\), \(Q=N-1\), every initial graph is a spanning path. For any \(2\le R<N/2\), the fraction of first creations outside the pool that lie within distance \(R\) in the **supplied ring reference** is

\[
\frac{N(R-1)}{M-N}=\frac{2(R-1)}{N-3}\longrightarrow0.
\]

This survives arbitrary record visibility, and a charged bond held in every branch also leaves the ratio unchanged. `coherent_ring.py` checks the integer deletion counts, including that locked-bond control. A ring reference and the initial path coherence are supplied resources, not assembled spacetime.

## 3. Full matter dynamics breaks the uniformity later

There is no corresponding general equality at every order. An independent unregulated integer Hamiltonian builder gives an explicit counterexample. At \(N=4,Q=3\), prepare the neutral state proportional to the sum of the graphs \(\{01,12,23\}\) and \(\{02,12,23\}\). Edges 03 and 13 are absent in both. Their occupation probabilities agree through \(T^4\), but

\[
P_{13}(T)-P_{03}(T)=\frac{T^6}{72000}+O(T^8).
\]

With matter hopping disabled as a declared diagnostic, all common-vacant probabilities are equal at every tested order. In fact the rewiring-only operator has the exact symmetry exchanging the vacant edge slots, so that comparator equality holds at all times.

The primary distance control uses one sector \(N=5,Q=4\), with supplied pool \(F=\{01,12,23,03,34\}\). Its initial state is the neutral equal superposition of the five graphs obtained by deleting one pool edge. The candidates 04 and 14 have equal endpoint-degree distributions in this preparation, but reference-pool distances two and three. Exact full-law arithmetic gives

\[
\boxed{P_{04}(T)-P_{14}(T)=-\frac{1103}{3240000000}T^{10}+O(T^{12}).}
\]

All coefficients through \(T^8\) cancel. The longer destination is initially favored in this particular comparison. Under complete initial graph dephasing the coefficient is half as large, and partial record visibility multiplies it by \((1+\eta)/2\). The leading normalized near/far difference is therefore independent of that visibility: coherence is not an independent leading locality filter here.

The exact sixth-order differences among the five destinations in this control track endpoint-degree sums. The two equal-degree destinations were chosen to separate that effect from the later difference. The statement about degree sums is an executed finite-case identity, not a universal claim for arbitrary graph preparations.

Full finite evolution confirms unequal destinations and a later reversal of the near/far sign. No monotone distance law follows. The references are properties of \(F\), not operational distances derived from the evolving state; no finite-speed signal is inferred.

## 4. Unordered initial coherence changes global shape without changing edge means

Prepare the neutral pure state with equal amplitude over **every** labelled \(Q=N-1\) graph. No spatial pool, loop reward, or ordering is supplied in its amplitudes. Its full-H energy is \(\kappa Q\), and its variance is \(2t^2Q\). Its density matrix is invariant under simultaneous event permutations, so each labelled edge has mean \(Q/M\) at every time.

Global permutation-invariant shape observables need not be stationary. The matter term alone preserves graph occupation, while the neutral rewiring term alone leaves this initial state stationary. Their coupled dynamics changes the graph distribution.

The exact signed invariant-sector calculation retains all physical configurations at its cutoff. It is a restriction to the symmetry of the chosen state, not a claim that relabelling has been derived as a gauge redundancy. A global neutral-background character is included so the prepared orbit sums have a consistent phase. The density and all measured shape probabilities do not depend on that overall character convention.

An exact integer/rational expansion at \(N=4,5,6,7\) finds

\[
P_G(T)=\frac1{\binom MQ}\left[1+c_N(\overline W-W_G)T^6\right]+O_N(T^8),
\qquad W_G=\sum_i\binom{d_i}{2},
\]

for every graph isomorphism class, with the probability per labelled member understood in this expression. The second- and fourth-order probability changes vanish. The verified coefficients are

| N | Exact coefficient \(c_N\) |
| --- | --- |
| 4 | \(5/1296\) |
| 5 | \(13/3125\) |
| 6 | \(8729/2025000\) |
| 7 | \(349/79380\) |

All four match \(\ell_N(2\Delta-\ell_N)/90\), where \(\ell_N=\Gamma(M-1)/M\) and \(\Delta=2m+\epsilon_E/2\), at the frozen couplings. That closed expression is recorded as a candidate general identity; the script proves its finite cases, not a uniform-in-\(N\) expansion theorem. These are short-time coefficients of the full frozen law, not the weak-hopping ground-energy expansion.

The positive coefficient biases against degree concentration. In particular the leading change of the mean wedge count is \(-c_N\operatorname{Var}_{\rm uniform}(W)T^6\). This distinguishes a genuine dynamically generated correlation from a geometric preference inserted into the action. It does not distinguish all spatial graphs from nonspatial graphs sharing their degree statistics.

The independent full-builder intertwiners at four and five events agree with the signed compression, including time evolution. Larger calculations retain only the invariant state sector but include all graph, matter, and electric configurations within it. The seven-event cutoff-five space has 487,273,213 full physical states and 139,097 invariant states; only the latter Hamiltonian was evolved. It is not a full-spectrum solve.

## 5. Relation to the entropy result

The finite geometry change is real, but its pure coherent preparation has zero entropy. `GIBBS_AND_PREPARATION.md` derives why low entropy matters, and why a high-entropy neutral mixture cannot overcome narrow graph-family rarity in \(\nu>1\) through isolated unitary evolution.

The formation controls have \(Q=N-1\). The thermal rarity theorem's leading \(N\log N\) gap requires \(\nu>1\). No tensor pole, Lorentzian continuum, universal coupling, or spatial thermodynamic phase is established by combining these separate statements.
