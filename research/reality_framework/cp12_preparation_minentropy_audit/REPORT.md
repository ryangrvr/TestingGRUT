# CP-12: an all-time min-entropy strengthening of the loop-rich preparation obstruction

**2026-10-14 — independent analytical extension of Work commit 63bbd4af; partner audit recommended.** Additive branch only. Frozen CP-12, canonical GRUT-RAI, and previous results untouched. All statements concern the complete **labelled** physical Hilbert space at fixed Q=nu N+O(1), and the event of at least alpha N vertex-disjoint cycles of microscopic length <=R, where R, alpha, nu are fixed positive constants. Does not claim emergent spacetime, geometry selection or experimental gravity.

## Main result: stronger than logarithmic suppression for the explicit incoherent mixture

Let M=binom(N,2), D=binom(M,Q), and |G,0> denote the physical neutral sea on labelled Q-edge graph G. Prepare

\[
\rho_0=\frac1D\sum_{G:|E(G)|=Q}|G,0\rangle\langle G,0|
\]
and evolve with the **complete frozen CP-12** Hamiltonian H (fermion hopping, positive rewiring Laplacian including diagonal, U(1) rotors with unbounded integer fields), i.e. rho_t=e^{-itH}rho_0e^{itH}. Write Pi_S for the graph-diagonal projector onto S={G:K_R(G)>=ceil(alpha N)}.

**Theorem:** For every fixed exponent 0<delta<alpha/nu,
\[
\boxed{\sup_{t\in\mathbb R}\operatorname{Tr}(\rho_t\Pi_S)
\le (2|t_{\rm hop}|+\Gamma)\,Q\,N^{-(1+\delta)}
+\exp[-(\alpha-\nu\delta)N\log N+O(N)].}
\]
Hence
\[
\boxed{\sup_t P_t(S)=O(N^{-\delta})\quad\text{for every }\delta<\alpha/\nu.}
\]

This is a strict strengthening of the parent paper's general \(O(1/\log N)\) result for **this explicit high-entropy preparation**. It does not improve the full class \(D(\rho_0\|\sigma_\beta)=O(N)\) without an additional min-entropy condition: a small eigenvalue spike may be consistent with small relative entropy.

For conventional prepared square and cubic tori with alpha=1/4, representative valid choices are delta=1/16 for Q=2N, and delta=1/24 for Q=3N (neither is optimal). Their finite-size bounds are extremely loose; this is an asymptotic theorem, not a laboratory estimate.

## Independent proof (no relative entropy is needed)

1. **Conserved maximum eigenvalue.** rho_0 is a normalized rank-D projection, so
\[
\|\rho_t\|_\infty=D^{-1}\quad\text{for ALL times}.
\]
Closed unitary evolution preserves the entire spectrum; this is a min-entropy fact and needs no thermal state.

2. **Conserved extensive energy controls electric fields.** Let
\[
F=\frac{\epsilon_E}{2}\sum_{e,i\in e}E_{ie}^2\ge0,\quad H=\kappa Q+2mN_v+F+T_{\rm hop}+R.
\]
Because m>=0, R>=0 and T_hop>=-2|t_hop|Q,
\[
F\le H-\kappa Q+2|t_{\rm hop}|Q
\]
as a quadratic-form inequality. On each neutral sea, the exact **full H diagonal** is \(\kappa Q+\Gamma Q(M-Q)/M\); thus
\[
\langle F\rangle_t\le 2|t_{\rm hop}|Q+\Gamma Q(M-Q)/M\le(2|t_{\rm hop}|+\Gamma)Q.
\]
Let \(P_>=1_{\{F>N^{1+\delta}\}}\). Markov's inequality on the positive observable F gives, for every t,
\[
\operatorname{Tr}(\rho_tP_>)\le (2|t_{\rm hop}|+\Gamma)Q/N^{1+\delta}.
\]

3. **Few low-electric graph states are loop-rich.** The graph event Pi_S commutes with the electric-energy cutoff \(P_\le=1-P_>\), although neither generally commutes with H. On a fixed Q-edge graph, the number of physical fermion configurations is at most binom(N+Q,Q). There are 2Q electric half-legs. For any s>0, the number of integer vectors obeying F<=N^{1+delta} satisfies
\[
\#\le \exp[(2s/\epsilon_E)N^{1+\delta}]\Theta(s)^{2Q},
\qquad\Theta(s)=\sum_{z\in\mathbb Z}e^{-sz^2}.
\]
Choose \(s=N^{-\delta}\). Using \(\Theta(s)\le1+\sqrt{\pi/s}\), this gives
\[
\log\operatorname{rank}(\Pi_S P_\le)
\le\log|S|+\log\binom{N+Q}{Q}+\delta Q\log N+O(N).
\]
The exact parent counting bound gives
\[
|S|/D\le A_{N,Q,R}^{\,\lceil\alpha N\rceil}/\lceil\alpha N\rceil!
=\exp[-\alpha N\log N+O(N)].
\]
Consequently,
\[
\operatorname{rank}(\Pi_S P_\le)/D
\le\exp[-(\alpha-\nu\delta)N\log N+O(N)].
\]

4. **Split the final probability.** Since Pi_S and P_\le commute,
\[
\operatorname{Tr}(\rho_t\Pi_S)
=\operatorname{Tr}(\rho_t\Pi_S P_\le)+\operatorname{Tr}(\rho_t\Pi_S P_>)
\le\frac{\operatorname{rank}(\Pi_S P_\le)}{D}+
\operatorname{Tr}(\rho_tP_>).
\]
The result follows. It holds **uniformly for arbitrarily late t**, without assuming the graph sector or electric cutoff is invariant, and without modifying the Hamiltonian.

## General min-entropy/energy variant

For any initial density operator satisfying (i)
\[
\|\rho_0\|_\infty\le e^{cN}/\binom M Q
\]
with fixed c, and (ii) \(\langle H-\kappa Q\rangle_0=O(N)\), the same derivation gives \(\sup_t P_t(S)=O(N^{-\delta})\) for every delta<alpha/nu. This is a **sufficient condition**, not a conclusion for all O(N)-relative-entropy preparations. In particular, a single pure graph state or a pure uniform coherent graph state does not satisfy (i).

The result translates the operational resource requirement into a bound on initial *min-entropy* (largest eigenvalue) as well as energy, without claiming either can produce spacetime.

## Numerical/independent controls and scope

The standalone checker reproduces explicit N=6,Q=6 short-cycle graph fractions by exhaustive labelled enumeration; it tests the finite-N union bound, electric theta estimate and limiting exponent at Q/N=2,3. Small-N inequalities are generally vacuous, and even the refined polynomial estimate may become nonvacuous only at astronomically large N because all states and flux configurations were generously overcounted.

Other essential limits:

- The theorem is an **all-time upper bound on a specific graph projector**. It does not prohibit nonlocal, coarse-grained, dynamically generated relational space or a topology with cycle lengths growing with N.
- It **does not apply to the pure coherent neutral sea**; a larger quantum information resource is necessary to evade this kind of bound, but is not sufficient to create spacetime.
- It depends on the complete labelled Hilbert space, energy conservation, positive mass and rewiring, and a fixed density Q/N.
- No claim about a permutation-gauge quotient, open-system cooling, measurement/postselection, infinite-temperature thermodynamics, or the higher-density ground state is made.
- No new conceptual principle is asserted: the proof combines unitary invariance of operator norm, Markov's inequality and the parent's cycle counting. Independent scientific novelty is not established.

**Next Work test:** In the same Q=3N state, calculate how candidate low-min-entropy and coherent preparations differ in *dynamically generated* loop packing, graph distance/volume growth, and relational source–detector response; do not infer formation from the mere ability to evade an entropy obstruction.
