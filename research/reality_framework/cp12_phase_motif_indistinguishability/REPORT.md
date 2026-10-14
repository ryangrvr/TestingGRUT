# CP-12: topologically different phase masks indistinguishable to second order

**14 October 2026, independent adversarial audit of [Work checkpoint 86f1d77](https://github.com/ryangrvr/TestingGRUT/commit/86f1d77c544e001dc0fc119b069442e655e625d0).** Additive research only, full frozen CP-12. No changes to canonical GRUT-RAI, Stage-3 Card 3 or CP-12. The explicit six-event exhaustive checker in this directory was run locally. Parent's many-body finite-time eigensolves were inspected in the report but **not independently rerun**.

## Verdict on the parent

The parent derives the exact full-H four-event square-probability coefficients \(-1/810\) versus \(+1/1215\) for two pure preparations with equal graph marginals, purity, neutral-basis relative-entropy coherence, total energy and energy variance. Its Taylor coefficients are unregulated short-time results. Finite time \(T=20\) uses an electric cutoff; \(Q=3N\), \(N\le256\), is a short-time coefficient calculation, not a many-body evolution at those N. These distinctions are maintained.

## A sharper topology-blindness corollary

Fix N,Q and let S be a set of a marked potential edges out of M=binom(N,2). Prepare a uniform neutral-sea superposition with phase \(z^{k_S(G)}\), where \(k_S(G)=|E(G)\cap S|\), \(|z|=1\). Choose the same phase z and marked-edge count a for two preparations. If they also share
\[
W(S)=\sum_v\binom{\deg_S(v)}2,
\]
then for **every** permutation-invariant graph-diagonal observable O, full frozen CP-12 gives
\[
\boxed{\langle O(T)\rangle_A-\langle O(T)\rangle_B=O_N(T^3).}
\]

The equality includes their \(T^0,T^1,T^2\) Taylor coefficients. No restriction to cycle counters, regularity or connectivity is needed. It does not claim equality at higher orders, of labelled observables, or at finite T. At each fixed finite N the neutral initial vectors have finite energy moments and the derivatives are well defined despite unbounded electric rotors.

**Proof:** [Parent DERIVATION.md, sections 1–4](https://github.com/ryangrvr/TestingGRUT/blob/86f1d77c544e001dc0fc119b069442e655e625d0/research/reality_framework/cp12_energy_matched_preparation/DERIVATION.md) derives that all neutral-preparation graph-observable response derivatives through \(T^2\) equal those of the neutral Johnson-graph Laplacian R. Its two-color quotient gives a graph probability \(P_T(G)=D^{-1}[1+q_1(k)T+q_2(k)T^2+O_N(T^3)]\), with \(q_1\) linear and \(q_2\) **quadratic** in \(k=k_S(G)\) (given fixed M,Q,a,z). Therefore equality of \(\mathbb E_{\rm unif}[O k]\) and \(\mathbb E_{\rm unif}[O k^2]\) suffices. The permutation action on any one potential edge is transitive. Its action on distinct unordered pairs of potential edges has exactly two orbits: incident and vertex-disjoint. Thus the averages \(\mathbb E[O n_e]\) are the same for all e, and \(\mathbb E[O n_en_f]\) depend only on whether e and f share a vertex. The number of incident pairs in S is exactly W(S), and all other pairs number binom(a,2)-W(S). Hence two masks with the same (a,W) have identical 0th–2nd full-H graph-observable response derivatives for **every** vertex-permutation-invariant O.

This is a direct exact corollary of the parent's previously proved quadratic reduction. The new contribution is a particularly hostile **topologically distinct**, exactly computed pair and a broader all-observable indistinguishability interpretation. It is not a new fundamental physical principle.

## Explicit six-event counterexample pair

Use N=Q=6, M=15, a=7 and b=8. Both masks have marked vertex-degree multiset (3,3,2,2,2,2), W=10; neither is a vertex relabelling of the other.

- Mask A is C6 plus diameter (0,3): edges (0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,3). It contains **zero triangles and two simple four-cycles**.
- Mask B consists of two triangles linked by one bridge: edges (0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3). It contains **two triangles and zero simple four-cycles**.

Choose the parent's energy-matching phase \(z=\cos\theta+i\sin\theta=-7/8+i\sqrt{15}/8\), because
\[
\cos\theta=1-\frac{M(M-1)}{2a(M-a)}=-7/8.
\]
Both **full-H** initial states have exactly the same graph probabilities 1/5005, pure-state min-entropy zero, neutral-basis coherence log2(5005), and the same two energy moments:
\[
\boxed{\langle H\rangle=84/25=3.36,\quad \operatorname{Var}(H)=97527/8125.}
\]
Variance equality is the parent's permutation-of-potential-edge-coordinates argument for the neutral Johnson Laplacian combined with \(\|T_{\rm hop}\psi\|^2=2Q\), not a claim about higher moments.

An independent symbolic evaluation of the neutral Johnson amplitude gives, for each graph with k marked edges, **the exact full-H T² graph probability correction**
\[
\boxed{q_2(k)=\frac{203k^2-1133k+1386}{24000}.}
\]
Note \(D^{-1}\sum_{G:|E|=Q}q_2(k_S(G))=0\) exactly.

Independent exhaustive enumeration of all C(15,6)=5005 graphs, with exact rational q₂, produces for **both** masks:

| Invariant graph observable O | Initial expectation | Full-H T² coefficient |
|---|---:|---:|
| triangle count | 80/91 | -29/35750 |
| simple four-cycle count (count each distinct cycle) | 45/91 | -87/357500 |
| vertex degree-wedge count \(\sum_v\binom{d_v}2\) | 60/7 | -29/13000 |
| connected-graph projector | 732/1001 | 29/66000 |
| exact degree-two projector | 2/143 | 203/4290000 |

The enumeration checker uses a distinct small exact rational representation \((a,b)\mapsto a+i\sqrt{15}b\) and constructs q₂ from **two applications of the Johnson graph generator** to the phase amplitudes, then enumerates graph observables. It is not a full many-body time-evolution simulation. The theorem transferring these second derivatives to full CP-12 is explicitly proved in the parent DERIVATION and was separately inspected here.

## Interpretation: state selection needs more than simple correlation statistics

- Amount of purity, coherence, initial graph distribution and the first two energy moments do not fix geometric response (as Work proved).
- Even a **different topology of the external phase mask**, while matching its degree statistics, cannot affect *any unlabeled graph observable's first two time derivatives* in this family. This is a rigorous finite-time **leading-order blindness** of the specified controller/preparation family, not universal blindness of CP-12 to higher correlations.
- The first dynamically distinguishing invariant graph coefficient, if any, must occur at order **three or higher**. Proving which coefficient is first in the **full Hamiltonian**, and whether it remains finite or relevant as N grows, is a targeted next test.
- A controller selecting phases based on supplied edge labels is not an internally selected state. Independent topology selection, emergence of locality, long-range dimension and an experimental prediction all remain open.
- The illustrative pair lies in **Q=N**, where conventional 3D local gauge capacity is obstructed; it is an adversarial test of the *phase-control mathematics*, not a prospective 3D physical phase. The high-density \(Q=3N\) question requires its own same-state and finite-time investigation.

**Status:** Exact finite combinatorial certificate + full-law second-order reduction carried from the published parent, independently checked. No unregulated finite-time claim or resolved ToR selection mechanism.
