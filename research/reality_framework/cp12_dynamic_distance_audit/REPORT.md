# CP-12 independent audit: charge-defined distance, neutral shortcuts, and a conditional locality timescale

**2026-10-11.** Independent additive research, separate \`codex/grut-cp12-dynamic-distance-audit\` branch forked from frozen [CP-12 commit 0014d615](https://github.com/ryangrvr/TestingGRUT/commit/0014d6153a03a27cd1212ac176b77ceb67a07830). This does not edit CP-12, earlier scouts, GRUT-RAI or Stage-3 Card 3, and does not open CP-13. The original CP-12 authors independently tested up to 2,866,255 states and derived a graph-departure bound. Our new contributions are the **explicit charged-state neutral-shortcut family**, a **bounded-detector corollary**, and a sharper **uniform \(Q=N\) ring-counting baseline**. Their ingredients are established graph/gauge algebra; no scientific priority or empirical novelty claim is justified.

## 1. What the frozen law actually supplies

CP-12 has Q conserved connections among M=N(N-1)/2 possibilities; vertex and edge fermions; two compact U(1) half-link rotors on each occupied connection; and the exact Gauss constraints
\[
G_i=\sum_{e\ni i}E_{ie}-n_i=0,\qquad
G_e=-E_{ie}-E_{je}-n_e+1=0.
\]
A present edge can be removed by the graph-reconnection Laplacian only if it is neutral: its edge fermion is filled and both half-link electric fields vanish. An absent edge is likewise neutral. Each allowed undirected old/new bond swap has matrix element \(-\Gamma/M\), irrespective of the old graph's shortest-path separation.

The frozen short-time neutral-edge removal probability is \(O(T^2)\), versus \(O(T^4)\) for a specifically prepared single charge–electric excitation, due to the extra fermionic de-excitation step. The fixed-charge **bare diagonal** energy minimum on graph G is
\[
E_{\min,G}(h,i)=\kappa Q+2m+\frac{\epsilon_E}{2}d_{\rm inc}^G(h,i),
\]
for an occupied edge hole h and occupied vertex i. This is a precise, conditional distance-sensitive bare electric energy, not a full screened force law or metric.

## 2. An explicit family: Gauss-respecting neutral shortcut lowers available distance

For even N>=8 take the connected cycle G with edges (j,j+1 mod N), choose one hole at edge h=(0,1), one vertex charge at i=N/2, and an integral unit electric flux along the shortest incidence path from the edge center of h to i. The shortest path has length \(d_{\rm old}=N-1\), and its signed electric fields satisfy **all CP-12 Gauss constraints**. In particular they vanish on e=(N/2+1,N/2+2), whose edge fermion remains filled.

Now remove this neutral occupied edge e, and add the absent chord f=(1,N/2). The new graph G' is still connected, has **exactly Q=N links**, and is accessible by one *permitted* CP-12 reconnection. All fermionic occupations and all nonzero E-fields remain unchanged; the new chord has zero electric flux and a filled edge mode. Consequently every Gauss constraint stays satisfied after the rewire.

The state-to-state Hamiltonian matrix element is **exactly \(-\Gamma/M\)** for the neutral graph swap, and both computational-basis **diagonal** matrix elements are equal. Their \(\kappa Q\), mass, electric \(\sum E^2\) and positive graph-Laplacian degree terms are equal because the number of neutral removable edges is unchanged. The transition does **not** instantly release electric energy.

Yet G' admits a *different* flux configuration with only a three-incidence-link path h -> vertex 1 -> chord f -> charge vertex i. Thus
\[
\boxed{d_{\rm inc}^{G}(h,i)=N-1,\qquad
       d_{\rm inc}^{G'}(h,i)=3.}
\]
The **minimum available** bare electric energy drops by \(\epsilon_E(N-4)/2\). Reaching that lower minimum requires real electric/fermion dynamics and is not guaranteed by the single rewire. The result is that the *instantaneous-distance energy landscape* can change by O(N) even when a permitted neutral topology move is diagonal-energy-degenerate.

| N | Initial distance | Distance after one allowed neutral rewire | New minimum bare-electric saving (\(\epsilon_E=1/4\)) |
|---:|---:|---:|---:|
| 8 | 7 | 3 | 0.5 |
| 12 | 11 | 3 | 1 |
| 24 | 23 | 3 | 2.5 |
| 48 | 47 | 3 | 5.5 |
| 96 | 95 | 3 | 11.5 |

Independent [check.py](check.py) computes the integer flux, validates every vertex/edge Gauss constraint before and after each rewire, checks original/new graphs are connected, equality of the complete frozen-Hamiltonian diagonal for the same flux state, and both shortest distances. This is **not** an exact many-body transition rate or demonstration of rapid remote signalling. Indeed the next result shows why the matrix element suppression matters.

### Physical interpretation

The original flux-carrying links on the charged path are locked until discharged, as CP-12 intended. A **different neutral edge** can still be removed to introduce a distant shortcut. Gauss neutrality protects the *excited links*, not a fixed graph distance between charges. If one wants a stable emergent metric, either the rates of such shortcut creation must vanish in the relevant operational limit, or a physical mechanism must penalize their creation or suppress their observable consequences. An explicit distance-penalty inserted by hand would assume the very geometry the theory seeks to derive.

## 3. The balancing exact bound: prepared-graph dynamics becomes frozen at fixed time

CP-12's own [PHYSICAL_DERIVATIONS.md](https://github.com/ryangrvr/TestingGRUT/blob/0014d6153a03a27cd1212ac176b77ceb67a07830/research/constructive_physics/attempt12/PHYSICAL_DERIVATIONS.md) §4 derived a bound on graph-departure probability. Here we extract a useful **uniform bounded-observable corollary** that accounts for *virtual return paths*, not just final graph leakage.

Let P_G project onto one graph G, but retain **all fermion and permitted electric states** on it. Decompose the full Hamiltonian H as \(H_0=P_GHP_G+P_G^\perp HP_G^\perp\) plus the off-diagonal graph-changing block \(B=P_G^\perp H P_G\). Rewiring is the only graph-changing term and the distinct destination graphs are mutually orthogonal. If r(G,\chi) of the Q occupied edges are neutral in an input matter/flux basis state \(\chi\), the outgoing squared amplitude is \(r(G,\chi)(M-Q)(\Gamma/M)^2\). Hence
\[
\|B\|\le b_N:=\frac{\Gamma}{M}\sqrt{Q(M-Q)}
\]
independently of the electric cutoff. In fact it can be tightened to the maximum removable-edge count in the initially selected graph sector, but the universal Q bound suffices. Because \(\|B+B^\dagger\|=\|B\|\), Duhamel's formula for the self-adjoint operators H,H0 yields for **any** initial \(|\psi_G\rangle\in P_G\mathcal H\)
\[
\boxed{\left\|e^{-iHT}|\psi_G\rangle-
                 e^{-iH_0T}|\psi_G\rangle\right\|\le |T|b_N.}
\]
Thus, for **any bounded observer/detector operator** O with \(\|O\|\le1\),
\[
\boxed{|\langle O\rangle_H(T)-\langle O\rangle_{H_0}(T)|
\le 2|T|\frac{\Gamma}{M}\sqrt{Q(M-Q)}.}
\]
The comparison evolution H0 is the CP-12 physical within-G fermion/gauge Hamiltonian including the diagonal contribution of the positive rewiring Laplacian; it merely sets **off-block** graph hops to zero. No commute, perturbation-series truncation, numerical gauge cutoff or chosen geometric detector is assumed.

For Q=νN and M~N²/2,
\[
b_N\sim\Gamma\sqrt{\frac{2\nu}{N}},\qquad
|\Delta\langle O\rangle|\lesssim 2\Gamma |T|\sqrt{2\nu/N}.
\]
Therefore **the effect of all possible graph swaps on any bounded observable tends to zero at fixed external T for any initially definite graph**, even though the Hamiltonian has nonzero old-distance-independent shortcut matrix elements.

The same derivation bounds the probability of leaving the initial graph by \(T^2b_N^2\). Achieving a graph departure probability \(p_*>0\) *requires*, though is not guaranteed at, a time
\[
T\ge \sqrt{p_*}/b_N\sim \frac{\sqrt{p_*}}{\Gamma}
                  \sqrt{\frac{N}{2\nu}}.
\]
This is an exact **conditional locality–geometrogenesis time-scale tradeoff** for this 1/M kinetic normalization: on fixed microscopic times the prepared graph freezes relative to graph change, while preparing a different global geometry with O(1) probability takes a diverging time at least O(sqrt N).

**Boundaries:** the statement is for normalized states initially confined to *one graph*, bounded norm-one observables, fixed external time measured in CP-12's supplied units, fixed Γ, and Q=νN. It is **not** a theorem of microscopic locality for the full fluctuating ground state or thermal Gibbs state, nor for coherent initial superpositions over distinct graphs or times T~sqrt N or larger. The within-G graph can itself be a hubbed nonspatial graph; its propagation speed is not assumed finite. The result *does* constrain how much virtual topology-changing paths can influence a prepared spatial state at finite time. No emergent Lorentz invariance or physical clock has been earned.

## 4. The overlooked \(\nu=1\) counting baseline is still hostile to a single ring

CP-12's extensive quantum relative-entropy thermal bound excludes a thin spatial-family class for ν>1, but gives no exclusion theorem at ν=1. This doesn't make the ν=1 simple N-cycle typical. For the **uniform classical graph measure** at exact Q=N, the probability of a single connected N-cycle is
\[
p_N=\frac{(N-1)!/2}{\binom{\binom N2}{N}}.
\]
A direct Stirling expansion, retaining the fact that M=N(N-1)/2 rather than N²/2, gives the sharp asymptotic
\[
\boxed{p_N=\pi e^2\left(\frac{2}{e^2}\right)^N
(1+O(1/N)).}
\]
Thus a single spanning ring has an **exponential** counting disadvantage even at ν=1, although the leading N ln N terms cancel. Verified exactly for N=4,5,6 against CP-12's 0.2, 0.047619, 0.011988 values and tested numerically through N=1,000.

More broadly, a uniform Q=N graph has many isolated vertices. For M=N(N-1)/2 and Q=N,
\[
p_1=\frac{\binom{M-(N-1)}{Q}}{\binom{M}{Q}},\quad
p_2=\frac{\binom{M-(2N-3)}{Q}}{\binom{M}{Q}}.
\]
If I is the number of isolated vertices, \(E I=Np_1\sim Ne^{-2}\), and
\[
\mathrm{Var}(I)=Np_1(1-p_1)+N(N-1)(p_2-p_1^2)=O(N).
\]
Chebyshev therefore yields
\[
P_{\rm uniform}(G\ {\rm connected})
\le P(I=0)
\le \mathrm{Var}(I)/(E I)^2=O(1/N).
\]
This is a conservative bound; the actual uniform connected probability is smaller. For N=1000, \(E I=134.79\) and the bound is below 0.00442.

**Scientific scope:** These are exact *reference graph counts*, not probability statements about the **interacting CP-12 quantum ground state or Gibbs state**. At ν=1, a sufficiently strong physical energy advantage O(N) could in principle overcome the O(N) entropy cost at sufficiently low temperature. The original quantum thermal exclusion theorem remains silent there. CP-12's computed N-cycle probabilities 0.28956 (N=4), 0.09487 (N=5) and 0.03220 (N=6 invariant sector) exceed uniform counts but decrease over the tested sizes; no thermal thermodynamic extrapolation follows.

## 5. Physical decision and next discriminating test

CP-12 supplies a genuine finite model of matter-gauge back-reaction and an electric-path energy diagnostic. These achievements survive this independent critique.

It currently **cannot identify the instantaneous shortest-path electric cost with a stable spacetime metric**: permitted graph swaps can generate new shortcuts at no bare diagonal energy premium proportional to their old-graph distance. Conversely, because of the 1/M scaling, all topology-changing effects vanish uniformly on bounded prepared-graph observables over fixed microscopic time when Q∝N. That opens a potential *conditional semiclassical locality regime*, but **postpones geometrogenesis to increasingly long times** rather than solving it.

The decisive next physics experiment is one properly prepared protocol using an explicit charged reference and detector on an initially spatial-looking graph: compute their finite-time gauge-invariant response, compare to the **within-G** (graph-frozen) response as N increases with Q/N fixed, and then separately measure graph assembly over T scaled as sqrt N. These are two distinct hypotheses and should not be merged. Include adverse hub-rich graphs and superposed-graph initial states. Treat the absence or persistence of long-distance detector tails as the observed outcome. Do not infer Einstein dynamics from graph distance or add an ad hoc local rewiring metric penalty.

No original physical law or new empirical prediction is claimed. The uniform counting and graph-freezing corollaries are consequences of known combinatorics and quantum evolution inequalities specialized to frozen CP-12. A strict novelty review of these model-specific results remains outstanding.

## Reproduce

[check.py](check.py) uses only Python's standard library and NO CP-12 solver internals. It verifies Gauss-valid integral flows on explicitly constructed rings, connectivity, neutral rewiring legality, *equal full Hamiltonian diagonal values* of the two graph basis states, geodesic distance and available bare energy shift, the exact N-cycle combinatorics and uniform isolated-vertex bound, and numeric time-scale coefficients through N=1000. It writes results.json when run. **The complete 2,866,255-state frozen CP-12 eigensolver is not rerun.** The results here are analytical finite-graph corollaries, not claimed independent reproduction of every CP-12 numerical assertion. The exact committed script remains available for verbatim partner execution.
