# CP-12 — Connectedness, neutral bridges, and a gauge-cost trade-off

**13 October 2026 — independent additive research.** Frozen CP-12 law, complete matter–gauge Hamiltonian, unbounded integer U(1) fields; no change to GRUT-RAI or Stage-3 Card 3. Full standalone audit packet: GRUT_CP12_Connected_Bridge_Audit_2026-10-13.zip (provided in the chat).

## Exact connectivity-exit theorem

Let N>=5, M=N(N−1)/2 and Q=N. Every connected simple N-vertex N-edge graph G is unicyclic, with one cycle of length ell(G) and b(G)=N−ell(G) bridges. Removing a bridge e partitions G into components of sizes s_e,N−s_e. Exactly s_e(N−s_e)−1 **originally absent** edge destinations reconnect those components. The remaining originally absent destinations disconnect G:

\[
\boxed{d_e(G)=M-N+1-s_e(N-s_e).}
\]

Nonbridge removals preserve connectivity regardless of the absent destination. Thus, in the fully neutral sea, the number of one-rewire exits from connected to disconnected graphs is exactly

\[
X(G)=\sum_{e\in\operatorname{Bridges}(G)}d_e(G).
\]

In the **full CP-12 physical Hilbert space** let \(P_e^{neutral}\) project onto the removable occupied-edge configuration (center fermion filled, both electric incidence legs zero). All absent destination edges already have these neutral values. For any physical initial state supported on a single definite connected graph G, the short-time probability of disconnected geometry is

\[
\boxed{
P_{\rm disconnect}(T)=\frac{\Gamma^2T^2}{M^2}
 \sum_{e\in\operatorname{Bridges}(G)}d_e(G)
 \langle P_e^{neutral}\rangle +O_G(T^3).
}
\]

The frozen hopping operator cannot change graph occupation in one step, and distinct removed–created pairs reach orthogonal final graphs. This is not a long-time or uniform-N Taylor bound and does not extend without modification to initially coherent superpositions over several graphs.

## Positive rewiring boundary for every physical state supported on connected graphs

Denote by \(P_C\) the projector onto all connected graphs, with every allowed matter/electric state, within Q=N. The frozen rewiring operator R is a positive Laplacian on each fixed fermion/electric fiber. Its quadratic form on a vector supported on connected graphs is the sum of positive internal graph-difference squares **plus** the squared boundary amplitudes at connected-to-disconnected rewires. Therefore

\[
\boxed{
P_CRP_C\succeq\frac{\Gamma}{M}
 \sum_{G\in C}|G\rangle\langle G|\otimes
 \sum_{e\in\mathrm{Bridges}(G)}d_e(G)P_e^{neutral}.
}
\]

Set
\[
A_N=\max(0,M-N+1-\lfloor N^2/4\rfloor),\qquad
w_N=\Gamma A_N/M\longrightarrow \Gamma/2=0.05.
\]
Since \(s_e(N-s_e)\le\lfloor N^2/4\rfloor\), the connected-sector rewiring penalty is at least \(w_NB_{\rm neutral}\), where \(B_{\rm neutral}\) counts removable bridges.

## Locking bridges consumes actual matter/electric energy

Define \(L_{\rm bridge}\) as the number of bridges that are NOT neutral/removable. For an occupied edge with an edge-center hole \(h_e\in\{0,1\}\) and integer incidence fields, the diagonal operator inequality is

\[
1-P_e^{neutral}\le h_e+E_{ie}^2+E_{je}^2.
\]

The conserved total fermion number implies \(\sum_{\text{occupied }e} h_e=N_v\), the total vertex fermion occupation, so

\[
L_{\rm bridge}\le N_v+\sum_{e,i\in e}E_{ie}^2.
\]

At frozen m=1, epsilon_E=1/4,

\[
H_{\rm mass}+H_E=2mN_v+\frac{\epsilon_E}{2}\sum E^2
 \succeq\frac18 L_{\rm bridge}.
\]

Since \(w_N\le\Gamma=0.1<1/8\), combining this with the positive graph-Laplacian boundary gives the exact **connected compression** inequality

\[
\boxed{
P_C(H_{\rm mass}+H_E+R)P_C
 \succeq w_N\sum_{G\in C}b(G)P_G.
}
\]

**IMPORTANT:** The full Hamiltonian also includes the fermion-hopping operator, which can contribute NEGATIVE energy. This inequality alone cannot decide whether the actual full-H ground state is connected. It constrains the cost of creating a connected correlated phase in the frozen law.

## Topology/regularity relation

For every connected unicyclic G, \(D_2(G)=\sum_i|d_i(G)-2|\le2b(G)\). Proof: start from the unique cycle and attach each of the b tree edges and its new vertex; each step increases the absolute degree defect by at most two. Therefore states with extensive D_2 necessarily have extensive bridges; long cycles with attached paths can conversely have many bridges but only two degree defects.

This relates the previous Q=N near-degree-two **ground-state exclusion candidate** to a new distinct constraint. **Their combination does not yet yield a ground-state connectedness no-go:** the boundary connected-to-disconnected coupling is not known to be sufficiently small.

## Exact independent controls

A standard-library-only checker enumerates all occupied-to-absent swaps and tests the bridge count against brute-force graph connectivity for cycles, triangles with tails, and squares with tails at N=5,...,10. The checker also constructs a **Gauss-law-valid** occupied bridge locked by an edge-center hole, endpoint vertex fermion, and one electric incidence unit, verifying removal suppression. Its recorded status is CP12_CONNECTED_BRIDGE_BOUNDARY_EXACT_PASS.

Examples: at N=6, a pure cycle has 0 outgoing disconnection swaps and a triangle with a tail has 8; at N=9 the counts are 0 and 70. For N=9's triangle-tail case the exact T² coefficient is \(\Gamma^2(70)/36^2=0.000540123456790...\).

## Next genuinely open step

Use this connected-sector inequality together with the exact degree-two exclusion and full-H coherent variational methods to attempt a **controlled projector bound for actual ground-state connectedness**, including the negative hopping term and non-negligible graph mixing. If no bound can be earned, study gauge-invariant relational geometry and propagation in a *specified same-sector state*, distinguishing ground states from supplied preparations. Do not infer that an energy witness is the ground state or claim emergence of Einstein gravity from topology alone.
