# CP-12 destination-uniformity and first-rewire locality audit

**13 October 2026. Independent additive research; partner review pending.** Frozen CP-12 unchanged; branch forked from GRUT III adjudication `694463f`. No change to GRUT-RAI or Stage-3 Card 3.

## Exact single-host theorem (full frozen CP-12)

Fix one labelled graph G with N events, M=N(N−1)/2 possible connections and Q present connections. Let |ψ_G⟩ be **any physical matter/U(1) state supported on G**, including arbitrary coherent combinations of gauge/fermion states. For each present edge e, write P_e^neutral for the projector onto removable-edge configurations (edge fermion occupied, both electric legs zero), and r_G=Σ_{e∈G}⟨ψ_G|P_e^neutral|ψ_G⟩.

For any particular absent edge f, the full Hamiltonian's *only* first-order graph-changing term removes a neutral edge e and creates f at constant amplitude −Γ/M. Different removed edges lead to orthogonal final graphs, and the graph swap preserves the physical matter/gauge configuration. Consequently

`||P_(f present) H_12 |ψ_G>||^2 = Γ² r_G/M²`

and, for T→0 at **fixed N**,

`P(f occupied; T) = Γ² r_G T²/M² + O(T³)`.

For a real initial state, real H and graph readout, time reversal improves O(T³) to O(T⁴). **Every initially absent f has the same leading coefficient**. Therefore, for any chosen subset S of absent edges, the fraction of leading occupation creation in S equals exactly `|S|/(M−Q)` when r_G>0. Charges and electric flux change r_G (what can be removed) but do not preferentially choose a destination at this order. Matter hopping is included, not switched off; it cannot change graph occupation in one Hamiltonian action.

Summing over all absent f, the first departure probability has coefficient `Γ² r_G (M−Q)/M²`. At sparse Q=νN and r_G=O(N), it is O(Γ²/N): graph *preservation* at fixed microscopic time is fully consistent with *nonlocal creation when a rewire occurs*.

## Exact ring-locality count

Let G=C_N, Q=N. For a fixed graph distance R with 2≤R<N/2, the candidate absent edges joining sites at initial distance 2,...,R number N(R−1). All absent edges number N(N−3)/2. The share of first successful rewires depositing an edge locally is

`p_local^(leading)=2(R−1)/(N−3) → 0`.

For any fixed-degree graph family with fixed R, local candidate pairs are O(N), while absent pairs are O(N²). Hence a first neutral rewire is overwhelmingly nonlocal **conditional on rewiring**, no matter which physical matter/gauge state lives on the single initial graph.

## Scope, critical exceptions, and next test

1. **Not a universal geometrogenesis no-go.** A coherent superposition over distinct graphs can interfere between transition paths. It must be analyzed separately, including coherence preparation and energy cost.
2. This is the first nonzero T² coefficient at fixed N; the error remainder is not proved uniformly in N and no long-time conclusion is made.
3. 'Local' here is relative to a supplied graph distance, **not** an operationally derived metric; this is an adversarial baseline.
4. No graviton pole, stress-energy Ward identity, Lorentz limit, Einstein equations, or new empirical prediction follows.
5. This is compatible with Work's exact cubic *source-dependent retarded* connection response: the two calculations measure different experiments and orders of time.

**Research recommendation:** register this baseline in the same-sector/same-state CP-12 charter. Before asserting formation of localized spacetime, show that correlated quantum dynamics overturns this leading random-destination behavior on physically generated relational shells. In parallel require independent source/detector momentum, z=1 two-TT propagation, stress-energy universality, Weinberg–Witten/Marolf applicability analysis, and GW170817-compatible speeds.

The standalone stdlib `check.py` uses exact first-swap counting and verifies Gauss-law-valid neutral and one-charge/one-flux locked states; it writes `results.json` and emits `CP12_DESTINATION_UNIFORMITY_EXACT_COUNT_PASS`. No full finite-T quantum propagation is claimed.
