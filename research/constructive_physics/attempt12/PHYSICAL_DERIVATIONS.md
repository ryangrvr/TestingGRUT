# CP12: constraints, energy, physical records and distance

All results use LAW.md unchanged. These are author derivations and numerical controls, not independent review or established gravity.

## 1. Gauge invariance and stability

An elementary hop a_i† U_ie b_e raises E_ie and the vertex occupation together and lowers the edge occupation. Both vertex and edge-center Gauss operators therefore commute with it. Electric and mass terms commute with Gauss. A permitted rewiring acts only on a neutral, zero-electric-field edge and preserves every constraint. Summing Gauss over a connected component gives its charge neutrality. No charged excitation exists on an isolated vertex in this sector.

For fixed classical phases, the fermion hopping has magnetic incidence matrix B, Tr(BB†)=2Q. Its fixed-M ground energy, including kappa Q, is

 U_phi = kappa Q - Tr[sqrt(m² I+t² BB†)-m I].

Concavity of sqrt yields, with nu=Q/N,

 U_phi/N >= kappa nu - sqrt(m²+2t²nu)+m.

The rotor phase representation is a direct integral of these fermion operators. Restricting to Gauss and inactive-edge constraints cannot lower this bound. The electric and rewiring terms are positive, so the same bound holds for the full CP12 Hamiltonian and its electric-cutoff projections. At the frozen couplings its minimum over nu is -1/4, reached by the bound at nu=3/2. This does not say any graph saturates it. The neutral filled-edge state has energy kappa Q when coherently uniform over graphs, giving a finite extensive trial upper bound.

## 2. Matter excitation delays edge removal

Let a fixed graph G have Q edges and prepare its neutral filled-edge vacuum |0;G>. Select a present e incident on i. The normalized state |r> = a_i† U_ie b_e |0;G> is gauge invariant relative to the vacuum, with a vertex particle, an edge hole and electric field one on their common half-link. Its bare mass/electric cost is 2m+epsilon_E/2. It is an explicit physical distinction, without computing global minimal/maximal causal roles.

The frozen rewiring rule permits removal of e in |0;G> but not in |r>. At short external time T,

 P_remove,e(T;vac) = (M-Q)(Gamma/M)² T² + O(T³),
 P_remove,e(T;record) = (M-Q)(Gamma/M)² t_hop² T⁴/4 + O(T⁵).

For the recorded state, the unique leading two-step path first annihilates the particle-hole/electric excitation and then removes e. These coefficients were checked at four shrinking times, two electric cutoffs, with norm, energy and reversal controls in record_dynamics.py. A coherent partial preparation restores a cos²(angle) fraction of the T² term.

This is reversible excitation back-reaction, not a permanently stable observer memory or measurement outcome. The selection rule also depends on the SUPPLIED rewiring operator: gauge invariance alone would allow more elaborate rewiring that transports flux and charges. It is not a universal gauge-protection theorem. Identifying the initially addressed edge and vertex is supplied; a complete relational apparatus has not been constructed.

## 3. A bare relational separation diagnostic

Fix a graph and the occupation sector with one vertex particle at i and one hole at an active edge center e, all other vertex modes empty and edge modes filled. Gauss requires an integer electric flow from that edge center to i in the bipartite incidence graph. Let d_inc(e,i) be its shortest path length.

Any integral flow of unit divergence decomposes into one or more source-to-sink paths plus cycles. Since integer E has E²>=|E|, sum E_l²>=d_inc. A unit flow on a shortest path attains equality. Thus

 min_E [E_diagonal - kappa Q] = 2m + (epsilon_E/2) d_inc(e,i).

If no path exists, no state with those fixed charges satisfies Gauss. This result uses existing physical charges and does not change Gauss to insert external sources. distance_control.py independently compares all available charge placements with breadth-first distances for the small graph ensembles.

The separation is a property of the instantaneous connection structure and constrained bare energy. The full interacting system may polarize, create pairs and screen; this is NOT a derived interacting static potential, Newtonian attraction, continuum distance or light cone. It identifies an operationally motivated quantity for a future spectroscopy calculation while exposing a confining bare electric energy, not gravity.

## 4. Large-N geometry-change rate

Let P_G project on one initial graph, retaining all allowed matter/gauge states. Only rewiring changes G. Different single swaps from G have orthogonal final graphs. Each is a partial isometry of norm at most Gamma/M. Therefore

 ||(I-P_G) H P_G|| <= (Gamma/M) sqrt[Q(M-Q)].

Compare H with its block-diagonal part relative to P_G. Duhamel's formula for their unitary groups, whose bounded off-block difference has this norm, gives for ANY normalized initial state supported on G

 P_leave,G(T) <= Gamma² T² Q(M-Q)/M².

For Q=nu N this is O(T²/N). The empty-electric vacuum saturates the short-time coefficient. This is a consequence of the supplied 1/M rewiring rate, not emergent relativistic locality. A prepared graph changes little on fixed external time scales as N grows; generating another geometry from it can require diverging times. The bound does not apply to an initial superposition over many graphs and does not establish a thermodynamic relaxation time.

## 5. Quantum preparation is separate from a ground-state solution

quench.py starts from the physical, neutral, equal-amplitude graph state, not the computed ground state. Its energy remains kappa Q. For N=4,Q=4 the probability of a four-cycle changes from 0.2 to approximately 0.262418715 at T=20, compared with the ground-state value 0.289560821. Cutoffs 4 and 5 agree on that probability within 1.6e-8; cutoff 3 is less accurate. The Wilson expectation also evolves and can change sign. There is no bath, irreversible cooling or demonstrated convergence to a selected geometry. Unitarity preserves the initial energy distribution.
