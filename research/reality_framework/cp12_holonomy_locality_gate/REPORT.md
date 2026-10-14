# CP-12: graph connectivity is not yet spatial gauge geometry

**Additive independent mathematical audit — 14 October 2026 UTC.** Parent: [full-ground connectivity c3c3de3](https://github.com/ryangrvr/TestingGRUT/commit/c3c3de35729b175aab01fe45624ae65a10ad0e9d). No modification to frozen CP-12, canonical GRUT-RAI, Stage-3 Card 3, or the parent's numerical records.

## Verdict on the pinned ground-connectivity packet

The parent's all-irreducible-representation regulator-L=3 calculation gives P(connected)=0.94903651 at N=5 and 0.86697123 at N=6, versus *uniform labelled Q-edge graph* baselines 0.88095238 and 0.73126873. This is compelling **finite-system numerical evidence** that the true small regulated ground is *more* connected than uniform graphs. The previous disconnected coherent matching energy witness does not identify the actual ground-state graph distribution. The parent independently constructs signed S_N irreps (11 at N=6), reports a 36,801,265-dimensional physical Hilbert space, numerically locates the all-irrep lowest energy, and distinguishes the N=3 exact infinite-rotor certificate from N=5/6 numerical regulator checks. It does **not** prove thermodynamic connectivity or a 3D spacetime phase.

## Theorem A: single-cycle and local-holonomy obstruction, independent of dynamics

For every connected undirected simple graph with N vertices and Q occupied edges, the graph cycle-space rank is

\[
b_1(G)=Q-N+1.
\]

Therefore **every** connected Q=N graph, regardless of whether it is a spanning ring or a short-cycle-plus-trees graph, has exactly **one independent cycle**.

Define the graph-diagonal, permutation-invariant bounded observable

\[
\Lambda_R(G)=\frac{1}{N}\#\{v:\ v\ \text{lies on a simple graph cycle of length}\le R\}.
\]

A connected Q=N graph is unicyclic, with a single simple cycle of length ell(G), so for fixed R

\[
\boxed{\Lambda_R(G)\le R/N\qquad\text{(every connected Q=N graph).}}
\]

Consequently, for **every** quantum state rho of frozen CP-12, including arbitrary fermion/electric correlations and coherent graph superpositions,

\[
\boxed{\mathrm{Tr}(\rho P_C\Lambda_R P_C)\le\frac{R}{N}\mathrm{Tr}(\rho P_C).}
\]

If P_C has nonzero probability, the expected fraction of vertices on short microscopic cycles **conditional on connectedness** is <=R/N and tends to zero at fixed R. No Hamiltonian, electric regulator, spectral gap or equilibrium assumption enters this theorem.

This excludes an interpretation of the connected Q=N microscopic **graph itself** as a homogeneous short-loop-rich spatial lattice at large N (for example an ordinary cubic or diamond lattice where a positive fraction of sites lie on elementary plaquette loops of bounded length). It does NOT exclude effective loop observables emerging through graph superpositions, relational coarse-graining, matter correlations, or spatial organizations whose microscopic graph does not supply local plaquette cycles.

## Theorem B: exact count of independent CP-12 electric circulations

For each present edge e=(i,j) with i<j, write its edge-center hole number as h_e=1-n_{b_e}, and set x_e=E_{ie}. Edge Gauss law fixes E_{je}=h_e-x_e.

For each vertex i, vertex Gauss law gives

\[
\sum_{e=(i,j)} x_e-\sum_{e=(j,i)}x_e
=n_{a_i}-\sum_{e=(j,i)}h_e .
\]

Writing B_G for the oriented N-by-Q incidence matrix, this is B_G x=b. A necessary and sufficient condition for an integer solution is Gauss neutrality on **each** connected component (sum of b_i in every component zero). If nonempty, the entire integer solution set is an affine lattice of rank

\[
\boxed{\operatorname{rank}_{\mathbb Z}\ker B_G=Q-N+C(G),}
\]

where C(G) is the number of components. This follows from a spanning forest: choose one integer value for every chord not in the forest, then solve the tree flows uniquely over integers. A connected Q=N graph therefore has **one freely specifiable integer circulation variable at fixed matter occupations**, even though there are 2Q electric half-leg variables before constraints. In Q=2N connected graphs there are N+1; in Q=3N connected graphs there are 2N+1.

These numbers concern *independent pure electric circulations at fixed graph and fermion occupations*, not the total dimension of the full coupled quantum Hilbert space. A superposition over graphs and matter arrangements need not be reducible to one collective gauge oscillator; new physical mechanisms would have to demonstrate local field emergence.

## Corollary C: connection to the previous near-degree-two ground exclusion

For each connected unicyclic G, let b(G)=N-ell(G) be the number of bridges and D_2(G)=sum_i |d_i(G)-2|. Building G by attaching b(G) tree edges to its unique cycle proves

\[
\boxed{D_2(G)\le 2b(G).}
\]

Conditionally on the earlier full-H **theorem candidate** asserting ground probability P(D_2<=K_N)->0 for all K_N=o(N), it follows immediately for any K_N=o(N):

\[
\boxed{\Pr_{\rm ground}(C(G)=1,\ b(G)\le K_N)
\le\Pr_{\rm ground}(D_2(G)\le2K_N)\to0.}
\]

Thus if the connected ground-state weight persists at nonzero probability in the thermodynamic limit, it cannot be concentrated on connected graphs with only subextensive bridges. The asymptotic density and true connectivity remain unproved. This corollary is contingent on final acceptance of the earlier two-star/near-degree-two proof, not on the N=5/6 numerical fit.

## Finite observation: mild bias toward longer cycles, but no continuum

For connected labelled unicyclic graphs, the exact number with a unique cycle of length ell is

\[
\boxed{U_{N,\ell}=\binom N\ell\frac{(\ell-1)!}{2}\,\ell N^{N-\ell-1}
=\frac{N!}{2(N-\ell)!}\,N^{N-\ell-1},\quad 3\le\ell\le N.}
\]

Derivation: select the ell cycle vertices, cyclically arrange them, then attach an oriented rooted forest on the remaining vertices with these ell specified roots (Cayley's generalized forest count). The final integer expression is well defined even for ell=N.

This exactly reproduces enumerated connected-graph cycle distributions: N=5: {3:150, 4:60, 5:12}; N=6: {3:2160, 4:1080, 5:360, 6:60}; N=7: {3:36015, 4:20580, 5:8820, 6:2520, 7:360}.

Comparison with the parent's ground-state analysis (Q=N, L=3):

| N | Uniform connected mean cycle | Quantum connected ground mean cycle | Conditional ring probability / uniform baseline |
|--:|--:|--:|--:|
| 5 | 3.37837838 | 3.51412513 | 1.84936 |
| 6 | 3.54098361 | 3.70865083 | 2.26572 |

Thus although short cycles dominate the finite quantum ground, it is **mildly biased to longer cycles relative to the uniform labelled-graph baseline**, not purely toward the shortest possible cycle. This does not isolate the physics of selection from matter/electric state-count factors; the uniform graph baseline is NOT a uniform physical-Hilbert-space ensemble.

At N=6, the connected ground's mean bridge count conditional on connectedness equals 6-3.70865083=2.29134917. Multiplying by connected probability 0.86697123 gives 1.98653381, independently matching the saved \`bridges_on_connected_graph\` graph-projector measurement.

## Scientific implications and next tests

1. Connectivity and an extensive largest component are distinct from spatial dimensionality.
2. The connected Q=N graph-as-geometry route has exactly one independent microscopic cycle/holonomy, and zero density of bounded-length local cycles. This is a rigorous structural obstruction to conventional homogeneous lattice gauge geometry, **not** a general no-go for state-derived relational spatial geometry.
3. The next spatial/gauge-state investigation should study same-state, same-Q connected and irregular competitors with **extensive cycle rank**, especially Q=2N (where the degree-four near-regular ground result is negative but irregular phases remain possible) and Q=3N (the usual cubic-lattice edge density).
4. Apply an operational gate: joint volume growth, spectral dimension, independent gauge-invariant detector causality, and low-energy transport. Do not fit dimension from one graph silhouette; do not conflate the permutation-singlet property with gauge redundancy.
5. The complete N=5/6 eigensolves were **read and compared from published results**, not independently rerun here. New independent controls are purely graph combinatorics and Gauss linear algebra; they do not change the quantum spectrum.

Run \`python check.py\` from this directory (Python standard library). It checks combinatorial cycle counts against independent exhaustive graph enumeration and verifies integer Gauss incidence-rank on connected and disconnected example graphs.
