# CP-12 general regular-sector audit — new conditional theorem

**Date:** 2026-10-13. **Scope:** independent additive math and numerical adversaries only. Based on the *unchanged* CP-12 Hamiltonian; no modification to canonical GRUT-RAI, Stage-3 Card 3 or other frozen branches. Parent: `0463ec1` (classical degree-two fragmentation theorem); companion comparison: `43d706b` (exact ring-sector graph cut). **Not peer-reviewed or partner-audited.**

## Result A — exact switching bound for *any* fixed degree k

Let labelled simple graphs have N vertices, Q=kN/2 edges, with 2 <= k <= N-2 and kN even. Define W(G)=sum_i binom(d_i,2), and a separate classical Gibbs ensemble p(G) proportional to exp[-a W(G)] with a=beta lambda finite and fixed. All k-regular graphs minimize W.

From each k-regular graph choose an oriented edge (u,v) in kN ways and choose w not adjacent to u in N-1-k ways. Replace (u,v) with (u,w). The degrees at v,w become k-1,k+1, and W increases by **exactly one**. Each defective graph has at most k+1 regular predecessors (the possible edges incident on its unique degree-(k+1) vertex). Consequently the number of distinct one-pair defects is >= kN(N-1-k)/(k+1) times the number of regular graphs. Their Boltzmann weights prove

\[
\boxed{\Pr_a(\text{k-regular})\le\left[1+\frac{kN(N-1-k)}{k+1}e^{-a}\right]^{-1}.}
\]

For fixed k and fixed finite a, exact regularity tends to zero at least as N^-2. This is a **classical surrogate theorem**, not CP-12 quantum thermodynamics. The CP-12 fourth-order coefficient is uncontrolled at the frozen hopping t=1, so it cannot be substituted for a in a physical Gibbs prediction.

## Result B — *actual neutral graph-swap* quantum block-cut bound

Let P_k project onto all k-regular graph configurations, retaining all their fermion/U(1) gauge states, and C_k=I-P_k within Q=kN/2. Under CP-12's frozen neutral occupied-to-empty edge-swap off-diagonal `-Gamma/M`, with M=binom(N,2), every k-regular graph has D=Q(M-Q) distinct external graph-swap neighbors. A single swap **never** maps between distinct k-regular graphs. Every external neighbor has at most k+1 k-regular predecessors. The unrestricted scalar graph-cut adjacency therefore has squared norm <=(k+1)D by a nonnegative matrix row-sum bound.

Matter and electric neutrality constraints delete graph-changing matrix elements; graph-preserving matter hopping adds none. In the common universal occupation basis, this bounds the physical cross-sector operator at every cutoff:

\[
\boxed{b_{N,k}:=\|C_kH_{12}P_k\|\le\frac{\Gamma}{M}\sqrt{(k+1)Q(M-Q)}
=\Gamma\sqrt{\frac{k(k+1)}{N}}(1+O(N^{-1})),}
\]
for fixed k.

Duhamel's formula, using the block-diagonal part of the **same** full Hamiltonian, then gives for any initial physical state supported on k-regular graphs:

\[
\boxed{\Pr(\text{outside k-regular at }T)\le \min(1,b_{N,k}^2T^2).}
\]

Thus disconnected regular graphs can be *conditionally dynamically persistent* for T=o(sqrt N), just like regular spatial graphs. Stability of supplied geometry is **not** selection.

If a physical trial state outside P_k has energy below the exact lowest P_k block eigenvalue by delta_N>0, the global ground-state overlap satisfies
\[
\langle P_k\rangle_{\rm ground}\le b_{N,k}^2/(b_{N,k}^2+\delta_N^2).
\]
No such large-N gap premise has been established for CP-12.

## Result C — exact adversarial equal-degree geometry examples

Periodic d-dimensional L^d cubic lattices with k=2d (L>=4) and disjoint copies of K_(2d+1) have identical vertex degrees, Q and W whenever (2d+1) divides L^d:

| Candidate spatial lattice | Fragmented competitor | (N,Q,k,W) |
|---|---|---|
| 1D ring, L=6 | 2 K3 | (6,6,2,6) |
| 2D 5x5 torus | 5 K5 | (25,50,4,150) |
| 3D 7x7x7 torus | 49 K7 | (343,1029,6,5145) |

**These examples match the degree-regularizing contribution, not CP-12's complete energies.** High-order virtual matter/gauge loops can distinguish them; at weak hopping the established C6 versus 2 C3 comparison actually favors the fragmented case.

## Result D — the important entropy competitor at k>=3 is often CONNECTED

For fixed k>=3, standard results of Wormald show nearly all uniformly sampled labelled k-regular simple graphs are connected. Therefore enforcing connectedness as well as degree k is not enough to earn finite-dimensional geometry.

McKay–Wormald's regular-graph enumeration has log(count)=(k/2)N log N+O(N). Any *single* N-vertex lattice isomorphism class has at most N! relabellings, giving the classical uniform bound
\[
\Pr(\text{one specified lattice class})\le\exp[-(k/2-1)N\log N+O(N)].
\]
For k=6 the leading exponent is -2N log N. This bound is only against a specified lattice class in the uniformly weighted regular-graph ensemble, **not against all conceivable spatial geometries or full interacting CP-12**.

A separate controlled geometry comparison at N=343 and k=6, all Q=1029: the 7x7x7 periodic cubic lattice has graph-ball sizes r=0..6 equal to [1,7,25,63,123,195,263], diameter 9. Five seeded connected uniform random 6-regular controls show radius-three balls of size 140–153 and diameter 5. These are graph-metric controls, **not quantum detector measurements**.

## Verification, interpretation and required next test

`check.py` exhaustively enumerates all labelled graphs at (N,k)=(6,2),(6,3),(6,4),(7,2),(7,4), confirms defect counts and the <=k+1 reverse bound; verifies the three exact spatial adversaries and the fixed-k asymptotic formula's leading log coefficient. It does **not** diagonalize full CP-12, independently prove the McKay–Wormald enumeration, or check gauge fiber block matrices. The latter require a separate partner code audit.

**Primary next calculation:** compare the full *frozen-coupling*, rewiring-diagonal-retaining, gauge-constrained conditional energy density of connected low-dimensional lattice families against equal-degree disconnected cliques **and connected nonspatial high-girth regular graphs**; track spectral dimension and gauge-invariant relational responses only on backgrounds whose geometry is actually supported by CP-12. Pay the energy cost of state preparation. Preserve all regulators, controls and negative results. Do not insert a geometric selection reward.

Prior art: N. Wormald, JCT B 31 (1981), https://doi.org/10.1016/S0095-8956(81)80021-4 ; McKay and Wormald, Combinatorica 11 (1991), https://doi.org/10.1007/BF01275671 ; Wilkinson and Greentree, PRD 92 (2015), https://doi.org/10.1103/PhysRevD.92.084007 .

**Verdict: new conditional, model-specific analytical bounds and exact counting controls; no selection of physical spacetime, no Einstein gravity, no observational prediction.**
