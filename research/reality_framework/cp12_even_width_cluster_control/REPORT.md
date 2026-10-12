# CP-12 independent even-width and graph-mixing control: the first winding coefficient is not an odd-cycle artefact

**2026-10-12.** Additive partner review branched from the frozen-coupling CP-12 cluster-expansion result [27010b8](https://github.com/ryangrvr/TestingGRUT/commit/27010b8951cfde6ad06912dd7975ecb79af04fa1). Original files, CP-12 Hamiltonian, GRUT-RAI, Stage-3 Card 3 remain unchanged. This packet reuses the owner's already-computed finite-cluster energies **without any new large block diagonalization** and independently reconstructs their combinatorial host embeddings. It does not claim a new law of physics.

## 1. A physically meaningful confound in the odd-width tube

The author's explicit width-nine comparison is between a bipartite extended square grid and a **non-bipartite** width-nine periodic tube. The latter contains noncontractible odd 9-cycles; the former has no odd cycle. Therefore the first width-sensitive energy contribution tests a mixture of finite circumference and odd-cycle topology. In a fermion-and-gauge model this distinction should not be silently labelled “dimensionality selection”.

The owner already enumerated both the ninth-edge winding loop **and** the lost ninth-edge open path, correctly obtaining
\[
\Delta e_9=\frac{E(C_9)}9-[E(P_{10})-E(P_9)].
\]
At frozen coupling t=1, a=.1 the stored L4 data give +1.0414431046e-9, but the owner then computed an edge-ten correction making the cumulative value **−5.85430687e−9**. This is a robust reversal of **truncated** differences, not a proof of the final host ranking.

To control the *parity* issue at less computational cost, compare square 18×18 even torus with 8×36 periodic tube. Both hosts are bipartite, degree four and have the same local neighbor graph. The width-eight tube first differs at 8-edge winding clusters. Independent control 10×36 verifies no eight-edge differences at width ten.

## 2. Independent exhaustive rooted-injective embedding computation

Use all stored connected bipartite degree≤4 edge clusters of size≤8 in the owner's pinned \`clusters_k8_L3.json\` (136 rows) and \`clusters_k8_L4.json\` (124 embedded-contributing rows). Unlike reusing \`clusters.py:rooted_count\`, [check.py](check.py) implements its own standard-library-only graph adjacency, depth-first rooted injective non-induced graph embedding enumeration, and graph automorphism counter.

For **every row**, the independently recomputed 18×18 square embedding density matches the owner's stored count exactly. For **every row**, width ten and square agree. At **width eight**, exactly two unequal embedding densities occur, both at edge order eight:

| Topology | Rooted maps on square | Rooted maps on width-eight tube | Automorphisms | \(\Delta m\) per vertex |
|---|---:|---:|---:|---:|
| Eight-cycle \(C_8\) | 112 | 114 | 16 | \(+1/8\) |
| Eight-edge path \(P_9\) | 5,916 | 5,914 | 2 | \(-1\) |

There are **no** extra unequal branched clusters through eight edges. Thus the first width-eight-sensitive relative-energy contribution is
\[
\boxed{
\Delta e_8=\frac{w(C_8)}8-w(P_9)
          =\frac{E(C_8)}8-[E(P_9)-E(P_8)].
}
\]
The identity between the two expressions follows from the same exact Boolean linked-cluster subtraction used in the frozen owner report. The relative-energy baseline cancels identically because both expressions carry one bond of energy per site.

## 3. Computed frozen-coupling control

The finite-block numbers, taken directly from the author's **pinned** recorded ground energies (frozen matter hopping t=1, inherited graph-block rewiring diagonal a=.1), yield:

| Electric regulator | \(E(C_8)\) | \(E(P_8)\) | \(E(P_9)\) | \(\Delta e_8\) |
|---:|---:|---:|---:|---:|
| L=3 | −5.383379590363609 | −4.769196950386977 | −5.442119407729160 | **+8.54673143e−9** |
| L=4 | −5.383379590363596 | −4.769196950386970 | −5.442119407729161 | **+8.54674109e−9** |

The independently calculated *weight times new embedding density* result at L4 is **+8.54674031547e−9**, agreeing with the independently telescoped finite-block expression within 8e−16. The two regulator values differ by 9.66e−15, orders of magnitude below the positive first width-sensitive contribution.

**Physical deduction:** the first positive width-sensitive coefficient does NOT require an odd-cycle/frustration difference, because it also occurs for a bipartite tube. But the coefficients are of edge order eight in one host and nine in another; they are **not comparable final energy rankings**. The width-nine owner's tenth-edge sign reversal is positive evidence against inferring an order-independent preference from the new width-eight sign. The later width-eight terms are uncomputed here.

## 4. A separate unearned ingredient: coherent graph mixing

The linked-cluster ground energies are energies of graph-compressed CP-12 blocks. The **same frozen Hamiltonian** also includes off-block neutral graph swaps. At fixed N,Q, at most D=Q(M−Q) off-block swaps depart any basis state and their amplitudes are bounded by \(\Gamma/M\), where \(M=\binom N2\). The off-block matrix V is Hermitian. Its uniform maximum absolute row-sum bound gives by the Schur/Gershgorin norm inequality
\[
\boxed{\|V\|\le \frac{\Gamma Q(M-Q)}{M}.}
\]
Let \(H_{\rm block}\) retain ALL complete fixed-graph blocks and the positive rewiring diagonal, but omit only their offdiagonal graph changes, and let \(E_B\) denote its **minimum over all allowed graph blocks**. A best block eigenstate is a variational trial state for the full H and has zero V expectation; Weyl's bound yields the exact (but loose) statement
\[
\boxed{E_B-\frac{\Gamma Q(M-Q)}M\le E_0(H)\le E_B.}
\]
For Q=2N, \(\Gamma=.1\), N→∞, the uncertainty of this universal bound per vertex approaches \(2\Gamma=.2\). That is astronomically larger than \(\sim10^{-9}\) differences between edge-truncated fixed-graph energies. **The loose bound does NOT mean graph mixing really changes the energy density by .2**; it means we have not certified mixing to be negligible at the desired precision.

Note that CP-12's *fixed-time prepared-graph* detector error O(T²/N) is a response theorem for special initial states. It does NOT bound the energy correction of the full quantum *ground state*, which can be a coherent superposition of all graphs. Conflating these two limits would be a serious inference error.

## 5. A precise next requirement rather than another ambiguous winner

For any two candidate phases A and B, a rank claim must eventually provide, on the **same** thermodynamic prescription:
\[
\Delta e_{\rm total}=\Delta e_{\le K}^{\rm block}
                   +\Delta e_{>K}^{\rm block}
                   +\Delta e_{\rm mix}.
\]
The available reported calculations give \(\Delta e_{\le K}^{\rm block}\) and converged finite-block electric-cutoff diagnostics, **not** bounds on the latter two terms. A decisive sign requires
\[
\boxed{
|\Delta e_{>K}^{\rm block}+\Delta e_{\rm mix}|
<
|\Delta e_{\le K}^{\rm block}|,
}
\]
or an alternative independent bound on the sign of the entire remainder. Neither condition has been established.

The next calculation should first target a **tail bound / empirical reliability check**, not another interpretation of a tiny signed partial sum: record absolute edge-order contributions and embedded-cluster counts at multiple orders, finite-size cross-validation against complete graph-block energies on small periodic hosts (where all edge subsets fit), and if attainable a provable convergence criterion (gapped locality or a controlled parameter). Then separately test quantum off-block mixing with an invariant subspace or variational/response estimate. Preserving t=1 and a=.1 is essential when auditing CP-12's physical ranking; varying parameters may test mathematical convergence but is not a result about frozen CP-12.

The width-eight control is valuable because it is a **same-parity adversarial comparison**, derived entirely from already recorded cluster energies. Its positive first term is *not* evidence of a completed spatial selection mechanism.

## Reproduction and verification

Run from repository checkout:
\`\`\`sh
python research/reality_framework/cp12_even_width_cluster_control/check.py
\`\`\`
The script uses Python standard library only, reads the pinned existing cluster-result JSON files, exhaustively recomputes rooted embeddings and automorphisms, checks every owner's square-grid count, and verifies two numerical expressions and both electric cutoffs. It writes \`results.json\` locally when executed.

The same independent combinatorial calculations and L3/L4 energy identities were executed in a separate JavaScript implementation against source files fetched directly at the pinned commit. The **exact committed Python script was not executed by this assistant in a checkout**; a partner should run it verbatim. No new finite Hamiltonians were diagonalized here. All scientific results inherit the scope and numerical precision of the source's precomputed block eigenvalues.

## Prior research framing

Linked-cluster expansion/inclusion-exclusion is an established many-body method; see Tang, Khatami & Rigol, *Computer Physics Communications* 184, 557 (2013), arXiv:1207.3366, already identified in the owner's analysis. Cyclic bipartiteness is elementary graph theory. No literature novelty or experimental prediction is asserted.
