# CP-12: an exact ring-sector coupling bound before attempting large-N emergent spacetime

2026-10-13. **Independent additive audit** forked from frozen CP-12 full coherent follow-up [7be83fc](https://github.com/ryangrvr/TestingGRUT/commit/7be83fc1f3d6a5799f11a587d3df5abc023dda4f). No CP-12 Hamiltonian change, no CP-13 or Stage-3 card, and no merge to canonical GRUT-RAI. This is a strict finite quantum spectral statement about the Q=N graph sector, not a derived physical universe or a claim of novelty.

## Motivation and existing evidence

The frozen law has N vertices, M=N(N-1)/2 possible edges, Q=N conserved occupied edges, fermions, compact U(1) fields subject to Gauss constraints, and a positive neutral-only graph rewiring Laplacian with off-diagonal -Gamma/M between graph states differing by a single occupied/empty edge swap. The full matter/gauge Hamiltonian changes no graph occupation in one action. Rewiring preserves all fermion/field coordinates and only acts when the removed edge is neutral.

The owner's latest fixed-graph (including the diagonal graph-rewiring contribution) energies at frozen couplings for N=6,Q=6 are
\[
E_{\rm ring}(C_6)=-0.635203537847,\qquad
E_{\rm comp}(C_3\sqcup C_3)=-0.637262028038.
\]
The competitor wins by delta_6=0.002058490191. These are **regulator-converged numerical block results at two tested cutoffs**, not infinite-cutoff spectral proof or the global full-graph ground state. A further weak-hopping x^12 triangle advantage is independently reported by the owner. Can we quantify whether full graph mixing could make a ring-dominated ground state anyway without simulating the huge joint Hilbert space?

## The exact spanning-ring / non-ring rewiring cut

Let R be the set of every **labelled simple connected N-cycle**, N>=4, and P_R the projector onto all physical gauge/fermion states whose graph belongs to R. Let C be all other Q=N graph configurations.

**Lemma 1 (no direct ring-to-ring rewires):** Two distinct 2-regular simple graphs cannot differ by a single occupied-edge/empty-edge exchange: removing edge (a,b) changes degrees at a,b by -1, while adding a different edge (c,d) increases degrees at c,d; for every degree to remain two the unordered endpoint sets would have to coincide. Thus \(P_R H P_R\) is the direct sum of complete fixed-graph CP-12 blocks over spanning rings. All their lowest energies coincide by simultaneous permutations, including matter and rotors.

Write \(K_{CR}=P_C H P_R\), the full gauge/fermion-coupled graph-changing block. On the *scalar unrestricted* matching graph, let A_{CR} be the 0/1 adjacency matrix between spanning rings and outside graphs under one bond exchange. Every ring has
\[
D=N(M-N)
\]
outside neighboring graphs.

**Lemma 2 (outside preimage multiplicity):** Any outside graph reached in one swap from a spanning ring is a connected unicyclic graph. It is obtained either by a swap sharing one endpoint (degrees 1 and 3) or by a swap with disjoint endpoints (degrees 1,1,3,3). It has **at most two spanning-ring predecessors**. The degree-1/3 case has at most two repairs, obtained by removing one of the two cycle edges incident on the degree-3 vertex. The 1,1,3,3 case has at most one. The exact swap count for shared-endpoint replacements from each ring is \(2N(N-3)\); those outside graphs have two ring predecessors. All remaining outside neighbors have one ring predecessor.

Since every ring row is equivalent by event permutations, the nonnegative \(A_{RC}A_{CR}\) has a *constant* row sum:
\[
D+2N(N-3)=N(M+N-6).
\]
Its uniform positive ring vector attains the Perron eigenvalue, and a constant row sum bounds every eigenvalue in magnitude. Therefore the **exact scalar graph-cut operator norm** is
\[
\boxed{\|A_{CR}\|_2=\sqrt{N(M+N-6)}.}
\]

**Lemma 3 (full gauge/fermion bound):** Neutral-rewiring moves leave the full fermion occupations and electric values unchanged. In a universal matter/electric basis, the full cross-sector matrix is a direct sum of graph adjacency matrices with some entries *removed* by the gauge/neutrality conditions. Their individual spectral norms are dominated by that of the nonnegative unrestricted matrix A_{CR}. Consequently
\[
\boxed{b_N:=\|P_C H_{12}P_R\|
\le\frac{\Gamma}{M}\sqrt{N(M+N-6)}
\sim\Gamma\sqrt{2/N}.}
\]
This is cutoff-independent: restricting electric configurations cannot increase the universal upper bound. It uses the **supplied 1/M graph kinetics**, not a universally local microscopic principle.

## Ring-sector ground-state weight theorem

Let \(E_R=\inf\sigma(P_R H P_R|_{P_R\mathcal H})\) be the exact minimum of the ring-compressed full quantum block. Let \(E_{\rm comp}\) be any *valid full-Hilbert trial-state expectation* outside R, e.g. the lowest fixed-graph two-triangle eigenvector. By Rayleigh–Ritz the global ground energy \(E_0\le E_{\rm comp}\).

If the predeclared energetic advantage
\[
\boxed{\delta_N=E_R-E_{\rm comp}>0,}
\]
then for **every normalized full-Hamiltonian ground state** \(\psi\),
\[
(P_R H P_R-E_0)P_R\psi=-P_R H P_C\psi.
\]
The spectral gap on P_R is >=delta_N; taking norms gives
\[
\delta_N\|P_R\psi\|\le b_N\|P_C\psi\|.
\]
Hence the **exact conditional quantum probability bound** is
\[
\boxed{
\langle\psi|P_R|\psi\rangle
\le\frac{b_N^2}{b_N^2+\delta_N^2},\qquad
b_N=\frac{\Gamma}{M}\sqrt{N(M+N-6)}.}
\]
This includes coherent mixing with *every non-ring graph*, all fermions and U(1) electric fields: no Born–Oppenheimer elimination, Gibbs approximation or assumed graph selection.

At N=6 and the published regulated energy numbers, b_6=0.0632455532 and delta_6=0.002058490191, so the bound is **P_ring<=0.9989418**. This is **too weak to rule out a ring-dominated ground** at N=6. It must NOT be presented as negative evidence for the six-vertex ring probability; it is a correct, weak quantified obstruction.

For Q=N and fixed Gamma, \(b_N=O(N^{-1/2})\). If there is a fixed \(\delta_*>0\) for an infinite sequence of eligible system sizes, the bound would give \(P_R=O(N^{-1})\). If an *extensive* advantage \(\delta_N\ge cN\) is established, it gives \(P_R=O(N^{-3})\). **Neither energy-gap hypothesis has yet been proved.** They suggest a less costly decisive benchmark: estimate ring and adversarial fragment **fixed-graph energy densities** with controlled errors, including the CP-12 rewiring diagonal and gauge cutoffs, before attempting a full fluctuating N~1000 quantum system.

## Exact combinatorial checks

The independent [check.py](check.py) enumerates every labelled connected N-cycle and every graph reachable by one edge swap for N=4..8, without importing GRUT code, records outside predecessor multiplicities, verifies all identities and calculates the exact scalar norm. Results:
| N | Ring configurations | Max incoming ring predecessors | Scalar cut norm (Gamma=.1) |
|---:|---:|---:|---:|
| 4 | 3 | 2 | .066666667 |
| 5 | 12 | 2 | .067082039 |
| 6 | 60 | 2 | .063245553 |
| 7 | 360 | 2 | .059093684 |
| 8 | 2520 | 2 | .055328334 |

The fully gauge-constrained matter coupling might reduce that norm further. The script does **not** diagonalize CP-12's matter sector, recompute its N=6 block energies, or establish the large-N ring/triangle gap. All limits above are conditional mathematical consequences of the *unchanged* frozen law.

## Discussion with the GRUT III independent reviewer

**Lorentz protection should be an explicit fourth gate.** CP-12's exact gauge and event-permutation symmetries do not establish a common low-energy speed, and its Hamiltonian includes no bare magnetic plaquette term to demonstrate propagating gauge waves. The Collins et al. 2004 percent-level naturalness issue is a serious warning, not a CP-12-specific computed speed splitting or a species-independent 10^-20 bound. See https://doi.org/10.1103/PhysRevLett.93.191301 . Infrared emergent universal speed is possible in certain interacting theories but requires an explicitly derived symmetry/RG basin, e.g. Roy, Juricic, Herbut (2016), https://doi.org/10.1007/JHEP04(2016)018 .

**Computability:** A literal full Q~N, N~1000 fermion/electric/fluctuating-graph evolution is not credible at present. A realistic sequence is: (a) establish exact fixed-graph selection/fragmentation spectral bounds; (b) use prepared bounded-degree graph families as **comparators, not preselected winners**, to compute and cross-validate gauge-invariant charge-based \(R(r,T)\) via fixed-graph local Krylov/linked-cluster/tensor methods at sizes where they are controlled; (c) apply the proved CP-12 graph-preservation error bound \(O(\Gamma^2 T^2/N)\) for graph-preserving detectors to transfer only the **appropriate** prepared-graph results, and use the pair-referenced observable with gauge-invariant sources; (d) compare matter/gauge sector low-k velocities, test species-dependent speed differences, and ultimately compute radiative protection, not just tree-level dispersion. Graph selection and universal velocity protection are separate mandatory gates.

The exact full CP-12 preserved-geometry detector cone is conditional on prepared bounded-degree backgrounds; it does not prove those backgrounds dynamically occur, globally stationary states, universal speeds or Einstein gravity.

**Research ruling:** Work out **selection and spectral exclusions analytically first**, then targeted, fully validated retarded responses on bounded-degree adversarial comparisons; do not declare the original R(r,T) decisive merely because it can be numerically estimated at seven vertices. If the ring/fragment preference remains hostile under controlled growth, terminate that sector early without adding a hand-tuned geometric reward.

No novel empirical physical effect or theory of reality is established by this note.
