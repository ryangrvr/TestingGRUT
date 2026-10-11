# CP11 exact mechanism and adversarial large-N analysis

Frozen controls t=m=1, κ=1/2, Γ=1/10. All energies in supplied microscopic units. This is a conventional quantum graph–fermion construction, not derived quantum gravity.

## 1. Fermions generate the geometry energy

For fixed graph G with incidence B_G, the one-particle Hamiltonian is

 D_G = [[m I_N, t B_G], [t B_G^T, -m I_M]].

Its square has blocks m²I+t²B B^T and m²I+t²B^T B. With r=rank B, singular values sqrt(λ_a), it has r eigenvalue pairs ±sqrt(m²+t²λ_a), N-r unpaired +m, and M-r unpaired -m. Here λ_a are the nonzero graph-Laplacian eigenvalues. In the fixed-M fermion sector all M negative-energy levels are filled; m>0 makes this filling unambiguous even when the graph disconnects. Pauli exclusion and the prescribed particle inventory are inputs.

After the common mM offset, the conditional minimum is exactly

 U(G) = κ |E(G)| - Tr[sqrt(m² I_N+t² L_G)-m I_N].

Thus matter motion and geometry binding are determined by the SAME incidence coupling. This is substantially different from prescribing a geometric energy independently of propagation.

## 2. Extensive stability

Concavity of sqrt and Tr L=2E imply

 U(G)/N ≥ κ e - sqrt(m²+2t²e)+m,   e=E/N.

For κ<t²/m the scalar minimum is m-t²/(2κ)-κm²/(2t²). For κ≥t²/m it is zero, at e=0. At the frozen values it is -1/4. The graph kinetic operator is positive, so

 H11 ≥ -N/4.

The empty-graph filled-edge trial state has energy a_N=Γ M/N=Γ(N-1)/2. Therefore -N/4≤E0≤a_N. This proves bounded ground energy per vertex. It does not establish a continuum, finite energy density per physical volume, or renormalizability.

## 3. Sparse low-energy graphs and the actual size of quantum mixing

A stronger result addresses whether the quantum graph term could reverse the conditional spectral comparison extensively.

Let Q count occupied connections and x=<Q> in a ground state. The elementary inequality sqrt(m²+t²λ)-m≤t sqrt(λ), followed by Cauchy–Schwarz, gives the operator bound on each graph block

 H_matter + κQ ≥ κQ - t sqrt(2NQ).

Because the graph kinetic term is positive and E0≤a_N≤ΓN/2, Jensen gives

 κ x - t sqrt(2Nx) ≤ ΓN/2.

Writing x=N y² gives y≤Y where

 Y = [sqrt(2)t + sqrt(2t²+2κΓ)]/(2κ).

Thus x≤C N, C=Y², uniformly at the frozen controls (a deliberately loose bound). This proves only bounded mean degree; it does not impose a maximum degree or Lorentzian causal-link cutoff. These graphs were defined as connection graphs, not causal sets.

For each graph qubit, even when entangled with all fermions,

 |<X_e>|≤2 sqrt(p_e(1-p_e))≤2 sqrt(p_e), p_e=<n_e>.

Hence |Σ_e<X_e>|≤2 sqrt(Mx). With x≤CN and M≤N²/2,

 (Γ/N)|Σ_e<X_e>|≤Γ sqrt(2CN).

Let U_min(N) be the minimum conditional energy over ALL graphs. Variational bounds then give

 U_min+a_N-Γ sqrt(2CN) ≤ E0 ≤ U_min+a_N.

This controls the full quantum Hamiltonian without a Born–Oppenheimer approximation. At fixed Γ and the FROZEN 1/N kinetic normalization, graph mixing changes ground energy only subextensively, O(sqrt N), beyond the common a_N shift. It need not give a small state-vector correction, nor a small finite-N correction. The bound is conservative.

Also <U(G)>-U_min ≤ Γ sqrt(2CN). If every graph in a class A_N has U(G)≥U_min+δN for fixed δ>0, its ground-state probability obeys

 P(A_N) ≤ Γ sqrt(2C)/(δ sqrt N) -> 0.

This exact conditional implication lets an extensive spectral-energy disadvantage exclude a graph class even in the full quantum ground state. It does NOT prove which remaining class dominates.

## 4. Which spectra does the mechanism reward?

For d-regular G, L=dI-A. Let q=m²+t²d and r=t²/q, so r*d<1. The absolutely convergent expansion is

 Tr sqrt(qI-t²A)/N = sqrt(q)[1-Σ_{k≥1} c_k r^k Tr(A^k)/N],

where sqrt(1-z)=1-Σ c_k z^k and c_k>0. Tr A^k counts closed walks. The infinite d-regular tree has no more rooted closed walks of any length than a d-regular graph: projection of walks from the universal covering tree is injective for a fixed root. Therefore

 U(G)/N ≥ e_tree(d),

where e_tree is obtained from the Kesten–McKay adjacency spectral distribution. Graph sequences with locally tree-like neighborhoods approach this lower bound. Finite regular graphs have cycles, hence strict inequality at some order; periodic fixed-dimensional lattices retain a positive density of short cycles and incur a nonzero limiting penalty whenever compared at the same degree.

This is a mechanism for suppressing local loops, not for selecting Euclidean spatial geometry. High girth alone does not imply expansion. The conclusion is that nongeometric, locally tree-like competitors are admitted at lower energy, not that every minimizer is a proved expander.

At frozen controls, asymptotic conditional energies per vertex:

- ring / one-dimensional chain bulk (degree 2): -0.177609971862
- square lattice (degree 4): -0.184539293475
- cubic lattice (degree 6): -0.100067375815
- honeycomb lattice (degree 3): -0.197229928022
- locally tree-like degree 3: -0.197366056730
- locally tree-like degree 4: -0.186535041049
- locally tree-like degree 6: -0.101910317859

The degree-three nongeometric sequence also beats all listed lattices. This is not an exhaustive exclusion of every manifoldlike graph or weighted effective geometry. For each listed lattice sequence its positive energy-density difference from an available competing sequence implies an asymptotic suppression result by section 3. The honeycomb difference is only 0.000136128708, so the rigorous finite-N probability bound becomes useful only at enormous N; do not use it as a practical crossover estimate.

## 5. Propagation earned and propagation not earned

For a fixed graph, D²=m²+t²L in the vertex block. On a SUPPLIED cubic spatial lattice of spacing a,

 ω²=m²+4t² Σ_α sin²(p_α a/2)
     =m²+t²a²|p|²+O(a⁴p⁴).

This gives a conditional relativistic-looking dispersion with limiting coefficient c=t a, but does not select a, dimension, Lorentz invariance of interactions or physical spin. Incidence fermions are not automatically Standard Model Dirac spinors.

On a bounded-degree connection graph, the one-particle hopping has finite-range matrix support, so the Taylor series for propagation begins at the appropriate graph distance (two incidence hops per vertex edge). This operational distinction was absent for CP10's complete comparability graph. It defines graph locality, not a recovered physical light cone; high-girth expanding graphs can have logarithmic diameter and no fixed-dimensional spatial volume growth.

The full theory still permits allpairs creation of new connections and uses external time. Conditional graph locality is insufficient to establish dynamical spacetime causality.

## 6. What has and has not been established

Earned: one Hamiltonian yields matter-induced binding, connected small-N ground distributions, a derived spectral geometry energy, an extensive lower bound and bounded mean connection degree in the ground state. Its conditional propagation distinguishes graph separation.

Failure: the same spectral functional systematically favors locally tree-like competitors over the tested regular lattices. At the frozen kinetic scaling, quantum graph mixing cannot reverse a fixed extensive energy penalty asymptotically. There is no derived physical spacetime, graviton, universal gravitational coupling, measurement outcome law or experimentally new consequence.

This failure will not be repaired by changing κ, Γ, adding a loop reward, switching the particle filling or calling a lattice a vacuum. Such changes require a distinct accounted hypothesis. The useful unresolved opportunity is a microscopic interaction with an independently motivated preference for the loop structure that fermionic free binding alone suppresses. Naming gauge flux, interactions or topology does not establish that preference.

## Sources and priority

- Konopka, Markopoulou, Severini (2008), Quantum Graphity, https://arxiv.org/abs/0801.0861 : dynamical graphs and emergent-locality program.
- Spector and Schwartz (2018), https://arxiv.org/abs/1808.05632 : nonlattice degeneracies of simple graph Hamiltonians.
- Bianconi (2021), https://arxiv.org/abs/2106.02929 : topological Dirac dynamics on vertices and edges, established square-root-Laplacian framework.
- Sarid (2023), https://doi.org/10.1002/rsa.21150 : random-regular spectral results and stated Kesten–McKay density.
- Benincasa and Dowker (2010), https://arxiv.org/abs/1001.2725, and Aslanbeigi et al. (2014), https://doi.org/10.1007/JHEP06(2014)024 : causal-set continuum wave operators and stability qualifications. These do not supply a selector for this model and are not imported as a proof of its continuum.

The model-specific derivations above are author calculations pending independent review. No claim of mathematical priority follows from reproducing these standard spectral and variational methods.

## 7. Finite temperature exposes a second, independent failure

Ground-state sparsity does not imply a sparse thermal phase. At any fixed positive temperature T=1/β, let Z_N be the canonical partition function in the frozen M-fermion sector, including ALL 2^M graph configurations. Its fermionic dimension is d_F=binomial(N+M,M)=binomial(N+M,N), with log d_F=O(N log N).

Separate H=κ Q⊗I_F+R. The mass term after the common offset is 2m N_vertex, with norm≤2mN. The hopping norm in every graph block is ≤t Σ_a sqrt(λ_a)≤t N sqrt(N-1). The graph kinetic norm is Γ(N-1). Thus uniformly

 ||R|| ≤ 2mN + tN sqrt(N-1) + Γ(N-1) = o(N²).

Eigenvalue monotonicity gives

 |log Z_N - [log d_F + M log(1+exp(-βκ))]| ≤ β ||R||.

Consequently, for fixed β>0 and frozen finite couplings,

 lim (log Z_N)/M = log(1+exp(-βκ)).

Convexity in κ and differentiability of the limiting function imply

 lim <Q>_β/M = 1/(1+exp(βκ)) > 0.

The conclusion uses the full quantum Hamiltonian; it does not require commuting graph and fermion operators. Every fixed nonzero temperature has an asymptotically dense MEAN connection count. The density can be tiny at low T but is independent of N. The T->0 and N->infinity limits therefore do not commute in the desired sparse sense. This is not a statement that every individual thermal measurement has that density; concentration would require an additional argument.

The conditional ground-energy and full zero-temperature successes do not extend to a stable sparse finite-temperature thermodynamic limit under the supplied counting measure. A sparse thermal limit would need an N-dependent temperature or another physical mechanism. Neither has been earned here. No temperature scaling is adopted as a repair.

This closes two separate tests of this candidate: its ground energy fails to favor the listed local lattice geometries, and its fixed-positive-temperature ensemble has dense mean connectivity. It remains a mathematically defined finite quantum model, with a genuinely derived matter–geometry coupling, but fails the intended route to physical spatial geometry.
