# GRUT constructive physics: conserving bonds is not enough to generate space

Independent theorem and controls, 2026-10-11. **This is NOT CP-12 or a repaired CP-11 law.** Separate branch from frozen [CP-11](https://github.com/ryangrvr/TestingGRUT/blob/bbbe393ca36a449c45ed59e02c5a385857aa1b7f/research/constructive_physics/attempt11/LAW.md). All earlier packets, main GRUT-RAI, and Stage-3 Card 3 are untouched. The graph-counting ingredients are known mathematics. This is a proposed model-specific synthesis and falsification criterion, not independently established novel fundamental physics.

## 1. Exact finite-temperature quantum theorem

Let N labelled vertices have M=N(N-1)/2 possible edges. Impose a **new, expressly supplied** conservation rule fixing occupied bonds to K_N=alpha N+O(1), alpha>1. Thus only simple K-edge graphs are configurations; every such labelled graph is included and receives one orthogonal geometry basis state. Let every graph carry a matter sector of the same dimension d_F(N), as in CP-11's supplied N-vertex plus M-edge fermion modes at fixed total fermion count M.

Assume the FULL Hermitian Hamiltonian on this joint space can be written
\[
H_N=E_{\rm ref}(N)I+R_N,\qquad\|R_N\|\le CN
\]
for N-independent C. No diagonal/commutation assumption: R_N may coherently mix geometry and matter. At any fixed inverse temperature beta>0, with rho=e^{-beta H}/Tr(e^{-beta H}), let A_N be a prescribed graph class with joint projector Pi_A.

**THEOREM (quantum Gibbs class-probability bound):**
\[
\boxed{\operatorname{Tr}(\Pi_A\rho_\beta)
\le \frac{|A_N|}{\binom{M}{K_N}}\exp(2\beta C N).}
\]

**Proof:** Since each eigenvalue of R_N lies in [-CN,+CN], each eigenvalue of e^{-beta H_N} lies between e^{-beta(E_ref+CN)} and e^{-beta(E_ref-CN)}. Therefore Tr(Pi_A e^{-beta H}) <= |A| d_F e^{-beta E_ref+beta CN}, whereas Z >= choose(M,K)d_F e^{-beta E_ref-beta CN}. Divide. The argument uses scalar spectral bounds, not an invalid general operator-monotonicity assertion about exponentials. QED.

If the target graph class has at most exp(O(N)) unlabelled spatial templates, then |A_N|<=N!exp(O(N)). Meanwhile Stirling gives
\[
\log\binom{M}{K_N}=\alpha N\log N+O(N),\quad
\log|A_N|\le N\log N+O(N).
\]
Consequently
\[
\boxed{\log \Pr_\beta(A_N)\le-(\alpha-1)N\log N+O_\beta(N)\to-\infty.}
\]

**Result:** At every fixed positive temperature, extensive energy differences and extensive-norm quantum bond fluctuations cannot give thermodynamic dominance to this thin class of spatial graphs inside the full uniform K-edge Hilbert space. This does **not** forbid geometrical ground states at T=0, metastable nonergodic states, special constrained kinematics, graph ensembles outside the thin-template assumption, nonextensive interactions, or a physically derived alternative measure.

The count condition covers any finite collection of regular finite-dimensional lattice templates. Published graph enumeration also establishes N!exp(O(N)) growth for *all graphs embeddable in a fixed-genus surface*, not just a single lattice (Giménez–Noy; Chapuy et al.). This supplies a broad and independently studied geometric example, not a classification of all possible emergent spacetime geometries.

## 2. A complete conserved-bond quantum counterfactual with CP-11 matter

The actual frozen CP-11 graph flip \((\Gamma/N)\sum_e(I-X_e)\) changes K, so it cannot be quietly described as number-conserving. To test a conservation repair, declare instead a **positive fixed-K bond-swap Laplacian**:
\[
H_{\rm swap}=\frac{\Gamma}{N^2}
\sum_{\{G,G'\}:G\sim G'}
(|G\rangle-|G'\rangle)(\langle G|-\langle G'|)\otimes I_F,
\]
where neighbors differ by replacing one active edge with one inactive edge. Each undirected swap is included exactly once. This conserves K, is Hermitian, permutation covariant and nonnegative. Its underlying Johnson graph has degree K(M-K), hence
\[
\|H_{\rm swap}\|\le 2\Gamma K(M-K)/N^2=O(N)
\]
for K=alpha N. **Swaps are still globally nonlocal**, so this Hamiltonian is only an adversarial thermodynamic control, not a claimed physical microscopic mechanism.

Retain CP-11's original one-body Dirac incidence coupling, m,t fixed and canonical total fermion number M. Its common +mM offset converts the mass term to \(2m N_{\rm vertex}\), with norm <=2mN. For any graph B_G with K active edges, \(\operatorname{Tr} B_G B_G^\dagger=2K\), so singular-value Cauchy–Schwarz gives \(\sum_a s_a\le\sqrt{2NK}\). Its second-quantized hopping operator norm is <=t sum s_a. Therefore
\[
\boxed{\|H_N-\kappa K I\|\le 2mN+t\sqrt{2NK}+
       2\Gamma K(M-K)/N^2.}
\]
For K=alpha N, this supplies \(C=2m+t\sqrt{2\alpha}+\Gamma\alpha+o(1)\). Theorem 1 therefore controls the **fully noncommuting graph–fermion Gibbs state**, not just a classical Born–Oppenheimer conditional energy.

At fixed K the original CP-11 edge price kappa K becomes a common energy constant: it does not choose graph shape. Giving occupied bonds arbitrary unitary magnetic-incidence phases also preserves \(\operatorname{Tr} B_\phi B_\phi^\dagger=2K\), so the fermionic uniform bound survives. A new extensive gauge Hamiltonian or physically nonuniform graph-state measure would be an additional hypothesis requiring separate tests.

## 3. Fixed degree is also insufficient without more structure

Even restricting all vertices to degree d>=3 leaves \(\exp[(d/2)N\log N+O(N)]\) labelled d-regular graphs (the established fixed-degree enumeration). A thin spatial family with at most exp(O(N)) unlabelled templates has at most exp[N logN+O(N)] labels. Thus any uniformly extensive-norm Hamiltonian on the *full d-regular graph sector* suffers an analogous equilibrium probability bound:
\[
P_\beta(\mathrm{spatial\ templates})\lesssim
\exp[-(d/2-1)N\log N+O_\beta(N)].
\]
A four-valent port constraint therefore does not, by itself, select a square lattice. Any additional graph-space constraint must be derived and counted rather than smuggled in as geometry selection.

An independent uniform G(N,K) structural control is the exact expected number of simple ell-cycles:
\[
\mathbb E C_\ell=\frac{(N)_\ell}{2\ell}\frac{(K)_\ell}{(M)_\ell}
\longrightarrow\frac{(2\alpha)^\ell}{2\ell}.
\]
For K=2N and ell=4 this tends to **32**, independent of N, while a 2D square LxL torus has at least N elementary 4-cycles. This counting is for the uniform measure, not an interacting thermal cycle theorem.

## 4. Executable verification

Independent [check.py](check.py) reconstructs CP-11 canonical fermions **without importing GRUT code**: for every four-vertex graph with K=4 (15 total) and five-vertex graph with K=6 (210 total), it diagonalizes the exact one-body Dirac Hamiltonian and determines the canonical M-fermion ground and ceiling energies by filling the M lowest and highest single-particle levels. It also builds the full Johnson bond-swap Laplacian and verifies positivity, its zero uniform mode and the rigorous norm bounds. It checks simple four-cycle counting by exhaustive finite enumeration, plus large-N log counts through 1,000 vertices. Equivalent standalone controls were run before publication; execution of the exact committed script is a separate reproducibility task.

| Control | Exact maximum fermion energy magnitude | Rigorous bound |
|---|---:|---:|
| N=4, K=4 | 10.700169593 | 13.656854249 |
| N=5, K=6 | 13.913591358 | 17.745966692 |

For a **single unlabelled square-lattice template orbit**, the loose *uniform-ensemble* upper bound N!/choose(M,2N) is:

| N | log10 upper bound for orbit probability |
|---:|---:|
| 16 | -15.8163 |
| 25 | -32.3025 |
| 36 | -54.3623 |
| 100 | -204.2843 |
| 1,000 | -3092.2053 |

Those figures are raw counting at beta=0, not actual Gibbs probabilities under interacting CP-11. The Gibbs upper bound contains the large e^(2 beta C N) factor and can be numerically vacuous at accessible N. The analytic limit is what establishes fixed-beta suppression.

## 5. What this does and does not mean for GRUT

**Exact outcome:** The tempting successor idea, “conserve O(N) link carriers and let free fermions plus gauge flux select space,” does **not**, by those ingredients alone, ensure a stable equilibrium spatial phase. In the stipulated full fixed-inventory graph Hilbert space, extensive-norm dynamics cannot thermally select a thin target class for alpha>1 at fixed beta.

**Physical escape conditions:** At least one premise must change for that specific type of geometry to dominate: (a) the physical accessible topology/state-counting measure must be derived with enough correlated constraints to eliminate the O(N logN) surplus; (b) a physically justified energy/free-energy contrast O(N logN) or more must emerge, without arbitrarily inserting an N-dependent price; (c) the relevant state is not the unrestricted fixed-temperature Gibbs equilibrium; or (d) the actual spatial phase has far more configurations than the thin-template count. These alternatives are research hypotheses, not repairs achieved here.

Naively quotienting graph labels also does not magically fix the *count*: the number of unlabelled K-edge graphs is at least choose(M,K)/N!, since an isomorphism class has at most N! labellings. If the desired target contains only exp(O(N)) unlabelled templates, its fraction among uniformly counted unlabelled classes has the same large suppression. Actual gauge-constrained matter/reference Hilbert spaces can have different weights and must be treated explicitly.

Even a successful solution of the state-counting problem must separately earn operational light cones, Lorentzian propagation, Einstein dynamics and universal gravitational coupling. **No theory of reality has been completed by this theorem.**

## Relevant, already published work

- [Giménez & Noy, labelled planar graph enumeration (2009)](https://arxiv.org/abs/math/0501269).
- [Chapuy et al., fixed-genus graph enumeration (2011)](https://arxiv.org/abs/1001.3628).
- [Bender & Canfield, prescribed-degree graph enumeration (1978)](https://doi.org/10.1016/0097-3165(78)90059-6).
- [Konopka, Markopoulou & Severini, Quantum Graphity (2008)](https://arxiv.org/abs/0801.0861).
- [Chen & Plotkin, emergent manifold graph statistical mechanics (2013)](https://doi.org/10.1103/PhysRevD.87.084011), especially relevant because it **already studies a fixed number of vertices and edges** with additional valence and rotational-symmetry energy terms and reports finite-size two-dimensional-like structures.

The combinatorial and statistical principles are familiar. An equation-by-equation novelty assessment of the combined CP-11-conditional theorem is still required before claiming a publication-level original discovery. The independent mathematical contribution here is a rigorous **model-specific falsification criterion**, not a novel physical constant or experimental prediction.
