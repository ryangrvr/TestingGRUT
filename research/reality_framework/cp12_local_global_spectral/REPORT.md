# CP-12 partner physics audit: local degree/motifs do NOT identify macroscopic spatial geometry

**2026-10-12.** Independent, additive constructive/adversarial research branching from verified owner [5bbac5491f4ec506fe55050db796ba0faa9a8d38](https://github.com/ryangrvr/TestingGRUT/commit/5bbac5491f4ec506fe55050db796ba0faa9a8d38). This packet does **not** edit CP-12, GRUT-RAI, Stage-3 Card 3, or previous frozen models. It uses exact graph/Laplacian mathematics and finite random-regular graph samples **as comparators only**, not CP-12 many-body quantum evolution. It is not a fundamental-law discovery. The graph ideas have existing mathematical prior art.

## 1. Why the owner’s result matters

CP-12's full matter/electric theory produces a positive fourth-order conditional degree term: \(E_{\mathrm{degree}}=\lambda\sum_v\binom{d_v}{2}\) with \(\lambda=2t_{\rm hop}^4/\Delta^3\) in a **weak-hopping, fixed-graph** approximation. The owner verified the exact classical degree-switching bound and the exact 2-regular cycle-partition entropy: regularity has a finite-temperature defect obstruction, and even perfect 2-regularity does not select one spanning ring. These remain classical counterfactual statements; this packet does not extend them to the full interacting quantum Gibbs state.

The obvious escape is to prefer a larger mean degree and add rules that favor short spatial loops. This packet tests that idea more aggressively, with **identical connected graphs that even agree on every rooted neighborhood of a chosen fixed radius**.

## 2. Exact local indistinguishability theorem

Let \(C_a\square C_b\) be the Cartesian product of cycles with a,b>=5, so that it is an undirected simple 4-regular toroidal grid. For any fixed nonnegative integer r choose a fixed width w>2r+2, and any k large enough that L=wk>2r+2. Compare graphs on the **same \(N=L^2\) vertices**:
\[
G_N=C_L\square C_L,\qquad
H_N=C_w\square C_{N/w}
     = C_w\square C_{wk^2}.
\]
They have, exactly and at every N:
- 4 edges incident to each vertex; \(Q=2N\).
- \(W=\sum_i\binom{d_i}{2}=6N\), so *identical CP-12 fourth-order degree energy* in its explicit domain of applicability.
- One connected component each.
- No triangles and, for w and both side lengths >=5, **exactly N four-cycles**.
- Every rooted radius-r graph neighborhood is isomorphic; this follows because cycles shorter than \(2r+2\) cannot wrap within an r-ball and both products locally agree with the infinite square lattice.

Therefore **any classical diagonal geometry energy of the form** \(E(G)=\sum_{v\in V(G)} f(B_r(G,v))\), where f depends only on an isomorphism-invariant rooted r-ball, gives *exactly equal values* on \(G_N\) and \(H_N\), for every N. This also covers arbitrary fixed linear combinations of degree, short-cycle and other bounded-radius motif terms.

**This is a strict conditional non-selection theorem for bounded-radius diagonal geometry energies. It is NOT a no-go theorem for the full CP-12 Hamiltonian**, which has noncommuting geometry rewiring, dynamical matter/gauge modes and higher-order topology-sensitive correlations. Quantum spectra of local Hamiltonians can distinguish globally different graphs despite identical local geometric motifs.

## 3. Exact propagation spectrum distinguishes them

For the ordinary (unweighted) graph Laplacian on \(C_a\square C_b\),
\[
\lambda_{p,q}
=4\sin^2(\pi p/a)+4\sin^2(\pi q/b).
\]
Hence for a,b>=5,
\[
\boxed{\lambda_2(C_a\square C_b)
=4\min\{\sin^2(\pi/a),\sin^2(\pi/b)\}.}
\]
Here \(\lambda_2\) is the first nonzero Laplacian eigenvalue, not CP-12's physical many-body excitation gap. Taking \(N=(wk)^2\) and fixed w gives:
\[
\boxed{\lambda_2(G_N)\sim 4\pi^2/N,\qquad
\lambda_2(H_N)\sim4\pi^2w^2/N^2.}
\]

For numerical tests choose **w=9,r=3** and \(k=2,3,6\). Fully constructed network graphs agree on every radius-three neighborhood, component count, W, number of triangles and 4-cycles, while the independently computed exact eigenvalue formula gives:

| N | Square side | Thin rectangle | \(\lambda_2(G_N)\) | \(\lambda_2(H_N)\) | Gap ratio |
|---:|---|---|---:|---:|---:|
| 324 | 18×18 | 9×36 | 0.1206147584 | 0.0303844940 | 3.970 |
| 729 | 27×27 | 9×81 | 0.0539102588 | 0.0060141177 | 8.964 |
| 2,916 | 54×54 | 9×324 | 0.0135232845 | 0.0003760591 | 35.961 |

A genuine continuum candidate must therefore make the **long-wavelength response** of its matter/gauge fields depend on something beyond the local invariants they share.

## 4. An operational diffusion observable detects the dimension crossover

Use the same graph Laplacian's normalized heat return \(P_N(s)=N^{-1}\mathrm{Tr}(e^{-sL})\), and \(d_s(s)=-2\,d\ln P_N(s)/d\ln s\). For a Cartesian torus, the heat return factorizes exactly:
\[
P_{a,b}(s)=
 \Big[\frac1a\sum_{p=0}^{a-1}e^{-4s\sin^2(\pi p/a)}\Big]
 \Big[\frac1b\sum_{q=0}^{b-1}e^{-4s\sin^2(\pi q/b)}\Big].
\]
This gives the *exact finite-size* local spectral dimension
\[
d_s(s)=2s\left(
 \frac{\sum_p\lambda^{(a)}_p e^{-s\lambda^{(a)}_p}}
      {\sum_p e^{-s\lambda^{(a)}_p}}+
 \frac{\sum_q\lambda^{(b)}_q e^{-s\lambda^{(b)}_q}}
      {\sum_q e^{-s\lambda^{(b)}_q}}
\right).
\]
When \(1\ll s\ll N\) in the square grid, both axes support diffusion, so \(d_s\simeq2\). In the thin strip with fixed w and \(w^2\ll s\ll(N/w)^2\), the transverse axis saturates while the longitudinal axis diffuses, so \(d_s\simeq1\). (Numerical crossover starts earlier than the loose w² bound, controlled by the first transverse Laplacian gap.)

At the declared time **s=20**:

| N | Square spectral dimension | Thin spectral dimension |
|---:|---:|---:|
| 324 | 1.47175255 | 1.00962995 |
| 729 | 2.00365719 | 1.00964085 |
| 2,916 | **2.01282539** | **1.00964085** |

The N=324 square is already finite-size affected at time20. The N=729/2916 cases show the clean 2D-vs-1D crossover. At extremely late diffusion time on *any finite graph*, \(d_s\to0\) because the heat flow equilibrates; no asymptotic dimension is inferred from that limit.

The thin graph is mathematically a 2D torus at every finite N; **its fixed circumference collapses its large-scale effective dimension**. Thus even *topological torus structure* and local 2D appearance do not guarantee extended macroscopic 2D space.

## 5. Hostile connected high-valence control: random 4-regular graphs

Exact degree selection can also produce **connected expanders**, whose Laplacian gap does not vanish in the thermodynamic limit. For random fixed-degree 4-regular graphs the established Alon–Friedman spectral theorem gives (with high probability, up to arbitrarily small epsilon)
\[
\lambda_2(L)\ge4-2\sqrt3-\epsilon\simeq0.535898-\epsilon,
\]
in sharp contrast to the \(1/N\) grid gap. This is a known theorem, **not an original GRUT claim**. Source: Joel Friedman, [A Proof of Alon's Second Eigenvalue Conjecture and Related Problems](https://arxiv.org/abs/cs/0405020) and [AMS Memoir](https://doi.org/10.1090/memo/0910).

Three independently seeded NetworkX 4-regular samples per N give:

| N | Random-graph \(\lambda_2\) range | Random 4-cycle counts | Grid 4-cycle count |
|---:|---:|---:|---:|
| 324 | 0.541595–0.598334 | 7–11 | 324 |
| 729 | 0.529561–0.571999 | 10–11 | 729 |

All these sampled graphs are connected and have the exact same N,Q,W as the balanced grid. They generally *do not* share its short loops; their role is an additional adverse control showing why regularity and connectedness do not select geometric propagation.

A sparse eigensolver using only k=3 smallest-magnitude eigenvalues initially **missed degenerate low eigenmodes** in a different 4-regular clique-necklace diagnostic; that preliminary result was discarded. The published random-graph script requests k=10, checks residuals and never relies on sparse diagonalization to compute the torus gaps, which are analytic.

## 6. A frozen-CP-12 *physical* diagnostic distinguishes these otherwise locally identical geometries

The graph Laplacian above is a geometry diagnostic; it is not automatically CP-12's fermion/gauge propagation operator. However **CP-12 itself supplies a different conditional, gauge-respecting distance-energy measurement**: for a charge at vertex v and a hole on occupied edge e={u,w}, the least bare electric flux energy is
\[
E_{\min}(v,e;G)=\kappa Q+2m+\frac{\epsilon_E}{2}d_{\rm inc}^G(v,e),
\quad
\boxed{d_{\rm inc}^G(v,e)=1+2\min[d_G(v,u),d_G(v,w)]}.
\]
This follows from the bipartite incidence path and the owner's frozen-law Gauss constraints; it does NOT account for flux-screening or interacting-dressed charge energies.

For a Cartesian cycle graph \(C_a\square C_b\), average uniformly over vertices v and occupied edges e. Let \(D_v(\ell)\) be the mean cycle distance from a fixed vertex to a uniform vertex and \(D_e(\ell)\) the mean minimum distance from the endpoints of a fixed edge to a uniform vertex. An elementary exact counting identity for every cycle length \(\ell\ge5\) is
\[
D_v(\ell)+D_e(\ell)=\frac{\ell-1}{2}.
\]
Half the Cartesian-product edges run along each axis. Inserting the incidence distance formula therefore gives the **exact finite-N identity**
\[
\boxed{
\langle d_{\rm inc}\rangle_{v,e;C_a\square C_b}
=\frac{a+b}{2}.
}
\]
So CP-12's average bare electric-field excess for a *single fixed charge/hole pair* over uniform placements is
\[
\boxed{\langle E_{\rm electric}\rangle_{\min}
=\frac{\epsilon_E}{4}(a+b).}
\]

At frozen \(\epsilon_E=1/4\):

| N | Square mean incidence distance | Thin mean distance | Square bare electric excess | Thin bare electric excess |
|---:|---:|---:|---:|---:|
| 324 | 18 | 22.5 | 2.2500 | 2.8125 |
| 729 | 27 | 45 | 3.3750 | 5.6250 |
| 2,916 | 54 | 166.5 | 6.7500 | 20.8125 |

The exact values were computed independently by enumerating all vertex/edge orientations on both cycle factors, rather than inferred from the Laplacian. They demonstrate a **genuinely model-grounded physical energy diagnostic** distinguishing global aspect ratio despite identical local motif data.

**Boundary:** averaging *minima* of charged conditional configurations does not imply the interacting CP-12 ground state prefers the square, because its vacuum, kinetic, gauge and charged-sector state multiplicities also determine the free energy. Indeed cheap charge separations might favor a random expander with short paths over any finite-dimensional manifold. A dynamical selection calculation would have to compare the full dressed energies/free energies and geometrical probabilities, not merely this bare electric excess.

## 7. Cooling degree-four selection also fails an extensive spatial-family counting test

The owner's earlier two-regular theorem shows that temperature reduction does not choose a spanning ring from equally weighted cycle partitions. There is a stronger higher-valence version using established **Bender–Canfield fixed-degree graph enumeration**, without claiming an original counting theorem.

For N vertices and Q=2N simple undirected bonds, the degree-only energy \(\lambda W\), \(\lambda>0\), is minimized by **all** four-regular graphs, each with exactly \(W=6N\). The conditional distribution given perfect degree four is therefore uniform for *any* inverse-temperature sequence \(\beta_N\lambda_N\), including an idealized zero-temperature projection. The standard labelled count satisfies
\[
\boxed{
R_{N,4}\sim
\frac{(4N)!}{2^{2N}(2N)!(4!)^N}e^{-15/4}
=\exp[2N\log N+O(N)].
}
\]
If a specified spatial class contains at most \(e^{cN}\) unlabelled templates and counts every labelling of each, then \(|\mathcal S_N|\le N!e^{cN}\), while every member is four-regular. Therefore
\[
\boxed{
P_{\rm degree\ only}(\mathcal S_N)
\le P(\mathcal S_N\mid\mathrm{4\!-\!regular})
\le \frac{N!e^{cN}}{R_{N,4}}
=\exp[-N\log N+O(N)]\to0.
}
\]
The same asymptotic suppression persists if we *condition on connectedness*, because fixed d>=3 random regular graphs are connected with probability approaching one (Wormald's classical result). For one **single** square-torus isomorphism orbit, its number of labelled realizations is at most N!, giving the stronger explicit upper bound \(N!/R_{N,4}\).

This is only a **degree-energy Gibbs theorem** at Q=2N. Exact four-cycle density is a much stronger, separately imposed condition and could alter entropy at leading order; the local-indistinguishability counterexample above proves a distinct point, that bounded-radius diagnostics do not guarantee the macroscopic geometry even when some of them are fixed exactly. The full quantum CP-12 dynamics is not replaced by uniform labelled graph counting.

Prior literature: [Wormald, Models of Random Regular Graphs](https://web.williams.edu/Mathematics/sjmiller/public_html/ntprob19/handouts/graphs/Womald_ModelsRandomGraphs.pdf), Theorem 2.2 and §2; [Wormald, asymptotic connectivity of labelled regular graphs (1981)](https://doi.org/10.1016/S0095-8956(81)80021-4); [Friedman, Alon's second eigenvalue conjecture](https://arxiv.org/abs/cs/0405020).

## 8. The physics consequence

There are **three logically separate properties** that CP-12 must obtain, not infer from one another:

1. **Degree organization**: CP-12 has a genuine conditional Pauli-based candidate at fourth hopping order, with a restricted perturbative domain.
2. **Large connected structure**: the owner's sixth-site results show that identical degrees still allow lower-energy disconnected short cycles. Global structure has to beat component entropy and dynamics.
3. **Macroscopic dimension and causal propagation**: this packet proves that even exact degree, perfect connectedness, triangles, 4-cycles and arbitrarily prescribed finite local neighborhoods can match while the Laplacian's long-wave transport is qualitatively different.

A fourth obligation is that CP-12's physical fermion/gauge excitation propagator must genuinely reduce to an appropriate geometric operator; an arbitrary graph Laplacian is only a **geometric diagnostic**, not a derivation of the physical signal speed.

This gives a concrete **unfitted** large-N falsifier for frozen CP-12:

- At the same Q/N, measure *in the same full quantum state* the degree distribution, short loops, graph component sizes, graph ball growth, spectral return \(P_G(s)\) on sampled observed graph configurations, and a relational matter/charge detector response as a function of graph distance.
- Compare connected grids, thin grids, clique-like fragmentation and random regular expanders, with identical N,Q and degree statistics; do not let an arbitrarily selected graph be treated as dynamically preferred.
- If only local regularity/plaquettes improve while the physical detector sees expander-like or compactified geometry, **the spacetime-selection claim fails**.
- If a smooth manifoldlike sector develops appreciable physical probability *and* its dynamical response has a controlled Lorentzian continuum, proceed to gravitational degrees of freedom. Neither has been established.

**No new CP-12 law, metric penalty, flux preference, or geometry selection assumption was introduced here.** This is a mathematical counterexample that forces the *full* microscopic dynamics to demonstrate global propagation selection, rather than accept short local motifs as its proxy.

## Reproduction

[check.py](check.py), Python with NumPy/SciPy/NetworkX, builds the exact competing graphs at N=324,729,2916, verifies graph isomorphism of radius-3 neighborhoods and exact motif equality, computes all exact torus Laplacian gaps and heat-return dimensions, enumerates the **CP-12 bare charge/hole incidence-distance identity**, and supplies seeded 4-regular adversarial samples with checked sparse eigenpair residuals. It writes [results.json](results.json) when run. An equivalent independent local script with these same formulas passed; **verbatim execution of the committed script remains a separate audit step**.
