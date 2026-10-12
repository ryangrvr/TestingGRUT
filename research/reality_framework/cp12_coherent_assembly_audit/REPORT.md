# CP-12 coherence-enhanced connection formation: exact neutral-graph mechanism and full-Hamiltonian obstruction

**2026-10-12.** Independent additive mathematical investigation rooted at the completed [CP-12 operational follow-up, 65835d0](https://github.com/ryangrvr/TestingGRUT/commit/65835d0c974250314fa2832fd512ea850fac0961). No modification to its frozen Hamiltonian, previous independent branches, canonical GRUT-RAI, Stage-3 Card 3 or a new CP-13. The calculations combine established Johnson graph quantum walks and the already declared CP-12 fermion/gauge terms. No claim of new fundamental physics or literature priority.

## Exact question

The follow-up proved an error \(O(\Gamma^2 T^2/N)\) for bounded graph-preserving detector responses prepared on **one definite geometry**, throughout \(T=o(\sqrt N)\), while N=4,Q=3 dynamics joined a triangle to an isolated vertex at later times. Might coherent initial **disconnected** geometry states alter this formation timescale without changing the law?

**Finding:** Yes for the neutral-only rewiring dynamics and at second order in exact full CP-12; but the full matter/gauge Hamiltonian has rapidly proliferating fermion excitations and the neutral-only finite-time result cannot be promoted to a full formation theorem. Moreover, even the solvable neutral-only limit exhibits revivals and lacks spatial selection.

## 1. Isolate the unchanged graph-rewiring Laplacian

Let \(M=\binom N2\), fix conserved occupied-edge number \(Q=\nu N+O(1)\) with \(\nu>0\), and choose one **physically referenced** vertex \(v\). There are \(S=N-1\) potential connections incident to v and \(A=M-S\) potential nonincident connections. Assume \(Q\le A\). Work with electrically neutral, filled-edge matter basis states for each graph.

CP-12's unchanged positive all-pairs neutral-rewiring term, restricted as a **separate graph-only comparator**, is exactly the Johnson-graph Laplacian
\[
H_J=\frac{\Gamma}{M}\left[Q(M-Q)I-\mathcal A_{J(M,Q)}\right].
\]
The effective comparator omits \(H_{\rm hop}\), the electric interactions and mass dynamics. **The neutral-vacuum collection is not an invariant subspace of full CP-12**, so its finite-time evolution is not CP-12 evolution. The link-rewiring matrix elements themselves are the actual frozen \(-\Gamma/M\) amplitudes.

Define two initial preparations:
- \(|G\rangle\), a single neutral Q-link graph without edges incident on v.
- \(|\Phi_0\rangle = \binom AQ^{-1/2}\sum_{G\subseteq A,\,|E(G)|=Q}|G,\text{neutral}\rangle\), an **explicit coherent superposition** of every neutral graph having v isolated. This requires a selected physical vertex reference and nontrivial graph-sector coherence.

Measure **attachment of v**, i.e. whether one or more occupied edges touch v. This is weaker than the event that the entire graph is connected, and neither is equivalent to a spatial manifold.

## 2. Exact symmetry reduction from \(\binom MQ\) graph states to \(O(N)\)

The graph-only comparator is invariant under permutations of the S incident potential bonds and separately the A nonincident potential bonds. Let \(|k\rangle\) be the normalized sum of all Q-link graphs with exactly k connections incident to v. Their multiplicities are
\[
n_k=\binom Sk\binom A{Q-k},\qquad
0\le k\le\min(S,Q).
\]
\(|\Phi_0\rangle=|0\rangle\). The subspace spanned by \(|k\rangle\) is **exactly invariant** under \(H_J\), not a projection assumed without checking off-sector leakage. The reduced symmetric tridiagonal matrix is
\[
\boxed{
\begin{aligned}
(H_J)_{kk}&=\frac{\Gamma}{M}
 \left[Q(M-Q)-k(S-k)-(Q-k)(A-Q+k)\right],\\
(H_J)_{k,k+1}&=-\frac{\Gamma}{M}
\sqrt{(k+1)(S-k)(Q-k)(A-Q+k+1)}.
\end{aligned}}
\]
The off-diagonal expression follows by counting forward and reverse swaps; the coefficients are determined by the graph law, not fitted. Independent [check.py](check.py) builds the complete Johnson adjacency from scratch for N=4 (20 states), N=5 (210 states), N=6 (3,003 states), and verifies **the entire finite-time attachment probability**, at six times each, to maximum discrepancies under 1e-10.

## 3. Exact coherent-amplitude enhancement at short time

For any single graph with v isolated, exactly \(Q S\) possible one-swap outputs attach v, all mutually orthogonal. Therefore in the graph-only comparator
\[
\boxed{P_{\rm attach}^{\rm definite}(T)
=Q S(\Gamma/M)^2T^2+O(T^3).}
\]

For \(|\Phi_0\rangle\), each output in the k=1 sector has \(A-Q+1\) coherent predecessor graphs. Their amplitudes add before squaring, so the exact result is
\[
\boxed{
P_{\rm attach}^{\rm coherent}(T)
=Q S(A-Q+1)(\Gamma/M)^2T^2+O(T^3).
}
\]
The enhancement factor is **exactly \(A-Q+1\)**, growing as \(O(N^2)\). It is a genuine phase-coherence effect, not additional edges or a changed connection rate.

Because the full frozen CP-12 fermion–gauge hopping Hamiltonian **does not change the graph in one application**, both leading \(T^2\) attachment coefficients are also **exact for the corresponding full-law, neutral-filled-edge initial states**, at each finite N. This does **not** guarantee a uniform-in-N time interval on which the \(T^2\) term dominates.

For Q/N->ν>0,
\[
 P_{\rm attach}^{\rm definite}
 \sim4\nu\Gamma^2T^2/N^2,\quad
 P_{\rm attach}^{\rm coherent}
 \sim2\nu\Gamma^2T^2.
\]
Thus coherence evades the one-definite-graph departure bound's *preparation premise*; it does not falsify that theorem.

For the published N=4,Q=3 **actual assembly preparation**, a half-neutral sea plus half one-edge gauge excitation, the predicted initial connected probability is
\[
P_{\rm conn}(T)
=\frac12[3\times3+2\times3](\Gamma/6)^2T^2+\ldots
=\boxed{\frac1{480}T^2+\ldots}.
\]
The exact full follow-up at T=0.4 reports 0.000320885, close to the leading prediction 0.000333333. At this smallest size there is only one graph with a particular vertex isolated, so the coherent enhancement factor is exactly one; the genuine coherent-versus-definite distinction starts at larger N.

### Numerical trend at frozen \(\Gamma=0.1\), \(Q=N-1\)

| N | One-graph leading \(T^2\) coefficient | Coherent leading \(T^2\) coefficient | Neutral-only \(P_{\rm attach}(T=10)\) |
|---:|---:|---:|---:|
| 4 | .002500 | .002500 | .214226 |
| 5 | .001600 | .004800 | .382415 |
| 6 | .001111 | .006667 | .494971 |
| 12 | .000278 | .012500 | .721464 |
| 40 | .000025 | .017575 | .814675 |
| 80 | .00000625 | .018769 | .828668 |
| 160 | .0000015625 | .019380 | .835025 |
| 320 | .000000390625 | .019689 | .838055 |
| large-N | zero | .020000 | .840990 |

The last column describes the **neutral graph-only quantum walk**, not the full CP-12 time evolution.

## 4. Exact large-N oscillator contraction—and instability of neutral-only assembly

For each fixed k at Q/N->ν>0, the above Jacobi matrix converges entrywise to
\[
\boxed{
 H_\infty=\Gamma\big[
   2\nu+a^\dagger a-\sqrt{2\nu}(a+a^\dagger)\big]
 =\Gamma(a^\dagger-\sqrt{2\nu})(a-\sqrt{2\nu}),
}
\]
on the oscillator number basis \(|k\rangle\). The limiting Jacobi operator is essentially self-adjoint on finite-support sequences (off-diagonal coefficients grow as \(\sqrt{k}\)); the entrywise convergent finite Jacobi matrices give the expected unitary evolution from the vacuum on bounded time intervals. This is a standard bosonic/collective-spin contraction, not a new quantum oscillator.

The initial state is \(|0\rangle\). Under the displaced oscillator its k-distribution is Poisson with mean
\[
\mu(T)=2\nu|1-e^{-i\Gamma T}|^2
      =4\nu(1-\cos\Gamma T).
\]
Therefore the asymptotic quantum-walk attachment probability is
\[
\boxed{
P_{\rm attach}(T)\longrightarrow
1-e^{-4\nu(1-\cos\Gamma T)}.
}
\]
For ν=1 and Γ=.1: \(P_{\rm attach}(10)=0.840990412\), \(P_{\rm attach}(\pi/\Gamma)=0.999664537\), and \(P_{\rm attach}(2\pi/\Gamma)=0\). The finite-N reduced calculations approach these limits; the return at \(2\pi/\Gamma\) becomes sharper as N grows. **This is coherent reversible attachment with near-revival, not stable geometrogenesis or relaxation.**

## 5. Crucial stronger test: vertex attachment is NOT global connectivity

The neutral model's impressive single-vertex attachment must confront actual **full connectedness**, not merely k>=1. For Q=N-1, a globally connected graph is exactly a labelled spanning **tree**. Prüfer's theorem gives the number of such trees with reference vertex degree k:
\[
\#\mathrm{trees}(N,k)=
\binom{N-2}{k-1}(N-1)^{N-k-1},
\quad k\ge1.
\]
Because the quantum state is **uniform inside each |k> orbit**, its *exact* conditional global-connectivity probability is
\[
\boxed{
f_{N,k}
=\frac{\binom{N-2}{k-1}(N-1)^{N-k-1}}
{\binom{N-1}{k}
 \binom{(N-1)(N-2)/2}{N-1-k}},
\quad
P_{\rm conn}(T)=\sum_k |c_k(T)|^2 f_{N,k}.
}
\]
This is exact finite-N probability for the coherent neutral-only graph walk. An independent full enumeration of every Q=N-1 graph through N=7 verified the Prüfer sector counts, including the Cayley sum N^(N-2). It does **not** assume that attachment is connectivity.

| N | \(P_{\rm attach}(T=10)\) | \(P_{\rm globally\ connected}(T=10)\) |
|---:|---:|---:|
| 4 | 0.21422594 | 0.20652398 |
| 5 | 0.38241522 | 0.28773132 |
| 6 | 0.49497073 | 0.27228435 |
| 8 | 0.62131131 | 0.17427045 |
| 12 | 0.72146393 | 0.05184259 |
| 20 | 0.78086067 | 0.00391106 |
| 40 | 0.81467460 | 0.000006430 |
| 80 | 0.82866778 | 2.19032e-11 |
| 160 | 0.83502540 | 3.42751e-22 |
| 320 | 0.83805546 | 1.16168e-43 |

The apparent N-independent **attachment** channel does NOT assemble a connected macroscopic universe. At N=320 the reference vertex attaches with ~84% probability, but global connectedness is about \(10^{-43}\).

### Uniform-in-time **large-system no-assembly bound** for this prepared neutral walk

There is also a mathematical limit, independent of extrapolating the table. The reduced off-diagonal coefficients obey
\[
 |(H_J)_{k,k+1}|
 \le g_N\sqrt{k+1},\qquad
 g_N=\Gamma\sqrt{\frac{QS}{M}}.
\]
Writing \(m(T)=\langle k\rangle_T\), the Heisenberg commutator and Cauchy–Schwarz give
\[
 |\dot m(T)|\le 2g_N\sqrt{m(T)},\qquad m(0)=0.
\]
Hence for **all T at finite N**:
\[
\boxed{\langle k\rangle_T\le
   \Gamma^2\frac{Q(N-1)}{M}T^2.}
\]
This operator bound was tested on every reduced curve at the numerical times.

Conditioned on a reference degree k, the remaining Q-k edges are **uniformly sampled** among the A potential edges joining the other N-1 vertices. For Q/N->ν>0, the number I of vertices isolated within that remaining graph has
\[
E[I\mid k]=(N-1)
\frac{\binom{A-(N-2)}{Q-k}}{\binom A{Q-k}}
\ge c_\nu N,\quad
\operatorname{Var}(I\mid k)=O_\nu(N)
\]
uniformly in allowable k, with \(c_\nu>0\). (The variance follows from the analogous two-vertex hypergeometric probability.) If k is less than \(c_\nu N/2\), the probability of having I<=k is O(1/N) by Chebyshev. If k is greater, Markov and the above moment bound give probability O(Γ²T²/N). Every globally connected graph must have I<=k, since the new edges to the reference vertex can cover at most k otherwise isolated vertices.

Therefore:
\[
\boxed{
P_{\rm conn}^{\rm neutral\ walker}(T)
=O_\nu\left(\frac{1+\Gamma^2T^2}{N}\right),
\qquad
T=o(\sqrt N)
\Longrightarrow P_{\rm conn}(T)\to0.
}
\]

**Scientific consequence:** Within this specific fully coherent initial state and neutral rewiring comparator, there is **NO macroscopic globally connected formation during the growing-time propagation window** of the previous CP-12 operational theorem. Coherence can improve attachment but does not solve spatial selection. The result does **not** extend to the full matter/gauge Hamiltonian, to an alternative physically prepared coherent state, or to arbitrary initial graph classes without fresh proof.

The independent runnable test is [global_connectivity.py](global_connectivity.py).

## 6. Full-CP-12 microscopic matter coupling makes the neutral-only approximation nonuniform

Now restore the exact frozen CP-12 matter/gauge Hamiltonian, with t_hop=1 and both U(1) half-links. On each neutral filled-edge graph, one matter hop moves an edge fermion to either endpoint and creates a distinct gauge-invariant edge hole, vertex charge, and unit half-link electric field. The \(2Q\) outputs are mutually orthogonal and orthogonal to all neutral configurations, for any coherent superposition over different neutral graphs. Thus exactly
\[
\boxed{
\|P_{\rm excited}H_{12}|\Phi_0\rangle\|^2
=2Qt_{\rm hop}^2.
}
\]
In particular, the full-law short-time **matter excitation probability** is
\(2Qt_{\rm hop}^2T^2+O_N(T^3)\), whereas the graph-only/coherent attachment curvature is \(\sim2\nu\Gamma^2T^2\). Their ratio scales as \(N(t_{\rm hop}/\Gamma)^2\), so at the frozen values it diverges with N.

One can express the same issue as an **exact spectral-variance identity** for the coherent neutral state:
\[
\boxed{
\begin{split}
\langle H_{12}\rangle_{\Phi_0}
 &= \kappa Q+\frac{\Gamma Q S}{M},\\
\mathrm{Var}_{\Phi_0}(H_{12})
 &=2Qt_{\rm hop}^2
  +\frac{\Gamma^2}{M^2}QS(A-Q+1).
\end{split}}
\]
The microscopic fermion sector has norm of order \(\sqrt N\) on that state, even though the neutral-only graph attachment scale is O(1). This is a **nonuniform expansion warning**, not a demonstration that the full theory suppresses formation.

**Independent actual-law operator verification:** [full_law_moments.py](full_law_moments.py) constructs \(H_{12}|\Phi_0\rangle\) *directly* from the frozen Gauss-law fermion hopping and neutral-rewiring matrix elements, without constructing the full many-body matrix or importing its solver. It checks Gauss at **every generated output**, keeps the exact CAR fermion signs and incidence orientations, and compares \(\langle H\rangle\), \(\mathrm{Var}(H)\), the attachment coefficient, and the matter-excitation coefficient with the above formulas.

| N | Coherent initial graph states | Distinct \(H_{12}|\Phi_0\rangle\) output states | Exact full variance | Coherent attachment curvature |
|---:|---:|---:|---:|---:|
| 4 | 1 | 16 | 6.0025000000 | 0.0025000000 |
| 5 | 15 | 215 | 8.0048000000 | 0.0048000000 |
| 6 | 252 | 3,822 | 10.0066666667 | 0.0066666667 |
| 7 | 5,005 | 83,083 | 12.0081632653 | 0.0081632653 |

All algebraic moment and Gauss checks passed in a standalone local implementation. This independent full-law result is still a **one-step Hamiltonian/spectral moment computation**, not later-time full coupled quantum evolution. The exact committed script should be rerun verbatim by a partner. It means the favorable neutral oscillator curve has **not** been validated as a finite-time result of the full coupled theory.

The graph-only coherent preparation \(|\Phi_0\rangle\) has reconnection energy \(\Gamma QS/M\sim2\nu\Gamma\), while dephasing it into the equal classical mixture over the same isolated-vertex graphs increases the initial energy by
\[
\boxed{\Delta E_{\rm dephase}
=\frac{\Gamma Q(A-Q)}{M}
\sim \nu\Gamma N.}
\]
This holds for the same **neutral** full-CP-12 states because fermion/electric diagonal expectations do not change. A long-range graph-coherent state is not an innocuous substitute for an incoherent ensemble: its relative phases and referencing of v are physical resources, and its preparation is not supplied by CP-12's small-N assembly result.

## 7. Scientific verdict

**What is earned:** the exact neutral-reconnection dynamics has an analytically solvable **coherence-assisted attachment channel**; the effect has an \(N\)-independent leading coefficient and nontrivial finite-time large-N neutral-only limit. The short-time coefficient and full-law spectral variance follow from frozen CP-12 without modifying its Hamiltonian. This makes the distinction between *formed by dynamics* and *prepared with extensive graph coherence* explicit and quantitative.

**What is NOT earned:** this preparation does not select realistic space, a preferred dimension or even a globally connected graph. Most neutral graphs used in \(|\Phi_0\rangle\) have additional disconnected components and are generally non-manifoldlike. Fixing a particular vertex v presupposes a physical relational reference. Neutral-only revival is not thermalization. And the full coupled fermion–gauge dynamics **cannot** be approximated over its formation timescale by the neutral-only quantum walk on the basis of this calculation.

This does not weaken the owner's prepared-geometry detector cone: that result genuinely holds for one definite graph or an incoherent graph mixture under its specified detector assumptions. It only identifies an inequivalent **coherent** preparation route, with its own physical costs and a much harder large-N full-Hamiltonian verification obligation.

### Recommended next decisive calculation (still using frozen CP-12)

In the *same* \(Q=\nu N\) sequence, begin with \(|\Phi_0\rangle\) and retain **all fermion and electric excitations**. Measure: (i) graph attachment and global connectedness separately; (ii) graph coherence in a physically referenced basis; (iii) charge/flux back-reaction; (iv) occupation of spatially acceptable graph families, not merely k>=1; and (v) source–detector responses as a function of graph distance. Compare at physical times O(1), O(1/Γ), and O(√N)/Γ. The full result must be checked for electric cutoff convergence and independent geometry-space orbit validity. No extra preferred lattice energy or ad hoc metric penalty is authorized.

If connectivity accelerates but graph structure remains random/nonspatial, the fast route has failed its true purpose. If the fermion/gauge dynamics destroys the coherence-assisted channel, the mechanism fails operationally. If spatial graphs emerge with appreciable probability **and** their detector cone remains controlled, that would be the first valid foundation for a stronger continuum investigation—not gravity by itself.

## Reproducibility and scientific priority

[check.py](check.py) is self-contained (NumPy/SciPy), does not import the CP-12 code, constructs the EXACT Johnson graph for N<=6, compares entire time evolution with the exact 1D coherent quotient, and records N up to 320. [full_law_moments.py](full_law_moments.py) independently verifies the full-law Taylor coefficients and spectral moments by explicit one-application state enumeration without requiring SciPy. Independent local equivalent algorithms passed; the verbatim committed script should also be run by a partner. The exact tree-degree formula and global-connectivity probability are independently checked by [global_connectivity.py](global_connectivity.py), including direct enumeration through N=7 and the finite-N expected-reference-degree bound. The graph-only full-sector comparison and the exact full-law derivative/variance identities should not be mislabeled as a full-CP-12 numerical evolution. Standard quantum walks on Johnson graphs, collective spin contractions and coherent-state oscillator revivals have prior art. A publication-level originality search is still needed; no discovery of a new force, gravity or experimental prediction is asserted.
