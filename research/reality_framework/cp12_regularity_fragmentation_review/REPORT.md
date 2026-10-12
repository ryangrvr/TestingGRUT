# CP-12 independent review: degree selection is not connected geometry

**2026-10-12.** Additive independent analysis forked from frozen full-coupled CP-12 follow-up [7be83fc](https://github.com/ryangrvr/TestingGRUT/commit/7be83fc1f3d6a5799f11a587d3df5abc023dda4f). No modifications to the CP-12 Hamiltonian, original follow-up, canonical GRUT-RAI, Stage-3 Card 3 or earlier scouts. The calculations are **two precise classical Gibbs surrogate theorems** that expose what the author's emergent degree interaction alone does not earn. They are NOT the exact CP-12 joint fermion/gauge/geometry thermal state, not an explanation of real spacetime, and their graph enumeration has extensive prior art.

## 0. Owner result respected, not silently extended

CP-12's frozen matter+U(1) gauge hopping has exact neutral-vacuum fourth moment
\[
\langle V^4\rangle=12Q^2-8Q-4W(G),\qquad
W(G)=\sum_i\binom{d_i}{2}.
\]
In the **declared bare fixed-graph, weak-hopping comparator only**, with \(\Delta=2m+\epsilon_E/2\),
\[
E_G=\kappa Q-\frac{2Qt_{\rm hop}^2}{\Delta}
+\left[4Q+2W(G)\right]\frac{t_{\rm hop}^4}{\Delta^3}
+O_G(t_{\rm hop}^6).
\]
It is a physically derived positive coefficient of degree irregularity. The frozen \(t_{\rm hop}=1\) lies outside the proven small-hopping perturbative limit, and the *full graph compression includes the rewiring diagonal*. Accordingly, we do NOT replace full CP-12 thermodynamics by that truncation. Instead we test the **precise independent statistical hypothesis** that positive local degree energies alone could make spatial geometry typical.

## 1. An exact finite-N injection theorem: local degree energy cannot enforce perfect regularity at fixed temperature

Declare a separate classical Gibbs measure on ALL labelled simple graphs G with N vertices and exactly Q=N edges:
\[
P_{N,a}(G)=Z^{-1}\exp[-a W(G)],\quad
a=\beta\lambda\in(0,\infty),\quad W=\sum_i d_i(d_i-1)/2.
\]
In the CP-12 fourth-order comparator, \(\lambda=2t_{\rm hop}^4/\Delta^3\). Here \(\lambda\) is held fixed and the declared Gibbs model is studied on its own terms.

For Q=N, convexity gives \(W(G)\ge N\), with equality if and only if **every vertex has degree 2**. Let \(\mathcal R_N\) be that 2-regular family.

**THEOREM:**
\[
\boxed{P_{N,a}(\mathcal R_N)
\le \bigl[1+\tfrac23 N(N-3)e^{-a}\bigr]^{-1}
\longrightarrow0\qquad (a\ \mathrm{fixed}),\quad N\ge4.}
\]

**Elementary proof:** For each G in \(\mathcal R_N\), choose an edge (a,b), designate one end b (2N oriented choices), and choose one vertex c different from a,b and from the other neighbor of a (N−3 choices). Replace (a,b) by (a,c), maintaining Q=N and simplicity. Exactly one vertex (b) now has degree1, one (c) has degree3, and the other N−2 vertices have degree2, increasing W from N to N+1. There are exactly \(2N(N-3)\) directed construction choices per G. Any output graph H has uniquely identifiable defective vertices b (degree1) and c (degree3). Its only possible 2-regular precursors remove one of the three edges (a,c) touching c and insert (a,b); hence **at most three preimages**. Thus \(\#\mathcal D_N\ge [2N(N-3)/3]\#\mathcal R_N\) for the class of these one-pair defects. Every G in \(\mathcal R_N\) has weight e^{-aN}; every H in \(\mathcal D_N\) weight e^{-a(N+1)}. Bounding the partition function below by these two disjoint classes proves the inequality. QED.

A necessary (not sufficient) size-dependent condition for \(P_{N,a}(\mathcal R_N)\ge p_*>0\) is
\[
\boxed{a\ge \log\!\left[\frac{2N(N-3)}{3}\frac{p_*}{1-p_*}\right]
=2\log N+O(1).}
\]
This follows directly by requiring our *upper bound* to be at least \(p_*\). The frozen CP-12 coefficient does not come with an \(N\)-dependent reinforcement of this kind.

This theorem does **not** mean a spatial phase can tolerate no defects. Real extended phases may be robust to a positive defect density. It establishes that *perfect* degree regularity cannot become typical at fixed positive temperature solely from any N-independent local energy coefficient of this form.

### Exhaustive controls

[check.py](check.py) enumerates ALL \(\binom{N(N-1)/2}{N}\) graphs through N=7, constructs all forward defect swaps, and verifies at most three reverse predecessors for each, then computes actual Gibbs probabilities and checks the bound.

| N | Number of 2-regular graphs | Distinct degree-1/3 defects reached | Largest reverse multiplicity | Exact P(2-regular), a=4 |
|---:|---:|---:|---:|---:|
| 4 | 3 | 12 | 2 | 0.931738 |
| 5 | 12 | 120 | 2 | 0.843978 |
| 6 | 70 | 1,080 | 3 | 0.775186 |
| 7 | 465 | 10,500 | 3 | 0.697556 |

The declining probability at substantial finite a is already visible. The analytic bound, not the small N fit, establishes the thermodynamic result.

A supplemental **standard sparse-graph entropy saddle** for this declared degree-weighted model gives
\[
p_d\propto\frac{z^d}{d!}e^{-ad(d-1)/2},\quad
\sum_d d\,p_d=2.
\]
For every finite a, \(p_1,p_3>0\). Illustrative asymptotic defect fractions \(1-p_2\) from this standard approximation are 0.525 (a=1), 0.181 (a=4), 0.0752 (a=6), 0.01088 (a=10). This saddle is **not** needed for the rigorous injection theorem and is not a separately proven CP-12 thermodynamic law. The small fraction at large a can still represent extensive O(N) defects in an N-vertex system.

## 2. Even PERFECT regularity leaves a separate fragmentation entropy

Now make the deliberately *more favorable* assumption that all degrees are exactly 2. Every simple graph is then a vertex-disjoint union of simple cycles of lengths \(\ell\ge3\), including the single spanning N-cycle as one option. Exact labelled counting uses the exponential generating function
\[
\frac{\#\mathcal R_N}{N!}=[z^N] A(z),\qquad
A(z)=\exp\left[\sum_{\ell\ge3}\frac{z^\ell}{2\ell}\right]
=(1-z)^{-1/2}\exp[-z/2-z^2/4].
\]
The number of labelled spanning N-cycles is (N−1)!/2. Consequently the **EXACT finite-N uniform probability**, conditional on perfect 2-regularity, is
\[
P_N(\mathrm{one\ ring}\mid\mathcal R_N)
=\frac{1}{2N[z^N]A(z)}.
\]
From the standard singularity asymptotic,
\[
\boxed{
P_N(\mathrm{one\ ring}\mid\mathcal R_N)
\sim \frac{\sqrt\pi e^{3/4}}{2\sqrt N}
\approx \frac{1.8761424159}{\sqrt N}\longrightarrow0.
}
\]

| N | Conditional probability of one ring |
|---:|---:|
| 6 | 0.857143 |
| 8 | 0.718563 |
| 12 | 0.571152 |
| 20 | 0.433017 |
| 40 | 0.301349 |
| 80 | 0.211410 |
| 160 | 0.148904 |
| 320 | 0.105085 |
| 1,000 | 0.059366 |

Thus **enforcing exact two-regularity is not enough**: the unrestricted cyclic orderings still generally consist of multiple components. A selected degree sequence and a selected macroscopic connected geometry are different results. This asymptotic is standard enumerative combinatorics and must not be claimed as GRUT novelty.

## 3. Even a fixed energy penalty per disconnected component is insufficient in this classical model

As a second hostile control, suppose hypothetically that **a positive extra cost Δ_c per cycle** has somehow been earned, on top of the exact degree-2 requirement. (CP-12 does not currently derive this term, so it is *not* part of its law.) At fixed inverse temperature β, define
\[
\theta=\tfrac12 e^{-\beta\Delta_c}>0,\qquad
A_\theta(z)=\exp\left[\theta\sum_{\ell\ge3}\frac{z^\ell}{\ell}\right]
=(1-z)^{-\theta}\exp[-\theta z-\theta z^2/2].
\]
The exact probability of one ring is \(\theta/[N[z^N]A_\theta(z)]\). Therefore
\[
\boxed{
P_{\beta,\Delta_c}(\mathrm{one\ ring})
\sim\Gamma(1+\theta)e^{3\theta/2}N^{-\theta}
\longrightarrow0
\qquad(\beta\Delta_c\ \mathrm{fixed}).
}
\]
There is also a simple **uniform finite-\(N\) bound** requiring no singularity analysis. Keeping just one-cycle and two-cycle terms in the positive partition sum, the ratio of *two-cycle weight* to *one-cycle weight* is exactly
\[
\frac{Z_{N,\mathrm{two\ cycles}}}{Z_{N,\mathrm{one\ cycle}}}
=\theta\sum_{\ell=3}^{N-3}\frac1\ell,\qquad N\ge6.
\]
Therefore
\[
\boxed{P_{\beta,\Delta_c}(\mathrm{one\ ring})
\le\left[1+\theta\sum_{\ell=3}^{N-3}\frac1\ell\right]^{-1}.}
\]
To keep one-ring probability at least \(p_*>0\), this *necessitates*
\(\theta\lesssim 1/\log N\), or \(\beta\Delta_c\gtrsim\log\log N+O(1)\). This is a rigorous necessary condition for the **hypothetically supplied** cycle-penalty model, not a viable CP-12 repair. The exact coefficient asymptotic above proves the stronger \(N^{-\theta}\) decline at any fixed \(\theta>0\).

Finite energy penalties per component only **slow** the algebraic decline; they do not make one connected ring thermodynamically typical at arbitrarily large N. The limit is subtle: with βΔ_c=4, one ring remains ~94.7% probable even at N=1000, but its asymptotic probability must eventually fall. This helps explain why strong finite-system loop preference can mislead.

This classical no-go can be evaded by a growing effective penalty, a fundamentally restricted configuration space, non-Gibbs states, or richer nonadditive many-body correlations. In particular it is **not** a no-go for full noncommuting CP-12.

### A third exact finite-graph reminder: Q=N−1

When Q=N−1, the minimal value of W is N−2, attained by degree multiset (1,1,2,...,2). A **spanning path** and a **disjoint triangle plus an independent path** (for N>=6) have precisely this multiset and identical W, although only one is connected. So the fourth-order degree term cannot distinguish them even at zero temperature in its declared perturbative comparator.

## 4. Relation to the owner 12th-order loop result

The owner's rigorous rational recursion establishes for bare fixed-graph energy at weak hopping:
\[
E(C_6)-2E(C_3)=
\frac{47080992982591078400}{310998849352074430281}t_{\rm hop}^{12}
+O(t_{\rm hop}^{14})
=0.1513863896303\,t_{\rm hop}^{12}+O(t_{\rm hop}^{14}).
\]
This **favors disconnected triangles at weak hopping**, despite identical regularity. The coefficient is source-derived from the pinned follow-up and has **not been independently recomputed by this present checker**. At frozen t_hop=1 the owner separately finds \(E_{\rm block}(C_6)=-0.635203537847\) and \(E_{\rm block}(C_3\sqcup C_3)=-0.637262028038\) (including the rewiring diagonal, regulator-controlled), again favoring disconnected geometry. These are not global ground-state comparisons.

Our additional result is independent of the detailed 12th-order coefficient. It shows that **(a) degree regularity, (b) global connectivity, and (c) extended spatial dimension** are three genuinely separate physical selection obligations. Two ports and degree two lead at most to cyclic one-dimensional neighborhoods. A 3D geometry, if it exists in CP-12, must be selected in a different coordination/interaction regime or through a richer relational notion of locality.

## 5. Concrete next discriminating physics experiment — no new law

Keep CP-12 frozen and test, within the *same state*, a **fragmentation-sensitive quantum observable**, not solely the degree variance or connected probability:

- At Q=N−1, prepare and distinguish a spanning path versus triangle+path, both with identical W; at Q=N, a spanning N-cycle versus a partition into smaller loops. Measure **graph class support** after full evolution, keeping matter/electric correlation and graph coherence.
- Test length-resolved connected-component counts \(C_\ell\), the largest-component size and heat-kernel/spectral dimension on actual observed components, not just one indicator.
- Compare conditioned relational source/detector responses for competing graph classes, retaining full rewiring diagonal and physically charged reference setup. Track whether local propagation survives for components with comparable degrees but different global topology.
- Extend exact linked-cycle matter-energy comparisons through l=8 or 12 using the original code, and report the first coefficient distinguishing cycle lengths; compare weak hopping and the frozen t_hop=1 *separately*.
- Preserve all state-preparation energy differences and electric cutoff convergence. A structural energy term may be dynamically earned, but cannot be inserted retrospectively because it selects the answer we want.

The crucial question is whether CP-12's **full** quantum geometry/matter state develops long-range correlations that override these local-energy and combinatorial tendencies. The present results supply exact adversarial baselines and early rejection criteria—not a completed theory of reality or a novel experimentally tested prediction.

## 6. Prior mathematical physics

- [Enumeration of labelled 2-regular graphs; Wilf's asymptotic, cited in this combinatorics discussion](https://math.stackexchange.com/questions/4118524/is-there-a-closed-form-formula-for-counting-2-regular-labelled-graphs).
- [Wormald, asymptotic connectivity of labelled regular graphs (1981)](https://doi.org/10.1016/S0095-8956(81)80021-4): notably fixed r>=3 differs sharply from r=2.
- [Bauer & Bernard, maximal-entropy random networks with given degree distribution (2002)](https://arxiv.org/abs/cond-mat/0206150).
- [Quantum Graphity, Konopka et al. (2008)](https://arxiv.org/abs/0801.0861).

All source methods above have significant prior art. Only the explicit application to CP-12's derived regularizer and the combined adversarial bounds belong to this project's model-specific work; novelty has not been externally established. Independent executable controls: [check.py](check.py); it writes results.json when run. Equivalent standalone controls were run locally before publication. The exact committed script is awaiting verbatim partner execution.
