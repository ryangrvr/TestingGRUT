# CP12 operational scaling: a conditional growing-time propagation bound

2026-10-12. Additive analysis of frozen CP12, rooted at partner audit `c98a3fd1df5b15951c79eefab3a6bab2a835a321`. **No CP13, changed coupling, or modification to CP12's law.**

## Main result

The distance-shortcut audit is physically correct: neutral rewiring can lower the available bare electric-distance minimum without immediately redistributing the actual flux. Nevertheless, shortcuts do not invalidate every operational locality limit.

For a source B and detector A of norm at most one that preserve geometry, and any initial state on one definite graph, the full detector response obeys

\[
 |R_{\rm full}(T)-R_{\rm block}(T)|
 \le 4\Gamma^2T^2\frac{Q(M-Q)}{M^2},\qquad M=\binom N2.
\]

Here the block Hamiltonian is **P_G H P_G**, retaining the diagonal rewiring energy. This improves the general bounded-detector O(T/sqrt(N)) estimate to **O(T²/N)** for these physically relevant geometry-preserving operations at Q proportional to N. It controls virtual excursions and return paths, not only final graph-departure probability. The proof is in [THEOREM.md](THEOREM.md).

For prepared graphs of uniformly bounded degree, the compressed Hamiltonian has bounded nearest-neighbor hopping on its incidence graph. Its unbounded electric energy is on-site. The standard even-fermion Lieb–Robinson argument, with those on-site terms treated exactly, consequently gives

\[
 |R_{\rm full}(T)|\le C_0e^{-a(d_G-v_{\rm LR}|T|)}
 +4\Gamma^2T^2\frac{Q(M-Q)}{M^2}.
\]

This establishes a **conditional operational propagation regime of the unchanged law**: for T_N=o(sqrt(N)) and detector separation outside the bounded-speed front by a growing margin, the response tends to zero. The theorem assumes the prepared graph family; it does not derive its emergence, Lorentz invariance or gravity. It is a response upper bound, not proof of massless or ballistic waves.

The proof extends by convexity to arbitrary incoherent mixtures of graphs. For imperfect geometry selection and residual coherence, a useful sufficient bound is

\[
 |R|\le C_0e^{-a(d-v_{\rm LR}|T|)}+2p_{\rm bad}
 +4\Gamma^2T^2Q(M-Q)/M^2+2\|\rho-\mathcal D(\rho)\|_1.
\]

Here p_bad counts initial weight outside the specified bounded-degree/separated-detector family; D dephases graph sectors. This separates geometry quality, coherence and topology-change rate. Decoherence can affect one term without selecting the geometry or solving the other terms.

## Exact audit execution and two implementation repairs

The committed partner `check.py` was executed verbatim and failed at the small-N isolated-vertex calculation: it passed a negative factorial argument to `lgamma` when the binomial coefficient should be zero. A corrected copy then exposed a second edge case, division by zero when the expected number of isolated vertices is zero.

[audit_corrected.py](audit_corrected.py) makes only the following mathematical boundary repairs (plus a distinct output filename):

- Return negative infinity from logchoose(n,k) when k is outside 0..n.
- Return the vacuous connectedness bound 1 when E[I]=0, rather than dividing by zero.

The corrected exact script passes all shortcut, Gauss, energy, counting and asymptotic controls; results are in [audit_results.json](audit_results.json). The original published audit remains unchanged. These execution defects do not invalidate its mathematical arguments. No claim that the original file itself passed is made.

## What the numerical experiment actually measures

A prepared state on a ring or a hostile hub-rich graph is a coherent superposition of the neutral filled-edge state and one gauge-invariant particle-hole/electric excitation on edge (0,1). The source is B=n_0 and the detectors are bounded vertex occupations n_j. All states obey Gauss. Source and detector identities are supplied physical-reference labels; this is not a constructed autonomous relational apparatus.

We evolved the complete Q=N graph–matter–gauge Hamiltonian for N=4 and 5 at electric cutoffs 2,3,4, comparing to its exact initial-graph compression and, as an adversarial comparator, to the bare fixed-graph Hamiltonian that deletes the entire graph Laplacian. Both fixed microscopic times and T=sqrt(N),2sqrt(N) were tested. The largest dynamical run has 622,194 physical basis states. Norm, energy, graph leakage and both response bounds are checked during evolution.

For the detector j=floor(N/2), at T=1:

| Preparation | Full response | Difference from exact graph block | Difference between block and bare comparator |
|---|---:|---:|---:|
| N=4 ring | 0.000745431290 | -0.000001399204 | -0.000006654192 |
| N=5 ring | 0.000789742066 | -0.000000478135 | -0.000011430665 |
| N=4 hub-rich | -0.012346679240 | 0.000001725620 | 0.000130125592 |
| N=5 hub-rich | -0.010933232284 | 0.000001077914 | 0.000169381511 |

The last column matters: on a fixed graph the Laplacian contributes

\[
 \frac{\Gamma(M-Q)}M\sum_{e\in G}R_e,
\]

where R_e projects onto a removable neutral edge. Apart from a common constant, this is a negative local potential for nonremovable excitations. Its coefficient tends to Gamma. It cannot be discarded as a negligible finite-N correction. The complete partner report already retains it; the numerical control confirms why that distinction is necessary.

Cutoffs 3 and 4 agree on the T=1 responses at approximately 1e-15 in these runs. At 2sqrt(N), their largest detector difference is below 1.7e-6, and cutoff-4 boundary probability is below 4e-6. These are numerical convergence checks, not rigorous infinite-cutoff error certificates. Cutoff 2 alone has large late-time boundary populations and would be an inadequate control.

N=4 was additionally evolved to T=sqrt(N)/Gamma=20 at cutoffs 5 and 6. Ring graph-departure probability is 0.154854; the hub-rich preparation gives 0.147718. The ring's selected response differs from its block comparator by about -0.0004701; the hub-rich case by -0.0084607. The two highest cutoffs agree on the detector responses within 1.4e-12. This is real long-time topology change, not proof of geometrogenesis or equilibration. The analytic small-error bound is no longer strong at this scale.

Two sizes cannot establish a scaling exponent. Indeed b_N is not monotone from N=4 to 5 in this charge sector. The asymptotic claim comes from the theorem, not a fit to these points.

## Distance-resolved conditional response

The exact compressed ring Hamiltonian was also evolved at N=6,8 and three electric cutoffs. For N=8,T=1, the response at vertex distances 1,2,3,4 is approximately

\[
 (5.75974\times10^{-2},\;7.85240\times10^{-4},\;
 1.14480\times10^{-6},\;4.43268\times10^{-10}).
\]

This is an actual hierarchy of density responses from the same matter/gauge dynamics, not an unexplained small excitation gap. It belongs to the prepared-graph compression, not a full eight-vertex fluctuating-graph simulation. Finite responses are not exactly zero outside a sharp cone. The numbers do not determine a Lorentzian limiting velocity or dispersion relation.

## Independent controls and coherent counterexample

`impulse_control.py` checks the retarded response using independently propagated positive and negative source kicks; the discrepancy is below 6e-12. `operator_controls.py` tests the quadratic operator bound on 320 independently generated dense block-Hamiltonian/time combinations and verifies the CP12 diagonal identity directly.

For two neighboring graphs prepared with relative phase i, the graph probability changes at first order, with derivative magnitude Gamma/M=1/60 at N=4. An incoherent mixture of the same two graphs has only a quadratic change and obeys the stronger bound. [mixture_results.json](mixture_results.json) records both. This counterexample prevents extending the quadratic claim to arbitrary graph superpositions. Dephasing the uniform neutral graph state also raises the system energy by exactly Gamma Q(M-Q)/M, asymptotically Gamma nu N. The N=4,Q=4 value 2/15 was checked. Thus the favorable incoherent-preparation hypothesis is not a free replacement for a coherent low-energy state; its apparatus and energy source remain to be supplied.

## Actual joining experiment, separately from preservation

A further full-Hamiltonian test starts from a triangle plus an isolated vertex (N=4,Q=3), with the same local physical charge preparation. The detector is placed at the initially isolated vertex. Nothing postselects connected graphs.

| External time | Connected probability | Initially isolated detector response |
|---|---:|---:|
| 1 | 0.00141932 | 0.000009962 |
| 4 | 0.01416601 | 0.000299061 |
| 20 | 0.14037727 | 0.007315869 |
| 40 | 0.24251626 | 0.010050325 |

Electric cutoffs 4 and 5 agree on these connected probabilities within 9e-11, while the cutoff-5 boundary weight is below 8e-11. Total energy and norm are conserved. The exact fixed initial-graph detector response is zero because its vertex is isolated and Gauss forces its occupation to vanish. Joining and a detector response develop together under the unchanged law.

This demonstrates actual finite-system connectivity creation on longer times, not just a lower bound on a hypothetical formation time. It does not prove monotone convergence, a selected spatial geometry, a macroscopic assembly rate, or a stationary connected phase. Its N=4,Q=3 sector is an explicit finite control, not evidence for the Q/N=1 thermodynamic limit. See `assembly_control.py` and `assembly_results.json`.

## What the growing-time result does and does not buy

For a supplied D-dimensional bounded-degree family, linear size scales as N^(1/D). At times proportional to that diameter the topology correction scales as N^(2/D-1). It therefore vanishes for D>2, including D=3; D=2 is marginal for this bound, and a full-ring D=1 traversal is outside the guaranteed window. **This does not select three dimensions**: every D>2 shares the favorable exponent, and actual ballistic traversal has not been derived.

The result makes the proposed two-timescale possibility mathematically concrete on the preservation side. A prepared three-dimensional-type geometry can admit uniform detector control over increasing diameter-scale times even while the law permits microscopic rewiring. The large-system formation side remains open: substantial graph departure requires at least an order-sqrt(N)/Gamma time from a definite graph, with no guarantee of producing spatial geometry at that time. The thermal atypicality result also remains in force under its original assumptions.

CP12 now earns more than an instantaneous electric-distance identity: it has a conditional detector-propagation bound robust to its own virtual topology changes. It has not earned cosmological initial conditions, a selected dimension, a universal signal speed, Einstein dynamics or a new empirical prediction.

## Literature and reproduction

The fixed-geometry locality ingredient is established: [Nachtergaele–Sims, arXiv:1410.8174](https://arxiv.org/abs/1410.8174), particularly Theorem 3.1 and its interaction-picture treatment, and [Nachtergaele–Sims–Young, arXiv:1705.08553](https://arxiv.org/abs/1705.08553) for even fermionic observables. Their full papers were inspected. The CP12-specific block estimate and its combination with these bounds are author deductions pending independent review; no broad priority claim is made.

Run `assembly_control.py`, `operator_controls.py`, `mixture_control.py`, `impulse_control.py`, `audit_corrected.py`; `scaling.py --n 4 --L 4` and `--n 5 --L 4`; `scaling.py --n 4 --L 6 --long`; and `fixed_propagation.py`. NumPy/SciPy required except for the repaired audit. Set OPENBLAS_NUM_THREADS=1. Source imports the unchanged model from `research/constructive_physics/attempt12/model.py`. Outputs distinguish full, compressed, bare and corrected-audit calculations. The follow-up branch is additive to c98a3fd; CP12, all earlier candidates, canonical GRUT-RAI and Stage-3 Card 3 remain unchanged.
