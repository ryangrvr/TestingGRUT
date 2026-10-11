# CP11 scientific checkpoint: binding generates connectivity, but fails spatial selection

## Outcome

One new candidate was frozen before execution at local commit ee63567. It supplies a graph–fermion Hamiltonian in which the same hopping interaction produces propagation and an energy preference among graphs. It improves on CP10 in two ways: connections can lower matter energy, and conditional propagation resolves graph distance. Exhaustive conditional enumeration and full finite quantum calculations favor connected small systems without imposing connectivity.

The completed adversarial analysis finds two failures of the intended physical interpretation. The derived free-fermion energy rewards locally tree-like competitors over the regular spatial lattices tested. A full-Hamiltonian bound shows quantum connection mixing cannot reverse a fixed extensive energy disadvantage asymptotically at the frozen normalization. Separately, every fixed positive temperature gives asymptotically dense mean connectivity. No physical spacetime or gravitational theory has been established.

Verdict: **CONDITIONAL FINITE MECHANISM; SPATIAL-SELECTION TARGET FAILED for the tested lattice families and fixed-positive-temperature sparse phase.** This is not a universal no-go against interacting emergent geometry or an exhaustive classification of all CP11 ground graphs.

## Supplied physics

See LAW.md for the exact unmodified Hamiltonian. Unlike CP09/10, these are undirected connection graphs, not posets. Canonical fermions live on vertices AND every possible edge; there are O(N²) supplied edge modes and exactly M particles. The mass, hopping, edge cost, particle inventory, quantum formalism, external time, and Γ/N graph-flip scaling are inputs. No dimension or preferred graph is supplied. Counting labelled graphs and summing the canonical fermion sector define the finite ensemble.

The common mM subtraction makes the empty filled-edge state have zero matter energy. It cannot be credited as a solution to gravitational vacuum energy. The graph Dirac dispersion is established conventional structure, not evidence of physical Standard Model spinors. Previous relational clocks, causal-order interpretation and transport-holonomy results are not automatically inherited by this different action.

## Derived results

For fixed graph G:

 U(G)=κE(G)-Tr[sqrt(m²I+t²L_G)-mI].

This is obtained by filling the negative one-particle spectrum, independently checked against the complete incidence matrix spectrum. The full theory conserves fermion number and is covariant under simultaneous vertex permutations including edge-orientation and fermionic wedge signs.

At m=t=1, κ=.5, Γ=.1:

 -N/4 ≤ E0(N) ≤ .05(N-1).

The ground mean edge count is O(N), with an explicit loose bound in DERIVATION.md. These facts earn an extensive lower bound and sparse mean connectivity at zero temperature; they do not earn a physical spatial continuum.

## Finite calculations

Every labelled graph was included through N=6 in conditional enumeration: 33,866 graphs in total. At N=6 the ten conditional minimizers are connected degree-three graphs. N=5 minimizers have degrees 2,3,3,3,3. These are small-system facts, not continuum extrapolations.

Full graph–fermion ground calculations:

| N | Hilbert dimension | E0 | Connected probability | Mean edges |
|---|---:|---:|---:|---:|
| 2 | 6 | -.190489394 | .966898402 | .966898402 |
| 3 | 160 | -.427613105 | .990966602 | 2.770290245 |
| 4 | 13,440 | -.644784918 | .994697315 | 4.850662019 |

The full quantum ground differs from the static minimizer: at N=4 the complete graph has probability only .2073 although it is the conditional energy minimizer. This is why the full-Hamiltonian analysis was necessary.

## Larger competitors and the quantum loophole

| Conditional family | Energy per vertex in its large-size limit |
|---|---:|
| 1D ring | -.177609971862 |
| 2D square lattice | -.184539293475 |
| 3D cubic lattice | -.100067375815 |
| 2D honeycomb lattice | -.197229928022 |
| Locally tree-like degree-three graphs | -.197366056730 |

The closed-walk expansion proves a fixed-degree lower bound attained in the locally tree-like limit. This is stronger than selecting a favorable numerical competitor. It does not prove every minimizer is an expander or exclude every conceivable manifold discretization.

For the full quantum problem, letting U_min be the conditional minimum and a_N=Γ(N-1)/2:

 U_min+a_N-O(sqrt N) ≤ E0 ≤ U_min+a_N.

Ground weight in any graph class lying at least δN above U_min is O(1/sqrt N) for fixed δ>0. Thus the positive energy-density deficits of the listed lattice sequences cannot be erased by the frozen quantum mixing term. The conservative bound has no useful small-N crossover prediction.

## Thermal check

For the full canonical quantum ensemble at any fixed β>0:

 lim_{N->infinity} <E(G)>_β / binomial(N,2) = 1/(1+exp(βκ)).

This is a positive constant. Graph entropy therefore yields dense mean connectivity at every fixed nonzero temperature, despite sparse ground-state bounds. The proof controls the full operator norm of the matter and graph-mixing remainder relative to the O(N²) edge-count term; it is not a classical approximation. Exact small-system thermal traces verify the noncommuting thermodynamic derivative used in the argument.

## Propagation and physical limits

On a fixed path, the vertex-to-vertex amplitude at graph distance r first appears at order time^(2r), because each vertex edge takes two incidence hops. Direct matrix exponentiation agrees with that expansion. On a supplied spatial lattice the dispersion is ω²=m²+4t²Σ sin²(p_a a/2), with the familiar quadratic continuum limit. This is conditional graph propagation, not a derived Lorentzian light cone. The full geometry still has allpairs connection-changing interactions and external time.

No graviton, gravitational universal coupling, observed dimensionality, quantum outcome selection, or empirical prediction beyond conventional graph–fermion theory was derived. No experimental comparison can be claimed from the finite tests. These are explicit missing sectors, rather than additional unexplored numerical tolerances.

## Executed controls and reproducibility

Run with Python, NumPy, SciPy; OPENBLAS_NUM_THREADS=1:

1. `python research/constructive_physics/attempt11/check.py`
2. `python research/constructive_physics/attempt11/large_graph_controls.py`
3. `python research/constructive_physics/attempt11/adversarial_checks.py`
4. `python research/constructive_physics/attempt11/thermal_control.py`

Checks completed: exact conditional enumeration; filled one-particle versus spectral energy agreement; full canonical joint spectra and eigenvector residuals; independent dense Jordan–Wigner Hamiltonian equality; signed permutation covariance at N=3,4; path propagation order; extensive and mixing bounds; lattice integrals at three resolutions; Kesten–McKay normalized spectral integration; declared-seed regular graph controls; exact noncommuting finite-temperature derivative and partition bound.

The random regular graphs are samples, with seeds and rejection counts recorded; the theorem does not rely on a lucky sample. The full N=5 quantum system has 3,075,072 basis states and was not diagonalized. The analytical large-N statements do not require that calculation. No stability under all possible added EFT operators, interacting renormalization flow or continuum completion has been certified. All derivations and scripts remain author-level pending independent review.

## Research consequence

The construction has earned matter-induced connectivity and exposed why that achievement falls short. Free fermionic binding penalizes a density of short loops; familiar spatial lattices contain them. A successor would need an independently motivated interacting mechanism that changes this competition and controls finite-temperature graph entropy while also producing usable physical propagation. Arbitrarily adding a favored-loop reward would merely supply the requested geometry. CP11's parameters and particle sector will not be changed to hide these failures.

Primary literature and proofs are in DERIVATION.md. The GRUT research mission remains active; this candidate is not a theory of reality.
