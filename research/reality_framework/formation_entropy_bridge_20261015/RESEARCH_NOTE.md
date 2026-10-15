# GRUT constructive continuation: graph entropy versus induced regularization

**Status:** new analytic benchmark and research design, not a claim of geometric selection.  
**Parent:** constructive sprint `c2120bc719ba81f35d01e24b4dc4bbd13bddbd2a`. Frozen CP-12 Hamiltonian and original neutral ray are unchanged.

## Question

The constructive sprint derived a positive fourth-order coordination term in the *controlled static heavy-matter limit*, proportional to `sum_i d_i^2` at fixed graph-edge count. Does this energetic preference overcome the enormous number of irregular graphs, and can it explain dynamical geometry formation?

The first half can be quantified as a **combinatorial benchmark**. The second half is not implied by an energy calculation: the original CP-12 state evolves unitarily, not by an assumed Gibbs distribution.

## Exact graph-combinatorial facts

Fix `Q=2N` undirected edges among `M=N(N-1)/2` possible edges, with `N>=5`. The degree sum is `sum_i d_i=4N`. Hence

```
sum_i d_i^2 = 16N + sum_i (d_i-4)^2.
```

The equality follows by expanding the square and using the fixed degree sum. Thus the static fourth-order term has a unique *degree-sequence* minimum at `d_i=4` for all i, whenever a simple 4-regular graph exists. This does **not** distinguish the different 4-regular topologies or select intrinsic dimension.

For a uniformly chosen labelled Q-edge graph, the degree of any fixed vertex is hypergeometric with population `M`, marked edges `N-1`, and draws `Q`. Consequently

```
E[d_i] = 4,
Var(d_i) = Q (2/N)(1-2/N) (M-Q)/(M-1),
E[sum_i(d_i-4)^2] = N Var(d_i).
```

At `N=6`: `M=15,Q=12`, `Var(d_i)=4/7`, and exactly 15 of 455 labelled graphs are 4-regular (the octahedral class). Their fraction is `15/455 = 0.03296703297...`.

At `N=7`: `M=21,Q=14`, `Var(d_i)=1`, and 465 of 116280 labelled graphs are 4-regular (360 complements of C7 plus 105 complements of C3 union C4). Their fraction is `465/116280 = 0.003998968...`.

As `N -> infinity`, `Var(d_i) -> 4`. Thus a typical uniform graph has a degree-irregularity cost of order `4N`, even though all graphs share the same mean degree.

For fixed degree 4, the standard configuration-model asymptotic count of simple labelled regular graphs yields

```
log [ # labelled simple 4-regular graphs / binom(M,2N) ]
  = N (4 log 4 - log(4!) - 4) + o(N)
  = -1.632876... N + o(N).
```

This asymptotic invokes the standard fixed-degree regular-graph enumeration result and is not a new GRUT theorem. The finite N=6,7 counts above are exact.

## Induced energy scale and the central missing step

In the *controlled* heavy-matter static-band expansion of the sprint, the fourth-order graph-dependent term is

```
Delta E_4(g) = alpha sum_i (d_i-4)^2,
alpha = t^4/(8m^3) > 0,
```

after removing graph-independent constants at fixed Q. The uniform-graph mean penalty is therefore `alpha N Var(d_i)`, asymptotically `4 alpha N`.

Both the degree-selection energy difference and the log-number disadvantage of regular graphs are extensive in N. Neither automatically dominates the other. More importantly, an equilibrium free-energy argument `E - T S` **cannot** be used to claim CP-12 formation: no thermalizing reduced graph dynamics, bath, effective temperature, or detailed balance has been derived from the original coherent neutral ray. At frozen m=t=1, the fourth-order static-band expansion is not controlled; the exact Feshbach calculation establishes an interaction but does not establish the large-N coordination preference.

**Conclusion:** the static mechanism can favor degree regularity in an auxiliary limit, but the number of irregular alternatives is exponentially larger. A dynamically justified graph-selection mechanism must be computed from the actual initial state and frozen Hamiltonian.

## Constructive next calculation: a charged-sector kinetic theory with a controlled error

1. Define `P` as a physically gauge-invariant charged-sector projection using circulation, matter density and graph structure. Do not select a target dimension in P.
2. Derive an exact finite-time projected evolution identity with memory kernel, e.g. `P exp(-iHt) P` or Nakajima-Zwanzig / Feshbach, retaining `PH(1-P)` and initial `(1-P)psi_0`. State clearly that a static Feshbach resolvent alone is not a Markovian dynamical law.
3. Use the exact charged formation coefficients from N5–N7 as short-time boundary data; compare the reduced kernel against full unfiltered finite evolution where feasible. Test for finite-time occupation and return, not just nonzero Taylor coefficients.
4. Measure in the **same unfiltered state**: charged-network weight, graph persistence, degree variance, intrinsic ball growth, and source-to-detector response, with total event probability. Include controls matched in Q, matter density, degree sequence and spectral resources.
5. If a thermal or stochastic approximation is proposed, derive its regime, memory-time and error from the unitary model, and compare the resulting graph weights against the exact uniform-graph combinatorial baseline.
6. Demand an operational geometry discriminator that rejects 1D, 2D, higher-dimensional and random 4-regular alternatives. Regular degree is not three-dimensionality.

## Verdict

**New analytic benchmark:** degree regularity is an extensive energy-vs-multiplicity competition; uniform 4-regular graphs are exponentially rare among Q=2N labelled graphs.  
**Not established:** equilibration, effective temperature, 3D selection, dominant charged phase, continuum gravity, or a physical prediction beyond the input quantum Hamiltonian.  
**Decisive missing bridge:** a controlled, state-dependent, finite-time charged-network kinetic mechanism derived from the original neutral ray.
