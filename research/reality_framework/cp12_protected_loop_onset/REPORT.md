# CP-12 protected-loop onset: exact hopping-order bounds

**14 October 2026.** Additive analytical follow-up to 9d0f2ef. Frozen CP-12 Hamiltonian unchanged. This proves short-time order bounds, not spatial selection.

## Theorem A: protected core containing k edges

Let psi_0 be the frozen equal-amplitude neutral ray over all Q-edge graphs, with empty vertex fermions, filled active edge fermions and E=0. Let Pi_{\ge k} project onto configurations with at least k protected present edges (the electric-support protected core). Every hopping matrix element changes the electric support on **at most one** potential edge. Neutral rewiring commutes with every protected-edge projector, and all diagonal Hamiltonian terms do likewise. Consequently,
\[
\Pi_{\ge k}H^j\psi_0=0\quad\text{for every integer }0\le j<k.
\]
For finite N,Q and initial psi_0 in the domain of H^k, Taylor's integral remainder for the self-adjoint unbounded frozen Hamiltonian gives the **unregulated** bound
\[
\boxed{\Pr_T(|P|\ge k)=\|\Pi_{\ge k}e^{-iHT}\psi_0\|^2
\le \frac{|T|^{2k}}{(k!)^2}\|H^k\psi_0\|^2.}
\]
The bound holds for any real T but is informative only where its RHS is below 1. The same statement holds for any event requiring at least k protected edges, including a protected cycle of length at least k. In particular a protected triangle has probability O(T^6) as T tends to zero at fixed N,Q. This does not prove the leading coefficient is nonzero in every sector.

## Theorem B: fully occupied-edge flux cycle

Let Pi^{flux}_{\ge k} project onto configurations containing at least k present edges with occupied edge fermion (b=1) **and** nonzero electric field. For each such edge, its initial edge fermion b=1 must have undergone an even, positive number of hopping flips, hence **at least two hops on that particular edge**. Rewiring cannot change b or E, nor can it rewire an edge carrying flux. Thus
\[
\Pi^{flux}_{\ge k}H^j\psi_0=0\quad (0\le j<2k),
\]
and
\[
\boxed{\Pr_T(k\text{ occupied-edge flux edges})
\le \frac{|T|^{4k}}{((2k)!)^2}\|H^{2k}\psi_0\|^2.}
\]
In particular, an occupied-edge-flux triangle requires at least six hopping events, and its probability is O(T^12) at fixed N,Q. A flux-only *closed circulation* with n_v=0 is a subset of this event and obeys the same bound.

These are **short-time fixed-system** upper bounds; the norms \|H^k psi_0\| may grow rapidly with N. They are not uniform-in-N suppression results and cannot be used to infer a macroscopic onset time or thermodynamic spatial phase without further norm estimates. The first-order flux-protected fraction k/Q in 8ece3e2 need not obey T^12: it can arise from ordinary hole protection, whereas a fully occupied-edge flux triangle is a much more demanding event.

## Exact three-event reachability check

The companion checker enumerates all fermion and electric configurations reachable from the N=3, Q=3 neutral triangle by up to six gauge-compatible matter hops, with no electric cutoff and no rewiring (M=Q). It finds 44 distinct reachable states, the first fully protected triangle at hopping order 3 and the first flux-only triangle at hopping order 6. This verifies **reachability orders**, not quantum amplitudes or probabilities; amplitudes can in principle cancel. The general operator theorems above are independent of this finite check.

## Research consequence

The 9d0f2ef protected-core loop probabilities mix two physically different routes: (1) holes/vertex charges making edges nonremovable, and (2) genuine closed electric circulation on occupied edge modes. Their dynamical minimum orders differ by a factor of two in hopping count for a fixed loop. The next empirical gate should resolve protected-core cycles by **hole content versus fully occupied-edge flux content**, with energy-matched coherent controls, over several early times. A genuinely gauge-supported spatial loop phase should not be inferred from the easier hole-protection channel alone.

No dimension, causal cone, state-selection law or gravity is derived here.
