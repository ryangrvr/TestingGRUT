# GRUT construction: dimension-neutral circulation density and its exact memory balance

**Status:** analytical construction and executable next-test specification. No new simulation claimed.  
**Parent:** published charged-memory bridge `5d634060aa5c0d25740f1c7db6e9f1f1d41b9c1a`. CP-12 Hamiltonian and neutral ray unchanged.

## 1. A physically meaningful target without a preferred lattice

At fixed N,Q=2N, for each physical basis configuration define its active graph g, its integer half-leg fields E, and vertex occupations. Define the directed **through-flow subgraph** F(g,E) to contain an active edge e=(i,j) if E_ie=-E_je is nonzero, oriented according to its sign. This is a graph-and-field diagnostic; no Euclidean embedding or chosen dimension enters.

Define C(g,E) as the subset of through-flow edges belonging to at least one directed cycle. A concrete finite algorithm computes it by strongly connected components: retain an edge (u,v) if u and v lie in the same nontrivial strongly connected component of F. Let `c=|C|`, `rho_c=c/Q`, `rho_v=n_v/N`. These are diagonal, gauge-invariant physical observables; the cycle test depends on electric flux, not on a choice of gauge potential.

**Caution:** with nonzero vertex charges, directed paths may terminate at matter, so through-flow and circulation must be distinguished. A charged loop event need not require all Q edges to circulate. Use thresholds `P_(a,b)=1{rho_c>=a and rho_v>=b}` with declared a,b in (0,1], and report sensitivity to thresholds. A spanning-cycle or connectedness requirement is an additional diagnostic, not part of the primary loop-density definition.

This choice is *not* the previous full-circulation event and its short-time coefficients cannot be silently reused. The N5-N7 full-circulation coefficients are exact boundary checks for the corresponding stricter projector only.

## 2. Exact unitary flux into an event

For any orthogonal event projector P, let Z=I-P, A=PHP, C=PHZ, D=ZHZ, and x=P psi, z=Z psi. Then

```
d/dt ||x||^2 = 2 Im <x,C z>
```

with the convention inner product conjugate-linear in the first argument and Schrödinger equation i dot psi=H psi. Check: dot x=-iAx-iCz; 2 Re<x,-iCz>=2 Im<x,Cz>.

For initial neutral state x(0)=0, z(0)=psi_0, exact elimination gives

```
dot x = -i A x -i C exp(-iDt) psi_0
        - integral_0^t C exp(-iD(t-s)) C† x(s) ds.
```

The memory kernel K(t)=C exp(-iDt) C† is a matrix-valued positive-type correlation function, but its entries and instantaneous contribution to d||x||²/dt need not be nonnegative. In particular, **no autonomous positive formation rate** follows from the nonzero early-time amplitude.

By substituting the exact complement solution, the probability current splits into an initial-source contribution and a coherent return contribution:

```
dot p(t) = 2 Im <x(t), C exp(-iDt) psi_0>
           - 2 Re <x(t), integral_0^t K(t-s)x(s) ds>.
```

The second sign follows from z=e^{-iDt}psi_0 -i integral e^{-iD(t-s)}C†x(s)ds, and 2 Im[-i q]=-2 Re[q]. These terms can each change sign at later times; the decomposition is exact and does not assume Markovianity.

## 3. A useful exact conservation / current bound

Since C is bounded for fixed N,Q, Cauchy-Schwarz gives

```
|dot p(t)| <= 2 ||C|| sqrt(p(t)(1-p(t))).
```

Set theta(t)=arcsin sqrt(p(t)) (interpreted by continuity at p=0,1). Then `|dot theta| <= ||C||` almost everywhere, yielding

```
p(t) <= sin²(min(||C|| t, pi/2))
```

for initial p(0)=0. This is a generic event-coupling speed limit; it is not an improved CP-12 large-N bound unless ||C|| can be bounded sharply for the chosen circulation-density event.

For a time-independent H and any initial state, the spectral time average of P is governed by energy-block coherences:

```
lim_(T->infinity) (1/T) integral_0^T p(t) dt
= sum_lambda <psi_0, Pi_lambda P Pi_lambda psi_0>,
```

where Pi_lambda projects onto the exact energy eigenspace. In the unbounded-rotor finite-N model this identity needs the appropriate pure-point spectral condition; do not claim it unconditionally in a possible continuous-spectrum setting. For the finite electric-cutoff models it is exact. This is a practical way to distinguish transient formation from persistent spectral support without fitting a dissipative law.

## 4. Geometric observables that cannot be replaced by degree regularity

On every basis graph g, measure:
- connected component fraction and cycle rank `Q-N+components(g)`;
- degree variance `N^-1 sum_i (d_i-4)^2` as a diagnostic, never a filter;
- intrinsic ball count `V_g(v,r)=|{w:d_g(v,w)<=r}|`, averaged over vertices and (where physically justified) conditioned on an actual record;
- heat-kernel return `h_g(s)=N^-1 Tr exp(-s L_g)`, `L_g` the combinatorial graph Laplacian; where a scaling window exists define running spectral dimension `d_s(s)=-2 d log h_g/d log s`;
- gauge-invariant source-to-detector density response with receiver selected independently of the final graph.

Do not interpret a diameter-two N4 graph as evidence for a spatial dimension. A finite system with no scaling window has no established continuum dimension. For coherent graph superpositions, these diagonal geometry observables have well-defined expectation values, but a claim of *one classical persistent geometry* additionally requires an operational record or a decoherence mechanism, with its cost included.

## 5. The next executable benchmark

1. First implement the circulation-density projector on existing N4,Q4 full-basis data and verify that threshold a=1 with b>0 recovers the published charged-circulation event. This is an equivalence check, not new large-N physics.
2. On N5,Q10, compare full-circulation and intermediate loop-density thresholds. Record p(t), source current, return current and probability conservation over t=0..10; verify current by independent finite difference.
3. For N5,Q10 and N6,Q12, use orbit symmetry and charged-port Krylov spaces seeded by both psi_0 and C†; check all outgoing electric-budget residuals. Do not infer scaling speedup from seven ports at N4.
4. Use early-time exact coefficients as checks only for matching target definitions. Benchmark a=1 and b>=1/N against prior N5/N6 results.
5. Define a persistent-event correlator `S_P(t,tau)=||P exp(-iH tau) P exp(-iHt) psi_0||²` and the conditional ratio S_P/p(t), explicitly noting this is a two-time sequential-projection diagnostic and changes the unmeasured history. Separately measure the unmeasured one-time probability p(t+tau). Never conflate the two.
6. Compare dimensionally different graph families and random degree-matched controls, and report formation probabilities rather than only conditioned responses.
7. Acceptance: publish a reproducible code path, finite-size and electric-cutoff convergence, conservative full-law residuals, and at least one nontrivial finite-time lower bound for an extensive loop-density event. If not achieved, publish an explicit numerical obstruction, not a selected spatial phase.

## 6. Scientific decision

The exact sourced memory equation already proves why irreversible formation-rate extrapolations fail in the N4 control. The new contribution here is a dimension-neutral charged-loop observable, an exact source-versus-return probability-current identity, and an operational specification that keeps probability, geometry and propagation attached to the **same evolving state**. No claim of macroscopic spatial emergence, three dimensions, Lorentz symmetry or Einstein gravity is made.
