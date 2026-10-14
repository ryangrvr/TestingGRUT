# CP-12 unmasked formation: integrated benchmark and precision gate

Status: **FORMATION NOT ESTABLISHED — SIZE AND ELECTRIC-REGULATOR GATES UNRESOLVED.**

Parent: `22ebf4f6a7202b107fa04e4a90dd6e33e9c96dc5`. Frozen Hamiltonian and canonical GRUT-RAI are unchanged. This is an additive computational research packet, not an independently audited theorem of emergent space.

## What changed

The test now starts from the equal-amplitude neutral superposition over **all** Q-edge graphs, with Q=2N. It contains no spatial graph, topology mask, graph-dependent phases, or geometry classifier. The ray still supplies a special pure boundary condition. Its energy is Q/2, its energy variance is 2Q, its min-entropy is zero, and its neutral-reference information deficit is log₂ binom(M,Q). A preparation resource construction and its missing physical controller are explicit in DERIVATION.md. This packet does not derive that boundary state from frozen CP-12.

This neutral state is stationary under neutral rewiring alone. Unlike the preceding phase-mask experiment, any change in its graph distribution under the full law requires coupled matter/electric dynamics. That observation identifies the relevant mechanism but does not prove that the mechanism forms space.

The numerical solver retains a finite total electric-square budget, keeps every frozen-law diagonal, and records all outgoing full-law transitions in a separate matrix. An analytic Duhamel bound, evaluated numerically, measures omitted-state error. It does not infer accuracy merely from conserved energy or a stable norm. The signed quotient fixes a source vertex for response measurements, allowing representations excluded from the fully invariant formation sector.

The diagnostics are graph-by-graph intrinsic volume growth, exact packing of vertex-disjoint triangles/squares, sequential loop-class survival, and a gauge-invariant linear density response. The fully integrated five-event control evaluates them together. At six events, source and formation calculations use different exact symmetry compressions of the same preparation and electric-budget rule. Source response is measured at wait=lag=0.1; formation diagnostics include times 0.1 and 0.2. The separate six-event persistence control uses wait=lag=0.05. These times are not conflated as one measurement record.

## The finite-size trap is exact

For every graph with M−Q<N−1, connectedness and diameter at most two are automatic. At Q=2N:

| N | Q | Geometry available in the entire sector |
|---:|---:|---|
| 5 | 10 | Only the complete graph K₅ |
| 6 | 12 | Every graph connected, diameter two |
| 7 | 14 | Nontrivial distance structure becomes possible, but remains very short |

Root-averaged B(1)=1+2Q/N=5 is fixed in all these sectors. At six events B(r)=6 for every r≥2, regardless of Hamiltonian or preparation. Neither this ball curve nor high connectedness can evidence emergent three-dimensional space.

The seven-event initial state already has B(2)=6.9793602 and B(3)=6.9984520. Nearly the entire graph is inside two steps before evolution. No dimension is fitted to these saturated curves.

## Computed results and their error controls

**Five-event response control.** With electric budget 10, a local n₀ kick at waiting time 0.1 gives, after lag 0.1,

    χ₁ = 6.5829686884×10⁻⁶,
    conservative absolute error upper = 1.1100494431×10⁻⁶.

The analytic bound is evaluated with 1,024 midpoint intervals and ordinary floating point. The sign is resolved under that bound; interval arithmetic is not claimed. The budget-eight pilot gives 6.5829694898×10⁻⁶, consistent with this value, but its looser error estimate did not resolve the sign. The complete graph's loop event is certain and persistent before the experiment. This is a gauge-invariant RESPONSE CONTROL, not propagation over an extended geometry or formation of loops.

**Six-event candidate response fails the error gate.** At the same waiting time and lag, electric budget 6 gives

| Final-distance bin | Projected susceptibility | Conservative absolute error upper |
|---:|---:|---:|
| 1 | +6.4336206014×10⁻⁶ | 1.1923224960×10⁻³ |
| 2 | −1.2447905321×10⁻⁶ | 1.1923224960×10⁻³ |

Both uncertainty intervals contain zero. We cannot credit the distance-two number as an established full-law signal. The source calculation has 147,153 retained signed-quotient basis vectors and 180,850 outgoing quotient orbits. A large matrix and accurate energy conservation do not repair the error gate.

The fully invariant six-event calculation has 24,830 retained basis vectors and 30,201 outgoing orbits. At time 0.2 its total graph-measurement variation error is bounded above numerically by approximately 3.8440×10⁻⁴. Its very small triangle and square changes lie below conservative observable error bars; they are not certified formation effects. Intrinsic volume saturation is exact independently of numerical errors.

For the separate six-event loop-class control at wait=lag=0.05, the initial event probability is approximately 0.9560439560, raw conditional survival 0.9999935250, and joint-probability absolute error bound 1.5890×10⁻⁵. That class was already highly probable in the initial graph distribution; this does not demonstrate assembly or persistence of an individually identified loop. The apparent departure of the survival ratio from one is not resolved by the joint-probability bound.

**Seven-event long-time pilot is inconclusive.** At budget 3, the invariant closure has 59,041 basis vectors representing 65 graph isomorphism classes. At time 2, projected norm error is about 9.1×10⁻¹⁵ and energy error about 6.6×10⁻¹⁴. Yet the omitted-state norm bound is approximately 6.41, giving only the trivial graph total-variation bound one. This calculation cannot decide finite-time formation, even though its internal conservation checks look excellent. The raw record is preserved to expose this failure, not to support a physical conclusion.

The five-, six-, and seven-event runs are a declared COMPUTATIONAL size series with uneven regulator adequacy. They are not a controlled thermodynamic extrapolation. In particular, the failed seven-event bound cannot be repaired by combining it with the five-event control.

## Source-detector interpretation

The source and detector commute with Gauss constraints. The response is a derivative with respect to a physical local density kick, not an initial correlation mislabeled as propagation. However, detector shells use the final graph distance and can respond through geometry redistribution as well as matter changes. A resolved susceptibility alone would not establish a propagating photon or a light cone. The phase kick is an external diagnostic intervention, not a controller that prepares the formation state.

## Verification

Fresh execution of the parent verifier passed its parent checker and four cubic controls. The new independent controls compare the new quotient with complete regulated physical matrices at N=3 and N=4, both with and without a marked source. They verify:

- the projected Hamiltonian intertwiner (errors below 2.2×10⁻¹⁴);
- the outgoing-leakage Gram matrix (errors below 1.1×10⁻¹⁴);
- actual direct full-matrix evolution errors below the analytic error upper sums;
- the source susceptibility against an independent central finite pulse (difference about 1.1×10⁻¹³);
- the exact fixed-Q one-step ball identity;
- the user's rational coefficient reconciliation, 3(181/63,000,000)=2353/273,000,000.

The comparison full matrices have a half-leg cutoff two. For these small checks that cutoff includes every outgoing transition from the retained budget-three space. Direct-evolution comparisons test the bound on that full regulated comparison matrix; the bound's unbounded-law applicability follows separately from its Duhamel proof, not from declaring the comparison matrix unregulated.

## Gate ruling and next work

This packet **does not pass** the proposed integrated spacetime-formation gate. It establishes a usable masked-free boundary-state benchmark, correct joint diagnostics, and an explicit way to reject numerically misleading results. It does not show a physical failure of CP-12, autonomous initial-state selection, a persistent spatial phase, dimensional scaling, or gravitational propagation.

The immediate numerical target is now a quantitative one: reduce the six-event source-response absolute uncertainty below roughly 10⁻⁶, and obtain meaningful electric-field error control at sizes whose intrinsic balls span a growing radial range. Raising the budget or developing an observable-specific bound may help the first target; merely adding time steps does not remove actual omitted-state leakage. A matrix-free or more efficient symmetry/occupation construction is needed before a controlled large-size formation claim is plausible. A successful abstract quantum preparation circuit still would not establish that frozen CP-12 autonomously selects its input state.

The scientific question remains open. The present obstacle is demonstrably inadequate size and numerical control, so this calculation cannot yet choose among all-to-all rewiring, the state family, or missing state selection as the physical cause of a formation failure.
