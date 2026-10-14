# CP-12 flux-loop gate: exact orders, circulation and protection transfer

14 October 2026. Parent `fe90e7f88735aa9159f69c248c0a3b76a5680f2b`. Frozen law and neutral initial preparation unchanged. The earlier report and executed reachability checker were fetched by immutable Git blobs and rerun.

## 1. Three different loop properties

For a fixed simple cycle C of length l, define physical diagonal projectors:

* P_C: every edge of C is present and protected (nonzero electric support).
* F_C: every edge of C additionally has its edge fermion occupied. This is the occupied-edge-flux support in the onset theorem. It is an undirected support criterion.
* G_C: a directed electric circulation can be extracted around C on the full vertex/edge incidence graph.

The edge Gauss constraint is E_ie + E_je = 1 - n_be. Thus an occupied edge has opposite fields (f,-f). A hole edge can also have opposite nonzero fields, for example (2,-1). The latter carries one unit of through-flow plus the local edge charge. A hole does not by itself exclude circulation.

For G_C, both half-leg fields on every cycle edge must be nonzero and of opposite sign, with their directions consistent around the cycle. Their integer magnitudes need not be equal. Subtracting the smallest common positive flow around the directed incidence cycle leaves every Gauss divergence unchanged. This establishes the existence of an integer circulation component; it does not establish a stationary current or an isolated subsystem.

G_C implies P_C. An undirected F_C need not imply G_C in a larger charged network: its directions can fail to close. Conversely G_C can contain holes and need not imply F_C. We report the directed occupied-edge subset as well as general incidence circulation.

For final configurations, cycle counts decompose exactly into entirely occupied-edge cycles and cycles containing at least one hole. At the level of events, a configuration may contain both kinds; disjoint hole-only, flux-only and both events sum to the protected-loop event. The probability of at least one circulation loop must not be equated with the fraction of individual cycles that circulate.

All these diagonal predicates are gauge invariant and invariant under relabelling when summed over cycles. For same-cycle persistence we use one named triangle, with the correct physical stabilizer representation. This is a diagnostic choice, not spontaneous preference for those labels.

## 2. Exact formation amplitudes

The parent theorem's unregulated hopping-count bounds hold. Starting from the neutral ray, each protected edge requires at least one hop, and an occupied edge with flux requires at least two flips on that edge. The initial finite-support ray belongs to every finite power domain of the self-adjoint frozen H, so the Taylor integral remainder is available without an electric regulator.

Reachability does not alone prove a nonzero leading amplitude. `onset_exact.py` independently applies the signed full Hamiltonian using integer/rational arithmetic. For N=Q=3 it finds:

\[
\boxed{\Pr(P_C;T)=2T^6+O(T^8),\qquad
\Pr(F_C;T)=\frac{2}{225}T^{12}+O(T^{14}).}
\]

At order three, the two protected target configurations have all three cycle edges as holes. Their H^3 amplitudes have magnitude 6. At order six, the two occupied-flux targets are the two unit-circulation orientations. Each H^6 amplitude is -48. Thus 2(48/6!)^2 = 2/225. Both minimum orders are realized rather than cancelled. Odd probability terms vanish for this real initial ray, real H and diagonal real observables.

A hole edge on a through-flow cycle has E_ie + E_je = 1 and both fields nonzero. Their absolute values sum to at least three, whereas an occupied through-flow edge requires at least two. A circulation cycle with h holes therefore requires at least 2l+h hops. In the three-vertex example, the first hole-bearing circulation amplitudes occur at seven hops; the probability coefficient at T^14 is exactly 32/33075. This does not change the T^12 onset of circulation itself, which is supplied by the occupied channel.

For the uniform fixed-Q graph preparation let M=binom(N,2), and let

\[
\bar C_3=\binom N3\frac{\binom Q3}{\binom M3}
\]

be its mean triangle count (zero if Q<3). Minimal paths stay within that triangle and cannot include a rewiring or diagonal event. Hence

\[
\Pr(\text{some protected triangle})=2\bar C_3 T^6+O(T^8),
\]

\[
\Pr(\text{some occupied-flux triangle})
=\frac{2\bar C_3}{225}T^{12}+O(T^{14}).
\]

At this leading order there cannot be two different completed triangles. The leading coefficient for any electric-circulation loop is the same T^12 coefficient: longer cycles and hole-bearing circulation require more hops. Existing N4, N5 and N6 matrices reproduce the triangle coefficients; minimal flux paths have electric square at most six, below every retained budget used in this check.

The one-edge mean flux fraction is a separate observable. No occupied flux is possible after only two total hops from empty vertex modes: a return hop then cancels its own electric shift. Three hops can create a flux edge using an adjacent edge as a matter supplier. Each ordered adjacent-edge pair gives amplitude magnitude two and probability coefficient 1/9. Thus

\[
\boxed{\frac{\langle k\rangle}{Q}
=\frac{2(N-2)(Q-1)}{9(M-1)}T^6+O(T^8).}
\]

Total protection instead starts as (Q-<r>)/Q = 2T^2+O(T^4). Here k counts occupied-edge flux; n_v counts holes. Ordinary hole protection and the k/Q flux diagnostic are not interchangeable. All coefficients use the frozen hopping t=1. The statements are fixed-system onsets, not uniform-in-N long-time estimates.

## 3. A protected loop can outlive its occupied-flux classification

Every hopping event on an occupied edge creates a hole. By the edge Gauss law a hole remains protected. All other terms preserve electric support and edge occupancy. Therefore the full physical operator identity is

\[
\boxed{(1-P_C)H F_C=0.}
\]

For a normalized vector phi supported in F_C and in Dom(H^2),

\[
\|(1-P_C)e^{-iH\tau}\phi\|^2
\le\frac{|\tau|^4}{(2!)^2}\|H^2\phi\|^2.
\]

Define V_C = F_C Pi_(n_v=0). All vertex modes are then empty and all edge modes occupied. In two hopping events, either two holes are created or the first particle returns to its original edge and cancels that hop's field shift. Neither process erases an initially protected cycle edge. Rewiring and diagonal insertions cannot change this conclusion. Consequently

\[
\boxed{(1-P_C)H^j V_C=0\quad(j=0,1,2),}
\]

and geometric loss from this vacuum-flux subspace is O(tau^6). V_C alone need not orient a particular undirected cycle in an arbitrary larger network; the protection statement does not require that extra condition. In our N4,Q4 named-triangle cases it is a circulation because the fourth edge is a bridge and cannot carry vacuum flux.

The occupied-flux label can nevertheless change immediately. For a vacuum state, hops on different cycle half-legs have orthogonal fermionic targets. Its conditional F_C loss curvature is 2l at t=1. The stricter V_C event loses probability with curvature 2Q, since a hop anywhere creates a vertex particle.

For a supplied three-edge unit-circulation basis control, exact unregulated arithmetic gives

\[
\begin{aligned}
\Pr(\text{lose }F_C;\tau)&=6\tau^2+O(\tau^4),\\
\Pr(\text{lose }G_C;\tau)&=3\tau^2+O(\tau^4),\\
\Pr(\text{lose }P_C;\tau)&=\tfrac13\tau^6+O(\tau^8).
\end{aligned}
\]

A hop on a positive unit-flux half-leg leaves through-flow, e.g. (1,-1) -> (2,-1); a hop on the negative half-leg gives (1,0), which retains protection but breaks circulation on that edge. A leading all-hole triangle basis control can instead lose protected geometry with coefficient 3tau^2. These supplied basis controls demonstrate different local mechanisms, not energy-matched autonomous phases.

## 4. Persistence, quantum correlations and diagnostic backaction

For the unmeasured main state at waiting time t, write its purification as X (one column for a pure state). Let I denote an initial loop projector and J a later detector on the same named triangle.

The sequential diagnostic joint probability is

\[
j_{J,I}(\tau)=\|J U(\tau) I X\|_F^2,
\qquad q_{J|I}=j_{J,I}/\|IX\|_F^2.
\]

This includes an ideal projective filter at the waiting time. It is not an unmeasured classical history, free preparation, or autonomous outcome selection. Conversion F -> hole can preserve P and sometimes G, so the three later detectors are reported separately.

We also compute the unforced two-time correlator

\[
c_{J,I}(\tau)=\operatorname{Tr}[(U X)^\dagger J U I X]
\]

without interpreting it as a positive joint probability. For J=I, the noninterventional projector-change score is

\[
\|I U X-U I X\|_F^2=p_I(t+\tau)+p_I(t)-2\operatorname{Re}c_{I,I}.
\]

It lies in [0,1]. Normalized real autocorrelation uses Re(c)/sqrt(p_I(t)p_I(t+tau)); rarity alone cannot be mistaken for a stable record.

The computational subgroup fixing vertex 0 permits the triangle (1,2,3) to be resolved. It changes no Hamiltonian term or physical initial state. The exact inclusion J from the fully invariant basis satisfies J^dagger J=1 and H_stabilizer J=J H_invariant. Both are numerically checked. An independent 6,075-state raw physical matrix checks the same probabilities and correlators, the outgoing Gram matrix and purification compression.

For mixed references, an SVD of IX compresses purification columns, retaining every singular value. This is an exact right-isometry representation, not a physical measurement beyond the declared initial filter or a rank approximation.

The binary nonselective measurement work is

\[
W_I=-2\operatorname{Re}\operatorname{Tr}[(IX)^\dagger H(1-I)X].
\]

Only hopping contributes because I commutes with diagonal terms and rewiring. Its error is controlled by a bounded hopping operator. For a pure input the measured system's pinching entropy is h_2(p_I); this does not create entropy in a larger unitary system plus apparatus or derive an outcome.

## 5. Coherent energy references and conditioning

The two references are exactly those constructed at 9d0f2ef in the retained N4,Q4 invariant sector: spectral dephasing of the initial ray and the maximum-entropy state at mean energy 2. They remain configuration coherent. Their full-H mean energies equal their retained-H means because they have finite retained support. Their global entropies, variances and instantaneous matter/flux resources need not match the evolved pure state. They are not full-unbounded-H equilibrium states or theory-selected preparations.

Mean energy is matched BEFORE filtering. A selected flux branch can have a different energy, which is explicitly computed. For the main state, a norm-only error estimate for this unbounded conditional energy would be invalid. Instead use conservation of the initial second moment: ||H psi(t)||=sqrt(12). Since [H,I]=[Hhop,I] is bounded,

\[
\|I H I\psi\|\le\|H\psi\|+2\|H_{\rm hop}\|.
\]

The same inequality for the projected vector v gives

\[
|\langle IHI\rangle_\psi-\langle IHI\rangle_v|
\le e\,[\sqrt{12}+\|Hv\|+4\|H_{\rm hop}\|].
\]

Here ||Hv||^2 = ||H_P v||^2+||Bv||^2 and ||Hhop||<=sqrt(32). Divide its signed numerator interval by the positive event-probability interval. Thus the main conditional-energy intervals use an energy-moment estimate, not an unsupported bounded-observable claim.

## 6. Full-law errors

Use the earlier complete reachable electric-square projections and their full outgoing residual B. If the purification error is at most e in Frobenius norm, any positive contraction O satisfies

\[
\boxed{\langle O\rangle\in
[(\sqrt p-e)_+^2,\,\min(1,(\sqrt p+e)^2)]},
\qquad p=\langle O\rangle_P.
\]

This follows by applying the norm triangle inequality to sqrt(O)X. It is often sharper for rare events than 2e+e^2. If the projected event vanishes, the upper bound is e^2, not zero.

For a filtered path, e_y <= e_wait + e_aux. A correlator with projected path norms x,y has complex error at most e_x||y||+e_y||x||+e_x e_y. Sequential probabilities use the positive-contraction interval and conditional ratios use positive denominator intervals.

For the energy references, projected purification columns are eigenvectors. Their full-law unfiltered-path error obeys e_x(tau)<=tau||BX||_F, because the residual Frobenius norm is constant under their column phases. Filtered paths are bounded separately by the existing leakage estimator.

A detector-specific dual refinement also gives |<Ov,psi-v>|<=e d for a positive contraction commuting with the retained projector. Cauchy-Schwarz in the O form yields

\[
\langle O\rangle_{\rm exact}\ge\frac{(p-ed)_+^2}{p}\quad(p>0),
\]

while its upper bound is p+2ed+e^2. We intersect it with the basic amplitude interval. In the dense T1 circulation test even this refinement does not exclude zero; the numerical center is not promoted to a detection claim.

These are analytic full-law bounds evaluated in ordinary floating point. They are not interval-arithmetic roundoff certificates. Exact rational onset coefficients are distinguished from floating-point dynamical calculations. Displayed decimal intervals are rounded outward.

## 7. A density response in the same diagnostic state

Normalize the named F_C branch of the N4,Q4 state at wait 10. Apply the gauge-invariant source exp(-i eta n_0), evolve by the unchanged H for lag 0.5, and detect D=(n_1+n_2+n_3)/3. D is independent of graph distance and geometry bins.

Because n_0 is a projector, for every linear observable A its exact derivative at eta=0 is [f_A(pi/2)-f_A(-pi/2)]/2. This is not an approximate finite difference. A separate commutator-vector calculation checks its sign and normalization.

To control reclassification of the circulation detector, also set g(eta)=<G_C> and f(eta)=<G_C D>. The conditional density derivative inside the continuing-circulation readout is

\[
\boxed{\frac{d}{d\eta}\frac{f}{g}\bigg|_0
=\frac{f'(0)g(0)-f(0)g'(0)}{g(0)^2}.}
\]

We bound f(0),g(0),f'(0),g'(0) and propagate their intervals through this formula. This separates a within-circulation density response from merely changing that bin's probability. Final conditioning remains a postselection diagnostic; no local observer or free selection procedure is supplied.

A normalized main filtered pure state has input error <=2 e_wait/sqrt(p_projected). Source kicks are unitary; their work comes only from hopping and is explicitly reported. Positive-detector path intervals then bound both density and same-loop retention.

If hopping is removed after the same preparation, vertex density commutes with every remaining term and the susceptibility is exactly zero. In a vacuum-flux branch n_0 rho=rho n_0=0, so this density-phase source acts as the identity and the response is exactly zero even with full hopping later. These are mechanism controls, not changes to the main law. Since total vertex-fermion number is not conserved, a density response does not by itself distinguish migration of a preexisting particle from edge/vertex conversion, and does not establish a causal cone.

## 8. Network capacity and the sparse-sector limitation

A connected spanning union of directed occupied-flux cycles has at least N occupied flux edges. Hence k>=N and n_v+r<=Q-N. At Q=N it must be a vacuum, zero-removability simple ring. Its density-phase source is inert for the reason above.

A spanning incidence-circulation network may include holes. It still needs at least N through-flow graph edges, but it can coexist with vertex matter. At Q=N it is again a simple ring and r=0, while n_v may be positive. Occupied through-flow edges cost E^2>=2; hole through-flow edges cost E^2>=5. Therefore, in this sector,

\[
\boxed{E^2\ge2N+3n_v\quad\text{on a spanning circulation ring}.}
\]

These two network criteria are different. Charged circulating rings are not ruled out by the occupied-edge vacuum observation. Both ring families have one-dimensional large-radius growth if made large. This Q=N sequence of spanning circulation networks therefore cannot establish three-dimensional extended geometry; the statement is about these circulating networks, not every possible sparse graph metric.

In N6,Q12 the existing budget10 projection excludes every spanning circulation network, which requires E^2>=12. Its projected zero is therefore not evidence of a full-law zero. The rigorous analytic upper bound from the inherited norm error is e(T)^2, approximately 0.01052034 at T1.
