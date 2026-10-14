# CP-12 distance-history gate: backaction, shortcuts, and a suppression condition

**Ruling: a resolved rewiring-mediated response survives between initially disconnected components in the tested preparation. A growing intrinsic causal region and state-selected spacetime are NOT established.**

Parent: `284da0a4af8603d9759a358c9f1b5c6c1a0b7c94`. The frozen Hamiltonian is unchanged. Main runs use a topology-neutral uniform superposition over all Q-edge graphs, without a lattice or topology mask. They are four-, five- and six-event AUXILIARY propagation tests with Q=N−1, waiting time 0.1 and lags 0.1 and 0.3. Connected graphs in this sparse sector are trees. This is not a replacement for the higher-density gauge-capacity/spacetime gate and does not overturn the earlier Q=12 matter-response result.

## 1. The adversarial test needs three protocols

We distinguish an individual source and receiver, and record their distance at the intervention time with orthogonal quantum record states. That operation can disturb a coherent geometry; it is not merely adding a classical column to an unchanged calculation.

We calculate (1) response classified by distance BEFORE intervention, (2) response classified by distance AFTER intervention in the SAME RECORDED experiment, and (3) final-distance response WITHOUT the initial record. The first comparison isolates classification changes. The second isolates recording backaction. The record's energy and information changes are reported below.

The record is an ideal gauge-invariant diagnostic requiring a global graph-distance computation. It is not an autonomous local apparatus derived from CP-12. Its use therefore does not complete the requested physical-preparation gate. The initial pure ray and Q sector also remain supplied.

## 2. Initially disconnected receivers respond

The table gives branch-weighted susceptibilities and conditional intervals after dividing by the initial-distance probability. The source is n₀, the receiver is n₁, and the source operation is exp(−iηn₀). Disconnected means there was no source–receiver path at the moment of the record.

| N | Lag | Disconnected response | Analytic truncation bound | Conditional response interval |
|---:|---:|---:|---:|---|
| 4 | 0.1 | +1.940393e-12 | 5.688e-17 | [+1.94034e-11, +1.94045e-11] |
| 4 | 0.3 | +6.583295e-10 | 3.547e-12 | [+6.54782e-09, +6.61877e-09] |
| 5 | 0.1 | -4.560624e-09 | 8.246e-17 | [-2.66036e-08, -2.66036e-08] |
| 5 | 0.3 | -2.178965e-07 | 5.647e-12 | [-1.27110e-06, -1.27103e-06] |
| 6 | 0.1 | -9.140767e-09 | 9.775e-13 | [-4.19765e-08, -4.19676e-08] |
| 6 | 0.3 | -4.322783e-07 | 8.089e-09 | [-2.02206e-06, -1.94777e-06] |

Bounds are analytically justified upper bounds on omitted-electric-state and quadrature effects, evaluated in ordinary floating point. They are not formal rounded interval-arithmetic certificates and do not bound all floating-point roundoff. Small complete-matrix controls and an independent finite-pulse calculation test numerical consistency. We do not assign physical significance to the tiniest responses solely because their truncation estimates are even smaller.

At six events and lag 0.1, the disconnected conditional magnitude is about 0.624% of the adjacent conditional response. It is much weaker than the nearest response, but larger than the distance-three and distance-four responses in this test. Thus the measured hierarchy is not controlled by initial distance alone. This is not evidence of equal-strength all-to-all transport; it is a resolved shortcut channel. Three small sizes and two times do not establish its thermodynamic survival or rule out state-dependent emergent locality elsewhere.

| Initial distance, N=6 and lag=0.1 | Initial class probability | Conditional interval |
|---|---:|---|
| 1 | 0.333333 | [+6.730292e-06, +6.730308e-06] |
| 2 | 0.273060 | [-1.061485e-07, -1.061363e-07] |
| 3 | 0.127872 | [-3.644628e-08, -3.643634e-08] |
| 4 | 0.039960 | [-8.348354e-09, -8.340707e-09] |
| 5 | 0.007992 | [+2.208891e-13, +7.571205e-12] |
| disconnected | 0.217782 | [-4.197655e-08, -4.196756e-08] |

Conditional source and receiver populations differ among distance classes and are included in the raw records. These are not otherwise identical local environments. A magnitude-versus-distance curve alone would therefore be insufficient evidence for a cone even without the disconnected channel.

## 3. Final-distance bins can greatly overstate the apparent reach

For the six-event disconnected class at lag 0.1:

| Protocol | Weighted response | Analytic truncation bound |
|---|---:|---:|
| Initial-distance classification, recorded | -9.14076748e-09 | 9.775e-13 |
| Final-distance classification, same recorded experiment | -8.56252958e-09 | 3.208e-12 |
| Final-distance classification, no initial record | -8.59177872e-07 | 3.208e-12 |

Both effects matter: reclassification changes the recorded response, while the initial record removes coherent cross-class contributions and changes the state. Treating the first and third rows as two distance labels on the same state would be incorrect. The longest finite-distance bin gives another warning: its early past-class response is extremely small, whereas its final-bin response can be much larger. No causal speed is inferred from that reassignment.

The initially-disconnected shortcut is a result for the INSTRUMENTED state. Recording injects rewiring energy and may itself activate channels by disturbing coherent cancellations. It cannot be assigned unchanged to the unrecorded input. A viable causal test in a candidate spatial phase must derive its physical metric records or establish that readout backaction is small for the relevant observables.

The joint signed history matrix and positive geometry-transition matrix are saved. Row sums reproduce past-class responses; column sums reproduce recorded final-class responses. Their total equals the recorded receiver-density response. Unrecorded total response differs by exactly the cross-record coherence terms.

## 4. Matched-state controls identify the mechanism

The controls retain the same full-law waiting state and distance record, and change terms only during diagnostic post-source evolutions. At five events and lag 0.1:

| Comparator for the initially disconnected class | Response |
|---|---:|
| Frozen full law | -4.56062371e-09 |
| Rewiring removed after the same preparation | -1.35525272e-20 |
| Hopping removed after the same preparation | +0.00000000e+00 |

The no-rewiring result has truncation bound 8.246e-17, compatible with its exact full-law null: without graph changes, disconnected component density algebras cannot signal to each other, even for initially entangled states. Removing hopping also makes the source response exactly zero, while distance still changes with probability 1.263729e-05 in this branch. Geometry changes alone are therefore not counted as a matter signal.

The numerical projection uses a GLOBAL electric-square budget, which could itself couple nominally independent components. The no-rewiring null is checked with its full-law error bound rather than assumed for the projected matrix. This excludes that artifact as the explanation of the reported resolved shortcut.

An independent operational check applies finite source phases ±π/2. Since n₀ is a projector, [f(η)−f(−η)]/[2 sin η] equals the linear susceptibility exactly. Its five-event result agrees with the bilinear computation within 7.115e-19. The finite pulses' conditional energy changes are -0.187946 and +0.233347; they are explicit diagnostic work, not resource-matched preparations.

## 5. Exact size dependence: one edge versus a pool

Let M=N(N−1)/2 and let Z project onto graphs containing none of L specified potential edges. For a nonempty zero-pool sector, the exact full-Hamiltonian cross norm is

    g_L² = ||(I−Z)HZ||²
         = Γ² Q L(M−L−Q+1)/M².

Hopping and diagonal terms preserve graph bits. The upper bound follows from exact source/target degree counts; uniform physical neutral strata attain it. This applies with UNBOUNDED integer electric fields. Nine finite incidence checks verify the counting identities, including a Q=2N case. Large-N rows are evaluations of the proved formula, not many-body simulations.

For a single specified edge and Q~νN,

    g₁² ~ 2νΓ²/N.

For a pool L~λN,

    g_L → Γ√(2νλ).

Thus Γ/M weakens one specified shortcut but leaves finite collective operator capacity for O(N) candidate edges. It cannot, by normalization alone, establish locality. This is a capacity result, not a prediction that every state realizes that norm or a fixed-N Taylor approximation extrapolated to macroscopic times.

A separate full-law speed bound gives P(a specified initially absent edge is present after τ)≤sin²(min(g₁|τ|,π/2)). The saved distance-transition weights satisfy the corresponding weighted bounds. This controls one edge, not the entire change of graph distance or an aggregate causal cone.

## 6. The tested state's coherence actually matters

For the initially disconnected class, the audit computes the instantaneous squared matrix-element norm for creating the direct source–receiver edge. It compares the actual coherent conditional state with complete PHYSICAL-CONFIGURATION dephasing, including graph, fermions and rotors. This is not a simulation of dephased finite-time transport and is stronger than graph-only dephasing.

| N | Coherent coupling square | Configuration-dephased coupling square | Ratio | Removable edge fraction |
|---:|---:|---:|---:|---:|
| 4 | 8.168945e-04 | 8.168945e-04 | 1.000 | 98.027% |
| 5 | 1.192675e-03 | 3.921126e-04 | 3.042 | 98.028% |
| 6 | 1.316644e-03 | 2.178419e-04 | 6.044 | 98.029% |

The dephased state has different energy and entropy; these are saved, so this is not presented as a free coherence benefit. In particular, at N=6 the conditional energies are approximately 2.662 and 3.036 in the projected resource diagnostic. The measured coherent coupling is substantial, rather than merely a worst-case norm unrelated to this state.

For a normalized state with no pool edges, a further exact Cauchy–Schwarz bound is

    ||(I−Z)HZψ||²
       ≤ [Γ² L(M−L−Q+1)/M²] ⟨r_removable⟩.

At L~λN, a sufficient condition for suppressing collective pool coupling is ⟨r_removable⟩=o(N). Gauge constraints make an occupied zero-electric-field edge removable, so electric/matter records can inhibit rewiring. The exact energetic requirement includes

    H_E ≥ (ε_E/2)(Q−r_removable).

At the frozen ε_E=1/4 this is at least 1/8 electric-energy unit per pinned edge. Such energy must be accounted for together with mass, graph cost and hopping. The early states tested here retain about 98% removable edges in this class; they have not reached the sufficient suppression regime. Pinning an arbitrary graph would still not select space.

There is a corresponding finite-time bound, rather than a Taylor extrapolation: starting with no pool edges, its creation probability is at most

    min(1, [Γ² L(M−L−Q+1)/M²] [∫₀ᵀ √⟨r_removable⟩_t dt]²).

Thus a uniformly subextensive removable-edge count over a fixed time window would suppress collective shortcut creation. This supplies a concrete mechanism to investigate inside the existing law; the required state behavior has not been demonstrated here.

## 7. Recording resources and preparation limits

| N | Exact initial mean energy | Record energy change | Record-work truncation bound | System entropy after ignoring record, bits |
|---:|---:|---:|---:|---:|
| 4 | 1.5 | +0.081688456 | 2.334e-11 | 1.685475 |
| 5 | 2.0 | +0.110350843 | 4.846e-11 | 1.985511 |
| 6 | 2.5 | +0.132253876 | 1.650e-08 | 2.139329 |

Initial variance is 2(N−1) and initial system min-entropy zero. Conditional branches remain pure; ignoring the orthogonal record gives eigenvalues p_r and the displayed entropy. A pure apparatus can retain total purity while storing this information. Its memory, graph-distance computation and work are external; no local physical apparatus or state-selection principle is derived.

The record-work bound uses a bounded difference of rewiring operators, avoiding an invalid inference from a norm error to an unbounded-energy error. The full state before recording has conserved mean energy Q/2.

## 8. Verification and the next constructive gate

Fresh checks include the parent chain, 16 dense Hermitian residual tests, an independently built 872-state physical matrix checking the ordered-pair quotient and recording purification, exact edge-pool incidence counts, matched-state term controls, and finite-pulse extraction. Source-code/dependency hashes and hashes of the six saved main evolutions are recorded. The verifier does not silently rerun those expensive records or the state-exposure evolution audit.

This packet does not demonstrate a state-selected spatial phase or a growing intrinsic causal cone. It also does not prove all states fail. It identifies the tested shortcut mechanism and a quantitative route that a successful state would have to satisfy: create spatial organization while suppressing collective far rewiring, potentially through dynamically generated persistent gauge/matter records, without merely freezing an arbitrary supplied graph. The higher-density gauge capacity and a derived physical preparation mechanism remain essential unresolved requirements.
