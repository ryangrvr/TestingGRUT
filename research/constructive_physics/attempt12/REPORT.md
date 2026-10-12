# CP12 — conserved connections, dynamical compact gauge fields and physical back-reaction

**Status:** a tested finite quantum gauge–matter–geometry construction; NOT a theory of gravity or a selected physical spacetime. The conserved resource removes the dense-edge pathology. The full thermal analysis still excludes the specified higher-density lattice-template classes at fixed positive temperature. The new partner scouts are incorporated as separate constructions in [PARTNER_SCOUT_BRIDGE.md](PARTNER_SCOUT_BRIDGE.md).

## What was actually constructed

One new law was frozen before execution, in local commit `74d1ca7098ce53a93df9d9af997b1435bd4a7759`. Remote publication of that identical law occurred later at `e48524cc40b446b5b27e82eea74e5d586d09840c`; this is not a claim of pre-execution remote preregistration. The research branch is rooted at the CP11 partner magnetic review `2bca657df1a3539449f0468d06280a9b56e079a8`.

[LAW.md](LAW.md) supplies a conserved number Q of undirected connections, vertex and edge fermions, compact U(1) half-link rotors, electric energy and Gauss constraints. Present connections determine fermion hopping. A positive connection-rewiring Laplacian changes geometry only where the removed connection has zero electric field and its neutral filled edge mode. Absent connections cannot carry gauge excitations. No lattice, valence cap, preferred dimension, optimized phase or magnetic plaquette reward is included.

The frozen numerical controls are m=t=1, kappa=1/2, epsilon_E=1/4, Gamma=1/10. External time, quantum mechanics, gauge symmetry, background edge charges, resource sector and all-pairs rewiring are supplied assumptions. The model does not explain those ingredients. Fixed Q is an actual conserved operator of its dynamics, but the existence and density of that resource are hypotheses, not discoveries about nature.

The exact partner magnetic script was executed and reproduced its static-phase results. CP12 then quantized the gauge sector instead of choosing pi flux. Static optimized-flux energies must not be interpreted as CP12's joint quantum energies.

## 1. Exact constraints and extensive stability survive

Every implemented basis state satisfies all vertex and edge Gauss constraints, fixed fermion number, fixed Q, and frozen inactive-edge conditions. A separate brute-force electric-field enumeration checked the constructed basis on four small controls. Every hop and allowed reconnection stays inside the physical subspace. Hermiticity and signed simultaneous-permutation covariance were checked independently.

The Hamiltonian obeys the same magnetic-incidence lower bound, now with positive electric and reconnection terms:

\[
 H_{12}/N\ge\kappa\nu-\sqrt{m^2+2t^2\nu}+m\ge-1/4,
 \qquad \nu=Q/N,
\]

at the frozen controls. See [PHYSICAL_DERIVATIONS.md](PHYSICAL_DERIVATIONS.md). No large-N negative-energy collapse is hidden in the new gauge dynamics.

Because [H,Q]=0, the connection fraction is exactly Q/binomial(N,2), at any temperature and in every state of that sector. For Q proportional to N it tends to zero. This fixes CP11's dense-edge thermal failure within a declared resource premise; it does not select dimension or connectivity.

## 2. Full joint gauge–matter–geometry ground states

All labelled graphs at each Q, all compatible fermion states and the regulated Gauss-law electric states are retained in the full calculations. The cutoff is |E|<=L; projected rotor raising is not unitary at its boundary. It is a numerical regulator, not a modified claim about the compact gauge group.

| N,Q | Full dimension at L=3 | Ground energy | P(connected) | Uniform-graph P(connected) | P(single cycle) | Uniform-graph P(single cycle) |
|---|---:|---:|---:|---:|---:|---:|
| 4,3 | 1,032 | -0.508218868 | 0.882530 | 0.800000 | 0 | 0 |
| 4,4 | 6,075 | -0.598552282 | 1 | 1 | 0.289561 | 0.200000 |
| 5,5 | 427,986 | -0.728493364 | 0.949037 | 0.880952 | 0.094871 | 0.047619 |
| 5,6 | 2,866,255 | -0.782952055 | 0.992195 | 0.976190 | 0 | 0 |

For N=4,Q=4, connectedness is automatic from graph counting; it is NOT an earned dynamical result. The nontrivial comparison there is cycle versus triangle-with-pendant. For Q!=N a single spanning degree-two cycle is impossible, so zero in that column is kinematical.

The N=5,Q=6 matrix has 29,878,495 nonzero entries. Changing L=2 to L=3 changes its ground energy by about 4.3e-10. Its L=3 boundary probability is 6.1e-11 and the recorded eigen residual is 8.3e-11. These are convergence controls, not rigorous infinite-cutoff error certificates.

A vacuum-relative permutation sector makes N=6,Q=6 accessible: 53,333 orbit states at L=3 represent a sector of a 36,801,265-state physical cutoff space. Its lowest energy is -0.858318930, connected probability 0.866971 (uniform 0.731269), and single-cycle probability 0.032202 (uniform 0.011988). **This is a sector result, not a proved full N=6 ground state.** The full and reduced N=4,N=5 ground energies agree where both were calculated.

[validate_reduction.py](validate_reduction.py) independently constructs the orbit isometry W and checks W†W=I, W†HW=H_reduced and HW=W H_reduced, with errors below 2e-14. The sector uses the signed fermion representation normalized relative to the filled-edge vacuum; treating event permutations as a physical gauge redundancy would be an additional interpretation, not a result.

## 3. Dynamical flux and hostile geometry controls

The gauge fields are not frozen at an optimized angle. Gauge-invariant Wilson-loop expectations in the ground state are weakly negative: the conditional triangle expectation is about -0.0456235 for N=3,Q=3; the selected square expectation in the full N=4,Q=4 state is about -0.0110196. This is a dynamically generated correlation favoring a negative cosine component, NOT a classical pi-flux state (which would give -1).

At common N=4,Q=4 and with reconnection disabled only for a declared conditional comparison, C4 has energy -0.622125877, versus -0.591167719 for the triangle-with-pendant. This compares full gauge and fermion dynamics on each graph, not a hand-chosen flux. Degree regularity and loop effects are both involved; the calculation does not isolate flux as the sole cause.

A hostile comparison remains: conditional long cycles approach energy per vertex -0.15551810, while disconnected triangles have -0.15576652 per vertex. Thus loops alone do not prove a single large spatial component is best. The full reconnection term is extensive and may alter that comparison; this is not a proof that the full CP12 ground fragments.

## 4. Physical records affect the geometry, without a global causal-role oracle

A gauge-invariant particle-hole/electric excitation on one present connection changes its removability. Independent evolution gives

\[
 P_{\rm remove}(T;\mathrm{neutral})
 =(M-Q)(\Gamma/M)^2T^2+\cdots,
\]
\[
 P_{\rm remove}(T;\mathrm{excited})
 =\tfrac14(M-Q)(\Gamma/M)^2t_{\rm hop}^2T^4+\cdots.
\]

The excitation must be erased before that connection can be removed by this law. This is a mechanical back-reaction effect, supported at two electric cutoffs with energy, norm and reversal checks. It does not require consciousness, irreversible collapse or an observer with access to a global geometry oracle. It also does not yet construct an autonomous observer, and it depends on the supplied neutral-edge rewiring rule.

The same electric constraints supply a bare separation diagnostic. For a vertex particle at i and an edge hole at e,

\[
 \min E_{\rm diagonal}-\kappa Q
 =2m+\frac{\epsilon_E}{2}d_{\rm incidence}(e,i).
\]

Integer flux along a shortest path saturates the lower bound. This was checked independently for **495 physical charge placements**. It is an instantaneous constrained energy relation; screening, full interacting spectroscopy, relativistic causality and Newtonian gravity have not been obtained.

## 5. Preparation and thermal tests do not equal ground-state selection

A neutral equal-amplitude graph preparation was evolved unitarily without encoding the ground state. At N=4,Q=4 its cycle probability rises from 0.2 to 0.262418715 at T=20. Cutoffs L=4,5 agree within 1.6e-8. Its energy remains 2, whereas the ground energy is about -0.59855. The simulation therefore demonstrates generated correlations, not cooling or convergence to a chosen ground geometry.

Full finite-temperature traces were calculated for N=4,Q=4 at L=2,3,4,5. At L=5, cycle probability is 0.217807417 for beta=5 and 0.296509879 for beta=20. These agree with L=4 to displayed precision. At beta=1 the L=4-to-5 change is 4.3e-7; at beta=0.5 it is 8.2e-5 and L=5 boundary probability remains 4.8e-4. The hottest case is **not claimed fully cutoff converged**. Good ground-state convergence did not justify treating the thermal cutoff as converged.

The analytic large-system result is stronger than a small thermal scan: [THERMAL_THEOREM.md](THERMAL_THEOREM.md) includes the full noncommuting quantum gauge and reconnection terms through relative entropy. For Q=nu N with nu>1, fixed beta>0, and a target class with at most N! exp[O(N)] labelled graphs,

\[
 P_\beta(\mathrm{target})=O_\beta(1/\ln N)\longrightarrow0.
\]

The class assumption covers fixed lattice sequences and a bounded exponential number of variants, not every imaginable curved spacetime discretization. This theorem does not exclude nu=1, zero temperature, non-Gibbs cosmological preparation or different microscopic state spaces. Nevertheless, CP12 does not solve higher-density equilibrium spacetime selection merely by imposing sparsity and adding gauge fields.

## 6. The partner advancements sharpen the locality obligation

The owner supplied SCOUT-RG01 and SCOUT-RG02 while this packet was being completed. Their exact laws were retrieved and the RG02 complete joint and shortcut scripts executed verbatim. The reproduced favorable-sign connected probabilities are 0.995021, 0.985453, 0.968941 at N=3,4,5. These are meaningful separate results; they are not CP12 outputs.

A new CP12 countertest shows that Gauss constraints alone do not fix RG02's underlying shortcut issue. One neutral reconnection followed by two gauge-invariant fermion hops relocates a physical charge, with

\[
 |(H^3)_{fi}|=\Gamma t_{\rm hop}^2/M,
\]

independent of the original separation in the corresponding family of configurations. The N=4 explicit matrix value is 1/60 at all three tested cutoffs. The amplitude is size-suppressed; no observable superluminal signal has been proved. But locality is not enforced by Gauss alone.

Conversely, a reversible graph move supported inside a single connected component on BOTH sides cannot change the component partition. All RG02 matching transitions at N=3,4,5 were enumerated: the condition excludes every merging/splitting transition. This is a conditional compatibility lemma, not a prohibition on emergent space. It exposes why simply restricting all reconnection to existing graph-local neighborhoods can prevent disconnected pieces from ever assembling.

See [PARTNER_SCOUT_BRIDGE.md](PARTNER_SCOUT_BRIDGE.md) for assumptions, exact tests and the distinction between its bounded-Hamiltonian entropy theorem and CP12's unbounded-rotor analysis.

## Disposition and scientifically useful next target

Retain CP12 as a consistent finite compact-gauge construction with conserved connection resources, explicit matter back-reaction, dynamical flux correlations and a bare physical-distance diagnostic. Its tested finite-size connectedness is encouraging. Its specified canonical thermal spatial-selection target fails under the theorem's assumptions, and its reconnection has not earned operational locality. No continuum, Lorentzian sector, graviton, Einstein equation or new empirical prediction is established.

The combined work now identifies an exact question for a successor: can a physically prepared matter/reference sector make graph-changing shortcuts irrelevant to detector responses in a controlled large-system limit, while retaining the process that builds connected structure? Merely imposing local rewiring or gauge symmetry cannot be credited with solving both. This is a prospective physical target, not an altered CP12 law or an opened CP13.

## Prior work and reproducibility

The use of dynamical compact gauge variables and Gauss constraints belongs to Hamiltonian lattice gauge theory: Kogut and Susskind, Phys. Rev. D 11, 395 (1975), https://doi.org/10.1103/PhysRevD.11.395. Physical Hilbert spaces and cutoff issues have established literature, e.g. Kaplan and Stryker, https://arxiv.org/abs/1806.08797. Lieb's flux-phase theorem, https://arxiv.org/abs/cond-mat/9410025, concerns a different Hamiltonian and does not select CP12's flux. The graph-entropy difficulty has important prior work by Si Chen and Steven S. Plotkin, https://arxiv.org/abs/1210.3372. No priority claim for the general mechanisms is made.

The proofs here are author derivations pending external review. Numerical checks are code/theory controls, not measurements. [REPRODUCE.md](REPRODUCE.md) lists exact commands and the different scopes of full, sector, conditional, thermal and dynamical calculations. All previous frozen packets and Stage-3 Card 3 remain unchanged.
