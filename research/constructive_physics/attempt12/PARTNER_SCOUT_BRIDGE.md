# Compatibility with the newly supplied SCOUT-RG01 and RG02 advances

2026-10-12. The owner supplied the partner's new work during CP12 completion. This note preserves the separate hypotheses; it does not rename RG02 as CP12, alter either frozen Hamiltonian, or claim a combined successful theory.

## Verified sources and executions

RG01: commit 8813abc0228ee572f9d2a580776eb0e55355af5b, `research/reality_framework/relational_matter_gap_scout/LAW.md`. Particle order on a supplied S1 defines adjacency; equilibrium spacing and acoustic waves follow from its declared interaction. It is a useful benchmark for connecting structure and transport, not a derivation of S1 or local gauge symmetry.

RG02: commit bbac3f19ddf7cca12393cc158aad40380b4e969d, `research/reality_framework/port_pairing_scout/LAW.md`, `FULL_JOINT.md`, `REPORT.md`, and associated code. The complete port-fermion Hamiltonian is distinct from its earlier Born–Oppenheimer graph Hamiltonian. Both retain two ports per vertex; the joint law has supplied relative hopping sign, external time and nonlocal reconnections.

The exact committed `shortcut_control.py` was fetched and executed. Its N=4 matrix result |(H²)_fi|=0.025 passes. The committed `full_joint.py` was also executed; its results are retained in PARTNER_EXECUTION.json. Independent random-start eigensolvers recover the N=3 and N=4 ground energies, and a dense N=3 diagonalization checks the first four energies including multiplicity. A single Krylov run may miss repeated eigenvalues; no full-shell counting claim is based on a four-vector sparse run.

## 1. What can be combined scientifically

RG02 provides a much more restricted microscopic connection space and strong finite-size connectivity in its favorable-sign ground state. CP12 provides actual electric-field degrees of freedom, Gauss constraints and back-reaction on allowed connection removal. Both couple the same connections to matter dynamics. These are complementary mechanisms, but their fields, occupations, measures, flux premises and reconnection operators are different. Their numerical successes cannot simply be added.

The earlier bounded-matter fixed-bond entropy theorem and CP12's result also have different hypotheses. CP12 has unbounded compact-rotor electric energy, so its full Hamiltonian is not bounded in operator norm by CN. The relative-entropy argument in THERMAL_THEOREM.md explicitly handles its trace-class electric reference ensemble, yielding a conservative O(1/ln N) exclusion for the specified target classes. Importing the stronger finite-dimensional exponential inequality without controlling rotor weights would be invalid.

## 2. Explicit countertest: Gauss alone does not remove long connections

In frozen CP12, prepare a physical vertex particle at 0 and a hole on edge (0,1), with E_(0,01)=1. Start on a four-cycle. Rewire another neutral edge to the previously absent (0,2). Hop its filled edge fermion to vertex 2, then fill its hole from vertex 0. The charge has relocated to 2, leaving the required electric flux on the new connection; every intermediate state obeys Gauss.

For the selected initial and final states,

 H_fi = (H²)_fi = 0,
 |(H³)_fi| = Gamma t_hop²/M.

The four-cycle test gives 1/60 and unitary T³ coefficient 1/360, identically for electric cutoffs 1,2,3. The construction extends to arbitrarily separated endpoints whenever the required other neutral removable bond exists. Its perturbative order does not grow with the original separation. The per-channel coefficient decreases with N; this is not a proof of macroscopic superluminality. It IS an explicit failure of Gauss constraints to enforce locality of the supplied all-pairs reconnection.

## 3. Reversible intrinsic locality has a component constraint

Consider a graph-changing transition with support S, which changes bonds only among vertices in S. Require S to lie inside one connected component of the PRE-transition graph. A stricter bounded-radius graph-locality rule implies this condition. Such a transition cannot join previously disconnected components: it never touches two of them.

For a Hermitian Hamiltonian the reverse transition exists. If the same intrinsic locality requirement applies to the reverse, the forward transition cannot split a component either; its reverse would join disconnected components. Hence, under these precisely stated conditions,

 **the partition of labelled vertices into connected components is conserved.**

This is a kinematic compatibility lemma, not a universal no-go theorem for local geometrogenesis. It excludes identifying locality solely with each instantaneous graph's connected components while also using those same reversible moves to assemble previously disconnected components. Extra mediating structure, different locality definitions, a prepared connected sector, changing vertex ontology or cosmological boundary conditions would change its premises and must be supplied and evaluated openly.

`local_reconnection_control.py` enumerates all two-bond moves for the RG02 matching kinematics at N=3,4,5. Requiring support in one component on BOTH sides excludes every component-partition-changing move. No new Hamiltonian or ground state is selected by this test.

## Constructive consequence

A successor cannot obtain all three properties merely by combining the current successes:

1. Joining initially disconnected components;
2. Reconnection local solely in the instantaneous graph;
3. Reversible Hamiltonian dynamics with the same locality rule in both directions.

At least one premise needs a physical explanation or refinement. Starting in a connected superselection sector is a legitimate hypothesis but supplies global connectivity. Introducing a mediator with its own accessible neighborhood is another hypothesis but must explain that neighborhood. Allowing microscopic shortcuts is also possible, but the actual detector response must show that they become irrelevant in a controlled continuum limit.

This sharpens the next physical test: establish an operational continuum causal bound for joint matter and reconnection, rather than assuming that gauge invariance or favorable connectedness already supplies it. Neither scout nor CP12 currently earns that bound, and none of these results establishes a gravitational spin-2 mode.
