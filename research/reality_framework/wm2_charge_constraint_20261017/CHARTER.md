# WM2-UC1: charge constraints and the full matter algebra

Base: `090c2ff2fd4b7cdf6428d0dbb8d7c4e73fc5b35e`.
Date: 2026-10-17. One separately named architectural construction.

## Frozen question

Can explicit preserved microscopic constraints exclude off-bond transport
for two distinct matter sectors, including dressed/virtual processes,
without silently assuming that every allowed interaction follows the graph?

Preserve WM2, S9's author-derived status and the closed numerical branch.
No formation solver, six-event run, gravity or dimensional search belongs
to this campaign. Do not claim independent review from Builder's own work.

## Candidate, supplied before diagnostics

WM2-UC1 keeps the single-end word grammar, storage and physical receiver.
It replaces the original neutral probe by two number-conserving sectors:
A with integer gauge charge 1 and two internal levels; B with charge 2
and three internal levels. There is at most one carrier of each species.
Explicit immobile root compensators cancel their total charge; their
occupations equal the conserved particle numbers.

Add one oriented electric register E_ij=-E_ji per event pair, with seven
levels -3,...,3. Orientations are conventions, not pair affinities. Let
B_ij be the occupied consecutive-bond projector of the generated word.
Impose commuting constraints

    G_i=sum_j E_ij - n_Ai - 2 n_Bi + R_i(f_A+2f_B)=0,
    C_ij=(1-B_ij)E_ij^2=0,
    N_A=f_A, N_B=f_B.

The exact physical space is their joint kernel. Zero electric energy is
a declared choice, not a gap claim. Gauge-invariant append/pop carries a
zero-flux projector on the changed bond. Gauge-dressed charged hopping
is permitted when it preserves these constraints. No global P h P
shortcut is admitted as a primitive implementation.

The vacuum seed adds E=0 and empty compensators/carriers to the original
WM2 preparation. Prove whether its formation restriction is exactly the
old law rather than assuming inheritance.

## Allowed operator class U

Every allowed Hamiltonian preserves the displayed constraints and matter
numbers, is covariant under event relabeling, and has bounded primitive
register support independent of N. Geometry-changing terms remain the
gauged WM2 single append/pop operations and their adjoints. Other terms
are geometry-diagonal and touch at most r electric registers per summand,
with fixed r. The total changing interaction norm is uniformly bounded
on the declared zero/one-per-species charge sectors. Configuration-preserving
diagonal energies and receiver motion are treated exactly.

The class includes internal, density and inter-species operators when
they satisfy these conditions. Do not exclude a constraint-preserving
product of distant internal observables merely because it is inconvenient.
Do not narrow U after seeing its counterexample. A theorem for charged
position transfer alone is not a theorem of universal matter causality.

## Mandatory attacks

1. Compute commutators/leakage for the old bare all-pairs and star terms,
   their single-link gauge dressings, and additions of local hopping.
2. Include Wilson strings and all complete primitive words. Count both
   electric-register changes and intermediate prefix removal. Pay the
   full norm tail, not the shortest path alone.
3. Exhibit nonzero nearby motion for both distinct sectors. Preserve all
   conditional receiver states and reverse geometry transitions.
4. Test internal matter interactions explicitly. An anticipated failure
   is a gauge-invariant, root-blind coupling of A's internal observable
   to B's distant internal observable, with bounded norm.
5. Price separately: occupied-bond rule; constraints enforcing it; any
   mechanism selecting the constraints/physical sector. No third item
   follows from declaring the first two.

## Outcomes and stop

Produce the explicit operators, a proof or failure of invariance, the
two-species Hamiltonian, complete word analysis, witness results and a
supplied/derived ledger. If a neutral/internal channel survives, retain
the charged-motion theorem at its exact scope and record the full gate
as failed. A commutant-closure obstruction is an acceptable negative
result; do not repair it by silently imposing connected graph support.

This freeze follows analytical design and anticipated failure modes,
before new matrix controls. It is not blind external preregistration.
All mathematical results remain author-derived pending independent review.
