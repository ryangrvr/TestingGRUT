# CP12 frozen law — conserved connections with compact U(1) transport

2026-10-11. Owner explicitly requests one combined resource, fermion and dynamical gauge construction and at least a sustained 40-minute investigation. Root: CP11 flux review 2bca657df1a3539449f0468d06280a9b56e079a8. This is one NEW candidate; no changes to CP01–11 or Stage3 Card3.

## Microscopic inputs and physical Hilbert space

N vertices; M=N(N-1)/2 possible undirected connections with occupation n_e∈{0,1}. Work in a sector Q=Σ n_e of conserved connection number. No local valence, dimension, embedding or lattice is prescribed. Q is a supplied resource/charge, not derived. Thermodynamic targets Q/N=ν=1 and 2 are distinct declared charge sectors of the same Hamiltonian, not a fitted density. Small controls N=3,Q=2,3 and N=4,Q=3,4; denser sectors only if feasible.

There are canonical vertex fermions a_i and possible-edge fermions b_e, total fermion number M. A present edge e=(i,j), i<j, has TWO compact U(1) rotors U_ie, U_je, with integer electric fields E_ie,E_je and [E,U]=U. Both legs are oriented from edge center to vertex. Absent edges are restricted to E_ie=E_je=0 and b_e†b_e=1; they carry no independent excited gauge states. Present edges have no such restriction.

Physical states satisfy Gauss constraints at each vertex and edge center:

 G_i=Σ_{e contains i} E_ie-a_i†a_i=0,
 G_e=-E_ie-E_je-b_e†b_e+1=0.

The background +1 at every edge center and total M particles are supplied, with zero total gauge charge. Vertex and edge gauge transformations, compact U(1), Fermi statistics, external time and canonical quantum mechanics are inputs. This is NOT derivation of gravitational constraints.

## Frozen Hamiltonian

 H12 = κQ + 2m Σ_i a_i†a_i + (ε_E/2) Σ_{e,i∈e} E_ie²
       + t Σ_{e=(i,j)} n_e[(a_i† U_ie-a_j† U_je)b_e + h.c.]
       + (Γ/M) Σ_{allowed undirected rewirings r} L_r†L_r.

Equivalently the mass is mΣ a†a -mΣ b†b +mM in the fixed-M sector. Controls m=t=1, κ=1/2, ε_E=1/4, Γ=1/10. M=0 sector has no rewiring term.

A rewiring removes one occupied edge e and activates one unoccupied edge f, leaving all fermion and electric variables unchanged. It is allowed only if e is in its removable state b_e†b_e=1 and E_ie=E_je=0; f already satisfies these conditions while absent. Each unordered pair of physical basis states is counted once, L_rψ=ψ_final-ψ_initial. Hence graph rewiring is a positive Laplacian and preserves Q, total fermion number and every Gauss constraint. The 1/M rate is supplied; allpairs rewiring remains a possible operational-nonlocality problem. It is not called relativistic transport.

The electric term provides genuine compact-gauge quantum dynamics. No magnetic plaquette term or favored-loop reward is supplied. Any flux preference must follow from fermions coupled to the electric dynamics. Link phases must not be replaced by a preselected π flux in the physical test.

## Execution and failure tests

First rerun the exact partner flux script. Freeze this law before CP12 executions. Derive gauge invariance, allowed rewiring, energy bounds, component charge neutrality, conditional propagation, and canonical finite-temperature counting. Show exactly what finite Q fixes and what it does not.

Construct the Gauss-law physical basis directly. Use electric cutoff |E|≤L only as a regulator, with L=1,2,3 and larger if needed; projected U must be acknowledged as nonunitary at the boundary. Retain exact gauge invariance. Compare cutoffs and boundary probabilities; no continuum claim from a single cutoff.

Solve joint connection+fermion+gauge systems for the declared small sectors. Measure connectedness, energy, flux Wilson loops, and actual electric excitations. Compare conditional cycle and hostile geometries at common Q. Test finite quantum electric dynamics against the supplied static-flux comparator without identifying their energies.

Analyze large-N fixed-Q counting, including graph entropy and physical gauge-state weights. Test whether a spatial-family advantage survives the common ensemble; do not declare geometry selected because the graph is sparse. Static lattice/flux integrals are benchmarks, not full dynamical solutions. Failure does not authorize changing couplings or silently adding a plaquette term. Report limits and continue useful in-scope derivations.
