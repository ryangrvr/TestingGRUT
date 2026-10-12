# What the frozen CP12 matter Hamiltonian can select

Author derivation, 2026-10-12. No new action, coupling, candidate, or physical law.
This separates fixed-graph energy information from full graph dynamics.

## 1. Neutral vacuum hopping moments

Fix a graph G with Q edges, degrees d_i, and neutral filled-edge vacuum |0_G>.
Write V for the hopping operator with its hopping coefficient removed. Let

    W(G) = sum_i binom(d_i,2),    Delta = 2m + epsilon_E/2.

V|0_G> contains 2Q orthogonal physical charge-hole/one-electric-unit states.
Each has energy Delta above the neutral vacuum. Thus <V^2>=2Q.

At the second application, hopping can undo its first action, or excite another
edge into a different vertex. Two edge fermions cannot occupy the same vertex;
nor can one filled edge be emptied twice. The number of compatible unordered
two-hop choices is

    K = binom(2Q,2) - Q - W(G).

The two orders of compatible even-fermion hopping operators add coherently.
The final electric fields distinguish different incidence choices. Each of K
two-excitation states therefore has amplitude of magnitude 2. The vacuum
return amplitude is 2Q. This proves the exact moment

    <0_G| V^4 |0_G> = (2Q)^2 + 4K
                    = 12Q^2 - 8Q - 4W(G).

This is an exact finite-graph algebraic statement for the hopping part, not a
claim that the full Hamiltonian has this fourth moment. It survives a hard
electric cutoff L>=1 because every intermediate field needed here is 0 or 1.

## 2. A derived degree regularizer, with a controlled domain

For this subsection only, set graph rewiring aside as a declared fixed-graph
comparator. Retain the CP12 mass and electric terms. Let the hopping strength
be t_hop and use perturbation theory about t_hop=0. The neutral vacuum is the
unique zero-excitation ground state at fixed finite graph when m,epsilon_E>0.
All one-excitation intermediates above have gap Delta; the compatible
two-excitation intermediates have gap 2Delta. The second-order and fourth-order
Rayleigh-Schrödinger corrections are therefore

    E(G) = kappa Q - 2Q t_hop^2/Delta
         + [4Q + 2W(G)] t_hop^4/Delta^3 + O_G(t_hop^6).

The fourth-order coefficient is obtained as (2Q)^2/Delta^3 minus
4K/(2Delta^3). The subscript G is deliberate: this note proves no uniform
thermodynamic bound on the remainder. At the frozen t_hop=1, Delta=2.125,
the fourth-order truncation is quantitatively inaccurate and is not used
as a replacement for exact diagonalization.

At fixed Q, W=(sum_i d_i^2 -2Q)/2, so the derived positive term penalizes
degree variance. This interaction arises from Pauli exclusion and virtual
electric/matter excitations. It was not added as a geometric energy penalty.
It does not pick a spatial dimension or graph connectivity.

Example: at Q=N, both a spanning ring and a disjoint union of cycles have
d_i=2 and W=N. They remain degenerate through this order. At Q=N-1, degree
regularization alone also allows a path component accompanied by cycles,
not just one spanning path. Regular degree is insufficient for space.

`matter_selection.py` checks the hopping moment, evaluates the perturbative
coefficient directly using the diagonal reduced resolvent, and compares
with exact ground energies at small hopping. Varying hopping here is a
mathematical coefficient check, not a retuning of the frozen construction.

## 3. The complete graph-block comparator

When comparing with full CP12, the graph compression is

    P_G H P_G = H_matter,G + a_N sum_{e in G} R_e,
    a_N = Gamma (M-Q)/M,

where R_e projects onto a neutral removable edge. This diagonal is retained
by `compressed_competitors.py`. It changes the actual energy hierarchy; one
must not substitute the no-rewiring comparator from section 2 for this block.

At N=6,Q=6, a_N=.06, electric cutoffs 3 and 4 give:

    one spanning ring:  E_block = -0.635203537847...
    two triangles:     E_block = -0.637262028038...

The disconnected competitor has lower block energy by about .00205849.
With the rewiring term entirely omitted, the respective energies are
-.933108967404 and -.934599106378; fragmentation wins in that comparator too.
These are converged finite-matrix results, not rigorous infinite-cutoff error
bounds. Cutoff differences and eigenvector residuals are recorded.

There is a narrow variational consequence. No nontrivial single-edge swap
maps one 2-regular simple graph to another: preserving all degrees would
require the replacement edge to have the same unordered endpoints. Thus
compression to the collection of spanning rings is block diagonal in graph.
At N=6 all rings have the same block spectrum by relabelling/gauge convention.
The two-triangle trial state therefore has lower energy than **every state
supported entirely on spanning rings**, in the tested regulated Hamiltonian.
This does not bound the ground-state ring probability away from a value close
to one: virtual mixing with other graphs can lower the energy, and a full
ground-state computation would be required. It also says nothing by itself
about a nonequilibrium spatial phase.

## 4. Full dynamics breaks the neutral quotient

The exact full evolution in `evolve.py` uses all fixed-Q graphs, all Gauss-law
fermion states and all electric states within the stated cutoff. It does not
reuse the neutral-only k-sector reduction.

At time T, let p_G be the graph probability and p_k the probability of reference
degree k. Define the neutral-quotient-matched comparator

    u_G = p_{k(G)}/n_{k(G)},   D_TV = (1/2) sum_G |p_G-u_G|.

Neutral-only evolution from the supplied initial coherent state has D_TV=0
exactly. Full CP12 develops D_TV>0. Comparing a target probability against
sum_{G in target} u_G holds the reference-degree distribution fixed and tests
shape information beyond attachment. u_G is a diagnostic distribution, not
a second quantum evolution or a physical preparation.

For N=5,Q=4 the directly calculated graph-probability Taylor coefficients
remain k-uniform through fourth order and first break that uniformity at
sixth order. The L1 norm of the sixth-order coefficient outside the uniform
sector is .00087702222223 (double precision control). This is a finite-system
computed coefficient, not an arbitrary-N theorem. Odd orders vanish for the
real Hamiltonian, real initial amplitudes, and diagonal graph readout.

The fourth-order vacuum energy preference and sixth-order probability change
are different quantities. The former offers a microscopic interaction to
investigate; it is not a derivation of the finite-time distribution from a
Gibbs state or an adiabatic approximation.

## 5. Physical boundary

Full matter back-reaction can remove the uniform graph weighting responsible
for the neutral-only counting obstruction. That removes a premise of that
particular theorem; it does not invalidate it or prove spatial assembly.
Paths and rings are transparent finite-size targets, not sufficient tests of
relativistic spacetime. The preparation still supplies graph coherence and a
reference vertex. No large-N geometry selection, Lorentz limit, gravity, or
autonomous preparation is proved here.


## Exact first topology-sensitive energy coefficient

A further conditional calculation compares the bare fixed-graph vacuum energy
of a six-cycle with twice that of a triangle. Rewiring is omitted in this
comparator only; it is not a change to CP12. Write x for the hopping coefficient.
The extensive zero-hopping terms cancel. Exact rational perturbation theory gives

E(C6)-2E(C3) = (47080992982591078400/310998849352074430281) x^12
 - (899838970108087191238151766016/662798998370613420396196869375) x^14
 + O(x^16).

All coefficients through order ten vanish identically. The leading coefficient
is positive, approximately .15138638963030954. Thus sufficiently weak hopping
favors the disconnected triangles, despite identical degree sequences.
The negative fourteenth-order coefficient is about -1.357634776636958: do not
extrapolate this truncated series to the frozen hopping x=1.

`loop_series.py` uses intermediate normalization with vacuum |0>, diagonal
excitation energy D, and R=D^{-1} on its orthogonal complement:

    E_n = <0|V|psi_(n-1)>
    |psi_n> = -R [ V|psi_(n-1)> - sum_(a=1)^n E_a|psi_(n-a)> ].

All operations in the recursion are rational; the frozen diagonal gaps are
integer multiples of 1/8 and hopping matrix elements are integers. Electric
cutoffs L=2 and L=3 give identical differences through order14. This explicitly
checks regulator independence across those two cutoffs, not every cutoff.
`loop_numeric_control.py` separately compares sparse ground-state eigenvalues
at small hopping against the series. A closed triangle can support a six-hop
flux excursion and its return; this suggests why the first energy difference
appears at order12, but the calculation does not isolate all contributing paths.
