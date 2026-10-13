# CP-12: connection back-reaction, physical sectors, and the missing gravitational response

Research branch: `codex/grut-cp12-gravity-response`
Parent: `553aeffc10faef6928a7646d4a9173019273e2d7`.

## Result

The unchanged CP-12 Hamiltonian produces an exact matter-driven connection response. Larger symmetry-reduced calculations also identify collective connection excitations. Neither result establishes a gravitational field: their size scaling, symmetry, and operational interpretation must be distinguished from a propagating spin-two metric.

This packet supplies derivations, independent executable controls, and a conventional gravitational benchmark. It changes no frozen law and opens no new constructive-physics candidate.

## An exact local response

Prepare the neutral sea on one definite graph, intervene with the frozen vertex hopping operator normalized by the square root of its initial degree, and read a connection occupation. Its first nonzero retarded response is cubic in time. For an occupied edge incident on the source it is

`-Gamma^2 t_hop^2 (M-Q) T^3 / (M^2 sqrt(d_i))`.

For any absent edge it is

`+Gamma^2 t_hop^2 sqrt(d_i) T^3 / M^2`.

The second coefficient does not depend on the absent edge's distance from the source. The sum over all connections vanishes, as required by conserved connection number. An independent integer-state builder verifies these identities and Gauss constraints without an electric cutoff. Complete finite-system evolutions and independent source impulses confirm the response.

At extensive connection number and bounded source degree these coefficients decrease with system size. Stronger full-Hamiltonian bounds control individual connection readouts beyond this short-time expansion. They constrain this particular metric interpretation; they do not exclude every collective encoding of geometry.

## Preparation and ground state are different physical regimes

The neutral rewiring operator has an exactly solvable Johnson-graph spectrum. A definite sparse graph has almost all its spectral weight in the highest graph harmonic, yielding an all-time neutral-sector departure bound of order `1/N`. This statement excludes matter hopping; it must not be applied to full CP-12 evolution.

The full-law harmonic calculation sharpens that distinction. Locking an occupied neutral connection preserves the highest harmonic; unlocking leaks a precisely calculable fraction. Combining this algebra with the parent fixed-block lower bound and a full-Hamiltonian variational upper bound gives a ground-state highest-harmonic weight bounded by order `1/N`. A definite graph instead has weight `1-O(1/N)` there. Stable prepared geometry and the ground-state ensemble therefore occupy very different sectors.

An ideal graph record also has different energetics in these sectors. Within the highest harmonic, dephasing the graph changes the system energy by `-Gamma <r>/M`, a small nonpositive shift. This is a system-energy statement; the apparatus and total energy exchange are not modeled. Projecting operators into this harmonic generates noncommutation, but the controlled local commutator bound vanishes with size. It is not an emergent gravitational observable algebra.

## Tensor response exists, but permutation symmetry is not Lorentz spin

A supplied tetrahedral propagation comparator shows that four mean bond weights miss two independent shear components. Existing electric-activity pair correlations provide those components without adding operators to the Hamiltonian. An exact, unregulated four-arm graph block has nonzero doublet and triplet responses, with finite excitation gaps. Complete fluctuating five-event calculations retain both responses and agree between electric cutoffs two and three.

Their lowest common fivefold shell transforms as the event-pair representation `[N-2,2]`. Its multiplicity is an event-permutation property, not evidence for five polarizations of a graviton. The larger calculation makes this distinction especially clear:

| Events | Connections | Pair-shell multiplicity | Pair gap | First permutation-singlet gap |
| --- | --- | --- | --- | --- |
| 5 | 4 | 5 | 0.06847112 | 0.15723255 |
| 6 | 5 | 9 | 0.06938808 | 0.16837800 |
| 7 | 6 | 14 | 0.07012838 | 0.17603917 |

The multiplicity follows `N(N-3)/2`. These finite gaps do not establish their thermodynamic limit. The connection density here is `Q=N-1`, not the proposed four-regular spatial sector.

The seven-event result is an exact signed symmetry reduction of a cutoff-two physical space containing 161,804,041 basis states to 761,244 states. We did not diagonalize the entire space. The subgroup admits singlet, standard, and pair representations; other representations were not searched for a lower ground energy. The five- and six-event electric-cutoff comparisons are excellent: the six-event cutoff-two/three ground energies differ by approximately `8.1e-11`, and pair gaps by `1.4e-11`. Seven events has only one completed cutoff. Its cutoff-three attempt built the basis but did not produce a completed spectrum and is not counted as a passing test.

Independent intertwiners against the frozen Hamiltonian builder validate the signed quotient at four and five events. Complete integer symmetry checks are used for moderate matrices; the largest result uses 64 deterministic columns plus numerical spectral residuals. These are numerical controls, not interval-certified thermodynamic results.

## The physical target

A separate conventional spin-two benchmark verifies four first-class constraints, two transverse-traceless propagating modes, and the universal-coupling condition of the soft-emission gauge identity under its stated Lorentzian assumptions. Those structures are supplied in the benchmark. CP-12 has not derived them.

The decisive remaining calculation is a relational, low-energy response with a controlled spatial continuum, an appropriate energy constraint, two tensor propagation modes, and universal matter coupling. A permutation multiplet, an electric Gauss law, or a finite tensor-shaped response cannot substitute for those requirements. Geometry selection against thin, fragmented, and nonspatial competitors also remains unresolved.

## Verification and partner entry points

`verify.py` passes ten groups. It validates recorded controls and recomputes the outward ground-sector certificates; it is not a rerun of every large eigensolve. `DERIVATION.md` contains the analytical proofs and primary-source distinctions. `MANIFEST.json` records input hashes, file hashes, numerical scope, and reproduction commands. Individual scripts and JSON outputs are retained for adversarial review.

The achieved result is an exact connection back-reaction and a clearer separation between prepared geometry, collective label symmetry, and gravitational response. A working theory of reality has not been established. The frozen CP-12 mechanism now has concrete physical tests it must satisfy rather than another inferred continuum claim.
