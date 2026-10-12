# CP12: full coherent evolution, derived degree bias, and a fragmentation control

2026-10-12. Additive follow-up based on `7452a462ba787e5b04e01a4ca4b92dc68a84fb6e`.
The frozen CP12 law is unchanged. No CP13, Stage-3 card, canonical branch or PR
was modified. These are author calculations and reproducibility controls, not
independent experimental validation or a completed theory of spacetime.

## Outcome

The full fermion/electric dynamics does **not** preserve the uniform graph
weighting used by the neutral-only no-assembly theorem. It generates additional
shape information. Coherence-assisted attachment also survives in small full
systems, although it is markedly slower than the neutral comparator.

A microscopic source of graph dependence can be derived: the fourth-order
fixed-graph vacuum energy penalizes unequal degrees through Pauli exclusion.
But that term does not distinguish a single ring from disconnected cycles.
At the frozen couplings an explicitly tested disconnected competitor actually
has lower graph-block energy than the spanning ring.

Thus there is a genuine matter-induced geometric bias, together with a concrete
reason not to interpret it as macroscopic spatial formation.

There is also a target correction: the supplied initial distribution already
typically has a giant nonspatial random-graph component. Failure to connect
every vertex is not failure to have a macroscopic component. See
[PREPARATION_TARGET.md](PREPARATION_TARGET.md) for the primary literature and
an independent N=64–4096 sampling control. The missing achievement is its
spatial organization and operational propagation, not just its existence.
This packet does not compute a relational retarded detector cone for the
coherent initial ensemble. The earlier definite-graph propagation theorem
cannot automatically be transferred to it.

## Verbatim partner audit

All three committed scripts from the parent audit were fetched and executed
without edits: `check.py`, `full_law_moments.py`, and `global_connectivity.py`.
All passed. Source SHA256 hashes and their actual JSON outputs are recorded in
`partner_verbatim_results.json`. In particular the full-law one-step moment
identity is reproduced; it is not silently treated as a later-time simulation.

The imported `attempt12/model.py` was also compared byte-for-byte with its file
at the pinned parent commit and matches exactly.

## Full dynamics, not the neutral reduction

`evolve.py` constructs the complete fixed-Q graph Hilbert space, retaining every
Gauss-law matter occupation and electric configuration within the cutoff L.
It starts with the partner's equal-phase neutral state over graphs in which
vertex 0 is isolated. Electric cutoffs regulate the U(1) fields; they are not
modifications to the physical law. There is no graph symmetry quotient.

Measured targets are reference attachment, global connectedness, a spanning
path for Q=N-1, a spanning ring for Q=N, isolated-vertex count, largest component,
maximum degree, charge density and electric excitation. Paths/rings are simple
shape controls; they are **not** evidence for three-dimensional spacetime.

At five vertices, Q=4, time 20 (L=5):

| Quantity | Neutral-only walk | Full CP12 |
|---|---:|---:|
| Reference attached | .84446707 | .31311136 |
| Globally connected | .53521517 | .23860901 |
| Spanning path | .23969148 | .11703030 |

The full dynamics slows this assembly route substantially. Nevertheless,
the conditional distribution within fixed reference degree k ceases to be
uniform. Replacing it by a uniform distribution **at the actual full-model k
probabilities** gives connected probability .21682710 and path probability
.09633119. The actual values are larger by about .02178 and .02070 respectively.
This diagnostic controls for attachment rather than confusing more incident
edges with geometric organization.

At time 40 the full connected/path probabilities are .34209591/.08588467.
The full distribution's total variation distance from its k-uniform comparator
is .13408029. The path probability is not monotone: this is unitary evolution,
not a relaxation or stable-phase demonstration.

## An explicit coherence comparison in the full theory

`dephased_control.py` evolves all 15 initial neutral graph basis states for
N=5,Q=4, then averages their probabilities. It compares this exact mixture to
the coherent state with the same initial graph probabilities.

| Time | Connected, dephased | Connected, coherent | Path, dephased | Path, coherent |
|---:|---:|---:|---:|---:|
| 10 | .04111161 | .11035697 | .01558914 | .04463197 |
| 20 | .10450128 | .23860901 | .04320794 | .11703030 |
| 40 | .14483507 | .34209592 | .05452973 | .08588465 |

These are L=4 full-law results, checked against L=3 and the coherent single-state
calculation. Dephasing raises the initial system energy from 2.16 to 2.24. This
is an inequivalent physical preparation, not an equal-energy experiment or a
proof that coherence can be prepared for free. The advantage survives at these
sizes; its large-N survival is open.

## Q=N: a stronger geometric target

At N=5,Q=5,time20, L=5, the full system has:

    attachment .22895268
    connectedness .22259310
    spanning ring .0071175813

The neutral comparator gives .74804909, .69397134, .027311388 respectively.
The k-matched uniform distribution gives ring probability .0059359481.
The full theory modestly favors rings relative to that matched distribution,
but fewer than one percent of outcomes are rings at this time. The completely
uniform fixed-Q graph baseline is 12/252 = 4.7619% rings. These different
comparators answer different questions; none justifies a large-N selection
claim.

## A derived interaction, and why it is insufficient

The full derivation is in [DERIVATION.md](DERIVATION.md). For a fixed graph with
neutral filled-edge vacuum and hopping operator V stripped of its coefficient,

    <V^4> = 12Q^2 - 8Q - 4 sum_i binom(d_i,2).

Compatible two-hop processes are reduced by Pauli exclusion at shared vertices.
For the declared comparator omitting graph rewiring, perturbation about zero
hopping gives, with Delta=2m+epsilon_E/2,

    E_G = kappa Q - 2Q t_hop^2/Delta
        + [4Q+2 sum_i binom(d_i,2)] t_hop^4/Delta^3 + O_G(t_hop^6).

The degree penalty is derived from CP12's own matter/gauge physics. At fixed Q
it favors more even degrees. The coefficient is checked both by matrix moments
and the perturbative resolvent, then by small-hopping exact diagonalization.
The expansion is **not controlled at the frozen t_hop=1 by this calculation**;
all frozen-coupling comparisons use direct eigensolvers instead.

Regularity still allows fragmentation. At N=6,Q=6, retaining the actual rewiring
diagonal in the graph compression, cutoff-converged block energies are:

| Graph | E_block |
|---|---:|
| One spanning ring | -.635203537847 |
| Two disconnected triangles | -.637262028038 |

No single-edge swap connects distinct 2-regular graphs. Consequently the
triangle trial state lies below every state confined entirely to spanning
rings in this regulated six-vertex system. This is not a full ground-state
distribution calculation: admixture of other geometries can lower energies.

Separate fixed-graph matter calculations for rings N=3 through 8 also favor
the triangle's energy per vertex over each larger tested ring. That is a
finite-family diagnostic, not a proven thermodynamic limit.

## Exact loop discrimination beyond degree regularity

The bare fixed-graph energies of one six-cycle and two triangles agree exactly
through hopping order ten. At order twelve their difference is positive:

    E(C6)-2E(C3) = .15138638963030954 t_hop^12 + O(t_hop^14).

Exact rational recursion and an increased electric cutoff agree. A separate
small-hopping eigensolver control is included. This is a conditional
weak-hopping result; the full rewiring dynamics and the frozen t_hop=1 cannot
be inferred from this series. It explains why degree regularity alone is not
the relevant selection test: topology-dependent matter energy can prefer
small disconnected loops. See DERIVATION.md and loop_series_results.json.

## Numerical scope and error control

Every Hamiltonian construction checks Gauss's law and Hermiticity. Every time
point checks state norm and energy conservation. Sparse matrix exponentials
are used, without a manual time-integration step. Independent cutoff scans,
neutral-walk comparisons, analytic initial variances, batch versus single-state
propagation, and small-time Taylor coefficients provide separate checks.

N=4,Q=3 is run through 2pi/Gamma at L=4,5. N=5,Q=4 is run through T=40 at
L=3,4,5; N=5,Q=5 through T=20 at L=3,4,5. The largest five-vertex basis has
839,082 states. Differences between L=4 and 5 in the selected N=5 geometric
probabilities are small, but are empirical convergence evidence, not rigorous
infinite-cutoff error bounds. Electric observables converge more slowly.

`verification.json` records the actual cross-cutoff maxima and all checks.
The six-vertex partial execution and cutoff limitations are recorded separately in
`LARGER_SYSTEM.md`; it is not assumed converged from five-vertex behavior.

An initial six-vertex run was interrupted to avoid excessive peak memory.
The implementation was improved to release the basis dictionary/state tuples
after constructing the sparse Hamiltonian and readout arrays. It now also
writes a checkpoint after each completed time. This memory change does not
alter matrix elements or evolution. Interrupted attempts are not counted as
successful late-time calculations.
The first memory-reduced N=6,L=3 run still ended with exit137 during the
10-to-20 advance while another large run was active. Its valid time-10 checkpoint
was preserved. The final implementation composes intervals of duration at
most one under the identical Hamiltonian; a five-vertex comparison checks this
execution change against the original direct interval calculation.

## Research decision

The neutral-only no-assembly theorem remains correct in its stated domain.
Its uniform conditional graph premise fails in full CP12, and finite-system
coherence enhancement is real. What survives is a concrete microscopic
selection question: can degree regularization and gauge/matter dynamics beat
fragmentation and graph entropy in the same large-system state?

This packet does not answer yes. The hostile energy comparison explains why
another favorable attachment curve is insufficient. A successor argument must
distinguish one extended geometry from many locally similar fragments, track
the same state's operational propagation, and pay for preparation/coherence.
No additional geometric penalty has been inserted to obtain a desired answer.

## Reproduce

From this directory, with NumPy/SciPy installed:

```sh
OPENBLAS_NUM_THREADS=1 python evolve.py --n 5 --q 4 --L 5 --end 40
OPENBLAS_NUM_THREADS=1 python evolve.py --n 5 --q 5 --L 5 --end 20
OPENBLAS_NUM_THREADS=1 python dephased_control.py
OPENBLAS_NUM_THREADS=1 python matter_selection.py
OPENBLAS_NUM_THREADS=1 python cycle_competitors.py --L 4 --maxn 8
OPENBLAS_NUM_THREADS=1 python compressed_competitors.py
OPENBLAS_NUM_THREADS=1 python short_time.py
OPENBLAS_NUM_THREADS=1 python loop_series.py
OPENBLAS_NUM_THREADS=1 python loop_numeric_control.py
python verify.py
```

`verify.py` validates the recorded controls; it does not rerun the expensive
evolutions. Actual numerical outputs and all author scripts are included.
