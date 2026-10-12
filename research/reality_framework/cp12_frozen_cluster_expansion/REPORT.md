# CP12 at frozen hopping: edge-cluster geometry comparison

2026-10-12. Additive research based on 9758c60a7af62f76da9e706fb96c0c548a1cebab.
No frozen law, canonical GRUT-RAI branch, PR12 or Stage-3 card was changed.

## Scientific outcome

We have moved beyond extrapolating a weak-hopping Wilson coefficient. The
finite-cluster ground energies here use the frozen hopping t=1 and retain
fermions, dynamical integer electric fields and exact Gauss constraints.
The true rewiring diagonal is retained with the common thermodynamic host
coefficient a=Gamma=.1. Off-diagonal changes between graph configurations are
not included in this declared graph-block comparison.

The ninth-edge width-sensitive contribution is positive, but the tenth-edge
correction reverses it. The apparent preferred geometry changes with cluster
truncation even after the electric cutoff has converged numerically. This is
the central result of the frozen-coupling test. An isolated flux channel, a
single linked weight and the full thermodynamic energy are different objects;
none can be substituted for another.

That result does not establish a spatial phase. The square-versus-locally-tree-like
energy difference repeatedly changes sign as cluster order increases, and the
tube comparison requires control of the later corrections. Numerical eigensolver precision
and electric-cutoff convergence do not provide a thermodynamic remainder bound.

## Incoming review was executed verbatim

The pinned check.py from the winding-dimensionality review was fetched and
executed without source edits. It passed; partner_verbatim_results.json records
its actual output and source SHA256. The preliminary copy had an extra final
blank line; it was replaced with the exact fetched source before the recorded
verbatim run. The content correction did not change its calculation.

## The nonperturbative finite-cluster calculation

For each connected edge cluster c, define its weight by subtracting its
connected proper edge subclusters:

    w(c) = E(c) - sum_s n(s,c) w(s).

Every E(c) is obtained by numerical diagonalization of the physical CP12 graph
block at t=1. It contains all hopping orders for that cluster. Edge-cluster
order is not hopping order. The inclusion-exclusion identity is exact for a
finite host if all its edge subsets are included. Truncating at finite edge
order is an approximation, with no proven tail bound in this packet.

DERIVATION.md specifies the Hilbert space, inherited compression term,
embedding convention and the exact combinatorial identities. Numerical
linked-cluster expansions are an established method; see Tang, Khatami and
Rigol, Computer Physics Communications 184, 557 (2013), arXiv:1207.3366.
This is a CP12 implementation and test, not a new general method.

All connected bipartite clusters of maximum degree four through eight edges
were generated: 136 topologies. Their union of nonzero square/tree embeddings
contains 124. The L=3 run computes all 136; the L=4 run computes the 124 that
contribute. Both give the same targeted sums. None of the omitted zero-embedding
clusters can be a subcluster of a contributing cluster.

The tree densities are the exact finite-cluster embedding densities in the
infinite four-regular tree, representing the locally tree-like limit of random
four-regular graphs. These are not finite random-graph exact ground energies.
The balanced, fixed-aspect-ratio and width-nine tori have identical embedding
densities through this order. They need longer clusters to be distinguished.

## Square versus locally tree-like competitor

The values below are the difference of truncated relative energy densities;
common neutral baselines cancel. Positive means the tree comparator is cheaper.
These are L=4 results at a=.1, t=1:

| Maximum cluster edges | Square minus tree |
|---:|---:|
| 4 | +0.000102700545 |
| 5 | -0.000303132076 |
| 6 | +0.000209311497 |
| 7 | +0.000045608382 |
| 8 | -0.000056288597 |

The absolute sum of the eight-edge square terms is .02150727, while their
signed increment is only .00009437. That cancellation is another reason a
small last signed increment cannot be treated as a tail bound. The analogous
tree absolute sum is .02930523. These are diagnostics, not a convergence proof.

At order 8 the relative densities are approximately -1.178929998489 (square)
and -1.178873709891 (tree). Adding the common 1.2 baseline gives the corresponding
truncated compressed-energy densities. These are not variational bounds.
The change in their difference from order 7 to8 is about1.02e-4, larger than
the final5.63e-5 apparent preference. A winner is therefore not established.
The sign changes are a direct convergence warning, not a fitted asymptotic
claim or proof that no convergence is possible at higher order.

## First width-nine-sensitive term: both cycle and path matter

The tube has extra nine-cycle embeddings. It also loses nine-edge path
embeddings that wrap back onto their starting vertex, violating injectivity.
The exhaustive relevant catalogue through nine edges has exactly two unequal
embedding densities:

    Delta m(C9)=+1/9, Delta m(P10)=-1.

Thus the first density difference is

    Delta e9 = w(C9)/9 - w(P10)
             = E(C9)/9 - [E(P10)-E(P9)].

At L=4, the cycle-only contribution is 9.6887204432e-9; subtracting the open-path
weight 8.6472773386e-9 leaves

    Delta e9 = +1.0414431046e-9.

L=3 gives +1.0414449303e-9. The difference is about1.83e-15. Independent random
starting vectors reproduce the critical ground energies within 1e-14.
The positive sign is reproducible for this first edge-truncated contribution.
It is not the full tube-minus-plane energy density.

The intermediate27x108 and balanced54x54 families also tie through this order;
their shortest winding cycles have length27 or54. No claim about their complete
relative energies follows from the width-nine calculation.

## The tenth-edge correction reverses the apparent preference

At edge order 10 the embedding differences occur for four branched trees,
one open path, and a nine-cycle with one pendant edge. `tenth_boundary.py`
enumerates them. `tenth_energy_expression.py` performs independent Boolean
inclusion-exclusion to express their entire increment as a linear combination
of thirteen finite-block ground energies. The coefficient-weighted bond count
vanishes exactly; no hidden neutral-baseline shift appears.

Both cutoffs completed. The largest block has 1,240,618 physical states.
The cumulative tube-minus-plane density differences are:

| Maximum relevant cluster edges | L=3 | L=4 |
|---:|---:|---:|
| 9 | +1.04144493e-9 | +1.04144310e-9 |
| 10 | -5.85431770e-9 | -5.85430687e-9 |

The tenth-edge increment alone is about -6.89575e-9, over six times the
magnitude of the positive ninth-edge term. The cutoffs agree on the cumulative
result within about 1.1e-14. Recorded weighted residual indicators are below
1.3e-12; they do not certify the infinite-cutoff or thermodynamic limit.
The positive ninth-edge sign therefore cannot establish that the extended
plane wins. Nor does the negative tenth-edge sign establish that the tube
wins: no bound on the later terms has been proved.

`tenth_results_L3.json`, `tenth_results_L4.json` and `verification.json` contain
the actual outputs. All thirteen finite energies were evaluated at each
cutoff; the critical path/cycle energies were reused only after graph
isomorphism checks. New blocks use seeded random eigensolver starts.

## Numerical error and remaining physical scope

Ground-state residuals and cutoff-boundary probabilities are saved for every
cluster. L=3-to4 differences in the common eight-edge cluster energies are
below 1.6e-14. The square partial-sum difference is below 1e-11. Those errors
are far smaller than the observed change with cluster order. Electric cutoff
control and cluster truncation control are separate requirements.

A Ritz residual is a numerical consistency diagnostic, not an infinite-field
error certificate. The random-start audit checks the critical ground energies.
The initial symmetric vector missed an excitation sector in one path block;
all stored gap entries are consequently marked as uncertified diagnostics.
No propagation speed or physical massless mode is inferred from them.

This packet does not compute a full quantum Gibbs state, a thermal free energy,
or the effects of off-diagonal quantum graph mixing. It does not compare the
complete fragmented four-regular family. The previously established full-CP12
thermal entropy theorem remains applicable in its own stated domain.

## Research consequence

At the actual hopping, graph-dependent matter energy can be computed without
using a small-t polynomial. The first two width-sensitive linked contributions are
now specific physical calculations, including the companion open paths and
branched clusters that a loop-only argument misses. This is a stronger test than an isolated loop
amplitude or a static spectral-dimension diagnostic.

The unresolved obligation is quantitative: control the omitted linked weights
and incorporate graph mixing before identifying a selected macroscopic
geometry. The current square/tree oscillations show why that obligation is
substantive. Local organization, winding response, extended geometric order,
and operational propagation remain distinct achievements.

## Reproduce

With numpy, scipy and networkx installed, run from this directory:

    OPENBLAS_NUM_THREADS=1 python clusters.py --L 3 --maxedges 8
    OPENBLAS_NUM_THREADS=1 python clusters.py --L 4 --maxedges 8 --embedded-only
    OPENBLAS_NUM_THREADS=1 python winding_weight.py --L 3
    OPENBLAS_NUM_THREADS=1 python winding_weight.py --L 4
    python embedding_boundary.py
    python tenth_boundary.py
    python tenth_energy_expression.py
    OPENBLAS_NUM_THREADS=1 python tenth_energy.py --L 3
    OPENBLAS_NUM_THREADS=1 python tenth_energy.py --L 4
    python verify.py

The scripts require the frozen attempt12/model.py already in this repository.
`verify.py` checks saved numerical records; it does not rerun the expensive
basis construction and eigensolves. Tenth-order runs checkpoint each computed
graph in separate JSON files. They can be resumed without inventing late-time
or missing-cutoff results.
