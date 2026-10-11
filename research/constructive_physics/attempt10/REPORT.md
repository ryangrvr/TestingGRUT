# CP10: joint causal geometry, matter and dynamical reference

2026-10-11. One new constructive candidate, frozen in local commit f0a0832 before execution. Owner instruction work it out. Separate from Stage3; no Card3, canonical integration or earlier frozen-law changes. Verdict: functioning finite positive joint dynamics, no physical continuum established.

## Explicit mechanism

CP09 plus one N-position reference, reference hopping via the positive Laplacian of the current undirected comparability graph, and interaction g sum_i(I-Z_i)|i><i|_R. k=g=1/4. The same Hamiltonian evolves geometry, matter and reference. No instantaneous geometry-role oracle reads a prepared full identity record. Matter interaction occurs at reference occupancy; its mathematical form is explicitly supplied.

All three contributions are positive. At fixed N, total energy is conserved. A total-order system trial with all matter X+ and uniform reference has reference hopping energy zero, interaction expectation g, CP09 energy J(N-1). Thus 0<=E0_10(N)<=J(N-1)+g. This is a ground-energy bound, not a proof of excitation scaling or a thermodynamic continuum.

Simultaneous event permutations act on reference as well as matter and geometry. Hamiltonian covariance was verified exactly for an adjacent label swap at N2,3. The initial CP09 ground times uniform reference is invariant; the full evolution therefore remains invariant. This actual reference has nontrivial motion and interaction energy rather than an assumed zero Hamiltonian compensating representation.

## Executed construction results

N2 joint dimension24, groundenergy .260465360127; initial totalenergy .344755283698. N3 dimension456, groundenergy .482868701245; initial totalenergy .560979226644. Initial geometry/matter groundstate and reference are a product; no precomputed geometry-reference correlations supplied.

At t1:

| N | reference purity | geometry/reference classical mutual information (nats) | mean relational orientation |
|---|---:|---:|---:|
|2|.94468135|.00036222169|.02217897|
|3|.93343162|.00074708852|.01912571|

Orientation is the same geometrical D_i, evaluated at the reference's occupied site, an invariant joint observable. Mutual information uses the full joint geometry/reference endpoint distribution, not selected records. It is small but nonzero; purity alone would show correlations with the combined system, not specifically geometry, which is why both diagnostics are retained. These are endpoint measurement statistics, not a stabilized memory or definite outcome mechanism.

Norm, inverse evolution and energy checks pass within 1e-10; energy differences observed near machine precision. Reversing the complete state restores the initial state, including reference. No fundamental irreversibility.

## Supplied structures and failure boundaries

The comparability graph is not derived spatial adjacency. It connects events separated by long causal paths as well as cover-related events, and erases edge orientation for reference hopping. This is an explicitly postulated motion rule, not Lorentz locality, a physical causal cone or an explanation of distance. It respects event relabelling but is not spacetime diffeomorphism invariance. H09's geometry transitions remain globally prescribed, with external time.

The labelled reference basis and finite tensor factors are supplied. Operationally reading the invariant orientation still needs a measurement interaction or instrument; the current calculation predicts its mathematical distribution and mean, not a local laboratory implementation. Reference correlations do not show that a physical observer knows its position in an emergent metric. No local light cone, four-dimensional geometry selection, gravitational constraints or spin2 sector has been derived.

The new law therefore earns joint reference dynamics and correlation generation. It does not yet earn propagation locality. A subsequent analysis must test its response against independently specified reference separation and determine whether those quantities are stable across fluctuating configurations, instead of defining distance to fit the response. No parameter scan, new penalty or preferred dimension has been inserted after tests.

## Reproduce

python research/constructive_physics/attempt10/check.py, NumPy/SciPy. Uses the published independent CP09 sparse builder. All N2,3 labelled orders included. Terminal CP10_JOINT_DYNAMICAL_REFERENCE_PASS. Numerical probes evaluate total energy, covariance, inverse evolution, reference reduced purity, geometry-reference joint distribution and relational orientation. Standard supplied quantum dynamics; not a novel empirically established law. Independent review outstanding.
