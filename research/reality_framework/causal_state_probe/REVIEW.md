# Causal response and state probe: reproduced results and adversarial boundary

2026-10-11. Additive analysis of supplied notes; no new action or candidate law. CP-05 and other frozen packets remain unchanged.

## Reproduction

The uploaded script is preserved verbatim as `supplied_probe.py`; its notes are preserved as `SUPPLIED_NOTES.txt`. All supplied assertions pass. Full outputs are in `supplied_results.txt`.

| 180-event order | Relation count | Ordering fraction | Height | Commutator rank | SJ positive modes |
|---|---:|---:|---:|---:|---:|
| Supplied 4D sprinkling | 1413 | 0.0877094972 | 6 | 156 | 78 |
| Complete 37-by-38 layers plus isolates | 1406 | 0.0872749845 | 2 | 2 | 1 |
| New random 90-by-90 layered control | 1413 | 0.0877094972 | 2 | 180 | 90 |

The last row is a new conventional counter-control, executed by `adversarial_rank.py`. Results are preserved in `adversarial_results.json`. It uses the identical unnormalized link prescription, no within-layer relations, fixed seed 71011 and exactly the sprinkled sample's relation count.

## Exact explanation and stronger boundary

For any two-layer order with cross-layer link block B and response normalization gamma,

    Delta = gamma * [[0, -B.T], [B, 0]].

Each nonzero singular value sigma_j(B) produces eigenvalues +gamma*sigma_j and -gamma*sigma_j of i*Delta. Thus

    rank(Delta) = 2*rank(B),
    number of positive SJ modes = rank(B).

Complete bipartite coupling has rank(B)=1. That explains the supplied rank-two result. Arbitrary layered coupling need not have low rank. In the added binary control B is full rank 90; its minimum singular value is 0.0024505, comfortably above the numerical tolerance.

Consequently a large commutator rank cannot be promoted into a criterion selecting continuum histories. Height already distinguishes these specific fake orders from the sample. The quantum response remains useful for testing propagation on a specified order, but this example does not establish a stronger general manifold classifier than order statistics.

These height-two orders are non-continuum controls in the intended large-system sense; no claim is made that an individual finite order can never embed in any specially chosen spacetime region.

## State construction: correct, conditional, and not a history amplitude

For a finite real antisymmetric Delta, H=i*Delta is Hermitian and W=H_+ is positive, with W-W.T=H and N=|H|/2. The supplied alternative state (1+lambda)*W+lambda*W.T has the same commutator, preserves its kernel and multiplies symmetric fluctuations by 1+2*lambda. These identities reproduce numerically.

A field algebra still needs its zero-mode relations and dynamics specified; a positive matrix alone does not construct an interacting matter theory. The SJ prescription depends on the chosen region, propagator and inner product/volume measure. Relabelling covariance does not imply independence from changing the region or physical measure.

The 4D script uses a bare link matrix, omitting physical density-dependent normalization. This is harmless for its rank comparison but prevents interpreting its noise magnitudes as calibrated detector fluctuations. Its introductory docstring describes the 2D toy; the actual script also includes the 4D benchmark, as reported in the notes.

Most importantly, W is a state covariance, not a derived relative amplitude for competing causal histories. It must not be inserted as a history-selection weight without another physical assumption. A finite retarded matrix is strictly triangular after topological sorting and therefore singular; writing a Gaussian history weight as det(G_R)^(-1/2) would not supply a defined action or measure. A legitimate determinant construction needs an actual kinetic operator, boundary/zero-mode prescription and integration measure. It also needs a common normalization across histories.

## Literature and physical next step

Johnston's 2008 propagator construction checks correspondence on sprinkled 1+1 and 3+1 Minkowski spacetimes; it does not derive a microscopic wave equation on every arbitrary order. Johnston's 2009 work establishes the causal-set free-field construction. Fewster and Verch's 2013 work shows generic SJ states on ultrastatic slabs can yield infinite fluctuations for derivative Wick squares: finite matrix positivity cannot substitute for continuum state regularity.

The notes' coupled matter/causal-history direction is scientifically meaningful. The outstanding calculation must specify a physical history dynamics and matter measure, obtain their joint amplitude rather than fit a geometric preference, and test it against the abundance of non-continuum competitors. The present state prescription cannot provide that amplitude by itself.

For detector records, the operational next calculation would use physically specified source smearing and response/noise kernels in an influence functional. The order-defined probes in the supplied code are intentionally chosen after inspecting future profiles; they are diagnostic witnesses, not independently calibrated detectors. No detector experiment, spin-2 sector, universal coupling, quantum completion or new GRUT law is established here.

## Reproduce

From this directory, with NumPy installed:

    python supplied_probe.py
    python adversarial_rank.py

Primary sources:

- Steven Johnston, Particle propagators on discrete spacetime (2008): https://arxiv.org/abs/0806.3083
- Steven Johnston, Feynman propagator for a free scalar field on a causal set (2009): https://arxiv.org/abs/0909.0944
- Christopher J. Fewster and Rainer Verch, The necessity of the Hadamard condition (2013): https://doi.org/10.1088/0264-9381/30/23/235027
- Sumati Surya, The causal set approach to quantum gravity (2019), section 5.2: https://doi.org/10.1007/s41114-019-0023-1
