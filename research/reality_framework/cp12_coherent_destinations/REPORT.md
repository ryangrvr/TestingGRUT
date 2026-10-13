# CP-12: coherent formation, destination locality, and preparation entropy

13 October 2026. Additive follow-up on `codex/grut-cp12-coherent-destination-response`, based on the destination audit at `8cd95db`. Reviewed inputs include the gravity-response packet `11fbaec`, Claude's adjudication `694463f`, and the new unregulated Gibbs audit `4bb0fa5`. Frozen CP-12, canonical GRUT-RAI, prior constructions, and Stage-3 Card 3 are unchanged.

## What moved forward

The full frozen Hamiltonian can generate nonuniform graph correlations from an unordered pure preparation. This is stronger evidence of a mechanism than a favorable energy assigned to a supplied lattice. It is still finite-system organization, not a demonstrated spatial phase.

The calculation also identifies its resource: a globally pure coherent initial state. A separate unregulated proof shows that isolated unitary dynamics cannot make a sufficiently narrow spatial target typical from a high-entropy neutral graph mixture, at bounded energy density in sectors with `Q/N>1`. Nonequilibrium formation is therefore not an automatic escape from the Gibbs obstruction.

## Formation without a supplied spatial ordering

Prepare equal neutral amplitudes over **every** labelled graph with `Q=N-1`. No spatial pool, dimensional embedding, loop penalty, or graph classifier is inserted in the Hamiltonian or these amplitudes. The initial energy is `kappa Q`; the state is permutation symmetric.

We evolved its complete signed permutation sector, retaining matter, Gauss constraints, rotors, and off-diagonal graph mixing. Global connectedness and path probability increase:

| Events | Connected initially | Connected at T=16 | Spanning path initially | Spanning path at T=16 |
| --- | --- | --- | --- | --- |
| 5 | 59.5238% | 74.2115% | 28.5714% | 43.9348% |
| 6 | 43.1568% | 58.3277% | 11.9880% | 23.1970% |
| 7 | 30.9727% | 47.1043% | 4.6440% | 12.0642% |

The absolute probabilities decrease across these sizes. Their enhancement relative to the initial ensemble is not a proof that a spatial class becomes typical. A spanning path is also not three-dimensional space.

For five events, every labelled edge maintains mean occupation 0.4 while these shape probabilities change. Uniform labelled responses do not imply a stationary relational geometry. Permutation symmetry constrains one- and two-point label readouts; it does not eliminate higher correlations.

An exact rational expansion at four through seven events identifies the first graph-probability bias, at order `T^6`, against concentrated vertex degree. It vanishes at orders two and four. This is full-law short-time dynamics, not an extrapolation of the weak-hopping ground-energy expansion. Degree organization still cannot distinguish all spatial and nonspatial graphs.

The seven-event cutoff-five Hilbert space has 487,273,213 physical states. Its chosen invariant sector has 139,097 states; only that sector was evolved. The compression is exact for the prepared state. Complete intertwiners against the frozen full builder at four and five events validate the Hamiltonian and evolution. No claim is made to a full-spectrum solution at seven events.

## Coherence does not fix the first destination rule

The definite-graph first-rewiring theorem extends to arbitrary coherent graph/matter states whenever the candidate edges are absent in **every** initial branch. For the finite-support preparations tested here, the common leading coefficient is

\[
P_f(T)=\frac{\Gamma^2}{M^2}\operatorname{Tr}(D\rho D^\dagger)T^2+O_N(T^3),
\]

where `D` coherently deletes a removable edge. Coherence and records change this deletion norm but cannot favor one common-vacant destination over another at that order. Coherent spanning-path preparations inside a supplied ring retain the same vanishing fixed-radius local-destination fraction as the original ring audit.

Full matter dynamics does break the equality later. In one five-event, four-connection control with a supplied square-tail reference, two equal-degree candidate destinations have reference distances two and three. Exact unregulated arithmetic gives

\[
P_{\rm near}(T)-P_{\rm far}(T)
=-\frac{1103}{3240000000}T^{10}+O(T^{12}).
\]

The longer destination is initially favored; the sign reverses at later tested times. Thus microscopic destination blindness does not prove all-order response blindness, but the actual generated bias here does not produce a monotone locality law. These probabilities are not the earlier cubic source-response experiment, and the reference distance is supplied rather than operationally derived.

## The Gibbs result and the nonequilibrium qualification

The unregulated Gibbs argument at `4bb0fa5` is consistent with the frozen law. `GIBBS_AND_PREPARATION.md` gives a shorter proof by directly bounding the graph probability in the positive electric/mass reference Gibbs state, then applying quantum data processing. It yields `O(1/log N)` rarity without introducing an intermediate cutoff.

The same file strengthens the electric-energy control and extends sufficient suppression to a slow logarithmic cooling window. These are bounds, not phase-transition temperatures.

It also proves an all-time result for the full neutral graph mixture under isolated unitary evolution. Conservation of energy, invariance of the density-matrix spectrum, and positive electric energy give

\[
P_T(\mathcal S_N)\le O(1/\log N)
+\exp[-(\nu-1)N\log N+O(N\log\log N)]
\]

for `Q=nu N`, `nu>1`, and the explicitly narrow graph-family count. More generally, bounded energy density and an initial entropy within `O(N)` of the neutral mixture cannot generate a nonvanishing probability for that target through closed unitary evolution.

Our pure unordered preparation has an entropy deficit of `log binom(M,Q)` relative to that mixture. It supplies this resource without supplying a spatial ordering, but no mechanism preparing its purity and phases has been derived. The formation calculations have `Q=N-1`; they must not be combined with the `nu>1` theorem as one realized phase.

A conditional extension rebuilds the Gibbs argument in the specified signed permutation sector. It still excludes exponentially many target templates in `nu>1`; simply calling event labels gauge does not rescue that particular ensemble. Other gauge prescriptions require their own definitions and counts.

For a genuinely broad operational quantum geometry, the appropriate entropy diagnostic is the reference probability of its physical projector. Broad support over graph labels alone is insufficient if the physical code subspace remains small. A sufficiently large fluctuating phase or a justified low-entropy/open-system preparation remains possible.

## Numerical and reproduction limits

Five-event full evolutions use electric cutoffs two through six. Cutoffs five and six agree on the reported destination occupations to below `1e-15`, and on unordered shape statistics to about `1e-13`. Six-event invariant evolutions use cutoffs two through six; seven-event calculations use two through five. Seven-event connected/path probabilities differ by less than `2e-10` between cutoffs four and five at the reported times. Electric-energy observables converge less tightly. These comparisons are numerical evidence, not interval-certified removal of the regulator or a thermodynamic limit.

Exact Taylor calculations use an independent integer reachable-state builder without an electric regulator. The unregulated entropy claims are analytical estimates; small matrices validate their algebraic steps, not the asymptotic theorem. Main source files are hashed in `MANIFEST.json`.

One correction to the adjudication: SciPy 1.17.0 is already listed in the pinned gravity-response manifest's `final_validation_versions`. Installation instructions were missing. This packet supplies `requirements.txt`, explicit commands, and executable imports that also work under `python -I`.

## What this changes about the next test

The next physical test must distinguish a broad operational spatial phase from a restricted list of geometric templates, and account for its preparation entropy. It must then obtain relational distance, tensor propagation, and universal coupling in that same state and sector. The current finite formation result does not supply a physical clock, a continuum momentum, a graviton, or a Lorentzian gravitational constraint.

The earned advance is a full-law correlation-generating mechanism from an unordered preparation, together with a sharper accounting of when it can compete with entropy. It supplies a concrete foundation for further testing; it does not establish a working theory of reality.

## Reproduce

Use Python 3.12 and install the pinned requirements. From the repository root:

```sh
python -m pip install -r research/reality_framework/cp12_coherent_destinations/requirements.txt
python -I research/reality_framework/cp12_coherent_destinations/verify.py --execute-small
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_coherent_destinations/finite_evolution.py --cutoff 6
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_coherent_destinations/invariant_evolution.py --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_coherent_destinations/invariant_evolution.py --exact
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_coherent_destinations/invariant_evolution.py --exact-n 7
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I research/reality_framework/cp12_coherent_destinations/invariant_evolution.py --n 7 --cutoff 5
```

The verifier checks records and selected fresh small calculations; it does not silently rerun all large evolutions. The source copy `frozen_model.py` is byte-for-byte the pinned frozen model, included because the adjudication branch does not itself contain Work's frozen packet. It is a reproduction input, not a changed law.
