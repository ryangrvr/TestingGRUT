# Causal geometry and quantum response: notes integration

Date: 2026-10-11. Additive research analysis, not a new generating-law attempt. Frozen CP packets unchanged. No claim of an established theory.

## Executed controls

Both supplied Python scripts were executed unchanged and passed.

`grut_causal_geometry_reconstruction.py`: ordering fractions for dimensions 2, 3 and 4 were 0.488653, 0.228241 and 0.094128; the discrete estimator selected the correct supplied dimensions. Poisson proper-time estimator mean was 0.99999508, standard deviation 0.00250483. The nonmanifold two-layer-plus-isolates control had fraction 0.100191 and chain height 2; a separate four-dimensional sprinkling had height 12.

The dimension-4 fraction and height 12 belong to different generated samples. This does not invalidate the counterexample but should not be presented as a paired measurement on one sample. Pair comparability indicators are correlated; their uncertainty is not the binomial error for independent pairs. The density and continuum interval used to calibrate the clock are supplied. The script demonstrates reconstruction within a known model, not dynamical emergence or a general manifold certificate.

`grut_reality_bridge_recon_test.py`: the rotated outgoing momenta produce different transverse-traceless soft factors despite equal total momentum; Ward residual was 2.51e-16. The supplied Newtonian branch interaction gave phase 0.133632359172 and concurrence 0.0667664749696. The separate Gaussian radiative benchmark gave decoherence exponent 8.71e-65. Their ratio is not a universal ratio for one complete apparatus protocol.

This usefully separates conservative correlations from outgoing radiative records. Endpoint conservation and a soft-factor check do not replace a conserved stress-energy history including the preparation apparatus.

## Literature check

Hiroki Matsui, *Graviton-induced which-path decoherence in matter-wave interferometry*, arXiv:2607.20867, July 2026, and *Gravitational decoherence from soft-graviton S-matrix*, arXiv:2609.08584, September 2026, are real preprints. The September bibliographic record and abstract were retrieved; its full derivation was not independently reproduced here. These calculations assume linearized quantum gravity, rather than generate it.

The 2024 causal-set article is by **S. Carlip**, not Peter Carlip: *Causal sets and an emerging continuum*, General Relativity and Gravitation 56, 95, DOI:10.1007/s10714-024-03281-1. Relevant sections were inspected. It discusses suppression of dominant nonmanifold families, while explicitly leaving the remaining unsuppressed sets insufficiently understood. Its BDG action uses dimension-dependent coefficients and an Einstein-Hilbert continuum target. This is important progress, but does not establish dimension selection from unstructured quantum histories.

Causal-order reconstruction theorems concern suitably distinguishing Lorentzian spacetimes. They constrain how much metric information is redundant once that class is assumed. They do not prove arbitrary finite orders possess such a continuum interpretation. Volume fixes the remaining conformal scale only with appropriate measure information.

## Additional exact interference control

Fix two labelled layers of sizes n1 and n2. Allow each possible cross-layer relation independently; no within-layer relations exist. Every such graph is a valid poset. There are L=n1*n2 possible links. For the supplied link phase theta, summing all configurations gives exactly

    Z(theta) = sum_m binomial(L,m) exp(i*theta*m)
             = (1+exp(i*theta))**L.

Consequently

    |Z| = [2*|cos(theta/2)|]**L.

This illustrates an important entropy-versus-interference distinction. Relative to the raw configuration count 2**L, the magnitude is suppressed whenever |cos(theta/2)|<1. But its absolute magnitude still grows exponentially when 2*|cos(theta/2)|>1. Absolute suppression requires |cos(theta/2)|<1/2. At equality it is order one. No claim about the entire causal-set path sum follows: other histories, measures, boundary states and observables must be included.

The formula can be verified directly by the binomial theorem; no numerical precision campaign is needed. It is a conventional link-action control, not a proposed GRUT law. A small normalized contribution alone is not proof that manifoldlike histories dominate.

## Research implication

The supplied work makes causal quantum dynamics a credible foundation question, complementary to the response/noise and gravitational-dressing work in DEEP_RESEARCH_2026_10_11.md. A meaningful mechanism must do more than fit a dimension statistic: it must yield a controlled continuum regime, physical probabilities, propagation and a gravitational response under its own stated dynamics.

The strongest synthesis is to treat causal structure, physical response and records together: which histories support stable causal propagation; what correlations their dynamics transfers; and which relational probes can recover those correlations. Using measured quantum rules or established low-energy gravity as correspondence requirements is legitimate. Calling supplied geometry or amplitudes derived would not be.

A falsifiable first target for a future authorized model is a response or fluctuation relation that remains fixed after independent low-energy calibration and differs from complete gravitational EFT. Existing scripts do not provide that relation. They provide controls against false emergence and false identification of information transfer with radiation.

No CP-07, Stage-3 card, canonical change, remote push or publication claim is made.

## Primary sources

- https://arxiv.org/abs/2607.20867
- https://arxiv.org/abs/2609.08584
- https://doi.org/10.1007/s10714-024-03281-1
- https://doi.org/10.1007/s41114-019-0023-1
