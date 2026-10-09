# Operational discriminator boundaries

**DRAFT — pending external review.** These are specification/audit results,
not proposed GRUT candidate laws, certified Stage-3 locks, or new empirical
findings. No experiment was performed. Sources: [SOURCES.json](SOURCES.json).

## Exact discrimination furnished by this package

1. **Memory-capacity test.** T4 specifies four two-bit preparations and a
   common two-step probe. The exact response matrix I₄ excludes all
   classical retained-state models with N≤3. The comparator resource bound
   and common interface are independent requirements. A fixed four-state
   register passes. Conservative finite-data singular-value bounds are
   given in PROOFS.md; a noisy nonzero determinant alone is insufficient.
   A classical digital implementation is straightforward. It is a standard
   control, not a test of a new law of reality.
2. **Historical reinforcement versus a static latent parameter.** T6 fixes
   α=β=1. Passive records of the urn and beta mixture coincide at every
   length. Force the first record to success and increment the urn's
   actual counter; the predicted next-success frequencies are 2/3 for the
   urn and 1/2 for the static latent comparator whose parameter is unaffected
   by forcing. Standard conditioning gives 2/3 in both. This is an exact
   causal contrast. It requires the forcing operation to alter the physical
   reinforcement state as specified; simply relabelling a readout does not
   implement it. A fixed count-state dynamics duplicates the urn under the
   same intervention, so the contrast excludes only the stated static
   comparator. A software urn is feasible; no universal precedent reservoir
   or corresponding quantum intervention has been independently certified.

## Plausible physical signatures and their comparators

| Signature class | Observable / controlled comparison | Fixed-law comparator | Falsification and scope |
|---|---|---|---|
| Dimensionless-constant temporal drift | Calibrated frequency ratios of multiple clocks with distinct sensitivity vectors, same physical scale, time-dependent environmental nuisance controls | Frozen-parameter Standard Model plus GR and apparatus noise; additionally compare a fixed scalar-field EFT | Significant reproducible drift rejects the frozen-parameter model. A **specified** drift law is refuted by incompatible amplitude/phase/correlation bounds. Drift alone does not reject fixed scalar-field dynamics |
| Spatial/branch-history correlations | Independent spectral/clock ratios conditioned on calibrated backgrounds and relic histories | Scalar backgrounds, environment/systematics, fixed phase-transition/vacuum sectors | Reject only if quantitative joint correlations lie outside the preregistered conventional comparator family. Arbitrary history correlations can be hidden-state effects |
| Novelty/repetition response | Held-out identical operational preparations, randomized repetition/conditioning protocols, stable apparatus, predicted dependence on prior encounters | Quantum instruments with specified stationary or aging environment, learning apparatus and latent noise | A numerical novelty law needs an operational sameness criterion, count reservoir and reset/forcing definition. Bounds inconsistent with that law refute it; no such full law is certified here [M02] |
| Nonergodic parameter distributions | Repeated trajectories or independent patches with specified sampling/selection, occupation and switching statistics | Quenched disorder, mixtures, conserved/topological sectors, slow fixed dynamics | Nonergodicity rejects an ergodic comparator. It does not by itself establish evolving fundamental laws or a population measure |
| Parameter attractor distributions | Independently controlled initial distributions and common dynamics; approach rate and terminal distribution | Fixed potential, mutation/replication or scalar dynamics | Exclude a **specified** basin/convergence/distribution prediction with incompatible trajectories. A single observed parameter cannot identify its unique attractor or distinguish an encoded target |
| Departures from EFT running | Compare couplings at matched energy/resolution and chronology, including thresholds and matching uncertainties | Fixed beta functions plus known thresholds; then explicit new-field EFTs | An incompatible running curve excludes the stated EFT, not fixed-law architecture. Temporal drift must be separated from scale dependence |
| Bounded-history/memory violation | Response matrices under independently justified capacity/locality/reset bounds | Finite-state or finite-dimensional specified process class | Certified rank/causal witness outside that bound rejects it. Unbounded hidden-history models remain; quantum and classical bounds differ |

These rows specify signal classes and falsification obligations. They are
not predictions with freely fitted amplitudes dressed as discoveries. The
architecture alone gives none of their quantitative amplitudes, parameter
values, stochastic spectra or distributions. T2 proves the absence of a
discriminator against its unrestricted fixed lift.

## Why varying constants are insufficient to identify law evolution

A conventional fixed extended action can have a scalar-dependent
electromagnetic kinetic factor. In natural units and metric signature
(−,+,+,+), for example, the usual scalar/Maxwell sector contains

$$S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2}{2}R-\frac12g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi
-V(\phi)-\frac14B_F(\phi)F_{\mu\nu}F^{\mu\nu}+\mathcal L_{\rm matter}
\right].$$

With supplied B_F>0, fixed charge normalization and canonically normalized
electromagnetism on a slowly varying background, the local effective
fine-structure parameter scales as α(φ)=α₀/B_F(φ). V, B_F, matter
couplings, gravitational/quantum structure, initial conditions and validity
regime are supplied. This is a conventional comparator family [M14,M15],
not a fully specified new candidate or an assertion that any arbitrary
drift can satisfy all physical bounds. Such models can impose correlated
equivalence-principle and matter-coupling constraints. A physically useful
selector must beat the applicable constraints of an explicit comparator,
not just a constant-parameter fit.

For a clock ratio R_ij=ν_i/ν_j, the small-change relation is

$$\delta\ln R_{ij}=(K_{\alpha,i}-K_{\alpha,j})\delta\ln\alpha
+(K_{\mu,i}-K_{\mu,j})\delta\ln\mu+\text{other couplings and nuisances}.$$

Sensitivity coefficients and nuisance channels need calibration. One ratio
normally cannot identify two independent varying constants; the sensitivity
matrix must have sufficient rank. μ conventions must be fixed (proton/
electron versus electron/proton), since reciprocal conventions reverse
the logarithmic change. Published drift and short-period analyses already
implement restricted tests of this kind [M23,M24]. The dated 2021 Yb+
benchmark reports α drift 1.0(1.1)×10^-18 per year; this is a quoted
historical benchmark, **not** the strongest 2026 bound or a detected drift.
The 2025 Th-229 sensitivity paper is an opportunity/calibration reference,
not evidence of variation [M25]; its full calibration analysis was not
independently accessed here.

## Physical-baseline limitation

The abstract eight-state and register countermodels obey normalized
probabilities and discrete causal updates. They are conventional stochastic/
computational systems. They are not represented as equilibrium Gibbs models
obeying every KMS/FDT condition, relativistic completions, or apparatus-free
laws generating memory/readout from first principles. Thermodynamic costs,
energy bookkeeping and physical apparatus would have to be included for
those stronger comparisons. No such completion is claimed as tested.

This matters in both directions: an arbitrary history kernel is not
automatically conventional local physics, and a hidden-state objection
cannot excuse a proposed physical law from conservation/locality or
experimental bounds. T1 is exact at its stated architecture/operation level.
