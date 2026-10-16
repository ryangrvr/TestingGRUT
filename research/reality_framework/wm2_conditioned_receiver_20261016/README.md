# WM2-S6: history-conditioned receiver calculation

Read REPORT.md, THEORY.md, COST_AND_SCOPE.md and NUMERICAL_LIMITS.md.
The two numerical protocols have separate public pre-trajectory charters:
CHARTER.md and GLOBAL_RESIDUAL_CHARTER.md. The physical law and preparation
are unchanged from S4/S5. No prior research packet is modified.

The exact content is a complete cross-history Fock transition calculus
and a physical-H state-error bound using the full projected-emission
residual. ACTION_ENRICHMENT.md additionally derives exact coherent
action densities and a conditional Krylov-enclosure error-scaling bound.
Its algebra controls are separate; it was not used to change either
primary trajectory protocol. This must not be conflated with efficient ranks,
all-size formation, an isolated real-energy sector or a spatial phase.

## Saved-evidence verification

From the repository root, with Python, NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/wm2_conditioned_receiver_20261016/reproduce.py
```

This checks immutable charter hashes, every packet file, both complete
attempt lists, original root weights and sign reuse, numerical gates,
unknown target-time populations, all reference comparisons and the
unchanged S5 packet. It does not rerun the long trajectories.

## Independent execution

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/wm2_conditioned_receiver_20261016/reproduce.py --recompute --output /tmp/grut_s6_rerun.json
```

This makes a temporary copy of S6 and its S5/S4 dependencies, executes
the original local-quota campaign, additional controls, then the separately
frozen cumulative-budget campaign. The archive is unchanged. The run can
take tens of minutes. It exports new evidence, not a byte-identical
wall-time transcript. Different hardware or exact arithmetic reuse can
change the endpoint of a computational stop without changing the test.

## Evidence and code

* RESULTS.json: original S6 protocol, controls and every primary stop.
* GLOBAL_RESULTS.json: G1 protocol and its SHA-256 of the unchanged
  original S6 evidence. A G1 outcome is not an original-quota outcome.
* ADDITIONAL_CONTROLS.json: larger matrix-free maps, bit-identical
  rotation caching, the common-frame S5 limit and exact storage counts.
* ACTION_CONTROLS.json and action_enrichment.py: coherent action-density
  and exact Krylov-support identities, tested separately after the
  primary runs; no new formation result is claimed from these controls.
* CONTROLS.json: initial algebra/wavefunction checks; the campaign
  repeats these checks and records them with its physical references.
* fock.py: normalized rectangular Fock lifts, dense and matrix free.
* conditioned.py: exact histories, per-history frames, full-action
  residuals, projected exponentials and all reduction costs.
* campaign.py and global_residual.py: the two predeclared physical
  test lists and distinct numerical allocation policies.
* verify.py, reproduce.py, manifest.json and VERIFICATION.json: read-only
  evidence audit and archive-preserving rerun entry point.

State-error accounts are evaluated in ordinary double precision.
NUMERICAL_LIMITS.md records cancellation, Gram defects, exact-arithmetic
scope, the caching implementation change and remaining certification
debt. No original charged-matter universality, dimension selection or
distinctive empirical prediction follows from these numerical methods.
