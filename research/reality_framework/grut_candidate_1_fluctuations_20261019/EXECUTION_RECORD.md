# Execution and provenance

Parent state: `7c06320abb66c8ba44f02a397abac4062eabc6be` in `ryangrvr/TestingGRUT`. Work continues the owner's request to concentrate on the rotating finite-charge fluctuation problem at unchanged action and matching. No original candidate, classical follow-up or preceding scalar-loop packet was edited.

## Actually executed

| Calculation | Source/protocol | Evidence |
|---|---|---|
| Rotating scalar spectra and stationary backgrounds | `modes.py`, `mode_protocol.json` | 25 initial spectra, all completed, about16.56s |
| Matched-spacing volume and tightened phase controls | `modes.py --box-check`, `mode_protocol_amendment_01.json` | 11 additional spectra, about6.49s; initial phase split retained |
| Saved-profile angular comparison | `angular_bound.py`, `angular_bound_protocol.json` | Three grid refinements and spline maximum; whole-half-line bound explicitly uncertified |
| Additional ell5 controls | `modes.py --ell5-check`, `mode_protocol_amendment_02.json` | Four additional spectra, about3.13s; total40 |
| Independent charged vacuum pole threshold | `pole.py`, `pole_protocol.json` | Ten acceptance checks; analytic bubble identity, quadratures, derivative/scale/curvature comparisons |
| All scalar vacuum poles by full vertex tensor | `neutral_poles.py`, `neutral_pole_protocol.json` | Six checks; reproduces X and adds phi,D |
| Integrating linear-algebra, numerical and boundary checks | `verify.py` | 42/42 checks in `VERIFICATION.json` |
| Fresh-directory replay | `reproduce.py` | Eight commands succeed;90 substantive comparisons pass in `REPRODUCTION.json` |
| Scientific plot | `plot_results.py` | `fluctuation_diagnostics.png`, generated from recorded data and visually inspected |

The scalar modes use NumPy/SciPy and one BLAS thread. Source outputs contain the exact runtime versions. All workloads finished; no ongoing job, package installation, paid compute or empirical dataset was used. The eight-command replay preserves originals and uses the same runtime with numerical tolerances; it is not an external replication or bitwise identity of sparse eigenvectors.

## Preserved negative/unfavorable diagnostics and corrections

- The initial N180 background tolerance produced an apparent real phase-zero split, approximately±3.40e−7. Original roots, phase overlaps and files were preserved. A pre-run amendment tightened that stationary solve; the roots changed to an imaginary pair near±1.65e−8 with phase overlap essentially one. The report describes this as numerical zero/Jordan conditioning, not as evidence of a physical instability that was removed.
- The unconstrained radial amplitude operator's negative eigenvalue, about−.188834, is retained. The actual generator has the original gyroscopic term; the positive fixed-charge rank-one term appears only in the separate constrained-energy diagnostic. No ad hoc dynamical stabilization was added.
- Discrete stationary preparations are at fixed original Omega, with charges converging to the original BVP value. They are not falsely labeled as exactly equal-Q backgrounds on every grid.
- The first X-pole run reached a NumPy-boolean JSON serialization error before saving its results. The explicit correction record changes serialization and an equivalent expression only. No physical input, counterterm, formula or acceptance threshold changed.
- Internal review corrected two report phrases before posting: the calculation covers bounded parts of the prior specification, and the conditional stability implication is proved rather than its continuum hypotheses being established. No scientific rerun was required by those wording changes.

No spatially varying quantum energy or raw frequency sum was computed. The absence of a subgap mode in bounded probes, or a small elementary-particle pole correction, is not relabeled as quantum binding or lifetime.

## Input integrity

The original profile SHA256 is `0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3`. Other imported inputs are the preceding packet's `flow_check.py`, `effective_potential.py` and `potential_results.json`. Their hashes are recorded in the pole outputs and packet manifest. The replay copies these parent inputs into their original sibling paths.

The new packet is posted on the separately named branch `codex/grut-candidate-1-fluctuations-20261019`. The only out-of-packet change is the coordination claim/log. The frozen default/historical branches, benchmark coefficients and existing reports are preserved.

## Review and unresolved work

Constructor, calculator, challenger and integrator worked on complementary parts and inspected actual files. This is internal adversarial and numerical cross-implementation review. It is not independent validation, an experiment, a complete continuum spectrum or a quantum theory-of-reality demonstration.

The remaining decisive calculation is the renormalized rotating fixed-charge quantum energy with continuum normalization, collective modes and ultraviolet counterterms, followed by like-order particle-emission/fragmentation comparisons. Physical total-vacuum matching, apparatus implementation, empirical parameter constraints and Universal GRUT remain separate obligations.
