# GRUT interacting detector histories — connected construction and results

7 November 2026 · Draft research packet · Zenodo paused

The pass constructs an independent massive bosonic reference for the same intended Schwinger-model detector histories as the lattice calculation. It fixes the boundary-induced mean, separates detector compression from Hamiltonian truncation, and tests which coherences a reduced description must retain. **Interacting continuum matching remains unresolved at the registered accuracy.** No missed gate has been converted into a pass by changing a tolerance or substituting an extrapolation.

| Frozen calculation | Physical question | Result |
|---|---|---|
| [M0-01](RESULT_M0_01.md) | Do fixed-length lattice massless histories match the independent interval vacuum? | Maximum finest TV 0.010011; 0.001 window missed. |
| [M0-02](RESULT_M0_02.md) | Does an analytic endpoint correction improve that comparison? | Maximum TV 0.001678; substantial improvement, still a miss. |
| [M1-01](RESULT_M1_01.md) | Do the independently constructed massive and lattice histories agree with resolved errors? | Maximum TV 0.002472; reference-cutoff and overflow prerequisites also miss. |
| [D1-01](RESULT_D1_01.md) | Which information can be dropped from the interacting finite parent? | Energy-population comparator differs by 0.0221–0.0375. The frozen 32-level coherent reduction also misses. |

TV is total variation on the joint outcome distribution. Percent versions of these numbers are percentage points of absolute probability distance, not relative errors on a rare event.

[Theory connection](THEORY_CONNECTION.md) derives the common instrument description, an autonomous clock-and-record embedding, and the role of effective-state sufficiency. It distinguishes an analytic finite construction from a physical implementation. The new bosonic code is independent of the lattice representation, but has not received an independent external implementation review. Its algebraic and massless checks were performed by the same author.

The prior spectrum result remains a successful benchmark. The prior Gaussian comparison used the lattice's own mean/covariance and therefore did not establish an independent continuum bridge. Neither the spectrum nor the new near-agreement certifies the full equation-(14) claim, causal realization of the ideal extended detector, chiral matter, dynamical gravity, or new GRUT laws.

## Reproduction from a clean repository checkout

Use Python 3.12; the recorded runs used 3.12.14. From the repository root:

```sh
export GRUT_REPLAY_ROOT="$PWD"
export GRUT_OUTPUT="$(mktemp -d)"
python -m venv "$GRUT_OUTPUT/venv"
. "$GRUT_OUTPUT/venv/bin/activate"
python -m pip install -r "$GRUT_REPLAY_ROOT/research/reality_framework/grut_detector_matching_20261107/requirements.txt"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export GRUT_MATCH="$GRUT_REPLAY_ROOT/research/reality_framework/grut_detector_matching_20261107"
python "$GRUT_MATCH/verify_packet.py" --output "$GRUT_OUTPUT/verification.json"
```

`verify_packet.py` reconstructs the stored analyses, checks source hashes, distributions and finite-parent bounds, and explicitly verifies that the misses remain misses. It exits nonzero if preservation/analysis checks fail. It does not rerun production or certify a physical theory.

For an independent reference check and one small interacting point:

```sh
python "$GRUT_MATCH/bosonized.py" --verify --output "$GRUT_OUTPUT/algebra.json"
python "$GRUT_MATCH/bosonized.py" --mass 0.5 --cutoff 4 --output "$GRUT_OUTPUT/massive_E4.json"
```

To reproduce the newly timed finite-parent reduction (one reconstruction, originally 36 seconds):

```sh
timeout 300 python "$GRUT_MATCH/reduce_d1.py" --output-dir "$GRUT_OUTPUT/d1"
```

To rerun the entire bounded massive comparison:

```sh
python "$GRUT_MATCH/run_m1_01.py" --output-dir "$GRUT_OUTPUT/m1"
```

The execution wrapper preserves completed points, exceptions and resource caps. Its zero exit means the result was recorded, **not that the scientific gates passed**; read the emitted verdict and `analysis.json`. It requires a clean committed tree and refuses existing output directories. The standalone workers also refuse overwriting their output. A timeout exit is a resource miss and must be retained.

The massless comparisons are independently reproducible with:

```sh
python "$GRUT_MATCH/run_m0_01.py" --output-dir "$GRUT_OUTPUT/m0_01"
python "$GRUT_MATCH/run_m0_02.py" --output-dir "$GRUT_OUTPUT/m0_02"
```

They preserve their original protocol definitions, including the stale descriptive JSON string disclosed in RESULT_M0_02.md. The code and task list determine the actually executed controls. Do not replace a protocol's original verdict with a rerun using a changed dictionary.

## Optional original DMRG validation

This uses Perplexity's existing solver and is separate from the independent bosonic reference. Install its recorded dependency and invoke the corrected command:

```sh
python -m pip install physics-tenpy==1.1.1
python "$GRUT_REPLAY_ROOT/research/reality_framework/grut_qed_detector_histories_20261107/mps/validate_vs_ed.py" --output "$GRUT_OUTPUT/dmrg_validation.json"
```

The original four-case gate now returns 0 for pass and 1 for fail, writes outside the archive and records reference hashes. [Validation repair](VALIDATION_REPAIR.md) preserves the first wrapper exception, corrected clean-checkout replay and deliberate rejection of a modified reference. This is evaluator-entry-point verification, not a frozen external test of an evaluator on unfamiliar theories.

## Compact archive and provenance

From a clean committed tree, build the portable source/results packet:

```sh
python "$GRUT_MATCH/build_replication.py" --output "$GRUT_OUTPUT/GRUT_interacting_history_packet.zip"
```

The ZIP includes the new packet and the original QED sources/caches it depends on, a SHA-256 manifest and extraction instructions. It includes no Python packages or Git history. Its local bootstrap commit is a new replica identity; the manifest and individual raw files preserve the original freezes and source hashes. Put environments and outputs outside the extracted tree.

[Standalone validation](PACKET_VALIDATION.md) records a fresh virtual-environment run using only an extracted packet plus its declared dependencies. The small massive point and full D1 calculation reproduce the archived probabilities exactly in that environment, including the failed scientific gates. The optional `diagnose_field_moments.py` performs a labeled post-result analysis of archived ground-state checkpoints; it does not alter a registered comparison.

| Registration | Published before its new production |
|---|---|
| M0-01 | `2b13287ae398d45bed7b4dcdcbb083747aaa4a00` |
| M0-02 | `c5e107a9920c57b13482d63f290fa9f7bc1fc7b5` |
| M1-01 | `479eaed2f906e323dac491f8f308c0c05718d63b` |
| D1-01 | `90745ef0d7aaf4c6d8fe3a11edc7ee0d12f7d339` |

Original failures, raw points, resource records, ground-state checkpoints and exact analysis code remain alongside each report. Registration establishes the sequence of this work; it is not evidence that the model describes nature. No mass scan, new cosmological abundance campaign, merge or publication is part of this pass.
