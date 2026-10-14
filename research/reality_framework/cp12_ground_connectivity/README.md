# Running the CP-12 ground-connectivity packet

The scientific entry point is REPORT.md; proofs and symmetry conventions are
in DERIVATION.md. Frozen law and all physical constants are retained.

## Dependencies and ordinary verification

Python 3.11 or newer, NumPy and SciPy. Tested versions are recorded in
VERIFIED.json. From this directory:

~~~sh
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify.py
~~~

The verifier recomputes complete small spectra, exact rational
unbounded-rotor certificates, and uniform graph counts. It audits saved
large eigensolves and their independent-builder agreement; it does not
rerun every expensive solve. This distinction is recorded in VERIFIED.json.

The verifier explicitly inserts its own directory before local imports,
so it also supports Python isolated mode if NumPy and SciPy are installed
in that interpreter's normal site packages:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I verify.py
~~~

## Recomputing the large calculations

Run the following from this directory, with one BLAS thread. Direct solver
entry points below use ordinary Python, not isolated mode.

~~~sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python ground.py --N 4 --L 3 --full
python ground.py --N 5 --L 3 --full
python ground.py --N 4 --L 3
python ground.py --N 4 --L 4
python ground.py --N 5 --L 3
python ground.py --N 5 --L 4
python ground.py --N 6 --L 3
python ground.py --N 6 --L 4
python representation_ground.py --N 5 --L 3 --fresh
python representation_ground.py --N 6 --L 3 --fresh
python cutoff_comparison.py --N 5 --L 3 --z -0.728494 --fresh
python cutoff_comparison.py --N 6 --L 3 --z -0.85832 --fresh
python cutoff_certificate.py
python analyze.py
python verify.py
~~~

The all-representation scripts normally resume already completed sectors.
Use --fresh as shown to actually recompute every requested eigenproblem.
Basis caches are kept under a temporary-directory cache, or a directory
set by CP12_CACHE_DIR. Use a new cache directory after changing physical
basis/transition code. Caches contain only reproducible intermediates and
are not part of the committed scientific record.

The six-event all-sector calculation reduces a 36,801,265-state Hilbert
space to multiplicity blocks; its largest block has 817,489 entries.
This remains a substantial CPU/memory calculation.

The retained N=7,L=2 checkpoint can be regenerated separately:

~~~sh
python ground.py --N 7 --L 2
~~~

It is provisional: no all-sector check or cutoff refinement was run for
N=7. It must not be presented as a verified full unregulated ground state.

## Files and provenance

- ground.py: fixed-Q=N ground state, graph projectors, connected/disconnected
  matrix compressions retaining the boundary diagonal.
- representation_ground.py: every permutation representation, without
  assuming a singlet ground; fresh small full-spectrum checks.
- cutoff_certificate.py: exact rational small certificates for unbounded
  integer rotors, using the CP-12-specific tail bound.
- cutoff_comparison.py: larger floating all-sector Schur-comparison checks,
  explicitly distinguished from exact certificates.
- analyze.py: summaries from saved eigensolves plus fresh uniform graph counts.
- verify.py: declared verification scope and dependency versions.
- frozen_model.py, exact_destination.py, invariant_evolution.py: exact
  support-module copies from the previously published coherent-destination
  packet, commit 67fa7d4da7c9d6e0bc2330b6a764e1997957d455. Their older
  standalone entry points are not required for this packet.
- ground_*.json, representations_*.json, cutoff_*.json, analysis.json:
  raw and derived numerical records with explicit sector and cutoff scope.
- VERIFIED.json: freshly executed verification groups and limits.
- MANIFEST.json: file hashes and pinned provenance for this publication.

Only this new research directory is added on a separate TestingGRUT branch.
