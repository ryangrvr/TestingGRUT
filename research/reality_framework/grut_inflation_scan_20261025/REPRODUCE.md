# Reproduce the scan

Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Matplotlib and GCC were used.
The C solver generates its DOP853 coefficient header from local SciPy.
No raw arrays from old collapse runs or reviewer code are needed.

From this packet's code directory, run in order:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GRUT_SCAN_WORKDIR=/tmp/grut_inflation_scan_20261025
python build.py
python controls.py
python scan.py --workers 8
python refine.py --workers 8
python ensemble.py
python validate_candidates.py
python local_window.py
python analyze.py
python report.py
python verify.py
python manifest.py
```

The scan is computationally uneven: low H requires resolving millions of
massive-field oscillations. Eight CPU workers were used. Temporary full
trajectories and ensemble outputs are cached outside the repository;
physics-source hashes distinguish trajectory caches. Reproduction can use
an empty work directory. Runtime fields vary and are not scientific data.
The literal first calculation used scratch build/cache paths; these were
made configurable after the runs without changing evolution equations,
coefficients or scientific values.

The first local-map summary serialization failed on a NumPy integer;
it was corrected to a native integer and regenerated from the completed
trajectories. No physical run or point was discarded because of it.

The reported H/phase grid and local map are in SCAN_RESULTS.csv,
REFINEMENT_RESULTS.csv, REFINEMENT_SHIFT_RESULTS.csv and
LOCAL_WINDOW_RESULTS.csv. The 2,048 phase-ensemble members and canonical
endpoints are in ENSEMBLE_MEMBERS.csv; their grouped means and configuration
are in ENSEMBLE_SUMMARY.json. Original full temporary outputs have hashes
in that summary. FINDINGS.json separates numerical evidence from physical
status. CONTROL_RESULTS.json and CONTROL_HISTORIES.json preserve the
mathematical checks. Y_SOURCE_RESULTS.json gives every assumption in the
late-source estimate.

Independent reimplementation targets: the explicit H=1.1211167304e13,
c_H=1, delta=pi/2 point; the full inflationary amplitudes returned by
model.minimum; late N_X≈-238.649, E_X≈741.650 and J_Y≈13.269; and the
2e-12/2e-14 tolerance comparison. The 0.03-radian source-phase ensemble
uses N=128 and seed 2026102600. The actual random sequence and supplied
distribution are part of the experiment, not an inferred inflationary state.
