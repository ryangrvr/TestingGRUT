# Reproduction

Use a complete checkout of this research branch. Dependencies: Python 3, NumPy, SciPy. Frozen sibling modules are read and byte-checked; no parent physics is edited.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/cp12_unmasked_formation/verify.py
```

This runs the parent verifier and fresh independent small-matrix controls, then checks saved run records and source/dependency hashes. It does not silently rerun expensive pilot evolutions.

To rerun the measured cases:

```bash
P=research/reality_framework/cp12_unmasked_formation
OPENBLAS_NUM_THREADS=1 python "$P/run.py" --n 5 --budget 10 --probe --wait .1 --lag .1 --steps 1024
OPENBLAS_NUM_THREADS=1 python "$P/run.py" --n 6 --budget 6 --probe --wait .1 --lag .1 --steps 64
OPENBLAS_NUM_THREADS=1 python "$P/run.py" --n 6 --budget 6 --end .2 --steps 128
OPENBLAS_NUM_THREADS=1 python "$P/run.py" --n 7 --budget 3 --end 2 --steps 64
```

The six-event probe takes several minutes on the recorded host and roughly 1 GB RAM; resource use varies. The current probe code also records joint geometric/persistence diagnostics; the original six-event and budget-eight five-event probe pilots predate those extra output fields. The source-free seven-event pilot likewise predates the added sequential persistence fields. Their unchanged response/formation kernels and declared settings are reproducible; timestamps, runtimes, roundoff, and newly added diagnostics may differ. These pilot schemas are identified by verify.py, not represented as fresh runs of later-added diagnostics.

An optional CP12_FORMATION_CACHE environment variable names a LOCAL trusted scratch directory for intermediate Python pickle caches. Never load caches from an untrusted source. Omit it for uncached reproduction; caches are not research deliverables. All JSON output goes beside run.py. The verifier writes CHECKS.json and VERIFICATION.json locally; --record-manifest additionally refreshes the code/dependency manifest. Update saved run records intentionally when changing numerical settings.

REPORT.md states the failed precision and size gates. DERIVATION.md supplies the error bound, preparation accounting, exact finite-size saturation argument, and diagnostic definitions. No result here certifies emergent dimension or a causal propagation cone.
