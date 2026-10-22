# Reproduce the proposal calculation

From the TestingGRUT repository root, with Python, NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_higgs_portal_feasibility_20261022/code/calculate.py
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_higgs_portal_feasibility_20261022/code/verify.py --replay
```

The calculator reads the parent RESULTS.json and three parent profile CSVs.
It writes this packet's RESULTS.json and one raw reference recoil spectrum.
It does not solve a new finite-portal saddle, run the old 51 checks or fit
experimental data. `verify.py --replay` copies the specified inputs and
calculator into an isolated temporary directory, executes it afresh, and
compares every scientific numeric leaf and the raw spectrum. The temporary
copy is removed when done; the archived outputs are not modified by verify.

If intentionally changing this proposal, refresh its manifest only after
all intended outputs and documents are final:

```sh
python research/reality_framework/grut_higgs_portal_feasibility_20261022/code/make_manifest.py
```

Fresh replay is same-author reproducibility, not an independent review,
parameter freeze, observed signal, allowed-region determination or canonical
theory adoption. All extrapolated couplings are formal until their validity
and the stated phenomenological requirements have been checked.
