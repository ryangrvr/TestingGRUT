# Reproduction

From the TestingGRUT root at this packet's commit, using the existing
Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8:

```sh
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_formation_wall_audit_20261022/code/calculate.py
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_formation_wall_audit_20261022/code/verify.py --replay
OPENBLAS_NUM_THREADS=1 python research/reality_framework/grut_formation_wall_audit_20261022/code/plot.py
```

The calculator reads only the frozen parent's DETECTOR_RESULTS.json
and HISTORY_RESULTS.json. It does not require the large external toy
ensemble, ACROPOLIS or previous transient source caches. It writes
RESULTS.json and two raw CSV arrays. The plot reads these outputs.

The verifier uses separately organized full linear matrix/ODE controls
and physical peak/thermodynamic checks, then performs a fresh isolated
replay when requested. Its denser control densities improve numerical
conditioning for independent matrix diagonalization; the physical
dilute rate is computed directly from the exact scaled determinant.
The temporary replay tree is removed automatically. Neither program
simulates nonlinear collapse, injects guessed formation radiation or
reruns an older experimental likelihood.

After deliberate changes to calculations/outputs/documents, refresh
and verify the manifest:

```sh
python research/reality_framework/grut_formation_wall_audit_20261022/code/make_manifest.py
python research/reality_framework/grut_formation_wall_audit_20261022/code/make_manifest.py --check
```

The neutral wall comparison uses the stated 0.81 MeV primary bound
and its stable-network scope. No likelihood data are fetched or fitted.
Source credit and exact access scope are in SOURCES.md. Reproduction
of a conditional onset does not establish formation of the target
objects, empirical evidence or a vacuum-energy solution.
