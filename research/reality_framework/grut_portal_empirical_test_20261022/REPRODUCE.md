# Reproduction

Run from the TestingGRUT root at this packet's commit. Python 3.12.14
was used. The tested numerical environment is NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8, tqdm 4.67.1, ACROPOLIS 1.3.1, numba 0.62.1,
llvmlite 0.45.1. Use an isolated environment rather than changing an
existing project installation. For example:

```sh
python -m venv /tmp/grut-portal-env
/tmp/grut-portal-env/bin/python -m pip install numpy==2.3.5 scipy==1.17.0 matplotlib==3.10.8 tqdm==4.67.1 acropolis==1.3.1 numba==0.62.1 llvmlite==0.45.1
```

The official 33 MB compressed / 235 MB uncompressed toy likelihood
is external. The program's memory use is appreciably larger than its
compressed file size. Fetch the exact pinned data, then calculate:

```sh
TASK_PY=/tmp/grut-portal-env/bin/python
TASK_PACKET=research/reality_framework/grut_portal_empirical_test_20261022
TASK_CACHE=/tmp/grut-portal-inputs
export OPENBLAS_NUM_THREADS=1
"$TASK_PY" "$TASK_PACKET/code/fetch_inputs.py" "$TASK_CACHE"
"$TASK_PY" "$TASK_PACKET/code/prepare_inputs.py" "$TASK_CACHE"
"$TASK_PY" "$TASK_PACKET/code/detector.py" "$TASK_CACHE"
"$TASK_PY" "$TASK_PACKET/code/history.py"
"$TASK_PY" "$TASK_PACKET/code/shape.py"
"$TASK_PY" "$TASK_PACKET/code/photodisintegration.py"
"$TASK_PY" "$TASK_PACKET/code/verify.py" "$TASK_CACHE"
"$TASK_PY" "$TASK_PACKET/code/verify.py" "$TASK_CACHE" --replay
"$TASK_PY" "$TASK_PACKET/code/plot.py"
```

`fetch_inputs.py` checks the official tag's exact commit plus all nine
source byte counts and SHA-256 hashes before accepting them. It never
silently replaces an unexpected cache. It requires Git and curl.
The HEPData endpoints are versioned; if server serialization changes,
a hash failure is an input change requiring inspection, not permission
to weaken the check. Our source cache was temporary and is not needed
if these public files can be fetched.

`detector.py` writes the twelve-case detector results, full toy-coupling
table, reconstructed arrays and one observed likelihood curve. It
uses first-order portal amplitudes and the existing parent profiles,
not a fresh finite-kappa soliton solve. `history.py` writes the
twelve-case hot/cold results and reference trajectories. `shape.py`
writes the point-WIMP fits and reference spectrum comparison.
`photodisintegration.py` runs the original package four times, including
an effective-grid refinement; each run takes tens of seconds to a few
minutes depending on hardware. `--only same_envelope_grid_refinement`
replaces that one case in an existing result. Nuclear runtime fields
are metadata and will differ.

`verify.py` executes original official detector-code controls, distinct
formula controls and the source-tail test. Its `--replay` mode creates
a temporary isolated packet, reruns detector/history/shape and compares
all recorded numerical leaves plus raw tables/arrays. It does not
repeat ACROPOLIS, fit observational cosmological parameters or rerun
the old 51/99-check packets.

If intentionally modifying the delivered packet, regenerate its manifest
after outputs and documentation are final, then check it:

```sh
"$TASK_PY" "$TASK_PACKET/code/make_manifest.py"
"$TASK_PY" "$TASK_PACKET/code/make_manifest.py" --check
```

Source/data hashes and a clean replay establish reproducibility under
the stated assumptions. They do not establish an allowed particle model,
independent verification, novelty, blind predictive success or a solution
for the vacuum density.
