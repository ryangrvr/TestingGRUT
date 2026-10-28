# Reproduce the fixed M1 stress tests

From the TestingGRUT repository root, with Python3, NumPy, SciPy, Matplotlib,
g++17 and OpenMP available:

```bash
export GRUT_M1_DYNAMICS_WORK=/tmp/grut_m1_dynamics_20261027
python research/reality_framework/grut_m1_dynamics_20261027/code/verify_engine.py
python research/reality_framework/grut_m1_dynamics_20261027/code/run.py
python research/reality_framework/grut_m1_dynamics_20261027/code/analyze.py
python research/reality_framework/grut_m1_dynamics_20261027/code/exterior.py
python research/reality_framework/grut_m1_dynamics_20261027/code/bulk.py
python research/reality_framework/grut_m1_dynamics_20261027/code/radial.py
python research/reality_framework/grut_m1_dynamics_20261027/code/radial_ladder.py
python research/reality_framework/grut_m1_dynamics_20261027/code/verify_radial.py
python research/reality_framework/grut_m1_dynamics_20261027/code/summarize.py
python research/reality_framework/grut_m1_dynamics_20261027/code/verify.py
python research/reality_framework/grut_m1_dynamics_20261027/code/package.py
python research/reality_framework/grut_m1_dynamics_20261027/code/figure.py
python research/reality_framework/grut_m1_dynamics_20261027/code/audit.py
```

Use an empty raw directory for a fresh replication. The initial-state constructor
and actual compiled evolution are in code/common.py and code/field3d.cpp. Individual
case names can be passed to run.py. The production histories have separate raw
directories; analyze.py collects their records into RUNS.json. Do not run multiple
collectors against the same analysis output simultaneously. The cache rejects a
source mismatch rather than silently accepting a stale physical run.

The 13 full 3D histories took about an hour of summed process wall times while
several ran concurrently on this host; per-case timings are in RUNS.json. The
nine spherical histories are much cheaper. Hardware, thread scheduling and runtime
columns vary. The short DOP853 comparisons are separate numerical implementations
by the same author, not independent-author replication.

TRACE_DATA.npz preserves every scalar diagnostic trace. RADIAL_DATA.npz and
EXCITATION_DATA.npz preserve radial histories and final profiles. FINAL_PROJECTIONS.npz
and OBJECT_FINAL_RADIAL.csv preserve figure inputs. The full 3D raw arrays are
reproducible intermediates, not durable exported full-volume checkpoints; their
initial/final/snapshot sizes and hashes are in RAW_MANIFEST.json. No lossy projection
is represented as an exact full-field restart. Analysis of full arrays requires
regenerating the named history in GRUT_M1_DYNAMICS_WORK.

SOURCES.json pins every inherited file consumed. The original M1 packet and older
reviewed engines are unchanged. No current observational fit, new inflationary
state, quantum collision operator or vacuum-density prediction is reconstructed
by these commands. Hot N64/N96 outputs are explicitly not spatially converged.
