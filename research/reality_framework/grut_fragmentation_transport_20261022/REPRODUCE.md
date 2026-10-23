# Reproduce the bounded calculations

Parent: bd2cac32d224f70c18359722646983e5b758c7e0.
Run from the repository root. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8 were used. No network downloads are required to calculate.
The parent formation and empirical-test directories must be present.

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
python research/reality_framework/grut_fragmentation_transport_20261022/code/fragment.py --n 64
python research/reality_framework/grut_fragmentation_transport_20261022/code/fragment.py --n 96
python research/reality_framework/grut_fragmentation_transport_20261022/code/fragment.py --n 128
python research/reality_framework/grut_fragmentation_transport_20261022/code/fragment.py --n 96 --max-dt .01
python research/reality_framework/grut_fragmentation_transport_20261022/code/fragment.py --n 96 --seed 29
python research/reality_framework/grut_fragmentation_transport_20261022/code/transport.py
python research/reality_framework/grut_fragmentation_transport_20261022/code/verify.py
python research/reality_framework/grut_fragmentation_transport_20261022/code/plot.py
```

The runs took roughly 20--243 seconds each on the shared machine, with
some overlap/CPU contention. This is not a runtime benchmark. Adaptive
event positions and roundoff may change slightly across platforms;
compare resolved common states and stated tolerances rather than requiring
bitwise-identical compressed NPZ containers or timing metadata.

`fragment.py --output /absolute/path/tag` writes tag.json, tag.csv and
tag.npz. `--end 6` is useful for a common-state rerun; do not call such a
truncated run a completed dense-core simulation. The standard outputs
overwrite this packet's local generated files. Use an alternate output
prefix if preserving the published snapshot during replication.

The transport CSV contains one row per conditional source/response kernel;
JSON also contains particle spectrum, energy/abundance bookkeeping,
counterfactual massless-receiver coefficients and channel-only examples.
NR JSON component charges depend on the selected density threshold and
must not be imported as a final soliton mass function.
