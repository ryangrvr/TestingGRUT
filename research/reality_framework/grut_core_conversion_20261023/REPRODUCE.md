# Reproduce the packet

Use a checkout containing the published continuation and its parent tree.
Python dependencies are NumPy, SciPy and Matplotlib; C++ uses g++ with
C++17. No new package installation was used for this packet.

This is a compact publication: raw arrays are inventoried in
RAW_DATA_MANIFEST.json and retained in the builder workspace, but are not
included in this GitHub commit. A fresh checkout must regenerate the
18 core cases with the production commands below before executing the
data-dependent convergence, exterior and full verification commands.
Recorded JSON outputs and all code are included for immediate review.

After generating the raw arrays, from the repository root:

```bash
python3 research/reality_framework/grut_core_conversion_20261023/code/convergence.py
python3 research/reality_framework/grut_core_conversion_20261023/code/exterior.py
python3 research/reality_framework/grut_core_conversion_20261023/code/verify.py
python3 research/reality_framework/grut_core_conversion_20261023/code/make_manifest.py --check
```

The verifier reruns selected scattering solutions, action/force variations
and optical controls, and reconstructs every new raw trace. It does not
rerun all 18 production evolutions. The parent baseline's 13 checks were
executed in a fresh copy and recorded in BASELINE_REPLAY.json; no second
complete replay of the old eight production evolutions is claimed.

To recompute all new rates and conditional transport, in that order:

```bash
python3 research/reality_framework/grut_core_conversion_20261023/code/rates.py
python3 research/reality_framework/grut_core_conversion_20261023/code/transport.py
```

Do not run transport while the rate table is being regenerated. Each
transport batch holds one table in memory and records its numeric hash.
Timing fields can differ between reruns. These are same-author controls.

To rerun the 18 new core cases (several minutes to tens of minutes,
depending on the computer), keep separate case names as below:

```bash
python3 research/reality_framework/grut_core_conversion_20261023/code/refinement.py H025 H025_time H0125 H0125_time H00625 H003125 H0125_box H0125_abs H0125_late
python3 research/reality_framework/grut_core_conversion_20261023/code/refinement.py N05 N025 N0125 N00625 N0125_time N0125_box N025_time N0125_abs N0125_far
```

Regenerate convergence and exterior results after the last core finishes,
then run the verifier, plot script and manifest generator. Raw traces are
lossless compressed float64 arrays with their full-array SHA256 in the
case summaries. Histories and final fields are preserved for every case.

```bash
python3 research/reality_framework/grut_core_conversion_20261023/code/plot.py
python3 research/reality_framework/grut_core_conversion_20261023/code/make_manifest.py
```

The manifest pins compact packet files and consumed parent inputs. Passing it
checks integrity; it does not establish a physical completeness gate.
