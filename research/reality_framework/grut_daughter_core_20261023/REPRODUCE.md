# Reproduce the point checks and radial prototypes

Run from the TestingGRUT commit containing this packet. The tested runtime
uses Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8 and a C++17
compiler invoked as g++. No additional package or system installation was
performed. The parent d71cc01 inputs must be present at their relative paths.

First verify the delivered bytes, before overwriting any outputs:

```sh
TASK_PACKET=research/reality_framework/grut_daughter_core_20261023
python "$TASK_PACKET/code/make_manifest.py" --check
```

The manifest check also loads every raw trace shard and verifies all eight
reconstructed float64 arrays against their originally recorded hashes.

In a copy of the packet, run the calculation:

```sh
export OPENBLAS_NUM_THREADS=1
python "$TASK_PACKET/code/daughters.py"
python "$TASK_PACKET/code/run_core.py" A A_time A_fine A_box B A_far A_finer A_remote
python "$TASK_PACKET/code/spectra.py"
python "$TASK_PACKET/code/verify.py"
python "$TASK_PACKET/code/plot.py"
```

`daughters.py` solves neighboring Q=199999/200001 profiles, computes the
annihilation/streaming screens and writes the new point record. The exact
charge yield and v calibration provenance are retained. Its separate fresh
temporary replay compared 346 numeric leaves; the result is recorded in
POINT_REPLAY_RESULTS.json. That is not an independent author or a new
rerun of all older source packages.

`run_core.py` compiles the included C++ in a temporary directory and
performs each case sequentially. It solves the supplied initial neutral
profile, evolves the full radial fields, and saves raw history/final fields
plus lossless float64 time-trace shards of at most 8192 samples each.
The most refined and longest runs each took a few minutes in this runtime;
timings and compiler/platform roundoff can differ. About 50 MB of compressed
trace/data artifacts are retained. No actual primordial-density 3D run is
performed. The numerical absorber has a measured energy/charge receiver.

`spectra.py` reconstructs exact traces with `run_core.load_trace`, measures
time/Fourier fluxes, and writes restricted particle-spectrum tables. Its
incomplete soft-field coverage is an output, not an exception to ignore.
No duplicate unsharded trace or regenerable full Fourier array is required.

`verify.py` distinguishes thirteen implementation/scoped-rate checks from
five physical-completeness gates. It intentionally prints OPEN for the
population gates; a successful program exit does not mean those gates
passed. The force control compiles the original C++ and compares its
first drift against a separate finite variation of the discrete energy.
All stored traces are rehashed; selected traces regenerate measured flux.
Free-wave, phase-space and decay-budget controls retain both charge signs
and the real-scalar normalization. There is no outside replication.

The source calculation was run in the eight stated production/control
cases, not independently replayed a second time in full. A full rerun
will reproduce the same mathematical problem, but floating-point sensitive
late populations should not be declared identical without comparison.
The recorded spatial convergence failure is retained as scientific evidence.
One additional A_box execution restored a trace shard found truncated in
the final delivery audit. Its complete trace hash and all scientific
summary fields exactly matched the original run; original production
timing remains in CORE_RESULTS.json, with both timings and the comparison
in TRACE_REPAIR_RESULTS.json. This was not an independent replication.

After intentional changes, regenerate/check the manifest only after all
outputs and documentation are final:

```sh
python "$TASK_PACKET/code/make_manifest.py"
python "$TASK_PACKET/code/make_manifest.py" --check
```

Absolute vacuum selection, in-core quantum capture/visible leakage,
asymptotic soft-particle extraction, a converged 3D charge population,
Lyman-alpha fitting and new detector likelihoods are outside this run.
