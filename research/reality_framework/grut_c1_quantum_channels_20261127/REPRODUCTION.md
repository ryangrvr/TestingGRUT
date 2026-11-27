# Reproduce or verify C1

Use the C1 draft branch with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and a C++17 compiler with OpenMP. From this directory run the archive verifier with the Python standard library.

~~~sh
python3 verify_packet.py
~~~

The verifier checks the frozen and archived source hashes, the result flags, exact log integrity, signed intervals, the absolute subtraction and the workload caps. The replay uses 1,932 continuum Jost points and 34 finite-difference batches in 163.7 seconds, within the 2,400 point, 40 batch and 1,680 second caps.

Regenerate the six spline records from the frozen branch solver using a fresh filename.

~~~sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 verify_packet.py --regenerate-profiles /tmp/grut-c1-profiles.npz
~~~

Each coarse or fine charge group contains radial knots, local descending-power PPoly coefficients, and the scalar frequency. The trailing coefficient index enumerates $(F,F',S,S')$. The frozen solver interpolates the derivative variables as independent first-order BVP components. The profiles are numerical splines on their recorded finite domains.

A complete same-protocol replay uses a fresh output directory.

~~~sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 -u run.py --output /tmp/grut-c1-production
~~~

The wall times and order of parallel log entries vary between invocations. Compare the channels, registered gates and signed allowances. The archive verifier checks the exact bytes of the retained execution log after decompressing it.

The freeze commit is [a3ff9f4](https://github.com/ryangrvr/TestingGRUT/commit/a3ff9f40a1588501a24afc87de996668ef237374). [Protocol](PROTOCOL.json) · [Raw results](RESULTS.json) · [Profiles](PROFILES.npz) · [Execution log](EXECUTION.jsonl.gz)
