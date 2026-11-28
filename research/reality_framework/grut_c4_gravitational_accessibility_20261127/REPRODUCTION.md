# Reproduce the C4 gate

Use the C4 draft on C3 head 0d1ff1e428fe49ffbdb55e5083da5f9e1df0d13c. Keep the C1 profile archive and results and the C2 and C3 protocols and results in their sibling directories. The local run used Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and the C2 extended-precision interval code. It inherits the host round-to-nearest and outward nextafter assumptions from C2. Hashes of the pinned inputs are in the results.

From this directory use a fresh output path.

~~~sh
timeout 45 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONDONTWRITEBYTECODE=1 python3 -u verify_gate.py \
  --output /tmp/grut-c4-replay.json > /tmp/grut-c4-replay.log
cmp RESULTS.json /tmp/grut-c4-replay.json
cmp EXECUTION.log /tmp/grut-c4-replay.log
~~~

The retained result and execution log replay byte-for-byte. The code reconstructs polynomial input data from the earlier archive and runs bounded quadratures and outward bounds. It solves zero new background, spectral or tidal systems. A sparse checkout needs all four C-series packet paths present. [Frozen protocol](PROTOCOL.json) · [Proof and units](ANALYTICAL_GATE.md) · [Initial non-scientific failures](NEGATIVE_RESULTS.md)
