# Replay X9

Run from the TestingGRUT root with Python, NumPy, SciPy and GCC
with OpenMP. Use one BLAS thread per simulation process.
The coevolved H55 runs also supply the early scalar
trajectories for the driven bath. The [input hashes](results/inputs.json)
pin the X5 background model and the X7 archive used to
select the six bath modes.

~~~bash
export OPENBLAS_NUM_THREADS=1
X9=research/reality_framework/grut_working_theory_v2_20261101/tracks/x9_memory
X9_RUN=$(mktemp -d /tmp/grut-x9.XXXXXX)
python -B "$X9/code/production.py" "$X9_RUN/production" "$X9_RUN/production.json"
python -B "$X9/code/envelope.py" "$X9_RUN/envelope" "$X9_RUN/envelope.json"
python -B "$X9/code/bath.py" "$X9_RUN/production" "$X9_RUN/bath" "$X9_RUN/bath.json"
python -B "$X9/code/bath.py" "$X9_RUN/production" "$X9_RUN/refined" "$X9_RUN/bath_refined.json" --refined
python -B "$X9/code/bath.py" "$X9_RUN/production" "$X9_RUN/dense" "$X9_RUN/bath_dense.json" --dense
python -B "$X9/code/live.py" "$X9_RUN/live" "$X9_RUN/live.json"
python -B "$X9/code/verify.py" --replay "$X9_RUN"
~~~

The first two bath archives include legacy whole-file
digests captured while the later part of the scalar
histories was still being written. All three archives
pin the early numerical columns actually used. The
verifier compares that trace digest during a replay.

A spot replay uses a freshly repeated coarse production
pair and the dense bath and live controls. The bath
runs can use the production histories already generated
by the full run.

~~~bash
X9_SPOT=$(mktemp -d /tmp/grut-x9-spot.XXXXXX)
python -B "$X9/code/production.py" "$X9_SPOT/production" "$X9_SPOT/production.json" --only coarse --workers 2
python -B "$X9/code/bath.py" "$X9_RUN/production" "$X9_SPOT/bath" "$X9_SPOT/bath_dense.json" --dense
python -B "$X9/code/live.py" "$X9_SPOT/live" "$X9_SPOT/live.json"
python -B "$X9/code/verify.py" --spot-replay "$X9_SPOT"
~~~

The full histories stay in the specified work directory.
The tracked JSON preserves the numerical comparisons,
including the intermediate misses and remaining cutoff
sensitivity.
