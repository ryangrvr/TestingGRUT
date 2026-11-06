# Replay X10

Run from the TestingGRUT root. The recorded runs used
Python 3.12, NumPy 2.3, SciPy 1.17 and GCC with OpenMP.
Use one BLAS thread per simulation process.

~~~bash
export OPENBLAS_NUM_THREADS=1
TRACKS=research/reality_framework/grut_working_theory_v2_20261101/tracks
X10="$TRACKS/x10_handoff"
X10_RUN=$(mktemp -d /tmp/grut-x10.XXXXXX)

python -B "$TRACKS/x9_memory/code/production.py" "$X10_RUN/production" "$X10_RUN/production.json" --only half --workers 2
H55="$X10_RUN/production/half_base.history"
python -B "$X10/code/h55_seed.py" "$X10_RUN/h55_seed.json"
python -B "$X10/code/h55_probe.py" "$H55" "$X10_RUN/h55_probe.json"
python -B "$X10/code/h55_live.py" "$H55" "$X10_RUN/h55_live.json"
python -B "$X10/code/controls.py" "$X10_RUN/controls.json"
python -B "$X10/code/handoff.py" "$X10_RUN/handoff" "$X10_RUN/handoff.json" --workers 3
python -B "$X10/code/retention.py" "$X10_RUN/retention" "$X10_RUN/retention.json"
python -B "$X10/code/hybrid.py" "$X10_RUN/hybrid" "$X10_RUN/hybrid.json" --workers 3
python -B "$X10/code/verify.py" --replay "$X10_RUN"
~~~

The physical probes consume the early numerical columns
through $tH_i=0.103$ and compare their digest. A previously
replayed X9 half-step baseline history can supply the same
input. The reconstruction also checks the digest of the
archived scalar diagnostic. The remaining trajectory
data stay in the work directory.

For a spot replay, generate the two physical probes, seed
and controls as above and use the following commands for
the handoff and hybrid stages.

~~~bash
python -B "$X10/code/handoff.py" "$X10_RUN/handoff" "$X10_RUN/handoff.json" --only reference --workers 1
python -B "$X10/code/hybrid.py" "$X10_RUN/hybrid" "$X10_RUN/hybrid.json" --only reference bath_holdout --workers 2
python -B "$X10/code/verify.py" --spot-replay "$X10_RUN"
~~~

Running the verifier without a replay argument checks
the archived evidence and input hashes.
[Inputs](results/inputs.json) · [Verification](results/verification.json)
