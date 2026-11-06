# Replay X8

Run from the TestingGRUT root with Python, NumPy, SciPy and GCC
with OpenMP. The full replay uses the X7 fine and coarse
collisionless H55 histories generated from the archived Bloch
checkpoint. [Input hashes](results/inputs.json) pin that
checkpoint and the inherited source and stress implementation.

~~~bash
export OPENBLAS_NUM_THREADS=1
TRACKS=research/reality_framework/grut_working_theory_v2_20261101/tracks
X8="$TRACKS/x8_unitary"
X8_RUN=$(mktemp -d /tmp/grut-x8.XXXXXX)
python -B "$TRACKS/x7_regulator/code/overlap.py" "$X8_RUN/x7" "$X8_RUN/x7.json"
python -B "$X8/code/controls.py" "$X8_RUN/controls.json"
python -B "$X8/code/source.py" "$X8_RUN/source" "$X8_RUN/source.json" --workers 6
python -B "$X8/code/source.py" "$X8_RUN/paired" "$X8_RUN/paired_source.json" --workers 6 --paired
python -B "$X8/code/source.py" "$X8_RUN/uv" "$X8_RUN/preparation.json" --workers 6 --preparation
python -B "$X8/code/transport.py" "$X8_RUN/transport" "$X8_RUN/transport.json"
python -B "$X8/code/memory.py" "$X8_RUN/x7/coherent_half.history" "$X8_RUN/memory.json" \
  --coarse-history "$X8_RUN/x7/coherent.history"
python -B "$X8/code/verify.py" --replay "$X8_RUN"
~~~

For a shorter replay run the X7, controls, transport and memory
commands above, then run these two commands.

~~~bash
python -B "$X8/code/source.py" "$X8_RUN" "$X8_RUN/unused.json" --workers 2 \
  --only paired_reference uv_reference
python -B "$X8/code/verify.py" --spot-replay "$X8_RUN"
~~~

The solver saves full sampled histories in the supplied work
directories. The tracked JSON files record the case settings,
comparisons, conservation errors, and gate decisions. Running
the verifier with neither replay argument checks the archive.
