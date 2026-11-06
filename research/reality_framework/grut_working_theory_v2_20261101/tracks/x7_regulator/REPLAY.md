# Replay X7

Run from the TestingGRUT root with Python, NumPy, SciPy and GCC
with OpenMP. Limit BLAS to one thread per simulation process.
The preregistered parent is
08b7e792c89dded6242a521c3cce581693a07e67.
[Input hashes](results/inputs.json) identify the previous source
and the exact archived H55 checkpoint.

~~~bash
export OPENBLAS_NUM_THREADS=1
X7=research/reality_framework/grut_working_theory_v2_20261101/tracks/x7_regulator
X7_RUN=$(mktemp -d /tmp/grut-x7.XXXXXX)
python -B "$X7/code/regulator.py" "$X7_RUN/static.json"
python -B "$X7/code/controls.py" "$X7_RUN/controls.json"
python -B "$X7/code/evolution.py" "$X7_RUN/series" "$X7_RUN/evolution.json" --workers 6
python -B "$X7/code/evolution.py" "$X7_RUN/series" "$X7_RUN/radial.json" --workers 6 --radial
python -B "$X7/code/evolution.py" "$X7_RUN/pulse" "$X7_RUN/pulse.json" --pulse
python -B "$X7/code/overlap.py" "$X7_RUN/h55" "$X7_RUN/overlap.json"
python -B "$X7/code/verify.py" --replay "$X7_RUN"
~~~

The radial follow-up reuses the uniform 0.25-bin run from the
first slab command. The scratch runs save full sampled histories.
The tracked JSON records test settings, conservation errors,
comparisons and acceptance decisions. For a shorter spot replay
run the static, controls and overlap commands above, then run
these two commands.

~~~bash
python -B "$X7/code/evolution.py" "$X7_RUN" "$X7_RUN/unused.json" \
  --only gauss_refined gauss_cutoff gauss_heavier --workers 3
python -B "$X7/code/verify.py" --spot-replay "$X7_RUN"
~~~

Running the verifier without a replay argument checks the
archived evidence and hashes. The test counts measure replay
agreement and preservation of the recorded gate failures.
