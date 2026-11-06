# Replay X6

The X6 branch starts at X5 commit
`7a4c5aa9caf534d806564e89dced4cd66dda28c9`.
Run from the TestingGRUT root with Python, NumPy and SciPy.

```bash
export OPENBLAS_NUM_THREADS=1
X6=research/reality_framework/grut_working_theory_v2_20261101/tracks/x6_slab
X6_WORK=$(mktemp -d /tmp/grut-x6.XXXXXX)
python -B "$X6/code/source.py" "$X6_WORK/source" "$X6_WORK/source.json" --workers 6
python -B "$X6/code/source_followup.py" "$X6_WORK/source" "$X6_WORK/source.json" \
  "$X6_WORK/followup" "$X6_WORK/source_followup.json" --workers 4
python -B "$X6/code/transport.py" "$X6_WORK/transport.json"
python -B "$X6/code/bench.py" "$X6_WORK/controls.json"
python -B "$X6/code/edge.py" "$X6_WORK/edge.json"
python -B "$X6/code/verify.py" --replay "$X6_WORK"
```

The source runs save the full time series in the specified
working directory. The archived JSON files contain the
case settings, summary measures and all acceptance decisions.
The source replays take longer at the widest cutoff. The
verification can also be run without the replay flag to check
only the archived evidence.

The initial six-run comparison is stored in
`results/source_initial.json`. The intermediate published
checkpoint `5f91e6150447a0f8473281e525828ac969e9db7e` contains
that original solver and comparison. The later source tests
reallocate the known local subtraction between ledgers and
add cutoff controls, retaining the original outcome.
The [input hashes](results/inputs.json) pin the inherited
X5 stress and state evidence.
