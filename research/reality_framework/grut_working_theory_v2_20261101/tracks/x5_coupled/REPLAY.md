# Replay X5

Use the TestingGRUT X5 branch with pinned X4 parent
`83ce5341a3819585d250f6a6bd6f3da6390496f0`. Run from the repository
root with Python, NumPy, SciPy, GCC and OpenMP.

```bash
export OPENBLAS_NUM_THREADS=1
X5=research/reality_framework/grut_working_theory_v2_20261101/tracks/x5_coupled
X5_WORK=$(mktemp -d /tmp/grut-x5.XXXXXX)
python -B "$X5/code/run.py" "$X5_WORK/base" "$X5_WORK/h55.json" --workers 6
python -B "$X5/code/run.py" "$X5_WORK/refine" "$X5_WORK/h55_refined.json" --workers 4 --threads 2 \
  --cases four_bins eight_bins wide_four_bins four_bins_half_step
python -B "$X5/code/controls.py" "$X5_WORK/controls" "$X5_WORK/controls.json"
python -B "$X5/code/spatial.py" "$X5_WORK/spatial.json"
python -B "$X5/code/phase_trace.py" "$X5_WORK/phase" "$X5_WORK/phase_trace.json" --threads 2
python -B "$X5/code/diagnose.py" "$X5_WORK/base" "$X5_WORK/refine" "$X5_WORK/diagnosis.json" \
  --state "$X5_WORK/h55_terminal_state.npz"
python -B "$X5/code/verify.py"
```

The main and refined histories take several minutes on multiple CPU
cores. Each history writes its full time series, terminal spectrum and
summary to the working directory. The pressure-work and moment tests
operate on those outputs. The final command verifies the archived
evidence in the checkout. Rounding and OpenMP reduction order can
change low-order digits on a separate replay.

## Saved spectral state

The archived `results/h55_terminal_state.npz` comes from the
5,760-bin evolution at the archived interface time. It stores the
one-dimensional comoving grid `q`, weights `w` per helicity, and
`bloch[species, helicity, momentum, component]`, with species N1–N4
and helicities −1 and +1. The axes are in the rotating convention
of MODEL.md. Scalar fields, velocities, scale factor, energies and
initial model parameters are in the corresponding `eight_bins`
record in `results/h55_refined.json`. The verifier reconstructs
particle moments, pressure and force from the saved state.

The original model inputs and their hashes are listed in
`results/inputs.json`. The 64³ phase replay uses an instrumented
scratch copy of the old lattice engine. Its added diagnostics leave
the field update equations intact.
