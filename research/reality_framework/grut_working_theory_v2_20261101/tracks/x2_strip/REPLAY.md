# Portable replay of the blind M2-X scan

The acceptance rule is the one committed at `6d836ca5051b7bf3e9d677e26c3b7b2db0ae0f6c`.
Start from the paused results in `800c3bfd21757f9d9826b6276686f6c3362be611`.
The additional scripts only reconstruct inputs, adapt paths, shard the
prescribed keys, and apply the unchanged `fp_lf.nd` / `fp_lf.root` functions.
They never use the older `lf_chain` scan rows as input.

From this track's directory, with Python, NumPy, SciPy, SymPy, mpmath and GCC:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
scan_work=/tmp/grut-x2-blind
python -B code/rebuild_history.py tables "$scan_work"
python -B code/rebuild_history.py hLF84_mid "$scan_work"
python -B code/rebuild_history.py hLF_mid "$scan_work"
python -B code/resume_blind.py computed "$scan_work" --jobs 8 --shards 4
python -B code/analyze_blind.py "$scan_work"
```

For the separate extended reading, reconstruct `hLF55_mid` and `hLF270_mid`
with the same history command, run `resume_blind.py extended`, then add
`--extended` to `analyze_blind.py`. These histories use a scalar lattice block
from a different inflation scale. Their roots never enter the strict count.

Rebuilding a history uses the archived, full-precision hand-over state at
`t=3` in `lf_chain/results/handover_ho<NODE>.json`. This avoids reconstructing
the lost reference scratch record that originally produced that state. The
dephased solver, moment occupations, splice and anharmonic spectator-mass
correction are the original implementations. The fermion transition masses
are the ones in the committed table summary: `0.41667`, `1.6`, `0.8`, `0.48`.
In particular, `0.41667` is the original tabulated input, not the exact fraction
`250/600`. The table integrals and each reconstructed history are checked
before the engine runs. No couplings or rates are adjusted.

Each shard uses a separate engine cache and result file. Re-run with the same
`--shards` after interruption to skip complete keys. `resume_blind.py progress
"$scan_work"` reports only counts. The original checkpoint rows are kept
exactly in the merge. The reconstructed baseline is compared to the archived
baseline for the computed nodes and used only as the engine reference for
missing points. The analyzer refuses a partial grid.

The strict rule uses the numeric bulk threshold in the inclusive
`[1.03e-6, 4.12e-6]` band, branch self-consistency and `0 <= w <= 1`. The
Q-matter table separately lists all self-consistent bracketed Q roots and the
subset that would also satisfy the bulk band. An isolated passing c0 is
reported as an isolated point. Neither interpolated roots nor the extended
reading constitute a new lattice simulation at their root scale.

The exact rebuilt benchmark bytes, their SHA-256 digests and the replay
checks are also archived in `results/restored_histories/`. To restore those
verified inputs instead of recomputing the tables and kinetics, use
`python -B code/restore_histories.py "$scan_work"` in an empty work directory.
This is a byte-for-byte restore, not an independent reconstruction.

After both computed grids finish, `python -B code/analyze_blind.py "$scan_work"
--seal-only` writes their lossless archives before any root evaluation.
The ordinary analysis command then refuses a grid that differs from its
sealed archive. Runtime-dependent fields such as `seconds` are preserved,
so whole-file hashes of a separately timed repeat need not coincide.
