# Replay the X3 calculation

Base repository commit: `9f71d34f8a768e4ec7bc60a46a65a833a818c1de`.
Read `SCOPE.md`, `RUN_PROTOCOL.md` and `REPORT.md` first. The original tracks
are inputs; none is rewritten. Python 3, GCC and the packages in
`requirements.txt` are needed. All calculations ran on CPUs.

To verify the archived bytes and inherited input identities without rerunning
the physics, execute the following from the repository root:

    python -B research/reality_framework/grut_working_theory_v2_20261101/tracks/x3_actual/code/manifest.py

This also reconstructs both engine source variants and checks all seven recorded
source hashes. A fresh numerical replay may change runtime metadata; the
archived manifest deliberately detects that change. Do not reseal changed
evidence while presenting it as the original archive.

To rerun the physics, from the repository root:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
X3=research/reality_framework/grut_working_theory_v2_20261101/tracks/x3_actual
X3_WORK=/tmp/grut-x3-replay
python -B "$X3/code/prepare_tables.py" "$X3_WORK"
python -B "$X3/code/moments.py"
python -B "$X3/code/kinetic.py" "$X3_WORK" --workers 3
python -B "$X3/code/quantum_controls.py"
python -B "$X3/code/populations.py"
python -B "$X3/code/verify.py"
```

The table replay must match the original `e2e_kinetic/results/ptable_summary.json`
to the stated tolerance. N1 uses the archived rounded mass ratio **0.41667**,
not a newly rounded rational substitution. The new spectra are constructed
directly on each destination grid, with exact discrete energy matching.

Each kinetic run saves its full input, time series and final state as lossless
gzip files with deterministic gzip timestamps. The compiled source hash is in
each case's JSON. `kinetic.py` copies the original C engine into scratch;
the half-step build changes exactly three adaptive step ceilings. This changes
no physical rate. No old output file is overwritten. The source's early-stop
condition ends the runs before 1e8, once condensates and singlet energy are
below 1e-9 of the network total. Late K drift is reported separately.

Quantum controls include a homogeneous line calculation and a true three-space-
dimensional, four-component **single-mode** Dirac evolution on 16^3, 24^3 and
32^3 grids. The latter is not a filled Dirac sea, H55 background, backreaction
calculation or cosmological closure certification. The memory estimate is a
Nyquist lower bound for a direct uniform-grid method, not a no-go for other
numerical methods.

Execution issues preserved: the first legacy table-rebuilder command imported
an unrelated symbolic model and stopped because SymPy was absent. A direct
table-only replay avoided that dependency and reproduced all four summaries
exactly. The first kinetic launch stopped before evolving because mpmath was
absent; mpmath 1.3.0 was then installed (SymPy 1.14.0 was also installed during
dependency recovery, but the X3 scripts do not require it). All seven physics
runs completed normally. These were environment failures, not failed physical
grid points. No unreadable or unsuccessful physical run was removed.
