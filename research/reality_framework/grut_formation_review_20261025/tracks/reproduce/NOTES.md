# Track `reproduce`: reproduction audit of ff49e5d (causal receiver), 9636950 (soft receiver) and the retention numbers of bb211d1

25 October 2026 · Claude Code review track · read-only on Codex's branches; no commits.

Two kinds of numbers appear below, always labelled:

- **REPRODUCTION** means Codex's own code, unmodified, run from a `git archive` extraction of
  ff49e5d (which contains the unchanged 9636950, bb211d1 and c4b99f3 packet trees). The extraction is in
  `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/reproduce/tree`.
  The pristine published copies are next to it in `ff49e5d/`, `9636950/` and `bb211d1/`.
- **INDEPENDENT** means my own code in `code/` and a different method from Codex's.
- **NEW CONTROL** means Codex's old executable driven by my own driver at parameters Codex did not run.

## 1. Environment, compiler and deviations from REPRODUCE.md

| Item | Codex (ENVIRONMENT.json) | This reproduction |
|---|---|---|
| Compiler | g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0 | g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 |
| Flags | -O3 -std=c++17 (hard-coded in run.py, run_interface.py, refinement.py, operator_controls.py) | see below |
| Python / NumPy / SciPy | 3.12.14 / 2.3.5 / 1.17.0 | 3.13.16 / 2.5.3 / 1.18.1 (Matplotlib absent, so report.py and plot.py were not run) |
| CPU | not recorded | 4-core Intel Xeon 2.1 GHz, x86-64 baseline (SSE2, no FMA), shared with other review tracks |

**Flags [C].** The brief asks for -O2. `code/gxx_O2_wrapper.sh` sits first on PATH, rewrites Codex's `-O3` to
`-O2` and logs every compile line. -O2 builds were used for all of these:
- radial_nodal.cpp and radial4.cpp (bb211d1);
- operator_controls.cpp;
- the ff49e5d smoke run;
- the opening segments of P00625 (to τ=2810), P0125 (to τ=8060), P025_remote and P025_remote_time.

For this code the -O2 build is about 3.4 times slower than -O3 (19.4 s against 5.65 s on the same 150-unit case).
The machine was also shared with another review track, so at 02:31 the wrapper was switched to keep Codex's
`-O3` for radial_cfs.cpp and radial_flux.cpp only. The cut-short -O2 runs were rerun from scratch at -O3.
The switch rests on byte-identical outputs from the two builds in three checks:
- smoke: all four raw-array hashes equal the published hashes;
- a 150-unit nonlinear run at h=0.25: final, history and memory CSVs are byte-identical;
- the -O2 opening segment of P00625 (282 history rows and 12035 trace rows) is bitwise equal to the -O3 rerun.

This agrees with expectation: x86-64 has no FMA by default, ISO C++17 mode turns off contraction, and GCC never
reassociates floating-point sums without -ffast-math. The full compile log is in the scratch file `gxx.log`.

**No run length was reduced.** Every reproduced run used Codex's parameters and tmax.

**Reproduction gaps found [C].** These are not physics errors.
1. bb211d1's `refinement.py` never creates `data/`, and REPRODUCE.md does not say to create it. Both C++ programs
   write through unchecked `std::ofstream`. When the directory is missing, a run computes for its full length,
   writes nothing, and fails only when Python tries to read the files afterwards. I lost the first N00625/H003125
   attempt this way.
2. The nine old-spatial-method H runs of bb211d1 were **not rerun** (CPU budget). They are not inputs to any
   claim audited here. bb211d1's `convergence.py` needs all 18 runs, so its spatial and control statistics were
   recomputed from the nine N runs by `code/compare_reproduction.py` instead.
