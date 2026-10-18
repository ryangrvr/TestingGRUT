# Candidate 1 stability follow-up

Start with [REPORT.md](REPORT.md): one focused result connecting the unchanged scalar action, a charged-state family, specific fragmentation costs, three-dimensional finite-time persistence and conditional observables. This packet extends [Candidate 1](../grut_candidate_1_20261018/REPORT.md) without replacing its frozen evidence.

## Reproduce without overwriting the original experiment

From this directory, choose a new path that does not exist:

```bash
python reproduce.py --output /absolute/path/to/new-reproduction
```

The wrapper copies scripts, frozen minimization prescriptions and the three original background/result inputs into a fresh sibling-directory layout. It runs the finite branch scan and refinement, twelve independent minimizations, three bounded Newton checks, seven nonlinear trajectories, analysis, verification and figure sequentially. It does not install dependencies or contact accounts. On the original hardware the large three-dimensional runs take minutes; timings are recorded, and performance on other machines will differ.

Use `--prepare-only` to inspect the clean layout before executing calculations. The wrapper refuses an existing output directory. The independent optimization/refinement scripts also refuse existing result files so their recorded stops cannot be silently overwritten. Other component scripts produce outputs beside themselves; do not rerun them over the preserved archive unless deliberately working in a copy.

Tested dependencies: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8. The existing candidate's [requirements](../grut_candidate_1_20261018/code/requirements.txt) are sufficient. No new package was installed.

Reproduction means agreement at the stated numerical accuracy. Adaptive meshes, solver roundoff, compressed metadata and runtime fields need not be byte-identical. Fixed physical couplings and the actual stopping limits are part of the calculation.

## Evidence map

| Files | Purpose |
|---|---|
| `PROTOCOL.md`, `PROTOCOL_AMENDMENT_01.md` | Original protocol and pre-trajectory switch to a volume-preserving, analytically controlled strain |
| `DERIVATION.md` | Nonlinear equations, energy/charge, exact potential bound, preparation energy and conditional escape/fission statements |
| `branch_scan.py`, `branch_protocol.json`, `branch_refinement_protocol.json` | Fixed-frequency/fixed-charge BVP continuation, retained call-limit failure and bounded threshold/stability refinement |
| `branch_results.json`, `branch_states.csv`, `branch_profiles.npz`, `branch_thresholds_refined.json`, `branch_execution.jsonl` | Family, alternative starts, thresholds, every resolved component comparison and raw states |
| `branch_first_pass_results.json`, `branch_notes.md` | Preserved first-pass stopping limit and final interpretation |
| `independent_minimize.py`, `independent_protocol.json`, `independent_results.*`, `independent_profiles.npz`, `independent_execution.jsonl` | Separate finite-volume fixed-Q energy method, two starts, two grids, all original optimizer stops |
| `independent_refine.py`, `independent_refinement_*` | Separately frozen sparse-Newton stationarity check, preserving first-pass limits |
| `independent_analyze.py`, `independent_comparison.json`, `independent_notes.md` | BVP comparison and assumed-h² Richardson consistency, without global/interval claims |
| `evolve.py`, `evolution_prescription.json`, `evolution_execution.log` | Same-law 3D Hamiltonian Fourier/Klein–Gordon split evolution and run matrix |
| `evolution_CASE.json`, `evolution_CASE_slice.npz`, `evolution_run_summaries.json` | Time histories, conserved totals, localized-core diagnostics, numerical controls and final central planes |
| `analyze.py`, `analysis.json` | Matched-time error comparisons, exact input-energy controls and inherited direct observable calculations |
| `compare_slices.py`, `convergence_slices.json` | Direct final central-plane spatial/box comparisons, resampling and explicit band-truncation remainder |
| `verify.py`, `VERIFICATION.json` | Final acceptance evidence and explicit unverified claims |
| `OBSERVABLES.md` | Elastic s-wave contribution and its measurement assumptions; inherited full-stress weak-gravity values |
| `CHALLENGE.md` | Actual automated challenge, corrections and remaining limits |
| `plot.py`, `STABILITY_FIGURE.png/pdf` | Scientific figure generated from raw branch and evolution outputs |
| `reproduce.py` | Fresh-layout finite reproduction; preparation-only path was checked |
| `MANIFEST.json` | Hashes of this published snapshot and pinned input commit |

Original inputs come from `../grut_candidate_1_20261018/code/` at construction commit `c2e2cb6d7be1b68ce2b4370d5fb87b08412d5e8e`. The background hash is recorded in both the independent and evolution records. Those original inputs were not modified.

## Interpretation

Global charge conservation is a numerical symmetry check, not proof of localization. Local core retention is observed only over the finite sampled interval. Energetic comparisons cover found branches/components, not every field state. The finite-volume/Newton check reproduces stationary states through a separate implementation; it is not outside peer review. The reporter instrument and empirical particle realization remain open. No quantum or nonlinear-gravity evolution was executed.
