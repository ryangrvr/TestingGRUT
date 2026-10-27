# Reproduction

This additive packet reads the prior pilot, onset and kinetic interfaces without editing them. Use Python3 with NumPy, SciPy and Matplotlib. Its raw inputs are the existing double-precision pilot snapshots; regenerate those with the8159556 reproduction instructions if unavailable. Set GRUT_3D_WORKDIR to their work directory when different from/tmp/grut_nonlinear_pilot_20261026. RUNS.json records the original absolute snapshot paths, so an external relocation must resolve them to the selected work directory; the audit code performs that resolution from the manifest case/name rather than requiring the original host layout.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_core_wave_audit_20261026/code/verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_core_wave_audit_20261026/code/audit_fields.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_core_wave_audit_20261026/code/summarize.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_core_wave_audit_20261026/code/plot.py
python -B research/reality_framework/grut_core_wave_audit_20261026/code/package.py
```

For deliberate regeneration of the packet's manifest/source/snapshot records, the author used package.py--seal after inspecting the actual scientific figure. A normal package.py call checks the sealed outputs and pinned source files without rerunning the3D histories.

The tests are same-author identity and control checks. Production post-processing covers27 unique saved states, with no new dynamical trajectories. No physical seed, heat source, bound population or final abundance is inferred.

During construction a NumPy boolean in console JSON output caused a serialization failure after the verification file had already saved; printing that serialized file fixed the output and the checks were rerun. No input or criterion was retuned. The new helper's distinct cw_common module name avoids a collision with the inherited onset module.
