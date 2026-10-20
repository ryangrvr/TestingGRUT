# Post-collision persistence and azimuthal response

This is an additive follow-up on the isolated two-object research branch. The physical scalar couplings and original kernel are unchanged. A generated extension synchronizes linear m=1 and m=2 azimuthal perturbations with the same background Verlet integrator. Perturbations do not backreact on that background.

**Completed:** 31 additional evolutions cover persistence, assembly from Gaussian charge clouds, and an exploratory phase-to-charge record. All 89 source/numerical/decoder checks pass. The clouds of widths 4 and 6 form localized cores under the same full action; a rigid-field feedback control disperses. The selected m=1 response remains resolution-sensitive.

Read [REPORT.md](REPORT.md) for measured outcomes, the algebraic feedback mechanism, source-freeze timing and limitations. [SUMMARY.json](SUMMARY.json) contains comparisons and check outcomes. [RAW_OUTPUTS_MANIFEST.json](RAW_OUTPUTS_MANIFEST.json) indexes the exact [raw archive](raw_outputs.tar.gz).

![Conditional assembly and charge record](RESULTS.svg)

To inspect a fresh checkout, run `python3 unpack_raw.py`, then `python3 analyze.py`. `python3 plot.py` recreates the figure using Matplotlib. Run `python3 build_kernel.py` and `node controls.cjs` to replay implementation checks; `node formation/controls.cjs` checks the Gaussian preparation and rigid control. The original physical kernel and each study's frozen evolution sources remain unchanged.

`node run.cjs RUN_ID`, `node formation/run.cjs RUN_ID` and `node record/run.cjs RUN_ID` replay individual cases listed in the corresponding protocol.json files. Their run_suite.py scripts replay complete suites. The selected canonical perturbation norm can grow through translations, deflections and phase shifts; growth alone is not a nonlinear stability or chaos theorem. The record uses an ideal observable; it does not implement a detector or establish consciousness.
