# WM2-S7: positive residual bounds and shared-frame action capture

Continuation of frozen S6 (`b4c7902`), with a separately published
pre-trajectory charter (`4bf1448`). Only this folder is changed.

Read [REPORT.md](REPORT.md) for the result and its scope,
[ARITHMETIC.md](ARITHMETIC.md) for the arithmetic contract,
[MATHEMATICS.md](MATHEMATICS.md) for the action/projection analysis and
[HANDOFF.md](HANDOFF.md) for the next boundary. The charter is preserved
verbatim. No prior CP-12 or WM2 Hamiltonian is changed.

[EXECUTION.md](EXECUTION.md) discloses the exact run sequence and numerical
assumptions. [RESOURCE_MODEL.md](RESOURCE_MODEL.md) distinguishes exact
array counts from measured whole-process memory. The separate post-campaign
[weighted residual derivation](NEXT_RESIDUAL.md) was not applied to S7 runs.

## Reproduce

Python 3.12, NumPy, SciPy and a C++17 compiler are needed. No package install
or network access is used by the solver. From this directory:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py
```

This checks archived file hashes and numerical consistency; it does not
rerun trajectories. Explicit reruns write to `replay/`, preserving all
archived outcomes:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py --controls
OPENBLAS_NUM_THREADS=1 python reproduce.py --recovery --baselines --six
```

The second command enforces the five-event recovery gate before six events.
Each primary worker has the frozen 600-second wall budget. Runtime-dependent
stops can differ across machines, so an exact match of stopped populations
is not expected. Target times and physical instances cannot change with
runtime. Native source is compiled with strict flags into ignored `build/`;
compiled libraries and transient arrays are not archived.

`worker.py` also exposes individual fresh components, but the charter's order
and recovery gate govern any research interpretation. Original S5 code is
imported from the preserved sibling folder. The comparison engine is not
rewritten to make S7 look better; resource/error allocations differ and this
is not an equal-error speed comparison.

The default verifier checks consistency under declared arithmetic assumptions,
not independent validation of the mathematics, compiler, BLAS or physics.
