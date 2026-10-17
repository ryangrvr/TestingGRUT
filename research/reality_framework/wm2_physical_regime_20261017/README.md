# WM2-S9 physical regime

Start with [REPORT.md](REPORT.md), [PHYSICAL_TARGET.md](PHYSICAL_TARGET.md)
and the frozen [DECISION_GATE.md](DECISION_GATE.md).

The central proofs are [FORMATION_ANALYSIS.md](FORMATION_ANALYSIS.md),
[BALLISTIC_TRANSPORT.md](BALLISTIC_TRANSPORT.md), and their
[integrated protocol](INTEGRATED_PROTOCOL.md). Read the conservative
N>=2^28 boundary and supplied-resource limitations with every summary.
[HAMILTONIAN.md](HAMILTONIAN.md) states the complete unchanged microscopic
law and the size-uniform norm bound used by the transport comparison.

[SWEEP_CERTIFICATE.md](SWEEP_CERTIFICATE.md) is the frozen production
specification; [SWEEP_IMPLEMENTATION.md](SWEEP_IMPLEMENTATION.md) records
its implemented residual, checked hypotheses and executed controls.
The physical, mathematical and methods ledgers are separate.
[DEPENDENCIES.md](DEPENDENCIES.md) lists supplied inputs;
[OBSERVABLES.md](OBSERVABLES.md) states what can reject this model.

Sources are preserved; no earlier campaign is rewritten. From this folder,
with Python, NumPy, SciPy and a C++ compiler for the frozen positive kernel:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py
```

The default verifies saved hashes, source pins and consistency only. Explicit
replays preserve originals in `replay/`:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py --analytic-controls
OPENBLAS_NUM_THREADS=1 python reproduce.py --sweep-controls
OPENBLAS_NUM_THREADS=1 python reproduce.py --primary
```

The last command runs only the same two fixed workers, up to 1200 seconds
each; it refuses overwrites. No default command runs trajectories, and
there is no six-event or fallback command. Timing/stops may vary across
machines. [RESULT_TABLES.md](RESULT_TABLES.md) records the final finite decision.
