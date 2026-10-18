# Bounded exact controls and reproduction

The charter was frozen at 614ceb2efae88c206c5fd24c987fe942e14cb516.
Proof/code were frozen at d83205b98ffec972d2d8ba7989071804ddf20cf7,
before primary execution. EXECUTED_SOURCES.json records nine source hashes.
No frozen proof or source changed in response to the controls.

One `python run_controls.py primary` and one identical-source replay ran.
Both exited zero and emitted identical CONTROLS.json, SHA-256
1db491f00205a2503c8219afe2ac54f22729e4e0f81e6e084e15dd4ee43706c2.
Primary wall time: 0.075634995 seconds; replay: 0.074597183 seconds;
combined: 0.150232178 seconds. Peak child RSS: 11,648 KiB.

Frozen ceilings: 25 seconds CPU/wall per child, 60 combined seconds,
1-GiB address space, 25-MiB files. No fallback, relaxed gate, paid compute,
installation, persistent process or parameter/time search ran.

All 49 named checks passed using Fraction arithmetic. The channel controls
check exact quadratic maximizers, the conditional matrix formula and its
Gram rewrite, cubic coefficients, normalization costs, and illustrative
completely positive success maps. These illustrative channels are NOT
counterexamples or observations of the actual B2 moving dynamics.

The separate N=4 kinematic enumeration has 571 valid word/charge/fuel
configurations. It classifies 904 directed plus, 904 minus, 348 conversion
and 382 geometry transitions. Forty-four plus hops form exchange edges.
It checks their unique gate, reverse pairing and annihilation of the
fuel-empty outcomes. Storage/receiver details can remove kinematic edges
but cannot introduce a new color permutation, so no full receiver model
is substituted or simulated. Geometry is not clamped in the physical law.

There were zero physical time trajectories. No actual moving success,
retention, coherence or crossing-access function was evaluated. Input
weight test points are algebra checks, not a claimed time-window grid
certificate. The continuous-input statements rest on the analytic
quadratic/cubic proofs, not sampling.

Run `python verify.py` for final hashes, all six input pins, 17 entries
of the parent B2-RB manifest, replay identity and resource ceilings. The
saved runner refuses duplicate primary/replay metadata. A separate replay
can use a fresh checkout of the source-freeze commit with the unchanged
EXECUTED_SOURCES.json, preserving this historical record.

MANIFEST.json covers every final file except itself; commit identity
pins the manifest. Publication is on an isolated branch, with no canonical
merge or historical overwrite. Shared-context challenge and exact finite
checks are not independent mathematical review or empirical validation.
