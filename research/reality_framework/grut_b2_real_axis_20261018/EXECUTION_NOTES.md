# Execution disclosure

Charter: 233b27781f89a617b3315481da5f180a4eec2724.
Proof/code freeze: 5e77eb0e9d4d950db68d5682069e95801d00011b.
Both were published before primary control execution.

The primary command was `python run_controls.py primary`, followed by
one `python run_controls.py replay` from unchanged frozen sources.
run_controls.py enforces 25 s timeout, 1 GiB address-space ceiling,
25 MiB child file ceiling and one recorded execution per mode.
No failed primary, fallback, rescue variant or post-result source edit.

Primary: 0.037979281 s, 9472 KiB peak RSS, exit 0.
Replay: 0.026736376 s, 9472 KiB peak RSS, exit 0.
Combined: 0.064715657 s; maximum 9.25 MiB.
Both output SHA256:
3ef33828b42d8fbfe8ae8fdb0a720e40bce4534940d579a6a475ff31dfb5cd8e.

All 110 controls passed. Arithmetic is Python exact Fraction/integer
arithmetic; square-root upper bounds are rational ceilings proved by
integer squaring. JSON float values are display only. No outward-rounded
floating-point claim is needed for these constant comparisons.

No Hamiltonian evolution, operating trajectory, spectral quadrature,
time/length search, new installation, paid compute, experiment, external
reviewer contact or persistent/background work occurred. Preliminary
work consisted of analytic derivation, source reads and syntax parsing;
no unreported numerical tuning runs preceded the primary.

The existing actual-source Challenger read the actual operator and the
new proof without executing the controls. Shared context, model and
derivation ancestry are disclosed in CHALLENGE.md. Its agreement is an
internal adversarial check, not independent mathematical validation.

Exact arithmetic validates only the displayed finite calculations. The
operator-theoretic argument and inherited all-size comparison remain
author-derived pending independent review. The final report and manifest
were written after the successful checks; frozen proof/code hashes are
preserved in EXECUTED_SOURCES.json.
