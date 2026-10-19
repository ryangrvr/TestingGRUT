# GRUT-M1-W candidate research packet

Start with [REPORT.md](REPORT.md), then the complete law in
[GRUT_CANDIDATE_M1.md](GRUT_CANDIDATE_M1.md). The outcome is a narrower
working matter–memory model, not a universal or empirically validated theory.

- [WORKED_REGIME.md](WORKED_REGIME.md): derivation, prepared record,
  readout, energy and retained failed attacks.
- [COMPATIBILITY.md](COMPATIBILITY.md): what inhabits the same law and
  what is still missing, especially gravity.
- [PREDICTIONS_AND_TESTS.md](PREDICTIONS_AND_TESTS.md): calibration and
  rejection conditions.
- [OPEN_CLAIMS.md](OPEN_CLAIMS.md): unresolved claims and known failures.
- [SOURCES.md](SOURCES.md): primary-source context and actual read scope.
- [EXECUTION.md](EXECUTION.md): numerical methods, interruption and replay.

The core results, stronger attacks and writer each have separate protocols
written before the corresponding calculation. Numerical outputs preserve
failed claims; no universal Lyapunov exponent was fitted.

## Reproduce

Use the recorded Python/NumPy/SciPy environment. All dependencies were
already installed; no package was installed for this work. To preserve
execution records, calculation scripts refuse to overwrite their default
outputs. Run the following in a fresh copy of this packet containing only
its Python files and protocol/amendment JSON files, with a different output
directory from the archived results:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python calculate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python challenge.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python algebra.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python autonomous_writer.py --dt 0.0002
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python writer_replay.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python assemble.py
```

The single-step-size original writer run only recreates the coarse
calculation; its self-comparison is not accepted as refinement evidence.
`writer_replay.py` supplies the distinct refined-time and spatial-grid
comparisons. Wall time depends on hardware; the first unoptimized coarse
writer took about seven minutes. None of the commands runs indefinitely.

From the repository root, verify the archived packet without rerunning the
physical calculations:

```sh
python research/reality_framework/grut_memory_candidate_20261019/verify.py
```

The finalizer verifies source/protocol identities, raw acceptance outcomes,
relative links and file hashes. It preserves false scientific attack
conditions; an integrity pass does not turn them into successful physics.
