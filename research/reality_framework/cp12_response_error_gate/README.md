# Reproduce the observable error gate

Use a complete checkout containing the parent `cp12_unmasked_formation` and `cp12_ground_connectivity` packets. Requirements are Python 3, NumPy and SciPy.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/cp12_response_error_gate/verify.py
```

The verifier reruns the parent controls and the new dense/complete-physical-matrix controls, then validates saved response records. The expensive six-event pilot is NOT implicitly repeated by that command.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/cp12_response_error_gate/probe.py --n 5 --budget 10 --steps 64
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/reality_framework/cp12_response_error_gate/probe.py --n 6 --budget 6 --steps 64
```

The six-event source quotient has 147,153 retained vectors and over four million matrix entries. Building it takes several minutes on the recorded host. Resource use varies. An optional CP12_FORMATION_CACHE can point to a trusted LOCAL scratch cache from the parent solver; omit it for an uncached run. Caches are disposable intermediates and not deliverables. Never load third-party pickle files.

The numerical method, Hamiltonian, initial ray, electric-square budget, and waiting/lag times match the parent benchmark. Backward residuals tighten the detector-error estimate. The global detector norm bound is one on the entire unregulated space, not an inferred maximum over retained states.

DERIVATION.md proves the bound and states its applicability condition. REPORT.md distinguishes a density response from a causal propagation cone, formation of extended space, and autonomous preparation. Bounds are analytically justified but evaluated in ordinary floating point. No interval-arithmetic certification is claimed.
