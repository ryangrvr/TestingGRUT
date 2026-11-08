# E1-01: joint corrections improve histories; the bound still misses

7 November 2026 · Finite-parent method result · Zenodo paused

**Changing the Hamiltonian and state without also matching the detector leaves most of the tested history error.** The jointly derived correction reduces the maximum discrepancy by factors of 3.36 and 3.93 in the two declared detector regimes. It does not certify the interacting continuum limit or pass every gate within the finite parent.

The [derivation](DERIVATION.md), [frozen protocol](PROTOCOL_E1_01.md) and [implementation](experiment.py) were committed at `e6dc9b3b1222239c22695a1c03d589d39d66143c` before the dressed histories and weak-detector experiment were evaluated. [Raw summary](raw/summary.json) · [source provenance](raw/provenance.json).

## Physical claim tested

Can a derived 72-state description reproduce the existing 1280-state interacting parent's detector histories better when its preparation, Hamiltonian and detector are transformed consistently?

The parent is exactly the old E=12 Schwinger interval calculation. The retained space contains the same 72 free-energy states in all three schemes. The correction uses a block resolvent and an isometry J, with no fitted output or parent ground eigenvector supplied to the construction. It computes H_eff=J^dagger H J and W_eff=J^dagger W J. The Hamiltonian-only control uses the same H_eff and the same corrected preparation, but leaves the detector at PWP. It therefore has exactly the same spectrum as the jointly corrected model.

| Detector strength | Bare maximum TV | Hamiltonian-only maximum TV | Joint maximum TV | Joint bound |
|---|---:|---:|---:|---:|
| Strong, c=2 | 0.00769698 | 0.00706778 | 0.00229074 | 1.0 |
| Weak, c=0.5 | 0.00051905 | 0.00046705 | 0.00013216 | 0.12537856 |

TV is total variation of complete outcome histories, including overflow. The two histories are at times (0,0.7,1.9) and (0,0.6,1.5,2.4). The maximum is taken over both. Weak c=0.5 is a newly declared experiment; it does not replace the original strong detector.

## Frozen decisions

- **Method-improvement gate: passes in both regimes.** Joint TV is less than half the bare TV and below the Hamiltonian-only TV.
- **Finite-parent 0.001 accuracy gate:** strong misses; weak passes by measured comparison.
- **Finite-parent 0.001 bound gate:** both miss. The bound correctly dominates the measured error but is much too conservative. A measured weak-case pass is not relabeled a bound-certified result.
- Construction, normalization, Hermiticity, contraction and parent-replay checks pass. The archived strong parent is recovered within TV 2.22e-15.

All three descriptions use the same retained dimension, so the improvement is not obtained by silently increasing the basis. The eliminated block is still computed explicitly: this prototype is not a speedup and has not integrated out states above the original parent cutoff.

## What the state and spectrum show

The anchor iteration converges in four steps. Its distance from the eliminated block's spectrum stays above 6.38, and ||T|| is 0.14665. Those are diagnostics of this construction; they do not prove a uniformly separated excited-state band or a convergent continuum expansion.

The corrected preparation agrees with the parent's ground vector within 2.33e-15, and its energy within 3.22e-15. This follows from solving the self-consistent finite block problem, not from a new physical prediction about vacuum energy. The first-gap error improves from 0.01325 to 0.00464. The full retained-space dynamical residual is still 0.75167 in the declared units, even though its action on the ground state is negligible.

Thus exact ground-state matching and a common corrected spectrum do not fix later detector histories. The comparison between the joint and Hamiltonian-only cases isolates the effect of the detector dictionary within this specified reduction. It does not show that every observable or every reduction needs equally large corrections.

## Error accounting

The new bound includes three sources: the prepared vector's error, dynamical leakage out of the isometric image, and the mismatch between parent and reduced instrument actions. Its propagation contribution is bounded by an exactly integrated residual Gram matrix and Cauchy–Schwarz; an independent direct quadrature agrees with that matrix formula to 2.78e-17 in the declared verification case. There is no omitted time-quadrature error hidden in the certificate.

The bound discards cancellations when accumulating vector errors and summing branches. Its size is therefore not an estimate of the actual TV. It remains a valid but loose conditional finite-parent bound, evaluated with ordinary floating-point checks rather than interval arithmetic. It supplies no bound on the parent's missing continuum modes.

The overflow conventions are preserved. Ordinary outcomes are never renormalized. Detector dressing is not an adjustable missing-energy reservoir or a new unitary exponential chosen to hide leakage.

## Connection to GRUT

The integrated framework requires a common physical preparation, dynamics and operational access when comparing descriptions. This calculation makes that requirement explicit through a common isometry, rather than adjusting three independent fits. It also derives the exact eliminated-sector memory kernel, showing what an autonomous approximation discards.

These are applications of standard quantum and effective-Hamiltonian structure. They add no GRUT-specific law. M1-01's continuum matching remains unresolved, and the independent review of the bosonic boundary/mass dictionary remains pending. The result worth retaining is a measured improvement from a consistent correction, together with a recorded limit on what it certifies.

## Reproduce

Use a clean committed repository checkout and the existing packet's NumPy 2.3.5 / SciPy 1.17.0 environment. From the repository root:

```sh
export GRUT_EFFECTIVE_ROOT="$PWD"
export GRUT_EFFECTIVE_OUTPUT="$(mktemp -d)"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
timeout 600 python "$GRUT_EFFECTIVE_ROOT/research/reality_framework/grut_effective_instruments_20261107/experiment.py" --output-dir "$GRUT_EFFECTIVE_OUTPUT/e1"
```

The recorded run took 29.4 seconds with one BLAS thread. The output directory must not already exist. A zero exit reports a completed calculation; read the separate gates in `summary.json`. A timeout or exception is a preserved execution failure, not a scientific pass. The effective-Hamiltonian and embedding checkpoints, anchor iterations, histories, branch bounds and source hashes are included. Parent operators are regenerated from the frozen builder. The unchanged earlier packet still verifies its original decisions.

A second complete run from the clean frozen checkout in the isolated dependency environment took 39.6 seconds. All 16 compared distributions and all bounds reproduced at stored precision, with identical gate decisions. See [the replay record](REPLAY.json). This is replay of the same implementation, not independent external replication.
