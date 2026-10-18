# Adversarial review route — prepared, not dispatched

One assistant performed the analysis and challenge with shared code/context.
There is no independent reviewer, separate implementation, external clearance
or empirical evidence. Pinned B1/B2 source is reused read-only. This packet
does not validate B1's all-size formation/receiver extension.

## Check the decisive claims in this order

1. Reconstruct the unchanged B2 773-state matrix and the exact invariant
   F_L+F_R+N_+=F_L+F_R+N_-=2. Check that q=0 means BOTH fuels excited,
   while either single-fuel marginal is a sum of exclusive sectors.
2. Multiply the sector projector/current identities in `IDENTITIES.md`.
   Verify all complex coherences are kept until after evolution; no rate
   process or measurement is inserted.
3. Check [H_convert,L_e]=0 and that only hopping crosses the endpoint
   resource boundary. Recompute the theta=2 to 6 fuel/matter/resource
   balance. Try to account for the observed refill without spatial inflow.
4. Verify P's dimension 38, its <=4+2J energy bound, QHQ>=9/2-A-5/48-2J,
   and the initial variance 2J^2. Check the exact rational high-band bound.
   Do not confuse high spectral weight with bare charged-configuration weight.
5. Independently reconstruct source return from all spectral phases.
   Challenge the numerical block grouping using the archived residual
   account. Pinching is explanatory only; do not turn its favorable output
   or a finite-window average into a passed B2 protocol.
6. On the interior fixture, directly assemble the 72-state colored matrix.
   Check H=A tensor I+B tensor SWAP and triplet/singlet evolutions. Verify
   the cross-Gram channel, input-dependent return instrument, triplet code,
   and the claimed recovery operations and their extra resource costs.
7. Recompute gamma and its error account. Its departure from one is small;
   a nonzero exchange word was not evidence of catastrophic decoherence.
   No approximate identity assertion for the actual growing geometry follows.

## Numerical reproduction

From repository root, with Python/numpy/scipy already present:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python research/reality_framework/grut_b2_diagnosis_20261018/diagnose.py \
  --out-dir /tmp/b2dx-review
```

This writes the ten frozen data files and prints execution metadata. The
fixed law, initial states, times and outputs are pinned before execution.
`verify.py` checks file/source hashes and accounting without trajectories.
Raw files give complete population/coherence/current observations; the
law and code reconstruct all physical wavefunction coefficients.

`CHANNEL_ACCOUNT.json` is post-execution analytic processing of the frozen
triplet/singlet propagation errors: delta_gamma<=eT+eS+eT*eS+1e-12. It is
not a new primary evolution. All numerical error statements remain
conditional on the pinned B2 binary64 arithmetic contract. The spectral
pinching decomposition is a numerical diagnostic, not independently
certified decimal mathematics.

## What remains open after a successful reproduction

Reliable two-composite operation on increasing hosts, full growing-geometry
channel coefficients, practical decoder resources, extensive matter,
autonomous fuel/record preparation, selected dimension, gravity and a
distinctive empirical prediction. An infinite-volume trapping theorem is
not inferred from the discrete eight-event spectrum. A successful replay
does not settle any of those questions or constitute external validation.
