# Executed verification record

22 October 2026. One author, no subagents, outside replication or new installs.

The complete calculation passed **27 internal checks**. A fresh calculation in
a temporary directory then passed all 27 checks, the complete result leaf
schema and **290 scientific leaf comparisons**, at 2e-6 relative tolerance
with zero absolute allowance. The verifier also checked exact payload hashes
and independently reconstructed probability moments from five raw grid files.

| Check | Executed result |
|---|---:|
| Maximum canonical finite-derivative epsilon relative difference | 1.81e-11 |
| Maximum finite-derivative eta absolute difference | 6.74e-10 |
| Maximum probability normalization residual | 7.13e-13 |
| Minimum recorded probability | Positive |
| Largest boundary-cell mass | <1.2e-15 |
| Fine-grid positive second-moment change | <5e-5 relative |
| Larger-domain moment/probability change | <7e-13 relative |
| Mode ODE/Hankel power relative difference | <6.1e-7 |
| Wronskian relative error | <2.7e-8 |

Direct and eliminated densities agree under their common approximations.
Brownian half-line probabilities and second moments match the exact control.
The stochastic integrations use a sparse matrix exponential; grid/domain
controls are relevant, rather than a fictitious Monte Carlo time-step claim.

Commands run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/calculate.py
python3 code/plot.py
python3 code/make_manifest.py
python3 code/verify.py --rerun
```

Python/JSON syntax, relative links, git whitespace and the scientific figure's
render were checked. A final manifest regeneration/read-only check follows
text edits. Publication must match the exact locally staged tree.

These are numerical checks on a conditional slow-roll relation and a local
stochastic limit. They are not an all-order quantum theory, a converged relic
distribution through reheating, current-data likelihood exclusion, independent
scientific replication or positive identification of the proposed fields.
