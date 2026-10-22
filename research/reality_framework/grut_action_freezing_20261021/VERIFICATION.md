# Restored execution: actual verification record

22 October 2026. One author. Not outside scientific replication.

The restored build passed **7 analytic controls and 28 numerical gates**,
with four full-action integrations and six late scenarios each refined once.
All nine stored headline references from the original execution match at the
declared continuity tolerance. The old raw artifacts are not claimed to have
survived; this archive is the regenerated execution.

A fresh temporary-directory run also passed all 35 checks, the complete
result leaf schema, and **237 scientific leaf comparisons** at 2e-6 relative
tolerance with zero absolute allowance. Timing, solver counters and near-zero
residual magnitudes are not demanded to be bitwise identical; new residuals
must pass the original gates. This is same-author re-execution.

The verifier checked the archive's exact 16-member set, sizes and hashes and
reconstructed metric kinetic energy, positive potential, full Friedmann and
scalar expansion-work accounting, and late energy/endpoints from raw CSV.

| Diagnostic | Restored execution |
|---|---:|
| Worst full Friedmann relative residual | 3.464e-9 |
| Worst full scalar-work residual / initial scalar energy | 6.234e-10 |
| Worst late scalar-work relative residual | 6.873e-15 |
| Full endpoint refinement change h/v | <4.5e-10 |
| Full endpoint refinement change log(phi/phi0) | <2.6e-15 |
| Late H endpoint refinement (km/s/Mpc) | <1.5e-14 |

Commands executed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/build.py
python3 code/plot.py
python3 code/make_manifest.py
python3 code/verify.py --rerun
```

Python/JSON syntax, relative links and git whitespace were checked. The figure
was opened and visually inspected. The final manifest was regenerated after
text edits and followed by a read-only verification. The complete local staged
tree must match the published tree before completion is reported.

No new empirical likelihood, quantum/thermal simulation, reheating receiver,
false-vacuum lifetime, external scientific replication, canonical merge or
adoption of QV3 as the project's target is claimed. The physical uncertainties
listed in LEDGER dominate the numerical errors.
