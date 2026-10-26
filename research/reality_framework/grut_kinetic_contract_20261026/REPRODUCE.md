# Reproduce this input audit

From the repository root, with Python, NumPy and SciPy installed:

```bash
python -B research/reality_framework/grut_kinetic_contract_20261026/code/audit.py
python -B research/reality_framework/grut_kinetic_contract_20261026/code/verify.py
```

The scripts read the already committed 225028b source JSON and write only
this packet's derived JSON. The reviewer implementation is inspected and
credited; it is not imported, modified or executed by this packet.
No release, quantum-mode production or relic scan is rerun.

Actual work: two analytic source transforms with 161 samples each; two
fixed-energy equilibrium comparators; five temperature examples; one
independent time-domain decay ODE; 240 on-shell collision-product samples;
16 quadrature/Bessel species comparisons; and 40 finite-Q derivative checks.
All five verification groups pass. Numerical residuals are in
[VERIFICATION_RESULTS.json](VERIFICATION_RESULTS.json).

The first verification completed its mathematics but stopped while saving
a NumPy boolean to JSON. The save conversion was fixed; the subsequent
console print needed the same adjustment. A clean verification then exited
successfully. No scientific input or tolerance was changed to obtain a pass.
An initial coordination-table assertion also stopped before editing; the
claim was then added specifically to the live-claims block.

Only the source's negligible finite-start cutoff and nonrelativistic-parent
limit are used in deriving the Rayleigh spectrum. The displayed reference
momentum is collisionless redshifting to a defined density coordinate.
The equilibrium comparator uses fixed vacuum masses and an explicitly linear
condensed energy omega*Q. It is not a finite-temperature action, a formation
law or a relic prediction.
