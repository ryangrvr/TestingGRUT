# Reproduce

From the repository root, with Python, NumPy, SciPy and Matplotlib already installed:

```bash
python research/reality_framework/grut_five_outputs_20261027/code/calculate.py
python research/reality_framework/grut_five_outputs_20261027/code/verify.py
python research/reality_framework/grut_five_outputs_20261027/code/figure.py
python research/reality_framework/grut_five_outputs_20261027/code/audit.py
```

The calculation reads only the pinned heavy inputs, the preceding finite-history
matrices and the stored refined radial profile. It does not regenerate old histories,
run a kinetic formation solver, download detector data or install dependencies.

The primary recoil setup is frozen in PROTOCOL.md. The 27-point comparator-halo scan
is explicitly secondary and was added after viewing the primary result. No real
observations are fitted. The six verification groups check the legal state family,
complex matching, cold charge/mass ledger, gravitational constraint under constant
shifts, direct halo integral, and radial/recoil quadrature. They do not certify the
uncomputed formation, quantum matching, portal cosmology or a unique cosmic state.
