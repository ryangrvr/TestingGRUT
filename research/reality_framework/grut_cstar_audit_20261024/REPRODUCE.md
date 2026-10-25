# Reproduce the new checks

Requires the parent ff49e5d tree and Python, NumPy, SciPy and Matplotlib.
No new package installation was used. From the repository root:

```bash
python3 research/reality_framework/grut_cstar_audit_20261024/code/closure_checks.py
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_cstar_audit_20261024/code/charge_leakage.py
python3 research/reality_framework/grut_cstar_audit_20261024/code/verify.py
python3 research/reality_framework/grut_cstar_audit_20261024/code/report.py
python3 research/reality_framework/grut_cstar_audit_20261024/code/manifest.py --check
```

The forced radiation calculation uses six coupled production solves,
two free Green-function controls and two lower-resolution background
controls. It does not rerun the formation simulation. The pilot result is
retained as a labeled exploratory precursor, not an extra independent check.

All arrays/results needed for this compact calculation are included in the
parent and current commit; no bulk formation traces are required. Literature
read scopes are listed separately. Passing these controls does not promote
C*, the unspecified Wilson coefficients, the proposed gauge completion or
a vacuum-density selection into established results. Independent review
and a complete modified-action evolution are not claimed.
