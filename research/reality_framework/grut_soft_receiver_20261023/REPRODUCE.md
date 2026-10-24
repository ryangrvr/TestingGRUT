# Reproduce

Requires Python, NumPy, SciPy and Matplotlib. No new package installation was
used. From the repository root:

```bash
python3 research/reality_framework/grut_soft_receiver_20261023/code/boundary.py
```

That independent linear benchmark needs no raw formation traces.

The optional saved-trace audit requires six parent raw cases. If absent,
regenerate them using the parent builder (several minutes):

```bash
python3 research/reality_framework/grut_core_conversion_20261023/code/refinement.py N025 N0125 N00625 N0125_far N0125_box N0125_abs
python3 research/reality_framework/grut_soft_receiver_20261023/code/trace_audit.py
```

Those are the original preparations and are not a new receiver. The compact
packet includes every audit result and raw float64 hash, without duplicating
bulk traces. Review the recorded controls, then verify and render:

```bash
python3 research/reality_framework/grut_soft_receiver_20261023/code/verify.py
python3 research/reality_framework/grut_soft_receiver_20261023/code/report.py
python3 research/reality_framework/grut_soft_receiver_20261023/code/manifest.py --check
```

Use `verify.py --with-traces` when the parent raw files are available. That
option was used for the builder's recorded validation and also checks the
matched two-probe final fields and physical histories.

The verifier distinguishes accepted numerical controls from the failed
physical source/receiver gate. A passing verifier does not turn that gate
green. It does not claim a fresh replay of the parent full-field runs.
