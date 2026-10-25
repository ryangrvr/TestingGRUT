# Reproduce

Python, NumPy, SciPy and Matplotlib were already available. From the repo:

```bash
python3 research/reality_framework/grut_z6_gauge_bridge_20261025/code/matching.py
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_z6_gauge_bridge_20261025/code/phase.py
python3 research/reality_framework/grut_z6_gauge_bridge_20261025/code/verify.py
python3 research/reality_framework/grut_z6_gauge_bridge_20261025/code/report.py
python3 research/reality_framework/grut_z6_gauge_bridge_20261025/code/manifest.py --check
```

This executes 48 complex matching cases, 512 independent angular Fourier
samples, 100 unreduced gauge/kinetic checks, and 36 prescribed phase
histories plus one tighter-tolerance history. It does not rerun the parent
radiation BVP, a full modified-action soliton, formation or a cosmological
population. The late-rate number is an explicitly conditional rescaling
of the parent response, with the new matched coefficient.

The phase history's radial envelope, Higgs value and expansion are supplied.
No abundance, entropy, inflationary preparation or initial phase is selected.
The parent scalar parameters are retained as low-energy matching targets;
their UV quantum thresholds have not been computed.
