# Reproduce C3

Start from the C3 draft and retain the preceding results and profile files. The code uses Python 3.12, NumPy 2.3.5 and SciPy 1.17.0. Run from this directory with fresh outputs.

~~~sh
timeout 90 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -u scatter.py --output /tmp/grut-c3-results.json > /tmp/grut-c3-execution.log
cmp RESULTS.json /tmp/grut-c3-results.json
cmp EXECUTION.log /tmp/grut-c3-execution.log
~~~

The variable-phase solution uses the two registered radii and differential equation tolerances. The full angular result uses Numerov with the two registered radial meshes and angular truncations. To control the singular centrifugal term, an individual partial wave begins where $\rho/h$ is at least $2(\ell+1)$ and at least the registered minimum. The inner seed is regular in the constant core potential. Matched vacuum integrations subtract numerical vacuum phase drift. The two resolutions vary the origin handling, support and angular cutoff together. The result and log replay byte-for-byte in the retained environment.

The first invocation used the wrong C2 JSON key name and exited before a scattering calculation. Its traceback is retained. The reader correction singles out the archived certified radius under its existing key. No physical input or registered gate changed. [Initial execution](INITIAL_EXECUTION.log) · [Protocol](PROTOCOL.json) · [Results](RESULTS.json)
