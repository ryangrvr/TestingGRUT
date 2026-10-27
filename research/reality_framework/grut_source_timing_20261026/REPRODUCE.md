# Reproduction and input chain

Dependencies: Python3, NumPy, SciPy, Matplotlib and a C11 compiler. The actual run used the installed runtime and gcc-O3, with two ordinary numerical worker processes. No AI subagents, paid compute, external data or installation was used.

The packet imports the pinned mode/preparation code already in its parent tree. The kinetic bootstrap may need the completed Claude builds, which reside on a different commit. Before reproduction from a shallow or single-branch clone, fetch the named commits:

```bash
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
git fetch origin bca4927214695ef7c668ebb1d20f26a9cba91ca7
```

bootstrap.py from34f23f8 extracts the immutable kinetic interfaces when its cache is absent. This packet does not edit those interfaces or replace the original source populations. The reheating loader reads the existing heavy-vector INPUTS.json; the new helper's distinct timing_common name avoids colliding with its common module. Set GRUT_TIMING_WORKDIR to relocate this packet's regenerable raw files. The inherited short matrix control uses GRUT_MODE_WORKDIR independently.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/verify_response.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/run_response.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/quartic_scan.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/source_clocks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/check_late.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/summarize.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_source_timing_20261026/code/plot.py
```

package.py without arguments checks the supplied sealed outputs, source hashes, syntax, finite arrays and local links. After an intentional numerical regeneration, package.py--seal records that regenerated packet; it does not establish bitwise equality with the original across different platforms. Compare the meaningful response/constraint/clock errors as well. The author sealed only after rendering and viewing TIMING.png.

Public profile archives are double precision and contain the homogeneous background, responses on every sampled momentum, and tangent vectors at .02,2.13,5.8. All full raw vectors remain regenerable; their actual binary hashes are in the response metadata. The archived linear vectors past a supplied seed's breakdown are not that seed's physical evolved state.

Construction corrections: the new quartic scan initially saved successfully and then failed while printing a NumPy boolean as ordinary JSON; printing the serialized file fixed the console output. Adaptive source quadrature warned at interpolation knots; it was replaced by piecewise Gauss integration, with orders4 and8 agreeing at3.4e-16. The old/current charge-sign convention is converted explicitly. Physical inputs, selected seed amplitudes and diagnostic thresholds were not retuned.
