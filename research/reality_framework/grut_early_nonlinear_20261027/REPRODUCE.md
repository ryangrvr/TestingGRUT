# Reproduction

Python3, NumPy, SciPy, Matplotlib and g++ with OpenMP are required; gcc builds the inherited scalar/neutral linear engines. The installed runtime was used with no package installation, paid compute or AI subagents. Two ordinary numerical workers each use four OpenMP threads for the3D runs.

This branch contains all audited mode/preparation code from its parent. A shallow checkout may need the separate review/build objects:

```bash
git fetch origin ec1eb7f860f98f8d889f90ff94ebd53396f740bb
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
```

The kinetic bootstrap locates/extracts the pinned completed interfaces. None of the original packages is edited. The new nonlinear work directory defaults to/tmp/grut_early_nonlinear_20261027; override GRUT_EARLY_WORKDIR if needed. GRUT_MODE_WORKDIR and GRUT_TIMING_WORKDIR independently control the inherited, uniquely named scratch inputs.

Run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/dense_response.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/correct_clocks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/match.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/verify_engine.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/matching_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/production.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/linear_phi_control.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/analyze.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/export_covariance.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/summarize.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B research/reality_framework/grut_early_nonlinear_20261027/code/plot.py
```

package.py checks the provided sealed sources/outputs; package.py--seal records an intentional regenerated packet after figure inspection. A regenerated manifest does not prove bitwise equality across compilers/platforms. Check the reported physical/numerical comparisons rather than treating a hash mismatch alone as a scientific failure.

The public dense phase archives retain maxima and four selected momenta at.025tau. The nonlinear main48³ curve is retained at.05, the half-step curve at.05 from its.025 dynamics, and other curves are decimated for public viewing. All thresholds, conservation and work checks use the original dense traces, with their actual step limits in RUNS.json. Full fields and full traces are regenerable; their input/final/trace hashes are recorded. VIEW_FIELDS.npz is block-averaged and not used in quantitative analysis. SPECTRAL_EXPORT.npz is a measured classical covariance, not a physical particle distribution.

The32³ tau80 run in BENCHMARK.json is a performance/control run, not an eighth production history. The engine reuses the verified fourth-order spatial/temporal methods with the stated early potential and radiation clock. The matching code shares the prior audited linear equations; the3D short comparison uses a separate FFT/DOP853 implementation.

No failed physical point was retuned. Public arrays were compacted to selected phase vectors/curves after all scalar data were computed; raw production data remain regenerable. The finite-band state, reduced-model approximations and quantum matching limits are explicit in MATHEMATICS.md and PROTOCOL.md.
