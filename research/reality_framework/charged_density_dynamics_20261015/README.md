# Charged-density dynamics, 15 October 2026

Entry points: `REPORT.md`, `DERIVATION.md`, `HANDOFF.md`.

Requires Python, NumPy, SciPy, and a C++17 compiler for the accelerated basis.
The compiler produces a temporary shared library from the included source.
No compiled binary or private cache is required in the repository. The
legacy Python canonicalizer remains an independent reference.

Run from the repository root. Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1.
An optional GRUT_DENSITY_CACHE directory caches constructed bases/matrices;
omit it for a cache-free basis build. The reference-builder control optionally
uses CP12_FORMATION_CACHE. Every trajectory/current is computed anew.

```bash
python research/reality_framework/charged_density_dynamics_20261015/run.py --n 4 --q 4 --budget 256 --a 1 --b .25 --times .5,1,2,5 --wait 5 --lag .5 --source-current --history
python research/reality_framework/charged_density_dynamics_20261015/checks.py
python research/reality_framework/charged_density_dynamics_20261015/run.py --n 5 --budget 30 --times .25,.5,1,2,3 --wait 2 --lag .25 --source-current --history
python research/reality_framework/charged_density_dynamics_20261015/run.py --n 5 --budget 50 --times .25,.5,1,2,3,5,10 --wait 2 --lag .25 --source-current --history
python research/reality_framework/charged_density_dynamics_20261015/run.py --n 6 --budget 8 --reference-builder --times .25,.5,1 --wait .5 --lag .25 --source-current --history
python research/reality_framework/charged_density_dynamics_20261015/run.py --n 6 --budget 12 --times .25,.5,1,2 --wait .5 --lag .25 --source-current --history
python research/reality_framework/charged_density_dynamics_20261015/response.py --n 5 --budget 50 --wait 2 --lag .25
python research/reality_framework/charged_density_dynamics_20261015/onset.py --n 5 --budget 30
python research/reality_framework/charged_density_dynamics_20261015/onset.py --n 6 --budget 12
python research/reality_framework/charged_density_dynamics_20261015/exact_boundary.py --n 5 --python-reference
python research/reality_framework/charged_density_dynamics_20261015/exact_boundary.py --n 6
python research/reality_framework/charged_density_dynamics_20261015/local_onset.py
```

The N6 budget-twelve basis has 715,985 retained and 705,212 outgoing orbits;
its fresh build takes roughly eleven minutes in the recorded environment.
The source-fixing N5 response basis has 810,377 retained orbits. These are
larger computational tests, not constant-cost demonstrations.

`verify.py` checks the packaged evidence, source hashes and cross-calculation
identities. It does not replace rerunning the commands or reviewing the
analytical bounds. `VERIFICATION.json` records execution scope and limitations.
