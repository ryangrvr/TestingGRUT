# Reproduce the author calculations

Requires the parent 7737487 tree, Python, NumPy, SciPy and Matplotlib.
No new packages were installed. From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_gauge_point_test_20261025/code/pair_rate.py
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_gauge_point_test_20261025/code/point.py
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_gauge_point_test_20261025/code/full_release.py
OPENBLAS_NUM_THREADS=1 python3 research/reality_framework/grut_gauge_point_test_20261025/code/screening.py
python3 research/reality_framework/grut_gauge_point_test_20261025/code/verify.py
python3 research/reality_framework/grut_gauge_point_test_20261025/code/report.py
python3 research/reality_framework/grut_gauge_point_test_20261025/code/manifest.py --check
```

The pair calculation has three radial/phase quadratures for each power,
plus an analytic Gaussian control. It independently integrates the same
vacuum Born approximation as the reviewer, using the original profile.
It does not reproduce the review's new profile, distorted-wave coherent
solver or full nonlinear production campaign.

The high-field minimum uses nine starting points and an analytic-Hessian
polish. The full-field release uses three initial phases and one tighter
tolerance run, out to H_release t=100 (a about 14). It is deliberately a
short diagnostic, not a complete inflationary/fragmentation history or a
frozen asymmetry. The radiation metric, phi=0 regime and initial release
are stated inputs.

Screening uses three radial quadratures and a Fourier/Parseval control.
The extra energy is a variational binding bound, not a fluctuation spectrum.
The wall timing is estimated, not simulated. The added bias and curvature
coupling are proposed variants; they are not silently added to the old
point. The six gate dispositions are separate from passing author controls.
