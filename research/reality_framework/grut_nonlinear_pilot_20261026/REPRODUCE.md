# Reproduce the cold pilot

The existing environment used Python3.12, NumPy, SciPy, Matplotlib and g++ with OpenMP. No dependency installation was needed. Run from the TestingGRUT root; the onset and kinetic packets are immutable input dependencies at2335653 and34f23f8, with the e71730c interfaces supplied by their bootstrap.

```bash
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/verify_engine.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/verify_analysis.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/production.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/analyze.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/track.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/expanded.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/summarize.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/association.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/negative_candidates.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/plot.py
python -B research/reality_framework/grut_nonlinear_pilot_20261026/code/package.py --seal
```

GRUT_3D_WORKDIR selects scratch output; the default is /tmp/grut_nonlinear_pilot_20261026. The C++ binary is compiled under a source-hash name with -O3, -fopenmp and -ffp-contract=off. Two ordinary numerical workers use four threads each; these are not AI subagents. Physics evolution is identical in every thread count; only reduction rounding in diagnostics can vary.

The initialization file contains magic GRUT3D1 followed by newline, int32N, float64L, four float64 parameters(v,g²,lambda_X,lambda_phi), and six double-precision arrays(q1,q2,phi,p1,p2,p_phi), in C order. Snapshots contain int32N, float64L, float64time, int32precision_bits, then the same six arrays. Headers use little-endian native values on the tested host. C-order Fourier functions have fixed coefficients across resolution controls; they are not newly sampled per grid.

Production writes traces at approximately unit time and snapshots at approximately20/m_X, with exact integer-step times retained. The final240 snapshot is additional if the preceding stored time is slightly earlier. Late summary values interpolate onto times160,180,200,220,240 to avoid giving near-duplicate endpoints extra weight. Snapshot scatter is not an independent ensemble uncertainty.

analyze.py caches completed runs only while its source hash matches. Run it again after all seven histories finish to append newly completed cases. Tracking and expanded masks must follow completed analysis. Final package checks require all seven cases and nine numerical groups. The dataset and hashes are described in DATASET.md; visualization arrays are block averages of densities and are not raw field checkpoints.

Before production, a C++ long/long-long template ambiguity was fixed, and the discrete gradient-energy diagnostic was expressed using differences to avoid cancellation on a constant field. The performance-only benchmark is not a production result. No seed, coupling, conserved target or physical conclusion was tuned to make these checks pass.

Reproduction can differ pointwise after long chaotic evolution across compilers/hardware. Use conserved quantities, early-time convergence and the stated statistical diagnostics; do not certify late populations merely from conservation or one matching snapshot. The code does not include injected heat, quantum kinetics or physical final f_sym.
