# Reproduce GRUT-M1

From the repository root, with Python, NumPy, SciPy, Matplotlib and gcc available:

```bash
python research/reality_framework/grut_matched_candidate_20261027/code/match.py
python research/reality_framework/grut_matched_candidate_20261027/code/deformation.py
python research/reality_framework/grut_matched_candidate_20261027/code/rematch.py
python research/reality_framework/grut_matched_candidate_20261027/code/consequences.py
python research/reality_framework/grut_matched_candidate_20261027/code/energy_stability.py
python research/reality_framework/grut_matched_candidate_20261027/code/cold_modes.py
python research/reality_framework/grut_matched_candidate_20261027/code/cold_loss.py
python research/reality_framework/grut_matched_candidate_20261027/code/verify.py
python research/reality_framework/grut_matched_candidate_20261027/code/figure.py
python research/reality_framework/grut_matched_candidate_20261027/code/audit.py
```

Set GRUT_MATCHED_WORKDIR to a task-specific temporary directory if desired. The wrapper
directs the previous engine's raw files and compiler output there, without changing any
old input or report. The candidate's compact background exports retain complete early
rows through tau30 and selected later rows. They are not a dense spectral/particle
source; full raw histories regenerate with match.py/rematch.py. verify.py uses the
published compact early rows for its short independent-method comparison and needs
no raw-history cache.

match.py runs9 homogeneous controls for the old profile's normalization. deformation.py
solves10 changed stationary profiles. rematch.py adds7 homogeneous histories for the
actual changed mass. energy_stability.py runs4 grid/domain controls; cold_modes.py
scans6 supplied cold densities, without picking an onset population. cold_loss.py
recomputes2 classical-response and2 Born-pair quadratures with the new profile and
neutral curvature. Eight mathematical/numerical groups must pass.

This does not run microscopic reheating, a physical phase holder, full gauged/quantum
kinetics, a three-dimensional changed-potential formation history, a population relic
solver, a detector likelihood or a fundamental vacuum calculation. Those are not
represented by completed flags in the machine-readable point.
