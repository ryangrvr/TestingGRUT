# Reproduce the background, modes, coupled check and credited proxy replay

Requires Python, NumPy, SciPy, Matplotlib and GCC. Run from the repository
root. These are independent numerical processes, not additional AI agents.
Typical wall time is several minutes with four CPU workers; the refined
review-code relic replay took about six minutes in this environment.

```bash
python research/reality_framework/grut_heavy_vector_budget_20261025/code/background.py
python research/reality_framework/grut_heavy_vector_budget_20261025/code/vector_modes.py
python research/reality_framework/grut_heavy_vector_budget_20261025/code/budget.py
python research/reality_framework/grut_heavy_vector_budget_20261025/code/coupled.py
```

To replay the corrected benchmark, prepare the review source without
changing this checkout:

```bash
git fetch origin f7d641a2570f148b67e4ad05bd8acbea5f97c960
mkdir -p /tmp/grut_review_f7d641a
git archive f7d641a2570f148b67e4ad05bd8acbea5f97c960 research/reality_framework/grut_formation_review_20261025/tracks/ls_threshold | tar -x -C /tmp/grut_review_f7d641a
python -B research/reality_framework/grut_heavy_vector_budget_20261025/code/replay_ke.py /tmp/grut_review_f7d641a/research/reality_framework/grut_formation_review_20261025
python research/reality_framework/grut_heavy_vector_budget_20261025/code/verification.py
python research/reality_framework/grut_heavy_vector_budget_20261025/code/plot.py
```

The extracted review is read-only for this work; bytecode generation is
disabled in the replay. `replay_ke.py` imports the pinned review's lsrelic
and lsthermo. It does not represent a second independent relic solver.
The primary vector and coupled equations are the author's own code.

Temporary integrator builds and raw traces use
/tmp/grut_heavy_vector_budget_20261025, overridable through
GRUT_HEAVY_VECTOR_WORKDIR. Published raw-trace SHA-256 hashes in the
background summaries record the actual runs; the primary Higgs norm is
also saved in HIGGS_PROFILE.npz. Long-profile controls are reconstructed
by the first command.

Actual computation counts:

- Four standalone scalar releases: a=40 central/tight, a=60 and a=100.
- 232 Proca mode histories, including controls and one gauge sensitivity.
- Two coupled three-mode controls, one tighter; 17 three-mode band batches,
  49 unique band momenta and two padding slots.
- Three actual review-code proxy replays, including 600/1200-node ke+s.

The source integrator window boundary can emit a negligible-width final
window; it is excluded from averages. Adjacent output samples sometimes
map to the same floating-point ln(a); the interpolation collapses those
duplicates. The initial attempted interpolation stopped on this issue
before any mode results were saved. The corrected production mode run is
the one recorded here. The review's original-proxy Radau replay emitted
an internal overflow warning; it returned success and finite results.
Both ke replays had zero reaction-root fallbacks. Very small proxy f values
are interpreted as erased, not as meaningful 180-digit relic predictions.

The Coulomb-gauge derivation, independent algebra checks and late-mode
projection limits are in MATHEMATICS.md. No full all-mode interacting
quantum state, phase-holder action or kinetic formation calculation is
claimed. MANIFEST.json hashes final packet files except itself.
