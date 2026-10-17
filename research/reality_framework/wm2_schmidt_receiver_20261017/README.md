# WM2-S8 receiver compression

Read [REPORT.md](REPORT.md), the original [charter](CHARTER.md), and the
separately frozen [G1 allocation charter](GLOBAL_ALLOCATION_CHARTER.md).
[MATHEMATICS.md](MATHEMATICS.md), [ARITHMETIC.md](ARITHMETIC.md) and
[RANK_ALLOCATION.md](RANK_ALLOCATION.md) distinguish the numerical tools
from the unchanged microscopic model and the actual prepared state.

From this directory, with Python, NumPy, SciPy and a C++ compiler available:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py
```

The default checks saved hashes/records only. Explicit replays write into
`replay/`, preserving the archived outcomes:

```sh
OPENBLAS_NUM_THREADS=1 python reproduce.py --controls
OPENBLAS_NUM_THREADS=1 python reproduce.py --original
OPENBLAS_NUM_THREADS=1 python reproduce.py --global --baselines --six
```

The six-event replay requires the G1 recovery gate. Original S8 and G1
have separate gates; neither replaces the other's stops. Timed stops may
vary across machines. No default command launches primary trajectories.

Fresh S7 workers import the frozen S7 source and retain their 600-second
budgets; S8/G1 workers have 900 seconds, all fixed before their trajectories.
Report actual target completion and whole-process resources, not retained
coefficient counts alone. This is not an equal-error speed comparison.
