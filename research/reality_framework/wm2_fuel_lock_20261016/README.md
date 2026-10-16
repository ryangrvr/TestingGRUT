# WM2-S2: explicit energy storage and the resonance boundary

Separate construction campaign based on frozen WM2-S1
`7b1d9234deaeaddea3f34c419bdb35b99efd944b`.
WM1 remains `8741b6e74d571e643e9208187f5d93b5a8e2471c`.

S2's declared storage-assisted law forms a macroscopic path at exact
resonance with probability >=3/8 at one uniformly sampled time in a
specified O(N log N) window. All reverse terms remain. A newly derived
42.5 spectral gap and complete `(1/680)(6/91)^d` virtual bound hold in
the enlarged consistent sector. Full-H source responses and bulk
prefix persistence also have uniform bounds in that class.

The robust-phase gate **fails**: every fixed nonzero local storage
detuning suppresses macroscopic formation at all times as O(N^-2).
The stored energy is reversible and does not implement a one-way lock.
Dimension is supplied by the path grammar. A separate naive branching
proposal has a five-move shortcut that restores the full original tree.

Read REPORT.md for results, MODEL.md for the complete microscopic law
and preparation, MATHEMATICS.md for all proofs, and HANDOFF.md for the
review boundary. RESULTS.json and CONFIGURATION.json contain the saved
controls. ASSESSMENT.json records pass/fail status; VERIFICATION.json
records consistency and provenance, not physical acceptance.

From the repository root:

```bash
python research/reality_framework/wm2_fuel_lock_20261016/reproduce.py
python research/reality_framework/wm2_fuel_lock_20261016/reproduce.py --recompute
```

The first command checks saved evidence and frozen dependencies. The
second regenerates controls, report and scientific figure, then refreshes
the content manifest. Requirements: Python, NumPy, SciPy, matplotlib;
exact versions used are saved. No remote service is used by reproduction.

S1, WM1, frozen CP-12, SCOUT-1 and canonical GRUT documents are unchanged.
There is no original-charged-matter universality theorem, selected
three-dimensional phase, continuum gravity, earned quotient constant or
distinctive empirical prediction. No external validation is claimed.
