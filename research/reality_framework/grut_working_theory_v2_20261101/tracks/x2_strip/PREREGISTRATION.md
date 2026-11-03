# M2-X branch determination, re-run blind: pre-registration

Committed **before** any scan in this track is run. Claude, for the owner (issue #13).

**Why.** In lf_chain the M2-X cutoff (e_crit ≤ 4×10⁻⁶) was written after part of the data had been seen
(lf_chain NOTES §2). This re-run fixes the rule, the grid and the counting before any new output exists. Its result
is the final Q2 answer for M2-X under the lattice-fermion closure.

## The data

- **Fresh scans,** written to this track's `results/`. Nothing is reused from lf_chain's scan files. The engine is
  deterministic, so previously computed keys will reproduce; they are recomputed anyway.
- **Code unchanged:** e2e_chain's `scan.py` with lf_chain's redirection (`lfc_common`), run through `code/x2scan.py`.
- **Histories:** lf_chain's final LF histories (KFmid).
  - Computed nodes [C]: `hLF84_mid` (8.4×10¹⁰ GeV) and `hLF_mid` (1.3×10¹¹ GeV).
  - [A] nodes, which reuse a lattice block computed at another H: `hLF55_mid` (5.5×10¹⁰) and `hLF270_mid` (2.7×10¹¹).
- **Grid:** δ ∈ {−4, −2, 0, +2, +4}°, and c0 ∈ the union of
  - {1.00, 1.04, 1.08, 1.12, 1.16, 1.20},
  - {1.30, 1.31, …, 1.45},
  - {1.355, 1.365, 1.375, 1.385, 1.395}.

  That is 27 values per δ and 135 keys.

## The acceptance rule (fixed now)

- **Root.** For each key, take the free-X root with lf_chain's branch-aware rule `fp_lf.root(lo, hi, 'free X')`,
  linear in ln H between the two computed nodes, with the code unchanged.
- **M2-X point (strict, the final answer).** The root must be:
  1. self-consistent (the interpolated S/Q(10⁻⁶) is on the free-X side);
  2. not extrapolated: 0 ≤ w ≤ 1, i.e. H_* inside [8.4×10¹⁰, 1.3×10¹¹] GeV;
  3. within **1.03×10⁻⁶ ≤ e_crit (bulk) ≤ 4.12×10⁻⁶**. This is a factor of 2 either side of M1r's
     e_crit = 2.06×10⁻⁶, a value fixed before M2 existed.
- **M2-X strip.** A run of consecutive grid c0 values of strict M2-X points at one δ, reported with its end points.
  A single isolated point is reported as such. If there are no strict M2-X points, the final Q2 result is "**no M2-X
  strip** within the simulated range".
- **Extended reading (reported, not counted).**
  - Roots that use the [A] nodes or are extrapolated, with fp_lf's piecewise rule over all nodes.
  - The same e_crit band.
  - For comparison, the earlier cutoff (1×10⁻⁶, 4×10⁻⁶].
- **Also reported, under the same strict rule:** the M2-Q points (self-consistent Q-matter roots, not extrapolated).
