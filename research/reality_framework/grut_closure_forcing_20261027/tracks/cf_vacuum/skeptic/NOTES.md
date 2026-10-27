# cf_vacuum skeptic: working folder

The review itself is `../SKEPTIC.md`. This folder holds its code and results only.

- `code/sk_rge.py`: own 3-loop SM RGE, transcribed by hand from 1307.3536 App. A, with two matchings: Buttazzo 2013
  ("B13") and Bednyakov et al. 2015 ("B15").
- `code/sk_line.py` → `results/sk_line.json`: validation against the primary sources, the PDG 2026 line, the
  tension, Λ_λ, and the scan over the Planck-suppressed c₆|H|⁶ operator.
- `code/sk_checks.py` → `results/sk_checks.json`: the hidden-sector P3a bound by direct grid minimisation, the
  CW threshold (sympy), the mixing and valley identities, the portal loop stabilisation, the invisible-width bounds,
  the hidden-sector Landau poles and the inflation numbers.
- `code/sk_summary.py` → `results/sk_summary.json`: key numbers and verdict checks.
- `results/sk_sources.json`: primary sources, read levels and searches.

Run, in about 3 minutes: `cd code && python3 -B sk_line.py && python3 -B sk_checks.py && python3 -B sk_summary.py`.
