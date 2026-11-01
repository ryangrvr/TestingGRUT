# Track m1_vacuum_skeptic (pointer)

1 November 2026 · Claude · adversarial skeptic of `tracks/m1_vacuum/` (Codex d4747cb / ed73791, GRUT-M1 one-vacuum
potential and the Q = 200000 object), issue #13.

As instructed, the write-up lives next to the reviewed track:

- report: `../m1_vacuum/SKEPTIC.md`
- own code: `../m1_vacuum/skeptic/code/` (`sk_algebra.py`, `sk_qball.py`, `sk_hessian.py`, `sk_coldfluid.py`, `sk_wall.py`)
- results: `../m1_vacuum/skeptic/results/*.json` (`sk_checks.json` is the per-claim table)

Digest: no claim of the track is refuted. Every number was recomputed with a different method (exact rationals and
30-digit roots; a second-order cell-centred fixed-Ω solver with Richardson; a discrete-energy Hessian with a
secular-equation VK test; 6×6 dispersion matrix; planar-wall surface tension). Soliton: Ω and E/Q agree with Codex to
2.5e-14 and 5.8e-14, mass to 9e-13 relative. Corrections to the track: roll time 1.9059 → 1.9079/m_X (0.1 %); the
M1V-10 "CORRECTED" label is strong given Codex's own finite-grid caveat (interpretation confirmed: lowest fixed-charge
eigenvalue = free-box level to 2e-5 at L = 200, none below the edge, L₊ has one negative eigenvalue, f(0) = −21.5 < 0);
the Q^{-2/3} coefficient in the E/Q fit drifts 75 → 111 while the leading coefficient 8.164 (σ = 99.04 m_X³) confirms
the fitted 8.2. Nothing committed; no validation, exclusion or new theory claimed.
