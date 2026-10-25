# Track gp_critic: notes

The main write-up is [CONSOLIDATED.md](CONSOLIDATED.md). It covers the verdicts, the six requirements, cross-track
flags, corrections, open items and the next calculation. These notes record method, limits and reproduction only.

## Method

- **Read.** b5174e5 REPORT, CLAUDE_HANDOFF, MATHEMATICS, CORRECTIONS, PAIR_RESULTS; 7737487 CONSTRUCTION and
  MATHEMATICS (via `git archive` into the scratchpad); the NOTES.md of all five gp_* tracks; tracks/affleck_dine NOTES;
  REVIEW_XCHARGE.
- **Own code** (no Codex code run):
  - `code/critic_checks.py` [D, C, E] covers:
    - matched alpha6, beta6 and c6;
    - the survival rate and its robustness to higher operators;
    - remainder-law corrections to sibling figures;
    - wall T_ann for both g* choices;
    - the GMR re-evaluation and the e needed for a heavy vector during inflation;
    - the Y-condensate extrapolation from gp_matching's [C] release (long_run_pi_over_4);
    - separate-universe domain sizes;
    - the truncated-EFT existence edge.
  - `code/companion_ad.py` [C]: one-field matched-EFT release including the radial companion beta6 abs(X)^6. It starts
    at the Hubble-held minimum with a supplied phase and runs a 282-run scan over phase and eta at r = 0.0707, with
    controls and r dependence. DOP853, rtol 1e-10.
  - `code/regular_yields.py` [C]: yield and eps against H_inf for c_H = 1, 10 and 100.
  - `code/full_minimum_vs_H.py` [C]: the full radial X, Y, Sigma minimum against H_inf, with psi2 aligned and at pi/2.
- **Web.** The arXiv abstract of 2101.12248 was read: Gmu ≲ 4e-15 for one loop model [Q]. The per-model body values
  were not re-read.

## Limits

- Homogeneous and classical throughout.
- The Y-condensate estimate extrapolates a [C] linear growth over two decades in a. It assumes a uniform-sign charge
  and full or no thermalization.
- GMR is a constant-mass scaling.
- The companion study uses the truncated EFT. Only the H ≤ 0.5-edge rows are within its validity
  (lambda_Y abs(Y)^2/M_Y^2 ≈ 0.1).
- The full-minimum scan is static. No full three-field release was run between the regular and chaotic regimes.
- A side remark: tracks/affleck_dine quotes eps ≈ 4.5 eta. In the definition it shares with gp_phase, and with
  A_6 = 12.97, my control gives eps ≈ 6.4 eta (q/2 at a late comoving energy of 3). This is not load-bearing here.

## Reproduction

See CONSOLIDATED.md section 8. Outputs: `results/critic_checks.json`, `companion_ad.json`, `regular_yields.json`,
`full_minimum_vs_H.json`.

Note on `results/full_minimum_vs_H.json`: in the c_H = 10, H = 4.7e13 GeV, psi2 = pi/2 row, abs(Sigma) did not converge
(2.65e8 GeV). There the Sigma potential is about 1e-11 of the total, below the optimizer tolerance. abs(X), the
quantity used, is unaffected.
