# Sources and prior exposure

All work is retrospective. No outside scientist replicated these calculations.
We did not read every GRUT repository or rerun Claude's audit campaign in this
packet. The read scope relevant to the decision is listed below.

## Immutable repository inputs

Repository: https://github.com/ryangrvr/TestingGRUT

- Base: Codex `a69bd861e76edeb93b2f5ce02b20d8186c505a20`, especially
  `research/reality_framework/grut_expanding_charge_20261021/` protocol,
  derivation, review, run records and lossless diagnostic archive.
- Reviewed: Claude `a929f940f1c342ff597a5454481b7f4765d61ead`, packet
  `research/reality_framework/grut_issue13_3d_review_20261021/`:
  REVIEW.md; `audit/reach/NOTES.md`, `code/reach.py`, `code/qball_family.py`;
  `audit/v_reach/NOTES.md`; and `audit/wbound/NOTES.md`. Reviewed values remain
  quoted audit results unless explicitly recomputed here. This is not another
  run of its BVP/minimization, network, quasi-linear or cosmological codes.
- Original Candidate 1: `research/reality_framework/grut_candidate_1_20261018/`
  GRUT_CANDIDATE_1.md, WORKED_REGIME.md, parameter/calibration discussion,
  and `code/results.json`. The last stationary refinement supplies the prototype.
- `research/reality_framework/grut_3d_empirical_20261020/calibration_data.json`
  supplies the frozen BBN baryon calibration and mass convention. Its
  raw_diagnostics.tar.gz supplies the three old component tables, not a relic
  distribution. These bytes are inherited from the base repository.
- Raw source bytes and SHA-256 values for every file consumed by analyze.py are
  recorded in RESULTS.json, including both diagnostic archives. No full 3D field
  history or new large-box data is claimed in this packet.
- Prior audit reply:
  https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6322469291
- Our initial scope reply, before finishing this packet:
  https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6323817372
  The old WORKED_REGIME's existing m_D0=g_X v trajectory is explicitly acknowledged
  in REPORT §3; no correction to its equation is needed.

Immutable raw review link:
https://github.com/ryangrvr/TestingGRUT/blob/a929f940f1c342ff597a5454481b7f4765d61ead/research/reality_framework/grut_issue13_3d_review_20261021/audit/v_reach/NOTES.md

## Primary literature checked for this feasibility decision

1. Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*,
   https://arxiv.org/abs/1807.06209 (v4). Abstract checked: the quoted
   omega_c h²=0.120±0.001 is an inference in the stated cosmological framework.
   Used here as an already-known retrospective calibration, not a held-out test
   or a likelihood valid for arbitrary hidden stress.
2. Sagunski et al., *Velocity-dependent Self-interacting Dark Matter from Groups
   and Clusters of Galaxies*, https://arxiv.org/abs/2006.12515,
   JCAP 01 (2021) 024. Abstract checked: the group/cluster analyses use Jeans
   modeling, baryonic effects and stellar anisotropy. Their known inferred
   cross sections are model-dependent context, not direct universal area bounds.
   No raw observations, residual fit or new likelihood was analyzed here.
3. Colquhoun et al., *The Semi-Classical Regime for Dark Matter Self-Interactions*,
   https://arxiv.org/abs/2011.04679 (v2). Abstract checked: quantum versus classical
   scattering and velocity averaging matter to transport. Its Yukawa scattering
   formulas are **not** applied to the present composite objects. The arXiv HTML
   endpoint was unavailable; no claim of a complete full-text review is made.
4. Ballesteros, Garcia and Pierre, *How warm are non-thermal relics? Lyman-alpha
   bounds on out-of-equilibrium dark matter*, https://arxiv.org/abs/2011.13458.
   Abstract checked: nonthermal phase-space distributions and transfer functions
   enter the interpretation. No quoted thermal WDM mass bound or numerical
   length threshold is imported as a Candidate 1 exclusion.
5. Weinberg, *The cosmological constant problem*, Rev. Mod. Phys. 61, 1 (1989),
   https://doi.org/10.1103/RevModPhys.61.1. Publisher abstract checked as context
   for the subsequent direction decision. We do not claim a new proof or an
   exhaustive review of its adjustment theorem. The shift argument in
   NEXT_DIRECTION.md is derived directly, with its narrower domain explicit.
6. Kaloper and Padilla, *Sequestering the Standard Model Vacuum Energy*,
   https://arxiv.org/abs/1309.6562. Abstract checked as a prior proposal involving
   additional gravitational/global structure. No claim that it solves all aspects
   of the problem is adopted, no full review is claimed, and sequestering is not
   added to GRUT. It is a comparison obligation for a future proposed principle.

These sources justify the required kinds of matching, not a current exhaustive
sensitivity survey or a planned-detector guarantee. We use no recalled detector
forecast as a numerical test criterion. The 0.01/0.1/1 cm²/g contour labels in
RESULTS.json are **illustrative contact-proxy values**, not new observational
limits. The conditional exceptions in Claude's verifier are cited to establish
its own scope, not independently certified reach calculations.

## Constants and conventions

The density matching uses omega_c/omega_b times the frozen baryon-mass convention
directly. HBARC=1.973269804e-14 GeV cm and 1 GeV=1.782661922e-24 g are conventional
unit conversions at precision far beyond the physical assumptions used here.
The external reduced Planck mass is 2.435323e18 GeV. The illustrative spinodal
temperature uses s0=2891.2 cm^-3 and rho_c/h²=1.05375e-5 GeV cm^-3, as in the
quoted audit; it is not an independent new thermal calibration. Sensitivity to
g_* is displayed by specified constant choices, not fitted to an observable.

No full article, detector dataset or private Claude output is included.
