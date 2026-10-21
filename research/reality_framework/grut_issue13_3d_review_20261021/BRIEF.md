# Brief: Claude's audit of Codex's independent 3D charged-exit cross-check (a69bd86)

21 October 2026. Issue #13.

Codex published `a69bd861e76edeb93b2f5ce02b20d8186c505a20` on branch `codex/grut-expanding-charge-20261021`,
packet `research/reality_framework/grut_expanding_charge_20261021/` (C++ solver expanding.cpp, Gaussian wrapper,
run.py/gaussian_run.py/analyze.py/verify.py, RUN_TABLE.md, RESULTS.json, DERIVATION.md, REVIEW.md, ASSUMPTIONS.md,
raw_diagnostics.tar.gz). Claude's packet under review by Codex: `c7b8480`,
`research/reality_framework/grut_issue13_wall_history_20261020/` (REPORT.md and modules/*).

Codex's main results: three dense Gaussian-spectrum 3D preparations (two 192^3 seeds, one 288^3 larger box) form
both vacua (a in (1.54, 1.57]) and both genuine-vacuum phases wind all axes through a = 2.2; one of three sparse
64-mode seeds (s101) evolves to an almost single positive vacuum, robust to dt halving, 256^3/384^3 refinement and a
larger box; the expanding cold (sub-critical) control is wall-free; strict refinement convergence not met; two
large-box seeds unexecuted. Corrections to Claude: (1) 0.456 is rms speed, mean squared speed ~0.208, Q ~ 0.188
rho_w/t; (2) gamma(v_rms) = 1.124 is not an upper bound; energy-weighted E/(sigma A) lies in [(1-s)^-1/2, (1-s)^-1]
~ 1.124-1.263; no bound from an area-weighted mean; (3) core velocity proxies exceed one in charged states (estimator
failure); (4) instantaneous scalar w <= 1/3 is not a theorem even at small mean density; the cosmological argument
needs a justified coarse-grained state family; DeltaN_eff-equivalent is an energy comparison, not a Planck likelihood.

Rules for every audit agent:
- Get Codex's files READ-ONLY: `cd /home/user/TestingGRUT && git show a69bd86:<path>` or extract a copy with
  `git archive a69bd86 research/reality_framework/grut_expanding_charge_20261021 | tar -x -C <your scratch dir>`.
  Scratch dirs go under /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/.
  Never edit Codex's packet or Claude's earlier packets; no git commands other than read-only show/archive/log.
- Write ONLY inside your own directory research/reality_framework/grut_issue13_3d_review_20261021/audit/<key>/:
  code, results/*.json, NOTES.md (method, numbers, verdicts CONFIRMED / CORRECTED / REFUTED / NOT CHECKED, limits).
- numpy/scipy, g++ allowed. 4-CPU machine shared by up to 5 agents: use OMP_NUM_THREADS<=2 and keep each job
  <= 30 min wall; prefer reduced grids where a full rerun is too expensive, and say so.
- Label estimates [E] and preparation inputs [P]. No claims of exclusion, validation or identification of Candidate 1.
