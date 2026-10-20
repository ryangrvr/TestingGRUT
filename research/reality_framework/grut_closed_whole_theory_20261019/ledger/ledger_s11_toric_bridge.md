# Campaign ledger entry — S11: Toric memory / scrambling bridge (key `s11_toric_bridge`)

Sector: the separate, unmerged branch `claude/toric4d-memory-vs-scrambling-20261019` (pre-registration commit `52b8f11`, results commit `8e0f89f`), which tested the "bridge claim" B that a passive self-correcting memory must be locally non-thermalising / non-scrambling, in the classical Metropolis dynamics of the 4D toric code.

Checkout note. The task names commit `96a7e7e` as the checkout; the working tree is at `4a53fb2` (WIP checkpoint 3), a descendant of `96a7e7e` (`git log --oneline -3`: 4a53fb2, f0e4cf0, da8f925; `git show --stat 96a7e7e` adds NEXT_CHAT_HANDOFF.md, START_NEXT_CHAT.txt, AGENT_COORDINATION.md). `git diff --stat 96a7e7e HEAD -- research/reality_framework/grut_closed_whole_handoff_20261019/` is empty, so the handoff text read here is byte-identical to the handoff commit's.

Merge status. `git branch -a --contains 8e0f89f` returns only `remotes/origin/claude/toric4d-memory-vs-scrambling-20261019`; merge-base with HEAD is `f28b49d`; the two branch-only commits are `8e0f89f` and `52b8f11`. The branch is not merged. Neither branch commit touches `AGENT_COORDINATION.md` (`git show --stat 52b8f11`: 7 files, all under the packet directory; `git show --stat 8e0f89f`: 19 files, all under the packet directory).

Execution note for this entry. No branch code was run. Two result JSONs were parsed with a stdin-only `python3 -I -c 'import json,sys; ...'` one-liner to verify two printed numbers (standard library `json` only; no file from `code/` imported or executed).

Path abbreviations used in citations:
- `tb:` = file on branch `origin/claude/toric4d-memory-vs-scrambling-20261019` at `8e0f89f`, under `research/reality_framework/toric4d_memory_vs_scrambling_20261019/` (read via `git show`). Line numbers are from `cat -n` of `git show`.
- `tb@52b8f11:` = the same path at the pre-registration commit.
- `hand/` = `research/reality_framework/grut_closed_whole_handoff_20261019/NEXT_CHAT_HANDOFF.md` (working tree, identical to 96a7e7e).
- `coord/` = `AGENT_COORDINATION.md` (working tree).
- JSON results are cited by key path or by `cat -n` line.

---

## 1. Question asked in this sector

The packet's own question was not the parent closed-whole question. It isolated three claims in "the last surviving sharp reading of GRUT's founding intuition": **M** (late-universe state is a passive self-correcting quantum memory whose rigidity is geometry), **S** (fundamental dynamics does not scramble, λ_L < 2πT), and the unstated **B (bridge)**: "a passive self-correcting memory is necessarily non-scrambling / non-thermalising in its local degrees of freedom, so M implies S" (`tb:REPORT.md:7-11`; `tb:README.md:7-13`). It asked whether B holds "in the one translation-invariant model proven to be a passive self-correcting memory, the 4D toric code", by "measuring the logical memory time and the local thermalisation time in the same thermal trajectories" (`tb:REPORT.md:13`; `tb:README.md:15`).

Under the parent question ("what does a closed system, whose every subsystem's bath is the rest of the universe, force beyond known physics?", `hand/` title/`START_NEXT_CHAT.txt:3`), the sector's actual role was a side test of one link in the "memory" route: if B fails, "memory" and "non-scrambling" are independent properties and "the unification that motivated the sharp reading is gone" (`tb:README.md:15`).

The pre-registration states explicitly: "This is a classical stochastic test. Its scope and its limits are stated in §5 before running." (`tb:README.md:17`).

## 2. Assumptions / model class

- **4D toric code, Z-error sector only.** Z₂ variables on the faces of an L⁴ torus (6L⁴ faces), one edge stabiliser per edge (6 faces per edge), energy E = 2·|∂e| for an error 2-chain e; logical readout on the X̄ membrane of plane (01) (`tb:REPORT.md:17`; `tb:README.md:21`; `tb:code/toric4d.py:13-36`). The thermal model is the 4D Z₂ 2-form gauge theory, dual to the 4D Ising model, β_c ≈ 0.954 via tanh β* = e^{−2K}, K_c ≈ 0.1497 (`tb:README.md:21`; `tb:REPORT.md:19`).
- **Dynamics: single-face Metropolis at fixed inverse temperature β**, "Local, Markovian, detailed-balance; a Glauber/Davies-type bath" (`tb:README.md:25`; `tb:code/toric4d.py:44-56`). This is an open system coupled to an external ideal heat bath at fixed β, not a closed whole whose bath is its own complement.
- **Decoder:** zero-temperature local cooling (ΔE<0 accepted, ΔE=0 with probability ½) for L² sweeps; "a local decoder, not minimum-weight" (`tb:README.md:27`; `tb:code/toric4d.py:82-95,115`).
- **Control:** 2D toric code on an L² torus, MWPM readout via pymatching, "No ordered phase (T_c = 0)" (`tb:README.md:23`; `tb:code/toric2d.py:11-27`).
- **Quantum illustration (not decisive):** exact 18-qubit 2D toric code on a 3×3 torus with field h = 0.3 on X and Z, infinite-temperature OTOC C(t) = 1 − Re⟨W(t)VW(t)V⟩ by stochastic trace; "2D is not self-correcting, so it does not test B" (`tb:README.md:29`; `tb:code/toric2d_quantum_otoc.py:1-5,46-70`).
- **Grid:** 4D L ∈ {3,4,5}, β ∈ {0.60, 0.90, 1.00, 1.10}; 2D L ∈ {4,6,8,12}; t_max = 10⁵ sweeps; 18-point log time grid [1, 2, 4, 8, 15, 30, 58, 115, 225, 444, 873, 1719, 3384, 6661, 13111, 25809, 50802, 100000] (`tb:README.md:42`; `tb:results/toric4d_L5.json` key `t_grid`; `tb:code/toric4d.py:192`). Trajectories as run: 4D 48/24/16 at L=3/4/5, 2D 64 (`tb:code/run_all.py:5`; `tb:REPORT.md:21`).
- **Pre-registered readings** R1 (method control), R2 (test of B, outcomes R2-a / R2-b), R3 (illustration), with the rule "No reading is changed after the numbers are seen" (`tb:README.md:44-54`).
- **Pre-registered scope exclusions:** no quantum OTOC or Lyapunov exponent in the 4D code ("infeasible here: 96 qubits at L=2"); the spin-2 sector question untouched; "Nothing here concerns de Sitter space; the dS chaos results are taken from the literature brief"; damage spreading secondary, "O2 (autocorrelation) carries R2" (`tb:README.md:58-61`).
- **External input not in the repository:** the "literature brief of 19 Oct 2026 (research report, this conversation)" that supplies the de Sitter λ_L results (Aalsma–Shiu 2020 λ_L = 2π/β; Kolchmeyer–Liu 2024 and Narovlansky 2025 4π/β) (`tb:README.md:11`; `tb:REPORT.md:13`). `git grep -i -E "Aalsma|Kolchmeyer|Narovlansky|literature brief"` on HEAD and on the toric branch outside the packet finds only an unrelated line (`research/reality_framework/grut_closure_fork_20261019/SOURCES.md:50`). See §12.

## 3. Executed evidence

### 3a. Pre-registration record and what changed before the reported runs

- Pre-registration commit `52b8f11` (2026-10-19 12:16:40 −0500): README, the three model scripts, runner, requirements, and `energy_scan_prelim.json`; "No memory results are included in this commit; the only data is the preliminary energy scan used to locate the transition (beta_c ~ 0.95 ...), which corrected the first draft's beta grid before any run" (`git show --stat 52b8f11`, message lines 1-9; 7 files, 618 insertions).
- Results commit `8e0f89f` (2026-10-19 12:58:05 −0500): adds REPORT.md, analyze.py, SUMMARY.md, gauge_diffusion_check.json, logs and nine result JSONs; modifies run_all.py, toric2d.py, toric4d.py, toric2d_quantum_otoc.py (19 files, 4576 insertions, 29 deletions) (`git show --stat 8e0f89f`).
- **β_c amendment (pre-run, recorded):** the first README draft gave the 1-form value β_c ≈ 0.44; the smoke-test energy scan showed the transition between β = 0.85 and 0.95 at L = 3, 4; the grid became {0.6, 0.9, 1.0, 1.1}; "the readings in §4 were not" changed (`tb:README.md:21`; `tb:REPORT.md:19`). The scan grid actually contains β ∈ {0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4}, with no β = 0.85 point; e.g. L=4 energy per face heating/cooling 0.0883/0.0747 at β=0.900 and 0.00706/0.03245 at β=0.950 (`tb:results/energy_scan_prelim.json:65-68,61-64`). "Between 0.85 and 0.95" is therefore a reading between the 0.9 and 0.95 grid points, not a measured 0.85 point.
- **Estimator change (recorded in REPORT, not in README):** `fit_tau` changed from an exponential least-squares fit (points 0.05<P<0.95, lower bound if P.min()>0.9) to the half-life t½ at which P_log first falls to 0.5, log-interpolated, with t_max as a lower bound when P never reaches 0.5 (`git diff 52b8f11 8e0f89f -- code/`: toric4d.py hunk @@ -172,17 +172,19 @@; toric2d.py hunk @@ -110,16 +110,19 @@; `tb:code/toric4d.py:174-187`). REPORT: "After the first production attempt the τ estimator was changed ... this is an estimator change, not a reading change, and was made before the reported runs" (`tb:REPORT.md:21`). The 8e0f89f commit message says the change is "recorded in README/REPORT", but `git diff 52b8f11 8e0f89f -- README.md` is empty; README still defines O1 as "τ_log from P_log = e^{−t/τ_log}" (`tb:README.md:35`). The first production attempt's outputs are not in the branch.
- **Trajectory-count change (not in README):** run_all.py changed from `NTRAJ = "16"` for every job to `NTRAJ4 = {"3": "48", "4": "24", "5": "16"}; NTRAJ2 = "64"` (`git diff 52b8f11 8e0f89f -- code/run_all.py`, line 5). README still says "16 trajectories per point" (`tb:README.md:42`); REPORT states 48/24/16 and 64 (`tb:REPORT.md:21`).
- **Quantum-illustration change (not recorded in README or REPORT):** `main(h=0.3, times=None, n_states=2, ...)` with 12 time points became `n_states=1` with 6 points [0, 0.5, 1, 2, 4, 8] (`git diff 52b8f11 8e0f89f -- code/toric2d_quantum_otoc.py`, hunk @@ -43,9 +43,9 @@). An earlier run on the longer grid is preserved in `log_2d.txt`, ending `t=  4.0 OTOC=0.9651 ... (927s)` followed by `Terminated` (`tb:results/log_2d.txt:17-24`). REPORT reports "one Haar state" without mentioning the reduction (`tb:REPORT.md:95`).
- **Mechanical thresholds introduced at analysis time:** `analyze.py` is a new file in `8e0f89f` and encodes R1/R2 numerically: t½ "O(1)" as `x < 50`; stabiliser L-independence as `max(edge) - min(edge) <= max(2.0, 0.5 * max(edge))`; "not frozen" as `dm[-1] != 1`; two consecutive lower bounds "treated as consistent" (`tb:code/analyze.py:111-125,137-146`). The README pre-registered these only qualitatively ("O(1)", "independent of L", "does not stay at its initial value") (`tb:README.md:46-50`). analyze.py's docstring "No reading is decided here that was not written before the run" (`tb:code/analyze.py:1-2`) is accurate for the qualitative readings, not for the numeric thresholds.

### 3b. O4 — transition location (4D energy per face at the end of the memory runs)

| L | β=0.60 | β=0.90 | β=1.00 | β=1.10 |
|---|---|---|---|---|
| 3 | 0.2889 | 0.0334 | 0.0043 | 0.0029 |
| 4 | 0.2941 | 0.0500 | 0.0030 | 0.0023 |
| 5 | 0.2961 | 0.0732 | 0.0057 | 0.0021 |

(`tb:REPORT.md:31-35`; identical in `tb:results/SUMMARY.md:5-9`.) The REPORT reads this as "The transition sits between β = 0.85 and 0.95 ..., consistent with the 4D-Ising duality value 0.954" (`tb:REPORT.md:29`).

### 3c. O1 — 4D logical half-life t½ (sweeps), cooling decoder; 48/24/16 trajectories

| L | β=0.60 | β=0.90 | β=1.00 | β=1.10 |
|---|---|---|---|---|
| 3 | 2.64 | 89 | 977 | 1.31e+04 |
| 4 | 2.21 | 180 | 1.64e+04 | > 10⁵ |
| 5 | 2 | 76.3 | 7.98e+04 | > 10⁵ |

(`tb:REPORT.md:39-43`; `tb:results/SUMMARY.md:15-19` prints "> 100000".) Raw: L=3 β=1.00 t_half = 977.3674002731318 (P_before 0.5417, P_after 0.2917); L=4 β=1.00 t_half = 16431.640985444083 (P_before 0.6667, P_after 0.1667); L=5 β=1.00 t_half = 79792.16887869872 (P_before 0.75, P_after 0.375); L=3 β=1.10 t_half = 13110.999999999989; L=4 β=1.10 and L=5 β=1.10: `'t_half': None, 'lower_bound_t_half': 100000.0, 'note': 'P_log > 0.5 at t_max'` (`tb:results/log_4d.txt:3,7,11,4,8,12`). At L=5 β=1.10 all 18 P_log entries are 1.0 with P_log_se 0.0 (`tb:results/toric4d_L5.json` keys `betas.1.10.P_log`, `betas.1.10.P_log_se`, verified by stdin parse; REPORT "1.00 ± 0.00 at t = 10⁵", `tb:REPORT.md:37`). With 16 trajectories at L=5 the P_log resolution is 1/16; L=4 β=1.10 ends at P_log 0.917 at t=10⁵ (`tb:results/log_4d.txt:8`).

### 3d. O5 — 2D control t½ (sweeps), MWPM, 64 trajectories

| L | β=0.60 | β=0.90 | β=1.00 | β=1.10 |
|---|---|---|---|---|
| 4 | 1 | 1.5 | 1.56 | 4 |
| 6 | 1 | 1.37 | 1.59 | 2.21 |
| 8 | 1 | 1 | 1.78 | 1.59 |
| 12 | 1 | 1.37 | 1.64 | 2.11 |

(`tb:REPORT.md:47-52`; `tb:results/SUMMARY.md:23-28`; raw in `tb:results/log_2d.txt:1-16`, e.g. 2D L=4 β=1.10 `'t_half': 4.0, 'P_before': 0.8125, 'P_after': 0.5`.) Values of 1 are "decayed by first readout" (`tb:results/log_2d.txt:1,5,9,10,13`).

### 3e. O2 — local thermalisation time τ_loc (1/e crossing), 4D, stabiliser / face spin

| L | β=0.60 | β=0.90 | β=1.00 | β=1.10 |
|---|---|---|---|---|
| 3 | 0.67 / 2.67 | 0.85 / > 30 | 0.90 / > 30 | 0.77 / > 30 |
| 4 | 0.66 / 2.51 | 0.82 / > 30 | 0.88 / > 30 | 1.18 / > 30 |
| 5 | 0.64 / 2.55 | 0.85 / 24.47 | 0.97 / > 30 | 0.95 / > 30 |

(`tb:REPORT.md:56-60`.) SUMMARY prints the same cells with columns in the order "site / stabiliser" and `inf` where REPORT prints "> 30", e.g. 4D L=5 β=0.90 "24 / 0.85" (`tb:results/SUMMARY.md:32-36`). The autocorrelation is computed from 400 observation sweeps at 1-sweep resolution over 64 faces and 64 edges after 400 equilibration sweeps, but only to lag 30 (`acf(S, 30), acf(Y, 30)`, `tb:code/toric4d.py:120-137`), so `inf` means no crossing within 30 lags; REPORT's "> 30" is the accurate rendering. Verified at L=4 β=1.10: acf_edge = [1.0, 0.411, 0.168, ...] (31 entries) giving τ ≈ 1.2; acf_face = [1.0, 0.949, 0.916, 0.902, 0.897, ...] with 0.764 at lag 30 (`tb:results/toric4d_L4.json` keys `betas.1.10.acf_edge`, `betas.1.10.acf_face`, stdin parse). 2D O2 rows (both observables ≤ 4.6 sweeps) are in `tb:results/SUMMARY.md:37-40` and not reproduced in REPORT.

**Which trajectories.** The O2 run uses seed `seed + 7` and the O3 run `seed + 11`, each a separate equilibration from the zero-error state, whereas memory trajectories use `seed + 1000*k + int(beta*100)` (`tb:code/toric4d.py:198,202-203`; same structure `tb:code/toric2d.py:136,139-140`). The measurements share the model, update rule and (L, β), but are not literally "in the same thermal trajectories" as `tb:REPORT.md:13` and `tb:README.md:62` phrase it.

### 3f. Gauge check (post-hoc) — face spin is gauge-variant

β = 1.1, L = 4, single trajectory from the zero-error state: rows [t, flipped-face fraction, energy per face] = [0, 0.0, 0.0], [10, 0.0, 0.0], [100, 0.010416666666666666, 0.0], [1000, 0.068359375, 0.005208333333333333], [3000, 0.1875, 0.0], [10000, 0.3489583333333333, 0.0], [30000, 0.515625, 0.0]; note: "flipped-face fraction grows while energy stays ~0: zero-energy stabiliser (cube-boundary) moves accumulate; face spins are gauge-variant" (`tb:results/gauge_diffusion_check.json:2-41`; REPORT table `tb:REPORT.md:62-72`). REPORT: "flipping the six faces of any cube is a stabiliser (zero energy), so a single face spin is defined only modulo stabilisers and drifts at zero cost. This was recognised only after the production runs" (`tb:REPORT.md:23`); "Its slow autocorrelation is gauge drift, not thermalisation" (`tb:REPORT.md:74`). The script that produced this JSON is not among the `code/` files (listing: analyze.py, requirements.txt, run_all.py, toric2d.py, toric2d_quantum_otoc.py, toric4d.py; `git ls-tree -r -l 8e0f89f -- .../code/`).

### 3g. O3 — damage spreading (secondary), 4D

Hamming distance D(t) at t = 1, 2, 4, 8, 16, 32, 64, 128 between twin copies driven by one random stream after one face flip: L=3 β=0.60 [35, 95, 177, 250, 243, 252, 259, 242] (486 faces); L=3 β=0.90 [3, 3, 3, 6, 6, 6, 6, 6]; L=3 β=1.00 and 1.10 all zeros; L=4 β=0.60 [13, 81, 443, 738, 766, 764, 759, 771] (1536); L=4 β=0.90 [10, 46, 137, 262, 430, 557, 716, 737]; L=4 β=1.00 [1, 0, 0, 0, 0, 0, 0, 0]; L=4 β=1.10 zeros; L=5 β=0.60 [26, 229, 1268, 1792, 1859, 1896, 1863, 1958] (3750); L=5 β=0.90 [4, 28, 99, 358, 779, 1377, 1754, 1841]; L=5 β=1.00 and 1.10 zeros (`tb:REPORT.md:78-91`; `tb:results/SUMMARY.md:48-59`). REPORT: "At β ≥ 1.0 the injected flip heals within one sweep at every size ... In 2D the behaviour is seed-dependent and not used" (`tb:REPORT.md:93`). The 2D rows that were set aside show nonzero plateaus, e.g. 2D L=4 β=0.60 [4, 4, 4, 4, 4, 4, 4, 4] and 2D L=12 β=1.00 [3, 5, 6, 16, 47, 72, 104, 123] (`tb:results/SUMMARY.md:60-75`). Note the twin update picks the same face and the same uniform u in both copies (`tb:code/toric4d.py:59-79`); at low T the single flipped face carries four violated edges and reverts deterministically when selected, so "heals within one sweep" is the expected behaviour of this update rule, consistent with the README's "update-rule dependent; secondary" label (`tb:README.md:37,61`).

### 3h. Q1 — exact 18-qubit perturbed 2D toric code (illustration only)

One Haar state, h = 0.3, times in 1/J: t=0: C=0.000, ⟨Z_loc(t)Z_loc⟩=1.000, ⟨Z̄(t)Z̄⟩=1.000; t=0.5: 0.068, 0.297, 0.890; t=1: 0.447, 0.115, 0.718; t=2: 0.764, 0.139, 0.385; t=4: 0.965, 0.028, 0.085; t=8: 0.981, 0.018, 0.019 (`tb:REPORT.md:97-104`; full precision `tb:results/toric2d_quantum_otoc.json:11-34`, e.g. otoc[t=8] = 0.980894781972593, `n_states: 1`; `tb:results/log_q.txt:1-6`, 477 s total). REPORT: "A perturbed stabiliser code scrambles; this is the expected quantum behaviour and does not bear on B, which concerns the 4D phase" (`tb:REPORT.md:106`). This is an infinite-temperature OTOC of an 18-qubit perturbed stabiliser Hamiltonian on a torus; it is not a thermal OTOC, has no Lyapunov-exponent extraction, no MSS comparison and no gravitational content (`tb:code/toric2d_quantum_otoc.py:46-70`).

### 3i. Readings as applied (mechanical, `analyze.py`)

- **R1 (method control): PASS** — "β=1.00: t½(4D) by L = L3: 977, L4: 1.64e+04, L5: 7.98e+04 → increasing; β=1.10: ... L3: 1.31e+04, L4: >100000, L5: >100000 → increasing; β=0.6 (disordered): t½(4D) = [2.64, 2.21, 2.0]"; 2D flat at 1–4 sweeps (`tb:REPORT.md:110-115`; `tb:results/SUMMARY.md:90-95`).
- **R2 by the letter of the pre-registration:** "R2-a on the stabiliser observable and on damage; FAILS on the face-spin observable" — face / stabiliser τ_loc: β=1.00 L3 inf / 0.9, L4 inf / 0.88, L5 inf / 0.97; β=1.10 L3 inf / 0.77, L4 inf / 1.2, L5 inf / 0.95 (`tb:REPORT.md:117-120`; `tb:results/SUMMARY.md:97-100`).
- **R2 with the face spin removed as gauge-variant (labelled POST-HOC):** "R2-a: the gauge-invariant local observables (stabilisers) thermalise in about one sweep, independent of L, and a local perturbation heals within one sweep, while the logical bit persists for 10^3 to >10^5 sweeps → B is FALSE in the 4D toric code" (`tb:REPORT.md:122-126`; `tb:results/SUMMARY.md:102-106`).
- **The packet's own fallback:** "A reviewer who rejects post-hoc observable corrections on principle should read O2 as INCONCLUSIVE and the packet as resting on O1 + O3 (memory grows with L while a local perturbation heals within one sweep), which already decides B against the bridge in this model" (`tb:REPORT.md:128`). O3 was pre-registered as secondary and "update-rule dependent" with "O2 (autocorrelation) carries R2" (`tb:README.md:37,61`), so this fallback rests the verdict on the observable the pre-registration ranked second.

### 3j. The reviewer's / report author's broader claims (recorded separately from what was computed)

The REPORT is authored by "Claude (claude.ai), owner-directed" with "Independent review outstanding" (`tb:REPORT.md:3`); there is no separate reviewer document on the branch. Its broader claims, none of which is computed in the packet:
- "A passive self-correcting memory does not keep local records and does not slow local thermalisation; it keeps exactly one global, topological bit. B is false here. Since M → S needed B, the sharp reading of the intuition reduces to two independent claims: M ..., which is unformulated and collides with the CLPW empty-de-Sitter result and the 3D no-go theorems; and S (non-scrambling), which stands alone and is contradicted at the horizon." (`tb:REPORT.md:132`).
- "The test came out as holography predicted. What remains of the fracton/memory route is the open question alone, with no mechanism connecting it to non-scrambling. The honest classification is the one pre-registered for outcome R2-a: retire the sharp reading. The program's result stands as the map: the responsive-vacuum and closed-whole intuitions are each consistent with known physics and add no law, no state selection, and, now, no memory–scrambling link beyond it." (`tb:REPORT.md:140`).
- "AdS/CFT (a code that is maximally chaotic) already made B implausible on general grounds. This packet makes the same point in a proven passive memory with numbers." (`tb:REPORT.md:136`).
The de Sitter λ_L values, the CLPW result, the 3D no-go theorems and the AdS/CFT argument are all imported from the unavailable literature brief (`tb:README.md:11,15,60`), not from anything in the repository.

### 3k. The packet's own statement of what is not established

"This is a classical stochastic (Davies-type) model. It measures local thermalisation and classical information spreading, not a quantum OTOC or a Lyapunov exponent in the 4D code (96 qubits at L = 2 is out of reach here). The open question 'can a passive self-correcting phase carry a Lorentz-invariant spin-2 sector' is untouched. Nothing here concerns de Sitter. The cooling decoder is local and not minimum-weight; its possible failure mode ... would *shorten* the measured memory time and so cannot manufacture the size dependence seen." (`tb:REPORT.md:134`).

## 4. Answer

What was compared: in the classical single-face Metropolis dynamics of the 4D toric code Z-sector at fixed β, the decoded logical half-life (O1) against the equilibrium autocorrelation time of local observables (O2) and twin-copy damage spreading (O3), at the same (L, β) but in separate runs with separate seeds (`tb:code/toric4d.py:198,202-203`); with a 2D toric-code MWPM control (O5) and an unrelated exact 18-qubit infinite-temperature OTOC illustration (Q1).

What the numbers show in that model: at β ≥ 1.0 the logical bit's half-life grows with L (977 → 1.64e+04 → 7.98e+04 sweeps at β=1.0; 1.31e+04 → >10⁵ → >10⁵ at β=1.1) while the edge-stabiliser autocorrelation crosses 1/e in 0.77–1.18 sweeps at every (L, β) and the injected single-face flip heals to D=0 within one sweep (`tb:REPORT.md:39-43,56-60,78-91`). The 2D control stays at 1–4 sweeps for all L (`tb:REPORT.md:47-52`). The face-spin autocorrelation does not cross 1/e within 30 lags at β ≥ 0.9 (except L=5 β=0.9: 24.47) and was reclassified post hoc as a gauge-variant observable, verified by a zero-energy flipped-fraction drift to 0.515625 at t=3×10⁴ (`tb:REPORT.md:58-60,62-74`; `tb:results/gauge_diffusion_check.json:35-39`).

The estimator change was from an exponential least-squares fit to a log-interpolated half-life, made after a first production attempt whose outputs are not preserved, recorded in REPORT but not in README (`tb:REPORT.md:21`; code diff; empty README diff). The post-hoc correction was the removal of the face spin from O2 as not an observable of the code, labelled as such in REPORT and SUMMARY (`tb:REPORT.md:23,122-128`).

The quantum illustration shows an 18-qubit perturbed 2D stabiliser Hamiltonian at infinite temperature scrambling (OTOC 0.981 by t=8/J) with the logical autocorrelation also decaying (0.019 at t=8), which the packet itself says "does not bear on B" (`tb:REPORT.md:97-106`).

Confirmation that this is not a closed-unitary gravitational MSS-exponent computation: the dynamics is an open Davies-type Metropolis bath at fixed β (`tb:README.md:25`); no Lyapunov exponent, no thermal OTOC, no 2πT comparison and no gravitational degree of freedom appears anywhere in the code or results (`tb:code/toric4d.py` entire; `tb:code/toric2d_quantum_otoc.py:46-70`); the packet states this itself (`tb:REPORT.md:134`; `tb:README.md:58`).

## 5. Closure verdict

Nothing additional beyond known physics was forced in this sector, and the sector did not operate in the closed-whole class at all. The computation is of an open system (a lattice code in contact with an external ideal bath at fixed β, `tb:README.md:25`), so it cannot speak to what a closed whole with internal baths forces. Within its own class it recovers known physics: the thermal stability of the 4D toric code (DKLP 2002; Alicki et al. 2010), the L-independent 2D memory time (Alicki–Fannes–Horodecki), and fast decorrelation of local gauge-invariant observables in a gapped low-temperature phase; the packet's own priority statement says "The only addition is the joint measurement of memory and local thermalisation as a test of B" (`tb:REPORT.md:136`; `tb:README.md:62`). The computation that shows "nothing additional" is the O1/O2/O3 table set in §3c-3g: a proven passive memory with size-growing logical half-life and ~1-sweep local thermalisation, i.e. a counterexample to bridge claim B in the classical 4D toric code. What this removes from the GRUT program is a proposed internal link (memory ⇒ non-scrambling), not a law of physics; no new constraint, state-selection rule, chaos bound or gravitational statement is produced.

## 6. Classification

**mixed**:
- **recovered known physics** — thermal stability / size-growing memory of the 4D toric code, L-independent 2D memory, fast local thermalisation of stabilisers; acknowledged as prior literature (`tb:REPORT.md:136`; `tb:README.md:62`).
- **rejected mechanism (conditional)** — bridge claim B is rejected in the stated class: classical Metropolis dynamics, Z-sector 4D toric code, β ∈ {1.0, 1.1}, L ∈ {3,4,5}, local cooling decoder, t_max = 10⁵ sweeps, with the O2 verdict depending on a labelled post-hoc observable correction and the fallback resting on secondary O3 (`tb:REPORT.md:122-128`).
- **supplied model choice** — Metropolis bath at fixed β, cooling decoder, grid, seeds, half-life estimator, numeric thresholds in analyze.py (`tb:README.md:25,27,42`; `tb:code/analyze.py:111-125,137-146`).
- **conjecture** — the broader claims about de Sitter λ_L, CLPW, 3D no-go theorems, "as holography predicted", and "retire the sharp reading" rest on an out-of-repository literature brief and are not computed (`tb:REPORT.md:132,136,140`).
- **unexecuted task** — quantum OTOC / Lyapunov exponent in the 4D code; spin-2 sector question; independent review (`tb:REPORT.md:3,134`; `tb:README.md:58-59`).

## 7. Remaining freedom

- Any quantum (unitary or Lindbladian) version of the 4D toric code with a perturbation: no OTOC or λ_L was computed; R2-b's escalation path ("a quantum OTOC computation in a large-N self-correcting model") was never triggered or run (`tb:README.md:50,58`).
- The full code (X- and Z-sectors together), non-Metropolis baths, non-Markovian baths, and any closed-system (internal-bath) variant: all unexamined; the test used one sector and one Davies-type rule (`tb:README.md:21,25`).
- Temperatures other than β ∈ {0.6, 0.9, 1.0, 1.1}, sizes L > 5, times beyond 10⁵ sweeps, and more than 16 trajectories at L=5 (`tb:REPORT.md:21,37`).
- Decoder choice: the local cooling decoder is not minimum-weight; its failure mode is argued only to shorten memory time (`tb:REPORT.md:134`).
- Whether a passive self-correcting phase can carry a Lorentz-invariant spin-2 sector (claims M) — untouched (`tb:REPORT.md:134`).
- Everything about de Sitter, the MSS bound, and gravitational chaos — outside the packet (`tb:README.md:60`).

## 8. Failures retained

1. The pre-registered O2 face-spin observable was ill-posed ("not an observable of the code"), recognised only after the production runs; by the letter of the pre-registration R2 "FAILS on the face-spin observable" (`tb:REPORT.md:23,117-120`).
2. The first README draft's β_c ≈ 0.44 was wrong (1-form value); corrected by the energy scan before memory runs (`tb:README.md:21`; `tb:REPORT.md:19`).
3. A "first production attempt" with the exponential-fit estimator was superseded; its outputs are not preserved in the branch (`tb:REPORT.md:21`; branch tree listing).
4. The 8e0f89f commit message states the estimator change is "recorded in README/REPORT", but the README was not changed between 52b8f11 and 8e0f89f (empty `git diff`), so README §3 still describes the superseded exponential fit and "16 trajectories per point" (`tb:README.md:35,42`).
5. The quantum illustration was reduced from 2 Haar states × 12 times to 1 state × 6 times; an earlier run was `Terminated` at t=4.0 after 927 s (`tb:results/log_2d.txt:17-24`; code diff); not mentioned in REPORT.
6. 2D damage-spreading results were set aside as "seed-dependent and not used" (`tb:REPORT.md:93`) although present in SUMMARY (`tb:results/SUMMARY.md:60-75`).
7. "Independent review outstanding" (`tb:REPORT.md:3`); the branch is unmerged and no coordination entry exists for it.
8. REPORT/README phrase the design as "in the same thermal trajectories" (`tb:REPORT.md:13`; `tb:README.md:62`), whereas the code measures O2 and O3 in separate seeded runs at the same (L, β) (`tb:code/toric4d.py:198,202-203`).
9. The numerical thresholds that make R1/R2 "mechanical" were written in analyze.py at results time, not in the pre-registration (`tb:code/analyze.py` new file in 8e0f89f; `tb:README.md:46-50`).
10. L=4 and L=5 at β=1.1 are lower bounds only ("> 10⁵"), with 24 and 16 trajectories (`tb:REPORT.md:37,42-43`).

## 9. What must not be promoted

- "B is FALSE in the 4D toric code" must not become "memory and non-scrambling are independent in general" or "GRUT's memory route is refuted"; the result holds in one classical stochastic sector of one model at two temperatures and three sizes (`tb:REPORT.md:122,134`).
- The result must not be described as a computation of a Lyapunov exponent, an MSS-bound test, a de Sitter result, or a closed-system/internal-bath result; it is none of these (`tb:REPORT.md:134`; `tb:README.md:25,58,60`).
- The 18-qubit OTOC (0.981 at t=8/J) must not be cited as evidence about the 4D memory phase or about B (`tb:REPORT.md:106`; `tb:README.md:29,52`).
- "> 10⁵ sweeps" at L ≥ 4, β=1.1 is a lower bound from 24/16 trajectories, not a measured memory time (`tb:REPORT.md:37`).
- The de Sitter λ_L values (2π/β, 4π/β), the CLPW collision and the 3D no-go theorems are imported from an unavailable brief; they must not be entered into the ledger as repository-verified facts (`tb:README.md:11`; `tb:REPORT.md:132`).
- "The test came out as holography predicted" and "retire the sharp reading" are the report author's interpretation; the handoff instructs that the computed content be kept distinct from them (`tb:REPORT.md:140`; `hand/:194-195`).
- The post-hoc gauge correction must stay labelled post hoc; the by-the-letter reading fails on O2's face spin (`tb:REPORT.md:117-128`).
- O3 damage healing is update-rule dependent and pre-registered as secondary; it must not be promoted to a primary non-scrambling measure (`tb:README.md:37,61`).

## 10. Sources read

- `tb:REPORT.md` (all 144 lines), `tb:README.md` (66 lines), `tb:results/SUMMARY.md` (106 lines), `tb:results/log_q.txt` (6 lines), `tb:results/gauge_diffusion_check.json` (42 lines), `tb:results/log_4d.txt` (12 lines), `tb:results/log_2d.txt` (24 lines), `tb:results/energy_scan_prelim.json` (78 lines), `tb:results/toric2d_quantum_otoc.json` (36 lines).
- `tb:code/toric4d.py` (219 lines), `tb:code/toric2d_quantum_otoc.py` (79 lines), `tb:code/toric2d.py` (lines 1-60, 100-156), `tb:code/run_all.py` (21 lines), `tb:code/requirements.txt`, `tb:code/analyze.py` (via the 8e0f89f diff, 159 lines). Read only; not executed.
- `tb:results/toric4d_L5.json` and `tb:results/toric4d_L4.json`: keys `t_grid`, `n_faces`, `betas.1.10.{P_log,P_log_se,acf_edge,acf_face,tau_log}` via stdin-only JSON parse.
- `git show --stat 52b8f11`; `git show --stat 8e0f89f`; `git diff 52b8f11 8e0f89f -- README.md` (empty) and `-- code/`; `git ls-tree -r --name-only 8e0f89f -- <packet>`; `git branch -a --contains 8e0f89f`; `git merge-base`.
- `hand/:100-135` (campaign table incl. row at line 124), `hand/:150-170`, `hand/:180-205` (toric branch note at 190-195); `START_NEXT_CHAT.txt:1-9`; `coord/:5800-5832` (entries mentioning the toric read, line 5828); `git grep` on HEAD and the toric branch for toric/scrambling/MSS/Lyapunov/literature-brief terms.

## 11. Handoff-row discrepancies

Row (`hand/:124`): "Separate report on another branch; read before importing. It compares thermal Metropolis trajectories and decoded memory, with a post-hoc gauge-observable correction and a separate quantum illustration. It is not a closed-unitary gravitational MSS-exponent computation." Supporting note (`hand/:190-195`): branch and commit named; "The first 105 lines of its REPORT were read ... No toric code or trajectories were rerun here. Preserve that distinction from the reviewer's broader claims."

1. **Partial read acknowledged by the handoff; the row omits the packet's own verdict.** REPORT.md is 144 lines; lines 108-144 hold the applied readings, the "established / not established" section, the program consequence and reproduction. The row therefore does not record that R1 is PASS, that R2 by the letter FAILS on the face spin, that the post-hoc reading is R2-a ("B is FALSE in the 4D toric code"), or that the author's recommendation is "retire the sharp reading" (`tb:REPORT.md:110-128,140`). Checked against primary sources: the row's description is accurate as far as it goes; it is incomplete rather than wrong.
2. **"Independent review outstanding" is not carried into the row** (`tb:REPORT.md:3`); the row's "read before importing" is the only caution.
3. **"The reviewer's broader claims"**: there is no reviewer document on the branch; the broader claims are the REPORT author's (Claude, owner-directed, `tb:REPORT.md:3`). The handoff's label conflates author and reviewer.
4. **The estimator change is recorded only in REPORT**, not in README as the 8e0f89f commit message claims; the handoff mentions "its method, estimator change and gauge-observable correction" (`hand/:193-194`) without noting the README/REPORT mismatch (empty `git diff 52b8f11 8e0f89f -- README.md`; `tb:README.md:35,42`).
5. **"compares thermal Metropolis trajectories and decoded memory"** is accurate; the primary sources add that O2/O3 are measured in separate seeded runs at the same (L, β), not in the memory trajectories themselves (`tb:code/toric4d.py:198,202-203`), contrary to the REPORT's own "in the same thermal trajectories" wording (`tb:REPORT.md:13`).
6. **The row does not note that the model is an open Davies-type bath** (`tb:README.md:25`), i.e. outside the closed-whole class that the parent question names; "not a closed-unitary gravitational MSS-exponent computation" is correct but understates it.
7. No discrepancy on the branch name and commit hash: `claude/toric4d-memory-vs-scrambling-20261019`, `8e0f89f5d675f13fe8aa62e5ee17196b406dca12` (`hand/:191-192`; `git show --stat 8e0f89f`).

## 12. Unread or missing

- The "literature brief of 19 Oct 2026 (research report, this conversation)" cited by `tb:README.md:11` and `tb:REPORT.md:13` is not in the repository on HEAD or on the toric branch (`git grep` for Aalsma/Kolchmeyer/Narovlansky/"literature brief" finds only `research/reality_framework/grut_closure_fork_20261019/SOURCES.md:50`, an unrelated line). Every de Sitter, CLPW, 3D no-go and AdS/CFT statement in the packet depends on it and could not be checked.
- The "first production attempt" with the exponential-fit estimator (`tb:REPORT.md:21`): no outputs preserved.
- The script that produced `results/gauge_diffusion_check.json`: not among the `code/` files.
- The quantum run on the 12-point grid with 2 Haar states: only the terminated partial log survives (`tb:results/log_2d.txt:17-24`).
- `tb:code/toric2d.py` lines 61-99 (2D memory trajectory and autocorrelation bodies) were not read in full; `tb:results/toric2d_L*.json`, `toric4d_L3.json` and the non-verified keys of `toric4d_L4.json`/`toric4d_L5.json` were not read beyond SUMMARY/log renderings.
- No independent review document exists for this branch; none was read.
- Primary literature named by the packet (DKLP 2002; Alicki et al. 2010; Alicki–Fannes–Horodecki; Derrida–Weisbuch 1987; Aalsma–Shiu 2020; Kolchmeyer–Liu 2024; Narovlansky 2025) was not consulted here; attributions are reported as the packet prints them.
