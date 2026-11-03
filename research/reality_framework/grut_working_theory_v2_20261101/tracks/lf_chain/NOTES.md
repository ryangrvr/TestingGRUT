# lf_chain: the M2 windows on the computed (LF) reheating closure

3 November 2026 · Claude, lattice-fermion closure pass (BRIEF_LF.md is the pre-registration; GitHub issue #13 context).
Same-session AI work on one supplied model. Every history, splice and hand-over rule, closure, holder amplitude and
phase, relic rule, interpolation rule, window definition and criterion is a model choice; none follows from closure.
Nothing here is a validation, an exclusion, a new theory or a GRUT prediction. Labels: [D] derived, [C] computed here,
[Q] quoted from a reviewed track, [E] estimate (interpolation/extrapolation), [A] argued or assumed, [R] literature.

**Inputs read:** BRIEF_LF.md, BRIEF_E2E.md (M2 definitions), e2e_chain NOTES + SKEPTIC, pr_principle NOTES +
PREREGISTRATION + code, lf_fermion NOTES (§0, §4, §5) + the lf_fermion skeptic's digest and corrections.
**Reuse, read-only, unchanged:** e2e_chain `handover.py` (run as `__main__` via runpy), `splice.build`, `scan.py` (run as
`__main__`), `chain.py`, `fixpoint.solve/node`; through them mf_ideal's release engine, e2e_kinetic's network and
e2e_relic's `pmap`/`relic`; pr_principle's `edge.py` and `phase_soft.py` (byte-identical copies in `code/`, run from here
only so that they write into this track). Redirections only: scratch (`LFC_SCRATCH`) and `E.RES` → this track;
`splice.lattice_block` reads lf_fermion's `hist_HFa_LF.json` for the `KF*` keys (column logic re-typed unchanged).
Python ran with `-B`; no binaries or `__pycache__` in the tree. Scratch: `.../scratchpad/wt2/lf_chain/`.

---------------------------------------------------------------------------------------------------------------------

## 0. Result in brief

1. **K on the computed closure [C]: 5.53 at 1.3×10¹¹ GeV (brackets KFlo 5.85 / KFhi 5.01), 6.25 at 8.4×10¹⁰.** Old:
   P0 8.77 / 9.21, GM 15.28 / 16.42. The lf_fermion [E] estimate (5.4; skeptic 5.06–5.37) was 2–9% low.
   T_RH = 3.2×10¹² GeV (design node).
2. **M2-Q exists and widens.** Self-consistent Q-matter for c0 ≈ 1.377 → ≥ 1.60 at δ = 0 (upper edge not reached),
   δ ≈ −4…0° at c0 = 1.38–1.40 and ≥ −4…+4° at 1.45 (free X at ±8°). At the old holder (1.45, 0):
   H_i = 6.9×10¹⁰ GeV (old P0 9.9×10¹⁰), K 6.6, S/Q(10⁻⁶) = 1.98 (2.43), e_crit = 2.5×10⁻⁷ (4.3×10⁻⁷),
   r = 7.7×10⁻⁸ (1.6×10⁻⁷), D = 7.7 (7.9), s_* ≥ 9.2×10¹⁵ GeV (1.36×10¹⁶) [E: root at w = 0.52 between the
   5.5×10¹⁰ [A] node and 8.4×10¹⁰; −10% with linear-in-H interpolation].
3. **M2-X (free X, e_crit ≈ 2×10⁻⁶) survives only as a thin strip on the M2-Q boundary.** At δ = 0: c0 ≈ 1.36–1.377
   (e_crit 1.3×10⁻⁵ → 1.3×10⁻⁶). Representative point (1.37, 0): H_i = 7.9×10¹⁰ GeV, K 6.4, S/Q 3.82 / 3.30,
   **e_crit = 1.9×10⁻⁶**, M = 1.0002, r = 1.0×10⁻⁷, D ≈ 9.9, s_* ≥ 1.36×10¹⁶ GeV. Old GM M2-X (1.55, 0) is Q-matter
   on LF. Other proper-strip points: (1.40, +4°) e_crit 2.4×10⁻⁶, (1.40, −8°) 3.4×10⁻⁶ [E].
4. **Movement:** the windows move down in c0 (M2-X 1.55 → 1.37; M2-Q's lower edge 1.40 → 1.377) and down in H_i
   by 20–30% (r roughly halves). Branch assignments at (1.37, 1.40, 1.45; δ = 0) are unchanged across the KFlo/KFhi
   brackets; roots move by ±5–9%.
5. **Bistability (Q4-b):** one band, c0 = 1.3766–1.3778 at δ = 0 [E], where both a free-X and a Q-matter root exist.
6. **A second, low-c0 free-X pocket with an e_crit** sits at c0 ≈ 1.08–1.16 (δ = −4…0°), on the N_YL interference
   minimum (N_YL 1.5–2.6), with self-consistent H_i ≈ 2.1–2.4×10¹¹ GeV, which lies **outside** the computed lattice
   range (only an [A] node at 2.7×10¹¹ that reuses the 1.3×10¹¹ lattice block brackets it). e_crit ≤ 4×10⁻⁶ only
   for c0 ≈ 1.086–1.109 at δ = −4° [E/A].
7. **H1 (pre-registered): FAILS** under the pre-registered window definition (M2-X = e_crit ≈ 2×10⁻⁶, operationalised
   as 1–4×10⁻⁶; or M2-Q). The c0_soft points (1.08, −4°), (1.13, 0), (1.17, +4°) are self-consistent free X with
   e_crit 5.3×10⁻⁶, 9.6×10⁻⁶ and none. Closest miss: δ = −4°, Δc0 ≈ 0.006 below the proper pocket's edge. Under a
   looser "any e_crit" reading, 2 of 3 would count, but only at H_i ≈ 2.3–2.5×10¹¹, off the computed lattice range.
   Full statement in §5.
8. **Q3:** B (max CP torque, δ = −11…−14°) lands on free X without e_crit: miss. D (max ε): ε rises monotonically
   (0.09 → 0.83 over c0 1.0–1.6): no value. E (edge of stability): no tachyon for c0 1–2.6: no value. F: the S/Q
   minimum (1.87 at (1.45, 0)) is inside M2-Q, branch-defining by construction. **No non-tuned hit.**
9. **Q4 (lf_fermion §0):** Q4-a (e = e_crit) is the M2-X/M2-Q boundary line (c0 ≈ 1.377 at δ = 0), branch-defining.
   Q4-b (edge of bistability) gives c0 = 1.3766 / 1.3778 (δ = 0), branch-defining. Q4-c (t_release = t_coh,end): the
   release takes 1.42–4.7 H_i⁻¹ (90% of |N_X|), but the LF coherent stage ends at 0.85 H_i⁻¹; no c0 in 1.0–1.6
   satisfies it, so it gives no value. **No principle becomes non-trivially testable.**

---------------------------------------------------------------------------------------------------------------------

## 1. Final histories on LF (`code/build_hist.py`; `results/splice_checks.json`, `results/handover_*.json`) [C]

- e2e_chain `handover.py`, unchanged: lattice → kinetic moments at t_sw = 3; the degenerate-sphere fix fires for all
  four species (n̄ needed 1.2–1.5 at the design node); the dephased network starts at t_sw. The N split uses the kinetic
  reference's own species ratios (handover.py default). lf_fermion's per-species export was **not** used for the split; that variant was not run.
- `splice.build`, unchanged. Lattice-stage m_s² is the LF column (running mean, as e2e_chain), as lf_fermion requires.
- e2e_chain skeptic correction 8 applied: kinetic-block m_s² = 9(1 − 2x), x = 5×10⁻⁴ (step at t_sw 9.075 → 8.991).

| history | lattice block | H (GeV) | K | H jump at t_sw | lattice budget |
|---|---|---|---|---|---|
| hLF_mid | KFmid (computed) | 1.3e11 | **5.528** | 9.6e-4 | 6.2e-3 |
| hLF_lo / hLF_hi | KFlo_nl / KFhi_nl | 1.3e11 | 5.852 / 5.014 | 1.1e-3 / 9.7e-4 | 6.2e-3 |
| hLF84_mid | KFmid_H84_nl (computed) | 8.4e10 | **6.254** | 5.2e-4 | 9.7e-3 |
| hLF270_mid | KFmid, nearest-H rule [A] | 2.7e11 | 4.553 | 1.9e-3 | – |
| hLF55_mid | KFmid_H84_nl, nearest-H rule [A] | 5.5e10 | 7.017 | 3.2e-4 | – |

The two extra nodes follow e2e_chain's stated rescaling rule (lattice block nearest in H, kinetic network at that H).
lf_fermion found the early stage H-dependent (8.4×10¹⁰ does not fragment), so they are [A].

## 2. Scans and the fixpoint rule (`code/scan_run.py`, `code/fp_lf.py`; `results/scan_*.json`, `fixpoints_LF.json`)

- e2e_chain `scan.py`, unchanged: c0 = 1 base through mf_ideal's bench_run, then `chain.release_c0` + `chain.full`
  (corrected relic, c_L = 0, grav). 45 (c0, δ) keys on the two computed nodes, a subset on the [A] nodes, and
  5 keys on KFlo/KFhi; 150 engine runs in all, bases included (13–24 s each).
- **Branch-aware fixpoint** (e2e_chain skeptic, sk_static D, re-typed unchanged): f_b = ln[(H_req/H_i)(M_b/M_node)^(−2/3)],
  interpolated linearly in ln H. A root counts only if the interpolated S/Q(10⁻⁶) is on branch b. "primary" uses the
  computed pair (8.4×10¹⁰, 1.3×10¹¹), with extrapolation flagged by w. "all" is piecewise over every node of the key.
  The e2e_chain rule is recorded for comparison. The KFlo/KFhi brackets use the KFmid slopes [E], as the skeptic did
  for GMg03.
- **Window definitions** [model choice]:
  - M2-Q: a self-consistent Q-matter root at e = 10⁻⁶.
  - M2-X proper: a self-consistent free-X root with 1×10⁻⁶ < e_crit ≤ 4×10⁻⁶. This follows BRIEF_E2E's "e_crit ≈ 2×10⁻⁶";
    the e2e_chain skeptic gave the old window's actual range as 1.3–3.8×10⁻⁶.
  - M2-X loose: any e_crit ≤ 10⁻³ (e2e_chain's operational wording).
  - **Disclosure:** the 4×10⁻⁶ cutoff was written into `fp_lf.py` after the first design-node scan was seen. That scan
    included e_crit = 6.3×10⁻⁶ and 1.5×10⁻⁵ at two c0_soft points. The cutoff was taken from the skeptic's range, not
    tuned, but it was not blind.

## 3. Q2: new points and windows [C at nodes, E between, A where an [A] node is used]

| | M2-X proper (1.37, 0) | M2-Q (1.40, 0) | M2-Q (1.45, 0) | old GM M2-X (1.55, 0) | old P0 M2-Q (1.45, 0) |
|---|---|---|---|---|---|
| H_i (GeV); w; pair | 7.91e10; 0.86; 5.5–8.4e10 | 7.98e10; 0.88; 5.5–8.4e10 | 6.86e10; 0.52; 5.5–8.4e10 | 1.123e11 [Q] | 9.89e10 [Q] |
| primary rule (computed pair) | 7.89e10 (w −0.14) | 7.96e10 (w −0.12) | 6.81e10 (w −0.48) | | |
| linear in H | 7.74e10 | 7.83e10 | 6.19e10 | | |
| K | 6.36 | 6.35 | 6.62 | 15.7 | 9.05 |
| S/Q 10⁻⁶ / 3×10⁻⁶ | 3.82 / 3.30 | 2.64 / 2.12 | 1.98 / 1.45 | 4.03 / 3.05 | 2.43 / 1.75 |
| e_crit | **1.9e-6** | 4.2e-7 | 2.5e-7 | 1.72e-6 | 0.43e-6 |
| M | 1.0002 (free X) | 0.804 | 0.797 | 1.0003 | 0.802 |
| ε; N_YL | 0.50; 43 | 0.55; 15 | 0.65; 2.5 | 0.75; 35 | 0.61; 6.4 |
| r (cm.r_of) | 1.02e-7 | 1.04e-7 | 7.7e-8 | 2.05e-7 | 1.59e-7 |
| D = −(3/2) dln H_req/dln c0 (8.4e10 node) | 9.9 (one branch side) | 8.3 | 7.7 | 5.77 | 7.87 |
| s_* floor (β_iso < 0.038) | 1.36e16 GeV | 1.15e16 | 9.2e15 | 1.13e16 | 1.36e16 |

The table quotes μ_eff from the code convention u = μ + c0 e^{iδ}(10i − μ); it is in `windows.json`. N_* was not
re-solved here: cmb.py was not run.

**Windows on the grid** (`windows.json` → `grid`, per key):
- **δ = 0, by c0:**
  - 1.05–1.30: free X with no e_crit (M = 1.1–2.9). The exceptions are 1.13–1.16 (loose M2-X, see below).
  - 1.35–1.374: free X with an e_crit (1.3×10⁻⁴ → 1.3×10⁻⁶). M2-X proper for c0 ≳ 1.365.
  - 1.3766–1.3778: bistable.
  - 1.378–1.60: M2-Q.
- **δ scans:**
  - c0 = 1.38: −4° Q, +4° free X (e_crit at one node).
  - c0 = 1.40: −8° M2-X proper, −4/0° Q, +4° M2-X proper, +8° free X without e_crit.
  - c0 = 1.45: ±8° loose M2-X (root e_crit ≈ 4×10⁻⁵ at −8°, root extrapolated below 5.5×10¹⁰; 6.8×10⁻⁶ at +8°), −4…+4° Q.
- **Low-c0 pocket** (roots at 2.1–2.4×10¹¹ through the [A] node):
  - δ = −4°: 1.08 (e_crit 5.3×10⁻⁶), 1.09 (3.4×10⁻⁶), 1.10 (3.0×10⁻⁶), 1.12 (5.8×10⁻⁶). The proper pocket is 1.086–1.109.
  - δ = −2°: 1.10 (loose).
  - δ = 0: 1.13 (9.6×10⁻⁶), 1.14 (5.2×10⁻⁶), 1.16 (7.0×10⁻⁶); all loose.
  - δ = +4°: 1.15 / 1.17 have no e_crit; 1.19 has an e_crit at one node only.
- **Brackets:** at (1.37, 1.40, 1.45; δ = 0) the branch assignment is the same on KFlo and KFhi. The roots move as
  follows (free X; Q; Q):
  - KFlo: 8.3 / 8.4 / 7.2×10¹⁰;
  - KFhi: 7.2 / 7.3 / 6.3×10¹⁰.

  The lf_fermion skeptic's corrections (N1 κ_frag lower, N2's early Q lower) point toward KFlo (K +6%).
- **Systematics:**
  - Interpolation rule: −2% (M2-X) to −10% (extrapolated M2-Q) in H_i.
  - Every M2 root at δ = 0 except (1.378, 0) lies below 8.4×10¹⁰ and needs the 5.5×10¹⁰ [A] node or extrapolation.
  - Phase steps are 2–4°. The e2e_chain skeptic showed that 1° steps can open pockets (N5), so the δ edges are coarse.

**Do they move? Yes.**
- M2-Q is no longer marginal anywhere on the grid. Its lower edge moves 1.40 → 1.377, its upper edge is beyond 1.60, and
  its H_i falls by about 30%.
- M2-X proper moves from c0 1.55 to the Q boundary at about 1.37. Its H_i is −30% and its r is halved.
- The old GM M2-X holder is now Q-matter.

## 4. Q3: candidates B, D, E, F (same code, same definitions; `results/edge.json`, `phase_soft.json`)

- `edge.py` and `phase_soft.py` re-run unchanged. They are inflation-stage quantities and reproduce pr_principle's
  numbers: δ_B = −11.03 / −11.93 / −12.77 / −13.54 / −14.26° at c0 = 1.40–1.60; softest held mode 8.747 / 8.749 /
  8.748 H² at c0 = 1.08 / 1.13 / 1.17 (δ = −4 / 0 / +4°); no negative held eigenvalue for c0 1–2.6.
- **B (maximal CP torque).** On the closure, (1.37, −11°), (1.40, −11°), (1.45, −12°), (1.50, −13°) are self-consistent
  free X with S/Q 4.8–5.8 and no e_crit. (1.55, −14°) is loose (e_crit 0.7–2×10⁻⁵, root below 5.5×10¹⁰). It misses M2-Q
  by ≥ 7° and misses M2-X proper. **Miss.**
- **D (maximal yield).** ε at δ = 0 (1.3×10¹¹) rises monotonically: 0.092, 0.127 (1.10), 0.252 (1.25), 0.455 (1.35),
  0.554 (1.40), 0.650 (1.45), 0.789 (1.55), 0.832 (1.60). There is no interior maximum. **No value.**
- **E (edge of stability).** No edge. **No value.**
- **F (extremal S/Q).** The maximum is 25.8 at c0 = 1 (free X, no e_crit, outside the windows). The minimum is 1.87 at
  (1.45, 0), inside M2-Q, and **branch-defining by construction**.
- Chance baseline unchanged (pr_principle): about 0.15–0.25 expected chance hits. Observed non-tuned hits: 0.

## 5. H1, answered exactly as pre-registered

H1: the softest held mode sits at c0_soft = 1.08 / 1.13 / 1.17 for δ = −4 / 0 / +4°. It **succeeds** if c0_soft at the
window's own δ lies inside a self-consistent window, M2-X or M2-Q, on the computed closure; it **fails** otherwise.

| δ | c0_soft | self-consistent root (all-node rule) | e_crit | M2-Q? | M2-X proper? | loose? |
|---|---|---|---|---|---|---|
| −4° | 1.08 | free X, H_i 2.38e11 (1.3e11–2.7e11 [A] pair; primary: 2.52e11, w 2.5) | 5.3e-6 | no | **no** (proper pocket starts at 1.086) | yes |
| 0 | 1.13 | free X, 2.35e11 [A] (primary 2.52e11, w 2.5) | 9.6e-6 | no | no | yes |
| +4° | 1.17 | free X, 2.32e11 [A] (primary 2.57e11, w 2.6) | none | no | no | no |

- **Verdict: H1 FAILS.** No c0_soft point lies in M2-Q, which needs c0 ≥ 1.377, or in M2-X with e_crit ≈ 2×10⁻⁶.
- **The miss is close.** At δ = −4° the softest-mode point is 0.006 in c0 below a small M2-X-proper pocket
  (1.086–1.109). That pocket coincides with the N_YL interference minimum. This is a [C] observation; no mechanism is
  claimed.
- **Under a looser reading** (any e_crit), 2 of 3 points would lie in a free-X-with-e_crit region. But their
  self-consistent H_i ≈ 2.3–2.5×10¹¹ is not "on the computed closure": no LF lattice exists above 1.3×10¹¹, and the
  root is either extrapolated ×1.9 or obtained from a node that reuses the 1.3×10¹¹ lattice block [A].
- The looser reading is reported, not counted, and the strict cutoff was fixed with partial sight of the data (§2).

## 6. Q4 (lf_fermion NOTES §0 list, unchanged definitions)

- **Q4-a, e = e_crit.** This is the line where S/Q(10⁻⁶) = 3.47 at the root, i.e. the M2-X/M2-Q boundary: c0 ≈ 1.377
  at δ = 0, between 1.40 and 1.45 at +4°, below 1.38 at −4°, above 1.45 at −8° (grid resolution). It was declared
  branch-defining in advance, so it is **not counted**.
- **Q4-b, edge of bistability.** One bistable band at δ = 0: c0 = 1.3766–1.3778 [E, from c0-interpolation of the root
  S/Q at 1.374 / 1.378]. It is an artefact-prone feature of the sharp relic step at S/Q = 3.4716 (e2e_chain skeptic N3).
  Branch-defining, **not counted**. No other key on the grid has two roots.
- **Q4-c, release synchrony.**
  - t_release is defined here as the first time |N_X| reaches 90% of its final value, from the engine timelines (a model
    choice). On the 1.3×10¹¹ node it falls monotonically along δ = 0, from 4.67 (c0 = 1) to 1.47 H_i⁻¹ (1.60); over
    the whole grid the range is 1.42–4.7.
  - t_coh,end = 0.846 H_i⁻¹ is the first time f_σ + f_χ < 0.1 in the KFmid lattice. KFhi and the 8.4×10¹⁰ run never
    reach it.
  - No c0 in 1.0–1.6 gives synchrony. **No value.**
- **Answer to Q4:** the three listed principles are evaluable now. Two are branch-defining by construction, and the
  third gives no value. No principle becomes testable in a non-tuned way.

## 7. The lf_fermion skeptic's corrections, as applied here

- **K is taken from the kinetic hand-over, not from the [E] rule:** 5.53 (KFmid), brackets 5.01–5.85.
- **Direction of the Q-bias corrections** (N1 κ_frag lower; N2's early Q 9–19% high): both push toward the KFlo
  bracket. That bracket keeps every checked branch and raises the roots by about 5%. No re-calibrated lattice exists,
  so this is a bracket, not a correction.
- **Lattice cutoff:** κ_frag and the locality effect sit on hot spots 1–2 cells wide, and dx was not varied. This
  uncertainty is inherited by K and has not been quantified.
- **The 8.4×10¹⁰ LF node** has energy budget 0.97% and does not fragment. Its κ_frag was not checked with fermions.
  This is an inherited caveat on every root below 1.3×10¹¹, i.e. all M2 roots.

## 8. Limitations [A]

- No skeptic of this track.
- No LF lattice exists outside 8.4×10¹⁰–1.3×10¹¹. The M2 roots need the 5.5×10¹⁰ [A] node (and the H1 roots the
  2.7×10¹¹ [A] node), or extrapolation.
- Phase steps are ≥ 2°. HFb and Γ_th were not run. N_* and α_s were not re-solved.
- The relic's sharp branch step makes the bistable and empty bands model artefacts.
- The ray reduction and the hot-spot statistics of lf_fermion propagate unquantified into K.

## 9. Files

- **Code:** `code/lfc_common.py`, `build_hist.py`, `scan_run.py`, `fp_lf.py`, `post_lf.py`, and `edge.py` /
  `phase_soft.py` (byte-identical copies of pr_principle's).
- **Results:**
  - `results/windows.json`: the export, with the grid classes, representative points, brackets, Q3, H1, Q4 and a summary;
  - `fixpoints_LF.json`, `scan_hLF{_mid,84_mid,55_mid,270_mid,_lo,_hi}.json`, `splice_checks.json`, `handover_*.json`,
    `edge.json`, `phase_soft.json`.
- **Scratch:** `.../wt2/lf_chain/sp0/hist/bench_*.json` (final histories); `sp1`, `sp2` (engine files).

**Literature:** none newly verified this session (budget). β_iso < 0.038, Planck 2018 X (arXiv:1807.06211), is quoted
through the e2e_chain skeptic's verification [Q/R]. The method literature is lf_fermion's.
