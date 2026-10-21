# Rerun audit of Codex's 3D charged-exit cross-check (a69bd86)

21 October 2026, issue #13. Track: `rerun`. Codex packet read-only from
`git archive a69bd86 research/reality_framework/grut_expanding_charge_20261021` into scratch;
nothing in Codex's or Claude's packets was edited. Machine: 4 CPUs shared with other audit agents;
every solver job ran with `OMP_NUM_THREADS=2 --threads 2`, sequentially, each under 30 min wall.

Labels: [E] = estimate/diagnostic of a finite numerical run; [P] = preparation input chosen by hand.
Nothing here is evidence of exclusion, validation or identification of Candidate 1.

## Method

1. **Own reader** (`reader.py`; numpy/json only, no Codex analysis code imported or run). Recomputes
   from `raw_diagnostics.tar.gz` for every RUN_TABLE row: a_form (first scalar output with
   `fplus_05 >= 0.10` and `fminus_05 >= 0.10`), end f+/f-, last-record winding masks, max relative
   charge drift `max|Q(t)/Q(0)-1|`, and the conformal work residual `H(t)-H(0)-W(t)` under four
   normalisations. Output: `results/archived_recompute.json`.
2. **Independent preparation code** (`make_modes.py`, `make_noise.py`): re-implemented the 64-mode
   and Gaussian recipes from their written description. All 5 archived mode tables and all 3 Gaussian
   noise files regenerate **bit-for-bit** (sha256 equal to the archived/config hashes).
3. **Build**: `g++ -std=c++17 -O3 -fopenmp -ffp-contract=off` (ENVIRONMENT.json flags, no fast-math,
   no -march=native), g++ 13.3.0-6ubuntu2~24.04.1. Source sha256 e62701e5... (equals SOURCE_FREEZE);
   executable hash differs from Codex's (different build host), as REPRODUCE.md anticipates.
   Launch via own `run_case.sh` / `queue*.sh` (not run.py / gaussian_run.py).
4. **Own topology reconstruction** (`topo.py`): periodic winding from saved final states with a
   different algorithm from Codex's cell-BFS (scipy `ndimage.label` on the open box, periodic face
   gluing as a label graph with lattice-translation edges, inconsistent offsets = non-contractible
   cycle). Unit controls `test_topo.py` (slab 6, slab complement 6, corner blob 0, its complement 7,
   tube 1, diagonal line 3, contractible ring across a boundary 0): ALL PASS
   (`results/test_topo.txt`). Also evaluated at stricter thresholds |S|>0.7/F2<0.3, 0.8/0.25, 0.9/0.1.
5. `compare.py` compares rerun CSV/topology with the archive (bitwise per column);
   `aggregate.py` writes `results/summary.json`.

## (1) Archived diagnostics vs RUN_TABLE.md: row-by-row

| run | N | L | a_form (mine/table) | end f+ | end f- | winding +/- (any; largest) | budget pointwise (mine/table) | budget / H0 | budget / final abs work | max charge drift |
|---|---|---|---|---|---|---|---|---|---|---|
| main_s101 | 192 | 76.8 | none / none | 0.977672 | 2.11928e-06 | 7/0; 7/0 | 0.001166 / 0.00117 | 0.006861 | 0.001383 | 1.77e-14 |
| main_s202 | 192 | 76.8 | 1.60 / 1.60 | 0.442538 | 0.473921 | 7/7; 7/7 | 0.0007326 / 0.000733 | 0.004591 | 0.0008709 | 1.3e-14 |
| main_s303 | 192 | 76.8 | 1.63 / 1.63 | 0.687723 | 0.231248 | 7/7; 7/7 | 0.0008938 / 0.000894 | 0.006114 | 0.001046 | 1.98e-14 |
| cold_s101 | 192 | 76.8 | none / none | 1 | 0 | 7/0; 7/0 | 0.0007007 / 0.000701 | 0.001659 | 0.001211 | 2.56e-14 |
| kappa4_s101 | 192 | 76.8 | 1.93 / 1.93 | 0.117289 | 0.743885 | 0/7; 0/7 | 0.0007568 / 0.000757 | 0.004947 | 0.0008927 | 2.84e-14 |
| dt_half_s101 | 192 | 76.8 | none / none | 0.97772 | 1.97799e-06 | 7/0; 7/0 | 0.0002907 / 0.000291 | 0.001711 | 0.0003452 | 1.77e-14 |
| grid256_s101 | 256 | 76.8 | none / none | 0.978352 | 1.66893e-06 | 7/0; 7/0 | 0.001182 / 0.00118 | 0.006954 | 0.001402 | 4.11e-14 |
| grid384_s101 | 384 | 76.8 | none / none | 0.978905 | 8.3005e-07 | 7/0; 7/0 | 0.0005303 / 0.00053 | 0.003123 | 0.0006295 | 9.64e-14 |
| large_s101 | 288 | 115.2 | none / none | 0.981402 | 0 | 7/0; 7/0 | 0.001022 / 0.00102 | 0.00611 | 0.001225 | 6.06e-14 |
| gaussian_s101 | 192 | 76.8 | 1.57 / 1.57 | 0.427757 | 0.516638 | 7/7; 7/7 | 0.0008995 / 0.0009 | 0.005252 | 0.001084 | 3.02e-14 |
| gaussian_s202 | 192 | 76.8 | 1.57 / 1.57 | 0.470218 | 0.469507 | 7/7; 7/7 | 0.0008578 / 0.000858 | 0.005061 | 0.001032 | 3.6e-14 |
| gaussian_large_s101 | 288 | 115.2 | 1.57 / 1.57 | 0.459575 | 0.487886 | 7/7; 7/7 | 0.0009232 / 0.000923 | 0.005361 | 0.001114 | 4.13e-14 |

- a_form, end f+/f- (6 significant digits), winding masks, N, L: **exact agreement, 12/12 rows.**
  Gaussian formation interval: previous output a=1.54, first output a=1.57 in all three runs
  (Codex's "a in (1.54, 1.57]": CONFIRMED). Max charge drift 9.64e-14 (grid384): CONFIRMED.
- Budget: my first-pass definition `max|res| / max(H)` disagreed with the table in the 3rd digit for
  4 rows (main_s101 0.001152 vs 0.00117; dt_half 0.000287 vs 0.000291; grid256 0.001168 vs 0.00118;
  grid384 0.000524 vs 0.00053). With the pointwise denominator `max(|H0|,|H(t)|,|W_net(t)|)`
  (found in analyze.py only after the mismatch) all 12 rows agree to 3 digits; max 0.0011816
  (grid256) = the README's "<= 0.1182 percent": CONFIRMED. PROTOCOL.md specifies
  `max(H0, H(t), accumulated absolute work)`; |W_net| <= accumulated |work|, so Codex's number is an
  upper bound on the protocol quantity (conservative). Normalised by H0 alone the residual is
  0.17-0.70 % (H grows ~6x over the window) - still below the 1 % gate. The residual is a timestep
  error: dt_half/main_s101 = 0.2493 (second order: 4.01x).

## (2) End-to-end reruns with the compiled frozen solver

### Full-size reruns of archived rows (exact archived arguments, 2 threads instead of 8)

| row | count columns (10) bitwise identical at all 41 outputs | topology JSONL identical (all records) | max rel. diff of reduction columns | charge drift rerun / archive | budget rerun / archive | wall s |
|---|---|---|---|---|---|---|
| cold_s101 | True | True (11) | 9.5e-11 | 9.28e-14 / 2.56e-14 | 0.0007007 / 0.0007007 | 235 |
| gaussian_large_s101 | True | True (11) | 1.7e-12 | 1.78e-13 / 4.13e-14 | 0.0009232 / 0.0009232 | 864 |
| gaussian_s101 | True | True (11) | 7.7e-13 | 8.43e-14 / 3.02e-14 | 0.0008995 / 0.0008995 | 234 |
| gaussian_s202 | True | True (11) | 6.9e-13 | 1.19e-13 / 3.60e-14 | 0.0008578 / 0.0008578 | 213 |
| kappa4_s101 | True | True (11) | 4.9e-13 | 1.13e-13 / 2.84e-14 | 0.0007568 / 0.0007568 | 231 |
| large_s101 | True | True (11) | 2.0e-11 | 1.40e-13 / 6.06e-14 | 0.001022 / 0.001022 | 590 |
| main_s101 | True | True (11) | 1.3e-12 | 1.05e-13 / 1.77e-14 | 0.001166 / 0.001166 | 194 |
| main_s202 | True | True (11) | 1.2e-12 | 1.28e-13 / 1.30e-14 | 0.0007326 / 0.0007326 | 254 |
| main_s303 | True | True (11) | 2.4e-13 | 1.40e-13 / 1.98e-14 | 0.0008938 / 0.0008938 | 211 |

Field evolution is per-site deterministic (only diagnostics use parallel reductions), so with a
different thread count every count-based column (f+-, 0.4/0.5/0.6 cuts, genuine+-, min/max S) and every
topology record reproduce **bit-for-bit**; energy/charge integrals differ only at <=1e-11 relative
(summation order; mean_S and the linear projection are near-zero sums and differ more in relative terms).
Charge drift at 2 threads is 6e-14 to 2e-13 vs Codex's 1e-14 to 1e-13: same roundoff class, both
<< the 1e-9 gate.

### Reduced-cost reruns (what changed is exact)

All: dt=0.06, T=48 (a 1 -> 2.2), rate 1/40, kappa=2, sample 1.2, topology every 4.8.
(a) Same box L=76.8, coarser grid: N=128 dx=0.6 and N=96 dx=0.8. Sparse s101 uses the identical
mode table (modes are L-defined). Gaussian s202 = Fourier restriction of the identical N=192 draw
(coarse field equals the fine field at coincident points to 3.5e-6 for dx=0.6, 1.7e-3 for dx=0.8,
unit-rms field) [P]. Walls are under-resolved at dx=0.6/0.8 (comoving width ~sqrt2/a ~ 0.64 at a=2.2).
(b) Same dx=0.4 (same L/N), smaller box N=128 L=51.2: sparse modes redrawn for L=51.2 (same seed,
same recipe; a different preparation), Gaussian native N=128 draw (a different realization).

| run | a_form | end f+ | end f- | max f- (a) | end winding +/- | 7/7 at every post-formation record | charge drift | budget |
|---|---|---|---|---|---|---|---|---|
| red128dx06_main_s101 | none | 0.9737 | 4.77e-07 | 0.0837 (1.69) | 7/0 | n/a | 3.60e-14 | 1.113e-03 |
| red96dx08_main_s101 | none | 0.9704 | 0 | 0.0865 (1.69) | 7/0 | n/a | 3.93e-14 | 1.009e-03 |
| red128dx06_gaussian_s202 | 1.57 | 0.4693 | 0.47 | 0.4698 (2.20) | 7/7 | True | 8.66e-14 | 8.176e-04 |
| red96dx08_gaussian_s202 | 1.57 | 0.4698 | 0.47 | 0.4701 (2.20) | 7/7 | True | 2.94e-14 | 7.669e-04 |
| red128dx04_main_s101 | none | 0.9680 | 0.00184 | 0.0353 (1.75) | 7/0 | n/a | 4.76e-14 | 1.133e-03 |
| red128dx04_gaussian_s202 | 1.57 | 0.4921 | 0.449 | 0.4488 (2.20) | 7/7 | True | 6.44e-14 | 8.696e-04 |

**Dichotomy persists**: s101 is one-sign at every resolution of the same preparation (dx 0.2..0.8,
archived + mine: end f- <= 2.2e-6, formation threshold never reached) and in its L=51.2 redraw
(peak f- 3.5 %, end 0.18 %); Gaussian s202 forms at a=1.57 and stays
two-sign ~47/47 with 7/7 winding at dx 0.6 and 0.8 (end f+/f- within 0.001 of the 192^3 values) and
in a fresh L=51.2 realization.

### Additional seeds at reduced cost (N=128, dx=0.4, L=51.2) [E]

| run | a_form | end f+ | end f- | max f- (a) | end winding +/- | 7/7 at every post-formation record | charge drift | budget |
|---|---|---|---|---|---|---|---|---|
| seed128dx04_gaussian_s303 | 1.57 | 0.4604 | 0.48 | 0.4797 (2.20) | 7/7 | True | 5.93e-14 | 8.819e-04 |
| seed128dx04_gaussian_s404 | 1.57 | 0.6277 | 0.319 | 0.4156 (1.72) | 7/5 | False | 6.71e-14 | 9.495e-04 |
| seed128dx04_gaussian_s505 | 1.57 | 0.4799 | 0.466 | 0.4663 (2.20) | 7/6 | False | 5.98e-14 | 9.067e-04 |
| red128dx04_gaussian_s202 | 1.57 | 0.4921 | 0.449 | 0.4488 (2.20) | 7/7 | True | 6.44e-14 | 8.696e-04 |
| red128dx04_main_s101 | none | 0.9680 | 0.00184 | 0.0353 (1.75) | 7/0 | n/a | 4.76e-14 | 1.133e-03 |
| red128dx04_main_s202 | 1.63 | 0.2590 | 0.677 | 0.6773 (2.20) | 7/7 | True | 5.22e-14 | 1.126e-03 |
| seed128dx04_main_s303 | 1.60 | 0.4848 | 0.404 | 0.4223 (1.66) | 7/7 | True | 5.88e-14 | 6.540e-04 |
| seed128dx04_main_s404 | 1.54 | 0.4628 | 0.478 | 0.4779 (2.20) | 7/7 | True | 8.66e-14 | 7.692e-04 |
| seed128dx04_main_s505 | 1.63 | 0.9168 | 0.0362 | 0.3268 (1.66) | 7/4 | False | 5.97e-14 | 1.344e-03 |

Formation (both raw fractions >= 10 %) occurred in 4/4 Gaussian preparations at L=51.2, all at
a=1.57. Two-sign fractions survived to a=2.2 in all 4 (minimum end f- 0.319), but **all-axis winding
of both phases at every post-formation record held in only 2/4** (s202-native, s303): Gaussian s404's
negative phase (38-41 % of volume) winds only x,z (mask 5) from a=1.60 on, while its fractions drift
0.45/0.38 -> 0.63/0.32; Gaussian s505's negative phase drops from 7 to 6 at a=2.08. Sparse seeds at
L=51.2: s202, s303, s404 two-sign 7/7; s505 forms at a=1.63, both phases wind 7/7 at a=1.72 and 1.84
(negative phase 24.5 % of volume, one component), then the negative phase shrinks to 3.6 % by a=2.2.

### Large-box rows Codex left unexecuted (N=288, dx=0.4, L=115.2; frozen source, archived argument set) [E]

Codex's RUN_STATUS lists large_s202 and large_s303 as `not_launched_cost_cap`. Run here with the same
arguments as large_s101 (modes regenerated by the same recipe for L=115.2: s202 sha256 8f8e746f...,
s303 e87bb60d...). New evidence, not a reproduction.

| run | a_form | end f+ | end f- | max f- (a) | end winding +/- | 7/7 at every post-formation record | charge drift | budget |
|---|---|---|---|---|---|---|---|---|
| new_large_s202 | 1.60 | 0.4119 | 0.502 | 0.5019 (2.20) | 7/7 | True | 2.08e-13 | 7.434e-04 |
| new_large_s303 | 1.60 | 0.5193 | 0.39 | 0.4186 (1.66) | 7/7 | True | 2.53e-13 | 7.295e-04 |

## (3) Topology criterion on my own reconstruction

- At the packet threshold (|S|>0.5, F2<0.5) my independent algorithm reproduces the solver's own
  topology record **exactly** (components, occupied and largest fractions, any/largest winding) on
  all 9 saved states (full main_s101, gaussian_s101, gaussian_s202 at a=2.2; four reduced runs at
  a=2.2; gaussian_s202 at a=1.60; main_s101 at a=1.69). The archived JSONL itself contains only
  summaries, so the archive alone cannot be re-segmented; the reruns reproduce it bit-for-bit, which
  makes the reconstruction a check on the archived records.
- Threshold robustness at a=2.2: Gaussian s101/s202 both phases still wind 7/7 at |S|>0.9 & F2<0.1
  (largest components 0.31/0.38 and 0.34/0.34 of the volume). At a=1.60 (first post-formation record)
  Gaussian s202 phases wind 7/7 at every threshold, even at |S|>0.9 where each occupies only ~8.7 %.
- **Winding is not a persistence indicator.** main_s101 at a=1.69 (between archived topology records
  1.60 and 1.72): the minority negative phase occupies 8.17 % and its largest component (3.6 %) winds
  all three axes by both the solver's diagnostic and my reconstruction; it is non-winding at the next
  archived record (1.72) and extinct (2e-6) at a=2.2. At |S|>0.7 it does not wind. Codex's own archive
  has a second case: kappa4_s101 has both phases 7/7 at a=1.96 and 2.08, and the positive phase loses
  all winding by a=2.2 (occupancy 0.286 -> 0.116). Third case: sparse s505 at L=51.2 above. Winding
  sets in at low filling (main_s303 positive phase 7 at 3.8 % occupancy, a=1.60).
- Bit semantics: an axis bit is set when some non-contractible cycle has a nonzero component along
  it; a single diagonal (1,1,0) cycle sets x and y (unit test diag_line -> 3). "7" therefore does not
  by itself certify three independent box-spanning cycles. For the ~45 % bicontinuous Gaussian phases
  this does not change the conclusion.

## Other observations [E]

- s101 is not a "nothing forms" history: at a=1.69 it is an 81:8 two-sign state (f- peak 8.17 %
  at dx=0.4; 8.09 % at 256, 8.04 % at 384, 8.37 % at dx=0.6, 8.65 % at dx=0.8; 9.24 % with the |S|>0.4
  cut) that then becomes one-sign. "none" means the 10 % formation threshold was missed by ~2 points;
  the outcome (minority extinction) is robust across all resolutions and dt. In the L=115.2 box the
  peak is only 2.9 %.
- Core velocity proxy v2_core at a=2.2 for s101: 1.26 (dx 0.4), 1.59 (0.6), 2.36 (0.8): strongly
  resolution dependent and >1, consistent with Codex's statement that the proxy fails as a wall speed.
- Charge location (own reconstruction, a=2.2): only 38 % (Gaussian s101/s202) of the total charge sits
  in |S|<0.5 cores and ~60 % in cells with F2<0.5. Not part of this track; flagged for the charge/
  estimator audit.

## Verdicts

See `results/summary.json` and the structured return. Short form:
- 12/12 archived rows recomputed by own reader: CONFIRMED (budget only under Codex's pointwise
  normalisation; first-pass normalisation differs in the 3rd digit for 4 rows - definitional).
- Formation interval (1.54, 1.57] for the 3 Gaussian runs: CONFIRMED.
- Max charge drift 9.64e-14, work residual <= 0.1182 % (pointwise denominator): CONFIRMED;
  residual is 0.70 % if normalised by H0 alone (still under the 1 % gate).
- 9 archived rows rerun end-to-end at full size (main_s101/s202/s303, cold_s101, kappa4_s101,
  large_s101, gaussian_s101/s202, gaussian_large_s101): CONFIRMED bit-for-bit in all count and
  topology data; integrals to <=2e-11 relative.
- s101 one-sign vs Gaussian two-sign dichotomy at reduced resolution (dx 0.6, 0.8 same preparation;
  dx 0.4 in L=51.2): CONFIRMED.
- The two unexecuted large-box rows (large_s202, large_s303), run here with the frozen arguments:
  both form at a=1.60 and keep 7/7 winding at all 6 post-formation records (end 0.412/0.502 and
  0.519/0.390). New evidence consistent with Codex's L=76.8 pattern; not Codex's claim.
- "Both phases wind all axes at every post-formation record to a=2.2" (3 Gaussian runs, L>=76.8):
  CONFIRMED for those runs; CORRECTED in scope: in a smaller L=51.2 box it held for 2/4 Gaussian seeds,
  and box winding at a sampled time does not imply survival (three counterexamples).

## Limits

- Same frozen C++ source and same preparations: this is a reproduction, not an independent
  integrator. Determinism makes the full-size agreement a check of archive integrity/provenance,
  not of the physics or discretisation.
- grid384_s101, grid256_s101 and dt_half_s101 not rerun end-to-end (grid384 ~45+ min at 2 threads);
  they are checked only through the archived-diagnostic recomputation.
- Reduced runs at dx=0.6/0.8 under-resolve walls; extra seeds use a smaller box (L=51.2) where
  finite-box effects are larger. Seed counts are tiny; no probabilities are inferred.
- Topology checked only at saved states (final, a=1.60, a=1.69); sampled cadence Delta a = 0.12.
- No physical calibration, no cosmological persistence, no observational comparison.

## Files

Code: `reader.py`, `make_modes.py`, `make_noise.py`, `run_case.sh`, `queue1.sh`..`queue4.sh`,
`topo.py`, `test_topo.py`, `compare.py`, `compare_full_all.sh`, `aggregate.py`, `fill_notes.py`.
Results: `results/archived_recompute.json`, `results/rerun_*.json`, `results/topo_*.json`,
`results/noise_*.json`, `results/test_topo.txt`, `results/summary.json`.
Raw rerun outputs (CSV, topology JSONL, meta, slices) are in scratch only; the saved f64 states
(2 GB) were deleted after the topology analysis (rerunning the queue recreates them):
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun/runs/`.
