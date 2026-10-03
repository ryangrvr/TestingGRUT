# DA0 · C3 — PRIMARY MODEL 2 (M2-D): STAGE A RESULT (deterministic)

**Status:** Stage A was run exactly on the preregistered grid. **Numerical illustration** (independent code path, not
independent reviewer). Exact statements are marked **INTERNALLY PROVED / NOT EXTERNALLY REVIEWED**. No novelty is
claimed. **C3-B is closed. No noise, no fluctuations, no PM2-BIS.** This document stops for owner review.

**Frozen inputs (unchanged after results):**
- Law: PM2-F0, with γ = 1/2, ν = κ = L_e = 1.
- Drive: PM2R-01, h_i = +1 and h₀ = −(N − 1).
- Geometry: L × L open square lattice, outlet at corner node 0.
- Initial ensembles:
  - uniform C ≡ 1;
  - seeds 1 – 20, C⁰ = 1 + δU(−1, 1) with δ = 0.01.
- Statistics: PM2R-04 / 05 / 06.
- Comparators: R-TREE (Wilson, seeds 1 – 20), SP-TREE, A0, A1, A2.
- Sizes: L ∈ {16, 24, 32, 48, 64}.

**Implementation** (ledger PM2-I0 … I5, committed and pushed as `c6de0fa` **before** the grid ran):
- integrator: explicit Euler in y = √C, the same ODE, at max_rel_change 0.0125;
- exact spanning-tree early stop;
- other tolerances as registered;
- grading on ensembles only;
- a preregistered robustness pass at max_rel_change 0.05 (PM2-I4).

**Code and logs:**

| file | content |
|---|---|
| `c3/pm2_stageA.py` | Stage-A code |
| `c3/pm2_stageA.log` | production run (step 0.0125) |
| `c3/pm2_stageA_rel005.log` | robustness pass (step 0.05) |
| `c3/pm2_stageA_summary.py` | aggregation |
| `c3/pm2_stageA_summary.log` | aggregation output |
| `c3/pm2_earlystop_check.*`, `c3/pm2_yvar_check.*`, `c3/pm2_dt_check.*` | implementation validation |

## 0. Two exact statements used in reading the results (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED)

**PROP PM2-U (the uniform start cannot reach a tree under exact dynamics, for even L).**
1. The diagonal reflection σ: (x, y) ↦ (y, x) fixes the outlet, the drive and C ≡ 1. The exact flow therefore preserves
   σ-symmetric states.
2. σ fixes no lattice edge. An edge (x, y)–(x+1, y) maps to (y, x)–(y, x+1), which is a different edge. So every edge has
   a distinct partner.
3. Pruning is σ-equivariant, so a symmetric active set has **even** size.
4. A spanning tree has N − 1 = L² − 1 edges, which is **odd** for even L.
5. Hence the exact uniform-start trajectory never reaches a spanning tree at any L in the grid.

**Consequence.** Every uniform-start tree in the logs is selected by **round-off symmetry breaking**. This agrees with
PM2-I3/I4: its topology does not converge under step refinement. Uniform-start trees are reported but used for nothing.

**PROP (from PM2R-06) n_reroute(f, e) = S_e.** This holds for every tree. Therefore **P3 grows with N for every tree
family whose subtree sizes grow**, R-TREE, SP-TREE and A2 included. ρ > 0 is generic tree combinatorics unless it
exceeds the comparators. This was already registered as "structural only" (PM2R-06).

## 1. Results (production run, step 0.0125)

**Convergence.** 100 / 100 perturbed runs and 5 / 5 uniform-start runs end on a connected spanning tree (N − 1 active
edges). There are no non-tree or step-capped runs.

### 1.1 Per-size ensemble statistics

Each cell gives the mean and 95% t-interval over 20 seeds. SP and A2 are single deterministic values.

**τ (P1):**

| L | PM2 | R-TREE | SP | A2 |
|---|---|---|---|---|
| 16 | 0.992 [0.948, 1.037] | 0.439 [0.386, 0.491] | 1.479 | 0.394 |
| 24 | 0.873 [0.821, 0.925] | 0.416 [0.378, 0.454] | 1.688 | 0.388 |
| 32 | 0.786 [0.749, 0.823] | 0.384 [0.355, 0.413] | 1.816 | 0.360 |
| 48 | 0.668 [0.634, 0.702] | 0.397 [0.376, 0.417] | 2.040 | 0.342 |
| 64 | 0.599 [0.577, 0.621] | 0.385 [0.362, 0.407] | 2.206 | 0.334 |

**η_H (P2):**

| L | PM2 | R-TREE | SP | A2 |
|---|---|---|---|---|
| 16 | 0.669 [0.644, 0.693] | 0.672 [0.633, 0.711] | 0.871 | 0.661 |
| 24 | 0.660 [0.643, 0.676] | 0.679 [0.661, 0.697] | 0.945 | 0.660 |
| 32 | 0.646 [0.625, 0.666] | 0.674 [0.654, 0.694] | 0.965 | 0.665 |
| 48 | 0.615 [0.592, 0.637] | 0.656 [0.630, 0.682] | 0.984 | 0.691 |
| 64 | 0.614 [0.600, 0.628] | 0.654 [0.640, 0.668] | 0.991 | 0.669 |

**P3 median:**

| L | PM2 | R-TREE | SP | A2 |
|---|---|---|---|---|
| 16 | 19.1 [18.0, 20.2] | 36.0 [29.4, 42.5] | 11 | 33 |
| 24 | 31.7 [29.6, 33.7] | 61.2 [48.8, 73.5] | 17 | 49 |
| 32 | 46.5 [44.3, 48.7] | 104.8 [86.3, 123.2] | 23 | 86 |
| 48 | 82.4 [78.2, 86.6] | 191.1 [159.9, 222.2] | 34 | 147 |
| 64 | 122.2 [115.7, 128.7] | 307.7 [243.7, 371.6] | 45 | 260 |

**Strahler order of the outlet:**

| L | PM2 | R-TREE | SP | A2 |
|---|---|---|---|---|
| 16 | 4.00 [4.00, 4.00] | 4.15 [3.92, 4.38] | 2 | 4 |
| 24 | 4.35 [4.12, 4.58] | 4.90 [4.69, 5.11] | 2 | 5 |
| 32 | 4.95 [4.85, 5.05] | 5.10 [4.89, 5.31] | 2 | 5 |
| 48 | 5.65 [5.42, 5.88] | 5.65 [5.42, 5.88] | 2 | 6 |
| 64 | 5.90 [5.76, 6.04] | 6.00 [6.00, 6.00] | 2 | 6 |

### 1.2 Fits across L

| fit | PM2 | R-TREE | SP-TREE | A2 |
|---|---|---|---|---|
| ρ: log(P3 median) vs log N | **0.674 [0.651, 0.698]** | 0.784 [0.705, 0.863] | 0.507 | 0.750 |
| Strahler slope vs log₂N | 0.509 [0.459, 0.560] | 0.447 [0.389, 0.506] | 0.000 | 0.508 |

PM2 and R-TREE intervals are pooled OLS fits with 95% t-intervals; SP-TREE and A2 are 5-point OLS fits.

### 1.3 Other ensemble diagnostics

**Energies and tree counts:**

| L | Ē perturbed, min / median / max | spread (max − min)/median | Ē uniform start (round-off) | distinct trees, seeds 1 – 20 |
|---|---|---|---|---|
| 16 | 59.38 / 60.94 / 61.90 | 4.1% | 61.52 | 20 |
| 24 | 99.41 / 101.54 / 102.81 | 3.3% | 102.87 | 20 |
| 32 | 140.68 / 145.14 / 149.29 | 5.9% | 149.22 | 20 |
| 48 | 233.23 / 236.97 / 243.27 | 4.2% | 249.27 | 20 |
| 64 | 330.29 / 335.06 / 342.27 | 3.6% | 360.12 | 20 |

**Static straight-path barrier diagnostic** (200 swaps per L, PM2R-06):

| L | median raw | median reduced | max raw | fraction > 0 |
|---|---|---|---|---|
| 16 | 9.44 | 0.235 | 1.2×10³ | 0.975 |
| 24 | 12.76 | 0.185 | 3.3×10³ | 0.970 |
| 32 | 14.35 | 0.141 | 4.2×10³ | 0.945 |
| 48 | 8.54 | 0.049 | 1.3×10⁴ | 0.895 |
| 64 | 20.18 | 0.079 | 8.4×10³ | 0.925 |

**Uniform-start trees** (round-off-selected; PROP PM2-U):
- They are SP-like: τ runs from 1.54 to 2.04, η_H from 0.85 to 0.96, Strahler is 2 – 3.
- They are never among the seeded trees.
- Their Ē is above every perturbed seed at L ≥ 48.

**Controls:**
- **A0** (C ≡ 1): all edges are active, so there is no selection. Ē is 2.6×10³ at L = 16 and 1.6×10⁵ at L = 64.
- **A1** (C relaxed under the fixed uniform-network flow Q⁰): **all edges stay active at every L** (8064 / 8064 at
  L = 64). Without conductance → flow feedback, no tree is selected.
- **A2** (supplied recursive-bisection hierarchy): τ, η_H and Strahler are **within or near the R-TREE intervals** at
  every L.

## 2. Preregistered criteria

| item (Revised Stage-A positive) | rule | production (0.0125) | robustness pass (0.05) |
|---|---|---|---|
| 1. sparse / tree architecture from the fixed homogeneous law | trees reached | **YES** (105 / 105; A0, A1: no) | YES |
| 2a. P1 | at L = 48 and 64: PM2 τ-interval disjoint from R-TREE, and τ_SP outside it | **PASS** (both sizes) | PASS |
| 2b. P2 | same rule for η_H | **FAIL**: L = 48 overlaps ([0.592, 0.637] vs [0.630, 0.682]); L = 64 is disjoint | FAIL |
| 3. P3 | ρ > 0 (declared reading: 95% lower bound > 0) | **PASS** (0.674 [0.651, 0.698]) | PASS (0.670 [0.647, 0.693]) |
| 4. Strahler growing beyond the controls | declared conservative reading (`pm2_stageA_summary.py` header) | **FAIL**: PM2 slope overlaps R-TREE; at L = 48 the intervals are identical; at L = 64 PM2 (5.90) is not above R-TREE (6.00) | FAIL |
| 5. ensemble robustness | identical verdicts at both step sizes (PM2-I4) | **MET**: all verdicts identical (§2.1) | |

**Disclosure on item 4's reading.** The charter gives no numeric rule for item 4. The reading used here was written into
the aggregation script after the L = 16 numbers were visible and before L ≥ 24 were. Under any reasonable reading the
outcome does not change: PM2's Strahler order never exceeds R-TREE's at any L, and the slopes are statistically
indistinguishable.

### 2.1 Robustness pass (PM2-I4)

Re-run of the full grid with y-Euler at max_rel_change = 0.05 (`c3/pm2_stageA_rel005.log`). **Every verdict is
identical**, so item 5 is **MET** for all criteria:
- 1 = YES (all runs are trees);
- P1 = PASS;
- P2 = FAIL;
- P3 = PASS;
- Strahler = FAIL.

| quantity at L = 64 | step 0.0125 | step 0.05 |
|---|---|---|
| τ | 0.599 [0.577, 0.621] | 0.604 [0.579, 0.629] |
| η_H | 0.614 [0.600, 0.628] | 0.610 [0.592, 0.629] |
| P3 median | 122.2 | 121.3 |
| Strahler | 5.90 [5.76, 6.04] | 5.95 [5.85, 6.05] |
| ρ (all L) | 0.674 [0.651, 0.698] | 0.670 [0.647, 0.693] |
| Strahler slope (all L) | 0.509 [0.459, 0.560] | 0.510 [0.461, 0.559] |

- **τ drift at 0.05:** 1.017 → 0.908 → 0.798 → 0.680 → 0.604, the same trend.
- **Static barriers at 0.05:** the median reduced barrier is 0.227, 0.178, 0.155, 0.062 and 0.103 for L = 16 … 64, again
  with no growth.
- **Conclusion:** individual trees are discretisation-sensitive (PM2-I3), but every **ensemble** statistic and verdict is
  stable to well within its own interval.

## 3. Answers to the Stage-A questions (exactly as registered)

**1. Convergence to sparse / tree architecture from homogeneous initial conductance, without a supplied hierarchy?**
**Yes**, for ensemble (ii). Every perturbed run ends on a spanning tree. For ensemble (i), the trees are round-off-selected
(PROP PM2-U). The selection requires flow feedback: A1 keeps the full lattice. This is the **KNOWN** concave-cost tree
property (PM2-F6), so "rediscovering multiplicity earns nothing".

**2. An asymptotic collective invariant (P1 / P2) that differs from R-TREE and SP-TREE?**
**Not established.**
- **P1** passes the preregistered finite-size rule at L = 48 and 64. However, τ_PM2 **decreases monotonically with L**
  (0.99 → 0.87 → 0.79 → 0.67 → 0.60) toward the R-TREE / A2 band (≈ 0.33 – 0.40), which itself is flat. The separation
  is shrinking: the gap is 0.55 at L = 16 and 0.21 at L = 64. **Asymptotic distinctness is not supported by the data.**
  No extrapolation rule was registered, so none is applied.
- **P2** fails the rule (overlap at L = 48).
- The P1 pass is therefore a **finite-size** difference: the deterministic dynamics at these sizes yield trees with
  steeper subtree-size tails than uniform spanning trees. It is not shown to be an asymptotic invariant.

**3. Does the number or depth of collective symmetry-inequivalent structures grow with N?**
- **Count:** 20 / 20 seeds give distinct trees at every L, so there are ≥ 10 classes modulo the order-2 reflection (an
  exact lower bound). This is a sampling diagnostic bounded by 20 and **not a growth criterion** (PM2R-03).
- **Depth:** outlet Strahler order grows ∝ log N **at the same rate as R-TREE and the supplied A2** (0.51 vs 0.45 vs
  0.51 per doubling of N, with overlapping intervals). It never exceeds R-TREE. **No growth beyond generic tree
  combinatorics.**

**4. Do elementary reroutings involve source mass growing with N (P3)?**
**Yes structurally (ρ = 0.67 > 0), but generically.**
- By n_reroute = S_e, every growing tree family has ρ > 0. PM2's ρ is **below** R-TREE's (0.78) and A2's (0.75), and its
  P3 medians are about 2.5× **smaller** than R-TREE's at L = 64.
- The swaps of PM2 trees are less collective than those of random trees.
- Per PM2R-06, this is structural only; it establishes neither metastability nor transitions.

**5. Do barriers between collective alternatives grow with N?**
**Not established; the static diagnostic does not indicate growth.**
- The median raw straight-path barrier is O(10) and roughly flat (9 – 20).
- The median reduced barrier falls from 0.23 to 0.05 – 0.08.
- Pruning is absorbing (PM2R-02), so dynamical barriers are outside Stage A. Barrier diagnostics never upgrade Stage A.

## 4. Terminal

**The Stage-A positive is NOT met.** Item 4 fails. Item 2 is met only by a finite-size P1 difference that is drifting
toward the random-tree band, and item 3 is generic. **PM2-A PARTIAL is not reached.**

**Proposed PM2 Stage-A terminal (owner to rule): PM2-B — GENERATED TREE / NETWORK — GENERIC TOPOLOGICAL HIERARCHY
ONLY.**

Reasons:
- The fixed homogeneous law generates spanning trees endogenously, and the selection is flow-feedback-dependent (A1
  negative). This is KNOWN behaviour of concave-cost transport networks.
- Hierarchy depth (Strahler) and rerouting collectivity (P3) are those of generic random spanning trees or less.
- The one surviving difference, P1 at finite L, is not shown asymptotic.

**Residuals preserved (not resolved by Stage A):**
- **(R1)** The asymptotic value of τ_PM2: does it converge to the R-TREE value, to an OCN-like value, or to a distinct
  value? This is a larger-L question; no extrapolation was registered.
- **(R2)** Whether trees dynamically selected under a different but preregistered physical mechanism (e.g. fluctuating
  sources, PM2-BIS) behave differently. **Not opened.**
- **(R3)** Detector limitation: **A2, a supplied hierarchy, is also indistinguishable from R-TREE on τ, η_H and
  Strahler.** The registered invariants have weak power to separate hierarchical from random trees on the lattice. A
  negative on them is therefore a negative for *generated hierarchy as measured by these invariants*, not a universal
  no-go.

**Not claimed:** C3-A1, PM2-A, a novel network law, or any statement about OCNs or rivers beyond the registered
comparators.

## 5. Information accounting

**Supplied (priced; unchanged from the charter):**
- the lattice and embedding;
- the outlet at corner 0 and the open boundary;
- the homogeneous source density h = +1 and the outlet sink;
- the adaptation ODE and γ, ν, κ, L_e;
- the initial ensembles:
  - uniform: 0 random numbers;
  - perturbed: |E| i.i.d. uniform numbers per seed, i.e. **O(N) random numbers per run, which select the instance**
    (PM2R-03);
- the drive power (dissipation) and the material cost.

**Implementation (not model) parameters:** max_rel_change, prune floor 10⁻¹², stationarity tolerance 10⁻⁶, 2000 stable
steps, step cap 4×10⁵ (never reached), y_ref = 10⁻³, and the exact early stop. Declared in PM2-I0 … I5 before the grid
ran.

**Derived (generated by the law, as measured):**
- spanning-tree topology: yes, a known property;
- a branch hierarchy with Strahler ∝ log N: **generic**, matching R-TREE;
- finite-size P1 steepness: yes, **not shown asymptotic**;
- collective rerouting: generic (ρ below the controls);
- metastable landscape and barriers: **not established**.

**Information bookkeeping for the selected architecture.** Per PM2R-03, the seed (O(N) random bits) selects which tree is
reached. The law contributes the *ensemble* statistics in §1, and these do not exceed generic tree combinatorics except
for the finite-size τ.

**TRUE COMPRESSION: 0.** No supplied architecture primitive was replaced by a law-generated one beyond generic tree
statistics.

**C2 architecture primitive.** Stage A does **not** show the C2 primitive to be replaceable by this bounded adaptive
mechanism, at the level of growing collective architecture beyond random trees.

## 6. Stop

Stage A is complete. Implementation changes are in the ledger (PM2-I0 … I5, before the grid). Nothing was varied after
results: γ, ν, κ, δ, geometry and drive are unchanged. **No PM2-BIS, no noise, no C3-B. Stopping for owner review.**
