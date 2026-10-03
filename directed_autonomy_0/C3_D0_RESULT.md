# DA0 · C3-D0 — ARCHITECTURE DETECTOR QUALIFICATION: RESULT

**Status:** D0 was run exactly as chartered (`C3_D0_ARCHITECTURE_DETECTOR.md` plus REPAIR 01 §R, pushed at `06cce52`
before any control value existed). This is a **numerical illustration** (independent code path, not independent
reviewer). It is an instrument result, **not GRUT physics**. No PM2 tree was scored. No PM2-BIS, no C3-B. **Stopping for
owner review.**

**Code and logs (`c3/`):**

| file | content |
|---|---|
| `d0_detector.py` | the Φ detector, the P1 / ε constructions and the secondary diagnostics |
| `d0_tests.py`, `d0_tests.log` | implementation-firewall tests |
| `d0_run.py`, `d0_run.log` | the full run (294 trees) |
| `d0_results.json` | per-tree values |

## 0. Implementation firewall (D0-O6)

**All 36 unit tests passed on the first execution** (`d0_tests.log`):
- relabelling and child-order invariance on synthetic trees;
- isomorphic trees built in different orders give equal Φ;
- **N1:** Φ = 0 with |U| = 0, and the exact class count is 2L − 1, at all 7 grid sizes (LEMMA D0-N1 confirmed);
- path: Φ = 0;
- **complete binary tree of depth 6:** Φ = 6/7 with |U| = 14 and 2 classes, as hand-derived;
- P0, P1 and the ε-trees are valid lattice spanning trees;
- P0-ε with b = 1 reproduces Stage-A P0 exactly.

**Inspection record:**
- No implementation change followed any look at a control value. **No debugging inspection of control Φ occurred.**
- Before the run, one structural check confirmed that the ε-trees differ from P0 / P1 and between seeds. It compared edge
  counts only; no Φ was computed.
- Per-tree timing was measured on N0 (L = 243, seed 1), P1-ε (L = 243, seed 1) and N0 (L = 128, seed 1), printing only
  wall time.

## 1. Primary Φ

N0 and the ε-families: mean ± 95% t half-width over seeds 1 – 20. N1, P0 and P1 are deterministic. |U| is the number of
units, and "cls" is the number of distinct skeleton classes (an ensemble mean where applicable).

**Grid G2:**

| L | N0 Φ | N0 \|U\| / cls | N1 Φ | P0 Φ | P0 \|U\| / cls | P0-ε Φ | P0-ε \|U\| / cls |
|---|---|---|---|---|---|---|---|
| 16 | 0.250 ± 0.115 | 3.1 / 2.2 | 0 | 0.667 | 3 / 1 | 0.567 ± 0.073 | 3 / 1.3 |
| 32 | 0.461 ± 0.061 | 7.3 / 3.9 | 0 | 0.714 | 7 / 2 | 0.586 ± 0.021 | 7 / 2.9 |
| 64 | 0.511 ± 0.041 | 15.9 / 7.7 | 0 | 0.857 | 14 / 2 | 0.636 ± 0.031 | 14 / 5.1 |
| 128 | 0.481 ± 0.037 | 33.0 / 17.0 | 0 | 0.897 | 29 / 3 | 0.691 ± 0.018 | 29 / 8.9 |

**Grid G3:**

| L | N0 Φ | N0 \|U\| / cls | N1 Φ | P1 Φ | P1 \|U\| / cls | P1-ε Φ | P1-ε \|U\| / cls |
|---|---|---|---|---|---|---|---|
| 27 | 0.446 ± 0.066 | 6.5 / 3.5 | 0 | 0.875 | 16 / 2 | 0.739 ± 0.034 | 13.2 / 3.5 |
| 81 | 0.524 ± 0.033 | 19.9 / 9.4 | 0 | **0.800** | 25 / 5 | **0.631 ± 0.025** | 22.1 / 8.2 |
| 243 | 0.416 ± 0.018 | 62.2 / 36.3 | 0 | 0.924 | 144 / 11 | 0.842 ± 0.005 | 121.2 / 19.1 |

## 2. Qualification (exact §R algebra)

**Gaps:**

| family | quantity | values by L | check |
|---|---|---|---|
| P0 | Δ_j | 0.417, 0.253, 0.346, 0.416 | adjacent allowances 0.176 / 0.102 / 0.078; 0.253 ≥ 0.241 ✓ |
| P0-ε | Δlow_j | 0.129, 0.042, 0.053, 0.156 | all > 0 |
| P0-ε | Δε_j | 0.317, 0.124, 0.124, 0.211 | 0.124 ≥ 0.317 − 0.270 ✓ |
| P1 | Δ_j | 0.429, **0.276**, 0.508 | required Δ₂ ≥ 0.429 − (0.066 + 0.033) = **0.330** ✗ |
| P1-ε | Δlow_j | 0.192, 0.049, 0.403 | all > 0 |
| P1-ε | Δε_j | 0.293, **0.106**, 0.426 | required Δε₂ ≥ 0.293 − 0.157 = **0.135** ✗ |

**Verdicts:**

| rule | P0 | P1 | P0-ε (Q6) | P1-ε (Q6) |
|---|---|---|---|---|
| Q1: above up_R at every L | **PASS** | **PASS** | — | — |
| Q2: above N1 at every L | **PASS** | **PASS** | — | — |
| Q3: adjacent non-shrinking | **PASS** | **FAIL** (27 → 81) | — | — |
| Q4: parent array only | PASS | PASS | — | — |
| Q5: defined on N0 / N1 | PASS | PASS | — | — |
| Q6.1: Δlow > 0 at every L | — | — | **PASS** | **PASS** |
| Q6.2: lo above N1 | — | — | **PASS** | **PASS** |
| Q6.3: adjacent non-shrinking | — | — | **PASS** | **FAIL** (27 → 81) |
| **Q6** | — | — | **PASS** | **FAIL** |

## 3. Terminal

**D0-B — REGISTERED TREE OBSERVABLES INSUFFICIENT — HIERARCHY NOT OPERATIONALLY IDENTIFIED.**

The primary fails Q1 – Q5 for the deterministic positive P1, on Q3. Per §R D0-O7, this is D0-B: neither D0-A-T nor D0-A-N
is available. The secondaries cannot rescue Φ, and nothing was varied after the results: no exponents, sizes, seeds or
rules.

**Scope (D0-O2).**
- This rejects **Φ, as chartered, under the preregistered finite-grid non-shrinking rule**.
- It does **not** show that hierarchy is undetectable or impossible.
- It does not reverse or modify the PM2 Stage-A terminal.

## 4. What the data show (descriptive; no grade changes)

1. **P0 family: complete pass, including noise tolerance.**
   - Φ separates the binary bisection hierarchy from N0 and N1 at every G2 size.
   - It does so with the fine scale randomised over about one-third of the log-mass range (P0-ε), and the noise-tolerant
     gap is non-shrinking over G2.
   - Taken alone, P0 would have met the D0-A-N conditions. **This is why D0 required two structurally different
     positives: P0 alone would have over-qualified the detector.**
2. **P1 family: separated everywhere, but not monotonically.**
   - P1 and P1-ε lie above the N0 interval at **every** size (Q1, Q2, Q6.1, Q6.2 pass). The gap is largest at L = 243.
   - The gap **dips at L = 81** by more than the registered sampling allowance.
   - At L = 81, P1 has 5 skeleton classes (2 at L = 27, 11 at L = 243), and N0's Φ peaks at 0.524 on this grid.
   - With base-3 recursion and fixed √N / √S scales, the unit threshold and the skeleton resolution fall at different
     positions relative to the recursion levels at each L. This is a discreteness effect of the kind the non-shrinking
     rule is designed to penalise.
   - Whether it would wash out at larger L is **not tested**. No sizes are added (D0-O5).
3. **The conjecture Φ(N0) → 0 is not supported on these grids.**
   - Uniform spanning trees show substantial coarse-skeleton recurrence: Φ = 0.25 – 0.52, with about half of the units
     sharing a skeleton class.
   - Only at L = 243 is a decline visible (0.416).
   - At the sizes reachable here, units are few (3 – 62 per tree) and their √S-resolution skeletons are small. Chance
     recurrence among random trees is therefore large. This is the main reason the separation margins are moderate rather
     than qualitative.
   - The charter's audit expectation for N0 ("CONJECTURE Φ → 0") is recorded as **not supported at tested sizes**.
4. **Conjecture "D(P0), D(P1) = O(log N)" is false at tested sizes.**
   - The exact class count grows as N^{0.70} (P0) and N^{0.53} (P1), versus N^{0.94} (N0) and N^{0.51} (N1). Block
     subtrees occur in many exit-position / shape variants.
   - Exact recurrence does **not** separate P1 from the SP comb N1 (both ≈ N^{1/2}).
   - It is destroyed by fine-scale noise: P0-ε gives 0.95 and P1-ε 0.92, both ≈ N0. This confirms the audit's
     "template detector, brittle" assessment of family 5a.
5. **Families 1 – 4 are offset class, as the audit predicted.**
   - P0 is not separated from N0 by λ, τ, η_H, Strahler, R_B or Tokunaga c (P0 Strahler = N0 Strahler).
   - The top-window light-fraction KS distances are small for P0 / P0-ε (0.05 – 0.14) against N1's 0.72 – 0.82.
   - P1 differs more (λ ≈ 0.21 vs 0.16; KS 0.17 – 0.24). These remain secondaries and cannot rescue the primary.

## 5. Information accounting

**Supplied:**
- the lattice and corner root;
- unit source mass;
- the control constructions (N0 UST and the ε-UST: 20 seeds each; N1 BFS; P0 bisection; P1 ternary nesting);
- the grids;
- the algebraic scales √N, √S and ⌊N^{1/3}⌋;
- the preregistered rules.

**Derived:** a measured, preregistered failure of Φ on the second positive family's non-shrinking criterion. Also
measured:
- Φ's noise-tolerant success on the first positive family;
- quantitative random-tree recurrence baselines;
- the falsification of two audit conjectures (Φ(N0) → 0 and D(P) = O(log N)) at the tested sizes.

**Not earned:** a qualified detector; any statement about PM2, GRUT physics, C3-A1 or C3-B. **TRUE COMPRESSION: 0.**

## 6. Consequence for C3 (as preregistered for D0-B)

C3 must **repair its operational definition of architecture before any further physical-model search**. PM2-BIS remains
unopened.

**What D0 established for that repair:**
- (i) Mass-law families (1 – 4) do not separate even a supplied hierarchy (P0) from random trees.
- (ii) Exact recurrence is a brittle template detector.
- (iii) Coarse-skeleton recurrence with fixed √N / √S scales recognises supplied hierarchies above the random baseline at
  every tested size and is noise-tolerant on P0. However, its gap is **not scale-stable across recursion bases** on the
  preregistered grid, and random trees carry a large chance-recurrence baseline at these sizes.

**No repair is proposed or applied here.** The choice of next step belongs to the owner.

**Stop.** No PM2-BIS, no C3-B, no PR, no merge.
