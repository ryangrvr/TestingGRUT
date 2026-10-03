# DA0 · C3-D0 — ARCHITECTURE DETECTOR QUALIFICATION (charter and candidate audit only; NOT RUN)

**Status:** CHARTER / AUDIT ONLY. **No D0 control computation has been made.** No adaptive-flow simulation, no PM2-BIS,
no C3-B.

> **D0 CHARTER REPAIR 01 (§R at the end) supersedes the text below wherever they differ.**
> - **Owner ruling:** charter CONDITIONAL PASS. After Repair 01 is committed, the D0 control computation is approved to
>   run without another owner stop.
> - Superseded by §R: the "vanishing fraction" wording (§D0-F2, §D0-F2.2), the Q3 rule and its interpretation, the Q6
>   rule, and the D0-A grade labels.

**Owner ruling that opens D0:**
- PM2 deterministic Stage A is accepted. Terminal: **PM2-B — GENERATED TREE / NETWORK; REGISTERED HIERARCHY IS
  GENERIC / RANDOM-TREE-LIKE; NO DISTINCT GROWING COLLECTIVE ARCHITECTURE ESTABLISHED** (ledger PM2-O1 … O5).
- Stage-A boundary: `grut-directed-autonomy-0 @ 2139cc436e79782c231286f4437bdf62402140e0`.

**Preserved:**
- PM1 terminal;
- all previous boundaries;
- RA0-frozen `ab4fd860ebbc72a6a3c0c5dd1ff2cef52616b7f5`.

**Constraints:** no PR, no merge.

**Labels:**
- KNOWN (literature; metadata / abstracts only);
- PROVED HERE (internally; INTERNALLY PROVED / NOT EXTERNALLY REVIEWED);
- CONJECTURE;
- UNKNOWN (to be measured).

## Purpose

**Central question.**
> Can we define a canonical, scalable observable of growing collective hierarchy that distinguishes known supplied
> hierarchy from generic random-tree geometry, without clustering, fitted macro-units, a chosen hierarchy depth, or a
> post-hoc objective?

**Why this gate exists.** The PM2 Stage-A lesson is that *generating a tree is easy; proving that a generated tree
contains non-generic growing architecture is the hard part.* The registered PM2 suite (τ, η_H, Strahler) failed its own
supplied-hierarchy positive control: A2 ≈ R-TREE.

**What D0 is.** D0 certifies an **audit instrument**, not a physical law. A D0-A outcome is **not** a result about GRUT
physics.

**The ladder D0 serves:**

  adaptive variables (PM1) → generated topology (PM2) → **generated non-generic architecture** (not yet shown; needs a
  qualified detector)

## D0-F0 — Firewall

1. **PM2 Stage A is frozen.** Its outputs (logs, trees, statistics) are **not inputs** to D0. No D0 statistic, scale or
   threshold may be chosen with reference to them.
2. **Calibration uses controls only** (§D0-F2). PM2 trees are not recomputed or re-scored in D0.
3. **No retroactive change.** A qualified detector may later be applied to a newly chartered generator (e.g. PM2-BIS, or
   another model). It **never** changes the accepted PM2 Stage-A terminal.
4. **Primary-first.** The primary detector family (§D0-F4) is chosen **here, before any control is computed**.
   Secondary diagnostics are reported but **cannot rescue** a failed primary.

## D0-F1 — Inherited C0 discipline

**Forbidden in any detector:**
- a clustering objective;
- a selected number of levels;
- a fitted number of modules;
- hand-selected subtree boundaries;
- a post-hoc threshold chosen because it separates the controls.

**Allowed (supplied by the problem; priced):**
- the rooted tree (root = outlet);
- unit source mass per node, hence the subtree mass S_v = |T_v|;
- the size N;
- fixed **algebraic** scale conventions declared **here, before computation**, in the same family as the √N convention
  already registered in PM2R-04.

**Not used by the primary:** the lattice embedding (coordinates), the construction history, any block labels.

## D0-F2 — Control set (fixed now)

All trees are spanning trees of the L × L open square lattice, rooted at corner node 0 (index x·L + y), exactly as in
Stage A. **Every supplied positive is graded ARCHITECTURE SUPPLIED and never counts as a C3 construction result.**

| ID | control | definition | role |
|---|---|---|---|
| **N0** | R-TREE | uniform spanning tree, Wilson's algorithm rooted at 0; seeds 1 – 20 per L (the Stage-A code) | generic random-tree geometry |
| **N1** | SP-TREE | the registered BFS tree (neighbours visited in increasing index) | simple geometric routing |
| **P0** | A2 | the Stage-A recursive-bisection tree, `bisection_tree(L)`, used at **L = 2^k** only, so every split is exact | supplied hierarchy #1 (binary, balanced) |
| **P1** | **ternary nested block tree** (new; defined in §D0-F2.1) | deterministic, at **L = 3^k** only | supplied hierarchy #2 (9-ary nesting, non-binary and asymmetric block branching) |
| **P0-ε, P1-ε** | **fine-scale-randomised supplied hierarchies** (proposed; §D0-F2.2) | P0 / P1 recursion stopped at blocks of mass ≤ ⌊N^{1/3}⌋; inside each terminal block, a Wilson UST of that block rooted at the block's exit cell; seeds 1 – 20 | stress test: hierarchy above a vanishing fraction of the scale range, randomness below it |

**Size grids (fixed now):**
- **G2** = {16, 32, 64, 128} for P0, P0-ε, N0, N1;
- **G3** = {27, 81, 243} for P1, P1-ε, N0, N1. The largest is N = 59 049.

N0 and N1 are computed on both grids. There are no other sizes.

### D0-F2.1 — P1 (second supplied hierarchy), exact definition

`build(Q, r)` for a square block Q of side s = 3^j and exit cell r ∈ Q:
1. If s = 1, stop.
2. Partition Q into 3 × 3 sub-squares of side s/3, indexed q = 3a + b, where a is the x-block index and b the y-block
   index (mirroring idx = x·L + y). Let q_r be the sub-square containing r.
3. **Block-level tree:** BFS on the 3 × 3 sub-square grid from q_r, neighbours visited in increasing q. Each non-root
   sub-square's parent is the sub-square it was discovered from.
4. **Junctions:** for each non-root sub-square Q′ with block-parent P′, let e = the cell of Q′ at the **midpoint** of the
   side Q′ shares with P′ (unique, because s/3 is odd). Set parent[e] = the cell of P′ adjacent to e across that side.
5. Recurse: `build(Q′, e)` for each non-root sub-square, and `build(Q_{q_r}, r)` for the root sub-square.

The top-level call is `build(whole lattice, 0)`.

**Properties:**
- It is deterministic, with no random component.
- Depth k = log₃ L grows.
- Block-level branching is **non-binary**: 2, 3 or 4 children depending on where the exit lies.
- Block-level masses are **asymmetric**. Example: for a corner exit, the block tree is
  (0,0) → {(0,1) → {(0,2) → (1,2) → (2,2), (1,1) → (2,1)}, (1,0) → (2,0)}, with branch masses 6 : 2 in sub-square
  units.
- It is not a relabelled A2: 9-ary rather than binary nesting, a different junction geometry, and different branching
  numbers.

**Disclosure.** The P1 block-level rule is BFS, i.e. SP-TREE at the block level. P1's hierarchy comes **from the
recursion**, not from the routing rule. A detector that merely recognises "SP-like" structure would therefore separate
P1 from N0 for the wrong reason. **This is why N1 must also be rejected** (Q2).

### D0-F2.2 — Why the ε-controls are proposed (owner may strike)

Any physically generated tree will carry fine-scale randomness. A detector that recognises only **exact** deterministic
templates would therefore fire on P0 / P1 and on essentially no generated tree. It would be a determinism detector, not a
hierarchy detector, which defeats D0's purpose (condition 5: "meaningful for a generated tree").

P0-ε / P1-ε keep every supplied level above mass ⌊N^{1/3}⌋. That is a growing number of levels, about (2/3)·log N, while
the randomised part is a vanishing fraction of the scale range. Q6 below grades this. **The owner may strike Q6.** Its
result then only labels the qualification (§D0-T).

### D0-F2.3 — Exact control facts (PROVED HERE, elementary)

**LEMMA D0-N1 (SP-TREE is a comb).**
- With index x·L + y and sorted neighbour visits, BFS dequeues each layer x + y = d in increasing x (induction on d).
- Hence parent(x, y) = (x − 1, y) for x ≥ 1, and parent(0, y) = (0, y − 1).
- N1 is a **comb**: a spine along x = 0 with L teeth, each a straight path of L − 1 nodes.
- Subtree masses:
  - tooth node (x, y): S = L − x;
  - spine node (0, y): S = (L − y)·L.

**COROLLARY.**
- Every **light** child of N1 (§D0-F4) is a tooth top, with S = L − 1 < √N = L.
- Hence N1 has **no primary units**, and **Φ(N1) = 0 identically** (Q2 is analytic for N1; it is still computed as a code
  check).
- The exact subtree-isomorphism class count is **D(N1) = 2L − 1**: L path lengths plus L − 1 combs, i.e. Θ(√N).

**LEMMA D0-P0 (A2 subtrees are blocks).** In `bisection_tree`:
- Every non-root cell is chosen as `bcell` exactly once, because N − 1 splits choose cells in disjoint new halves, and a
  chosen cell is thereafter only an exit r.
- The subtree of a cell is **exactly** the block for which it is the exit: cells of B drain to `bcell` by recursion, and
  B is attached through `bcell` only.
- So P0's subtrees are precisely its supplied blocks. The analogous statement holds for P1: the subtree of a junction
  cell e = the sub-square Q′ ∪ its block-level descendants.

## D0-F3 — Candidate detector audit (before any computation)

Columns: invariance under root-preserving tree isomorphism; supplied scale or threshold; expected behaviour on N0 / N1 /
P0 / P1; risk that a random tree passes for purely combinatorial reasons.

### Family 1 — branch-mass / child-mass ratio statistics across scales

**Definition.** For each node v with S_v ≥ 2:
- the heavy child h(v) = the unique child of maximum mass;
- the light fraction q_v = 1 − S_{h(v)}/(S_v − 1), with q_v = 1 if the maximum is tied.

Consider the law ℒ_k of q over nodes in dyadic window W_k = [2^k, 2^{k+1}).

**Invariant / supplied:** root-isomorphism invariant. Supplied: dyadic windows (base 2).

**Expected:**
- N0: a scale-stationary, non-degenerate ℒ_k (CONJECTURE, from the statistical self-similarity of UST).
- N1: degenerate; teeth have q = 0, and spine q is a deterministic function of S.
- P0 / P1: **UNKNOWN.** Node-level splits inside blocks depend on internal junction paths; they are not simple atoms.

**Random-pass risk: HIGH.** Any statistically self-similar tree gives a scale-stationary law. Separation could only come
from the **shape** of the limiting law (offset class).

### Family 2 — scale-to-scale recurrence of rooted subtree structure

**Definition.** Compare a unit's structure at scale 2s with the composition of units at scale s.

**Invariant / supplied:** invariant; the scale step is supplied.

**Expected:**
- N0: recurs **in distribution** (self-similarity), so it passes in distribution.
- P0 / P1: recur **pathwise** (identical blocks).

**Random-pass risk: HIGH in the distributional form.** Only the **pathwise** (isomorphism) form separates, and that form
is Family 5.

### Family 3 — Horton / Tokunaga self-similarity beyond raw Strahler

**Definition.** Strahler branches; counts N_k; side-branch counts T_{i,j}; Horton ratio and Tokunaga (a, c).

**Invariant / supplied:** invariant; no supplied scale.

**Expected:**
- N0: Tokunaga-type. **KNOWN** for critical binary Galton–Watson trees (Burd, Waymire & Winn, *Bernoulli* 6 (2000));
  CONJECTURE for lattice UST.
- N1: Strahler 2, so it is degenerate.
- P0 / P1: UNKNOWN.

**Random-pass risk: HIGH.**
- Raw Strahler already failed in Stage A: A2 equals R-TREE.
- Self-similarity itself is generic. Separation can come only through the (a, c) values (offset class).

### Family 4 — canonical source-mass filtration

**Definition.** T_m = {edges with S_e ≥ m} for dyadic m (upward closed, so it is a rooted tree). Use the branching
transition law between T_{2m} and T_m: leaf counts and merge statistics.

**Invariant / supplied:** invariant; the dyadic m is supplied.

**Expected:**
- **Scale stationarity is generic:** for a self-similar tree, leaves(T_m) ~ N/m, and edges(T_m) ~ N m^{−τ}.
- This reduces to τ-type exponents and constants. The PM2 τ, a Family-4 quantity, separated A2 from R-TREE only weakly
  (0.33 vs 0.38).

**Random-pass risk: HIGH** (offset class).

### Family 5 — description / recurrence on rooted subtree isomorphism classes

**(5a) exact.** D(T) = the number of distinct rooted-isomorphism classes among {T_v}, i.e. the size of the minimal DAG;
canonical AHU labels.
- **Invariant / supplied:** invariant; no supplied scale.
- **Expected:**
  - N0: D close to linear. **KNOWN** for random binary trees: Θ(n/√log n) (Flajolet, Sipala & Steyaert, ICALP 1990);
    CONJECTURE for lattice UST.
  - N1: D = 2L − 1 = Θ(√N) (PROVED HERE).
  - P0 / P1: small. CONJECTURE: O(log N), from the block lemmas.
- **Random-pass risk: LOW.** But it is **brittle**: any fine-scale randomness makes every large subtree distinct, so it is a
  **template detector** (it fails P0-ε / P1-ε by construction).

**(5b) coarse.** Recurrence of **coarse branching skeletons** of canonical large units. This is the **PRIMARY** (§D0-F4).
- **Invariant / supplied:** invariant. Supplied: √N (units) and √S (resolution), both algebraic.
- **Expected:**
  - N0: **CONJECTURE** Φ → 0: coarse skeletons have ~N^{1/4} or more leaves and are random, so collisions vanish.
  - N1: Φ = 0 (PROVED).
  - P0 / P1: CONJECTURE Φ high, from block recurrence.
  - ε-controls: **UNKNOWN**.
- **Random-pass risk: LOW to MODERATE.** Small skeletons collide by chance at small N, which Q3 must see through.

### Rejected outright

**Embedding compactness** (subtree perimeter exponent):
- not root-isomorphism invariant;
- **N1's spine subtrees are rectangles**, so compactness would grade the SP comb as "hierarchical". Compactness is not
  hierarchy.

**Light depth** λ = ⟨# light edges to the root⟩ / log₂N: retained as secondary only. CONJECTURE: Θ(1) for both random and
hierarchical trees, so it is offset class.

### Audit finding (central; PROVED HERE at the level of the argument, not as a theorem)

Families 1 – 4 are functionals of the **mass / branching law**. A uniform spanning tree is itself **statistically
self-similar**, so these functionals become scale-stationary on N0 as well. They can separate a supplied hierarchy from
N0 only through **different limiting constants or law shapes** (the "offset class"). This is exactly the class in which
PM2's τ separation was shown to drift.

Only **pathwise recurrence** (Family 5) offers a qualitatively different asymptotic on N0, where recurrence vanishes, and
on supplied hierarchies, where it persists. **Exact** recurrence is a template detector. **Coarse** recurrence (5b) is the
one candidate that could, in principle, be both qualitatively separating and tolerant of fine-scale randomness.

**Hence the D0 question sharpens to a definitional fork for C3:**
- **"architecture = a distinct universality class of the branching law"** (offset class); or
- **"architecture = pathwise multi-scale recurrence of collective motifs"** (5b).

D0 tests the second as primary. A D0-B outcome would mean that C3's operational definition of architecture must be
repaired, as anticipated in the ruling.

## D0-F4 — PRIMARY detector (preregistered): coarse-skeleton recurrence Φ

Input: a rooted tree T (parent array only), root o, N = |T|, and S_v = |T_v|.

1. **Heavy child.** h(v) = the unique child of maximum S. If the maximum is tied, v has **no** heavy child.
2. **Light nodes.** A light node is a non-root node that is not the heavy child of its parent. This is the canonical
   heavy-light decomposition: light nodes are the tops of "streams".
3. **Units.** U = {v light : S_v ≥ ⌈√N⌉}.
4. **Coarse skeleton of v.**
   - μ_v = ⌈√S_v⌉.
   - Take the node set {w ∈ T_v : S_w ≥ μ_v}. It is upward closed in T_v, so it is a rooted subtree at v.
   - Apply **homeomorphic reduction**: repeatedly suppress non-root nodes with exactly one child.
   - The result K(v) is an unlabelled, unordered rooted tree.
5. **Class.** The AHU canonical string of K(v).
6. **Statistic.** Φ(T) = 1 − |{class(K(v)) : v ∈ U}| / |U|. Set Φ := 0 if |U| ≤ 1.

**Reading.**
- Φ is the fraction of large collective units whose coarse branching skeleton **recurs** elsewhere in the tree.
- Φ → 0: every large unit is coarse-unique (generic).
- Φ bounded away from 0: the large-scale structure is built from **recurring collective motifs**.

**Properties:**
- Canonical and root-isomorphism invariant.
- No embedding, no labels, no construction history.
- Defined for every rooted tree.
- Supplied scales: √N and √S, algebraic and declared now.
- **Cost:** O(Σ_{v∈U} |T_v|) plus hashing. Each node lies in at most log₂N light subtrees, so this is O(N log² N), which
  is feasible at N = 59 049.

**Implementation (to be written only after owner review):** a new `c3/d0_detector.py` that reuses the Stage-A tree code
(`lattice`, `wilson_ust`, `sp_tree`, `bisection_tree`, `subtree_stats`) **unchanged**. P1 and the ε-controls are new
functions implementing §D0-F2.1 / F2.2 verbatim.

## D0-F5 — Qualification rules (preregistered)

For F ∈ {P0 on G2, P1 on G3}:

Notation:
- m_R(L), h_R(L) = the N0 mean and the 95% t half-width over seeds 1 – 20;
- Δ(L) = Φ_F(L) − m_R(L).

| ID | rule |
|---|---|
| **Q1** (vs R-TREE) | at **every** L in F's grid: Φ_F > m_R + h_R |
| **Q2** (vs SP-TREE) | at every L: Φ_F > Φ_{N1} (= 0 by LEMMA D0-N1) |
| **Q3** (scale-structural, non-shrinking) | Δ(L_max) ≥ Δ(L_min), **and** Δ(L_max) ≥ Δ(L_{max−1}) − h_R(L_max). A separation that shrinks across the grid (the PM2-τ pattern) fails |
| **Q4** (no hidden labels) | the detector reads only the parent array: verified by code inspection and by applying it to N0 / N1 with the identical call |
| **Q5** (meaningful on non-family trees) | Φ is defined and computed for N0 and N1 with no special-casing beyond Φ := 0 for \|U\| ≤ 1 |
| **Q6** (proposed; noise tolerance) | P0-ε on G2 and P1-ε on G3 satisfy Q1 – Q3, using the lower 95% bound of their own 20-seed interval in place of Φ_F |

**Secondary diagnostics** (reported; **cannot rescue** a failed primary):
- 5a exact D(T) and its log-log slope;
- Tokunaga (a, c) estimates;
- the top-window law ℒ_k distance (P vs N0, two-sample KS);
- light depth λ;
- τ, η_H and Strahler on G2 / G3.

## D0-T — Terminals

| ID | terminal | requirement |
|---|---|---|
| **D0-A** | **CANONICAL HIERARCHY DETECTOR QUALIFIED** | Q1 – Q5 for **both** P0 and P1. Labelled **(noise-tolerant)** if Q6 also passes, and **(template-only)** if Q6 fails. A template-only D0-A does **not** by itself license applying Φ to a stochastic generator; the owner rules |
| **D0-B** | **REGISTERED TREE OBSERVABLES INSUFFICIENT — HIERARCHY NOT OPERATIONALLY IDENTIFIED** | the primary fails Q1 – Q5 for P0 or P1. C3 must then repair its operational definition of architecture (the definitional fork above) before any further physical-model search |
| **D0-INDETERMINATE** | stated technical obstruction only | e.g. the primary is infeasible at G3's largest N |

**Overfitting guard.** Success on P0 alone does not qualify. A detector that requires the construction history of P0 / P1
fails.

## Information accounting

**Supplied (priced):**
- the lattice and the corner root;
- unit source mass;
- the control constructions:
  - N0: UST, 20 seeds;
  - N1: BFS;
  - P0: bisection;
  - P1: ternary nesting;
  - ε-controls: recursion plus UST, 20 seeds;
- the size grids;
- the algebraic scales √N, √S and ⌊N^{1/3}⌋.

**Derived only if earned:** a qualified canonical detector of growing collective hierarchy.

**Not earned in any outcome:** any statement about PM2, GRUT physics, C3-A1, C3-B, consciousness or TRUE COMPRESSION
(which stays 0).

**Novelty:** none claimed. Isomorphism-class compression, heavy-light decomposition and Horton–Tokunaga analysis are
standard tools.

## Hard stop

Charter and audit only. **No D0 computation, no detector code, no PM2-BIS, no C3-B.** For owner review:
1. Accept P1 (§D0-F2.1) as the second supplied hierarchy?
2. Accept 5b coarse-skeleton recurrence Φ as the **primary** detector, with its √N / √S scales?
3. Keep or strike **Q6** (the ε-controls), and the noise-tolerant / template-only labelling of D0-A?
4. Accept the Q3 non-shrinking rule as the operational meaning of "asymptotic / scale-structural"?
5. Accept the size grids G2 = {16, 32, 64, 128} and G3 = {27, 81, 243}?

**Sources** (metadata / abstracts; full texts not re-read):
- Burd, Waymire & Winn, "A self-similar invariance of critical binary Galton–Watson trees", *Bernoulli* 6 (2000);
- Flajolet, Sipala & Steyaert, "Analytic variations on the common subexpression problem", ICALP 1990;
- Aho, Hopcroft & Ullman (1974): the rooted-tree isomorphism canonical form.


## §R — D0 CHARTER REPAIR 01 (owner review of `e52a276`; before any D0 computation)

### D0-O1 — P1 accepted

The deterministic nested 3 × 3 block tree (§D0-F2.1) is the second supplied hierarchy. Its BFS block-level routing does
not disqualify it:
- N1 controls directly for the BFS / comb confound;
- P1's supplied architectural content is the **growing recursive nesting**;
- P1 differs from P0 in recursion arity, branching multiplicity, mass ratios and junction geometry.

P1 is graded **ARCHITECTURE SUPPLIED**. It can never count as a physical C3 result.

### D0-O2 — Φ accepted as the primary; scales frozen; scope

The primary is **Φ** as in §D0-F4:
- units: light nodes with S_v ≥ ⌈√N⌉;
- skeleton resolution μ_v = ⌈√S_v⌉;
- canonical heavy-light decomposition (a tied maximum means no heavy child);
- homeomorphic reduction;
- AHU rooted-isomorphism classes.

**The exponents ½ (in √N and √S) are frozen** and will not be varied after results.

**Scope.** Φ is an **operational detector of pathwise recurrence of collective coarse branching motifs across scales**.
- "Φ qualified" does **not** mean all possible mathematical notions of hierarchy are identified.
- A Φ failure rejects **this operational definition**. It does not show that hierarchy is impossible or undetectable.

### D0-O3 — ε-controls and Q6 kept; scale-range wording corrected

**Correction.** The earlier statement that randomising below N^{1/3} affects a "vanishing fraction" of the scale range is
**false** on a logarithmic mass scale: log(N^{1/3}) / log N = 1/3.

**Correct statement.** P0-ε / P1-ε preserve a **growing majority** (asymptotically about **two-thirds**, modulo discrete
recursion levels) of the supplied recursive log-mass range, while replacing a **non-vanishing fine-scale portion** (about
**one-third** of the log-mass range) with random-tree structure.
- The randomised portion is not vanishing.
- This makes Q6 a **substantial** robustness test.

### D0-O4 — Q3 and Q6 repaired (exact algebra, fixed before computing)

**Notation.** For each family F and its ordered grid L_1 < … < L_m (G2: m = 4; G3: m = 3):
- **m_R(L), h_R(L):** the N0 sample mean and 95% t half-width over seeds 1 – 20, with h = t_{0.975, 19} · s / √20. A
  zero sample s.d. gives h = 0.
- **lo_X(L) / up_X(L):** m_X ∓ h_X for any 20-seed ensemble X.

**Deterministic positives (F = P0 on G2, F = P1 on G3):**

| ID | exact rule |
|---|---|
| **Q1** | for all j: Φ_F(L_j) > up_R(L_j) |
| **Q2** | for all j: Φ_F(L_j) > Φ_{N1}(L_j) |
| **Q3** | let Δ_j = Φ_F(L_j) − m_R(L_j). For every adjacent pair j = 1 … m − 1: **Δ_{j+1} ≥ Δ_j − [h_R(L_j) + h_R(L_{j+1})]** (and Q1 at every size) |
| **Q4** | the detector reads only the parent array (code inspection, plus the relabelling unit test, D0-O6) |
| **Q5** | Φ is defined and computed for N0 and N1 with no special-casing beyond Φ := 0 for \|U\| ≤ 1 |

**ε-positives (Fε = P0-ε on G2, Fε = P1-ε on G3; seeds 1 – 20):**

| ID | exact rule |
|---|---|
| **Q6.1** | for all j: **Δlow_j = lo_{Fε}(L_j) − up_R(L_j) > 0** |
| **Q6.2** | for all j: lo_{Fε}(L_j) > Φ_{N1}(L_j) |
| **Q6.3** | let Δε_j = m_{Fε}(L_j) − m_R(L_j). For every adjacent pair: **Δε_{j+1} ≥ Δε_j − [h_{Fε}(L_j) + h_{Fε}(L_{j+1}) + h_R(L_j) + h_R(L_{j+1})]** |

**Q6** passes for a family iff Q6.1 – Q6.3 all hold.

**Interpretation (renamed).** A pass of Q3 / Q6.3 means **SCALE-STRUCTURAL, NON-SHRINKING OVER THE PREREGISTERED FINITE
GRID**. It does not mean "asymptotically proved". A D0-A pass establishes empirical scale-structural separation over the
tested growing sequences, not a thermodynamic-limit theorem.

### D0-O5 — Grids fixed

Exactly G2 = {16, 32, 64, 128} and G3 = {27, 81, 243}. N0 and N1 are evaluated on both. No sizes are added after results.

### D0-O6 — Implementation firewall

**Before** any control value is computed:
1. Implement Φ exactly as chartered (`c3/d0_detector.py`), reusing the Stage-A tree code unchanged.
2. Write unit tests (`c3/d0_tests.py`) for:
   - (a) **relabelling invariance:** a random permutation of node labels leaves Φ unchanged;
   - (b) **rooted-isomorphism invariance:** children-order permutations leave Φ unchanged, and isomorphic trees built
     differently give equal Φ;
   - (c) **analytic N1:** Φ(N1) = 0 with \|U\| = 0, at all grid L (LEMMA D0-N1);
   - (d) **analytic hand cases:** a path gives Φ = 0; the complete binary tree of depth 6 gives Φ = 6/7. Derivation: 14
     units, at subtree depths 3 / 4 / 5 (8 / 4 / 2 of them). Their skeletons are a cherry, the depth-2 complete tree, and
     the depth-2 complete tree, so there are 2 classes;
   - (e) **construction validity** (not Φ): P1 and the ε-trees are spanning trees using only lattice edges.
3. **Do not inspect** P0 / P1 / N0 / ε Φ values while modifying the implementation, except to debug a violation of these
   exact tests. Any such inspection is recorded in the additive ledger before the final computation is rerun.

**Secondary diagnostics** (cannot rescue Φ), with implementation fixed here:
- exact class count D(T) and its log-log slope;
- Horton ratio R_B (geometric mean of N_k / N_{k+1}), with Tokunaga T_1, T_2 and c = T_2 / T_1;
- light depth λ = mean light-edge count to the root / log₂N;
- a two-sample KS statistic of the light fraction q over nodes with S ∈ [√N, N), for P vs pooled N0;
- τ, η_H and Strahler (Stage-A definitions).

### D0-O7 — Qualification grades (supersede the D0-A row of §D0-T)

| ID | grade | requirement |
|---|---|---|
| **D0-A-N** | **CANONICAL HIERARCHY DETECTOR QUALIFIED — NOISE-TOLERANT** | Q1 – Q5 for P0 **and** P1, **and** Q6 for P0-ε **and** P1-ε. **This is the only grade that, by itself, may support a later owner decision to apply Φ to a stochastic / physically generated tree family** |
| **D0-A-T** | **CANONICAL TEMPLATE-HIERARCHY DETECTOR QUALIFIED — TEMPLATE-ONLY** | Q1 – Q5 for P0 and P1, but Q6 fails for at least one ε-family. Φ recognises supplied deterministic recursive motifs. It is **not** shown fit for PM2-BIS or another noisy generator. This is an instrument result, **not a failure** |
| **D0-B** | **REGISTERED TREE OBSERVABLES INSUFFICIENT — HIERARCHY NOT OPERATIONALLY IDENTIFIED** | the primary fails Q1 – Q5 for P0 or P1 |
| **D0-INDETERMINATE** | stated technical obstruction only | |

### D0-O8 — Run scope

Run only N0, N1, P0, P1, P0-ε and P1-ε, at exactly the preregistered grids and seeds. **PM2 Stage-A trees are not
scored.** No PM2-BIS, no C3-B.
