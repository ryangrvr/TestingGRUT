# NR-H — independent-code review of the S2-H2 priorities

**Files:** `review/nr_h2.py` (+ `.log`). No H2 code is used.

**Independent routines:**
- an own fixed-step RK4 (the original used `solve_ivp`);
- random-sequential zero-T Glauber (the original used checkerboard);
- exact rational arithmetic;
- fresh models and seeds.

**Independence level:** *independent code path, not independent reviewer.*

| ID | frozen claim | reproduction | grade |
|---|---|---|---|
| **NR-H1** | structural attractors: a symmetric contraction → 0; a doubly stochastic primitive chain → uniform from every start; uniformity is tuned | F = 0.7R(1.3)x: \|x₃₀₀\| ≤ 1.8·10⁻⁴¹ from every start, including 5·10⁵. A fresh 6-state Birkhoff chain → uniform to 1.7·10⁻¹⁶ from 11 starts. A 0.03 perturbation gives max \|π − 1/6\| = 0.018 | **REPRODUCED** |
| **NR-H2** | explicit-target laws relocate the state value into D | the affine x* equals (I − cR)⁻¹b to 10⁻⁶ for two fresh b. Gradient-flow endpoints from −8 / 0 / 8 equal the root of x³ + x − a for a = 0.3 and −1.1 (to 10⁻⁶) | **REPRODUCED** |
| **NR-H3** | all-to-all alignment synchronizes; a local ring keeps twisted basins | complete graph: r = 1.000000 for 20 / 20 starts. Ring: **19 / 30** twisted (original 21 / 30); windings −3 … 2 | **REPRODUCED** |
| **NR-H4** | 2D zero-T Ising freezes a finite fraction in straight stripes; the complete graph aligns every start mod Z₂ | 2D (L = 24, sequential): non-uniform **0.40**, all of them straight stripes (original 0.44 at L = 32, checkerboard; SKR-type, SECONDARY). Complete graph N = 201: 100 / 100 aligned, '+' 0.53 (gauge) | **REVIEW-CONFIRMED-WITH-SCOPE**: numerics cover random starts. The **every-state** statement rests on the majority argument for odd N: a minority spin sees h ≥ 2 toward the majority, a majority spin never sees h against it, so the majority only grows. This argument was checked by the reviewer, not simulated over all 2²⁰¹ states |
| **NR-H5** | bijective / finite-time / unitary / invertible-Markov laws preserve fine-grained information | an exact-rational affine contraction, 150 steps forward then back: **x₀ recovered exactly**. Gradient flow, T = 2, RK4: 2.50000000. Unitary dilation: global distinguishability **0.70710678 at every collision**. Invertible P³: p₀ recovered to 6.9·10⁻¹⁶ | **REPRODUCED** |

**S2-H2 numerical verdict: REVIEW-CONFIRMED** (H4 with the scope noted above).
