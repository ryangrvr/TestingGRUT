> **REPAIR 03:** the arrow / H_env statements below are scoped to the tested collision and channel dilations. Whether
> every effective arrow in a closed unitary universe needs a special state is OPEN (S2-D-arrow). See ZOOM_OUT_05,
> REPAIR 03 and Y-10.

# S2-H2 RESULT — can dynamics select the universe's state / basin?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-H2. Pre-registered at `90a133c`, including the admissible state spaces.
No admissible state was removed after the run.

**Files:** `probes/S2-H2/s2_h2_attractors.py`, with log `probes/S2-H2/s2_h2_attractors.log`.

**Wall attacked:** H.

**Notions kept separate:**
- **A** attractor uniqueness;
- **B** global reachability, from EVERY admissible state;
- **C** information erasure, judged from the exact final microstate including every degree of freedom in the model.

## 0. Headline

> **Dynamics can make the initial basin irrelevant without a measure. It does so in three ways, each with its own
> price, and fine-grained information is never erased unless the law is non-injective.**
>
> 1. **Explicit parameter** (contraction with offset b, potential with tilt a, a generic Markov chain). The attractor
>    *value* is written into the law. → **H → D RELOCATION**.
> 2. **Symmetry / class** (rotation-equivariant contraction, symmetric potential, doubly stochastic chain, all-to-all
>    alignment). The attractor is fixed structurally. → **TRUE H COMPRESSION INTO D**, but only while the symmetry is
>    exact (**TUNING-PRICED**). For the frame field, it holds only with **non-local, all-to-all** coupling.
> 3. **Dissipation.** Its unitary dilation (H2-1′) shows the attractor value equals the **prepared state of an
>    environment**, and global distinguishability is conserved exactly. → H is **relocated into the environment's
>    preparation**, not eliminated.
>
> **With local coupling (ring, 2D lattice), basin data survives.** Twisted states (21/30 starts) and frozen straight
> stripes (44% of starts) are stable alternatives.

## 1. Model by model

### H2-1 Deterministic contraction F(x) = 0.8 R(φ)x + b on ℝ²

- **A, B hold for every start**, including |x₀| = 10⁶. The max distance from x* after 400 steps is ≤ 1.1·10⁻¹⁶. No
  measure and no exclusions are needed (Banach).
- **The value of x\***:
  - with b = (1.3, −0.4), x* = (1.4009, 0.9507), set by b → **H → D RELOCATION**;
  - with b = 0, x* = 0, forced by rotation equivariance → **compression by symmetry**. A generic perturbation b ≠ 0
    moves x*, so the compression is **TUNING-PRICED**.
- **C (erasure): NO.**
  - The affine contraction is a bijection. In exact rational arithmetic the separation after 200 steps is 4.6·10⁻²⁶
    (from 10⁻⁶), and exact inversion recovers x₀ (True).
  - In float64 the recovery error grows from 2.6·10⁻¹⁴ (n = 20) to 1.1·10⁻⁴ (n = 120) to 1.1·10⁴ (n = 200).
  - Forgetting is **resolution- (A-) priced**.
- **Genuine erasure** needs a **non-injective** law. F(x) = 0.8|x| + 0.5 merges ±0.37 into 0.796 in one step.

### H2-1′ The contraction as a unitary dilation

A qubit system collides with fresh |0⟩ ancillas through a partial swap (θ = 0.6).

| collision k | 0 | 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| system TD | 0.7071 | 0.4121 | 0.2558 | 0.1657 | 0.1101 |
| global TD | 0.707107 | 0.707107 | 0.707107 | 0.707107 | 0.707107 |

- The system converges to the **ancilla state**.
- The attractor value is the bath's **preparation**: H moves into the environment.
- The information moves into the ancillas, and global distinguishability is conserved exactly.

### H2-2 Gradient flows (9 starts, from −50 to 50, including the stationary point 0)

| potential | outcome |
|---|---|
| x⁴/4 + x²/2 − 0.7x | every start → 0.5414. **H-SELECTED-FROM-D, H → D RELOCATION** (via a) |
| a = 0 | every start → 0. **Symmetry compression, TUNING-PRICED** |
| double well, tilt 0.05 | → −1.0241 or +0.9740 by basin. With the tilt, x₀ = 0 is no longer stationary: the unstable point moves to ≈ +0.05, so 0 flows to −1.0241. **BASIN DATA SURVIVES**; bistability is structurally stable |
| symmetric double well | ±1, gauge-related, **but x₀ = 0 stays at 0 forever**. **H-SELECTED-MOD-GAUGE only on the complement of the unstable point** (exclusion priced), and only at the tuned ε = 0 |

- **C:** the finite-time flow is a diffeomorphism. Backward integration recovers x₀ = 3.000000 (T = 2) and 3.000009
  (T = 8), but not at T = 20 (5·10⁵). → **asymptotic / resolution-priced forgetting.**

### H2-3 Primitive Markov chains and a primitive qubit channel

- **Doubly stochastic (Birkhoff mixture):** primitive, with π = uniform (0.2 × 5) reached from every initial
  distribution.
  - The uniform state is forced by double stochasticity, with no state-valued parameter: **TRUE COMPRESSION INTO D**.
  - **But the attractor is itself a measure** (counting measure). This parallels S2-3b: H_measure → D.
  - A 0.02 perturbation keeps uniqueness (primitivity is open) but breaks uniformity: π = (0.193, 0.198, 0.207, 0.199,
    0.204). → **TUNING-PRICED**.
- **Generic chain:** π = (0.233, 0.228, 0.153, 0.250, 0.136), encoded in the transition entries. → **RELOCATION**.
- **C:** P^t is injective (det P ≠ 0), so the exact distribution still encodes p₀. Recovery error grows 10⁻¹⁴ → 10⁻⁷
  → O(1) as cond grows to 10³³. → **resolution-priced.**
- **Probability is built into the transition law:** **H SELECTED, PROBABILITY D-PRICED**. This does not derive
  probability.
- **Amplitude-damping + dephasing channel:** TD 0.7071 → 0.0168 → 0.0004 → 0. The unique fixed point |0⟩⟨0| is the
  zero-temperature **bath state** (Stinespring = H2-1′).

### H2-4 Unitary hostile (random H, 6 qubits)

- |eigenvalues of U| = 1.000000000000. There is no attracting fixed point.
- Global TD = 0.992178987987 at t = 0, 1, 5, 20, 100, constant to 12 digits.
- Qubit-0 TD only fluctuates: 0.06 – 0.16.
- → **NO FINE-GRAINED H SELECTION** (scope: unitary / Hamiltonian / isometric). Cross-reference S2-3b.

### H2-7 Shared reference frames (the S2-5 / S2-6 follow-up)

| dynamics | result | classification |
|---|---|---|
| U(1) alignment, **complete graph**, identical frames | 20/20 random starts → r = 1.000000. The antipodal 23/17 configuration is **exactly stationary** (velocity ≤ 7·10⁻¹⁷, the float sin π) and unstable (Jacobian +1.026). Without a kick it does not synchronize in exact arithmetic | **H-SELECTED-MOD-GAUGE (U(1))** on the complement of a closed, nowhere-dense, Lebesgue-null set of unstable equilibria. Here Baire and Lebesgue genericity **agree** (unlike S2-3), but "every state" fails, so the exclusion is priced |
| same, with frequency disorder σ | r = 1.000 / 0.985 / 0.753 / 0.274 for σ = 0 / 0.2 / 0.5 / 1.0 | order becomes a **regime of D** (K > K_c): **TUNING- / REGIME-PRICED** |
| U(1) alignment, **ring** (local) | final windings over 30 starts: {−3:1, −2:7, −1:8, 0:9, 1:4, 2:1}. 21/30 end twisted | **BASIN DATA SURVIVES** (topological sectors) |
| Z₂ frame field (J ↔ −J), **2D lattice**, zero-T Glauber (pure energy descent) | 56% reach uniform alignment; **44% freeze in straight stripes**. All 53 remain non-uniform after 3000 more sweeps | **BASIN DATA SURVIVES** (SECONDARY: Spirin–Krapivsky–Redner frozen stripes) |
| Z₂ frame field, **complete graph**, zero-T, N = 201 | **200/200** reach uniform alignment; orientation + in 0.49 (gauge) | **H-SELECTED-MOD-GAUGE for every admissible state (odd N): TRUE H COMPRESSION INTO D**, priced by the all-to-all graph, the dissipative arrow and finite N |
| 2D, finite T (stationary measure of a primitive chain) | T = 1.8: ⟨\|m\|⟩ = 0.959, 2⟨f_A f_B⟩ = 1.896 at distance 16. T = 3.2: 0.068, −0.013 | needs T < T_c **supplied**. Per the charter this is **not counted** as dynamical selection: **REGIME-PRICED** |

**Shared-frame verdict.** Selection by D of the *existence* of a shared reference (orientation gauge) is achieved
**only** for all-to-all dissipative alignment. With local coupling, the kind S2-1b / S2-Σ locality criteria favour,
defects (twists, stripes) keep basin data alive.

→ **Shared-frame existence selected by D mod gauge: YES, conditional on non-local coupling + dissipation. NO for
local coupling.**

## 2. H2-8 Information accounting

| model | metric / topology | arrow | noise | exclusions | measure | access | final state unique (mod gauge)? | fine-grained erased? | verdict |
|---|---|---|---|---|---|---|---|---|---|
| contraction, b ≠ 0 | ℝ², Euclidean | dissipative | none | none | none | none | yes | **no** (bijection) | H-SELECTED-FROM-D; **H → D RELOCATION** |
| contraction, b = 0 | ℝ² with symmetry centre | dissipative | none | none | none | none | yes | no | TRUE COMPRESSION (symmetry); TUNING-PRICED |
| \|x\|-contraction | ℝ | dissipative | none | none | none | none | yes | **yes** (non-injective) | erasure only by a non-injective law |
| dilation | qubit + ancillas | unitary + fresh bath | none | none | none | system only | system → bath state | **no** (global TD const) | **H → H_env RELOCATION** |
| gradient, unique min | ℝ | dissipative | none | none | none | none | yes | no at finite T | RELOCATION (a) / symmetry compression (a = 0, tuned) |
| double well | ℝ | dissipative | none | none | none | none | **no** | no | BASIN DATA SURVIVES |
| symmetric well | ℝ | dissipative | none | **x₀ = 0** | none | none | mod gauge, except x₀ = 0 | no | H-SELECTED-MOD-GAUGE (exclusion priced); TUNING-PRICED |
| doubly stochastic Markov | 5 states | stochastic | **supplied** | none | counting measure = attractor | none | yes (uniform) | no at distribution level | TRUE COMPRESSION of a measure; PROBABILITY D-PRICED; TUNING-PRICED |
| generic Markov | 5 states | stochastic | supplied | none | none | none | yes | no | RELOCATION; PROBABILITY D-PRICED |
| qubit channel | Bloch ball | dissipative | — | none | none | none | yes (bath state) | no (Stinespring) | H → H_env RELOCATION |
| unitary | ℂ⁶⁴ | reversible | none | — | — | reduced only | **no** | **no** | NO FINE-GRAINED H SELECTION |
| Kuramoto, complete | T⁴⁰ | dissipative | none | unstable equilibria (closed, null, nowhere dense) | none needed | none | mod U(1) | no at finite T | H-SELECTED-MOD-GAUGE (exclusion priced); disorder → REGIME-PRICED |
| Kuramoto, ring | T⁴⁰ | dissipative | none | — | — | — | **no** (twists) | no | BASIN DATA SURVIVES |
| Ising zero-T, 2D | {±1}^1024 | dissipative (T = 0 bath) | tie-break | — | — | — | **no** (stripes) | yes (many-to-one) | BASIN DATA SURVIVES |
| Ising zero-T, complete | {±1}^201 | dissipative | none (odd N) | none | none | none | **mod Z₂, every state** | yes (many-to-one) | **TRUE H COMPRESSION INTO D mod gauge**; priced: non-local graph + arrow |

## 3. H2-9 Central question

**Did dynamics select the universe's state, or did we buy it with a contractive dissipative law that already contains
the answer?** Both, in different places.

- **Explicit-parameter laws** (b, a, generic π): the state value is written into D. → **H → D RELOCATION**, not
  elimination.
- **Structural laws** (symmetric contraction, doubly stochastic, all-to-all alignment): the attractor's identity
  follows from symmetry and class, with no state-valued parameter. → **TRUE H COMPRESSION INTO D**. This is the
  campaign's first one. Every instance is priced, by:
  - (i) exact symmetry (perturbations break it: TUNING-PRICED); or
  - (ii) non-local all-to-all coupling, which conflicts with the locality that Σ-criteria favour; and
  - (iii) dissipation.
- **Dissipation itself.** In every dilated case the attractor equals a prepared environment state, and global
  fine-grained information is conserved. The arrow that makes compression possible is paid for by an H-item: a
  low-entropy environment preparation.

**Status: S2-H2 COMPLETE.**
- H is **compressible into D** in structural cases (symmetric / doubly stochastic / all-to-all alignment, mod gauge).
- Each compression is **TUNING-, NON-LOCALITY- or ARROW-priced**.
- Explicit-parameter cases are **RELOCATION**.
- Dissipation **relocates H into the environment's preparation**.
- Unitary dynamics gives **NO FINE-GRAINED H SELECTION**.
- Fine-grained erasure occurs only for non-injective laws.
