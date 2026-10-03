# DA0 · C2 CHARTER — GROWING ENDOGENOUS ARCHITECTURE (preregistration; not started)

**Branch:** `grut-directed-autonomy-0`. C2 is part of the same programme as C1; DA0 is not frozen.

**Owner ruling.** C1 was accepted after C1 REPAIR 01, and **C2 is APPROVED TO CHARTER**. This file is that charter. **No
computation has been run.** C2 runs only after owner approval of this charter.

**Inherited prohibitions** (DA0 §0, C0):
- no consciousness interpretation;
- no quantum, collapse or Born-rule work;
- no frozen-branch or canonical-GRUT modification;
- no PR, no merge.

**Labels, as in DA0:** KNOWN / REDERIVED, PROVED HERE (internally, not externally reviewed), NUMERICAL ILLUSTRATION,
CONJECTURE, INTERPRETATION.

## 1. Why C2, stated from C1's boundary

At fixed K, C1's output is a **finite directed Markov macro-machine**: a K1-type loop over K2-type predictive states.
- K1 and K2 are undefeated.
- Making K larger by construction is trivial and teaches nothing.

**Central question.**

  Can dynamics derive an **unbounded endogenous state architecture** without supplying its capacity scale?

There are two routes. They are **not** equivalent, and they are tested separately:

| route | object | question |
|---|---|---|
| **C2-A** | flat growth: K_N = dim 𝒜_N → ∞ | Can blind recovery survive a growing number of macrostates, and is the growth forced rather than enumerated? |
| **C2-B** | depth growth: 𝒜_N^(1) ⊂ 𝒜_N^(2) ⊂ … ⊂ 𝒜_N^(d_N), with d_N → ∞ | Can the dynamics generate a nested filtration whose depth grows, without choosing the levels? |

A flat system with 10⁶ metastable states may have no hierarchy at all. A logarithmically deep filtration can support
multiscale organisation with small local branching.

## 2. Capacity firewall (C2-F; extends C0)

| ID | rule |
|---|---|
| **C2-F1** | **No capacity measure is chosen first.** No preregistered objective of the form "maximise entropy / mutual information / predictive information / statistical complexity / integrated information / number of states / depth". Nothing is optimised |
| **C2-F2** | **Outputs are structural and produced by the dynamics:** K_N = dim 𝒜_N; the certified cut ranks; the nested filtration and its depth d_N; branching ratios K^(j+1)/K^(j); the canonical macro current J_Q at each level and its cycle-space dimension. Any information-theoretic quantity may be computed **only afterwards**, as a description, and only alongside all four comparators |
| **C2-F3** | **Bounded-description rule (the owner's negative test, made operational).** A family whose definition **enumerates** its basins or levels (parameters that grow with K_N or d_N) has a **supplied architecture**. Its K_N → ∞ or d_N → ∞ is "forced by construction" and is graded **ARCHITECTURE SUPPLIED BY FAMILY**. A C2 positive of the strong kind needs at least one family given by an **N-independent finite rule** (plus N) whose architecture is not listed in its definition. This firewall is on **inputs**; it is not an objective and is never optimised |
| **C2-F4** | **No level chosen.** In C2-B the output is the whole filtration. No level is preferred, and no level index is supplied |
| **C2-F5** | **Cross-N discipline (inherited from RA2-02).** A growing rank sequence K_N is admissible only if each K_N is fixed by spectral data at that N, via a family-level certificate. No eigenvector or basin is matched across N |

## 3. The mathematical obstruction (recorded before any computation)

**G5's constants deteriorate with K.**
- The global idempotent threshold is ε_glob = 1/(36·K^{3/2}·Λ²), with Λ = p_min^{−1/2}.
- The tensor bound is ε_bd = s²(11Λ + 2C_V).

**Unavoidable floor.** Since p_min ≤ 1/K and Σ_x π_x Σ_i u_i(x)² = K:

  **Λ ≥ √K and C_V ≥ √K.**

**Anticipated corollary** (to be proved in C2-A1, not yet claimed). For balanced blocks (Λ = O(√K)) and C_V = O(√K),
G5's existing lemmas suffice for blind recovery whenever

  **s_N · K_N^{3/2} → 0**,

up to constants. So fixed-K recovery **does not** extend automatically to K_N → ∞.

**A further hypothesis fails.** H3′ (p_min ≥ p_* > 0) is **violated** for K_N → ∞. A growing-K theorem must carry p_min
explicitly.

**Recovery notion for growing K.** Vanishing *total* misclassified mass is too weak once blocks shrink. The preregistered
primary notion is **uniform per-block recovery**:

  max_i π(B_i Δ B̂_i)/π(B_i) → 0.

Total-mass recovery is reported as secondary.

## 4. C2-A — flat growth K_N → ∞

| item | target |
|---|---|
| **C2-A1** | **K-explicit blind recovery theorem.** State and prove sufficient conditions, in terms of (s_N, K_N, p_min,N, C_{V_N}), for G5 / C1-R to hold uniformly with per-block recovery. Expected first form: s·K^{3/2} → 0 under balance |
| **C2-A2** | **Rate question (potentially theorem-worthy).** *How fast may K_N grow relative to the metastable separation η_N/g_N before blind recovery fails?* (a) Can the global-exclusion step (Lemma G5-3(b), the source of K^{3/2}) be improved? The deterministic odeco perturbation literature, e.g. Auddy–Yuan, has gap-free character, so a K-uniform threshold may exist; audit before use, and do not upgrade algorithmic results. (b) **Necessity:** construct a family where s·K^α stays bounded and recovery provably or demonstrably fails, to test sharpness |
| **C2-A3** | **Computability, separated from canonicality** (RA3-03 / G5-R3 discipline). R's definition uses all 2^K idempotents, which is infeasible to enumerate for large K. Use the canonical set of **primitive** idempotents directly. Completeness of a found set is checked as a proved consequence: k = dim V distinct primitive idempotents ⇒ no others, once Lemma G5-3 applies. Not by exhausting starts |
| **C2-A4** | **Directed structure at growing K** (with C1). Report the canonical macro current J_Q, its cycle-space dimension and the C1R-02 relative-current condition with K-dependence (N4b may degrade with K) |

## 5. C2-B — depth growth d_N → ∞

| item | target |
|---|---|
| **C2-B1** | **Fixed-depth nesting theorem.** In RA0 G4.7 nesting was numerical only. Expected route: the Riesz ranges are exactly nested (V^(j) ⊂ V^(j+1)) for consecutive certified cuts of the same operator; each V^(j) is close to a hidden partition algebra; G5 recovery at each level then gives blind partitions nested up to vanishing mass. Prove it for fixed d, in the reversible and the C1 non-reversible classes |
| **C2-B2** | **Depth with d_N → ∞.** The constants accumulate across levels, so uniform-in-level conditions are needed. **Structural cost (to be stated precisely):** d levels, each with an adjacent decay-rate ratio ≥ R_N → ∞, require a total timescale spread ≥ R_N^{d_N}. Ask whether a **bounded-description** family (C2-F3) can supply this. **CONJECTURE C2-B (recorded, not assumed):** emergent unbounded depth with diverging ratios at every level is impossible for bounded-description local families without a supplied scaling of the separation (e.g. temperature or coupling scaled with N). Testing this conjecture is part of C2-B |
| **C2-B3** | **Directed hierarchy.** Report the macro currents per level, and whether coarser levels inherit or cancel finer currents. This is descriptive; no objective |

## 6. Controls (preregistered; each with a stated purpose; not to be multiplied)

| ID | control | purpose | expected if C2 is honest |
|---|---|---|---|
| **H1** | n independent metastable two-basin systems (product chain), K = 2^n | large K, structurally trivial | blind recovery may succeed. The recovered algebra is a **tensor product** of n bit-algebras (independent factors, G1-type), the slow spectrum k·ε gives **depth 1** and no further cut, and the macro dynamics are independent flips. Must be graded **LARGE K, TRIVIAL ARCHITECTURE** |
| **H2** | K_N metastable basins on a driven ring (R2 with K_N → ∞) | many macrostates, but a clock | blind recovery possibly positive. The macro cycle-space dimension stays 1, so the architecture is **one cycle**. Graded **CLOCK** (K1-equivalent at any K) |
| **H3** | an ε-machine with many causal states (a hidden-Markov output process; known theory) | predictive complexity with no endogenous physical boundary | K2 comparator. Its causal states are computed from a **supplied** observation process. Report it, and contrast it with C2 outputs, which come from D with no channel |
| **P1** | nested metastable hierarchy with fixed branching b and depth d_N → ∞ (ultrametric rates ε^j by level) | positive control for C2-B | nested blind recovery. Since the levels are enumerated by the family definition, it is graded **ARCHITECTURE SUPPLIED BY FAMILY** under C2-F3 |
| **E1** | **bounded-description candidate**: one N-independent local rule whose metastable architecture is *not* enumerated. Candidates for literature audit before selection: kinetically constrained models with hierarchical relaxation (East-type); Derrida-type hierarchical (GREM) landscapes under trap / Glauber dynamics; low-temperature local spin dynamics with many metastable configurations | the only control that can yield a **strong** C2 positive | outcome open. Any temperature or coupling scaling needed for diverging ratios is **priced** (C2-B2). Disorder, if used, is priced as a supplied random family with stated law |

Only H1, H2, H3, P1 and **one** E1 family are authorised. Adding more needs a stated purpose and owner approval.

## 7. C2 terminals (preregistered; per route, per family; not collapsed)

**C2-A:**
- **BLIND RECOVERY WITH K_N → ∞ PROVED** (state the growth condition, e.g. s·K^{3/2} → 0, and the recovery notion);
- **GROWTH-RATE OBSTRUCTION** (recovery fails beyond a stated rate, with a demonstrated or proved failure);
- **LARGE K, TRIVIAL ARCHITECTURE** (H1-type);
- **CLOCK** (H2-type).

**C2-B:**
- **CANONICAL DYNAMICAL FILTRATION OF UNBOUNDED DEPTH** (no level chosen);
- **FIXED-DEPTH NESTING ONLY**;
- **DEPTH REQUIRES SUPPLIED SCALING** (supports CONJECTURE C2-B).

**Either route:**
- **ARCHITECTURE SUPPLIED BY FAMILY** (C2-F3);
- **EMERGENT ARCHITECTURE FROM A BOUNDED-DESCRIPTION RULE** (strong positive; E1 only);
- **C2 INDETERMINATE** (stated technical obstruction).

**Information accounting** is required for every positive, as in DC1-04:
- **derived** — e.g. recovered K_N, the filtration, per-level currents;
- **supplied** — the family / rule, its scaling, reversibility or A, regularity, any disorder law;
- **not earned** — capacity in any information-theoretic sense unless separately shown; subsystem factorisation; inside
  / outside; consciousness; TRUE COMPRESSION.

## 8. Comparator obligations (C0 §2.3, at growing scale)

Every C2 positive carries a K1 – K4 table, plus explicit answers to:
1. Is the recovered architecture just a product of independent bits (H1)?
2. Is it just one cycle (H2)?
3. Does it exceed the causal-state description of its own derived macro-process (H3 / K2)?
4. Was the architecture enumerated by the family (C2-F3)?

## 9. Hard stop after C2

Answer:
1. Was blind recovery proved for K_N → ∞, and at what growth rate?
2. Is the K^{3/2} obstruction intrinsic or an artefact of Lemma G5-3?
3. Was fixed-depth nesting proved?
4. Did any family give d_N → ∞ without supplied scaling?
5. Did any architecture arise from a bounded-description rule?
6. How do H1 / H2 / H3 compare with the positives?
7. Was any capacity measure, level, K, ε or objective supplied?
8. Is C3 warranted, and with what starting object?

## 10. Deliverables (when C2 runs)

- `directed_autonomy_0/C2_GROWING_ARCHITECTURE.md`
- `directed_autonomy_0/c2/`
- `directed_autonomy_0/DA0_ZOOM_OUT_02.md`
- updates to `DA0_CORRECTION_LEDGER.md` and `DA0_STATUS.md`
