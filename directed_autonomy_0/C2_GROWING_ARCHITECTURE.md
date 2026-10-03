# DA0 · C2 — GROWING ENDOGENOUS ARCHITECTURE (result)

> **Repaired by C2 SCOPE REPAIR 01** (C2S-01, C2S-02; `DA0_CORRECTION_LEDGER.md`). The owner ruled that C2-A PASSES.
> Where wording differs, the ledger takes precedence.

**Preregistration.** C2 CHARTER + C2 CHARTER REPAIR 01 (C2R-01 … 10), with **E1 = East model** chosen and recorded
before any computation.

**Status of results:**
- Proofs: **INTERNALLY PROVED / NOT EXTERNALLY REVIEWED**.
- Numerics: `c2/c2_growth.py` + `c2/c2_growth.log` (independent code path, not independent reviewer; NUMERICAL
  ILLUSTRATION).
- Nothing is optimised (C2-F1).

## C2.1 Theory

**Notation.**
- A = span{1_{B_i}}, i = 1 … K, with p_i = π(B_i), Λ = p_min^{−1/2} and **κ = K·p_max ∈ [1, K]** (balance factor).
- V is the slow space, s = ‖E_V − E_A‖_π, and C_V = sup_{f ∈ V, ‖f‖ = 1} ‖f‖_∞ (≥ √K).
- ε_bd = s²(11Λ + 2C_V) is Lemma G5-2's tensor bound. It is already explicit in K: there are no hidden K factors.

### LEMMA G5-3′ (improved global exclusion; PROVED HERE)

**Statement.** In Lemma G5-3, the global threshold 1/(36K^{3/2}Λ²) can be replaced by

  **ε_glob′ = 1/(36·κ^{3/2}·Λ²)**, with κ = K·p_max.

**Where the old bound was crude.** G5-3(b) used λ_i ≥ 1 to bound idempotent norms. Using λ_min = p_max^{−1/2} instead:

**Proof.**
- **Norm bound.** For an idempotent x ≠ 0: 1/‖x‖ = ‖T(x̂, x̂)‖ ≥ (Σλ_i² x̂_i⁴)^{1/2} − ε ≥ λ_min/√K − ε = κ^{−1/2} − ε.
  So ‖x‖ ≤ 2√κ for ε ≤ 1/(2√κ).
- **Residual.** ‖r‖ = ‖E(x, x)‖ ≤ 4κε.
- **Per coordinate.** λ_i|x_i||x_i − 1/λ_i| = |r_i| gives t_i² ≤ |r_i|/λ_i ≤ |r_i|/λ_min. Hence
  ‖x − x_S‖² ≤ (√K/λ_min)‖r‖ = √κ‖r‖ ≤ 4κ^{3/2}ε.
- **Uniqueness.** This is < (1/(3Λ))², the Kantorovich uniqueness radius of G5-3(a), iff ε < 1/(36κ^{3/2}Λ²).
- **Remaining steps.** ε_glob′ ≤ 1/(8Λ) and ε_glob′ ≤ 1/(2√κ). The Weyl primitivity step needs (8Λ/3)ε + 2ε < 1/2, which
  holds since ε < 1/(36Λ²). ∎

**Consequence.** For balanced blocks (κ ≤ β), the threshold improves from **K^{−5/2}** to **K^{−1}**.

### THEOREM C2-A1a (K-explicit slow space → growing partition algebra; PROVED HERE)

Statement in terms of K, p_min, κ, C_V and s only (C2R-04).

**Hypothesis.** If ε_N := s_N²(11Λ_N + 2C_{V_N}) < 1/(36·κ_N^{3/2}·Λ_N²), then:

1. **Map defined.** R(V_N) is defined. The compressed product has exactly 2^K idempotents, and **exactly K_N primitive
   ones**. K_N is an output, not an input.
2. **Idempotent errors.** ‖ẽ_i − 1_{B_i}‖ ≤ d_i := (8/3)ε + √2·s·√p_i.
3. **Total misclassified mass.** m_tot ≤ 2Σd_i² ≤ (256/9)Kε² + 8s².
4. **Per-block recovery [C2R-05], proved here rather than inferred.**
   **m_rel ≤ ((256/9)Kε² + 8s²)/p_min.**
   - Every point of B_i Δ B̂_σ(i) is a misassigned point, and the misassigned set has mass ≤ 2Σd_l² (Chebyshev rounding).
   - So π(B_i Δ B̂_σ(i))/p_i ≤ m_tot/p_min.
   - σ is the G5 matching, chosen proof-side.
5. **Projector distance.** gap(V, A_Π) ≤ s + 2Λ√m_tot (when m_rel < 1, all blocks are non-empty).

**COROLLARY C2-A1a-bal.** Assume β-balance (p_max ≤ β p_min, hence κ ≤ β and Λ² ≤ βK) and C_V ≤ c√K. Then

  **s_N · K_N^{3/4} → 0**

suffices for 1 – 5 eventually, with m_rel ≤ (256/9)βc₁² s⁴K³ + 8βs²K → 0, where c₁ = 11√β + 2c.

**Status of the condition.**
- It is **sufficient, not shown necessary.**
- It **supersedes** the anticipated s·K^{3/2} → 0. That exponent was an artefact of the crude global bound.
- Whether K^{3/4} is itself intrinsic remains **open** (C2R-06).

### C2-A1b (microscopic family → accurate slow space): K-dependence of every constant

| step | constant | K-dependence |
|---|---|---|
| reversible, Davis–Kahan sin Θ | s ≤ η/(g − η) | **dimension-free** (no K) |
| Weyl cut certificate | ratio ≥ (g − η)/η | dimension-free |
| non-reversible, C1-T | e = ‖E₀‖, a = sup‖(z − G₀)^{−1}(I − E₀)‖ | E₀ and the fast resolvent are **block-diagonal**, so e and a are **maxima over blocks**; no K factor |
| tensor and recovery | Λ, κ, C_V | enter through Theorem C2-A1a only |

So the **dynamical step introduces no explicit K factor**. K enters through:
- the family's own η_N, g_N and per-block constants (e.g. η can grow with the number of neighbouring blocks);
- the **algebraic recovery step**.

**Sufficient, reversible, balanced, with C_V = O(√K):** (η_N/g_N)·K_N^{3/4} → 0. The assumption s = O(η/g) is **not**
hidden: for the reversible case it is Davis–Kahan, with constant 1.

### PROP C2-B1 (fixed-depth nesting; PROVED HERE) [C2R-07]

**Definition.** ν(coarse Ĉ, fine F̂) = Σ_f (π(f) − max_c π(f ∩ c)), the minimum π-mass to reassign.

**Statement.** Suppose the hidden partitions are nested (𝔅^(j+1) refines 𝔅^(j)). Then for blindly recovered partitions
at the two levels:

  **ν(Π^(j), Π^(j+1)) ≤ m_tot^(j) + m_tot^(j+1).**

**Proof.** Take a point x in a recovered fine block F̂ that lies outside Ĉ, the recovered coarse block matched to F̂'s
hidden parent. Then x is misassigned either at the fine level or at the coarse level. ∎

**Consequences.**
- At **fixed depth**, ν → 0 whenever each level satisfies Theorem C2-A1a (or G5 / C1-R at fixed K).
- With d_N → ∞, the total defect is ≤ Σ_j m_tot^(j), so **uniform-in-level** control is needed. That is not provided
  here.
- No post-processing, optimisation or choice of nesting is used. Exact numerical nesting is evidence only.

### CONJECTURE C2-B (class-restricted; unchanged in status) [C2R-10]

**Class:** finite-range local rates, bounded by a stated rate law, on a supplied lattice, with temperature / coupling
**held fixed** in N.

**Conjecture:** in this class, unbounded depth with every adjacent decay-rate ratio diverging does not occur without an
additional supplied scaling.

**Status:** **not proved.** The evidence below (P1, East) is consistent with it. It is **not** inferred from the product
bound.

## C2.2 Results (numerical illustration)

### C2-A2 — rate experiment

**Family FK.** Balanced, all-to-all; **ARCHITECTURE SUPPLIED BY FAMILY**. Coupling ζ_K = 0.3·K^{−α}, block size m = 6.
The columns compare the tensor bound with the old and new global thresholds, ε_glob (G5-3) and ε_glob′ (G5-3′).

| α | K = 4 … 64: s | s·K^{3/4} | s·K^{3/2} | cut ratio at K | primitive found | m_tot / m_rel (audit) | ε_bd < ε_glob (old)? | ε_bd < ε_glob′ (new)? |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.037 → 0.035 | 0.11 → 0.80 | 0.30 → 18 | 4.7 → 6.0 (**bounded: no certified cut**) | K / K at every K | **0 / 0** | no | no |
| 0.5 | 0.018 → 0.0044 | 0.05 → 0.10 | 0.15 → 2.3 | 8.9 → 42.9 | K / K | **0 / 0** | no | no |
| 1.0 | 0.0092 → 0.00055 | 0.026 → 0.013 | 0.073 → 0.28 | 17 → 339 | K / K | **0 / 0** | no (K ≥ 4) | **yes, at every K = 4 … 64** |

**Reading.**
- With the improved lemma, the α = 1 sequence is **rigorously certified at every tested finite size**. The certificate
  uses the proved upper bound ε_bd and an exactly computed s, both on the proof side. This is the first finite-size
  certified recovery in the programme. The old lemma certifies none of these sizes.
- Recovery is **exact** (m_rel = 0) even where no lemma applies (α = 0, 0.5), and even where no diverging cut exists
  (α = 0, the "algebraization without divergence" phenomenon seen in RA0 G4).

**Breakdown sweep A2b** (ζ up to 10; K = 8, 32, 64):
- Recovery stays exact up to ζ = 3 (s ≈ 0.32, cut ratio ≈ 1.2) **at every K**.
- It fails only at ζ = 10, where s ≈ 1 and the cut ratio is ≈ 1 (no metastability). There, at K = 32 and 64, R
  **refuses** (31/32 and 60/64 primitive found). At K = 8, R returns a partition with m_rel = 1.35.
- **No K-dependent breakdown was observed** in this family.

**GROWTH-RATE OBSTRUCTION: not found.** This is evidence **only**. FK is mean-field and balanced, so its perturbation is
spread evenly. A localised perturbation might behave differently. **The intrinsic threshold is open.**

### H1 — independent metastable bits

Dual grade **LARGE K, TRIVIAL ARCHITECTURE + ARCHITECTURE SUPPLIED BY FAMILY** [C2R-09].

| n (K = 2^n) | n·λ_s/g | cut ratio at K | slow eigenvalues / λ_s | max ratio inside slow sector | primitive | m_rel | macro Q vs Kronecker sum of single-bit generators |
|---|---|---|---|---|---|---|---|
| 1 (2) | 0.034 | 29.7 | 0, 1 | — | 2/2 | 0 | — |
| 2 (4) | 0.067 | 14.8 | 0, 1, 1, 2 | 2.0 | 4/4 | 0 | equal to 8e-17 |
| 3 (8) | 0.101 | 9.9 | 0, 1×3, 2×3, 3 | 2.0 | 8/8 | 0 | equal to 3e-16 |
| 4 (16) | 0.135 | 7.4 | 0, 1×4, 2×6, 3×4, 4 | 2.0 | 16/16 | 0 | equal to 3e-16 |

**Reading.**
- The whole 2^n slow sector separates only while n·λ_s/g is small. The cut ratio **falls** as n grows, as accounted for
  in C2R-09.
- The slow spectrum is k·λ_s, with internal ratios ≤ 2, so the **depth is 1**.
- The recovered macro-dynamics is **exactly a product of independent flips** (a G1-type factorisation). K is large; the
  architecture is trivial.

### H2 + H3 — driven ring of K basins; causal-state comparator

Dual grade **CLOCK + ARCHITECTURE SUPPLIED BY FAMILY**.

| K | decay-rate ratio at K | primitive | m_rel | macro cycle-space dim | ring current (uniform) | H3: distinct rows of uniformised macro kernel |
|---|---|---|---|---|---|---|
| 4 | 9.6 | 4/4 | 0 | **1** | 2.9e-3 | 4 |
| 8 | 8.6 | 8/8 | 0 | **1** | 1.6e-3 | 8 |
| 16 | 8.9 | 16/16 | 0 | **1** | 7.5e-4 | 16 |
| 32 | 8.6 | 32/32 | 0 | **1** | 3.9e-4 | 32 |

**Reading.**
- K grows, but the architecture is **one cycle**: a clock at any K.
- **H3 (comparator observation).** For the first-order Markov macro-process, distinct conditional future laws means
  distinct rows. So the K macrostates are the causal states. Predictive complexity grows with K, and adds nothing beyond
  the clock.

### P1 — nested ultrametric hierarchy

Grade **ARCHITECTURE SUPPLIED BY FAMILY**; b = 2.

**Results at ζ = 0.02:**

| depth d | states | cut ratios at ranks 2, 4, …, 2^d | blind recovery at every level | m_rel (audit) | nesting defect ν between consecutive levels (blind vs blind) |
|---|---|---|---|---|---|
| 1 | 8 | 59 | 2/2 | 0 | — |
| 2 | 16 | 28, 49 | yes | 0 | 0 |
| 3 | 32 | 53, 27, 33 | yes | 0 | 0, 0 |
| 4 | 64 | 50, 37, 27, 26 | yes | 0 | 0 |
| 5 | 128 | 50, 49, 36, 25, 20 | yes (up to 32/32) | 0 | 0 at all four nestings |

**Scaling with ζ (d = 3):** ratios (10.9, 5.7, 6.7) at ζ = 0.1 grow to (104, 54, 66) at ζ = 0.01, i.e. ∝ 1/ζ.

**Reading.**
- The recovered nested partition filtration is canonical and blind, with **no level chosen**.
- At **fixed ζ**, the level ratios stay **bounded as d grows** (and even decline at deep levels). So d_N → ∞ with every
  cut diverging requires **ζ_N → 0**: **DEPTH REQUIRES SUPPLIED SCALING**, on top of the supplied architecture.
- Exact nesting (ν = 0) is evidence; the guarantee is PROP C2-B1, at fixed depth.

### E1 — East model

Bounded description; no explicit hierarchy (C2-F3′ passed, C2 §6a).

**(a) Fixed q, L growing (no supplied scaling).**

| q | L = 4 | L = 6 | L = 8 | L = 10 | L = 12 |
|---|---|---|---|---|---|
| 0.1 (top ratio) | 10.0 | 6.7 | 5.1 | 4.1 | **3.4** |
| 0.3 (top ratio) | 3.4 | 2.3 | 1.8 | 1.6 | **1.3** |

The top adjacent ratios **decrease** with L, so there is **no diverging cut without scaling**. This is consistent with the
rigorous absence of separation at the equilibrium scale L ~ 1/q (Chleboun–Faggionato–Martinelli).

**(b) Fixed L, q → 0 (supplied scaling), q = 0.1 … 0.01.** The classification below is numerical only: ranks whose ratio
grows like q^{−a} with fitted a > 0.5. It is a diagnostic, not a certificate; the rigorous separation is the literature's,
for mesoscopic scales as q ↓ 0.

| L | diverging-cut ranks (fitted exponent a; ratio at q = 0.01) | depth d(L) |
|---|---|---|
| 3 | 2 (0.97; 76), 3 (1.01; 200) | 2 |
| 4 | 3 (0.97; 76), 5 (1.00; 100) | 2 |
| 5 | 2, 4, 8 | 3 |
| 6 | 3 (0.85), 6 (0.96), 13 (1.00) | 3 |
| 8 | 5, 13, 34 | 3 |
| 10 | 8, 28, 89 | 3 |

The depth grows only slowly with L in the tested range: **2 → 3 by L = 5, then 3 up to L = 10**. Every level exists only
because q → 0.

**Partition algebra in East's slow spaces (q = 0.01):**

| L | level (rank) | primitive found | Δ_HS | reading |
|---|---|---|---|---|
| 4 | 3 | 2/3 → **REFUSE** | 0.28 | |
| 4 | 5 | 4/5 → REFUSE | 3.7 | |
| 6 | 3 | **3/3**; block masses 0.0098, 0.0098, 0.98 | 0.32 | blocks = configurations with a **single isolated vacancy at site 5 or 6**, plus the vacancy-free bulk |
| 6 | 6 | 5/6 → REFUSE | 1.5 | |
| 6 | 13 | 6/13 → REFUSE | 46 | |
| 8 | 5 | **5/5**; masses ≈ q (four blocks), 0.96 | 0.51 | blocks = **single isolated vacancy at sites 5 – 8**, plus bulk |
| 8 | 13 | REFUSE | 5.9 | |
| 8 | 34 | REFUSE | 534 | |
| 10 | 8 | 7/8 → REFUSE | 2.2 | |
| 10 | 28 | REFUSE | 56 | |
| 10 | 89 | REFUSE | 6.1e3 | |

**Reading.**
- East's slow spaces are a **timescale filtration**, not a partition filtration.
- Only the first level is ever partition-like, and only at some L.
- Its recovered blocks are indexed by the **microscopic site of a stuck vacancy**. Which sites are slow is fixed by
  dynamics (the far sites), but the macrostate labels are supplied microscopic sites. Under **C2-F6** this is **not an
  emergent architecture**.
- Block masses are ~q → 0, so **G5's p_min regularity fails** in the very limit that creates the depth.
- Deeper levels have Δ_HS ≫ 1 and are refused.

## C2.3 Terminals (per route, per family; not collapsed)

| route / family | terminal |
|---|---|
| **C2-A**, FK (rate family) | **BLIND RECOVERY WITH K_N → ∞ PROVED** under Theorem C2-A1a (balanced: s·K^{3/4} → 0; certified at every tested size for α = 1), **+ ARCHITECTURE SUPPLIED BY FAMILY**. **GROWTH-RATE OBSTRUCTION: not found** (breakdown only at s ≈ 1, K-independent in FK); intrinsic threshold open |
| C2-A, H1 | **LARGE K, TRIVIAL ARCHITECTURE + ARCHITECTURE SUPPLIED BY FAMILY** (exact product of independent flips; depth 1) |
| C2-A, H2 | **CLOCK + ARCHITECTURE SUPPLIED BY FAMILY** (cycle-space dim 1 at every K) |
| C2-B, P1 | **FIXED-DEPTH NESTING** proved (PROP C2-B1); nested blind recovery to d = 5 numerically; **ARCHITECTURE SUPPLIED BY FAMILY + DEPTH REQUIRES SUPPLIED SCALING** |
| C2-B, **E1 East** | [C2S-02] **No endogenous growing partition architecture found; the tested spectral hierarchy requires q → 0, consistent with known East theory** (this is evidence at L ≤ 12, not a theorem that fixed-q East can never produce one). Original wording: **DEPTH REQUIRES SUPPLIED SCALING** (no diverging cut at fixed q; levels appear only as q → 0, depth ≤ 3 for L ≤ 10). **Timescale filtration without partition filtration.** The single partition-like level is a microscopic-site relabelling (C2-F6) with p_min → 0. **No EMERGENT ARCHITECTURE FROM A BOUNDED-DESCRIPTION RULE** |
| CONJECTURE C2-B | open; the evidence (P1, East) is consistent with it |

## C2.4 Comparator obligations (C2 §8)

| question | answer |
|---|---|
| 1. Product of independent bits? | H1: **yes, exactly** (Kronecker-sum test). FK and P1: no (all-to-all or ultrametric coupling) |
| 2. Just one cycle? | H2: **yes** at every K. FK / P1 / H1 are reversible, so they carry no macro current |
| 3. Beyond the causal-state description of the derived macro-process? | **No** for every family. The first-order Markov macro-processes have their macrostates as causal states (H3). K2 is not exceeded |
| 4. Architecture enumerated by the family? | FK, H1, H2, P1: **yes** (ARCHITECTURE SUPPLIED BY FAMILY). East: no, but East yields **no partition architecture** beyond a site relabelling, and needs q → 0 |

**K1 – K4:**
- K1 (bit / clock): matched by H1 and H2 at any K.
- K2: not exceeded.
- K3: no units or maximisation used.
- K4: extended to growing K by Theorem C2-A1a.

## C2.5 Information accounting

**[C2S-01] Blind recovery vs certification.** No K is supplied to the partition-recovery map once the canonical slow
space is certified: K = dim V. But **obtaining the correct slow space is certification, not blind discovery**.
- In the numerical controls, the known construction rank was used to extract the first K modes.
- Analytically, the family-level proof certifies that rank, and that proof may use the hidden architecture.
- This is the RA3-03 firewall: **blind recovery, conditional rank / cut certification.**

**Derived:**
- a **K-explicit blind recovery theorem**, Theorem C2-A1a with the improved Lemma G5-3′, replacing K^{3/2} by
  κ^{3/2}Λ² (balanced: s·K^{3/4});
- **per-block recovery** bounded explicitly;
- the K-free dynamical translation, C2-A1b;
- **fixed-depth nesting**, PROP C2-B1;
- **finite-size rigorous certification** in a growing-K sequence;
- blind recovery of nested filtrations in supplied hierarchies.

**Supplied:**
- every growing architecture that was recovered (FK, H1, H2, P1);
- the scalings (ζ → 0, q → 0) that create divergent level separations;
- reversibility, or A, for directed variants;
- balance and C_V regularity;
- the East lattice, locality and boundary.

**Not earned:**
- emergent unbounded architecture from a bounded-description rule;
- unbounded depth without supplied scaling;
- capacity in any information-theoretic sense;
- subsystem factorisation (H1's product structure is supplied);
- inside / outside;
- consciousness;
- **TRUE COMPRESSION relative to frozen GRUT: 0.**

## C2.6 Hard-stop answers

See `DA0_ZOOM_OUT_02.md`.
