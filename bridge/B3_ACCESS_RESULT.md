# B3 — A_res AGAINST CANONICAL GRUT ACCESS (result)

**Question (owner).** Once D_dyn and Σ are given, which components of A_res does GRUT genuinely derive, and which remain
independently supplied? "Access" is not treated as one object.

**Script:** `bridge/b3/b3_access.py` (+ `b3_access.log`), seed 20261003. It is an independent numpy re-implementation:
- the canonical `calc/p5_access.py` and `calc/p6_seed_selection.py` were **inspected for their closure rules only**, not
  run or imported;
- the label is **independent code path, not independent reviewer**.

**Classes used.**
- **Classical, earned Level-0 class:** K = 0.3·I + weighted Laplacian, ẋ = −Kx, linear readouts. Here the "closure" of a
  seed c is its observable (Krylov) subspace.
- **Auxiliary supplied quantum class**, used for the canonical P-5 / P-6 objects only. Per EA-0 these are **not defined**
  on the earned substrate without a lift. No auxiliary-class result is counted as Level-0 compression.

## B3-0 Canonical access decomposition (final owner-ruling scopes)

| component | canonical source and final scope | reviewed SCOUT A_res item |
|---|---|---|
| **A_seed** (what is coupled) | P-5 §2: "the seed itself remains a declaration". P-6: "access seed = NON-DERIVED INPUT", selected only up to block decomposition. Ledger A-13: "access seed … supplied" | what is coupled / accessible |
| **A_closure** (algebra generated from the seed under D) | P-5 §2: "canonical once the seed is declared" (standard bicommutant structure, NULL-REDUNDANT as a principle) | effect structure (partly) |
| **A_readout** (the actual exposed map) | L0 bridge ruling 02 §4: the earned end-site identity readout h = x₁ is globally observable (TRIVIAL/IDENTITY). §5: do not generalize. A-13: "declared readout" | effect / readout structure |
| **A_partition** (which sites / fragments count individually) | GS1: geometry is CLASS-SPLIT by access; "the access boundary itself is supplied". A-13: multi-site access sets are "declared … not earned" | fragment / locus grouping |
| **A_resolution** | P-6 L-D1: sufficiency is order-relative. GS1: the topology horizon ("invisible below the circumference order") | resolution / threshold |
| **A_time** | P-5: access is physical "exactly at boundary changes" (quench at t* = 3.0, declared) | time sequence / sampling, horizon |
| **A_coarse** | Layer 5 "environment / partition / coarse-graining" (SYN-0 §2); no canonical derivation | coarse-graining |
| (lift-dependent) | EA-0 = UNFORMULABLE at earned Level-0 scope; L0 bridge = NONUNIQUE-LIFT | — (category firewall) |

## Results

### B3-1 Closure given seed — **CONDITIONAL COMPRESSION FROM D + A_seed (+ a closure-rule convention)**

**Classical (earned class).** The closure is computed independently by Krylov iteration and by spectral projection. The
two agree to 1e-14, and the result is frame-covariant to 1e-14.

| case | closure dimension |
|---|---|
| generic net, generic seed | 24 |
| C1-a, earned end-site seed e₁ | 24 |
| C1-a, reflection-symmetric seed e₁ + e₂₄ | **12** (the even block) |
| decoupled abelian control | **2** |

**Auxiliary quantum class.** This reproduces P-5's dimension ladder of 16 / 8 / 4 (generic / parity-symmetric / abelian)
for both rules. Each rule is reproducible to 1e-15, and the ad_H rule is representation-covariant to 1.5e-15.

**New bridge finding:** the canonical record uses more than one closure rule.

| rule | where it is used | definition |
|---|---|---|
| **R_P5** | P-5 | alg{I, H, B} |
| **R_P6** | P-6 | unital algebra of the seed, ad_H-closed |
| **R_static** | P-6 L-D2 sufficiency test | ad_H step switched off |

They agree on the generic, parity and abelian cases. They **disagree when H carries structure outside the seed's reach.**
For site 1 decoupled from a hidden 2-spin sector:
- R_P5 gives dimension **16** with a 4-dimensional center, because it imports the hidden sector's energy projections into
  "access".
- R_P6 gives dimension **4**, the site-1 algebra, matching P-6 L-S.

So the compression is real but needs a declared closure rule.

> **A_closure = f(D, A_seed; closure rule)**: CONDITIONAL COMPRESSION. This is **not** "access derived".

**Redundant supply?** Ledger A-13 books "access seed / declared readout / access sets". In the canonical sources, the
"access sets" are **declared multi-site site sets** for geometry recovery (`L0_ACCESS_BRIDGE_CORRECTIONS_01`). Those are
A_partition / seed objects, **not closure algebras**.
- The ledger does **not** separately supply the closure, so **no REDUNDANT SUPPLY was found**.
- **CANONICAL UPDATE CANDIDATE B3-CUC-1** (bookkeeping clarification only): state in A-13 that the closure / effect
  algebra is downstream of (D, seed), and fix one closure rule. R_P6 is the natural candidate, since it does not count H's
  own hidden spectral projections as accessible. That choice is an owner decision.

### B3-2 Can D + Σ select A_seed? — **A_seed NOT REDUCIBLE TO D + Σ**

**Classical, earned class.** Same K, same net, single-site seeds.
- **C1-a path, N = 24:**
  - 16/24 sites are cyclic (informationally complete), and all of them have the **identical closure** ℝ^N with a trivial
    center.
  - Of the 120 pairs, 8 are transfer-function-identical. These are reflection pairs: **GAUGE**.
  - **112 pairs are physically distinct.**
- **Asymmetric weighted path:** 24/24 sites are cyclic, and **all 276 pairs are distinct**. This includes the two end
  sites e₁ vs e₂₄ (max |ΔG| = 0.143), which have the **same graph role in Σ**.

**Auxiliary quantum class.** This re-implements P-6 L-B with new Hamiltonians.
- Setup: a 3-spin chain, parity P = z₁z₂z₃ conserved, seeds σ_z¹ vs σ_z³.
- Over three random trials:
  - Both closures are exactly the 32-dimensional parity commutant (gap 1e-15), with center dimension 2 / 2. This holds
    under R_P6 **and** R_P5.
  - The probe coherence differs by **0.673 / 0.516 / 0.555**. The canonical value was 0.287 for its own H.
- **Σ does not resolve the pair.** The net is identical and both seeds are end sites. Any tie-break, such as "the end
  with the larger field", is an added criterion.

### B3-3 Readout completeness ≠ readout selection — **A_readout REMAINS SUPPLIED**

On C1-a (N = 24):
- 16 single-site readouts are informationally complete. The 8 that are not (0-based sites 1, 4, 7, …, 22) sit on nodes of
  an eigenvector.
- 200/200 random linear functionals are complete.
- e₁ vs e₂₄ is a gauge pair (ΔG = 5e-16). e₁ vs e₆ differ physically (ΔG = 0.51).

So **READOUT COMPLETENESS ≠ READOUT SELECTION.** State reconstruction *given* h stays an earned conditional result
(TRIVIAL/IDENTITY, L0 bridge §4). The readout itself is supplied.

### B3-4 GS1 inverse — **GS1 IS ONE-WAY: A → geometry, not geometry → A**

On one hidden weighted 3×3 grid:
- Full access reconstructs K to 1e-15.
- Every tested declaration is admissible on the same K: full, boundary {0, 2, 8}, another 3-set {1, 3, 5}, and the single
  sites {0} and {4}.
  - Their static response traces differ: 6.58 / 2.44 / 2.07 / 0.81 / 0.60.
  - **All 9 singletons and all 84 three-sets are informationally complete.**
- No candidate inverse selector returns a declared boundary:
  - max-information returns the full set (it is monotone);
  - min-cardinality gives 9 tied singletons;
  - GS1 selection power gives full access only.

**A_partition / A_seed REMAINS SUPPLIED.**

### B3-5 EA-0 category firewall

| statement | classification |
|---|---|
| Earned classical Level-0: the EA-0 quantum access object (P_ρ, Γ(ρ)) is undefined without a lift | **BLOCKED / CATEGORY MISMATCH** (carried exactly; not CLASS-SPLIT) |
| Auxiliary class: faithful states give P_ρ = I | TRIVIAL/IDENTITY (canonical; not re-run) |
| Auxiliary class: generic Γ(ρ) is non-selective | canonical L3′, at d = 4 and 8 only |
| Auxiliary class: non-trivial selective access only where D has exact invariant / block structure (P-6 "selected exactly up to its own block decomposition"; B3-1 parity block, dim 8, center 2) | **RELOCATION INTO D**: the access is recovered because D already contains the sectors. Theorem S is sufficient only; the retired "generator-carried iff" is **not** revived |

### B3-6 State-dependent values vs access structure — **H_corr is not counted as A compression**

The MI between site 1 and the rest at t = 1 is 0.071 for an uncorrelated preparation and 0.150 for a correlated one. The
closure stays at dimension 24, and the available effects, fragments and seed are unchanged. That is EA-0 case (a), values
on a fixed structure, not case (c), a change of structure. The same holds for SCOUT's MI / redundancy / distinguishability.

### B3-7 Resolution / hierarchy order — **A_resolution REMAINS SUPPLIED**

- **Order.** With the end-site seed fixed and one hidden edge changed at depth r = 3 / 6 / 10, the Markov moments agree
  through order **31 / 34 / 39** and first differ at **32 / 35 / 40**. The transfer differences are 3e-3 / 1e-4 / 2e-6.
- **Threshold.** The number of detectable edges (out of 23) at δ = 1e-2 / 1e-4 / 1e-6 / 1e-8 is **2 / 7 / 11 / 15**.
- The required order and δ are set by the unknown alternative, not by D. "Higher order eventually distinguishes" is not
  "the required order is derived". This is consistent with P-6 L-D1 and the GS1 horizon.

### B3-8 Path-dependent minimality — **MINIMALITY DOES NOT SELECT A_seed**

- On the uniform 5×5 grid there are 11 degenerate pairs, so no single readout is complete.
- In the pool {0, …, 4}, backward elimination gives **{1, 2, 3, 4}** forward and **{0, 1, 2, 3}** in reverse. Both are
  sufficient **and** minimal.
- All five 4-subsets are minimal.
- This classical re-implementation reproduces the P-6 L-D2 logic. The pool is the first lexicographic pool with differing
  outcomes (selection disclosed).

### B3-9 / B3-10 Grouping and temporal access — **A_partition, A_time REMAIN SUPPLIED**

**Record criterion:** on the C1-a net, a fragment "records" x₁₂(0) if the posterior variance from its noisy samples
(σ = 1e-3, prior variance 1) falls below the threshold. Counts are given as singles (of 24) / adjacent pairs (of 12) /
blocks of 6 (of 4).

| horizon T | cadence Δt | threshold < 0.6 | threshold < 0.3 |
|---|---|---|---|
| 1 | 0.5 | 1 / 2 / 2 | 1 / 1 / 2 |
| 5 | 0.5 | **3 / 5 / 2** | 1 / 3 / 2 |
| 20 | 2.0 | 1 / 2 / 2 | 0 / 0 / 0 |

- The "record redundancy" changes with grouping, horizon, cadence and threshold.
- D supplies **units** (relaxation rates 0.30 … 4.28), not a schedule.
- No canonical rule individuates fragments from D + Σ. The existence of the local net is not counted as a fragment
  selector.
- P-5's boundary-change event is where access becomes physical, but *when* it happens (t*) is declared.

### B3-11 / B3-12 Absorption tests

| test | evidence | verdict |
|---|---|---|
| **B3-11** (into Σ) | same Σ and D, different access, different interface data (B3-2, B3-3, B3-4) | **A_res ≠ Σ** |
| **B3-12** (into D) | same D, seeds that closure and center cannot distinguish, different physics (B3-2 classical and quantum) | **A_res ≠ D** |

## B3-13 Compression outcome (verdicts by component)

| component | verdict |
|---|---|
| A_closure / effect structure | **B: PARTIAL / CONDITIONAL ACCESS COMPRESSION**: f(D, A_seed; closure rule) |
| state reconstruction given the readout | **B** (earned, identity-grade; canonical) |
| block-level selection of access in the auxiliary class | **C: RELOCATION INTO D** |
| A_seed, A_readout, A_partition, A_resolution, A_time | **D: NONUNIQUE / SUPPLIED** |
| EA-0 quantum access objects on the earned substrate | **E: BLOCKED AT LEVEL-0** |
| — | **A (TRUE ACCESS COMPRESSION): not obtained** |

**Sharpened A_res (tested classes):**

    A_res  →  A_seed(+readout identification) ⊕ A_partition ⊕ A_resolution ⊕ A_time
    with  A_closure / effect structure = f(D, A_seed; closure rule)

- Seven listed SCOUT items reduce to **four supplied primitives plus one convention**:
  1. coupled/accessible + effect/readout → **A_seed**, with the effects downstream;
  2. grouping + coarse-graining → **A_partition**, but only where the coarse map is a fragment restriction. General
     aggregating coarse maps were **not tested**, so A_coarse ⊆ A_partition is **CONDITIONAL**;
  3. resolution / threshold → **A_resolution**;
  4. time + horizon → **A_time**.
- **A_partition is not reducible to A_seed.** The same total readout (all sites) grouped differently gives different
  record counts (B3-9).

## B3-15 Prediction firewall

No empirical payoff is claimed. The baseline stays at **zero confirmed distinctive GRUT quantitative predictions**. Only
B5 can test whether the removed freedom (the closure as a separate item) forces any observable relation. A priori it
should not, since the closure is a bookkeeping consequence.

## Scope and limits

- The numerics are a single agent's code path: an **independent code path, not an independent reviewer**.
- Small systems only: N = 24 chains, a 3×3 and a 5×5 grid, 2–3 qubits.
- The quantum results are auxiliary-class only.
- Interior-site and multi-site readouts are admissible declarations but **not earned** (L0 bridge §5). They are used here
  only as witnesses that alternatives exist.
- Linear classical class only. The nonlinear on-site case is not re-tested; L0 bridge Theorem 1 covers end-site
  observability there.
