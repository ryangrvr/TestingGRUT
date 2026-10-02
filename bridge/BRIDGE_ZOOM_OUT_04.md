# BRIDGE ZOOM-OUT 04 (after BRIDGE REPAIR 03 + B4) — HARD STOP for owner review

**Branch:** `grut-bridge-1`. **Canonical GRUT:** read-only at `master-w25bu9 @ b935099`, untouched.
**B4** is pure synthesis. No new numerical physics was run.

**Evidence:**
- `B4_NINE_LAYER_COMPRESSION_MATRIX.md`
- `B4_CANONICAL_UPDATE_CANDIDATES.md`
- `BRIDGE_CORRECTION_LEDGER.md` (BR3-01, BR3-02, B2-P1 review note)

## 1. Did any canonical supplied item become unnecessary?

**No.** The matrix has **0 TRUE COMPRESSION cells** (counted, not assumed). Every Ledger A entry that the bridge touched
(A-1 … A-4, A-9, A-11 … A-15) is still needed.

## 2. Did any layer become internally decomposable?

**Yes: four layers.** This is regrouping, not physical compression.

| layer | decomposes into |
|---|---|
| **2** | static-coupling-compatible generator class · optional drift · a **class-relative orientation carrier** (generator sign in G-D; preparation + two declared forwards in S6) |
| **3** | A_interface (= seed, readout) ⊕ A_partition ⊕ A_resolution ⊕ A_time, with A_closure = f(D, A_seed; R_closure) downstream |
| **4** | the sector label and statistics (**not** H_corr) vs H_marginals, H_cross, H_epoch and the S–B split |
| **5** | environment origin · bath state · noise law · partition listing · coarse map |

## 3. Which redundancies are genuine?

| ID | redundancy | grade |
|---|---|---|
| **R1** | Σ is re-encoded in the on-site drift and in non-uniform site noise | **genuine, class-scoped** (REDUNDANT SUPPLY / CONSISTENCY; recovery is exact, and the controls select whatever frame each was written in) |
| **R2** | the system–bath split, listed in layer 4 (A-14) and in the layer-5 name | **genuine duplicate listing** (bookkeeping only) |
| **R3** | the S6 bath state, listed under the preparation and under the environment | **genuine duplicate listing** |

**Not redundancies:**
- **R2′:** the split vs Σ vs A_partition only overlaps.
- **R3′:** S_ref is downstream, not duplicate.
- **R3″:** the noise stationary state vs the S6 bath is a different class.
- **B4-o1:** the A-5 price labels vs A-13 / A-14 are untested.

## 4. What is the reviewed dependency graph?

**Supplied roots:**
- K;
- the generator class with its class-relative orientation carrier;
- Σ;
- optionally, the drift and T_i, which re-encode Σ;
- A_seed, A_readout, A_partition, A_resolution, A_time and R_closure;
- the sector label;
- H_marginals (H_sys, plus the bath's T_b and Gibbs postulate), H_cross and H_epoch;
- the S–B split;
- the coarse map.

**Derived edges (all conditional):**

| inputs | → | output |
|---|---|---|
| (K, generator, A_seed, R_closure) | → | A_closure |
| (K, A_readout) | → | state reconstruction |
| (Σ, A, K) | → | geometry (one way only) |
| (K a.c., bath state) | → | S_ref, the local asymptotic state |
| (H_marginals, H_cross, K) | → | X_J(∞) |
| T_i | → | the noise stationary correlations |

**Fenced:** lift, ħ, outcome, gravity and cosmology.

**Non-edges:**
- K ↛ Σ;
- D + Σ ↛ any A block;
- D + Σ + A + environment law ↛ any H component;
- D ↛ orientation;
- geometry ↛ access.

(Diagrams are in the matrix file, §§B4-9 and B4-10.)

## 5. Is the nine-layer architecture physically smaller, bookkeeping-smaller, both or neither?

| question | answer |
|---|---|
| **Physically smaller?** | **No**: physical compression = 0 |
| **Smaller as a dependency structure?** | **Yes**: 2 items downstream (A_closure, S_ref), 2 duplicate listings mergeable, 1 redundancy edge, 4 layers decomposed |
| **Smaller as a count of explicit accounting items?** | **No: the explicit supplied / declarative accounting-item count grows** [BR4-01]. 6 items supplied in practice but not booked separately are made explicit: A_resolution, A_time (protocol declarations), R_closure (convention), H_epoch (coordinate value gauge; the special-form event is physical), the Gibbs-vs-GGE postulate (state-class choice) and H_cross (genuine physical boundary data). **Physical primitive count change: NOT ESTABLISHED** |

So the architecture is **bookkeeping-cleaner, not smaller in explicit accounting items**; no new count of physical primitives is claimed. Categories: **B (bookkeeping
compression only) + D (new redundancy / overdetermination).** Not A, not C.

## 6. Which CUCs survive B4?

| CUC | grade | status |
|---|---|---|
| **B1-CUC-1** | POTENTIAL CANONICAL DEPENDENCY UPDATE | survives |
| **B3-CUC-1** | BOOKKEEPING ONLY (plus the closure-rule disagreement for the owner) | survives |
| **B2-CUC-1** | SCIENTIFIC SCOPE CLARIFICATION, carrying the BR3-02 split of the two S6 forward notions | survives |
| **B4-CUC-1** (new) | BOOKKEEPING ONLY: merge duplicate listings; book the hidden items; mark S_ref downstream | new |
| **B4-CUC-2** (new) | SCIENTIFIC SCOPE CLARIFICATION: the orientation carrier is class-relative | new |

**None is applied.** No publication rewrite is recommended.

## 7. Does the SCOUT residual change after the full GRUT comparison?

**Not in size. Yes in resolution.**

    C5 → D_dyn ⊕ [Σ ⊗ (H_marginals, H_cross, H_epoch)] ⊕ [A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time]

where:
- D_dyn = static coupling K ⊕ generator class ⊕ a class-relative orientation carrier. The K vs generator split is a
  canonical GRUT contribution (S5).
- A_closure = f(D, A_seed; R_closure) is downstream, and R_closure is a supplied convention.
- The local asymptotic state is downstream of the bath part of H_marginals.
- Σ is redundantly re-encoded in site-local drift / noise when those are supplied.
- Orientation is **not selected**.
- Lorentz / causal cone is **supplied** in both programs (convergence, not compression).
- Lift, ħ, outcome and gravity / cosmology are outside the C5 bridge.

**GRUT's distinctive contributions to the residual:**
- the D_dyn split;
- an exact closed-form witness that the H price is correlation structure (B2);
- the R1 redundancy;
- the closure-rule finding.

## 8. Is there any plausible route for B5 to produce a distinctive observable?

**Plausible routes exist only as tests. Each looks non-distinctive on its face.**

| route | why it looks non-distinctive |
|---|---|
| **R1 overdetermination** (drift axes = noise eigenframe = net frame) | a checkable consistency relation, but generic to any site-local model |
| **The generalized relation** X_J(∞) = E₁(0) − E₁_G + Q12_G − Q12(0) | requires the supplied Q12(0), so it is not parameter-free |
| **The equal-T offset** ½T_b r² | a standard harmonic-chain number at the declared K₁₁ = 2.3 |

The expected B5 terminal is **NO DISTINCTIVE OBSERVABLE**, but B5 must test this rather than assume it. The baseline
stays at zero confirmed distinctive GRUT quantitative predictions.

## 9. Is Bridge-1 approaching saturation?

**Yes.** Every reviewed residual component has been audited (Σ: B1; A: B3; H: B2; D: B0 + B4), and the matrix closes
with 0 TC cells. The remaining open sub-items refine scope but cannot create TRUE COMPRESSION within the current
envelope:
- B1-o1, B1-o2;
- B3-o1 … o3;
- B2-o2, B2-o3;
- B4-o1.

After B5, Bridge-1's final product will most likely be **a more precise dependency theory**, with a sharper residual and
a short list of bookkeeping / scope CUCs, rather than new physics.

## Status

- **BRIDGE REPAIR 03:** applied. **B2:** accepted provisionally. **B4:** DONE.
- **Physical compression:** 0. **Dependency compression:** > 0. **Explicit accounting-item count:** ↑ (hidden items explicit). **Physical primitive count change:** NOT ESTABLISHED [BR4-01].
- **HARD STOP. Awaiting owner review before B5.**
