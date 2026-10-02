# BRIDGE ZOOM-OUT 03 (after BRIDGE REPAIR 02 + B2) — HARD STOP for owner review

**Branch:** `grut-bridge-1`. **Canonical GRUT:** read-only at `master-w25bu9 @ b935099`, untouched.
**Evidence:** `B2_HCORR_RESULT.md`, `B2_BOUNDARY_COMPONENT_LEDGER.md`, `b2/b2_hcorr.log`. The B3 repairs are in
`BRIDGE_CORRECTION_LEDGER.md` (BR2-01, BR2-02).

## 1. Is GRUT's S6 preparation merely one explicit H_corr|Σ boundary?

**Yes.** The S6 product state (H_sys, H_bath, H_cross = 0 at t = 0) is one admissible member of a large class. Under the
same D, Σ, A and bath law, other members include:
- correlated preparations;
- the marginal-matched product;
- a grafted correlation block;
- the time-reversed state;
- Gibbs itself.

They all behave differently (X_J(∞) ranges from −0.73 to +1.17; some relax, Gibbs does not; the reversed state
anti-relaxes). The S6 preparation plays exactly the role of SCOUT's H_corr|Σ. That is **RENAMING / IDENTICAL in role**,
reached independently.

## 2. Is temperature difference the real price, or correlation structure?

**Correlation structure.** The key evidence:

| state | temperatures | result |
|---|---|---|
| S6 factorized product | T_s = T_b | relaxes, with X_J(∞) = ½Tr² = 0.169426 (finite N 0.169425) |
| global Gibbs | same T | stationary (drift 9e-16) |
| marginal-matched product (Gibbs marginals, no cross-correlation) | every marginal at its equilibrium value | relaxes, with X_J(∞) = Tr² = 0.33885 |

- A and Gibbs differ only by the missing S–B correlation and the interaction-dressing of the marginals.
- With the S6 marginals fixed, changing only the cross block **flips the sign** of X_J(∞) at equal T (α > ½). The exact
  relation is X_J(∞) = (T_s − T_b) + T_b r²(½ − α).

> **TEMPERATURE DIFFERENCE NOT NECESSARY. CORRELATION / INTERACTION-ENERGY MISMATCH IS PHYSICALLY LOAD-BEARING.
> H_cross IS AN INDEPENDENT PREPARATION DATUM.**

## 3. Is productness necessary, sufficient, or one convenient preparation?

**Sufficient, not necessary, for local relaxation. It is one convenient member.**
- **Proposition B2-P1** (bridge theorem; reviewed at stated scope [BR3 review note]): the LS-1 mechanism (Gibbs invariance + a.c.
  spectrum + Riemann–Lebesgue) extends verbatim to **any trace-class** covariance perturbation, cross blocks included.
- Finite-N numerics confirm relaxation for rank-2 to rank-5 correlated deviations.
- The residual is therefore not "the preparation must be independent". It is **which correlation-boundary member was
  realized.**

## 4. Does the environment select any part of H_corr?

**Only partly, and conditionally.**
- Given a supplied Gibbs postulate and T_b, the bath state fixes the local asymptotic state (**PARTIAL RELOCATION INTO
  ENVIRONMENT**). Even this is conditional: the integrable chain also has stationary non-thermal GGEs, which give a
  different local limit (0.5622 vs 0.5821).
- The environment does **not** determine T_b, T_s, H_cross or the epoch.
- The infinite bath makes every trace-class H_corr locally forgettable, but **the integrated record keeps it**:
  X_J(∞) = E₁(0) − E₁_G + Q12_G − Q12(0). That is **D COMPRESSES MEMORY OF H_corr DOWNSTREAM**, not derivation.
- The stochastic noise layer fixes a correlated stationary state given T_i (RELOCATION) and erases the initial H_cross. It
  selects nothing. It was kept separate from the S6 bath.

## 5. Is the special epoch selected?

**By the preparation, not by D.**
- t = 0 is the zero-correlation point (C_SB = 7e-15).
- Conservative evolution generates S–B correlation in **both** time directions (C_SB(±t) → the Gibbs value 0.417).
- **SPECIAL MOMENT SELECTED BY PREPARATION.** Up to time translation of a history, the epoch is a gauge. It means only
  "where the declared form (productness) holds".

## 6. Is temporal orientation selected?

**No.**
- Σ₀ = RΣ₀R implies Σ(−t) = RΣ(t)R, so C_SB and D are even in t and J is odd: Janus-symmetric relaxation.
- A reversed preparation, with identical global and marginal entropies, runs back to the less-equilibrated state exactly
  (D_rev(t) = D_fwd(τ − t)).
- S6's NET-ARROW-CONFIRMED stands as stated. S6 contains two declared notions of forward [BR3-02]: a temperature-ordering
  sign for J (f_J = ±J on L1 / L2) and descent toward the supplied reference S_ref for D (σ = −Ḋ). Neither is selected
  by the time-symmetric dynamics. At T_s = T_b the J ordering disappears, yet the correlations set the sign of J.
- Proposed **B2-CUC-1** (bookkeeping clarification; not applied).
- **ORIENTATION NOT SELECTED; D + Σ + A DO NOT SELECT THE ARROW BOUNDARY.**

## 7. Which compression classes occurred?

| class | B2 instance |
|---|---|
| **TRUE COMPRESSION** | **none** |
| **CONDITIONAL COMPRESSION** | the local asymptotic state = f(bath state) via the a.c. relaxation mechanism, i.e. the memory of H_corr compressed downstream (B2-P1, bridge theorem reviewed at stated scope) |
| **RELOCATION** | reference / local-limit state into the bath state (partial, given Gibbs + T_b); stationary correlation structure into the noise parameters T_i |
| **RENAMING** | S6 declared product preparation ↔ SCOUT H_corr\|Σ |
| **SUPPLIED** | H_sys, H_cross, H_epoch, T_s, T_b, the Gibbs-vs-GGE choice, the orientation convention |
| **not the price** | temperature difference (equal-T transient); low entropy (the higher-entropy marginal-matched state relaxes; equal-entropy reversed states anti-relax) |

## 8. Did the reviewed residual shrink?

**No.** It is made explicit on the GRUT side:

    C5 → D_dyn ⊕ [Σ ⊗ H_corr]_coupled ⊕ A_res
       = D_dyn ⊕ [Σ ⊗ (H_marginals, H_cross) @ declared-form epoch] ⊕ (A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time)
    with A_closure = f(D, A_seed; R_closure),  A_interface := (A_seed, A_readout)

- Time orientation is not selected, and Lorentz structure is not derived.
- **TRUE COMPRESSION remains 0** across B1, B3 and B2.
- What Bridge-1 has gained is **convergence**. GRUT's canonical S6 algebra independently gives SCOUT's central sharpening:

> **The boundary price is correlation structure relative to the split, not low entropy or a temperature gradient.**

Empirical status: **zero confirmed distinctive GRUT quantitative predictions.** The generalized X_J relation is a model
consequence given a declared preparation, not a parameter-free observable.

## 9. Is Bridge-1 ready for B4?

**Yes.** Every reviewed residual component has now been audited against its GRUT counterpart:
- Σ (B1);
- A_res (B3);
- H_corr (B2);
- D_dyn (B0: IDENTICAL on "not derived").

The remaining layers (lift, ħ, Born rule, gravity, cosmology) are NO MAPPING / BLOCKED by rule. B4 can fill the nine-layer
× residual compression matrix from the existing ledgers without new physics. The open sub-items (B1-o1/o2, B3-o1/o2/o3,
B2-o1/o2/o3) refine scope but do not block it. **B2-o1, independent verification of Proposition B2-P1, is the most
valuable**, because it upgrades "productness is not necessary" from a bridge sketch.

## Status

- **BRIDGE REPAIR 02:** applied. **B2:** DONE.
- **TRUE COMPRESSION:** 0.
- **CANONICAL UPDATE CANDIDATES:** B1-CUC-1, B3-CUC-1, B2-CUC-1 (none applied).
- **HARD STOP. Awaiting owner review before B4.**
