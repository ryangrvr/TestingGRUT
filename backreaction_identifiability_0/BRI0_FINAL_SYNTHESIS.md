# BRI0 — FINAL SYNTHESIS

**Terminal:** **BRI0 COMPLETE — BACK-REACTION IDENTIFIABILITY EXISTS RELATIVE TO THE SHARED CAUSAL AFFINE EXOGENOUS CLASS;
UNIVERSAL EXOGENOUS NON-IDENTIFIABILITY REMAINS.**

**Branch:** `grut-backreaction-identifiability-0`. Frozen parent `scout-0 @ ab2da47` (unmodified). BRI0 is a
**new-premise campaign**: it does not reopen the SCOUT-0 saturation verdict.

**Reviewed boundaries:**

| boundary | content |
|---|---|
| `89b8236` | Charter 0 |
| `8df4b88` | Scope Repair 01 |
| `6220a9e` | Scope Repair 02, Candidate 1, preflight |
| `d5a0bdb` | PF4Q |
| `845b513` | BRI1-X1-THEOREM (accepted) |

**The ladder:** 𝓗 ⊂ E₁ ⊊ E₂± ⊊ E_univ. Every class is a **shared** model across the interventional clamp family 𝒳
(BRI0 §3: a single protocol never discriminates).

## 1. E₁ — additive exogenous forcing plus deterministic causal memory

- **P-17's class 𝓗** (harmonic bath, coupling linear in bath and system coordinates, system-independent free-force law)
  **lies in E₁** under clamping (PROP BRI-H).
- **Additive forcing plus arbitrary deterministic causal memory cannot identify bath ontology.**
- **Exact characterisation (BRI-C1):** E₁ membership holds exactly when the centred clamped force law is
  protocol-invariant.

## 2. E₂± — shared causal signed affine (location-scale) modulation of one exogenous process

- **C1** (harmonic bath, nonlinear system coupling A(q) = q + q³/3) **lies in E₂± exactly** (PROP BRI1-C1). Nonlinear
  system coupling and multiplicative noise **alone do not escape E₂±**; C1 is structurally BRI-E1.
- **X1** (a finite reciprocal Duffing bath with coupling 1/√N_B) **escapes E₂±** for every sufficiently large finite N_B
  (**THEOREM BRI1-X1; BRI-E2+O at class-theorem level**). The mechanism:
  - P0 has third cumulant exactly 0 (parity);
  - P1 has a rigorously negative leading coefficient on a small-time interval;
  - BRI1-R1 gives κ₃ = K/N_B + O(N_B⁻²);
  - the variances are positive;
  - |standardised skewness| is reflection-safe;
  - so there is an orbit violation.
- **The earned observable information** is a change in the **standardised response-law shape** that **cannot be removed
  by one shared causal affine modulation** (location, magnitude and sign).
- **Exact characterisation (BRI-C2±):** E₂± membership requires a common degeneracy set and one shared causal sign
  functional S_t[q_[0,t]] that aligns all standardised laws. Positive subtypes are BRI-E2+O (orbit / shape) and BRI-E2+C
  (causal coherence).

## 3. E_univ — the ceiling (BRI-UPPER)

- Every environment in the declared **deterministic classical causal** parent class is representable as
  **F_q = 𝔉[q, U]**: one exogenous random object U passed through a causal functional. Reciprocal and energy-absorbing
  back-reaction is included, **X1 included**.
- **Therefore BRI1 does not identify primitive randomness or a unique microscopic ontology.**
- The discrete-time stochastic extension is KNOWN (randomisation). No continuous-time general-stochastic claim is made,
  and quantum environments are out of scope.

## 4. Reservoir limit — the X1 affine escape is mesoscopic

- The non-affine skewness witness is **O(1/N_B)** and **vanishes as N_B grows**.
- **E₁-type reservoir limit for the frozen / pointwise-fixed protocol family:**
  - for P0 / P1 / P2, the centred force laws converge (fdd, via Lindeberg–Feller with the BRI1-R1 bounds) to the **same**
    Gaussian linear-response law;
  - the same pointwise argument applies to any separately fixed admissible bounded clamp satisfying the moment estimates.
- **Scope kept explicit:** no theorem gives one shared E₁ representation uniformly over the whole infinite clamp class 𝒳,
  and none is claimed.

## 5. Scientific result

> **Reciprocal anharmonic back-reaction can make information about environmental response visible in interventional
> reduced-force laws beyond deterministic location, scale and sign modulation, while still remaining universally
> representable by a causal exogenous random object.**

**Relation to SCOUT-0.** P-17 (class 𝓗) and P-18 (one-way drivers) stand exactly. Their non-identifiability mechanism
does **not** extend universally to finite reciprocal anharmonic baths under the affine competitor E₂±: X1 is an explicit
counterexample. The result locates the **boundary of the affine quotient**, not of universal exogenous representation.

## 6. What is NOT earned

- a GRUT-specific physical law;
- a GRUT empirical prediction;
- primitive randomness;
- a unique ontology;
- consciousness;
- TRUE COMPRESSION;
- an escape from E_univ.

## 7. Verification status

| item | status |
|---|---|
| Charter / class theorems (vacuity, the ladder and its strictness, BRI-C1, BRI-C2±, BRI-H, BRI-UPPER) | INTERNALLY PROVED / NOT EXTERNALLY REVIEWED |
| BRI1-C1, BRI1-R1, BRI1-T1, THEOREM BRI1-X1 | INTERNALLY PROVED / NOT EXTERNALLY REVIEWED |
| Exact symbolic series (C1 identity; small-t coefficients) | exact computer algebra, one code path |
| Frozen-τ K-coefficients (20 values) | **strong numerical evidence only — NOT CERTIFIED**. PF4Q-I / X1-PF-INDETERMINATE at the frozen τ, permanently |
| Certified-τ pipeline | documented possible future verification. Not built |
| Independent reproduction of the theorem | **future verification work**, if publication is contemplated |

## Stop rule

**BRI2 is not opened.** No second anharmonic candidate, no search for a larger effect, no parameter tuning, and no
attempt to make the result GRUT-specific. Lane 3 has answered its chartered question.
