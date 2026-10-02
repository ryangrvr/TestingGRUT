# B2 — H_corr VERSUS GRUT PREPARATION / ENVIRONMENT (result)

**Question (owner).** With Σ, A_interface, A_partition, A_resolution and A_time fixed, does any GRUT structure select
H_corr? Or does GRUT supply one member of the correlation-boundary class that SCOUT isolated?

**Script:** `bridge/b2/b2_hcorr.py` (+ `b2_hcorr.log`). It is an independent finite-N re-implementation of the canonical
S6 parent, read only from `S6_1_THEOREM_01.md` and `S6_1_COARSEGRAINED_ARROW_CHARTER_01.md`; `calc/s6_1_certify.py` was
not used.
- **Model:** N = 600, K = 2.3·I − T, q̈ = −Kq, classical Gaussian ensembles.
- **Finite-N role:** all times are far below the recurrence time. The finite-N values reproduce every S6 closed form to
  about 1e-5 (B2-0). They are **cross-checks, never adjudicating**. Theorem-grade statements are only the canonical LS-1 /
  LS-2 and the bridge algebra stated below.
- **Admissibility of a classical Gaussian state:** a positive-definite covariance, certified by its minimum eigenvalue
  for every state used. The quantum uncertainty condition does not apply (ħ is not reopened).

**Access held fixed throughout:**
- the retained 2×2 block of site 1;
- the coupling observable J = ⟨p₂q₁⟩;
- the same times and the same split (site 1 | bath {2, 3, …}).

**Label:** independent code path, not independent reviewer.

## B2-0 The canonical preparation, decomposed

S6 declares Σ₀ = Σ_G + Δ₀, with Σ_G = T_b·diag(K⁻¹, I), Δ_Q = (T_s/2.3)e₁e₁ᵀ − (T_b/r)uuᵀ and Δ_P = (T_s − T_b)e₁e₁ᵀ,
so rank Δ₀ ≤ 3. The re-implementation finds rank 3, and X_J(∞) matches the closed form at all four members.

| datum | canonical value | classification |
|---|---|---|
| **H_sys** (system marginal) | ⟨q₁²⟩ = T_s/2.3 (the *bare* site-1 Gibbs variance, not the dressed T_b·r), ⟨p₁²⟩ = T_s, ⟨q₁p₁⟩ = 0 | supplied (T_s and the bare form) |
| **H_bath** (bath marginal) | T_b·K_BB⁻¹ and T_b·I: the Gibbs state of the *uncoupled* bath | supplied (T_b and the Gibbs postulate); see B2-7 / 8 |
| **H_cross** (system–bath correlations) | Q_SB = 0, P_SB = 0, all mixed q–p blocks = 0 | supplied |
| **H_epoch** | the declared t = 0 at which the cross blocks vanish | supplied |
| **Σ_split** | site 1 vs bath {2, 3, …} | Σ (fixed) |
| reference | S_ref = T_b·diag(r, 1) and the forward-sign convention (B2-6) | convention |

None of these is "low entropy". B2-11 shows entropy is not the operative label.

## B2-1 Equal-temperature control (mandatory) — **CONFIRMED**

Same D, Σ and A, and T_s = T_b = T = 1.

| | State A (S6 factorized, equal T) | State B (global Gibbs, T·diag(K⁻¹, I)) |
|---|---|---|
| min eigenvalue | 0.233 | 0.233 |
| max local drift over t ≤ 40 | **0.339** (relaxes) | **9e-16** (stationary) |
| X_J(∞) | **+0.169425** (closed form ½Tr² = 0.169426) | 0 |
| D(0) | 0.01936 (> 0 because T/2.3 ≠ T·r) | 0 |

**What distinguishes A from B** (every temperature equal). rank(A − B) = 2, entirely in the q-sector; the p-sector is
identical.
- **S–B cross correlation:** ⟨q₁q₂⟩ is 0 in A and **T·r² = 0.33885** in B. The cross-block norm is 0 vs 0.417.
- **Interaction energy:** ⟨E_int⟩ is 0 in A and −T·r² in B.
- **Dressed marginals:** the site-1 q-variance is the bare T/2.3 = 0.435 in A and the dressed T·r = 0.582 in B. The bath
  q-block differs by a rank-1 dressing.

All of the difference is interaction-induced.

**Isolating H_cross.** Take the product of the **Gibbs marginals** with zero cross block ("marginal-matched", state M).
- S₁(0) equals S_ref exactly, so D(0) = 1e-16, and the bath marginal is exactly Gibbs.
- **The state is still not stationary:** D(t) rises to 4.95e-2 at t = 1 and decays back.
- X_J(∞) = **T·r² = 0.33885**, matching the closed form Q12_G − Q12(0). The marginals contribute nothing.

> **TEMPERATURE DIFFERENCE NOT NECESSARY FOR THE S6 TRANSIENT.**
> **CORRELATION / INTERACTION-ENERGY MISMATCH IS PHYSICALLY LOAD-BEARING.** With every marginal at its equilibrium value,
> the missing system–bath correlation alone drives the transient and a net J transfer.

## B2-2 Same marginals, different correlations — **H_cross IS AN INDEPENDENT PREPARATION DATUM**

The S6 marginals are held fixed and the S–B cross covariance is varied:
- a q-cross along the Gibbs direction, Q_{1B} = α·T_b·(K⁻¹)_{1B};
- a p-cross β;
- a mixed q₁–p₂ cross γ;
- a cross block grafted from the S6 state evolved to t = 3.

Admissibility at equal T is certified by the Schur complement T(1/2.3 − α²r³) > 0, i.e. |α| < 1.4847, and by the minimum
eigenvalue for every case. Every state has S₁(0) equal to the S6 value.

**Exact relation (bridge algebra).** This follows from the LS-2.1 identity, which holds for any state, plus return to
equilibrium (B2-3):

    X_J(∞) = E₁(0) − E₁_G + Q12_G − Q12(0),
    at the S6 marginals:  X_J(∞) = (T_s − T_b) + ½T_b r² − Q12(0) = (T_s − T_b) + T_b r²(½ − α)

| marginals | cross | X_J(∞), finite N | max \|S₁ − S₁^{S6}\| at t = 2 | X_σ(∞) |
|---|---|---|---|---|
| (1, 1) | none (S6) | +0.16942 | 0 | 0.01936 |
| (1, 1) | α = 0.3 | +0.06777 | 0.048 | 0.01936 |
| (1, 1) | α = 1.0 | **−0.16942** | 0.161 | 0.01936 |
| (1, 1) | α = 1.4 | **−0.30496** | 0.226 | 0.01936 |
| (1, 1) | β = 0.5 (p-cross) | +0.16943 | 0.166 | 0.01936 |
| (1, 1) | γ = 0.3 (q₁–p₂) | +0.16942 | 0.219 | 0.01936 |
| (2, 1) | α = 0.3 / 1.0 / 1.4 | +1.068 / +0.831 / +0.695 | 0.05 – 0.23 | 0.19967 |
| (2, 1) | graft from t = 3 | +1.007 (S6: +1.169) | 0.291 | 0.19967 |

**With every marginal fixed:**
- correlation changes alter the retained transient, J(t), and the **sign** of the integrated transfer (equal T, α > ½);
- the p-cross and mixed crosses change the transient and J(t) but not X_J(∞), since only Q12(0) enters the closed form;
- X_σ(∞) = D(0) is blind to all of them, because it sees only S₁(0) and S₁(∞).

This is the GRUT-side counterpart of SCOUT D3 / D6 (correlations, not marginals, carry the arrow-relevant boundary datum),
in a classical Gaussian infinite-bath model.

## B2-3 Productness is sufficient, not necessary — **bridge theorem (reviewed at stated scope) + numerics**

**Proposition B2-P1** (bridge; **BRIDGE-THEOREM REVIEWED AT STATED SCOPE** [BR3 review note]; not a canonical GRUT theorem).
- **Setup:**
  - K = K_∞ = 2.3·I − T on ℓ²(ℕ);
  - Σ_G = T_b·diag(K⁻¹, I);
  - Δ is any symmetric **trace-class** perturbation on ℓ² ⊕ ℓ² (in particular any finite rank) with Σ₀ = Σ_G + Δ ⪰ 0;
  - S–B cross blocks are allowed.
- **Claim:** for every fixed pair of local coordinates i, j, (Σ(t) − Σ_G)_{ij} → 0 as t → ∞.

*Sketch.*
1. Write Δ = Σ_k σ_k a_k a_kᵀ with Σ|σ_k| < ∞ and ‖a_k‖ = 1. Then (Σ(t) − Σ_G)_{ij} = Σ_k σ_k (M(t)a_k)_i (M(t)a_k)_j,
   using Gibbs invariance (LS-1.1).
2. Each (M(t)a)_i is a finite sum of ⟨e_i, h_t(K) b⟩, with h_t ∈ {cos √λt, λ^{−1/2} sin √λt, λ^{1/2} sin √λt} and b ∈ ℓ².
3. Under LS-0's unitary map ℓ² → L²(ρ) (e₁ is cyclic), this equals ∫ U_{i−1}(x/2) h_t(2.3 − x) b̂(x) ρ(x) dx.
   - The weight U_{i−1} b̂ ρ is in L¹ (Cauchy–Schwarz, finite measure).
   - After the smooth change of variables ω = √(2.3 − x) (with 2.3 − x ≥ 0.3), Riemann–Lebesgue gives → 0.
4. |(M(t)a)_i| ≤ c_i‖a‖ uniformly in t, because the h_t are bounded on [0.3, 4.3]. Dominated convergence over k then
   gives the claim. ∎

**Only three ingredients are used:** Gibbs invariance, trace-class Δ, and absolute continuity. **Productness is not
used.** The proof is the LS-1 mechanism, unchanged.

**Numerics (finite N).** |S₁(t) − S_ref| at t = 200 is below 1e-5 for every correlated case: rank-3, -4 and -5 deviations,
the marginal-matched product and the graft.

> **PRODUCTNESS IS SUFFICIENT, NOT NECESSARY, FOR LOCAL RELAXATION** (bridge theorem, reviewed at stated scope). The residual is not "the
> preparation must be independent". It is **which correlation-boundary member was realized.**

## B2-4 Correlation reversal / anti-arrow — **D + Σ + A DO NOT SELECT THE ARROW BOUNDARY**

**Construction.**
1. Take S6 (2, 1) and evolve it globally to τ = 8.
2. Form Σ_R = R·Σ(τ)·R, where R flips p → −p, and declare it the new t = 0. D, Σ and A are unchanged.
3. Σ_R is admissible (min eigenvalue 0.233).

**Result.**
- **D_rev(t) = D_fwd(τ − t) exactly** at t = 0, 2, 4, 6, 8, in the sequence 0.00028 → 0.00122 → 0.01409 → 0.14313 →
  **0.19967 = D(0)**. The retained state runs **back** to the less-equilibrated S6 preparation.
- X_J over [0, τ] is +1.137 forward and −1.137 reversed.
- **The global entropy and the retained-marginal entropies of Σ(τ) and Σ_R are identical** (difference 0.0). The two
  states differ only in the sign of their q–p correlations.

This is the direct GRUT analogue of SCOUT's D-arrow mirror (D1 / D7). The infinite-chain model is outside D0's
finite-dimensional scope, so **D0 is used only as a conceptual comparison.**

## B2-5 Special epoch — **SPECIAL MOMENT SELECTED BY PREPARATION; ORIENTATION NOT SELECTED**

Start from the S6 product state at t = 0.
- The cross-block norm is **C_SB(0) = 7e-15**. It is C_SB(±t) = 0.535, 0.810, 0.735, 0.510, 0.421 and 0.419 at
  t = 0.5, 1, 2, 5, 10, 20, approaching the stationary Gibbs value 0.417.
- **Symmetry:** C_SB(−t) = C_SB(+t) and D(−t) = D(+t) to all printed digits, while J(−t) = −J(+t).
- **Why:** Σ₀ = RΣ₀R (no q–p cross terms) and the dynamics are time-reversal symmetric, so Σ(−t) = RΣ(t)R. C_SB and D
  are R-invariant (even in t); J is R-odd.

Answers to the owner's three questions:
1. **Yes:** t = 0 is the zero-correlation point.
2. **Yes:** conservative evolution generates correlations in **both** time directions.
3. **No:** D does not distinguish + from −. The relaxation is Janus-symmetric.

The epoch is the point where the *declared form* (productness) holds, so it is selected by the preparation, not by D.
(That C_SB(t) > 0 for all t ≠ 0 is shown only at the sampled times.)

## B2-6 S6 orientation firewall (bookkeeping clarification; S6 is not downgraded)

**S6 contains two declared notions of forward. Neither is selected by the time-symmetric dynamics itself** [BR3-02]:
- **J channel:** f_J = +J on L1 (T_s > T_b) and −J on L2 (T_s < T_b). This sign is tied to the declared hot / cold
  ordering.
- **Entropy-reference channel:** σ = −Ḋ, i.e. descent toward the supplied reference S_ref = T_b·diag(r, 1). This direction
  is **reference-relative**, not the L1 / L2 temperature-sign convention.

So NET-ARROW-CONFIRMED means: given the supplied preparation, the J-channel sign convention and the supplied reference
state, the integrated observables keep a net forward value.

| earned (theorem-grade, canonical) | supplied / conventional |
|---|---|
| local return to equilibrium (LS-1) | which preparation is realized (H_sys, H_bath, H_cross) |
| closed-form integrated transfer (LS-2) | the reference state S_ref (= the bath's Gibbs marginal) |
| t⁻⁶ entropy tail with late sign reversal (LS-3) | the J-channel sign convention (L1 / L2); the σ-channel reference direction (descent to S_ref) |
| no erasure on the declared open members | the boundary epoch t = 0 |

**New bridge observation.** At **T_s = T_b** the L1 / L2 convention for J has **no ordering to refer to**, yet the raw J
transfer is non-zero (½Tr²) when the correlation boundary is mismatched. Its *sign* is set by the correlations (B2-2:
negative for α > ½). The entropy-reference construction keeps its supplied reference at T_s = T_b (D(0) = 0.019 > 0). In
neither channel does the convention carry the arrow; the preparation (and, for σ, the supplied reference) does.

**Proposed:** CANONICAL UPDATE CANDIDATE **B2-CUC-1** (bookkeeping clarification only; not applied). Annotate S6-1 [as split by BR3-02]:
"S6 has two declared forward notions: a temperature-ordering sign for J (L1 / L2), and descent toward the supplied
reference S_ref for D. Neither is selected by the time-symmetric dynamics. The equal-temperature offset ½T_b r² is a
correlation (interaction-energy) mismatch of the declared product preparation, and its sign depends on the declared
S–B correlations."

## B2-7 Temperature parameters — **THERMAL STATE CLASS CONDITIONALLY DEFINED; TEMPERATURE VALUE SUPPLIED**

**Members.** The S6 family (2,1), (10,1), (1/2,1), (1/10,1), plus (1,1), (1.05,1) and (1,1.05), all reproduce
X_J(∞) = (T_s − T_b) + ½T_b r² to 2e-5. Every member relaxes locally to T_b·diag(r, 1).

**Stationarity does not select T.** Gibbs is stationary at T = 0.5, 1 and 3 (drift ≤ 2e-15).

**Stationarity does not select Gibbs either.**
- The harmonic chain is integrable.
- For every positive function F(K), the zero-q–p-cross covariance Q = F(K), P = K·F(K) is stationary. Thus the harmonic
  chain admits an infinite non-Gibbs stationary family [BR3-01]. (No claim is made that this family exhausts all
  stationary covariances; cross blocks and spectral degeneracies are not characterized.)
- A non-thermal GGE, F(λ) = λ⁻¹(1 + 0.4·cos 2λ), is admissible (min eigenvalue 0.169) and stationary (drift 8e-16). Its
  retained marginal diag(0.562, 1.0015) is not a Gibbs marginal at any single T.

KMS / Gibbs parameterizes the bath equilibrium class **only after a Gibbs postulate is supplied**. It does not derive the
state or the value of β.

## B2-8 / B2-9 Does the environment carry H_corr? Mechanism vs boundary

**The bath state fixes the local limit.**
- A system Gibbs at T_s = 2 relaxes locally to the bath's state:

  | bath state | S₁(∞) |
  |---|---|
  | Gibbs, T_b = 1 | diag(0.582, 1) |
  | Gibbs, T_b = 1.5 | diag(0.873, 1.5) |
  | restriction of the GGE | diag(0.5622, 1.0015), the GGE marginal |

  So the **bath state determines the local asymptotic state.**
- With the same bath (T_b = 1), four different (H_sys, H_cross) preparations all reach S₁(∞) = diag(0.58211, 1). Their
  transfers differ:

  | preparation | X_J(∞) |
  |---|---|
  | S6 (2, 1) | +1.169 |
  | S6 (2, 1) + α = 1.4 | +0.695 |
  | S6 (0.1, 1) | −0.731 |
  | marginal-matched | +0.339 |

| H component | does layer 5 (environment) determine it? | classification |
|---|---|---|
| H_bath | given a supplied Gibbs postulate **and** T_b, the bath marginal and the local asymptotic state follow | **PARTIAL RELOCATION INTO ENVIRONMENT** (conditional); T_b and the Gibbs-vs-GGE choice stay supplied |
| H_sys | no; it is forgotten locally but retained in X_J | **SUPPLIED** |
| H_cross | no; the same as H_sys | **SUPPLIED** |
| H_epoch | no | **SUPPLIED** (selected by the preparation form, B2-5) |

**B2-9.** The infinite absolutely continuous bath is a **relaxation mechanism** (D / environment structure). It makes
every trace-class perturbation locally forgettable (B2-P1), but it **does not choose which perturbation occurred.**

> **D COMPRESSES MEMORY OF H_corr DOWNSTREAM** in the local state. **The integrated record keeps it**:
> X_J(∞) = E₁(0) − E₁_G + Q12_G − Q12(0). H_corr is not derived.

## B2-10 Noise-layer firewall (not identified with the S6 conservative bath)

**Setup.** The overdamped L0-1e process ẋ = −Kx + ξ, with Q = 2·diag(T_i) and n = 40.

**Uniform T.** The stationary covariance is T·K⁻¹, which is **cross-correlated** (⟨x₁x₂⟩ = 0.33885 = T·r²).

**Non-uniform T_i.** This gives a correlated non-equilibrium steady state (⟨x₁x₂⟩ = 0.177).

**Effect on an initial cross-correlation.** A product-of-marginals preparation is compared with the stationary
preparation. Their difference in ⟨x₁x₂⟩ falls from 0.339 to 1.4e-3 by t = 3.

**Reading.** The noise layer **fixes the stationary correlation structure conditional on T_i**. That is **RELOCATION** of
the reference state into supplied noise parameters, analogous to B1-4 for Σ. It **erases** the initial H_cross
exponentially and does **not select** which initial correlation was prepared.

## B2-11 Entropy firewall

| comparison | entropy | behaviour |
|---|---|---|
| equal-T factorized A vs Gibbs B | S(A) − S(B) = −0.1459 = −½ log(2.3r) | A relaxes, B stationary |
| **marginal-matched M vs Gibbs B** | **S(M) − S(B) = +0.1459** (= the Gibbs S–B mutual information = ½ log(2.3r), an exact identity via the Schur complement) | **the higher-entropy M relaxes**; B is stationary |
| forward-evolved Σ(τ) vs reversed Σ_R | **identical** global and marginal entropies | **opposite** local arrows |

**Reading.** Low entropy is neither the precise nor a sufficient label. For canonical S6 the exact supported claim is
narrower: **equal marginal temperature does not remove the transient when interaction correlations are absent.** Across
both programs the operative datum is **correlation structure relative to the split**. For reversal, that includes the
sign of the q–p correlations.

## B2-12 Residual sharpening (notation not forced)

    H_corr|Σ  →  H_boundary|Σ = (H_sys, H_bath, H_cross) @ H_epoch

**Independence:**
- **H_cross vs the marginals:** independent (B2-2: cross varied, marginals fixed).
- **The marginals vs H_cross:** independent (B2-7: marginals varied, cross fixed at 0).
- **H_bath:** partly conditionally structured by the environment (B2-8).
- **H_epoch is *not* independent of the others up to time-translation gauge.** (Σ₀, t₀) and (Σ(τ), t₀ + τ) describe the
  same history. H_epoch is meaningful only as "the epoch at which a declared *form* (productness) holds" (B2-5).

**Bridge notation:**

    [Σ ⊗ H_corr]_coupled  =  [Σ ⊗ (H_marginals, H_cross)]  at the declared-form epoch

This is **conceptual sharpening, not compression.**

## B2-13 Elimination test — **H_corr REMAINS SUPPLIED**

Fix D, Σ, A and the environment law (bath Gibbs at T_b = 1). The following admissible preparations all exist, behave
differently, and are not related by any gauge:
- S6 (2, 1);
- S6 (0.1, 1);
- the S6 marginals with α = 1.4;
- the marginal-matched product;
- the graft;
- the reversed state;
- the Gibbs state itself.

Their X_J(∞) values range from −0.73 to +1.17 (and −1.14 over [0, τ] for the reversed state). Their transients differ,
and some relax while others do not. No canonical GRUT layer reconstructs the realized member.

## B2-14 Prediction firewall

No new prediction is claimed. The generalized relation X_J(∞) = E₁(0) − E₁_G + Q12_G − Q12(0) is an exact consequence of
the model given a declared preparation, **not** a parameter-free observable. Baseline: **zero confirmed distinctive GRUT
quantitative predictions.** Only B5 can test any removed freedom.

## Scope and limits

- **Model scope:** one canonical model, the S6 infinite pinned chain (classical, Gaussian, harmonic), with finite-N
  numerics and an infinite-N bridge theorem (B2-P1, reviewed at stated scope; local covariance convergence only, not
  global trace-norm convergence).
- **Not covered:** non-Gaussian and nonlinear cases, and the quantum version (SCOUT's domain).
- **Verification:** the proposition received an independent logic review (owner) at its stated scope; grade
  **BRIDGE-THEOREM REVIEWED AT STATED SCOPE**, not a canonical GRUT theorem.
