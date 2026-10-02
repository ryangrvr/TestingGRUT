# QFT-SCOUT-1 · G2 — MODULAR FLOW, H_epoch AND ORIENTATION (result)

**Question (owner).** Given a standard pair (M, Ω) or (M, ω), Tomita–Takesaki supplies a canonical modular automorphism
group. Does this remove any frozen time / orientation information, or merely relocate it into the chosen algebra and
state?

**Script:** `g2/g2_modular.py` (+ `.log`).
- **Finite / type-I controls:** exact standard-form construction for (B(ℂ⁴), ρ).
- **Lattice illustration:** an 80-digit-precision Gaussian modular Hamiltonian for a free-field interval. This is an
  **illustration, not proof**.
- **Label:** independent code path, not independent reviewer.

## 0. Mandatory separations (charter firewall 3)

| level | statement | status |
|---|---|---|
| **(1)** modular flow exists | standard (M, Ω) ⇒ Δ, J; Δ^{it}MΔ^{−it} = M; JMJ = M′; ω is KMS at β = −1 for σ_t | **theorem** (Tomita–Takesaki; Takesaki 1970) |
| **(2)** modular flow is geometric | wedge algebra + vacuum + Poincaré covariance + spectrum condition ⇒ Δ^{it} = U(Λ_W(−2πt)) (boosts), J = CPT·rotation | **theorem under extra hypotheses** (Bisognano–Wichmann 1975 / 1976) |
| **(3)** modular flow is physical time | the thermal-time hypothesis | **interpretive hypothesis** (Connes–Rovelli 1994, gr-qc/9406019); **not adopted** |
| **(4)** orientation / sign | tested separately below | — |

**Source grades.**
- (1) and the imports below are **KNOWN-RESULT IMPORT — CITED FROM THE STANDARD LITERATURE, TEXT NOT RE-READ.**
  Primary-source fetches are blocked here; these are textbook / landmark results, and upgrading them is owner-reviewable.
- (3) uses the owner's pointer.

## G2-1 Finite / type-I control — **modular flow exists with no spacetime or physical-time content**

For (B(ℂ⁴), ρ) in standard form on ℂ⁴ ⊗ ℂ⁴:

| check | error |
|---|---|
| S = JΔ^{1/2} | 1.0e-14 |
| ΔΩ = Ω | 3.7e-15 |
| σ_t(A) = ρ^{it}Aρ^{−it} on M | 8.5e-14 (an **inner** flow) |
| Tomita-convention KMS ω(Aσ_{−i}(B)) = ω(BA) | 1.5e-14 |

So ω is KMS at **β = −1 for σ_t**, i.e. at **β = +1 for τ_t := σ_{−t}**.

**This is the baseline:** a canonical one-parameter group with no geometric or physical-time meaning.

## G2-2 State dependence — **type I: STATE-PRICED; type III₁: an injective outer modular homomorphism (CONDITIONALLY SELECTED)** [Q2R-01]

**Type I.**
- Different faithful states give different modular Hamiltonians (‖K₁ − K₂‖ = 3.73).
- **Every** Hermitian K is the modular Hamiltonian of ρ_K = e^{−K}/Z (error 3e-15).
- Flows of two states differ by the Connes cocycle u_t = ρ₂^{it}ρ₁^{−it} ∈ M (error 1e-14).
- **In type I the modular flow carries exactly the state's information: STATE-PRICED / RELOCATION.**

**Beyond type I** (imports):
- **Connes cocycle theorem (1973).** For any factor M, the modular automorphism groups of faithful normal states define
  one **canonical, state-independent modular homomorphism δ_M : ℝ → Out(M)**, whose kernel is Connes' T(M).
- **Semifinite factors (types I, II):** modular automorphisms are inner, so δ_M is trivial.
- **Type III₁:** T(M) = {0}, so δ_M is **injective / faithful.**
- **Other type III classes** (III_λ, III₀) can also have non-trivial outer modular structure, but **not** the III₁
  injectivity property.
- **Physical local algebras** [Q2R-02]: under the standard phase-space / nuclearity / scaling assumptions used in the
  canonical AQFT structural results (Buchholz–D'Antoni–Fredenhagen 1987), physically relevant local algebras fall into the
  **hyperfinite (injective) III₁ class.** Haagerup's uniqueness theorem (1987) identifies the injective III₁ factor
  uniquely up to isomorphism.

**New G2 structure, stated narrowly:** **TYPE III₁ CONDITIONALLY SELECTS AN INJECTIVE, STATE-INDEPENDENT OUTER MODULAR
HOMOMORPHISM δ_M.** It does **not** select a specific representative σ_t (those remain state-priced through Connes
cocycles), physical time, or orientation. Type III₁ does not uniquely create the idea of an outer modular structure; it
makes δ_M injective. This differs from G1's property: G1 needed only non-type-I, while G2's injectivity needs III₁.

**Limits.**
- δ_M is an outer **class**, not a specific automorphism group. Picking the group needs a state (STATE-PRICED).
- For a local algebra 𝒜(O), physical inertial time translations do not even preserve 𝒜(O). δ_M relates to physical
  dynamics only through (2) (wedges: boost / Rindler time) or (3) (an interpretation).

> **Classification:** the injective outer modular homomorphism is **CONDITIONALLY SELECTED** (given M type III₁). As **dynamics**
> it is interpretation-priced. The specific flow is STATE-PRICED.

## G2-3 Sign inversion — **ORIENTATION NOT SELECTED (CONVENTION- / Σ-PRICED)**

| operation | effect on the sign of modular time | evidence |
|---|---|---|
| algebra ↔ commutant (M ↔ M′) | **reverses**: Δ′ = Δ^{−1}. The commutant's flow is 1 ⊗ ρᵀ^{−it}Bρᵀ^{it} | finite check, error 9e-14. For wedges, the right and left wedges' modular flows are opposite boosts |
| modular conjugation J | commutes with Δ^{it} (antilinear) and **maps M → M′**: it carries the M-flow to the M′-flow | error 3e-14 |
| anti-automorphism α of M [Q2R-03] | σ^{φ∘α}_t = α^{−1}∘σ^φ_{−t}∘α: **reverses** | finite check with α = transpose, error 2e-15. For hyperfinite III₁, M ≅ M^op (by uniqueness), so δ_M is carried to its reverse by an anti-automorphism of the **same** algebra |
| Tomita parameter convention (Δ^{it} vs Δ^{−it}) | a definitional sign | convention |
| replacing the state | changes the inner flow, not the outer class or its sign | G2-2 |

**Reading.**
- The existence of σ_t is canonical; its **orientation is not an invariant of M up to anti-isomorphism.**
- The orientation is fixed only by:
  - declaring which algebra is "the system" (a Σ / A_partition choice);
  - the Tomita / KMS sign convention;
  - extra physical input (G2-4 / G2-5).
- **Scope [Q2R-02/03]:** the opposite algebra of an injective III₁ factor is again injective III₁, so Haagerup uniqueness
  gives R_∞ ≅ R_∞^op: an anti-isomorphism exists **for the hyperfinite III₁ class**. For that class, modular orientation
  is not invariant under algebraic structure considered up to anti-isomorphism. **This is not claimed for every factor
  or for arbitrary type III₁ factors.** Factors not anti-isomorphic to themselves (Connes 1975) are the control showing
  that anti-isomorphism cannot be assumed universally.

## G2-4 Bisognano–Wichmann — **CONSISTENCY / RELOCATION (orientation into the spectrum condition)**

**What is newly earned.** For (wedge algebra, vacuum, Poincaré covariance, spectrum condition), the state-defined modular
flow **coincides with** the supplied geometric boost flow, and J with CPT × rotation.

**Accounting.**
- **Nothing geometric is reconstructed.** Poincaré covariance is an input. The theorem is a consistency constraint
  linking the state to the supplied geometry.
- **The geometric sign / orientation in the BW identification is priced by the supplied relativistic structure** [Q2R-04]:
  - the Poincaré representation;
  - the wedge choice;
  - the spectrum / future-cone condition;
  - the modular sign convention.

  The modular flow is the boost by −2πt relative to these. The left wedge gets the reverse.

**Lattice illustration** (80 digits; vacuum of the free field with μ = ma = 0.02; a 30-site interval).
- Near the cut the modular Hamiltonian is boost-like: H_p[j,j] / 2π(j + ½) = **0.996, 0.976, 0.952** at j = 0, 1, 2.
- Across the interval it tracks the conformal two-ended profile 2πx(ℓ − x)/ℓ to within 1–16%. That profile is geometric
  only for conformal theories (Hislop–Longo); the mass causes the deviations.
- **A thermal state (β = 4) on the same algebra is not boost-like:** the ratios are 0.706, 0.394, 0.252, ….
- So the geometric identification belongs to **vacuum + covariance + spectrum condition**, not to the algebra alone.

> **Classification:** BW earns the **consistency relation "vacuum modular flow = wedge boost flow"** and does not by
> itself reconstruct the prior spacetime structure. The sign is priced by the supplied relativistic structure (Poincaré
> representation, wedge, spectrum condition, modular convention) [Q2R-04]. G3 tests the inverse direction.

## G2-5 KMS, passivity and orientation — **CONVENTION-PRICED / RELOCATION into passivity**

The four notions are separated as follows:

| notion | status |
|---|---|
| **Modular KMS** | the theorem: β = −1 for σ_t **in the code's modular-parameter convention** σ_t(A) = ρ^{it}Aρ^{−it}, i.e. β = +1 for σ_{−t}. **This sign is a modular-parameter convention, not physical evidence for either temporal orientation** [Q2R-06] |
| **Positive physical temperature** | the *declaration* that physical time τ_t = σ_{−βt} with β > 0 |
| **Passivity** | theorem-grade thermodynamics (Pusz–Woronowicz 1978) concerns **complete passivity** and KMS / ground states |
| **Orientation** | not fixed by any of the above without a further postulate |

**Finite check.** For H = +K (τ_t = σ_{−t}), ρ is passive: the maximum extractable cyclic work is **0.0000** (random
search and exact). For the reversed orientation H = −K, ρ is maximally active (**2.19** exact).

**Finite control:** for the chosen faithful state, H = +K is passive and H = −K is active. So **imposing** passivity
distinguishes the two candidate signs **in this control** [Q2R-05].

**Theorem-grade reading:** passivity / complete passivity can orient a dynamics relative to a state **once the Kelvin-type
thermodynamic condition is imposed.** This relocates the orientation into the thermodynamic postulate; it does not derive
that postulate. **"Passivity derives / selects the arrow" is not claimed. β > 0 is not a derived arrow.**

## G2-6 Split-buffer carryover — **RELOCATION into (d, 𝒩, ω)**

**Setup.** For the G1 buffered case (O₁ ⋐ O₂), the intermediate type I factor 𝒩 has a canonical choice once the vacuum
is given (Doplicher–Longo; bibliographic only).

**Consequences.**
- 𝒩 is **type I**, so its modular flow is **inner** and δ_𝒩 is trivial. 𝒩 carries **no intrinsic outer time.**
- Its modular Hamiltonian −log ρ_𝒩 exists but is state-priced (G2-2).
- **The buffer scale d is an input**: O₁ ⋐ O₂ is chosen before 𝒩 is constructed.
- **The special-form boundary event H_epoch is not fixed.** A one-parameter group has no distinguished origin, and no
  modular datum picks the instant at which the split product preparation is realized.

> **Classification: RELOCATION / CONDITIONAL.** (d, 𝒩, ω) must first be supplied.

## Classification against the frozen residual

| frozen item | G2 result | class |
|---|---|---|
| **time orientation** | flips under M ↔ M′, under anti-automorphisms (R_∞ ≅ R_∞^op for the hyperfinite III₁ class [Q2R-02/03]), and under the modular-parameter convention. Fixed only by Σ-choice, convention, the supplied relativistic structure (BW) or an imposed passivity criterion | **NOT SELECTED: CONVENTION-PRICED / RELOCATION** |
| **H_epoch** (special-form event) | no modular datum picks an origin; with a split buffer, (d, 𝒩, ω) are supplied | **RELOCATION / not selected** |
| **D_dyn (time-flow structure)** | type I: the flow is STATE-PRICED. **Type III₁: the canonical state-independent modular homomorphism δ_M is injective** (Connes T(M) = {0}), but it is physical dynamics only via BW (wedge boost time) or the thermal-time interpretation | injective δ_M **CONDITIONALLY SELECTED** (given M type III₁); as dynamics, interpretation-priced |
| **Lorentz / boost structure** | BW identifies the modular flow with supplied boosts | **CONSISTENCY** (a constraint), no reconstruction. The modular → geometry direction is G3 |

> **G2 TERMINAL (after QFT REPAIR 02):**
> - **ORIENTATION NOT SELECTED.**
> - **H_epoch NOT SELECTED.**
> - A type III₁ factor conditionally supplies an **injective, state-independent outer modular homomorphism δ_M**. Neither
>   a particular inner representative, a physical-time interpretation, nor orientation follows from that fact alone.
> - BW gives a geometric consistency relation under supplied relativistic hypotheses.
> - Thermodynamic orientation requires a supplied passivity / KMS criterion.
> - **TRUE COMPRESSION = 0.**

## Scope

- The finite checks are exact algebra. The lattice is one free-field illustration.
- The Connes / Haagerup / BDF / BW / Pusz–Woronowicz / Connes-1975 imports are cited from the standard literature
  (text not re-read).
- The σ^{φ∘α} reversal formula is verified for the finite transpose. For general anti-automorphisms it is an import /
  standard computation.
- Auxiliary to canonical GRUT (QP-5).
- No payoff claim; G4 has not been run.
