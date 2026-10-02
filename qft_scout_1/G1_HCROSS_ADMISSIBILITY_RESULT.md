# QFT-SCOUT-1 · G1 — IS H_cross = 0 ADMISSIBLE? (result)

**Target:** the frozen component **H_cross | Σ**. Can a state with zero system / complement correlations, i.e. a
**normal product state**, be realized for a given localization?

**Chartered correction (owner).** "Type III ⇒ no product state" is **not** claimed. The question is the existence of a
**normal** product state for a **chosen sharp complementary localization**, contrasted with split-buffer localization.

## 1. Theorem chain (as repaired by QFT REPAIR 01; see `QFT_CORRECTION_LEDGER.md`)

| ID | statement (hypotheses explicit) | role in G1 | grade |
|---|---|---|---|
| **S1 Uncorrelated-state criterion** | For von Neumann algebras A, B with A ∨ B a factor: **(A, B) is split iff there exists a normal A–B-uncorrelated state on A ∨ B.** For commuting A, B, a product state φ(ab) = φ₁(a)φ₂(b) is uncorrelated | **load-bearing** | KNOWN-RESULT IMPORT — PRIMARY / SOURCE-TEXT VERIFIED (Buchholz & Summers, *Quantum statistics and locality*, Phys. Lett. A 337 (2005) 17–21; Summers 2009, Thm. 5.4) [Q1R-06] |
| **S2 Split ⇒ tensor product** | If A ⊂ 𝒩 ⊂ B′ with 𝒩 a type I factor, then A ∨ B ≅ A ⊗̄ B spatially, and **arbitrary normal marginal states extend to normal product states** | buffered case | KNOWN-RESULT IMPORT — SOURCE-TEXT VERIFIED (Summers 2009, Thm. 4.1 and the product-extension construction). Doplicher–Longo, Invent. Math. 75 (1984) 493 is **bibliographically verified only** (paywalled body not inspected) |
| **L1 (proved here)** | If (ℳ, ℳ′) is split for a factor ℳ, i.e. ℳ ⊂ 𝒩 ⊂ (ℳ′)′ = ℳ with 𝒩 type I, then 𝒩 = ℳ, so **ℳ is type I.** *Proof:* bicommutant theorem. ∎ | load-bearing | proved (one line) |
| **G1-P1 Sharp normal-product obstruction** | **Let ℳ ⊂ B(H) be a non-type-I factor. Then the sharp commuting pair (ℳ, ℳ′) admits no normal product state.** *Proof:* (1) ℳ a factor ⇒ ℳ ∨ ℳ′ = B(H), a factor. (2) A normal product state is uncorrelated. (3) Normal uncorrelated state + factor join ⇒ split (S1). (4) Split (ℳ, ℳ′) ⇒ ℳ type I (L1). (5) Contradiction. ∎ **Geometric form:** for an AQFT region O with Haag duality 𝒜(O′) = 𝒜(O)′, if 𝒜(O) is a non-type-I factor, the same holds for the geometric sharp pair (𝒜(O), 𝒜(O′)) | **the G1 core** | **QFT-SCOUT ASSEMBLY OF SOURCE-TEXT-VERIFIED THEOREMS** |
| **S3 Non-normal products (sharp pair)** | For a type II or type III factor ℳ with N = ℳ′: normal product states are impossible; **product extensions of normal marginals can exist only non-normally (singularly)**; for type I, normal product states exist | control | KNOWN-RESULT IMPORT — SOURCE-TEXT VERIFIED (Rédei & Summers, as checked by the owner) |
| T2 Type (physical realization) | Local algebras of double cones / wedges are type III₁ factors under standard hypotheses: **one physically relevant route to non-type-I**, not the mathematical reason | realization | KNOWN-RESULT IMPORT — SOURCE LOCATED (Araki 1964; Fredenhagen 1985; Buchholz–D'Antoni–Fredenhagen 1987) |
| T3 Haag duality | 𝒜(O′) = 𝒜(O)′ (wedges; free-field double cones) | **geometric identification step only** [Q1R-03] | KNOWN-RESULT IMPORT — SOURCE LOCATED (Araki 1963; Bisognano–Wichmann 1975/76) |
| T1 Reeh–Schlieder | The vacuum is cyclic and separating for 𝒜(O) under the standard hypotheses | **not load-bearing for G1** [Q1R-04]; retained for vacuum structure and G2 | KNOWN-RESULT IMPORT — SOURCE LOCATED (Reeh & Schlieder 1961) |

**Load-bearing G1 inputs (priced):** **normality** (relative to the chosen representation) + **factoriality** +
**non-type-I** + **sharp complement (ℳ, ℳ′)**. For the spacetime reading, add **Haag duality**.

**Not required:** faithful marginals, standardness, Reeh–Schlieder, type III specifically [Q1R-01, Q1R-02, Q1R-04].
- **Type III is stronger than necessary.** The obstruction also covers type II factors. AQFT's type III local algebras are
  the physically relevant realization of the general non-type-I obstruction.
- **This is not a claim about all infinite-dimensional algebras.** Infinite-dimensional type I algebras (B(H) ⊗ B(K))
  remain the control: they admit normal product states.

**Normality is the load-bearing state-class restriction** [Q1R-05]. Normality is understood relative to the chosen
representation. Identifying it with "finite energy" or another physical property needs separate assumptions, which are
not made here.

## 2. Controls

| class | H_cross = 0 (product state) | evidence |
|---|---|---|
| **Finite or infinite type-I control** (B(H₁) ⊗ B(H₂); any lattice at fixed spacing) | **always admissible and normal** (ρ₁ ⊗ ρ₂) | elementary; consistent with L1, since type I sharp pairs split |
| **Sharp non-type-I factor / commutant** (algebraic; geometric with Haag duality) | **inadmissible in the normal class** (G1-P1) | S1, L1 (+ T3 for geometry) |
| **Split-buffer localization** (A ⊂ 𝒩 ⊂ B′, 𝒩 type I) | **admissible**; arbitrary normal marginals extend to normal product states | S2 |
| Non-normal state class (type II / III sharp pair) | product extensions of normal marginals exist **only non-normally** | S3 |

**Lattice illustration** (`g1/g1_hcross_lattice.py` + `.log`; free massive scalar in 1+1D, m = 1, box [−6, 6], a = 0.05 →
0.003125). This is **not proof.**

- **Sharp cut.**
  - The half-box entropy S rises by **0.1663 → 0.1667 per halving of a** (→ c/6 for c = 1), with S = 0.4997 → 0.9614,
    tracking (1/6)·ln(1/(ma)).
  - The mutual information I = 2S is the **minimum relative entropy from the vacuum to *any* product state**
    (min over σ_A ⊗ σ_B of S(ω‖σ_A ⊗ σ_B) = I(A:B)). It diverges like (1/3)·ln(1/a).
  - The energy cost of the product of vacuum marginals is 9.8 → 297.8, a ratio of about 2.3 per halving: roughly
    ln(1/a)/a.
- **Buffered cuts.** I(A:B) converges as a → 0. For d = 0.5, 1.0 and 2.0 it reaches **0.058398, 0.014897 and 0.001412**;
  the increments shrink to ≤ 4e-7.
- **The vacuum is correlated across every tested split** (I > 0), so it is not itself a product state. This is
  consistent with T1, and is a statement about one state, not about the class (T1 is not load-bearing for G1-P1).

**What non-type-I structure buys over the type-I control:** exactly the **sharp-split obstruction** (the AQFT type III
local algebras being the physical realization). At any
finite spacing the product state is admissible. Its minimal relative-entropy cost and its energy cost diverge in the
continuum, and the theorem chain makes that an inadmissibility statement in the normal class. A buffer removes the
obstruction.

## 3. Classification against the frozen residual

| item | result | class |
|---|---|---|
| **H_cross = 0 for a sharp complementary split Σ** | inadmissible in the normal-state class | **FORBIDDEN (in class)** |
| **H_cross overall** | the admissible set loses its product member for sharp splits; correlated normal states remain (the vacuum, its local excitations, …). **The realized correlated state is not selected** | **CONSTRAINED-NONUNIQUE** |
| **H_cross = 0 with a split buffer** | admissible; no constraint | no change |
| **Σ ⊗ H coupling** | **strengthened, not compressed.** Whether "zero correlation" is possible now depends on the **geometry of the split** (sharp vs buffered) | coupling made structural (CONDITIONAL, priced) |
| **New priced items exposed** | normality; factoriality; non-type-I; sharp complement; Haag duality (geometric reading only); for split preparations, the **buffer scale d** and an intermediate type I factor 𝒩. 𝒩 is not unique in general, though there is a canonical choice given the vacuum (Doplicher–Longo, import) | supplied / declarative accounting items (BR4-01 discipline: not automatically physical primitives) |

**Relevance to the Bridge-1 arrow finding.** The product ("fresh, independent") preparation that powers the S6 transient
and SCOUT's fresh-bath arrow (D4) is **not available as a sharp normal state in the AQFT envelope.** A fresh boundary
needs a buffer, i.e. a split scale. There is no contradiction with Bridge-1: S6 is a fixed-spacing, classical lattice
chain, which corresponds to the type-I control column. But the special-form boundary event H_epoch, read in this
envelope, must carry localization data (d, 𝒩), not just "product at t = 0".

> **G1 TERMINAL (after QFT REPAIR 01):** **H_cross = 0 INADMISSIBLE IN THE SHARP NORMAL-STATE CLASS**
> (FORBIDDEN-in-class, via G1-P1: non-type-I factor + sharp complement + normality; Haag duality for the geometric
> reading). **ADMISSIBLE WITH A SPLIT BUFFER.** Overall **H_cross: CONSTRAINED-NONUNIQUE.** Not selected; **not TRUE
> COMPRESSION.**
>
> **First genuine admissibility restriction on the frozen residual:**
> - finite / type-I envelope: H_cross = 0 allowed;
> - sharp non-type-I normal-state envelope: H_cross = 0 forbidden.
>
> This is a class restriction, not RENAMING. **The admissible H_cross depends on the localization class.**

## 4. Scope and limits

- **G1-P1** is an assembly of source-text-verified theorems (S1, owner-verified) and the one-line L1. It needs neither
  faithful marginals, standardness nor Reeh–Schlieder [QFT REPAIR 01].
- **The lattice** illustrates one model: a free scalar in 1+1D. The relative-entropy bound covers *all* product states;
  the energy figure covers only the product of vacuum marginals.
- **Bounded-energy and general-state entanglement claims** are not used. They are a separate import, reserved for later.
- **Auxiliary to canonical GRUT** (QP-5: the lift is supplied; canonical NONUNIQUE-LIFT stands).
- **No payoff claim.** G4 has not been run, and the baseline stays at zero confirmed distinctive GRUT predictions.
