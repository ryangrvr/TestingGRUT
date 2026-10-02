# QFT-SCOUT-1 · G1 — IS H_cross = 0 ADMISSIBLE? (result)

**Target:** the frozen component **H_cross | Σ**. Can a state with zero system / complement correlations, i.e. a
**normal product state**, be realized for a given localization?

**Chartered correction (owner).** "Type III ⇒ no product state" is **not** claimed. The question is the existence of a
**normal** product state for a **chosen sharp complementary localization**, contrasted with split-buffer localization.

## 1. Theorem chain (imports graded; QFT-SCOUT propositions marked)

| ID | statement (hypotheses explicit) | grade |
|---|---|---|
| **T1 Reeh–Schlieder** | In a Haag–Kastler / Wightman QFT with the spectrum condition and weak additivity, the vacuum Ω is cyclic and separating for 𝒜(O) whenever O has a non-empty causal complement. *(Bounded-energy-state extensions are NOT bundled here; they are a separate import, not used in G1)* | KNOWN-RESULT IMPORT — SOURCE LOCATED, TEXT NOT RE-READ (Reeh & Schlieder, Nuovo Cimento 22 (1961) 1051; owner pointer CMP 2026 for bounded-energy refinements) |
| **T2 Type** | Under standard hypotheses (e.g. free fields; more generally the scaling-limit / modular conditions), local algebras of double cones and wedges are factors of type III₁ | KNOWN-RESULT IMPORT — SOURCE LOCATED (Araki 1964 for the free field; Fredenhagen 1985; Buchholz–D'Antoni–Fredenhagen 1987) |
| **T3 Haag duality** | 𝒜(O′)′ = 𝒜(O) for wedges (Bisognano–Wichmann setting) and for double cones of the free field | KNOWN-RESULT IMPORT — SOURCE LOCATED (Araki 1963; Bisognano–Wichmann 1975 / 1976) |
| **T4 Split ⇔ normal product states** | For regions O₁ ⋐ O₂ (a buffer), the split property, i.e. a type I factor 𝒩 with 𝒜(O₁) ⊂ 𝒩 ⊂ 𝒜(O₂), gives 𝒜(O₁) ∨ 𝒜(O₂)′ ≅ 𝒜(O₁) ⊗̄ 𝒜(O₂)′. **Any normal states on the two algebras extend to a normal product state.** Conversely, for standard pairs the existence of a normal product state (faithful marginals) implies split. The free field and theories satisfying nuclearity have the split property for strictly separated regions | KNOWN-RESULT IMPORT — SOURCE LOCATED (Doplicher–Longo, Invent. Math. 75 (1984) 493; Buchholz 1974; Buchholz–Wichmann 1986; review: Summers 2009, arXiv:0812.1517; owner pointer CMP 2024 for split-property scope) |
| **L1 (QFT-SCOUT, elementary, proved here)** | Let ℳ be a factor and take the **sharp complementary pair** (ℳ, ℳ′), which Haag duality T3 identifies with (𝒜(O), 𝒜(O′)). If (ℳ, ℳ′) is split, i.e. some type I factor 𝒩 has ℳ ⊂ 𝒩 ⊂ (ℳ′)′ = ℳ, then 𝒩 = ℳ, so **ℳ is type I.** Contrapositive: **a non-type-I local factor with Haag duality gives a non-split sharp complementary pair.** *Proof:* (ℳ′)′ = ℳ by the bicommutant theorem, so the chain collapses. ∎ | **proved** (one line) |
| **C1 (combination)** | T2 + T3 + L1 + the converse direction of T4 imply: **for a sharp complementary localization with Haag duality and a type III local factor, there is no normal product state with faithful marginals on 𝒜(O) ∨ 𝒜(O′) = B(H).** | **QFT-SCOUT PROPOSITION, conditional on the T4-converse import.** Not independently peer reviewed |

**Hypotheses the conclusion actually uses (priced):**
- **normality:** the state class;
- **factoriality** of 𝒜(O);
- **Haag duality** for O (the sharpness of the complementary split);
- **non-type-I** (type III) of 𝒜(O);
- **faithful marginals / standardness** (via T4).

**Type III alone is not used and is not sufficient.** Without Haag duality, 𝒜(O′) can be strictly smaller than 𝒜(O)′
and the pair can split.

**What survives without normality.** Non-normal ("singular") product states on the C*-level always exist for
C*-independent commuting algebras. H_cross = 0 is therefore admissible **only if normality (locally finite energy /
regularity) is dropped.** Normality is the load-bearing price.

## 2. Controls

| class | H_cross = 0 (product state) | evidence |
|---|---|---|
| **Finite / type-I control** (B(H₁) ⊗ B(H₂); any lattice at fixed spacing) | **always admissible and normal** (ρ₁ ⊗ ρ₂) | elementary; L1 is consistent, since type I pairs can split |
| **Sharp AQFT localization** (Haag duality, type III) | **inadmissible in the normal class** (C1) | T1 – T4, L1 |
| **Split-buffer AQFT localization** (O₁ ⋐ O₂) | **admissible** (normal product states exist) | T4 |
| Non-normal state class | admissible | C*-level product states |

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
  consistent with T1, and is a statement about one state, not about the class.

**What infinite / type-III structure buys over the type-I control:** exactly the **sharp-split obstruction**. At any
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
| **New priced items exposed** | normality; Haag duality; type; for split preparations, the **buffer scale d** and an intermediate type I factor 𝒩. 𝒩 is not unique in general, though there is a canonical choice given the vacuum (Doplicher–Longo, import) | supplied / declarative accounting items (BR4-01 discipline: not automatically physical primitives) |

**Relevance to the Bridge-1 arrow finding.** The product ("fresh, independent") preparation that powers the S6 transient
and SCOUT's fresh-bath arrow (D4) is **not available as a sharp normal state in the AQFT envelope.** A fresh boundary
needs a buffer, i.e. a split scale. There is no contradiction with Bridge-1: S6 is a fixed-spacing, classical lattice
chain, which corresponds to the type-I control column. But the special-form boundary event H_epoch, read in this
envelope, must carry localization data (d, 𝒩), not just "product at t = 0".

> **G1 TERMINAL:** **H_cross = 0 INADMISSIBLE IN THE SHARP NORMAL-STATE CLASS** (FORBIDDEN-in-class, via
> Haag duality + type III + normality + the split ⇔ normal-product-state import). **ADMISSIBLE WITH A SPLIT BUFFER.**
> Overall **H_cross: CONSTRAINED-NONUNIQUE.** Not selected; **not TRUE COMPRESSION.**

## 4. Scope and limits

- **C1** depends on the T4-converse import. Its precise hypotheses (faithfulness of the marginals, standardness) are
  carried as stated; the primary text was **not re-read** (egress blocked).
- **The lattice** illustrates one model: a free scalar in 1+1D. The relative-entropy bound covers *all* product states;
  the energy figure covers only the product of vacuum marginals.
- **Bounded-energy and general-state entanglement claims** are not used. They are a separate import, reserved for later.
- **Auxiliary to canonical GRUT** (QP-5: the lift is supplied; canonical NONUNIQUE-LIFT stands).
- **No payoff claim.** G4 has not been run, and the baseline stays at zero confirmed distinctive GRUT predictions.
