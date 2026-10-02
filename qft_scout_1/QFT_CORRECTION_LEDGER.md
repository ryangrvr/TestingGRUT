# QFT-SCOUT-1 CORRECTION LEDGER

## QFT REPAIR 01 — sharpen the G1 theorem chain (owner, after `91401f8`; G1 accepted provisionally)

| ID | change |
|---|---|
| **Q1R-01** | C1 (which used the faithful-product converse) is replaced by **G1-P1**: for von Neumann algebras with a factor join, split ⇔ a normal uncorrelated state exists (Buchholz–Summers 2005; Summers 2009, Thm. 5.4). For a factor ℳ, ℳ ∨ ℳ′ = B(H), so a single normal product state would make (ℳ, ℳ′) split, and by L1 ℳ would be type I. **Faithful marginals, standardness and Reeh–Schlieder are removed** from the load-bearing hypotheses |
| **Q1R-02** | The mathematical core is **non-type-I factor**, not "type III factor". AQFT type III local algebras are the physically relevant realization. The obstruction also covers type II. It is **not** claimed for all infinite-dimensional algebras: type I remains the control |
| **Q1R-03** | **Haag duality scope.** The algebraic statement for the exact pair (ℳ, ℳ′) needs no Haag duality. The geometric AQFT statement needs Haag duality to identify 𝒜(O′) = 𝒜(O)′. It is priced only at that step |
| **Q1R-04** | **Reeh–Schlieder** is moved out of the G1 dependency chain. It is retained for vacuum structure and G2 |
| **Q1R-05** | **Non-normal products.** For a type II / III factor and its commutant, normal product states are impossible and product extensions of normal marginals exist **only non-normally / singularly** (Rédei–Summers, owner-checked). Normality is the load-bearing state-class restriction, **relative to the chosen representation**. It is **not** identified with finite energy without separate assumptions |
| **Q1R-06** | **Source grades.** Sharp normal-product obstruction: **PRIMARY / SOURCE-TEXT VERIFIED** (Buchholz & Summers, Phys. Lett. A 337 (2005) 17–21; Summers 2009, Thm. 5.4). Split ⇒ tensor product / normal product extensions: **SOURCE-TEXT VERIFIED** (Summers 2009, Thm. 4.1). Doplicher–Longo 1984: **bibliographically verified only**, not primary-text verified |

**G1-P1 grade:** QFT-SCOUT ASSEMBLY OF SOURCE-TEXT-VERIFIED THEOREMS.

**Terminal preserved:**
- H_cross = 0 INADMISSIBLE IN THE SHARP NORMAL-STATE CLASS;
- ADMISSIBLE WITH A SPLIT BUFFER;
- H_cross CONSTRAINED-NONUNIQUE;
- TRUE COMPRESSION not earned.

## QFT REPAIR 02 — sharpen G2 scope (owner, after `8f70a7a`; G2 accepted provisionally; terminal unchanged)

| ID | change |
|---|---|
| **Q2R-01** | **Modular-homomorphism scope.** For a factor M, the modular groups of faithful normal states define the canonical state-independent δ_M : ℝ → Out(M), with kernel Connes' T(M). Semifinite: δ_M trivial. **Type III₁: T(M) = {0}, δ_M injective.** Other type III classes can have non-trivial but non-injective outer modular structure. New G2 structure: **TYPE III₁ CONDITIONALLY SELECTS AN INJECTIVE STATE-INDEPENDENT OUTER MODULAR HOMOMORPHISM**. It selects no representative σ_t (state-priced via cocycles), no physical time and no orientation. The language that "type III₁ uniquely creates an intrinsic outer flow" is removed |
| **Q2R-02** | **Hyperfinite III₁ scope.** No unqualified "AQFT local algebras are the unique hyperfinite III₁ factor". Use: under the standard phase-space / nuclearity / scaling assumptions, physical local algebras fall into the hyperfinite (injective) III₁ class; Haagerup's uniqueness then identifies it up to isomorphism. The opposite of injective III₁ is injective III₁, so R_∞ ≅ R_∞^op. Not generalized to arbitrary III₁ factors |
| **Q2R-03** | **Orientation under anti-isomorphism.** The finite transpose calculation is kept. The structural conclusion is stated **only for the hyperfinite III₁ class**. Factors not anti-isomorphic to themselves remain the control |
| **Q2R-04** | **BW orientation pricing.** "Relocated into the spectrum condition" is replaced by: the geometric sign / orientation in BW is priced by the supplied relativistic structure (Poincaré representation, wedge choice, spectrum / future-cone condition, modular sign convention). BW earns the consistency relation "vacuum modular flow = wedge boost flow" and does not reconstruct prior spacetime structure in G2 |
| **Q2R-05** | **Passivity wording.** Finite control: H = +K passive, H = −K active, so **imposing** passivity distinguishes the signs in this control. Theorem grade (Pusz–Woronowicz): complete passivity ↔ KMS / ground states. Passivity can orient a dynamics relative to a state **once the Kelvin-type condition is imposed**; this relocates orientation into the postulate. "Passivity derives / selects the arrow" is not claimed |
| **Q2R-06** | **KMS sign is convention-dependent.** In the code convention σ_t(A) = ρ^{it}Aρ^{−it}, the finite result is β = −1 for σ_t and β = +1 for σ_{−t}. This sign is a modular-parameter convention, not physical evidence for either orientation |

**Revised G2 terminal:**
- ORIENTATION NOT SELECTED;
- H_epoch NOT SELECTED;
- type III₁ conditionally supplies an injective state-independent δ_M (no representative, no physical time, no orientation);
- BW gives a geometric consistency relation under supplied relativistic hypotheses;
- thermodynamic orientation requires a supplied passivity / KMS criterion;
- TRUE COMPRESSION = 0.

## G3 — source grades and notes (no repairs to earlier gates)

| item | grade |
|---|---|
| Borchers covariance Δ^{it}U(a)Δ^{−it} = U(e^{−2πt}a); HSMI ⇒ positive-energy affine representation; P = (log Δ_N − log Δ_M)/2π | used at the scope of the **owner-supplied modern restatement** (SOURCE-TEXT VERIFIED by the owner) |
| Wiesbrock, CMP 157 (1993) 83–92 (−hsm convention Δ_M^{−it}NΔ_M^{it} ⊆ N, t ≥ 0; correspondence with chiral CFT / III₁ subfactors); CMP 158 (1993) 537–543 | SOURCE LOCATED (bibliographic record + abstract via web search) |
| Araki–Zsidó, math/0412061 | SOURCE LOCATED |
| Kähler–Wiesbrock, JMP 42 (2001) 74–86 (3+1 Poincaré from algebras in specified modular position; net via BW) | SOURCE LOCATED (abstract) |
| Buchholz–Dreyer–Florig–Summers (2000), geometric modular action; Borchers–Yngvason thermal-state modular groups | CITED FROM THE STANDARD LITERATURE / MEMORY, TEXT NOT RE-READ. Not verdict-bearing |
| Lemmas F1, F2 (finite obstruction) | proved here |
| Antiunitary maps flip the hsm sign (Θ·Δ^{is}·Θ⁻¹ = Δ_Θ^{−is}) while P stays ≥ 0 | elementary QFT-SCOUT observation |

## QFT REPAIR 03 — sharpen G3 reconstruction scope (owner, after `80f250a`; G3 accepted provisionally)

Unchanged: H_cross CONSTRAINED-NONUNIQUE; orientation NOT SELECTED; H_epoch NOT SELECTED; TRUE COMPRESSION = 0.

| ID | change |
|---|---|
| **Q3R-01** | **HSMI ⇒ Borchers triple, not automatically a non-trivial local interval net.** Keep: HSMI (N ⊂ M, Ω) ↔ one-dimensional Borchers triple (M, U, Ω) (Wiesbrock / Araki–Zsidó scope), giving positive U(a), modular / dilation covariance and the ordered family M(a) = U(a)MU(−a). Candidate intervals A(a, b) = M(a) ∩ M(b)′ can be trivial: explicit HSMIs with N′ ∩ M = ℂ·1 exist (owner pointer, arXiv:2111.03172). New Σ grade: **PARTIAL CONDITIONAL COMPRESSION OF ONE-DIMENSIONAL ORDERED LOCALIZATION STRUCTURE.** Full frozen Σ remains uncompressed |
| **Q3R-02** | **Type III₁ scope.** F1 / F2 are kept (bounded generator + one-sided invariance ⇒ two-sided; finite standard proper HSMIs obstructed). "A non-trivial HSMI forces III₁" is narrowed to: **standard non-trivial HSMIs under Wiesbrock's additional relative-commutant hypothesis (factor + cyclicity / standardness of N′ ∩ M) force type III₁.** Not claimed from half-sidedness alone |
| **Q3R-03** | **Normalization ≠ physical scale.** N = U(1)MU(−1) fixes the translation parameter relative to (N, M). It does not derive a meter, a proper distance, an absolute null interval or a physical length scale. The toy's gaps used a supplied coordinate. **PARAMETER NORMALIZATION CONDITIONALLY FIXED; PHYSICAL SCALE SUPPLIED** |
| **Q3R-04** | **Logical interdefinability, not measured information equality.** Borchers + Wiesbrock interconvert a signed HSMI and the positive half-sided Borchers-triple data under the theorem hypotheses: **LOGICALLY INTERDEFINABLE AT THE THEOREM SCOPE**, hence **RELOCATION INTO HALF-SIDED MODULAR POSITION** of positivity and half-sided translation structure. No Shannon / Kolmogorov information is claimed equal. "Cheaper specification" is retained only informally |
| **Q3R-05** | **Source upgrades.** Wiesbrock CMP 157 (1993): **PRIMARY / SOURCE-TEXT VERIFIED** (author-uploaded full text inspected by the owner) for HSMI → U(a), the positive generator, the modular scaling relation, N = U(1)MU(−1), and the qualified III₁ statement. Kähler–Wiesbrock 2001: **PRIMARY-PUBLISHER ABSTRACT VERIFIED** for exactly a 3+1 Poincaré representation from the modular groups of algebras in specified modular position, plus a 3+1 net when combined with BW. Araki–Zsidó: SOURCE LOCATED / ABSTRACT VERIFIED |
| **Q3R-06** | **Revised G3 terminal** (nine points, as recorded in `G3_MODULAR_RECONSTRUCTION_RESULT.md`) |

## G4 — no repairs to earlier gates

- **The bar** is carried unchanged from Bridge-1 B5-0. The comparator is ordinary AQFT / relativistic QFT with the same
  premises.
- **New G4 numerics are illustrations:** lattice mass / regulator / buffer variation and the 80-digit 2π slope.
- **Q2a's "CFT universality" of the c/6 coefficient** is a KNOWN-RESULT IMPORT (Calabrese–Cardy-type), cited from the
  standard literature, text not re-read. The lattice reproduces it numerically.

## QFT REPAIR 04 — final G4 scope (owner freeze ruling, after `9a96734`; G4 terminal and prediction count unchanged)

| ID | change |
|---|---|
| **Q4R-01** | **The thermal modular profile does not disprove local 2π universality.** The finite thermal-interval lattice control shows only that the chosen interval's full modular Hamiltonian is not the vacuum BW boost profile. Replaced wording: "the 2π slope fails for a thermal state", "state-priced", "vacuum-only". Corrected reading: **the exact BW wedge identification is vacuum-specific; in the finite thermal-interval lattice control the full modular profile is not boost-like; this control does not adjudicate the more general continuum question of universal local 2π behavior arbitrarily near an entangling surface** (cf. the broader literature on universal local modular temperatures, e.g. arXiv:1611.08517, owner pointer, not re-read). Q3 stays **STANDARD AQFT / QFT STRUCTURE**: it fails the bar because it is standard, BW / modular-framework-dependent and not GRUT-distinctive, **not** because of thermal behavior. The script / log line "STATE-dependent" in `g4/g4_payoff.*` is kept as emitted, and this entry takes precedence |
| **Q4R-02** | **Saturation wording.** "No identified open item within AQFT + modular theory could produce TRUE COMPRESSION" is replaced by: **"No route identified or tested under the present QFT-SCOUT-1 charter yields TRUE COMPRESSION or a distinctive observable."** Added: "A further result would require either an untested theorem / structure within AQFT or a materially altered premise envelope." No impossibility over all AQFT / modular theory is claimed |

**Status:** G4 accepted. **QFT-SCOUT-1 FROZEN** by owner ruling, with Q4R-01 and Q4R-02 applied.
