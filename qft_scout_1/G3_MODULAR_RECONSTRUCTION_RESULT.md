# QFT-SCOUT-1 · G3 — MODULAR RECONSTRUCTION OF Σ + GEOMETRY (result)

**Central question (owner).** Can modular-position data reconstruct localization, translations and Lorentz / geometric
structure using **less supplied information** than the Σ / geometry being reconstructed? Or is modular position an
encoded form of the same information?

**Script:** `g3/g3_hsmi.py` (+ `.log`). It contains a finite / type-I obstruction (Lemmas F1, F2, proved, plus numerics)
and a light-ray geometric toy (supplied chiral BW / Hislop–Longo action, illustration).

**Source grades.** Borchers covariance Δ^{it}U(a)Δ^{−it} = U(e^{−2πt}a), HSMI ⇒ positive-energy affine-group
representation, and P = (log Δ_N − log Δ_M)/2π are used at the scope of the **owner-supplied modern restatement**
(SOURCE-TEXT VERIFIED by the owner).

| reference | grade |
|---|---|
| Wiesbrock, CMP 157 (1993) 83–92; CMP 158 (1993) 537–543 | SOURCE LOCATED |
| Araki–Zsidó, math/0412061 | SOURCE LOCATED |
| Kähler–Wiesbrock, JMP 42 (2001) 74–86 | SOURCE LOCATED; its abstract states a 3+1 Poincaré representation from algebras in specified modular position, plus a net via BW |
| Buchholz–Dreyer–Florig–Summers (2000), geometric modular action | CITED FROM THE STANDARD LITERATURE, TEXT NOT RE-READ |

## G3-0 Two directions (kept separate)

| direction | input | output | classification |
|---|---|---|---|
| **A: Borchers (1992)** | standard (M, Ω); one-parameter U(a) with a **positive** generator; U(a)MU(−a) ⊂ M for a ≥ 0 | Δ_M^{it}U(a)Δ_M^{−it} = U(e^{−2πt}a), J U(a) J = U(−a); and N := U(1)MU(−1) ⊂ M is a (−)HSMI | **CONSISTENCY / RELOCATION**: U and its positivity were supplied |
| **B: Wiesbrock converse (1993; Araki–Zsidó)** | N ⊂ M with common cyclic separating Ω; the **signed** half-sided condition Δ_M^{−it}NΔ_M^{it} ⊂ N for t ≥ 0 | a positive-energy representation of the affine (translation–dilation) group; U(a) with N = U(1)MU(−1); **P = (log Δ_N − log Δ_M)/2π ≥ 0** | a genuine **reconstruction**: positive energy is not separately supplied |

**Together, A and B make an equivalence:**

    (−)HSMI (N ⊂ M, Ω)   ⟺   (M, Ω) + a positive-energy half-sided translation group U, normalized by N = U(1)MU(−1)

**This equivalence is the key to the information audit (G3-1).**

## G3-1 Information audit of half-sidedness — **the smuggling located**

The supplied facts, item by item:

| # | supplied fact | what it encodes in the reconstructed structure |
|---|---|---|
| 1 | two specific algebras N, M | the **seed**: two members of the future translated family |
| 2 | the inclusion order N ⊂ M | the **ordering** along the reconstructed ray; which algebra is "further along" |
| 3 | a common Ω | the **state selector** (G3-5): the reconstruction is relative to Ω |
| 4 | Ω standard for N and M | **non-trivial infinite structure**: impossible in finite dimensions unless N = M (G3-2) |
| 5 | half-sidedness (invariance only for one half-line of the modular parameter) | **equivalent (by A ⇔ B) to the existence of a positive-energy translation with half-sided invariance.** It is the positivity + direction information in disguise |
| 6 | the sign (t ≥ 0 vs t ≤ 0) | whether N = U(+1)MU(−1) or N = U(−1)MU(+1): **the relation between the inclusion order and the positive-energy direction.** It flips under antiunitary maps (below) |

**Verdict on the owner's adversarial question.** Half-sided modular position is **not** cheaper *information* than "a
positive, oriented, half-sided translation structure along one ray." By the Borchers ⇔ Wiesbrock equivalence the two
are **logically interconvertible.**

It **is** a cheaper *specification*:
- the input is a finite datum: two algebras, one vector and one sign bit;
- the output is a continuum one-parameter group, the dilations, and the totally ordered continuum of algebras
  M(s) = U(s)MU(−s);
- the scale is fixed by the inclusion gap (toy: generator coefficient 1.00 for the gap 1, 2.50 for the gap 2.5).

> **Classification:** **PARTIAL CONDITIONAL COMPRESSION (specification-level)**. Information-level it is **RELOCATION
> into half-sidedness** (logical equivalence). The smuggling is located precisely in items 5 + 6 (and 2).

## G3-2 Finite / type-I hostile control — **obstruction proved**

**Lemma F1 (proved).**
- Suppose log Δ is bounded and σ_t(N) ⊆ N for t ∈ [0, ε).
- Then the derivation δ = i[log Δ, ·] maps N → N, by taking the derivative at 0⁺ of a curve that stays in the closed
  subspace N.
- So σ_t(N) = e^{tδ}(N) ⊆ N for **all** real t.
- **In finite dimensions (log Δ bounded), half-sided invariance implies two-sided invariance.**

**Lemma F2 (proved).**
- In finite dimensions, let Ω be separating for M and cyclic for N ⊆ M.
- Then dim N ≥ dim H ≥ dim M ≥ dim N, so **N = M.**

**Numerics** (`g3_hsmi.log`; B(ℂ⁴) with a random faithful ρ).
- **Random subalgebras U(B(ℂ²) ⊗ 1)U† are never half-sided.**
  - The σ_t-defect is present for both signs and symmetric to first order (relative asymmetry ≤ 2.2e-3 at |t| = 0.01).
  - Sample values: 3.4e-2 at ±0.01 and 3.3e-1 at ±0.1.
- **A modular-invariant subalgebra** (ρ = ρ_a ⊗ ρ_b) is invariant for **both** signs (defect ≤ 4e-15).
- **In the standard form, dim(NΩ) = 4 < 16 = dim H**, so a proper N is not standard.

**What the infinite / type-III setting buys.** An **unbounded modular generator** is necessary for a non-trivial HSMI
(F1). Moreover, a non-trivial HSMI forces type III₁ (Wiesbrock 1993, import). The finite and type-I control admits
**only** the trivial inclusion.

## G3-3 One translation line is not Σ — **PARTIAL CONDITIONAL COMPRESSION OF Σ (one-dimensional, chiral)**

| earned from one HSMI | still needed for Σ |
|---|---|
| one positive translation generator P; U(a) | transverse localization |
| the translated, totally ordered family M(s) = U(s)MU(−s) | causal complement structure |
| ordering along one coordinate | the full spacetime dimension |
| the dilation / modular relation (the affine group) | intersections producing bounded local regions |
| — | isotony / locality for a higher-dimensional net |
| — | identification of physical observables with these algebras |

One-dimensional localization is **not** promoted to full Σ.

## G3-4 Multiple HSMIs / modular intersections — **CONDITIONAL COMPRESSION with the geometry type RELOCATED into the relation pattern**

| imported theorem | modular-position assumptions (supplied) | output |
|---|---|---|
| Wiesbrock CMP 158 (1993) (conformal) | three algebras with a common Ω, pairwise in specified ± half-sided modular position | a positive-energy representation of the Möbius group PSL(2,ℝ); a chiral conformal net on S¹ |
| Kähler–Wiesbrock JMP 42 (2001) | a **finite set** of algebras in a **specified** modular position (pairwise HSMI / modular-intersection relations) | a representation of the **3+1 Poincaré group**; **combined with the Bisognano–Wichmann result**, a 3+1 local net |
| Buchholz–Dreyer–Florig–Summers (2000) | the condition of geometric modular action + modular stability, for an index set of wedge-like algebras | the Poincaré group as the symmetry group (in Minkowski space) |

**The key information test.**
- **The required relation table is isomorphic to an incidence pattern of the target geometry**, e.g. which half-lines /
  wedges are nested, which are related by ± half-sided position, and which intersect modularly.
- So the **geometry type**, meaning the dimension and causal / incidence structure, is **RELOCATED** into the relation
  pattern.
- But **a finite seed generates a continuum net and a Lie-group representation.** That is description-level
  **CONDITIONAL COMPRESSION.**
- **Caveat:** Kähler–Wiesbrock's net step itself consumes BW, i.e. the geometry → modular identification. So the full
  net reconstruction is **partly circular** with respect to G2's BW pricing.

## G3-5 State pricing — **VACUUM / STATE-SELECTOR SUPPLIED**

- The HSMI is relative to the common Ω. Changing the faithful state changes the modular groups by Connes cocycles
  u_t ∈ M, which generally lie outside N. So **the half-sided property is generically lost**.
- Where other states do carry HSMI-type structure, the reconstructed translations are relative to that state and need not
  be the physical ones. This is a subtle point in the thermal-state modular literature (e.g. Borchers–Yngvason),
  **cited from memory, not verified**, and not used for any verdict.
- Either way, **the reconstruction is not algebra-alone.**

## G3-6 Orientation firewall — **positive energy DERIVED-given-HSMI; orientation RELOCATED (G2 preserved)**

- **Positivity is convention-free.** Δ is canonical (S = JΔ^{1/2}, SAΩ = A*Ω). So P = (log Δ_N − log Δ_M)/2π is a
  canonical operator, and its positivity (for −HSMI) is a **theorem consequence**, not a separately supplied input.
- **The sign bit is relative to the complex structure.** For any antiunitary Θ, Θ·Δ^{is}·Θ⁻¹ = Δ_Θ^{−is}. So Θ maps a
  −HSMI to a +HSMI, while ΘPΘ⁻¹ stays **positive**. (This is an elementary observation of QFT-SCOUT.)
- **The toy confirms the structure.** The mirrored configurations, right rays (−hsm) and left rays (+hsm), reconstruct
  the **same** positive (rightward) translation generator. The sign records how the inclusion order sits relative to the
  positive-energy direction.
- **Reading.** Positive energy is **earned from the signed inclusion** (relocated into half-sidedness). Physical
  orientation, i.e. "future = the positive-energy direction", rests on the hsm sign + the inclusion order + the lift's
  complex-structure sign (QP-5), all supplied.
- **G2's ORIENTATION NOT SELECTED is preserved.** No thermodynamic arrow is involved.

## G3-7 Lorentz / Poincaré grades (separately)

| structure | grade |
|---|---|
| translations (one null direction) | **reconstructed** from one HSMI (CONDITIONAL; equivalent to a positive half-sided translation) |
| boosts / dilations (one) | the modular group itself acts as the dilation (CONDITIONAL; given HSMI) |
| full Lorentz / Poincaré group | **CONDITIONAL** on a finite set of algebras in a specified modular pattern (Kähler–Wiesbrock); the pattern encodes the group type |
| spacetime geometry | **RELOCATED** into the relation pattern (dimension, incidence) |
| local net | **CONDITIONAL**, and the 3+1 step consumes BW (partly circular) |

## G3-8 Σ comparison with the frozen residual

**Not asked:** "can geometry be reconstructed after a full net is supplied?" That is RELOCATION by definition.

**Asked:** can a small modular seed generate a larger localization net?
- One HSMI generates a **one-dimensional** net from two algebras + Ω + a sign. A finite modular pattern generates a
  3+1 Poincaré representation and, with BW, a net.
- **The seed encodes the geometry type** (the relation pattern) and **the positivity / direction** (half-sidedness), so
  **TRUE COMPRESSION is not earned.**

> **Terminal: PARTIAL CONDITIONAL COMPRESSION of Σ** (specification-level), with the geometry type and orientation
> RELOCATED into the modular-position pattern.

## G3-9 Dimension — **DIMENSION REMAINS SUPPLIED**

The same abstract one-HSMI data arise for light-ray (null-translated wedge) inclusions in any spacetime dimension
(Borchers' direction applies to wedge algebras and null translations generally). The 1D chiral, circle and 3+1 cases
need **different numbers and patterns** of algebras. **The dimension is encoded in the supplied relation pattern, not
selected.**

## G3-10 Payoff preview (no claim)

- The reconstructed structures are **abstract symmetry and localization** (affine, Möbius, Poincaré).
- The one fixed number is the universal modular factor **2π** in Δ^{it} ↔ boost(−2πt), the Unruh / Bisognano–Wichmann
  normalization. That is standard structure.
- No fixed spectrum or forced geometric ratio beyond group theory was found.
- G4 will apply the seven-criterion bar.

## Terminal

> **G3 TERMINAL:**
> - **One HSMI:** PARTIAL CONDITIONAL COMPRESSION (specification-level) of a one-dimensional localization + positive
>   translation + dilation structure. Information-level it is **RELOCATION into half-sidedness** (the Borchers ⇔
>   Wiesbrock equivalence).
> - **Positive energy:** DERIVED given the signed HSMI, i.e. relocated into it, not separately supplied.
> - **Orientation:** RELOCATED (hsm sign + inclusion order + complex structure); **NOT SELECTED** (G2 preserved).
> - **Σ:** PARTIAL (1D) CONDITIONAL COMPRESSION.
> - **Poincaré / 3+1 net:** CONDITIONAL on a finite modular pattern that encodes the geometry type; the net step is
>   partly BW-circular.
> - **Dimension:** SUPPLIED.
> - **Vacuum / state:** SUPPLIED.
> - **TRUE COMPRESSION = 0.**

## Scope

- **Exact controls:** the finite obstruction (F1, F2 proved; numerics exact).
- **Illustration only:** the light-ray toy uses the supplied chiral BW geometric action.
- **Imports:** Wiesbrock / Kähler–Wiesbrock / Araki–Zsidó are source-located; BDFS is literature-cited. Owner upgrade is
  welcome.
- Auxiliary to canonical GRUT (QP-5).
