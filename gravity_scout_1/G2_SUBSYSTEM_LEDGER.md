# G2 SUBSYSTEM LEDGER — what each subsystem structure buys, and what it is priced by

> **Repaired by GRAVITY REPAIR 02 (GR2-01 … GR2-06; `GRAVITY_CORRECTION_LEDGER.md`).** G2 is accepted provisionally
> after that repair. Where wording differs, the ledger takes precedence. Script and log outputs are kept as emitted.

## 1. New priced premises opened by G2 (added to QG-1 … QG-7)

| ID | premise | status |
|---|---|---|
| **QG-8** | asymptotic structure (asymptotically flat or AdS boundary), so that total Poincaré / ADM charges are defined | supplied |
| **QG-9** | perturbative order (O(κ), κ ∝ √G) | supplied (a truncation choice) |
| **QG-10** | dressing prescription (which gravitational "Wilson line" / dressing field renders an operator gauge-invariant) | supplied; non-unique |
| **QG-11** | observable class / resolution: which boundary or exterior observables count, including exact spectral projections of the boundary Hamiltonian | supplied |
| **QG-12** | the vacuum premises used in fine-grained split failure: a unique energy-minimizing state, plus a density (Reeh–Schlieder-type) property of boundary-generated states | supplied (per the scoped source) |
| **QG-13** | an observer / clock with energy bounded below (de Sitter static patch, CLPW) | supplied |

Carried: **U_ε** (the extended neighbourhood / collar), the **state class**, and **QP-5** (the quantum lift).

## 2. Structure-by-structure ledger

| structure | buys | priced by | residual reading |
|---|---|---|---|
| **S0. Flat AQFT split** A_in ⊂ 𝒩 ⊂ A_out, 𝒩 type I (G = 0) | normal product extensions of arbitrary normal marginals; independent state specification; spatial tensor-product implementation A_in ∨ A_out′ ≅ A_in ⊗̄ A_out′ | split / nuclearity (QP-4); collar d; the choice of 𝒩 (non-unique; canonical given a standard vector — Doplicher–Longo, bibliographically verified only) | H_cross = 0 admissible with collar (QFT-SCOUT-1) |
| **S1. Perturbative gravitational splitting** (Donnelly–Giddings, O(κ)) | localization of information **modulo total Poincaré charges**. Within subspaces of fixed charge matrix elements, exterior measurements at the tested order do not resolve the interior state | QG-8 (charges), QG-9 (order), QG-10 (dressing), U_ε, state class | independence survives **only charge-sector-wise**; the charges are shared (outside-measurable) data |
| **S2. Fine-grained split failure** (Raju, scoped examples) | no independent inside / outside state specification: boundary-near observables fix the full state | QG-8, QG-11 (exact boundary observables), QG-12 | **ordinary QFT split independence forbidden in class** at the source's scope; a Doplicher–Longo-style 𝒩 is blocked / not applicable there [GR2-02] |
| **S3a. Witten crossed product** | type III₁ → type II∞ (crossed product by the modular group); a semifinite trace; entropy defined up to a state-independent constant | specific large-N emergent black-hole setting; the boundary Hamiltonian fluctuation used in the crossing | type change, **no type-I interpolation reported** |
| **S3b. CLPW de Sitter static patch** | operators dressed to an observer worldline → type II₁; a finite trace; a maximum-entropy state (empty dS) | QG-13; the de Sitter static-patch setting | type change, **no type-I interpolation reported** |

## 3. Finite control (`g2/g2_subsystem.py`, `.log`): FINITE-DIMENSIONAL ILLUSTRATION ONLY

Inside: 4 levels with "charge" H_in = diag(0, 1, 1, 2), i.e. sectors of dims 1, 2, 1. Outside: 5 levels. Everything is
type I. The toy does **not** model gravity. It switches the priced ingredients on and off.

| case | outside algebra | dim A_out | dim A_in (= commutant) | center(A_in) | product states |
|---|---|---|---|---|---|
| **T0** (G = 0) | 1 ⊗ B(K) | 25 | 16 | 1 (factor) | any marginals |
| **T1** (charge dressing: outside also holds H_in) | alg{1 ⊗ B(K), H_in ⊗ 1} | 75 | 6 = Σ_E n_E² | 3 (= # sectors) | **only for charge-sharp interior marginals.** H_in lies in both algebras. The naive ρ ⊗ σ has ω(H_in·H_in) − ω(H_in)² = Var_ρ(H_in) = 0.584 ≠ 0. A sector-E restriction is a type-I factor |
| **T2** (outside holds P₀, the vacuum projector of H_tot = H_in + H_out + λV) | alg{1 ⊗ B(K), P₀} | see below | | | |

**T2 dimensions** (dim alg / dim A_in), at relative resolution ε:

| λ | Schmidt(Ω) | ε = 10⁻² | 10⁻⁴ | 10⁻⁶ | 10⁻¹⁰ |
|---|---|---|---|---|---|
| 0 | 1, 0, 0, 0 | 50 / 10 | 50 / 10 | 50 / 10 | 50 / 10 |
| 10⁻³ | 1, 9.1e-5, 6.9e-5, 2.4e-5 | 50 / 10 | 50 / 10 | **400 / 1** | 400 / 1 |
| 10⁻² | 1, 9.1e-4, 6.9e-4, 2.4e-4 | 50 / 10 | **75 / 6** | 400 / 1 | 400 / 1 |
| 10⁻¹ | 1, 9.1e-3, 6.9e-3, 2.4e-3 | 50 / 10 | 400 / 1 | 400 / 1 | 400 / 1 |

**Readings** (illustration grade):
1. With an **outside-cyclic vacuum** (λ ≠ 0, full Schmidt rank) and the **exact** vacuum projector, b₁P₀b₂ spans all
   matrix units. A_out = B(H) and A_in = ℂ: **no interior subsystem survives**. This is the algebraic skeleton of the
   fine-grained mechanism.
2. With a **product vacuum** (λ = 0), P₀ adds only |0⟩⟨0|_in ⊗ B(K): vacuum-sector data only.
3. **Resolution ladder.** As the resolution ε is coarsened, the computed algebra steps down: B(H) (400/1) → an algebra
   with the **same dimensions as the T1 charge-sector structure** (75/6, at λ = 10⁻², ε = 10⁻⁴) → vacuum-sector data
   (50/10).
   - The threshold tracks the Schmidt scale of Ω (~λ).
   - ε is a **numerical proxy** for observational resolution. This is an analogy, not a derivation of any gravitational
     resolution limit.
   - The equality with T1's 75 / 6 is a **dimension coincidence only, not an algebraic identification** [GR2-03].
   - The ε hierarchy is **A_resolution COUPLING CANDIDATE — ILLUSTRATION GRADE**: discarding numerically small
     generators is not a physically derived gravitational resolution [GR2-03].
4. **The T1 direct sum ⊕_E B(ℋ_E) is an algebraic caricature** of shared charge labels, not the continuum
   Donnelly–Giddings algebra [GR2-01].

## 4. Information ledger — where did the "independence" go?

| item | flat (S0) | perturbative gravity (S1) | fine-grained (S2) |
|---|---|---|---|
| independent interior state specification | yes (with collar) | modulo charges | no (in the scoped examples) |
| shared inside / outside data | none required | total charges (QG-8) | the full state, via the boundary algebra (QG-11, QG-12) |
| subsystem algebra | type-I-split factor pair | charge-sector direct sum (toy: non-factor with center = charges) | the interior commutant is trivial in the exact boundary algebra (toy) |
| what decides the answer | d, 𝒩 | order + dressing + charges | observable class + resolution + vacuum premises |
| supplied information removed? | — | **no** (relocated into QG-8 – QG-10) | **no** (relocated into QG-11 – QG-12) |
