> **AUDIT REPAIR 01 (external audit; see `AUDIT_REPAIR_01.md`) — TC-1 RENAMED AND NARROWED.**
> **TC-1 — EARNED-LAYER NON-SELECTION UNDER G.** A finite non-zero G-moved quotient component cannot be selected
> by the currently earned G-invariant predicate set. Any selector must introduce G-breaking structure. The current
> GRUT record contains no earned G-breaking reference (W2-TC1 census).
> **Withdrawn:** the slogan "any principle fixing such a component is information-equivalent to supplying its
> value". A selector could break G through some *other* independently selected reference, or through anomalous/RG
> structure. Part (ii)'s bijection only shows that the transported principles gP are equally compatible with the
> earned layer; it does not show that every G-breaking principle reduces to supplying Q.
> Status: CANONICAL-CANDIDATE, pending W2-DT.
> **W2-DT update:** TC1-SURVIVES-ANOMALY. Standing clarification: in a quantum theory with a scale anomaly, read G as
> the **RG-covariant** scaling action (rescaling + coupling flow). Dimensional transmutation generates an RG-invariant
> scale, but its absolute value still needs a coupling at a reference scale (a G-breaking boundary datum).

# SCOUT-1 ZOOM-OUT 01 (after W1-C, W1-A, W1-I, W1-P)

## Results so far

| Probe | Outcome | One line |
|---|---|---|
| W1-C | NON-SELECTION THEOREM (family) | Earned predicates are invariant under G = ⟨`S_λ`, `T_s`, U⟩, so only G-invariant targets are eligible. X-01: dimensionless ≠ eligible. |
| W1-A | CROSS-LAYER SELECTION (supplied layers) | U(1) × translation × ν fix the soft momenta `2πν·ℤ` across statistics and interaction; the velocity stays free. |
| W1-I | NON-SELECTION (values) | Chain-topology identifiability is Lanczos gauge fixing; interior readout breaks it (X-02); access-depth theorem. |
| W1-P | SELECTION PRINCIPLE (premise-priced) | Complete passivity ⇒ single β: `T(ω)` collapses to one T. A dense bath needs only single-copy passivity. |

## The eight questions

**1. Are we actually reducing freedom?**
- **Yes, conditionally, and for the first time in the program.**
  - W1-A fixes a component exactly (soft momenta).
  - W1-P reduces a function (`T(ω)`) to one number.
- **No** for anything the earned layer alone controls (W1-C proves this for G-moved components).
- What has shrunk is the *target space*. Only G-invariants are eligible: the affine shape class of `μ_r`,
  exponents, central charge, soft momenta in lattice units, `T/H`, `a/c`, `h(p)`, profile shapes.

**2. Are purported selectors just new supplied assumptions?**

| Probe | Selector's inputs | Supplied? |
|---|---|---|
| W1-A | sector, translation group, filling | yes, all |
| W1-I | the frame and the boundary readout | yes, pure gauge choice |
| W1-P | bath preparation (equilibrium as a no-resource condition) | yes, but different in kind |

W1-P's premise is the **first selector that is both scale-free (G-invariant) and physically compelled**:
the Kelvin–Planck second law imposed on preparations. It is still supplied relative to the record, which
admits non-equilibrium profiles on purpose. Adopting it would be an **owner premise decision**. It is recorded
as an option and is **not** requested now; nothing in Wave 1 depends on it.

**3. Is the same obstruction recurring?**
Yes. Two obstructions recur:
- **(a) The scale.** No G-moved component is ever fixed (W1-C; W1-P's T; W1-A's velocity).
- **(b) Relations between supplied data.** Where a value is fixed, it is a fixed **relation** among supplied
  layers (soft momentum ≡ 2π × filling; β profile ≡ constant). This is structurally like K1 in SCOUT-0.

**4. Have two layers begun constraining one another?**
- **Yes. W1-A is the explicit "two layers jointly constrain Q though neither does separately" form.** The
  hostile tests show each layer is needed:
  - translation alone does not suffice — the staggered cell with integer filling is gapped;
  - U(1) alone does not suffice — pairing gaps it, and the survivors sit at the k → −k fixed points.
- W1-P is the same form: bath Hamiltonian × preparation.

**5. Has a new quotient object appeared?**
Three refinements:
- the **affine shape class** of `μ_r` (W1-C);
- **access depth** — 2k moments ⟺ k Jacobi layers (W1-I). It explains why short-time earned criteria (E-15
  reads `a₀`) can never reach the IR edge data;
- the **β-profile shape** (W1-P).

**6. Is there a prediction path worth prioritizing?**
- None distinctive yet.
- The fixed-point corollary (W1-C) says where any scale-free selection must live: at G-fixed points
  (criticality `λ₀ = 0`, T = 0/∞) or on weight-0 universals.
- **Reprioritize Wave 1** to the rigidity probes, which target exactly those:
  1. W1-R (unitarity ⇒ discrete c);
  2. W1-S (spectral dimension ⇒ γ);
  3. W1-G (Gleason/Busch ⇒ `h(p)`);
  4. W1-L (local tomography ⇒ complex field);
  5. W1-F (prediction-first `T/H`).

**7. Should a route family be killed?**
- **KILL:** "scale/units selection from earned predicates" (W1-C theorem).
- **KILL:** "topology/access as a selector of K" (W1-I, X-02). It is gauge fixing.
- **KEEP:** "symmetry/anomaly cross-layer selection" (W1-A) as the positive template.
- **KEEP:** "thermodynamic preparation principles" (W1-P) as the only scale-free physically-compelled
  template.

**8. What result would most change our view?**

> A **scale-free principle that does not encode a symmetry/sector/preparation datum** and still fixes a
> weight-0 IR datum (γ, c, z, `h(p)`). Example: unitarity + scale invariance forcing a discrete `c` with no
> sector input (W1-R).

If W1-R/W1-S/W1-G all reduce to "supplied Hilbert space / supplied symmetry / supplied dimension", the
campaign's answer to "can anything make GRUT choose?" converges to:

> **only by composing supplied layers, or at fixed points.**

That answer would itself be a theorem-grade endpoint (below).

## Theorem recognition: TC-1 (triggered by W1-C + W1-A + W1-I + W1-P)

> **TC-1 (Selector–supply equivalence for G-moved components).**
>
> - **Setting.** Let 𝒫 be the earned layer, G its deformation group (W1-C), and Q a quotient component on
>   which `ρ(G)` acts freely on an orbit O.
> - **Claim.** For every principle P compatible with 𝒫 such that `P ⇒ Q = q* ∈ O`:
>   - (i) every transported principle `gP` is equally compatible with 𝒫 and gives `ρ(g)q*`;
>   - (ii) the map `gP ↦ ρ(g)q*` is a bijection of the G-orbit of P onto O.
> - **Conclusion.** Relative to the earned layer, **choosing P is information-equivalent to supplying the value
>   of Q.** This is the "any selector … is equivalent to supplying Q" form, proved for G-moved components.
>
> **Proof:**
> - 𝒫 is G-invariant, so `gP ∧ 𝒫 = g(P ∧ 𝒫)` is consistent iff `P ∧ 𝒫` is.
> - Equivariance gives `gP ⇒ Q = ρ(g)q*`.
> - Freeness gives injectivity.

**Grade and labels:** FAMILY THEOREM (exact given W1-C's invariance census); REDISCOVERED-KNOWN in form
(Curie/Buckingham); NEW-IN-GRUT as a statement about the selector problem.

**Status: CANONICAL-CANDIDATE**, flagged for hostile testing. The weakest point is the census itself: G must be
a symmetry of *all* earned structure, and E-11/E-22 involve the supplied H and ħ.

**Open half (the real remaining question):** for **G-invariant** targets, is there a principle that is
neither a supplied layer nor a gauge choice? Wave 1's remaining probes test exactly this.

## Queue changes

- New Wave-1 order: W1-R → W1-S → W1-G → W1-L → W1-F; then ZOOM_OUT_02.
- **Wave-2 spawn candidates:**
  - **W2-FP:** "fixed-point selection" — can a scale-free *dynamical* principle (self-organized criticality,
    RG attractor) select `λ₀ = 0` or a critical exponent without a tuned input? It tests the fixed-point
    corollary directly.
  - **W2-TC1:** a hostile test of TC-1 — search the record (E-11, E-22, the gravitational inputs) for any
    earned absolute reference.
- **Owner options logged (not requested):** adopting complete passivity (or local detailed balance) as an
  earned-layer premise. It would cut the record's non-equilibrium profile class.
