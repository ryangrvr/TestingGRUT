> **REPAIR 01 (owner audit):** the headline is scoped to a **FAMILY COMPOSITION DISCRIMINATOR**. In the tested finite
> GPT family, continuous reversible entangling interaction excludes the classical/simplex composite and the gbit
> min/max polytopic composites, while the real and complex quantum tensor composites pass. This is **not** a general
> theorem that continuous interaction forces tensor composition. Other non-polytopic GPT composites are untested.

# S2-2 RESULT — which composition rule is forced?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-2 (pre-registered). Local tomography is **not** used as an input.

**Files:** `probes/S2-2/s2_2_composition.py`, with log `probes/S2-2/s2_2_composition.log`.

**Labels:**
- **COMPOSITION PARTLY FORCED BY A DYNAMICAL PRIMITIVE** (continuous reversible interaction) + **NONUNIQUE at the
  field level**;
- KNOWN RESULT IMPORT: Gross–Müller–Colbeck–Dahlsten 2010 (boxworld has no non-trivial reversible dynamics); finite
  automorphism groups of polytopes — standard; core facts checked ✓.

## 0. Verdict

> **Independent preparability selects nothing.** Product states exist in every composite tested: Cartesian / simplex,
> gbit min tensor, gbit max tensor, real tensor, complex tensor.
>
> **A dynamical primitive, "continuous reversible interaction between the parts", separates the classes:**
>
> | Composite | Continuous interaction? | Reason |
> |---|---|---|
> | classical bits (simplex) | no | reversible maps = S₄ (finite; CNOT is a discrete interaction); Lie algebra 0 |
> | gbits, max tensor (boxworld / NS polytope) | no | finite automorphism group. No automorphism maps a product vertex (12 incident facets) to a PR vertex (8): an incidence invariant ✓ |
> | gbits, min tensor (local polytope) | no | finite group; every vertex is a product state |
> | **real rebits** (real tensor) | **yes** | `exp(t·iY⊗X)` is real orthogonal and generates entanglement 0 → 0.389 → 0.717 → 1 (t = 0, 0.2, 0.4, π/4) ✓ |
> | **complex qubits** (complex tensor) | **yes** | `exp(−itZ⊗Z)` gives the identical curve ✓ |
>
> So the criterion selects **non-polytopic (strictly-convex-type) systems composing by a linear tensor product with
> entangling continuous dynamics**. It does **not** fix the scalar field: ℝ passes as well as ℂ.

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| the single-system state-space shape (simplex, square, disc, ball) | **inserted** |
| "dynamics is reversible and continuous in time" | **inserted** (the primitive premise) |
| "composite dynamics can couple the parts" | **inserted** (interaction exists) |
| exclusion of polytope theories (classical, boxworld) | **forced** given the above |
| tensor-type composite vs min / max tensor | **forced** given continuous interaction, among the tested family |
| the field (ℝ vs ℂ) | **not forced**: still needs LT or something equivalent (SCOUT-1 D3; S2-5 target) |

## 2. Accounting

- **The C5-B choice is converted** into the dynamical premise "continuous reversible interacting dynamics" plus
  the single-system shape. This is essentially Hardy's / Masanes–Müller's continuous reversibility, recast as a
  statement about **dynamics** (a primitive-level property) rather than as an operational axiom.
- **That recasting is the useful output.** A **C5-B selector exists at the level of dynamics**, but it leaves the
  field and the dimension rule open.
- **GRUT note (deferred bridge, recorded only):** SCOUT-1 found the GRUT floor generator dissipative (SUPPLIED). A
  dissipative semigroup does not supply a continuous *reversible* group. So this selector is currently unavailable
  to GRUT.

## 3. Hostile notes

- Only polytopes, discs and balls were compared. Non-polytopic, non-ball state spaces (e.g. other strictly convex
  bodies) are not tested.
- The general result that only balls admit transitive continuous reversible dynamics with interaction (Masanes et
  al.) is SECONDARY here.
- A first run reported zero rebit entanglement because the probe state |++⟩ is mapped to a product state by
  `iY⊗X`. That was a test-state error, corrected (logged in the script and as Y-01).

**Status: S2-2 COMPLETE — composition is selected by a dynamical primitive (continuous reversible interaction) over
classical and boxworld composites; the field (ℝ vs ℂ) is not selected.**
