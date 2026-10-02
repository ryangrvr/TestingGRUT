> **OWNER REVIEW (after `59d77f6`): ZOOM_OUT_03 accepted provisionally.** The working decomposition is
> `C5 → D_dyn ⊕ H ⊕ Σ ⊕ A_res`. T2-2 stays a CLASSIFICATION ONLY.
> - **D_dyn** = dynamical structure once a decomposition is available.
> - **H** = state / basin / measure / preparation.
> - **Σ** = factorization + local dimension + grouping + relevant scale / objective structure.
> - **A_res** = residual access / readout choices after Σ.
>
> **S2-3b is upgraded.** It is an honest case of H splitting: **H_measure → D** (the invariant measure is fixed by the
> structure of the dynamics, compact case) and **H_preparation → A** (forgetting is relative to accessible readouts;
> fine-grained information is conserved). H is no longer monolithic.
>
> **Priority change:** S2-Σ runs together with S2-G1 (local dimension), with **no weighted objective first**. A
> zoom-out follows before S2-H2.

# SCOUT-2 ZOOM-OUT 03 (after REPAIR 01, S2-3b, S2-1b, S2-8, S2-7)

**Gate check.** T2-2 is allowed only after S2-1b, S2-3b and at least one access probe. All are done: S2-1b ✓,
S2-3b ✓, and two access probes, S2-8 ✓ and S2-7 ✓.

**Working frame:** T2-1′, C5 → D ⊕ H ⊕ A. A classification, not a derivation.

## The owner's question: has any of D, H or A been eliminated?

**No. None of the three walls is eliminated.** Each moved, and the motion was not toward elimination.

| Wall | Hostile probe | What happened | Eliminated? |
|---|---|---|---|
| **H** (measure / basin / preparation) | S2-3b | **Partially re-expressed**: the reference measure → *special* D (compact homogeneous geometry; fails on non-compact X₂); the preparation class → A (only spread preparations and coarse readouts forget; fine-grained contrast stays 1.0000) | **no**. The basin / phase data of S2-6 are untouched. The measure piece moved only into rigid, non-generic D |
| **D** (structure of the primitive dynamics) | S2-1b | **Sharpened**: locality needs an objective functional, a local dimension / number of factors, and a scale. Criteria disagree across (qubits vs ququarts) and within (Clifford frame) factorizations. They agree only when H has a dominant local frame | **no**. Its operational branch (Lieb–Robinson, autonomy) is **A** |
| **A** (access / readout / reference holding) | S2-8, S2-7 | **Reinforced**. Darwinism selects the pointer *given* the split, but not the split (R = 8 for every slot; R = 8 vs 0 for the same state across frames) nor the fragments (R = 8 / 4 / 2). Consistency filters but does not select (system and environment slots alike; the final basis is free; exactly consistent global sets exist) | **no** |

## The ten questions

1. **Derived or renamed?** Renamed, and in one place (S2-3b, measure → D) genuinely relocated, but only for rigid
   dynamics.
2. **Where did the measure enter?**
   - the stability of the TPS (a perturbation measure);
   - Haar on non-compact X₂;
   - the search ensembles (sampling only).

   On compact quotients, H's *measure* is earned. The *preparation class* is not.
3. **Where did access enter?** Everywhere:
   - the fragment grouping;
   - which slot is the system;
   - the history basis and time sequence;
   - the coarse readout of forgetting;
   - the autonomy / light-cone criterion for subsystems.
4. **Is the subsystem decomposition unique?** No. It is unique only given F, d, the number of factors and a scale, and
   only when a dominant local frame exists.
5. **Is composition unique?** Unchanged since Wave 1: a FAMILY COMPOSITION DISCRIMINATOR.
6. **Did a basin sneak in?** Yes:
   - the initial product state (Darwinism, histories);
   - ψ₀ itself;
   - the a.c. preparation class.
7. **Is an observer assumed?** Yes, as the holder of a split plus fragments plus a readout basis. The observer loop
   closes in S2-8: redundancy is meant to define objectivity for observers, but it presupposes the observer's
   fragmentation.
8. **Did information decrease?** No. One relocation (H-measure → special D). Several sharpenings: D gained an
   objective functional and a scale; A gained the split and the fragments.
9. **Has the C5 boundary moved?** Yes, inward onto a single object (below).
10. **The single hostile example that would destroy the candidate below:** an A-item that does **not** presuppose a
    factorization. Examples: a readout or record notion defined from (H, ψ) alone, or a probability rule whose access
    data are TPS-free. The Dowker–Kent sets are TPS-free, but they are not records. They are the closest hostile, and
    they fail on quasi-classicality.

## Theorem recognition — candidate T2-2 (classification; NOT a derivation)

> **T2-2 (split centrality).** In every toy universe tested, every A-item presupposes a choice of factorization Σ
> (TPS + grouping + scale). Examples:
> - the coarse partition (S2-3);
> - the hidden coin vs the tester (S2-4);
> - who holds the reference (S2-5);
> - the coarse readout that forgets (S2-3b);
> - the split and fragments (S2-8);
> - the history slot and basis (S2-7).
>
> The D-item "locality" that would supply Σ cannot do so without an objective functional, a local dimension and a
> scale. Its operational version (autonomy, light cones) is itself an A-criterion (S2-1b).
>
> Hence, in the tested universes:
>
>     C5 → D_dyn ⊕ H ⊕ Σ ⊕ A_res
>
> - D_dyn = dynamical structure given Σ (unique ergodicity, interaction type, non-scrambling);
> - H = state / basin data (the measure piece movable into rigid D);
> - **Σ = the factorization with its objective and scale: the common root of D-locality and A-access**;
> - A_res = residual access given Σ (fragment grouping, readout resolution, final basis, time sequence).

**Status:** a classification over 11 probes. It is labelled **CANONICAL-CANDIDATE: NO**; it is not a theorem.

What it claims: **D and A are not independent walls. They share Σ.**

What it does not claim: that Σ is derivable. The strong success condition (primitive structure → D, H, A jointly
forced) is **not met**. T2-2 is the precise statement of the obstruction.

**Hostile, already on record:** if Σ could be fixed from (H, ψ) alone by a single, uniquely justified joint
objective, the walls would collapse into D + H. Carroll–Singh-type "Hilbert-space fundamentalism" attempts this.
Stoica-type objections ("Hilbert-space fundamentalism is impossible") argue it cannot be unique (SECONDARY,
unverified here). S2-1b Case B is a concrete instance of non-uniqueness among natural objectives.

## Next (Wave 2)

| Probe | Question | Attacks |
|---|---|---|
| **S2-Σ** | joint selection: maximize a single combined objective (locality × quasi-classicality × redundancy) over TPSs for a fixed (H, ψ). Is the maximizer unique? What selects the weights of the combination? | Σ (T2-2) |
| **S2-G** (active) | dimension notions separately: local Hilbert d (now known to be a Σ input), graph / spectral, spacetime, capacity | Σ / D |
| S2-H2 | can the S2-6 basin data (ordered phase, preparation class) move into special D the way the measure did in S2-3b? | H |
