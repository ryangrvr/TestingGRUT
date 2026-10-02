# PROPOSED ERRATA TO FROZEN SCOUT-2 FILES (not applied; frozen branch untouched)

## E-01 (IR-01) — `results/S2-SigmaH_SH0_THEOREM.md`, Proposition 2

**Replace** "By CPR, when H is k-local at all, F_H is generically a finite union of LU × Sym orbits." **with:**

> Assumption (CPR-conditional): within the stated locality class (k, d, n), H has a unique local TPS class up to LU × Sym.
> CPR establish generic uniqueness conditionally — given at least one Hamiltonian in the class with a unique local TPS —
> and verify it explicitly for a restricted translation-invariant nearest-neighbour qubit class. No unconditional
> all-(k, d, n) theorem is attributed to them.

**Replace** case (a)'s "the joint principle adds no selection power beyond H-locality" **with:**

> If H's local TPS class is unique (the assumption above), the joint principle adds no selection power beyond
> H-locality. If H has several genuine dual local TPS classes, ψ-productness can distinguish among them.

## E-02 (IR-02) — same file, dichotomy (b), and `results/S2-SigmaH_RESULT.md` §ΣH-0 / §ΣH-16

**Replace** "(b) incompatible pair … there is a genuine trade-off: a Pareto family" **with:**

> (b) If F_H ∩ F_ψ = ∅, no frame achieves exact locality and exact productness together. This does **not** by itself
> imply several Pareto-incomparable optima: a single dominating compromise is not excluded by this proposition. Pareto
> nonuniqueness for incompatible pairs is a **numerical SCOUT classification** (S2-ΣH: H1, H5a, H7; S2-Σ), not a
> consequence of Prop. 2.

**In `S2-SigmaH_RESULT.md`:**
- replace "The numerics follow this dichotomy in every case" with "The numerics *classify* every tested case
  consistently with this dichotomy; the Pareto nonuniqueness in case (b) is numerical, not proved";
- in §ΣH-16, label the Pareto-family clause "(numerical, tested models)".

## E-03 (IR-03) — `results/S2-D-arrow_D0_THEOREM.md`, "Covered examples"

**Replace** "Rényi entropies;" **with:**

> continuous Rényi entropies S_α with α > 0 on finite-dimensional states (excluding the rank entropy S₀ = log rank ρ,
> which is discontinuous: diag(1−ε, ε) → diag(1, 0) jumps from log 2 to 0);

**Add** to the scope: "discontinuous functionals (rank, thresholded counts) are not covered."

## E-04 (IR-04) — same file

**(a)** In "Covered examples", **replace** "the S2-8 redundancy diagnostics built from these" **with:**

> the continuous mutual informations I(S:F) from which the S2-8 redundancy is built. The **thresholded** redundancy count
> R_δ (`I(S:F) >= (1-δ) H(S)`) is integer-valued and generally discontinuous at threshold crossings, so it is **not**
> directly covered. Record reversal for it was shown numerically (D7), not by D0.

**(b)** In Proposition 2's consequence, **replace** "The minimal case is Σ." **with:**

> Σ is one sufficient form of additional structure (used successfully in SCOUT-2). Its minimality is not proved.

## E-05 — `SCOUT_2_HANDOFF.md`

- §L: mark the Pareto-nonuniqueness bullet "(numerical)".
- §M scoped theorems: append "ΣH-0 Prop. 2 is conditional on CPR-type uniqueness of H's local class (IR-01); its
  Pareto clause is numerical (IR-02)".
- §J: "D0 covers continuous functionals (not rank entropy or thresholded counts)".
- §P: tick the theorem-review items, with IR-01 … IR-04 referenced.
