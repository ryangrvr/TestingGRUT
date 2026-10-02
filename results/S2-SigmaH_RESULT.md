# S2-ΣH RESULT — can one principle jointly select the factorization and its independence boundary condition?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-ΣH. Pre-registered at `bab78b5` together with the ΣH-0 theorem
(`results/S2-SigmaH_SH0_THEOREM.md`).

**Files:** `probes/S2-SigmaH/s2_sigmah.py` (+ `.log`). It reuses the S2-Σ machinery: 6 qubit slots, 202 groupings, and 7
frames including the ψ-product frame W_ψ. That gives 1414 candidates per case.

**Criteria:**
- L_H = −k_eff (H-locality);
- C_ψ = Σ_f S_f − S_global (total correlation).

## Headline

> **No joint principle on (H, ψ) selected [Σ, boundary] without a preference.** The cases split exactly along the
> pre-registered ΣH-0 dichotomy.
>
> **(a) Compatible pairs** (H local and ψ product in the same frame: PC, H3, H4, H6; Gibbs approximately).
> - The joint Pareto problem has a **dominant set**.
> - After quotienting by local-unitary equivalence it lies in **one frame**, the identity frame.
> - That set survives lexicographic reorderings A / B / D and every tolerance sweep: it is contained in all 8 winner sets.
> - **The frame is selected. But:**
>   - (i) the frame was guaranteed by construction, because the pair's compatibility *is* the boundary datum H_corr|Σ.
>     That datum is detected, not derived (ΣH-0 Prop. 2).
>   - (ii) the **grouping / local dimension is not selected**. The dominant set is 5 LU classes of contiguous bipartitions
>     (32×2, 16×4, 8×8). Rule C and MDL instead pick the finest (qubit) grouping. → **PREFERENCE-ORDER PRICED at the
>     d level.**
>
> **(b) Incompatible pairs** (the generic case: H1 conflict, eigenstate, Haar).
> - **PARETO NONUNIQUENESS.**
> - Different priority orders pick different frames. The tolerance sweeps share **no** common winner.
> - C-first rules pick the **trivial** ψ-product frame W_ψ, in which H is non-local (outcome E, as ΣH-0 Prop. 1 predicts).
>
> **Overall: outcome C/D, never A.**
> - C (boundary-compatible ⇒ frame selected, d not) for non-generic compatible pairs;
> - D (Pareto / preference nonuniqueness) for generic pairs;
> - E (trivial) whenever correlation is ranked first.

## ΣH-0 (theorem; see the separate file)

**Prop. 1.** Every pure ψ is a product state in a coset of frames of dimension (N−1)² + 1. **State independence alone
selects no TPS.**

**Prop. 2.** Joint compatibility F_H ∩ F_ψ ≠ ∅ is non-generic. When it holds, the joint selector reproduces CPR's class
(given k, d, n) and adds nothing. When it fails, there is a Pareto trade-off. **The numerics follow this dichotomy in
every case** (table below).

## ΣH-1 / 2 / 5 / 11 — Pareto fronts (L, −C), four epoch prescriptions, LU quotient

| case | pair type | C₀ / C_min front (bare → LU classes) | dominant? | frames on C_min front (LU) | local dims on front | C_avg / C_typ |
|---|---|---|---|---|---|---|
| **PC** MFI chain + product ψ (same frame) | compatible | 10 → **5** | **yes** | **id only** | 32×2, 16×4, 8×8 | **fail**: the front excludes the id frame (W_ψ, cliff1, cliff3, commutant) |
| **H1** H local in id, ψ product in cliff1 | **incompatible** | 8 → 3 | **NONE** | id **and** cliff1 | 32×2, 16×4 | no dominance |
| **H3** TFIM ring + translation-invariant product ψ | compatible | 25 → 15 (13 after translation + reflection) | yes | id only | 32×2, 16×4, 8×8 | wrong frame (cliff1) |
| **H4** record-forming star | compatible | 2 → **1** | yes | id: (01245)(3) | 32×2 | no dominance |
| **H5a** eigenstate | incompatible | 9 → 6 | NONE | W_ψ, cliff1, commutant, id | 32×2 | same (stationary) |
| **H5b** Gibbs β = 1 | ≈ compatible | 6 → 4 | yes | id, commutant | 32×2 | same (stationary) |
| **H6** Janus mean-field state | compatible | 10 → 5 | yes | id only | 32×2, 16×4, 8×8 | unique id (here it works) |
| **H7** Haar | **incompatible** | 5 → 4 | NONE | W_ψ, cliff1, commutant, id | 32×2 | no dominance |

**Epoch (ΣH-5).**
- C₀ (a supplied epoch) and C_min (infimum over the orbit, which *picks* an epoch) agree whenever ψ was given at its
  special moment.
- **C_avg / C_typ** (a time measure / a typicality measure) **fail the positive control**: the PC and ring fronts exclude
  the native frame.
- → **The boundary condition must be a special moment, not a time average: EPOCH-PRESCRIPTION PRICED.**

**Local dimension (ΣH-2).** On every compatible front, 2-factor cuts of different d (32×2 / 16×4 / 8×8) are exactly tied.
The joint core does **not** select d. → **d REMAINS Σ-PRICED** (G1 does not collapse).

## ΣH-3 / 4 — priority orders and tolerance sweeps

| case | A (max L; min C) | B (min C; max L) | C (min k; min C; min MDL) | D (min C; min k; max P) | ε / L_min sweeps (8 tolerances, full tie sets) |
|---|---|---|---|---|---|
| PC | 10 id-frame 2-cuts | same | **id qubits (unique)** | same as A | 5 sets; **all contain the 10-member dominant set** |
| H1 | id 2-cuts | **cliff1** | **cliff1** 8×4×2 | **cliff1** | 6 sets; **common: 0** |
| H3 | 25 id 2-cuts | same | id qubits | same | 5 sets; all contain the 25 |
| H4 | id (01245)(3) | same | id qubits | same | 6 sets; all contain the 2 |
| H5a | id 2-cuts | **W_ψ** (trivial) | 31 W_ψ groupings | W_ψ | 6 sets; common: 0 |
| H5b | id 2-cuts | same | same | same | 2 sets; all contain the 6 |
| H6 | id 2-cuts | same | id qubits | same | 5 sets; all contain the 10 |
| H7 | id 2-cuts (C₀) / **commutant** (C_min) | **W_ψ** | 31 W_ψ groupings | W_ψ | 6 – 7 sets; common: 0 |

**Readings:**
- Compatible pairs: the **frame** is robust to every priority order and tolerance. The **grouping** splits between
  "L first" (2-cuts) and "k first / MDL" (qubits).
- Incompatible pairs: **PREFERENCE-ORDER and THRESHOLD PRICED**. Ranking correlation first returns the trivial W_ψ frame:
  **DEFINITIONAL / TRIVIAL (E)**.

## ΣH-6 — time translation (s = 2.25)

| case | C₀ | C_min |
|---|---|---|
| PC, H6 | the dominant set is **lost** → **EPOCH-PRICED** | the same Σ class, with the selected epoch moving from 0.00 to −2.25 → **TPS COVARIANT / BOUNDARY EPOCH RELATIVE** (the partial positive) |
| H1 | no dominance either way | no dominance either way |

(The dominant-set listings in the log show LU-equivalent Clifford copies first; ΣH-11b identifies them with the id frame.)

## ΣH-7 — time reversal (Θ = K)

ψ and Θψ give the **same Σ** and the same epoch t* = 0. C(t) is symmetric about the special epoch:
- PC: 3.643, 3.691, **0.000**, 3.691, 3.643;
- H6: 3.056, 3.029, **0**, 3.029, 3.056.

→ the joint selector finds a **low-correlation middle condition** but **no orientation** (consistent with D-arrow).

## ΣH-8 — arrow quality as a third criterion

Adding R (forward relaxation from the C_min epoch) gives:
- dominance only for PC / H6 at horizon h = 2, and only for H5b, where it was already present without R; in each case
  the dominant set has 5 ties (2-cuts);
- **no dominance at h = 5**;
- fronts that **grow** in the incompatible cases (H1: 6 → 15; H7: 11 → 16).

→ The arrow criterion **does not break ties uniquely**. It adds a horizon. → **SCALE-PRICED.**

## ΣH-9 — records under one fragment rule

A common rule (each factor is the system; each other single factor is a fragment) **exists**, but the counts are not
comparable across factor numbers. Mean record count for H4 rises with the number of factors: 1.00 (2) → 1.69 → 1.98 →
2.11 → 2.29 (6).

→ **A_res / Σ COUPLING**: the record criterion's scale is set by the factorization it is meant to judge.

## ΣH-10 — MDL

Over 45 language / precision / cost settings:
- languages: computational, local-Hadamard, product-aware;
- δ ∈ {10⁻³, 10⁻⁶, 10⁻⁹};
- c_ψ ∈ {¼, 1, 4, 64, 1024}.

| case | distinct MDL winners |
|---|---|
| PC | 3: id qubits; cliff3 4×2⁴; **W_ψ** |
| H1 | 6 |
| H6 | 4 |
| H7 | 3 |

At low c_ψ, MDL is dominated by the H-description and picks id qubits. At high c_ψ it flips to W_ψ.
→ **DESCRIPTION-LANGUAGE / COST PRICED.** Combining H and ψ into one number does not make MDL fundamental.

## ΣH-15 — outcome

| outcome | holds? |
|---|---|
| A. JOINT BOUNDARY SELECTED | **NO** |
| B. Σ selected, boundary not | **no**: Σ is selected only when the boundary is already compatible |
| C. boundary selected given Σ | **yes, inverted form**: *Σ (frame) selected given a compatible boundary*; d / grouping still preference-priced |
| D. Pareto / preference nonuniqueness | **YES** for generic (incompatible) pairs |
| E. definitional / trivial | **YES** whenever correlation is ranked first (W_ψ) |

## ΣH-16 — scoped theorem (as supported)

> State independence cannot identify subsystem structure: every pure state is a product state in a (N−1)² + 1-dimensional
> family of frames (Prop. 1). Hamiltonian locality constrains subsystem structure only relative to a locality class
> (CPR). In the tested finite systems (n = 6, 1414 candidates per case, 8 cases):
> - jointly demanding locality and independence yields a **dominant frame only for non-generic compatible pairs**, where
>   the compatibility itself is the supplied boundary datum;
> - otherwise it yields a **Pareto family**.
>
> Selecting among the family, and selecting the grouping / local dimension even in the compatible case, requires a
> **preference order, a threshold, a description language, a horizon, or an epoch prescription**. Time-averaged or
> typical-time independence does not recover the native frame. The selected boundary is covariant under time translation
> only for the infimum prescription, and it is orientation-neutral.

**Status: S2-ΣH COMPLETE — outcome C/D(+E); A not obtained; d not selected; epoch moves but does not disappear;
orientation not selected.**
