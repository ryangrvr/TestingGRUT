> **AUDIT REPAIR 01:** the U(1) no-go is scoped. **A unitary 1+1D conformal IR theory with a non-trivial
> continuous U(1) current algebra has c ≥ 1.** So W1-A-style gapless U(1) Luttinger parents lie outside the c < 1
> FQS region.
> - This is **not** a theorem about every lattice system that merely has a U(1) conserved charge. Such a system
>   may be gapped, or its U(1) may not appear as a current algebra in the IR.
> - The mutual exclusion is between the **two declared gapless-conformal selector mechanisms** on the same parent.

# SCOUT-1 W1-R RESULT — unitarity rigidity and the IR central charge (C04)

**Charter:** `PROBE_CHARTERS.md` §W1-R. Preregistered outcomes: RIGIDITY SELECTION (premise-priced) / NOT FIXED.

**Files:** `w1r_unitarity_rigidity.py`, with log `w1r_unitarity_rigidity.log`.

**Labels:**
- RIGIDITY SELECTION (premise-priced);
- KNOWN RESULT IMPORT: Friedan–Qiu–Shenker 1984; GKO; Calabrese–Cardy entanglement; Bethe-ansatz XXZ —
  STANDARD-TEXTBOOK ✓;
- cross-layer NO-GO for GRUT's U(1) sectors;
- NON-DISTINCTIVE.

## 0. Verdict

> **RIGIDITY SELECTION — `P ⇒ Q ∈ S` with S discrete — but its premises price it out of GRUT's own sectors.**
>
> 1. **The rigidity is real and needs no sector datum.** Unitarity + Virasoro symmetry + `c < 1` gives
>    `c ∈ {1 − 6/(m(m+1))}`, with exponents in the finite Kac sets. The lattice lands on the quantized value
>    independent of microscopic couplings:
>
>    | Ising-class parent | c |
>    |---|---|
>    | integrable XY anisotropy g = 1.0 / 0.6 / 0.3 (N = 400) | 0.503 / 0.510 / 0.525 (crossover at small g) |
>    | TFIM | 0.4961 → 0.4988 |
>    | **non-integrable** TFIM + NNN (J2 = 0.3, `h_c = 1.4804` found by crossing) | 0.5052 → 0.5010 |
>
>    This is the closest any probe has come to "a scale-free principle fixes a weight-0 IR datum without
>    a sector input".
> 2. **But three supplied/tuned inputs remain:**
>    - **criticality** — a codimension-1 tuning of h. Away from `h_c` the effective c drifts (0.39 / 0.45 at
>      ±0.05) and → 0. Criticality is the W1-C G-fixed point `λ₀ = 0`: it is reached by **tuning**, never
>      selected;
>    - **`c < 1`** itself;
>    - **which member m** — that takes the symmetry class (Z₂ → Ising, m = 3), a supplied symmetry.
> 3. **Cross-layer NO-GO for GRUT's sectors.** A gapless 1D parent with a conserved U(1) charge carries a U(1)
>    current algebra, so `c ≥ 1`. That covers every SF-1/P-02/P-02b sector parent and the W1-A LSM systems.
>    - At `c = 1` the rigidity disappears: the XXZ `S^z = 1` scaling dimension varies **continuously** with
>      Δ (ED vs Bethe ansatz):
>
>      | Δ | −0.6 | 0 | +0.6 |
>      |---|---|---|---|
>      | ED, L = 20 | 0.1482 | 0.2505 | 0.3517 |
>      | Bethe | 0.1476 | 0.2500 | 0.3524 |
>
>    - **The very layer that W1-A used to fix the soft momenta (U(1)) is the layer that removes FQS rigidity.**

## 1. Results

**A. Kac table:**

| m | c | Primary weights h |
|---|---|---|
| 3 | 1/2 | 0, 1/16, 1/2 |
| 4 | 7/10 | 0, 3/80, 1/10, 7/16, 3/5, 3/2 |
| 5 | 4/5 | 0, 1/40, 1/15, 1/8, … |
| 6 | 6/7 | |
| 7 | 25/28 | |
| 8 | 11/12 | |

The values accumulate only at c = 1.

**B. Free-fermion XY chains (OBC, N = 400):**

| Chain | c (fit) |
|---|---|
| Ising line | 0.503 / 0.510 / 0.525 (g = 1, 0.6, 0.3; small-g crossover from the nearby XX line) |
| XX line | 0.992 / 0.992 (h = 0, 0.5) |
| off-critical | 0.000 |

**C. Non-integrable TFIM + NNN:**
- `h_c(J2 = 0) = 1.00030` (exact 1), which validates the crossing method.
- `h_c(J2 = 0.3) = 1.48044`.
- Pairwise c from the PBC `S(L/2)` for L = 8…18:
  - J2 = 0: 0.4961, 0.4974, 0.4980, 0.4985, 0.4988;
  - J2 = 0.3: 0.5052, 0.5032, 0.5021, 0.5014, 0.5010.
  - Both converge to 1/2 from opposite sides.

**D. XXZ scaling dimension `x₁` for L = 12/16/20 vs Bethe:**

| Δ | −0.6 | −0.3 | 0 | 0.3 | 0.6 | 0.9 |
|---|---|---|---|---|---|---|
| L = 20 | 0.1482 | 0.2021 | 0.2505 | 0.2988 | 0.3517 | 0.4159 |
| Bethe | 0.1476 | 0.2015 | 0.2500 | 0.2985 | 0.3524 | 0.4282 |

At Δ = 0.9 the convergence is slow (the marginal operator near the KT point Δ = 1); this is known.

## 2. Ten-point hostile test of "unitarity + Virasoro + c < 1 ⇒ c ∈ Kac"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | No: c is fixed to 1/2 across integrable/non-integrable couplings. |
| 2 | P silently contains Q? | No for the set. **Yes for the member** (symmetry class). |
| 3 | Representation-dependent? | No. c is a universal IR datum. |
| 4 | Physical or gauge? | Physical (entanglement, finite-size energy, specific heat). |
| 5 | Standard? | Yes (FQS 1984). |
| 6 | Parent variation? | Passes within the Ising class; the XX/XXZ hostile behaves as predicted. |
| 7 | Composition? | c is additive. A composite of rigid sectors is rigid **only while the total c < 1**. Two Ising copies (c = 1) leave the rigid region (the Ashkin–Teller line has continuous exponents) — standard. |
| 8 | Coarse-graining? | c is the RG-fixed-point datum: it survives coarse-graining by definition (c-theorem monotone). |
| 9 | Unique or stationary? | The set is discrete; the member is not unique without the symmetry. |
| 10 | Boundary condition selecting? | **The criticality tuning does the selecting of "being at a fixed point".** It is a codimension-1 choice, the G-fixed point of W1-C. |

## 3. What this means

- **Answer to ZOOM_OUT_01 Q8.**
  - **Partly yes.** A scale-free principle with **no sector datum** fixes a weight-0 IR datum to a
    discrete set.
  - **But only at a tuned fixed point and only below c = 1.**
  - **GRUT's U(1)-sector parents sit at c ≥ 1, where the principle is empty** — and continuous exponents
    (Luttinger K) are exactly what P-02b found to be supplied.
- **TP-1 refined.** The rigidity selectors that do not need a supplied sector live at fixed points (the W1-C
  corollary). Reaching the fixed point is a tuning. Choosing the member is a symmetry datum. The pattern
  becomes:

  > **selection = (fixed point, tuned) × (symmetry, supplied) × (rigidity theorem, imported).**

- **Cross-layer NO-GO (new placement):**
  - U(1) sector (needed for W1-A's selection) ⇒ c ≥ 1 ⇒ no unitarity rigidity ⇒ IR exponents free.
  - The two Wave-1 selectors are **mutually exclusive** on the same 1D parent: soft momenta can be fixed
    (W1-A), or the exponents can be quantized (W1-R), but not both by these theorems.
  - Exception: discrete-symmetry-broken U(1) at c ≥ 1 with additional structure (e.g. SU(2)₁ at the
    isotropic point fixes K = 1/2). That is again a supplied symmetry.
- **Empirical:** universal c and exponents are textbook. NON-DISTINCTIVE.

**Status: W1-R COMPLETE — RIGIDITY SELECTION (premise-priced: criticality tuned, c < 1, symmetry class).
Cross-layer NO-GO: U(1) sectors (c ≥ 1) are outside the rigid region; W1-A and W1-R selectors are
mutually exclusive on one 1D parent.**
