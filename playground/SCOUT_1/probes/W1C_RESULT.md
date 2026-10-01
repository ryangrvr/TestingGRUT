# SCOUT-1 W1-C RESULT — deformation-group non-selection (C06 + C07)

**Charter:** `PROBE_CHARTERS.md` §W1-C. Preregistered outcomes: NON-SELECTION THEOREM / COUNTEREXAMPLE FOUND /
UNFORMULABLE. The E-list is the record's (`SCOUT_0/BASELINE_MAP.md` Table 1, EDA-01 abridgements). No new
predicates were invented.

**Files:**
- script: `w1c_scaling_nonselection.py`
- log: `w1c_scaling_nonselection.log`

**Labels:**
- NON-SELECTION THEOREM (FAMILY THEOREM grade; row-checked over the E-list);
- REDISCOVERED-KNOWN — its logic is Curie's principle / Buckingham-π / RG fixed-point reasoning;
- NEW-IN-GRUT as a placement: it gives the *mechanism* for SCOUT-0's "classifies but does not fix".

## 0. Verdict

> **NON-SELECTION THEOREM (scoped).**
>
> Every E-entry is invariant under the deformation group G. In a weaker sense, its dimensionful identities are
> weighted-homogeneous. Therefore:
>
> 1. **The earned layer can only impose G-invariant constraint sets.** For a G-moved component these are
>    sign/class statements: `{0}`, `(0,∞)`, `{∞}` and their unions.
> 2. **The earned layer can fix values exactly only for G-invariant combinations.**
> 3. **Any principle that fixes a G-moved component at a finite non-zero value must break G.** That means it
>    supplies a reference datum: a scale, a reference pin, or a frame.
>
> The a-priori reconnaissance hypothesis ("only dimensionless components are eligible") is **necessary
> but not sufficient**. The dimensionless ratio `λ₀/W` is moved by the pin shift and is still not
> earned-fixable. This is recorded as **X-01**.

## 1. The deformation group (declared)

| Generator | Action | Admissible range | Physical reading |
|---|---|---|---|
| `S_λ` | `K → λK`, `t → t/λ`; every supplied rate (noise T, drift β·a², H, …) rescaled jointly by its weight | `λ > 0` | choice of time unit |
| `T_s` | `K → K + sI` (`μ_r` translated by s) | `s > −λ₀` (stay gapped) | pin shift inside the gapped class |
| `U` | `K → UKUᵀ`, `U e_r = e_r`, declared local frame transported with U | all such U | relabeling of the hidden sector |

## 2. Numerical checks (`w1c_scaling_nonselection.log`)

**A. Predicate invariance.** The pinned chain (n = 400, the S5-1/C1-a parent) and the pinned half-plane
(30×30, the EDA W-A witness) were each tested under `S_λ` (λ = 3, 0.2) and `T_s` (s = +0.5, −0.2).
- Every computable earned predicate is unchanged: symmetric/passive, gap > 0, real spectrum, local in the
  declared frame, CM kernel.
- The E-4-type asymmetric ring keeps its complex modes. Even the ratio `max|Im|/max|Re−λ₀| = 0.2` is unchanged.

**B. Hidden relabeling.** A random orthogonal U fixing `e_r` on the n = 120 chain:
- Moments 1…8 agree to `2·10⁻¹⁵` (relative), and the spectrum and retained weights are equal.
- The Lanczos coefficients from `e_r` agree to `2·10⁻¹⁵`.
- Bandwidth in the original frame goes from **1 to 119**. In the transported frame it is back to 1.

So **`μ_r` is the complete U-invariant of `(K, e_r)`**. This is the spectral theorem for a cyclic vector plus
Jacobi/Lanczos uniqueness (STANDARD-TEXTBOOK ✓). Locality (E-2) is U-invariant **only with the frame
transported**: the local frame is itself a supplied datum (it re-enters in W1-I).

**C. Component equivariance (chain | half-plane):**

| Component | `S_λ` | `T_s` | G-invariant? |
|---|---|---|---|
| `λ₀` (edge location) | ×λ | +s | no |
| W (support width) | ×λ | fixed | no |
| `m₁ = K_rr`, `m₂` | ×λ, ×λ² | shifted | no |
| `κ = 1/(2√m₁)` (E-14) | ×λ^{−1/2} | moves (0.3297 → 0.2988 / 0.3450) | no |
| `κ²m₁` | **0.25 fixed** | **0.25 fixed** | **yes** |
| `λ₀/W` (dimensionless) | fixed | **moves** (0.0750 → 0.2000 / 0.0250) | **no** (X-01) |
| edge exponent γ (fit) | fixed (0.4621 \| 0.3798) | fixed | yes |
| standardized cumulants (skew, excess kurtosis) | fixed (0, −1 \| 0, −0.778) | fixed | yes |

The G-invariant content of Q3 is the **affine shape class** of `μ_r`: `μ_r` modulo `λ ↦ aλ + b`. This class
carries γ, the standardized cumulants, and the visible sign. G-invariant does not mean earned-fixed. The chain
and the half-plane differ in γ and in kurtosis while satisfying the same predicates (EDA W-A). The class is
only the *eligible* target.

**D. Weighted homogeneity of the record's dimensionful identities (sympy).**
- E-15's `Δc₃ = 24βT₁a(44βa² + 5K₁₁)` is homogeneous iff `w_β + 2w_a = 1`, a consistent weight assignment.
  The record's pure numbers 12, 24, 44 and 5 have weight 0.
- E-14's `κ(λK₁₁)/κ(K₁₁) = λ^{−1/2}`. The relation `κ²K₁₁ = 1/4` has weight 0.
- These are **fixed relational coefficients**, not fixed quotient values. That is exactly what part 2 of the
  theorem permits.

**E. Fixed-point lemma.** For a weight-w quantity (w ≠ 0), `λ^w q = q` for generic λ only at `q ∈ {0, ∞}`.
`λ₀ + s = λ₀` has no solution for s ≠ 0.

## 3. Row-by-row invariance of the E-list

Class key:
- **INV** — the predicate is G-invariant;
- **EQV** — a dimensionful identity, weighted-homogeneous;
- **INV\*** — invariant only with the declared frame/sector transported;
- **N/A** — not a predicate on `(K, r)` (meta or non-implication).

| E | Content (abridged) | `S_λ` | `T_s` | `U` | Fixed values it contains | Weight |
|---|---|---|---|---|---|---|
| E-1 | gap ⇒ memory; passivity ⇒ positivity | INV | INV (open class) | INV | none (implications) | — |
| E-2 | locality ⇒ geometry | INV | INV | **INV\*** | none | — |
| E-3 | linearity ⇒ exact reduction | INV | INV | INV | none | — |
| E-4 | cycle affinity ⇒ ¬CM | INV (checked) | INV (checked) | INV\* | none | — |
| E-5 | seven non-implications | N/A | N/A | N/A | none | — |
| E-6 | no organizing pair | N/A | N/A | N/A | none | — |
| E-7 | identity readout: singleton classes | INV | INV | INV (r fixed by U) | none | — |
| E-8 | P-1/X2 intrinsic partition selection | INV | INV | INV\* | selects a **member** (which site), not a value | 0 |
| E-9 | influence cone exists | INV | INV | INV\* | cone speed `v ∝ λ` (READ) | 1 |
| E-10 | vocabulary demotion | N/A | N/A | N/A | none | — |
| E-11 | exact FRW kernel transport | EQV (H, T rescale jointly) | — | — | `T/H = 1/2π` is weight-0 (Q7; supplied background) | 0 |
| E-12 | counting asymptote `S_∞ ≈ 2d + 0.214` | INV (window in lattice units) | INV | INV\* | dimensionless number; lattice-dependent (READ) | 0 |
| E-13 | SF-1 sector ⇒ law class (z = 1 vs 2) | INV | INV | INV\* (sector transported) | z is weight-0; edge velocity weight 1 (READ) | 0 / 1 |
| E-14 | S5-1 Markov limit, `κ = 1/(2√2.3)`, `t^{−3/2}` | EQV | κ moves (W-C) | INV | 3/2 and `κ²K₁₁ = 1/4` are weight-0; κ is weight −½ (READ) | 0 / −½ |
| E-15 | S2-1 discriminator identities | EQV (checked) | — | INV | coefficients weight-0 | 0 |
| E-16 | realization-only | N/A | N/A | N/A | none | — |
| E-17 | formulable only with change | N/A | N/A | N/A | none | — |
| E-18 | no reversal | N/A | N/A | N/A | none | — |
| E-19 | net arrow; `D, Ḋ ~ t^{−6}` | EQV | INV (sign statements) | INV | exponent 6 is weight-0; prefactors move | 0 |
| E-20 | architecture | N/A | N/A | N/A | none | — |
| E-21 | sector selection (restated: counting does not select the exponent) | INV | INV | INV\* | none | — |
| E-22 | carrier elimination in a supplied quantum class | INV | INV | INV | none (ħ supplied with the class) | — |

**Tally:**
- G-moved values fixed by an E-entry: **0**.
- Every fixed number on the list (E-11 `1/2π`, E-12 `0.214`, E-13 z, E-14 `3/2` and `1/4`, E-15 coefficients,
  E-19 exponent 6) is **weight-0**.
- Every weight-0 number traces to a supplied class (background, lattice, sector, parent), per EDA-01.

`T_s` deserves a note. E-2's "pin-free fork (Outcome A)" is **held outside the floor** by the record, so the
gapless boundary is not an earned predicate. `T_s` is a symmetry of the floor as recorded.

## 4. Theorem (stated with its scope)

> **W1-C non-selection theorem.**
>
> **Setting:**
> - 𝓜 is the declared model space: linear, local, passive parents with a retained site, plus their supplied
>   rates (the C1-a/L0-1 family, S5-1, the SF-1 free parent, the S2 linear-response layer, the S6 parent).
> - G = ⟨`S_λ`, `T_s`, U⟩ acts on 𝓜 as in §1.
> - 𝒫 is any conjunction of E-list predicates.
> - Q : 𝓜 → V is any quotient component equivariant under a representation ρ of G.
>
> **Claim:**
> - (i) If 𝒫 ⇒ Q = q*, then `ρ(g)q* = q*` for every g that keeps the model admissible.
> - (ii) If 𝒫 ⇒ Q ∈ S, then S is a union of ρ(G)-orbits.
> - (iii) Consequently, no conjunction of earned predicates fixes:
>   - `λ₀`, W, any moment, κ, an edge or cone velocity, T, or the forcing amplitude (`S_λ` weight ≠ 0);
>   - `λ₀/W` (moved by `T_s`);
>   - any `K`-component beyond `μ_r` (moved by U).
>
>   The eligible targets are exactly the G-invariants: the affine shape class of `μ_r` (edge exponent γ,
>   standardized cumulants, visible sign), dynamical exponents, counting exponents, soft momenta in lattice
>   units, `T/H`, `a/c`, and relational coefficients.
>
> **Corollary (price of selection):** a principle that fixes a G-moved component at a finite non-zero value is
> not G-invariant. It must contain a G-non-invariant datum (a reference scale, pin or frame), and relative to
> the earned layer that datum is supplied.
>
> **Corollary (fixed points):** a G-invariant principle can select a G-moved component only at a G-fixed value.
> For a dimensionful quantity that means 0 or ∞ — gaplessness, criticality, zero velocity, infinite range.
> This is why scale-free (power-law, critical) data are the natural targets of any selection: RG fixed-point
> logic.

**Proof:** (i) If m is admissible, so is gm. Then `q* = Q(gm) = ρ(g)Q(m) = ρ(g)q*`. (ii) and (iii) follow from
(i) and the row checks of §3. The content is the invariance census, not the lemma.

**Grade:** FAMILY THEOREM (lemma exact; invariance row-checked on the recorded E-list, numerically on four
parents). It is not a statement about earned predicates the record does not contain.

## 5. Hostile test (ten points, applied to the non-selection claim)

| # | Question | Answer |
|---|---|---|
| 1 | Can Q vary while P holds? | Yes, along every G-orbit (exhibited numerically). |
| 2 | Did P silently contain Q? | No. The point is that P contains no scale. |
| 3 | Representation-dependent? | U-moved components are representation-dependent by construction. `μ_r` is the invariant. |
| 4 | Physical or gauge? | `S_λ` is a unit choice (gauge) **unless** an external clock/reference is supplied. Then "λ₀ in units of the reference" is weight-0 and fixing it needs the reference. |
| 5 | Already standard? | Yes: Curie's principle, Buckingham-π, RG fixed points. REDISCOVERED-KNOWN. |
| 6 | Parent variation? | Chain, half-plane, asymmetric ring, Lanczos chain: all pass. |
| 7 | Composition? | A conjunction of invariant predicates is invariant, so composing E-entries cannot escape. |
| 8 | Coarse-graining? | Decimation/RG maps commute with `S_λ` up to the flow. Fixed points of the flow are where selection can live (fixed-point corollary). |
| 9 | Unique or stationary? | Not applicable to a non-selection result. |
| 10 | Boundary doing the selecting? | Yes, and that is the loophole: the gapless boundary `λ₀ = 0` is `S_λ`-fixed. It is held **outside** the floor (E-2 Outcome A), so it is not earned. A future principle demanding criticality would be a supplied boundary choice. |

**Weakest point:** the theorem is only as strong as the claim that G is a symmetry of **all** earned structure.
- An earned result referring to an absolute reference would break it.
- Candidates are the cosmological background (E-11: H supplied, S-13) and ħ (supplied with the quantum class,
  A-6).
- Both are on the **supplied** list, so they do not count as counterexamples.
- If the owner ever admits H₀ or ħ as earned, part (iii)'s first bullet weakens to "fixed only in units of
  H₀/ħ".

## 6. What this does to SCOUT-1

- **Narrows the hunt.** Wave-1 probes are now judged only on G-invariant targets. W1-A (soft momenta in lattice
  units), W1-R (central charge), W1-S (edge exponent), W1-G (`h(p)`), W1-L (number field) and W1-F (`T/H`)
  are all weight-0, so all remain eligible. W1-P targets the *profile* `T(ω)` → constant, a weight-0 shape
  statement; T itself stays unfixed.
- **Sets the bar.** Any fix of a weight-0 component must also be checked against the rest of the supplied
  structure (graph, sector, symmetry). G-invariance only makes a component eligible.
- **SCOUT-0 relation.** This *explains* EDA-01's FIX = 0 for the dimensionful rows and does not overturn it.
  EDA-01's READ rows that carry weight-0 numbers (E-12, E-13, E-14) remain READ: the numbers depend on a
  supplied lattice/sector/parent.

**Status: W1-C COMPLETE — NON-SELECTION THEOREM (family/row-checked scope). Eligible targets = G-invariants.
X-01 recorded (dimensionless ≠ eligible).**
