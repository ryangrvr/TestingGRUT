# SCOUT-1 W1-P RESULT — does complete passivity fix the forcing-law temperature profile? (C12)

**Charter:** `PROBE_CHARTERS.md` §W1-P. Preregistered outcomes: SELECTION PRINCIPLE (premise-priced) / NOT FIXED.

**Files:** `w1p_complete_passivity.py`, with log `w1p_complete_passivity.log`. States are diagonal in the energy
basis, so ergotropies are exact by sorting (up to 2²² levels).

**Labels:**
- SELECTION PRINCIPLE (premise-priced);
- KNOWN RESULT IMPORT: Pusz–Woronowicz 1978; Lenard 1978 — STANDARD-TEXTBOOK ✓ (checked);
- REDISCOVERED-KNOWN: the "large bath acts as its own copies" refinement;
- NON-DISTINCTIVE.

## 0. Verdict

> **SELECTION PRINCIPLE (premise-priced): `P ⇒ Q ∈ S` with S drastically smaller.**
>
> **What gets fixed.** Complete passivity (no work extractable from any number of copies by cyclic unitaries)
> collapses Q2's mode-temperature profile `T(ω)` — a function — to a **single number T**. Equivalently, it
> fixes the **shape** of the noise/dissipation ratio, `coth(ω/2T)` as a function of `ω/T`, which is a
> weight-0 component.
>
> **What stays free.** The scale T is G-moved, so W1-C forbids any invariant principle from fixing it. The
> only selectable values are the G-fixed points:
> - `T = 0` (ground state, β = ∞);
> - `T = ∞` (β = 0).
>
> Both are completely passive, as the theorem requires.
>
> **New refinement.** For a **macroscopic bath with a dense mode spectrum, single-copy passivity already
> suffices.** The allowed profile window shrinks as m grows: 1.0 at m = 2, 0.15 at m = 4, 0.0073 at m = 20.
> In the dense limit, Cauchy additivity of `f(ω) = ω/T(ω)` forces f to be linear, i.e. KMS.
>
> **The price.** The premise is a statement about the **bath preparation** (S-11), which is supplied. Complete
> passivity is the Kelvin–Planck second law imposed as a state condition. It is physically defensible, but it
> is **not** a consequence of the dynamics: a non-equilibrium bath is an allowed preparation that is
> simply a work resource.

## 1. Numerical results

| State | single-copy `W₁` | `W_N/N` for N = 1…9 | asymptotic per-copy (`E − E_Gibbs` at equal S) |
|---|---|---|---|
| KMS T = 1 (ω = 1, 2) | 0 | 0 (≤ 2·10⁻¹⁶) | 0 |
| KMS T = 0.3 (three modes) | 0 | 0 | 0 |
| ground state | 0 | 0 | 0 |
| `T(ω)`: 1.0, 1.5 (`ω/T` = 1, 1.33) | **0 (passive)** | 0, 0, 1.6·10⁻³, …, 7.0·10⁻³ | 1.14·10⁻² |
| `T(ω)`: 1.0, 0.8 (`ω/T` = 1, 2.5) | **0 (passive)** | 0 up to N = 6; 5.7·10⁻⁶ at N = 7, 6.3·10⁻⁵ at N = 9 | 3.15·10⁻³ |
| three modes, T rising slowly | 5.4·10⁻⁴ (not passive) | rising to 1.05·10⁻³ | 2.03·10⁻³ |
| inverted (`ω/T` = 1, 0.67) | 7.0·10⁻² | ≈ 7.0–7.2·10⁻² | 7.48·10⁻² |

**Other checks:**
- **Two-qubit single-copy passivity** holds iff `ω₁/T₁ ≤ ω₂/T₂`. This matched 2000 random cases.
- **Many modes, one copy** (`T(ω) = T₀(1 + ε(ω−1))`, modes evenly spaced in [0.5, 1.5]). The maximal passive
  ε is:

  | m | 2 | 3 | 4 | 6 | 8 | 12 | 16 | 20 |
  |---|---|---|---|---|---|---|---|---|
  | max passive ε | 1.0 | 1.0 | 0.150 | 0.077 | 0.052 | 0.015 | 0.010 | 0.0073 |

  The window shrinks roughly as `m⁻²`.

## 2. Ten-point hostile test of "complete passivity ⇒ single T"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | The profile shape cannot vary. The scale T can. |
| 2 | P silently contains Q? | **Partly.** Complete passivity is equivalent to "Gibbs at one β" (that is the theorem). The principle is a restatement of equilibrium **as a no-resource condition**, which is the physically motivated form. |
| 3 | Representation-dependent? | It needs a declared bath Hamiltonian, i.e. the energy basis. The choice of H for the bath is supplied. |
| 4 | Physical or gauge? | Physical: ergotropy is measurable work. |
| 5 | Standard? | Yes: PW/Lenard. The dense-bath refinement is the infinite-system form of the same theorem. |
| 6 | Parent variation? | Qubits, three-mode and many-mode baths: all pass. Oscillator modes are untested but covered by the theorem. |
| 7 | Composition? | Complete passivity is defined **through** composition (copies). Single-copy passivity is **not** closed under composition: the passive two-mode states become active at N = 3 or N = 7. |
| 8 | Coarse-graining? | Coarse-graining modes into blocks keeps KMS (a Gibbs state restricts to Gibbs). Non-KMS passive states can turn active when blocked (the m-dependence above). |
| 9 | Unique or stationary? | Unique family (β ∈ [0, ∞]). |
| 10 | Boundary condition selecting? | No. But the endpoints β = 0, ∞ are the G-fixed points: the only values a scale-free principle could pick. |

## 3. What this means

- **Second genuine `P ⇒ Q ∈ S`** (after W1-A). It reduces a *function's worth* of freedom (Q2's `T(ω)`, the
  P-17 converse freedom) to one scale, **conditional on a supplied preparation principle**.
- **Pattern TP-1 again, now with a thermodynamic layer.** The selector is a supplied state-level layer (bath
  preparation = equilibrium) composed with a supplied bath Hamiltonian. Its output (a single β) is a
  G-invariant shape; the G-moved scale T stays free.
- **Record cross-check (conditional, classical-analog only).** S2 admits site noise profiles `T_i`. These
  include remote-only profiles such as G(∞) (`T₁ = 0`, remote `T > 0`), which are non-equilibrium.
  - A passivity principle carried over to the classical limit (uniform `T_i`) would **remove** such profiles.
  - It would **not** remove S2-1's leading discriminator `Δc₂ = −24βT₁a`, which needs only `T₁ > 0`.
  - The classical translation (site temperatures ↔ mode temperatures) is **not proven here**.
- **Empirical consequence:** FDT with one temperature. It is standard and non-distinctive.

**Status: W1-P COMPLETE — SELECTION PRINCIPLE (premise-priced: supplied preparation). Profile → single T; T
free (W1-C). Macroscopic-bath refinement (single-copy passivity suffices in the dense limit).**
