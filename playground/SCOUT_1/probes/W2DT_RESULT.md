# SCOUT-1 W2-DT RESULT — dimensional transmutation / anomalous scale generation vs TC-1

**Charter:** `PROBE_CHARTERS.md` §W2-DT (pre-registered).

**Files:** `w2dt_dimensional_transmutation.py`, with log `w2dt_dimensional_transmutation.log`.

**Labels:**
- TC1-SURVIVES-ANOMALY, with a stated reformulation of G;
- KNOWN RESULT IMPORT: Gross–Neveu large-N 1974; Coleman–Weinberg 1973; one-loop QCD running — STANDARD-TEXTBOOK ✓;
- a genuine **NET REDUCTION** of dimensionless input (see §3).

## 0. Verdict

> **TC1-SURVIVES-ANOMALY — and the probe found the first true information reduction of the campaign.**
>
> 1. **A scale is generated.** Large-N Gross–Neveu is classically scale-invariant, with one dimensionless
>    coupling λ. Quantum mechanically, the trace anomaly breaks scale symmetry and a mass
>    `m = Λ/√(e^{2π/λ} − 1)` appears.
> 2. **The scale is RG-invariant, and its absolute value is not fixed.**
>    - Cutoffs Λ = 10¹ … 10¹⁶ with λ(Λ) on one trajectory give m = 1.0000000 (gap equation and direct
>      minimization agree to 10⁻⁷).
>    - Which trajectory is chosen is one boundary datum: λ at a reference cutoff. That datum is a
>      **G-breaking reference** — a dimensionless coupling attached to a dimensionful reference scale.
>    - QCD shows the same thing: `Λ_QCD (1-loop, n_f = 5) = 87.8 MeV` needs `α_s = 0.118` **at μ = M_Z**.
>    - Coleman–Weinberg shows it too: `⟨φ⟩` is the renormalization point at which λ takes its CW value.
> 3. **But every dimensionless observable is fixed with no remaining dimensionless input.**
>    - `ΔV/m² = −1/4π` (−0.07957747) for all Λ from 10⁴ to 10¹⁶, and for all λ from 0.2 to 0.5 (λ = 1.0 shows a
>      finite-cutoff 10⁻⁴ deviation).
>    - The classical theory had **one** free dimensionless parameter. The quantum theory has **zero**, plus one
>      scale that is merely the unit.
>    - Coleman–Weinberg likewise eliminates λ: `m_S²/m_V² = 3e²/8π²` is fixed by the remaining coupling.

## 1. TC-1 status — survives, with G reinterpreted

**Is G anomalous?** The classical scaling G is **anomalous**: a pure rescaling is not a symmetry of the quantum
theory at fixed λ. It is a symmetry of the **RG-covariant family**: rescaling combined with the flow of λ
(Callan–Symanzik).

**Reformulation.** Replace G by its RG-covariant action on the family of quantum theories. Then:
- the whole family is a **single G-orbit**, coordinatized by m;
- the absolute m is G-moved, and no G-invariant predicate fixes it (TC-1 holds);
- G-invariants (all mass ratios, `ΔV/m²`, bound-state ratios) are **all fixed**.

**Outcome:** TC1-SURVIVES-ANOMALY.
- An absolute scale still requires a G-breaking boundary datum: a coupling at a reference scale, or a unit
  supplied by a second scale.
- No TC1-COUNTEREXAMPLE.
- The "needs reformulation" branch applies only to the extent that G must be read as the RG-covariant action.
  The theorem's logic is unchanged. **The TC-1 text gets the RG-covariant reading as a standing clarification.**

## 2. Where the absolute scale enters (the three distinctions)

| Distinction | Gross–Neveu / QCD / CW | Input needed |
|---|---|---|
| breaking classical scale symmetry | yes (anomaly) | none: automatic in the quantum theory |
| RG-invariant scale | yes (m, Λ_QCD, ⟨φ⟩) | none: it exists for every trajectory |
| absolute numerical scale in physical units | **no** | coupling at a reference scale (α_s at M_Z), **or** a second dimensionful scale to measure against |

**Hierarchy.** `Λ/M_Pl = exp(−2π/b₀α_UV)` = 1.6·10⁻¹⁸, 5.8·10⁻¹⁵, 1.4·10⁻¹² for α_UV = 0.020, 0.025, 0.030.
- A small *dimensionless* boundary datum produces a huge ratio.
- A 1 % change in α_UV moves the hierarchy by ×1.39.
- So transmutation explains *hierarchies* (ratios between two scales) from a dimensionless datum. It never
  produces the first scale.

## 3. Information accounting (input for ZOOM_OUT_03)

- **TRUE NET REDUCTION OF INPUT INFORMATION:** classical {λ} → quantum {} dimensionless parameters, plus one
  unit. This is not an exchange: the dimensionless coupling is **consumed** by the anomaly and reappears only as
  the choice of unit, which TC-1 says is unfixable anyway.
- **Price:** an interacting, asymptotically free (marginally relevant) sector with a non-zero β-function.
- **GRUT relevance:**
  - The current earned floor is **quadratic** (E-3 linearity, exact reductions). Free theories have no running,
    no anomaly and no transmutation.
  - So the one mechanism that genuinely removed input information is **absent from GRUT's earned layer by the
    same linearity that makes its reductions exact**. This is the second E-3 tension (after W2-ETH).
  - A GRUT route to "fewer inputs" would need a marginally relevant interacting sector. That is a premise change
    (owner level), recorded as the top Wave-3 target, not entered.

**Status: W2-DT COMPLETE — TC1-SURVIVES-ANOMALY (G read as RG-covariant scaling). Dimensional transmutation removes
a dimensionless input (net reduction), never an absolute scale. Absent from GRUT's quadratic earned layer.**
