# REALITY CHECK 04 — EXPERIMENT P INDEPENDENT ADJUDICATION

> Audit of `TestingGRUT @ 5febde9`. Canonical `GRUT-RAI` untouched.
> Source priority: executable code > raw results JSON > decision record > report.
> Independent verification performed: Controls A, D, H re-computed from scratch (numpy, exact finite-dimensional).

---

## 1. Method

- `calc/experiment_p_identifiability.py` read and re-executed under independent control.
- Control A re-verified numerically **and** analytically: the record-type environment (shift register, disjoint support) gives ⟨e₀|U₁|e₀⟩ = 0 by construction — decoherence is exact, not a numerical artifact.
- Control D scrutinized with a fresh independent calculation.
- Control H scrutinized with an independent audit of what observables are actually matched.

---

## 2. Status table

| Claim | Status |
|---|---|
| Closed unitary dynamics produce decoherence | **DERIVED** |
| Decoherence suppresses reduced coherence | **DERIVED** |
| Global unitary state selects one outcome | **NOT_DERIVED** |
| Pointer basis equals realized outcome | **NOT_DERIVED** |
| Outcome probabilities are dynamically derived | **NOT_DERIVED** |
| Born weights are derived | **NOT_DERIVED** (state-inherited) |
| Reduced state uniquely identifies global outcome info | **OUTSIDE_SCOPE as tested** (control D invalid) |
| Coarse-graining uniquely selects an outcome | **NOT_DERIVED** — control F: dephasing *manufactures* the outcome appearance |
| Additional outcome-selection structure required | **DERIVED** (control G) |
| Universal beyond tested model class | **OUTSIDE_SCOPE** (correctly not claimed) |

---

## 3. Control-by-control

| Control | Verdict | Note |
|---|---|---|
| A — decoherence | VERIFIED | cross term = 0.0 exactly; global purity 1.0; U unitary; analytic cross-check |
| B — definite outcome | VERIFIED / correct | global state remains pure entangled superposition; no selection op in code |
| C — weights vs amplitudes | VERIFIED | p₀ tracks |α|² exactly for 5 amplitudes → STATE_DEPENDENT |
| **D — purifications** | **INVALID AS EXECUTED** | see §4 |
| E — partitions | VERIFIED / weak | nested truncations, not independent partitions |
| F — coarse-grainings | VERIFIED | sharpest finding: dephased G₁ manufactures definite-outcome appearance |
| G — unitary vs collapse | VERIFIED / correctly labeled | MODEL C = ADDITIONAL_INPUT (selection + Born postulate) |
| H — identifiability | VERIFIED with scoping correction | see §5 |

---

## 4. CONTROL D IS INVALID AS EXECUTED — and the artifact contradicts itself

This is the most consequential audit finding.

Experiment P's own control D reports:

```
identical_reduced_decoherence: false
```

yet its verdict string simultaneously asserts:

```
"NON-IDENTIFIABLE at the reduced level; weights are STATE-DEPENDENT ... not dynamics-derived"
outcome_weight_identifiable_from_reduced_state_alone: true
```

That field is the **negation** of the verdict placed under it. Two internal contradictions.

**Mechanism (verified independently):** the second "purification" uses a generic environment superposition over record states. Under the same unitary, the reduced cross term is **0.163 ≠ 0** — i.e. the second state *does not decohere*. The pair is therefore **not an identical-reduced-state pair**, and the control never tests purificational indistinguishability at all.

**Consequence:** the theorem "reduced states do not determine purifications / outcome measures" is textbook-valid, but **this experiment does not demonstrate it**. All non-identifiability claims must be sourced from control H only, with the scoping correction below.

---

## 5. CONTROL H: scoping corrections

1. The compared reduced states have **different purities** (0.68 vs 0.5); only the pointer-basis decoherence class (cross = 0) is matched. Non-identifiability is therefore **relative to the declared decoherence-observable set**, not "identical observables" unqualified.
2. `p₀` is itself a reduced-state matrix element. If branch populations count as declared observables, the weights **are** identifiable from the reduced state — but as inherited initial-state information, not dynamics-derived.

**Correct precise statement:**

> Outcome weights are not identifiable from **decoherence structure** (coherence suppression) alone. They are identifiable from the full reduced state, but as initial-state input — hence STATE_DEPENDENT, not DERIVED.

---

## 6. Overreach list

| Location | Claim | Correction |
|---|---|---|
| control D verdict | "NON-IDENTIFIABLE at the reduced level" | Contradicts own computed flag; D invalid as executed |
| control H conclusion | "identical decoherence observables" | Only the decoherence class matched; purities differ |
| primary verdict | "non-identifiable from effective decoherence structure" | Accept as worded; reject generalized version elsewhere in artifact |

No universality overreach found — scope discipline is maintained throughout the original artifacts.

---

## 7. Relation to REALITY_CHECK_03

**CONSISTENT AND CLOSING.** RC03 identified: effective open-system structure plausibly derivable; outcome selection the sharpest gap. Experiment P confirms the first half (decoherence DERIVED) and establishes the second half as a scoped no-go (outcome selection NON-IDENTIFIABLE from decoherence structure in the tested class). The frontier is **closed at the tested model class**.

---

## 8. Final research status

> **Experiment P establishes decoherence but not outcome selection or Born-weight derivation within the tested model class** — with two audit corrections: (1) control D (purification non-identifiability) is invalid as executed and is removed from the evidence chain; (2) non-identifiability is scoped to the decoherence-observable set, not to the full reduced state.

`EXPERIMENT_P_FRONTIER = CLOSED_AT_TESTED_MODEL_CLASS`

---

## 9. UNITARY FORENSIC CORRECTION (APPLIED — supersedes conflicting claims above)

Verification artifact: `program/REALITY_CHECK_04_VERIFICATION.json`.

1. **Original model defect.** The original `calc/verify_experiment_p.py` U₁ was
   **not unitary** on its declared Hilbert space: dim = 16, `rank(U1) = 9`,
   `‖U1†U1 − I‖ = 1`; seven basis inputs had no image, seven outputs unoccupied.
2. **Restricted-subspace false positive.** The earlier reported zero unitarity
   error arose because the sampled initial states did not exercise the defective
   portions of U₁. It must not be cited as evidence of a valid closed-unitary model.
3. **Quarantine.** The original defective construction and outputs are preserved
   as quarantined history and carry no evidentiary weight.
4. **Cyclic repair.** Repaired U₁ is the complete permutation
   `U1|e_j> = |e_{(j+1) mod 2N}>` over the full record space
   (bijection verified: e₀→e₁, e₁→e₂, …, e₁₅→e₀).
5. **Unitarity verification.** Repaired U₁: `‖U1†U1 − I‖ = 0.0`,
   `‖U1U1† − I‖ = 0.0`, rank 16, all row/column norms 1. Full U (dim 32): rank 32,
   both unitarity errors 0.0 on the entire Hilbert space. Analytic argument: a
   bijection of an orthonormal basis is unitary.
6. **Regenerated A/B.** Control A: cross term exactly 0, global purity 1,
   branch populations (½, ½), record states E₀ = |e₁>, E₁ = |e₂> overlap exactly 0
   — **DERIVED**. Control B: global state remains the pure entangled superposition
   α|0>|e₀⟩ + β|1>|e₁⟩; no selection operation — **NOT_DERIVED**.
7. **Control D: INVALID_AS_EXECUTED** (repaired model: generic environment state,
   cross term 0.13849… ≠ 0). No evidentiary support for purificational
   non-identifiability. Original failure preserved.
8. **Control H: NARROWED.** Cross term = 0 for α ∈ {0.2, 0.5, 0.8} while
   p₀ = α² varies (0.04, 0.25, 0.64). Permitted claim only: *decoherence
   observables alone do not determine branch weights.* The stronger
   full-reduced-state claim is explicitly prohibited (p₀ is itself a
   reduced-state element).
9. **Weight origin.** p_i = |α_i|² classified STATE_DEPENDENT / INHERITED FROM
   INITIAL AMPLITUDES / NOT DYNAMICALLY DERIVED, within the tested model class.
   Not a universal disproof of the Born rule.

## 10. Corrected verdict

> **In the repaired model, closed unitary dynamics with orthogonal environmental
> records produce reduced-state decoherence while preserving a pure global
> entangled state. Definite outcome is therefore not derived by the tested
> unitary/coarse-graining route. Branch weights are inherited from the initial
> amplitudes and are not dynamically derived within this model class. Control D
> is invalid as executed and provides no evidence for purificational
> non-identifiability. Control H supports only the narrower claim that
> decoherence observables alone do not determine branch weights.**

Canonical GRUT-RAI: UNMODIFIED. All artifacts are TestingGRUT laboratory records.
