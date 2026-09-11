# REALITY_CHECK_05 — FROZEN INSTRUMENT VALIDATION PROTOCOL

> **Frozen before any calibration case is executed.** This file is the instrument
> specification. Once created, its substantive methodology is NOT to be modified
> in response to any calibration result. Ambiguities/defects discovered later are
> recorded as *protocol limitations*, not retroactively patched. Revisions belong
> to a future version (REALITY_CHECK_05.1+).
>
> Object under test: **the audit methodology**, not GRUT. Canonical GRUT-RAI untouched.

---

## 1. Purpose

REALITY_CHECK_02 validated the *taxonomy* (CLOSED/CONDITIONAL/PARTIAL/INPUT/
UNRESOLVED) on five external frameworks. It did NOT validate the full epistemic
machinery: input ledger → model reconstruction → control → identifiability →
provenance tracing → epistemic classification → adversarial attack → scoped claim.

RC-05 tests the missing layer, on theories the instrument did not originate.
The instrument must recognize **both** legitimate positive and legitimate negative
results. A methodology that rejects everything is broken; one that accepts
everything is broken.

## 2. Core epistemic vocabulary

Labels (NOT always mutually exclusive — multiple may apply simultaneously):

```
DERIVED / ASSUMED / EFFECTIVE / EMPIRICAL_INPUT /
NON_IDENTIFIABLE / UNRESOLVED / PREDICTIVE / CONSTRAINT
```

## 3. Required distinctions (must never be collapsed)

```
derived != assumed              derived != empirical input
effective != fundamental        consistent != derived
derived != predictive           prediction != postdiction
decoherence != outcome selection
pointer-basis selection != realized outcome
mathematical representation != physical ontology
model-class result != universal theorem
identifiability != mere calculability
agreement with experiment != derivation of the underlying quantity
```

## 4. Per-case procedure

1. **Input ledger** — every structure, constant, coupling, state, symmetry,
   regulator, truncation, fit parameter, probability postulate, interpretive
   assumption; each classified by provenance. Recursive: trace calculated
   quantities to authorized primitives, empirical inputs, or assumptions.
2. **Model reconstruction** — actual mathematical object, equations, inputs,
   calculated quantities, target observable, claimed scope. If not
   reconstructible from adequate sources: UNRESOLVED. No gap-filling from
   standard assumptions unless the case explicitly authorizes it.
3. **Provenance chain** — claim → equation → inputs → source. Classes:
   SOURCE-STATED / DERIVED-FROM-SOURCE / INFERRED / UNRESOLVED. AI interpretation
   is not evidence.
4. **Identifiability** — is the quantity uniquely determined by declared
   model+inputs? Non-identifiability requires an explicit admissible
   counterexample or rigorous argument — NOT mere multiplicity of
   representations/interpretations/difficulty/free parameters.
5. **Prediction test** — PREDICTIVE only if: fixed by model before comparison;
   inputs independently declared; not fitted to target; no hidden adjustable
   parameter chosen using the target; observable specified; quantitatively
   testable; derivation traceable. A derived quantity may be DERIVED but NOT
   PREDICTIVE if not independently fixed.
6. **Effective-theory test** — effective variables, underlying variables,
   reduction, regime, truncation, matching inputs. Do not penalize effectiveness.
   The correct result may be EFFECTIVE + DERIVED + PREDICTIVE WITHIN SCOPE.
7. **Adversarial attack** — attempt to falsify own classification. Change the
   verdict only if the evidence changes, not because a question exists.

## 5. Calibration set — FROZEN before execution

| # | Case | Pre-registered expectation |
|---|---|---|
| 01 | QED anomalous magnetic moment (a=(g−2)/2) | DERIVED + PREDICTIVE; renormalized inputs EMPIRICAL_INPUT. Instrument must NOT reject merely because QED has empirical inputs. |
| 02 | Standard finite decoherence construction | Reduced decoherence: DERIVED + EFFECTIVE. Definite outcome: NOT_DERIVED. Must reproduce the Experiment-P boundary WITHOUT originating from GRUT — critical cross-validation. Must not classify it as a GRUT result. |
| 03 | Born rule in Everett | POSTULATED / DERIVATION-CONTESTED. Do NOT force DERIVED; do NOT force FAILED; do not choose a philosophical winner. Report the dispute and the differing assumptions between competing derivations. If the vocabulary cannot represent this: record as instrument limitation. |
| 04 | GRW / objective collapse | Collapse dynamics: ADDITIONAL_INPUT / ASSUMED. Outcome selection: provided by the model's postulated dynamics. Probability law: explicitly locate where it enters. Do not classify collapse as derived merely because the model yields definite outcomes. Additional structure may solve the problem — it is still not derived from unitary QM. |
| 05 | Weinberg–Witten-constrained emergent gravity | CONSTRAINT: a scoped no-go ruling out a specified class under stated assumptions. Not a new positive theory, not an unscoped universal, not an empirical input. Record theorem assumptions and scope. |
| 06 | Chiral perturbation theory (specific quantitative observable) | EFFECTIVE + DERIVED within declared inputs/order + PREDICTIVE within regime/truncation (where prediction criteria are met). The instrument must not demand fundamental closure before recognizing a valid prediction. |

## 6. Success criteria (all required)

- **A. Positive-result recognition** — correctly recognizes QED g−2 (DERIVED+PREDICTIVE) and chiral EFT (EFFECTIVE+DERIVED+PREDICTIVE within scope).
- **B. External reproduction of the P boundary** — decoherence DERIVED/EFFECTIVE vs definite outcome NOT_DERIVED, on a non-GRUT construction.
- **C. Additional-input recognition** — GRW collapse recognized as postulated structure, not derived.
- **D. Dispute preservation** — Everett/Born disagreement preserved, not adjudicated without pre-registered criteria.
- **E. Constraint recognition** — scoped no-go recognized as a legitimate result.
- **F. No skepticism bias** — no systematic downgrading of established results because of declared empirical inputs. Rejecting g−2 or valid EFT predictions for having inputs = instrument failure.
- **G. No credulity bias** — no promotion of consistency, interpretive plausibility, postdiction, fitted agreement, or mathematical possibility into prediction/derivation.

## 7. Case-selection and source discipline

Do not substitute easier cases. Do not remove a contested case because its
literature is disputed. If a case is genuinely unreconstructable: UNRESOLVED,
with explanation — do not replace it after seeing the difficulty.

Secure authoritative primary/canonical sources BEFORE classification (for
disputed cases, multiple sources representing the actual positions). Record
source, exact section/page/equation, claim supported, scope/assumptions.

## 8. Case record format

Each case produces: case_id, subject, source_ledger, model_reconstruction,
input_ledger, provenance_map, identifiability_analysis, prediction_analysis,
effective_scope, adversarial_attack, final_classification, confidence,
limitations, evidence_locations. **Expected classifications are never used as
evidence.**

## 9. Instrument self-audit and verdict

After all six cases, blind-compare PRE-REGISTERED EXPECTATION vs OBSERVED
classification. For each mismatch determine: instrument failure / case
ambiguity / source ambiguity / legitimate disagreement / protocol limitation.
Do not silently correct the instrument.

**RC-05 = PASS** only if: all six cases reconstructed sufficiently to permit
judgment; required positive and negative discriminations correctly made; no
systematic skepticism or credulity bias; genuine interpretive disagreement
preserved; effective-but-predictive theory recognized as such.
Otherwise **FAIL** (a valid result) or **INCONCLUSIVE** if the evidence
cannot support a binary verdict.

## 10. Frozen

This protocol was written before any case research. It contains no
case-result information beyond the pre-registered expectations above.
No case may be executed until this protocol is reviewed and explicitly
authorized. No commit. No push.
