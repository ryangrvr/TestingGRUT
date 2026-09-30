# RC05_CASE_03 — Born-rule status in Everettian quantum mechanics (dispute-preservation calibration)

> REALITY_CHECK_05 calibration case; frozen protocol, unmodified.

## Sources

- **H. Everett III, Rev. Mod. Phys. 29, 454 (1957)** — universal wavefunction; branches; measure introduced `CANONICAL_PRIMARY`
- **D. Deutsch, Proc. R. Soc. A 455, 3123 (1999); D. Wallace, 'The Emergent Multiverse' (2012)** — decision-theoretic derivation of Born weights from rationality axioms `CANONICAL_PRIMARY`
- **W. H. Zurek, Phys. Rev. Lett. 106, 250402 (2011)** — envariance-based derivation of Born weights `CANONICAL_PRIMARY`
- **T. Maudlin, in 'Many Worlds?' (2010); A. Kent, Int. J. Mod. Phys. A 5, 1745 (1990)** — critiques: derivations rely on nontrivial assumptions (additivity, decision continuity, preferred factorization); Everett does not yield Born rule unaided `CANONICAL_PRIMARY`

## Input ledger

- universal unitary dynamics, no collapse → ASSUMED (Everett postulate)
- branch decomposition / preferred factorization → STRUCTURAL_SELECTION — arguably emergent via decoherence, contested as fundamental
- decision-theoretic rationality axioms (Deutsch-Wallace) → ASSUMED (explicitly declared in that derivation)
- envariance + Schmidt decomposition (Zurek) → ASSUMED (explicitly declared in that derivation)

## Model reconstruction

Everett: |Psi> branches after decoherence; question is whether p_i = |alpha_i|^2 is DERIVED or POSTULATED. Three published routes: (a) Everett 1957 measure postulate (explicit POSTULATE); (b) Deutsch-Wallace decision theory (DERIVATION under declared rationality axioms); (c) Zurek envariance (DERIVATION under declared symmetry assumptions).

## Provenance map

Every route's axioms are explicitly declared in its primary source. None is derived from unitary dynamics alone without declared additional axioms.

## Identifiability

The Born weights are not uniquely determined by unitary dynamics alone; each derivation requires a distinct additional axiom set. Different admissible axiom sets yield the same Born rule but the necessity of those axioms is disputed.

## Prediction analysis

No Everett-internal prediction distinguishing Born from non-Born weights is agreed upon; derivations are PREDICTIVE only relative to their declared axioms.

## Effective status



## Adversarial attack

Could the instrument pick a winner (e.g. 'Deutsch-Wallace derives it')? Refused: critics (Maudlin, Kent) contest the axiom necessity, and the protocol forbids adjudicating without pre-registered criteria. Could it force NOT_DERIVED? Refused: derivations exist under declared assumptions. Both refusal verdicts upheld on adversarial review.

## Final classification

```json
{
  "born_rule_status": "DERIVATION_CONTESTED / ASSUMPTION_DEPENDENT",
  "representable_in_frozen_vocabulary": "YES (ASSUMPTION_DEPENDENT via multi-label classification; flagged as a vocabulary stress point)"
}
```

Confidence: HIGH | Instrument pass: True

Limitations:
