# RC05_CASE_02 — Standard finite decoherence construction (external, non-GRUT)

> REALITY_CHECK_05 calibration case; frozen protocol, unmodified.

## Sources

- **E. Joos, D. Zeh, Z. Phys. B 59, 223 (1985)** — reduced coherence suppression under environmental coupling `CANONICAL_PRIMARY`
- **W. H. Zurek, Phys. Rev. D 26, 1862 (1982); Phys. Rev. D 24, 1516 (1981)** — pointer-basis selection and environment-induced superselection `CANONICAL_PRIMARY`
- **W. H. Zurek, Rev. Mod. Phys. 75, 715 (2003)** — decoherence selects pointer states; outcome problem explicitly NOT solved by decoherence alone `CANONICAL_PRIMARY`

## Input ledger

- H = H_S + H_E + H_int (system-environment model) → ASSUMED (declared construction)
- initial separable/weak-coupled state → ASSUMED
- interaction Hamiltonian (e.g. pointer coupling) → ASSUMED (declared)
- Born weights of initial amplitudes → INHERITED_STATE_DATA — NOT dynamically derived

## Model reconstruction

rho_S(t) = Tr_E[U rho(0) U^dagger]; off-diagonal elements decay as exp(-t/tau_D) in the pointer basis; diagonal populations = |alpha_i|^2 (inherited).

## Provenance map

claim -> master equation -> declared system/bath model -> primary sources. Decoherence is a THEOREM of the declared unitary model (partial trace of unitary evolution).

## Identifiability

Decoherence rate and pointer basis are uniquely determined by declared inputs. Definite realized outcome is NOT uniquely determined: global state remains an entangled superposition; multiple admissible purifications/global descriptions compatible with the same reduced state.

## Prediction analysis

Decoherence observable (coherence suppression) is DERIVED and testable. Definite outcome: NOT_DERIVED. No probability law derived beyond inherited amplitudes.

## Effective status

EFFECTIVE — reduced dynamics of an open subsystem of a closed unitary whole.

## Adversarial attack

Is decoherence 'fundamental'? No: it is derived for the reduced description; the global state remains pure. Does diagonal rho imply an outcome? Not within the model. Objection did not overturn the classification.

## Final classification

```json
{
  "decoherence": "DERIVED + EFFECTIVE",
  "pointer_basis": "DERIVED (basis selected by declared interaction)",
  "definite_outcome": "NOT_DERIVED",
  "outcome_weights": "INHERITED / NOT DYNAMICALLY DERIVED"
}
```

Confidence: HIGH | Instrument pass: True

Limitations: Sources canonical-primary, not re-fetched in this environment. Construction choice (pointer coupling) is declared, not derived — matches protocol expectation.
