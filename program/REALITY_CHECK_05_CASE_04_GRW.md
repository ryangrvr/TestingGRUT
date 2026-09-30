# RC05_CASE_04 — GRW / objective collapse (additional-input recognition calibration)

> REALITY_CHECK_05 calibration case; frozen protocol, unmodified.

## Sources

- **G. Ghirardi, A. Rimini, T. Weber, Phys. Rev. D 34, 470 (1986)** — spontaneous localization process with hit rate lambda and localization width a `CANONICAL_PRIMARY`
- **G. Ghirardi, P. Pearle, A. Rimini, Phys. Rev. A 42, 78 (1990)** — continuous spontaneous localization (CSL), relativistic-consistent stochastic master equation `CANONICAL_PRIMARY`
- **P. Pearle, Phys. Rev. A 39, 2277 (1989)** — parameterized stochastic collapse dynamics `CANONICAL_PRIMARY`

## Input ledger

- hit rate lambda (~10^-16 s^-1 per nucleon) → ASSUMED — phenomenological parameter, postulated, not derived from unitary QM
- localization width a (~10^-7 m) → ASSUMED — phenomenological parameter
- stochastic (Poisson) collapse law → ASSUMED — POSTULATED DYNAMICS (the model's defining additional input)
- probability measure over collapse trajectories → POSTULATED — the GRW stochastic law explicitly introduces it; not derived from unitary evolution

## Model reconstruction

Between hits: standard Schrodinger evolution. Hits: localization operator L_i(x) applied with Poisson-distributed rate. Master equation (CSL form) with Lindblad-type decoherence-plus-collapse term.

## Provenance map

outcome selection -> stochastic collapse law -> GRW postulates. The collapse term is the model's DEFINING ASSUMPTION, declared in the primary source.

## Identifiability

Outcome selection IS provided: the stochastic dynamics makes definite outcomes (to stochastic accuracy). Probability law: explicitly POSTULATED within the model (the Poisson/Itô measure), not derived.

## Prediction analysis

GRW/CSL produces testable deviations from unitary QM (spontaneous heating, macroscopic superposition bounds) — PREDICTIVE but with fitted/postulated parameters; experimental windows constrain lambda, a. Prediction status: PREDICTIVE_WITH_POSTULATED_PARAMETERS.

## Effective status

NOT an effective description of unitary QM; a rival fundamental stochastic dynamics.

## Adversarial attack

Would the instrument call GRW 'derived from QM' because its equations generate outcomes? No: the equations differ from unitary QM by postulated stochastic terms. Confirmed: classification survives attack.

## Final classification

```json
{
  "collapse_dynamics": "ASSUMED / ADDITIONAL_INPUT",
  "outcome_selection": "PROVIDED_BY_ADDITIONAL_MODEL_STRUCTURE",
  "probability_law": "POSTULATED (explicit stochastic measure)",
  "prediction": "PREDICTIVE_WITH_POSTULATED_PARAMETERS"
}
```

Confidence: HIGH | Instrument pass: True

Limitations: Sources canonical-primary, not re-fetched.
