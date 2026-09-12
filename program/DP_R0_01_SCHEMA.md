# DP-R0-01 SCHEMA — required record structures (frozen)

## Model-variant record
```
variant_id
name                    (exact literature name)
notation_map            {symbol_in_source -> canonical R0} with justification
primary_source          citation + arXiv/DOI
first_introduced        year/source
dynamics_class          ORIGINAL | MARKOVIAN_WHITE_NOISE | SMEARED_DENSITY | OTHER (explicit)
r0_definition           exact mathematical role; units; where it enters
r0_status               one of the R0 provenance taxonomy values
limit_behavior          is the model well-defined as R0 -> {0, infinity}?
observability           does changing R0 change predictions?
not_equivalent_to       variant_ids with similar notation but different physics
```

## Claimed-motivation record
```
motivation_id
variant_id
claim                   verbatim source proposition (substring-extractable)
source                  full citation; PRIMARY only
source_type             PAPER | REVIEW | PROCEEDINGS
fixes_or_constrains     FIXES | CONSTRAINS | NEITHER
argument_type           DERIVATION | HEURISTIC | EFFECTIVENESS_CRITERION |
                        CONJECTURE | POST_HOC | CONVENTION | NONE
framework_relative      yes/no
post_hoc_status         justified_before_data | after_data | unresolved
```

## Constraint record
```
constraint_id
observable              quantitative definition
experiment              citation; setup
applicable_variants     [variant_id,...]  (never "all DP")
parameter_constrained   R0 | collapse_rate | lambda | other (explicit)
bound_statement         numerical + units
bound_type              one of the frozen bound-type values
statistical_statement   confidence/limit type if available
assumptions             environmental/model/regularization (list)
direct_or_inferred      DIRECT | INFERRED
action                  EXCLUDES | BOUNDS | DISFAVORS | UNRESOLVED
```

## Region record (per variant, after all constraints)
```
region_id
variant_id
status                  one of the 7 permitted outcomes
excluded_region         interval + which constraints exclude it
open_region             interval + status of each endpoint (HARD vs CRITERION vs HEURISTIC)
unconstrained_region    interval + why not meaningfully constrained
endpoint_epistemics     per endpoint: source + bound_type
```

## Observable-discrimination record
```
observable_id
quantitative_signature  formula + expected magnitude in open region
comparison              vs Gamma_env (standard decoherence budget) and vs other
                        collapse models; degeneracy noted if any
channel                 HEATING | SPONT_RADIATION | FORCE_NOISE | COHERENCE_LOSS |
                        INTERFEROMETRIC_VISIBILITY | MECHANICAL_NOISE | OTHER (defined)
confound_budget         thermal, gas, blackbody, EM, charge, laser recoil,
                        vibration, ordinary decoherence — each compared
distinguishable         YES | NO | DEGENERATE_WITH_ENVIRONMENT | UNRESOLVED
```

## R₀ selector record (Phase 8 only)
```
selector_id
candidate_category      one of the frozen selector taxonomy values
claim                   what is claimed to select R0
evidence_level          PROVEN_THEOREM | RIGOROUS_PARTIAL | CONCEPTUAL_ARGUMENT |
                        HEURISTIC | SUPPORTING_CONTEXT | NONE
determines_or_narrows   DETERMINES | NARROWS | DOES_NOT_SELECT
                        (formula containing R0 ≠ determination)
status                  per selector taxonomy; default NO_DEMONSTRATED_SELECTOR
                        absent affirmative evidence
```

## Field discipline
`NOT_SPECIFIED` / `UNRESOLVED` preferred to invented values. Every positive
claim carries `source + exact location + verbatim proposition`. Primary sources
only for load-bearing items; reviews may corroborate, never substitute.
