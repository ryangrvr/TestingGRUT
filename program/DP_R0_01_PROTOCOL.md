# DP-R0-01 PROTOCOL — Diósi–Penrose R₀ parameter-space reconstruction (FROZEN)

> Frozen before execution. No DP analysis has been performed under this
> protocol. No commit/push in this run. Canonical GRUT-RAI untouched.

## Governing questions (strict order)

Q1 (empirical/definitional): **What DP parameter space, if any, remains after
variant-specific and provenance-clean constraints?**

Q2 (explanatory, only AFTER Q1): **Is R₀ independently selected, constrained,
or merely fitted?**

## Hard prohibitions

- Do NOT assume the window `4 Å < R₀ < 10⁶ Å` exists. Both endpoints are
  source-relative: the lower bound carries experimental/model assumptions;
  the upper figure is a phenomenological/effectiveness criterion, not data.
- Do NOT treat "DP" as one model. Variant conflation is a protocol failure.
- Do NOT rescore the comparative compression matrix in this run.
- Do NOT repair DP. Do NOT invent a new collapse model.
- Do NOT call an experimentally allowed range a prediction.
- No parameter may be called DERIVED because a numerical value is consistent
  with experiment; no parameter is PREDICTED if chosen after seeing data.

## Execution phases (may not be reversed)

1. Model-variant registry (see MODEL_REGISTRY; extend only by documented
   primary-source variants; never merge notation without recording the map).
2. Per-variant R₀ definition: symbol, units, exact mathematical role, where it
   enters the dynamics, regulator vs physical parameter, observability of
   changing it, well-definedness in limits.
3. Claimed-motivation reconstruction: exact claim, primary source, whether it
   FIXES vs CONSTRAINS R₀, heuristic vs derivation, post-hoc status.
4. Constraint registry, source-by-source: observable, experiment, applicable
   variant, statistical statement, assumptions, bound type.
5. Variant-specific allowed/excluded regions. Permitted outcomes:
   `NO_SURVIVING_REGION | CONDITIONAL_SURVIVING_REGION |
   EMPIRICALLY_OPEN_REGION | INDEPENDENTLY_CONSTRAINED_SCALE |
   FREE_PHENOMENOLOGICAL_PARAMETER | DISTINCTIVE_TEST_AVAILABLE |
   UNRESOLVED_PARAMETER_SPACE`.
6. Distinguishing-observable identification (quantitatively defined channel
   only: heating, spontaneous radiation, force noise, coherence loss,
   interferometric visibility, mechanical noise). "Quantum behavior changes"
   is not an observable.
7. Environmental-confound comparison: thermal, gas, blackbody, EM/charge
   noise, laser recoil, vibration, ordinary decoherence budget.
8. R₀ selector taxonomy (only now):
   `DYNAMICAL_SELECTION | MICROPHYSICAL_DERIVATION |
   SYMMETRY_OR_CONSISTENCY_CONSTRAINT | EFFECTIVENESS_SELECTION |
   BOUNDARY_OR_INITIAL_CONDITION | EMPIRICAL_CONSTRAINT |
   FRAMEWORK_RELATIVE_SELECTION | POST_HOC_FIT |
   NO_DEMONSTRATED_SELECTOR | UNRESOLVED`.
   A surviving interval is NOT a selected scale.
9. Prediction test: parameter fixed independently + observable pre-specified +
   quantitative consequence + distinguishable from QM/environment + no tuning.
10. Independent hostile audit (mandatory; must be capable of finding the DP
    branch weaker OR a viable region surviving).

## Bound-type taxonomy (never merged)

`HARD_EMPIRICAL_BOUND | CONDITIONAL_EMPIRICAL_BOUND |
THEORETICAL_CRITERION | EFFECTIVENESS_CRITERION | HEURISTIC_BOUND |
MODEL_CLASS_RESULT | UNKNOWN`

## Epistemic firewall (inherited from E&C-01 / RC-05)

`DERIVED ≠ ASSUMED ≠ EMPIRICAL_INPUT ≠ EFFECTIVE ≠ PREDICTIVE`
`CONSISTENT ≠ DERIVED` · `CONSTRAINT ≠ DETERMINATION`
`EXPERIMENTAL_AGREEMENT ≠ THEORETICAL_DERIVATION`
`MOTIVATION ≠ DERIVATION` · `CONSTRAINT ≠ SELECTION`
`POSTDICTION ≠ PREDICTION`

## R₀ provenance taxonomy

`DERIVED | EMPIRICAL_INPUT | PHENOMENOLOGICAL_PARAMETER | REGULATOR |
CONVENTION | CONJECTURAL_PHYSICAL_SCALE | EFFECTIVE_PARAMETER |
CONDITIONALLY_CONSTRAINED | UNRESOLVED`

## No-fitting firewall

If a proposed R₀ is chosen because it makes an existing experiment work, it is
`POST_HOC_FIT` unless an independent prior justification exists.

## Frozen-protocol discipline

This file is the instrument specification. Do not modify it after seeing DP
results. Ambiguities/defects get recorded as protocol limitations; revisions
belong to DP-R0-01.1 (future, separate). The execution must record this
protocol's SHA-256 and verify it unchanged at the hostile audit.

## Relation to prior artifacts

Inputs: `REALITY_CHECK_05_DP_AUDIT.*`, `REALITY_CHECK_05_PARAMETER_SPACE.*`,
`REALITY_CHECK_05_R0_DETERMINATION.*` (status: audit artifacts, not premises —
all their numbers are re-derived source-by-source under this protocol).
Compression matrix remains frozen and untouched.
