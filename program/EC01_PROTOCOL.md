# E&C-01 — EMERGENCE & CONSTRAINT RECONSTRUCTION PROTOCOL (FROZEN)

> Version 1.0 — frozen before any epoch population, case work, literature
> reconstruction, or substrate search. No research data are embedded in this
> document. Canonical GRUT-RAI untouched. All work in TestingGRUT.

---

## 0. CENTRAL AUDIT RULE

**A ⇏ B unless the mechanism connecting A and B is demonstrated.**

This is an EPISTEMIC/AUDIT RULE about what the program may claim. It is NOT a
metaphysical assertion that reality cannot contain unexplained transitions.
Unknown mechanism ⇒ record `UNRESOLVED_RELATION` or `NO_DEMONSTRATED_MECHANISM`.
Do not fill gaps by intuition, analogy, chronology, elegance, terminology, or
theoretical preference.

## 1. SCOPE

E&C-01 is an auditable reconstruction-and-failure program investigating how
mathematical, dynamical, and physical structures relate across physical
regimes. It is not a theory. The working hypothesis about organization —

    possibility → constraint → compatibility → realization → new possibility → ...

— is a hypothesis about the reconstruction's organization, not evidence that
reality works this way.

## 2. GENERAL TRANSITION FORM

Do NOT assume `M_{n+1} ⊆ M_n`. Use:

    (M_n, P_n, C_n) --F_n--> (M_{n+1}, P_{n+1}, C_{n+1})

- M_n = mathematical description/admissibility space
- P_n = physical regime / instantiated structure
- C_n = constraints relevant at that stage
- F_n = proposed transition map/mechanism (must be classified, not assumed)

## 3. TRANSITION TAXONOMY (FROZEN — no new types after freeze)

1. `HISTORICAL_SUCCESSION` — A occurred before B. Chronology only; never implies A ⇒ B.
2. `REPRESENTATIONAL_MAPPING` — descriptive/mathematical mapping; no physical generation claimed.
3. `EQUIVALENCE` — explicit appropriate equivalence map preserving relevant structure/observables.
4. `DUALITY` — demonstrated dual description under an established mapping. Analogy is insufficient.
5. `COMPATIBILITY` — at least one joint realization of A and B exists. Does not imply generation.
6. `CONSTRAINT_TRANSMISSION` — A restricts admissible possibilities for B:
   B(A) ⊆ B. Stronger than compatibility, weaker than generation.
7. `COARSE_GRAINING` — B obtained by specified elimination/averaging/tracing of degrees of freedom.
8. `SYMMETRY_BREAKING` — prior symmetry broken explicitly or spontaneously under specified dynamics/conditions.
9. `GENERATION` — specified mechanism operating on A produces B. Mechanism must be stated.
10. `DERIVATION` — B follows uniquely/necessarily from stated premises, equations, conditions, allowed operations.
11. `ONTOLOGICAL_CLAIM` — B is claimed physically real, not merely descriptive. Requires independent justification; never inferred from mathematical usefulness.
12. `NO_DEMONSTRATED_MECHANISM` — suggested/suspected connection, no supported mechanism.
13. `UNRESOLVED_RELATION` — relation cannot presently be classified.

Partial-order discipline: the taxonomy is a taxonomy of RELATION TYPE, not a
linear strength scale. An exact representational map can be stronger than a
vague generative claim; an empirical constraint can be decisive for one
question without being generative.

## 4. GATES (each independently passable / failable)

### 4.1 Object gate
Type endpoints before comparing: syntactic / semantic / interpretive /
unresolved typing / conflicting referents; distinguish mathematical object,
physical object, representation, observable, state-dependent vs dynamical
structure, effective variable, ontological candidate. Insufficiently typed
endpoint ⇒ **BLOCK THE RELATION**. Do not infer type from surrounding theory.

### 4.2 Realization gate (two separate questions)
- `M → P`: why is an admissible structure physically instantiated?
- `P → O`: what observable distinguishes the instantiated structure from alternatives?

Realization selectors (frozen, no additions after freeze):
`DYNAMICAL_SELECTION | STABILITY_SELECTION | COMPATIBILITY_SELECTION |
INITIAL_OR_BOUNDARY_CONDITION | EXTERNAL_PHYSICAL_PRINCIPLE |
OBSERVER_RELATIVE_SELECTION | NO_DEMONSTRATED_SELECTOR | UNRESOLVED`

## 5. RECONSTRUCTION (two independent maps)

- **Forward:** E0 → E1 → ... → En (chronological/historical/effective succession).
- **Reverse:** En → E_{n-1}? → E_{n-2}? → ... (what earlier structure is actually
  required to explain what is observed now?).

The reverse map must NOT be the forward map read backwards; construct it
independently. Compare only after both are built. If they disagree, DO NOT
REPAIR. Record as result: `HISTORICAL_CONTINUITY_WITHOUT_EXPLANATORY_CONTINUITY`,
`REVERSE_RECONSTRUCTION_GAP`, or appropriate unresolved status.

## 6. EQUIVALENCE / DUALITY DISCIPLINE

Classify as `SAME_CONSTRUCTION | EQUIVALENT_REPRESENTATION | DUAL_DESCRIPTION |
OBSERVABLY_EQUIVALENT | RELATED_BUT_NOT_EQUIVALENT | UNRESOLVED`.
Similar equations, single-point agreement, shared experiment, or verbal
translation do not suffice. Explicit mapping + domain of validity required.

## 7. PARAMETER PROVENANCE (frozen categories)

`PRIMITIVE | DERIVED | EFFECTIVE | EMPIRICAL_INPUT | FITTED | STATE_DEPENDENT |
BOUNDARY_CONDITION | ENVIRONMENT_DEPENDENT | CONVENTION_DEPENDENT | UNRESOLVED`.
Numerical expressibility is not derivation; track genuine determination.

## 8. ALTERNATIVE-TRANSITION ACCOUNTING

For every nontrivial explanatory edge, ask whether an alternative mechanism
produces the same result. Record alternatives. Possible mechanism ⇒ no
explanatory credit. Observationally equivalent mechanisms ⇒
`OBSERVATIONAL_DEGENERACY`.

## 9. SUBSTRATE HYPOTHESIS H1 (preregistered, downstream, NOT assumed)

**H1:** There exists a common underlying relational structure R from which major
physical regimes can be generated without independent primitive insertion at
every transition.

Permitted outcomes: `ONE_SUBSTRATE | MULTIPLE_RELATED_SUBSTRATES |
ONE_SUBSTRATE_PLUS_INDEPENDENT_BOUNDARY_CONDITIONS |
MULTIPLE_INDEPENDENT_PRIMITIVE_SECTORS | UNRESOLVED_SUBSTRATE_STRUCTURE`.

Sequence: reconstruct → identify recurring structure → test one substrate.
Never construct R first. Never reinterpret everything as 1space.

**H1 failure/weakening conditions:**
1. Every major transition needs a new unrelated primitive.
2. No common mechanism generates both quantum and geometric structure.
3. "Common substrate" is relabeling of target structures.
4. R ≡ Q + G + T + ... (target-collection renaming, no independent content).
5. No nontrivial shared constraints.
6. No reduction of independent inputs.
7. No consequence distinguishing it from disconnected effective descriptions.
8. Arbitrary parameter insertion at every transition.
9. Common structure only at representational level, no demonstrated physical content.
10. Reconstruction reveals genuinely independent primitive sectors with no demonstrated maps.

## 10. PROVISIONAL FUNDAMENTALITY

Status: `PROVISIONALLY_FUNDAMENTAL_UNDER_CURRENT_MAPS` — granted only when a
structure cannot currently be eliminated by demonstrated transition maps
without reintroducing an equivalent or independent primitive.

Safeguards: model-relative; not metaphysical proof; equivalent representations
are not elimination; "not yet derived" is insufficient; downgradable by later
demonstrated maps.

## 11. OUTPUT ARCHITECTURE (accounting, not ontological ladder)

    OBSERVED → FORCED → EFFECTIVE → GENERATED → ASSUMED
    → PROVISIONALLY_FUNDAMENTAL_UNDER_CURRENT_MAPS

Every entry preserves its actual evidence status.

## 12. RECORD SCHEMAS

Transition record: transition_id, source_object, target_object, source_type,
target_type, edge_type, premises, mechanism, mathematical_map,
physical_interpretation, parameter_inputs, parameter_provenance,
observable_consequence, alternative_mechanisms, evidence, failure_condition,
status, confidence, scope, unresolved_items.

Object record: object_id, name, mathematical_type, physical_role,
representation_status, state_dependence, dynamical_status, effective_status,
provenance, observables, alternatives, unresolved_identity, confidence, scope.

Missing field ⇒ `NOT_SPECIFIED` or `UNRESOLVED`. Never manufacture values.

## 13. EPOCH DISCIPLINE

No predefined number of epochs. Candidate areas for investigation only
(quantum structure; geometry/spacetime; gauge/matter; symmetry breaking;
QCD/confinement; nuclear; atoms/chemistry; condensed matter; thermodynamics;
classicality; cosmology; complexity; observers). Distinct emergence events are
not assumed; chronological transitions are not assumed explanatory.

## 14. STRONG EXPLANATORY ARROW — REQUIREMENTS

Typed endpoints; explicit premises; explicit mechanism; parameter provenance;
nontrivial constraint/generation; observable consequence; alternative analysis;
falsification/weakening condition. Historical association, conceptual analogy,
shared vocabulary, or successful effective theory alone do not qualify.

## 15. REALITY-FIRST OUTCOMES (all legitimate)

GRUT wrong; 1space wrong; emergence not fundamental; multiple primitives;
proposed transition nonexistent; only historical relation; effective
description sufficient; deeper structure currently inaccessible; substrate
question undecidable from current observables. The goal is not unification but
an increasingly accurate map of what reality forces versus what our
descriptions assume.

## 16. FREEZE VERIFICATION CHECKLIST (performed at freeze time)

1. Protocol contains no case/reconstruction data.
2. No GRUT conclusion encoded as a protocol premise.
3. H1 is a hypothesis, not an assumption.
4. Transition taxonomy complete for the intended protocol and frozen.
5. `UNRESOLVED_RELATION` and `NO_DEMONSTRATED_MECHANISM` are legal terminal outcomes.
6. `MULTIPLE_INDEPENDENT_PRIMITIVE_SECTORS` permitted.
7. Forward and reverse reconstruction explicitly independent.
8. Realization and observation are separate gates.
9. `PROVISIONALLY_FUNDAMENTAL_UNDER_CURRENT_MAPS` is not metaphysical FUNDAMENTAL.
10. No epochs populated.
11. No physics calculated.
12. No existing research artifact modified.

## 17. GOVERNING DISCIPLINE (per arrow)

claim → typed endpoints → edge type → premises → mechanism → observable
consequence → failure condition. Missing mechanism ⇒ the edge stops at
`UNRESOLVED_RELATION`. Follow every arrow; if the arrow stops, stop at the gap.
