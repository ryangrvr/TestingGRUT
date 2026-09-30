# E&C-01 ITERATION 2 — REPAIRED RECONSTRUCTION SCHEMA

> Instrument repair only. No new chains, no H1 test, no substrate construction.
> Historical baselines `5ab36d6` (protocol), `ccda3be` (pilot), `290ad5d` (hostile audit)
> are immutable. Every repair is traceable to a hostile-audit finding.

## Governing rule (unchanged, audit rule not axiom)

$$A \not\Rightarrow B \text{ unless the mechanism connecting } A \text{ and } B \text{ is demonstrated.}$$

## Repair 1 — Realization/Observation gates (audit finding: mp_o_gates)

Every Iteration 2 edge carries TWO independent gate records:

```json
"m_to_p": {
  "status": DEMONSTRATED|CONDITIONAL|EMPIRICAL_INPUT|EFFECTIVE|
            NO_DEMONSTRATED_SELECTOR|UNDERDETERMINED|UNRESOLVED|NOT_APPLICABLE,
  "mechanism": "<how math structure becomes physically instantiated>",
  "selector": DYNAMICAL_SELECTION|STABILITY_SELECTION|COMPATIBILITY_SELECTION|
              INITIAL_OR_BOUNDARY_CONDITION|EXTERNAL_PHYSICAL_PRINCIPLE|
              OBSERVER_RELATIVE_SELECTION|NO_DEMONSTRATED_SELECTOR|UNRESOLVED,
  "premises": [],
  "evidence": "",
  "scope": "",
  "unresolved": ""
},
"p_to_o": {
  "status": <same vocabulary>,
  "observable": "",
  "mapping": "<how physical structure maps to observable>",
  "measurement_conditions": [],
  "alternatives": [],
  "evidence": "",
  "scope": "",
  "unresolved": ""
}
```

Rule: the two gates are filled INDEPENDENTLY. `M→P demonstrated` never implies
`P→O demonstrated`, and vice versa. A missing gate is a result, not an omission.

## Repair 2 — Imported premises as first-class data (finding FB1)

Every edge carries:

```json
"imported_premises": {
  "state_assumptions": [],
  "interaction_assumptions": [],
  "environment_assumptions": [],
  "coarse_graining_assumptions": [],
  "boundary_conditions": [],
  "initial_conditions": [],
  "statistical_postulates": [],
  "matching_inputs": [],
  "parameter_inputs": []
}
```

Each imported premise is classified `NOT_DERIVED_FROM_SOURCE` where it is required
by the target but not produced by the source node. Empty lists are `[]`; unknown is
`["UNKNOWN"]`; not relevant is `["NOT_APPLICABLE"]`.

## Repair 3 — Mechanism scope (new field)

```json
"mechanism_scope": MODEL_SPECIFIC | FRAMEWORK_RELATIVE | REGIME_LIMITED |
                   GENERIC_UNDER_STATED_PREMISES | UNKNOWN
```

Use the narrowest justified value. GENERATION/DERIVATION labels never silently
expand beyond this scope.

## Repair 4 — A06 split (finding FA3 / object identity drift)

- **A06a** `effective classical statistics / robust records`: interference
  suppression, pointer stability, record robustness, classical statistics.
  Established by decoherence evidence (P + EID).
- **A06b** `classical trajectories`: approximate Newtonian trajectory recovery.
  Requires ADDITIONAL premises (semiclassical Hamiltonians, localized states,
  Ehrenfest/quasi-periodic dynamics) — carried explicitly, otherwise UNRESOLVED.
- **A07** `definite realized outcome` remains a separate node; never merged.

## Repair 5 — RC3 fresh re-anchoring (finding: INHERITED_BOUNDARY)

Reverse chain-C boundary re-derived from primary/authoritative sources
(non-uniqueness of UV completions given a low-energy EFT), no citation of RC-02
as evidence. See EC01_ITER2_MAPS.json RC3 record.

## Repair 6 — WITHIN_FRAMEWORK qualifier (finding RA1/RA3)

```json
"relation_scope": WITHIN_FRAMEWORK | FRAMEWORK_INDEPENDENT | UNKNOWN
```

Every reverse-map edge carries it. Framework-relative necessity never propagates
into substrate/ontology claims.

## Edge schema (full)

```json
{
  "edge_id": "", "prior_edge_id": null, "chain": "",
  "source": "", "target": "", "edge_type": "<frozen taxonomy>",
  "mechanism_scope": "", "relation_scope": "",
  "premises": [],
  "imported_premises": { ... },          // Repair 2
  "mechanism": "",
  "parameter_inputs": [], "parameter_provenance": [],
  "m_to_p": { ... }, "p_to_o": { ... },  // Repair 1
  "alternatives": [],
  "evidence": [], "failure_condition": "",
  "status": "", "confidence": "", "scope": "",
  "note": ""
}
```

Transition taxonomy: the frozen protocol's 13 types only. Qualifications
(`GENERATION_WITH_EXPLICIT_PREMISES` etc.) appear as notes + mechanism_scope,
never as new edge_types.
