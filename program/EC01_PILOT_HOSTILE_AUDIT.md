# E&C-01 — HOSTILE PILOT AUDIT

Baseline: protocol `5ab36d6` (hash verified unchanged), pilot `ccda3be`.
Original pilot artifacts: **unmodified**. Audit only; no expansion; no rescue.

## The audit question

> Did E&C-01 genuinely classify the relationships in the pilot according to its frozen protocol, or did it merely rename the existing explanatory narrative?

**Verdict: PARTIALLY → PASS_WITH_REVISIONS.** The pilot did real classification work, with three required scope revisions and one object split. Details below.

## 1. Edge-level reclassification (original labels quarantined during review)

| Edge | Original | Auditor label | Verdict | Key attack |
|---|---|---|---|---|
| FA1 | GENERATION | GENERATION_WITH_EXPLICIT_PREMISES | scope narrowed | entanglement also available as initial data; separability-preserving interactions exist |
| FA2 | GENERATION | GENERATION_WITH_EXPLICIT_PREMISES | confirmed | entanglement ALONE is not sufficient (needs environment+coupling+trace); pilot already refuses this collapse |
| FA3 | COARSE_GRAINING | COARSE_GRAINING (scope narrowed) | **OBJECT_REQUIRES_SPLIT** | node A06 "classicality" bundles interference suppression, pointer stability, robust records, and classical *trajectories*; the evidence establishes the first three; trajectories need extra premises (near-classical Hamiltonians) |
| FB1 | GENERATION | GENERATION_WITH_IMPORTED_POSTULATE | scope narrowed | the statistical postulate is an inserted primitive, not produced by dynamics |
| FB2 | DERIVATION | DERIVATION | confirmed | thermodynamic-limit construction is the derivation; premises stated |
| FB3 | NO_DEMONSTRATED_MECHANISM | NO_DEMONSTRATED_MECHANISM (+IMPORTED_CONDITION) | confirmed | coarse-graining gives a statistical tendency only; Past Hypothesis is the imported condition |
| FC1 | COARSE_GRAINING | COARSE_GRAINING (parameter payload flagged) | confirmed | EFT also *introduces* new independent inputs (LECs), not merely discards information |
| FC2 | NO_DEMONSTRATED_MECHANISM | NO_DEMONSTRATED_MECHANISM | confirmed | multiple inequivalent UV completions remain live; non-uniqueness not rescinded |

**Tally: 4 confirmed as classified, 2 scope-narrowed, 1 object split required, 0 downgraded, 0 upgraded.**

## 2. The primary attack surface: "classicality"

The word compresses at least four distinguishable phenomena:

1. suppression of interference (decoherence does establish this);
2. pointer-state stability / einselection (established given partition + coupling);
3. robustness of records (established);
4. recovery of classical **trajectories** (needs additional dynamical assumptions not stated in FA3's premises);
5. classical *ontology* (NOT established — Experiment P).

FA3's demonstrated scope covers (1)–(3). Recommendation (audit output only; frozen pilot untouched): split `A06` into `A06a` (decohered effective statistics / pointer stability — established) and `A06b` (classical trajectories — separate edge with its own premises or UNRESOLVED).

## 3. Reverse-map independence test

| Reverse edge | Independence |
|---|---|
| RA3 (decoherence→entanglement REQUIRED) | **PARTIALLY_DEPENDENT** — circularity risk: if decoherence is *defined* as entanglement-with-environment + trace, "requires entanglement" is partly imported from the forward definition. Survives as within-framework only. |
| RA4 (interference suppression→decoherence REQUIRED) | Independent — grounded in empirical observation. |
| RA2 (partition underdetermination) | Independent — genuine gap, correctly kept UNDERDETERMINED. |
| RA1 (→ unitary dynamics, "within framework") | Partially dependent; framework-dependence must be carried forward. |
| RB1 (→ low-entropy condition REQUIRED) | Independent — structural, not imported. |
| RB2 (substrate UNDERDETERMINED) | Independent. |
| RC3 (→ UV non-uniqueness) | **INHERITED_BOUNDARY** — leans on RC-02 without fresh primary-source anchoring in the pilot. Retained with flag. |
| RC4 (CONSTRAINT_WITHOUT_GENERATION) | Independent. |

**No full circularity found; one partial dependence (RA3) and one inherited boundary (RC3) flagged.**

## 4. Object-identity consistency

- `A06 classicality`: drift detected (see §2). **Split flagged.**
- `thermodynamics` vs `arrow of time`: correctly kept distinct (pilot does not repeat the state-space/entropy conflation).
- `EFT`: used consistently as a construction procedure with declared cutoff/order. No drift.

## 5. Historical vs explanatory

No forward edge upgrades chronology to causation. Reverse edges correctly distinguish REQUIRED / ONE_OF_MULTIPLE / UNDERDETERMINED. PASS.

## 6. Stopping-point audit

| Chain | Stop | Status |
|---|---|---|
| A | partition + definite outcome | **EVIDENCE_DRIVEN** (repaired, unitarity-verified Experiment P) |
| B | Past Hypothesis | **EVIDENCE_DRIVEN** (independently re-derived here) |
| C | UV completion | PARTIALLY_INHERITED (RC-02) — flag retained |

## 7. M→P / P→O gates

The pilot's edge records do **not** carry separate realization/observation gate fields, which the frozen protocol requires for major transitions. This is a protocol shortfall of the pilot execution — flagged for the next E&C iteration. No gate was asserted without support, so no wrong classification results from it.

## 8. Parameter provenance re-audit

Statistical postulate (B): POSTULATED, not derived — correctly carried. Wilson coefficients (C): EMPIRICAL_INPUT/FITTED — correctly carried. Coupling (A): model input. Branch weights (A): INHERITED per the P theorem — correctly never called dynamically derived.

## 9. H1 screen recheck

Agreed with the pilot: "coarse-graining with information loss recurs" is a property of our constructions, not evidence of a physical substrate. Each chain still imports genuine primitives at its stop (statistical postulate; Past Hypothesis; UV non-uniqueness). H1 remains **UNDETERMINED**.

## 10. Discrepancy matrix

```text
edges reviewed:            8
confirmed as classified:   4
scope narrowed:            2
object split required:     1
downgraded / upgraded:     0 / 0
reverse independence:      5 independent, 2 partially dependent, 1 inherited
overall verdict:           PASS_WITH_REVISIONS
```

## 11. Final answer to the audit question

**YES, with revisions.** E&C-01 did not merely rename the narrative: it refused the three canonical collapses (entanglement→decoherence alone, decoherence→classical ontology, decoherence→definite outcome), carried honest premises on every generative edge (including the inserted statistical postulate), and produced genuine reverse-map gaps. The revisions required are scope refinements — the A06 split, premise foregrounding on FA1/FB1, and fresh provenance for RC3 — not structural failures.

A discovered flaw is a successful result of the audit. Three were found.

---
Artifacts: `EC01_PILOT_HOSTILE_AUDIT.json` (machine-readable), this file.
Frozen protocol SHA verified: `6b012117…27a00a8`. Canonical GRUT-RAI: UNMODIFIED. Pilot artifacts: UNMODIFIED.
