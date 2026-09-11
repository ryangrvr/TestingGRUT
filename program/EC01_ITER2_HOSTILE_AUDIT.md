# E&C-01 — Iteration 2 Hostile Audit

> Independent hostile review of the instrument-repair pass. Historical baselines
> (`5ab36d6` / `ccda3be` / `290ad5d`) verified untouched; frozen protocol SHA verified
> `6b012117…27a00a8`. **No commit, no push** (per dispatch). Iteration-2 artifacts
> audited but NOT modified.

## 0. Integrity

| Check | Result |
|---|---|
| Protocol SHA-256 | MATCH |
| Baselines 5ab36d6 / ccda3be / 290ad5d | UNTOUCHED |
| Commit/push during audit | NONE |

## 1. Repair-by-repair audit (R1–R7)

| Repair | Verdict |
|---|---|
| R1 M→P / P→O gates | **REPAIRED_AND_CONFIRMED** — 14/14 edges, content inspected not just presence |
| R2 Imported premises first-class | **REPAIRED_AND_CONFIRMED** — statistical postulate, Past Hypothesis, separability, entangling-H, tracing, matching inputs all flagged `NOT_DERIVED_FROM_SOURCE` |
| R3 mechanism_scope | **REPAIRED_AND_CONFIRMED** — all edges, narrowest justified; no UNKNOWN abuse |
| R4 A06 split | **REPAIRED_AND_CONFIRMED** — A06a phenomena verified individually; A06b UNRESOLVED survives attack (no trajectory mechanism exists in chain evidence) |
| R5 RC3 fresh evidence | **REPAIRED_WITH_RESIDUAL_DEPENDENCY** — see Finding H2 (source scope, not conclusion) |
| R6 RA1 WITHIN_FRAMEWORK | **REPAIRED_AND_CONFIRMED** — qualifier-leak scan: none found |
| R7 statistical/boundary/matching inputs | **REPAIRED_AND_CONFIRMED** |

## 2. Finding H1 — FB1 label too strong

`NO_DEMONSTRATED_MECHANISM` is too strong **globally**: Sinai-type mixing/ergodic
theorems are demonstrated mechanisms by which dynamics determines the equilibrium
measure for specific system classes. The equal-a-priori postulate is imported only
for the generic claim. Independent reclassification:

$$
\text{FB1v2}: \quad \text{CONDITIONAL\_DERIVATION\_WITH\_IMPORTED\_POSTULATE}
$$

The imported-postulate fact is retained; the label is narrowed. (Iteration 2's own
text recorded the partial ergodic support but kept the stronger label.)

## 3. Finding H2 — RC3 source scope (the most important attack)

The φ⁴₄ triviality citation **overstates**: rigorous triviality is proven for
d > 4 (Aizenman 1982; Fröhlich); for d = 4 it is widely accepted but rests on strong
partial results and nonperturbative lattice evidence — not a fully proven continuum
theorem in the cited form. **However the RC3v2 conclusion survives independently:**

1. **Order-counting/matching argument** (direct, framework-independent): EFT data at
   order N cannot distinguish UV contributions first appearing at order N+1.
2. **RG information loss** — framework-level feature of Wilsonian effective descriptions.
3. **Landscape multiplicity** (Bousso–Polchinski; Douglas–Kachru) — explicit multiplicity examples.

φ⁴ triviality is demoted from establishing citation to supporting/partial status.
`CONSTRAINT_NOT_DETERMINATION` and `FRAMEWORK_INDEPENDENT` scope **retained**, justified
via (1)+(2).

## 4. Gate audit

- **M→P**: all 4 DEMONSTRATED claims survive attack (entangling interactions, decoherence
  experiments, EFT limits, equilibrium realization). No status asserted without mechanism.
- **P→O**: survives; decoherence→fringe-loss and EFT→scattering correctly regime-limited;
  UV and selector observables correctly `NO_DISTINCTIVE_OBSERVABLE`.

## 5. Firewalls and provenance

- **Decoherence ≠ outcome ≠ Born derivation ≠ classical ontology: INTACT.** A07 remains
  isolated; the A06 split actively protects the distinction.
- **Parameter provenance**: "9 parameters, zero derived without import" **CONFIRMED**.
- **Traceability 11/11**: CONFIRMED — every row maps to a visible repair.
- **RA3** remains PARTIALLY_DEPENDENT (definition-relative, correctly not called circular).

## 6. Blind reclassification result

**13/14 edges agree** with the Iteration-2 labels under independent reclassification.
One revision (FB1v2 per H1). No new substantive errors; minor vocabulary notes
(CT-edge `EFFECTIVE` vs `NOT_APPLICABLE`; RA3 type/status pairing — schema-v3 candidates).

## 7. Verdict

$$
\boxed{\text{PASS\_WITH\_REVISIONS}}
$$

**Final answers to the dispatch questions:**

- *Did Iteration 2 actually repair the instrument, or merely add fields?* — **Genuinely
  repaired.** The evidence is structural: an apparent derivational arrow (FB1) was
  actually downgraded, and the A06 split exposed a real explanatory boundary (A06b,
  classical trajectories) rather than papering over it. More fields did not produce
  more certainty; they produced more precise uncertainty.

- *Which repaired classifications survive hostile reclassification?* — All except FB1v2,
  which is revised to `CONDITIONAL_DERIVATION_WITH_IMPORTED_POSTULATE`, and RC3v2, whose
  conclusion is retained but with one source (φ⁴₄) scope-narrowed.

## HARD STOP

Audit complete. **No commit. No push. No repair of the findings** — H1 and H2 belong
to the next iteration per the standing principle: a correction enters the next version;
it does not rewrite the evidence that revealed it.

Canonical GRUT-RAI: UNMODIFIED.
