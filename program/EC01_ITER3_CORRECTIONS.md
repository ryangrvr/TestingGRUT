# File: program/EC01_ITER3_CORRECTIONS.md
# E&C-01 ITERATION 3 — TARGETED CORRECTIONS (H1 / H2)

> Scope-locked correction pass. Historical controls immutable:
> `5ab36d6` (protocol) / `ccda3be` (pilot) / `290ad5d` (audit 1) / `03a18e9` (iter 2 + audit 2).
> Protocol SHA-256 verified unchanged. Canonical GRUT-RAI untouched.

---

## Correction H1 — FB1

**Original (iter 2):** `NO_DEMONSTRATED_MECHANISM`
**Finding (03a18e9):** Too strong — Sinai-type mixing/ergodic results supply demonstrated
mechanisms for restricted system classes.
**Independent reclassification:** `CONDITIONAL_DERIVATION_WITH_IMPORTED_POSTULATE`
**Targeted re-audit verdict:** CONFIRMED.

What is now recorded:

- **Demonstrated mechanism (restricted classes):** chaotic mixing results (Sinai billiards;
  hard-sphere gas under stated conditions) yield microcanonical-type statistical behavior
  within those classes. Typicality (Goldstein–Lebowitz–Tumulka–Zanghì) justifies equal
  weights over the microcanonical shell for large systems.
- **Imported premise (still first-class):** the equal-a-priori/typicality measure remains
  necessary for the GENERIC claim. Dynamics alone does not fix the measure in general.
- **mechanism_scope:** `REGIME_LIMITED` — the mechanism is demonstrated only within mixing/
  typicality classes.
- **relation_scope:** `WITHIN_FRAMEWORK`.
- **M→P:** `CONDITIONAL` — chaotic many-body systems are physically realizable, but no
  selector is demonstrated for generic systems.
- **P→O:** `NOT_APPLICABLE` — measure provenance remains unresolved (unchanged).
- **Firewall preserved:** specific ergodic theorem ≠ generic thermodynamics; mixing ≠ universal
  equilibration; ergodicity ≠ Past-Hypothesis explanation; conditional theorem ≠ derivation
  from dynamics alone.

## Correction H2 — RC3 / φ⁴₄

**Finding (03a18e9):** φ⁴₄ citation overstated; rigorous triviality is a d>4 theorem.

Source-by-source audit (see `EC01_ITER3_SOURCE_AUDIT.json`):

| Level | Result |
|---|---|
| d > 4 | PROVEN_THEOREM (Aizenman 1982; Fröhlich 1983) |
| d = 4 | RIGOROUS_PARTIAL_RESULTS + NUMERICAL_LATTICE_EVIDENCE (Lüscher–Weisz 1987 et seq.) — NOT the same-level theorem |
| Final status of φ⁴₄ | **SUPPORTING_CONTEXT** (demoted from load-bearing) |

RC3 is retained as `CONSTRAINT_NOT_DETERMINATION` on grounds independent of the
overstated citation, via three separately audited lines:

- **Line A (structural):** order-counting/matching — finite-order EFT data cannot
  distinguish UV completions differing beyond that order.
- **Line B (framework-level):** Wilsonian RG information loss — multiple UV theories
  flow to the same IR EFT.
- **Line C (FRAMEWORK_RELATIVE):** landscape multiplicity (Bousso–Polchinski 2000;
  Douglas–Kachru 2007) — concrete multiplicity within string compactifications.

**Gates preserved:** M→P `UNRESOLVED`; P→O `NO_DISTINCTIVE_OBSERVABLE`;
mechanism_scope `GENERIC_UNDER_STATED_PREMISES`; relation_scope `FRAMEWORK_INDEPENDENT`
(structural) / `FRAMEWORK_RELATIVE` (line C).

## Locality & dependency results

| Edge | Depends on | Change |
|---|---|---|
| FB2v2 | statistical postulate (carried) | unchanged |
| CT2v2 | references FB1v2 | unchanged (form-constraint independent of measure status) |
| FC2v2 | references RC3v2 | unchanged (classification retained) |

**Corrections are locally contained. Zero unrelated edges altered. Zero out-of-scope
findings recorded.**

## Principle

> A correction enters the next version; it does not rewrite the evidence that revealed it.
