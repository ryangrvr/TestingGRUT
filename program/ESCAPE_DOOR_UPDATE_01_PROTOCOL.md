# ESCAPE_DOOR_UPDATE_01 — FROZEN PROTOCOL (method before execution)

> Status: FROZEN before any execution. This file must not be modified after commit-freeze.
> If a defect is found later, record it as a protocol limitation; corrections belong to version 2.
> Canonical GRUT-RAI untouched. TestingGRUT laboratory only. No Born-rule analysis authorized by this protocol.

## 1. Immutable inputs

- `program/REALITY_CHECK_05_MATRIX.json` / `.md` (committed RC-05 comparative results)
- `program/DP_RESCORE_01.json` / `.md` (commit 0ca550c)
- `program/DP_R0_01_EXECUTION.json` + hostile audit (commit 699fb78)
- Frozen RC-05 scoring rules (from REALITY_CHECK_05_PROTOCOL.md)

No other input may contribute scores. No candidate may be added or removed after freeze.

## 2. Immutability rules

1. **No score changes** unless the frozen RC-05 rules mechanically require them from committed evidence. Any such change must cite the rule and the evidence artifact.
2. **No promotion/demotion by judgment.** Every axis value must trace to a committed artifact.
3. Baseline (standard QM/QFT + decoherence) retains **unranked/reference** status. It is never forced into the ranking; it is the comparison floor.
4. Everett 0.82 is recorded as **highest frozen-heuristic compression score among compared alternatives — NOT empirical support**. This wording is mandatory wherever the value appears.
5. DP 0.63 is the audited rescored value; a future reversion is possible only via a new evidence-driven rescoring under the same rubric (not by this update).
6. The five malformed comparative-foundations JSON cases are NOT repaired, read as evidence inputs only to the extent already committed; their status is unchanged.

## 3. Axes

AXIS 1 — COMPRESSION SCORE: value from the frozen-rubric matrix (DP = 0.63 audited; others as committed).
AXIS 2 — EMPIRICAL DISCRIMINABILITY: assessed independently on four sub-fields, each classified separately:

- **D1 — test channel exists in principle** (quantitatively defined observable; YES/NO/UNRESOLVED)
- **D2 — parameter region experimentally constrained** (bounds exist that act on the candidate; YES/NO/UNRESOLVED)
- **D3 — clean discriminator demonstrated** (separation from environmental/systematic confounds actually demonstrated in the literature for a currently accessible region; YES/NO/UNRESOLVED)
- **D4 — successful novel prediction exists** (quantitative prediction fixed before comparison and confirmed; YES/NO/UNRESOLVED)

AXIS 2 DISCRIMINABILITY LEVEL (frozen operational definitions, assigned mechanically):

| Level | Definition |
|---|---|
| HIGH | D1 = YES **and** D3 = YES for a currently accessible parameter region |
| MEDIUM | D1 = YES and D2 = YES, but D3 = NO or region accessibility unresolved |
| LOW | D1 = YES but D2 = NO (no constraining evidence acting on the candidate) |
| NONE | D1 = NO |
| UNRESOLVED | evidence insufficient to assign any of the above |

Assign each candidate exactly one LEVEL + the four sub-field values.

## 4. Quadrant classification (frozen)

Using AXIS 1 relative to baseline and AXIS 2 LEVEL:

- **HIGH/HIGH** — compression ≥ baseline-heuristic comparison is favorable AND LEVEL = HIGH
- **HIGH/LOW** — favorable compression, LEVEL = LOW or NONE
- **LOW/HIGH** — unfavorable compression, LEVEL = HIGH
- **LOW/LOW** — unfavorable compression, LEVEL = LOW or NONE
- **UNRESOLVED** — either axis cannot be assigned under these rules

No candidate may be labeled a "winner" in the map. Quadrants describe position, not merit beyond the two axes.

## 5. Headline hypothesis (tested, not assumed)

Pre-register the claim to be mechanically tested:

> **H1-headline:** No audited candidate simultaneously exceeds baseline and possesses a demonstrated clean empirical discriminator (D3 = YES with accessible region).

The update must attempt to falsify H1-headline. It is CONFIRMED / REFUTED / UNRESOLVED strictly by the computed table. No wording may soften or strengthen the result beyond the classification.

## 6. Mandatory wording guards

- "test channel in principle" ≠ "demonstrated discriminator"
- "empirically constrained" ≠ "novel prediction confirmed"
- "highest compression score" ≠ "empirical support"
- Everett: mandatory qualifier attached (Sec. 2.4)
- DP: 0.63 audited; radiation channel = test vehicle; confound separation not yet demonstrated across surviving window

## 7. Outputs

- `program/ESCAPE_DOOR_UPDATE_01.json`
- `program/ESCAPE_DOOR_UPDATE_01.md`

Contents: per-candidate axis values + four D-subfields + quadrant + evidence citations; baseline row (unranked reference); H1-headline verdict; discrepancies vs prior map; limitations.

## 8. Prohibitions

- No Born-rule analysis (E&C-BORN-01 is a separate future decision).
- No malformed-case repair; no cosmetic-rewrite cleanup; no debris staging.
- No new candidates; no rescoring of other candidates beyond frozen-rule requirements.
- No canonical GRUT-RAI modification.

## 9. Stop condition

When both artifacts exist and every cell carries an evidence citation: STOP. Execution and any commit require separate authorization.
