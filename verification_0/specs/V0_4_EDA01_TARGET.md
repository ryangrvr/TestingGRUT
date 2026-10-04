# V0-4 TARGET SPEC — EDA-01 SECOND-READER AUDIT (blind, classification-first)

**What this spec contains:** the audit question, the classification rubric, the row source, the reading rules and the
acceptance tests. **It does NOT contain the original row classifications, the original verdict reasoning, or the
original witnesses.**

**Source of the target** (frozen `scout-0 @ ab2da47`): `playground/SCOUT_0/probes/EDGE_DATA_AUDIT_01.md` (blob
`a9acfd75…`).
- The orchestrator read lines 1 – 39 and the full classification table (lines 40 – 63), and grepped the headings.
- It also read `BASELINE_MAP.md` Table 1 row stubs.

**The orchestrator is therefore EXPOSED to the original classifications. The reproducer must be context-isolated.**

## 1. The audit question (as posed in the original)

For every **earned result E-1 … E-22** (`playground/SCOUT_0/BASELINE_MAP.md`, Table 1): does it **independently FIX** a
local spectral edge datum relevant to an effective law, or does it only state a consequence once that datum is supplied?

**Edge data in scope:**
- edge location / spectral bottom;
- edge exponent (local density-of-states exponent);
- edge velocity at an occupation edge;
- visible spectral sign;
- number and location of soft points.

The reader may add further data classes (e.g. low spectral moments) if it finds them relevant, but must say so
explicitly.

## 2. Classification rubric (labels from the original; operational definitions fixed here)

| label | name | meaning |
|---|---|---|
| **FIX** | FIXES-EDGE-DATA | the earned result, **with no supplied primitive varied**, determines the **value** of an edge datum relevant to an effective law |
| **CON** | CONSTRAINS-BUT-DOES-NOT-FIX | the earned result restricts the admissible edge data (existence, form, a selection among supplied options), but leaves values free |
| **READ** | READS-OUT-EDGE-DATA | the earned result is a statement whose content is **a function of** an edge / spectral datum that is supplied upstream |
| **IND** | INDEPENDENT | no spectral / edge datum is involved (selection terminals, non-implications, architecture, ontology) |
| **SUP** | SUPPLIED-UPSTREAM | the result transports or restates data that are inputs |
| **NF** | NOT-FORMULABLE | the question cannot be posed for that entry |

- A dual label (e.g. IND / SUP) is allowed if justified.
- **Decisive test for FIX vs READ / CON:** exhibit (or rule out) two models that satisfy the same earned predicate but
  have different values of the datum. If such a pair exists, the entry does **not** FIX it.
- For every READ / CON / SUP entry, trace the datum **upstream** to a supplied primitive (use the SUPPLIED list,
  `BASELINE_MAP.md` Table 3, labels S-1 …).

## 3. Mandatory special check (VER0 context)

**E-13 (SF-1 FORMATION-OF-LAW-CLASS)** must be classified against the **corrected** edge-data quotient banked by VER0
V0-3, **not** the old single-scalar formulation. The VER0 V0-3 replacement (owner-accepted) is:

> The relevant quotient is an edge-data set: positions and orientations of the relevant occupation edges, plus their local
> dispersion orders / jets. The dynamical exponent and the soft-point structure are separate invariants. Local support
> exponents depend on the first non-zero derivative order r_e at the relevant edges; soft momenta depend on the
> inequivalent allowed edge differences; multiple pockets cannot be reduced to one v_b; the scalar v_b = 0 / ≠ 0
> classification is only a special compression for the generic single-pocket cosine-like case.

**Headline question:** does **any** earned result genuinely **FIX** such data, or does it only **read / constrain /
import** them?

## 4. Reading rules (firewall)

**ALLOWED:**
- this spec;
- `playground/SCOUT_0/BASELINE_MAP.md` (entire file);
- primary record documents cited by Table 1 or Table 3 (owner rulings, verdicts, charters at the repository top level or
  under `playground/`);
- `playground/SCOUT_0/probes/P02_RESULT.md`, for E-13 only (note: its single-v_b theorem was refuted by VER0; use the
  §3 replacement);
- other `playground/SCOUT_0/probes/*_RESULT.md` files, **except** those on the deny-list.

**DENIED (they contain or quote EDA-01's classifications):**
- `playground/SCOUT_0/probes/EDGE_DATA_AUDIT_01.md`;
- `playground/SCOUT_0/probes/eda_witnesses.py` and `.log`;
- `playground/SCOUT_0/FRONTIER_QUEUE.md`;
- `playground/SCOUT_0/SATURATION_CERTIFICATION_01.md` and `SATURATION_CERTIFICATION_02.md`;
- `playground/SCOUT_0/SCOUT_0_HANDOFF.md`;
- `playground/SCOUT_0/ZOOM_OUT_03.md` and `ZOOM_OUT_04.md`;
- `playground/SCOUT_0/probes/P02B_RESULT.md`;
- **anything** under `verification_0/` other than this spec and the reader's own output files;
- any other branch or working tree (`/home/user/BRI0`, `/home/user/GRUT-SCOUT-2`, `/home/user/TestingGRUT`);
- git commands.

If an allowed file turns out to quote EDA-01 classifications, stop reading it, report it, and continue.

## 5. Acceptance outputs

1. **A row table for E-1 … E-22.** Each row gives: the label; the edge / spectral datum involved (if any); the upstream
   supplied primitive; and a one- or two-line justification **including the FIX-test reasoning** (a two-model pair, or a
   reason none can exist).
2. **At least three explicit two-model witnesses** for the most important READ / CON rows: same earned predicate,
   different supplied primitive, different edge datum, different effective law. Numerical where cheap; fresh code only.
3. **The headline verdict:** is any row FIX? If yes, give a full proof sketch.
4. **The E-13 classification** against the §3 quotient.
5. **A rubric-adequacy note:** are the labels well defined and mutually exclusive for these rows? Flag any row where the
   label depends on interpretation.

## 6. Grade rule (fixed before the reproduction)

| grade | condition (comparing the blind reader with the original matrix) |
|---|---|
| **V0-4-A** | the headline verdict agrees, and every row label agrees, or differs only by an explicitly justified dual-label / adjacent-label refinement that changes no FIX status |
| **V0-4-C** | the headline verdict agrees (no FIX), but ≥ 1 row is materially reclassified, or an upstream trace / witness is corrected |
| **V0-4-F** | the blind reader establishes a FIX row, so the original verdict fails |
| **V0-4-I** | the reader cannot classify enough rows from allowed sources (precise obstruction) |

The suffix "(orchestrator-exposed)" is appended to any VER-I1 grade.
