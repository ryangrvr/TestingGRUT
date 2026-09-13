# P1A_EDGE_ADJUDICATION_01 — PROTOCOL (frozen before execution)

**Commissioned:** owner directive 2026-09-12 (post-`dc110ba`): *"P1A edge adjudication →
graph recomputation → fixed-point determination → certification status update."* The
current certification `dc110ba` is **frozen as the historical audited state** and is not
amended by this operation (owner: "I would not amend dc110ba immediately").

## The question (owner's wording, the ONLY question this adjudication answers)

> Is the `background_time_translation_flow → rung1_inin_formalism / rung2_kms_gate`
> dependency edge genuinely present, and is its semantic representation correct?

Two sub-questions, answered separately and per node:
- **PRESENCE:** does the node genuinely presuppose the background time-translation flow —
  is the dependency real, load-bearing for the node's banked content?
- **REPRESENTATION:** is a `depends_on` edge — which the resident's tier-contradiction
  rule reads as "a result cannot cleanly stand on an open/assumed input (at most
  `derived-pending` on it)" — the semantically correct encoding of that presupposition,
  given that the edges were wired under OMISSION_STANDARD R5 ("wiring is part of banking
  an omission") whose stated purpose is computable multiplicity?

## Pre-registered outcomes and their pre-authorized consequences (owner's words)

- **EDGE_CORRECT** — "the eight-node demotion was legitimately reachable under the
  represented dependency graph; proceed to re-certification."
- **EDGE_INCORRECT** — "correct the graph, recompute the closure from scratch, and
  determine what classifications actually follow."
- **EDGE_SEMANTICALLY_AMBIGUOUS / UNRESOLVED** — "preserve the ambiguity and do not
  manufacture a final certification from it."

A mixed verdict (e.g. one node's edge correct, the other's not) routes each node down its
own branch; the certification update then reflects the mixed state honestly.

## Anti-optimization commitments (owner: "do not optimize the ruling for green, for fewer
demotions, or for preservation of |N| = 0")

Both directional bias risks are named in advance:
- **EDGE_INCORRECT is the framework-favorable direction** — it would re-promote nodes
  toward `shown` (the recorded bias pattern: framework-favorable findings, four collapsed
  in one session).
- **EDGE_CORRECT is the tidiness-favorable direction** — it discharges the certificate's
  conditional and closes the question cleanly (the comfortable-paperwork pattern).

Neither outcome is preferred. The verdict standard is: what the register's own booked law
(CHARTER tier vocabulary; OMISSION_STANDARD_v2 R1–R5; the resident's booked rule; the
register's own edge-typing precedents) plus the nodes' banked physics content forces.
The suite's color after propagation is NOT an input to the verdict.

## Method (two-audit family, adversarial-by-construction)

Three independent ADVOCATES, each assigned one verdict and charged to build the strongest
repository-grounded case for it: (1) EDGE_CORRECT, (2) EDGE_INCORRECT, (3)
AMBIGUOUS/UNRESOLVED. A provenance VERIFIER independently checks every factual claim the
advocates rely on. A mandatory REFEREE rules on evidence strength under the standard
above, answering presence and representation separately for each of the two nodes, and is
forbidden to split differences: a ruling with justification, or AMBIGUOUS with the
residual freedom named. The builder (this session) executes but does not vote.

## Evidence inventory (paths; advocates verify and may extend — this list does not fence them)

1. `provenance/claims.json` — the three nodes' full content: the flow node's booking
   tier_note (CHARTER tell run in full; R4 atomisation: only part (a), the background
   flow, is booked; part (b), state invariance, ABSORBED into rung2's KMS content; R1
   membership fence passed — named positions deny the flow and the denial changes what
   the framework says: no single-ω kernel), ledger_note (R5 wiring "because both
   presuppose it"; the tier question flagged there concerns rung3, not rung1/rung2);
   rung1's `laundering_ok`/`stance_justification` ("the LEDGER, not the tier, carries
   the inputs" — the tier grades four-leg literature verification); rung2's statement
   (KMS detailed balance — definitionally relative to a time-translation/modular flow),
   `edge_note` field (typed non-dependency annotations exist).
2. `provenance/OMISSION_STANDARD_v2.txt` R5 — "the nodes that PRESUPPOSE it MUST GAIN
   depends_on EDGES to it"; stated purpose: multiplicity must be computable, not prose.
3. `CHARTER.md` §1 — tier vocabulary: **derived-pending = "derived modulo a named open
   input"**; "Awaiting a computation ⇒ derived-pending, not derived"; §7 resident scope.
4. `provenance/resident.py` — the tier-contradiction rule's own comment: "a result cannot
   cleanly stand on an open/assumed input (at most 'derived-pending' on it)."
5. Edge-typing precedent: commit `d5e9a99` — "rung3 was never an actual dependency →
   edge-note; depends_on = rung1+rung2" (non-dependencies leave depends_on; genuine
   dependencies stay).
6. `PHYSICS_LEDGER/WALL_HELD_FLAGS_23_REVIEW.md` — F2 items 1–2: the live contradiction,
   with the tension named: "the omission was booked precisely to expose the
   presupposition those nodes rest on — arguably correct physics bookkeeping — while the
   resident's tier rule forbids `shown` resting on `assumed`. **Which wins is an owner
   call.**" (This adjudication, commissioned by the owner with pre-registered outcomes,
   is that call being taken.)
7. The demotion's provenance: `git diff d5e9a99..67059a8 -- provenance/claims.json` —
   eight nodes `shown → derived-pending`, tier strings only, no tier_note change, no
   commit-message mention ("R1/R2 complete + format pass"); sole prose rationale is the
   frozen RAI artifact's Phase-1 line. DISCLOSED here; the demotion's paperwork is NOT
   the question — the edge is.
8. `provenance/OPEN_PASSES.txt` — P1A stanza (question wording, MOOT history); P1B
   stanza (the prose-priced-inputs question is SEPARATE and stays open regardless of
   this verdict).
9. Physics sources booked on the nodes: `callen_welton1951`, `kubo1966` (KMS/FDT
   requires stationarity to define detailed balance), `schwinger1961`, `keldysh1964`,
   `feynman_vernon1963` (the influence functional exists on general backgrounds; the
   single-frequency kernel form K_R(ω,k) requires time-translation invariance).

## Out of scope (not decided here)

P1B (whether `shown` is right for prose-priced inputs — the edge-vs-prose asymmetry is
P1B's subject, not a refutation of edge semantics); P5; P6; the |N| = 0 classification
content; any tier change to any node other than what the ruled branch mechanically
forces; Born 2A/2B (unauthorized, unchanged).

## Post-verdict mechanics (pre-committed)

1. Record the verdict in `P1A_EDGE_ADJUDICATION_01.md` + `.json` (separate result files;
   this protocol is not edited after freeze).
2. Update the P1A stanza: append the ruling history (MOOT → CLOSED with the ruling, or
   MOOT preserved with the ambiguity recorded). History is appended, never rewritten.
3. Graph recomputation + fixed-point determination per branch: EDGE_CORRECT → verify the
   propagated register is the fixed point of the booked tier rule (zero
   tier-contradiction flags; no over- or under-propagation); EDGE_INCORRECT → correct
   the representation, recompute the closure from scratch, report what classifications
   follow (whatever they are); AMBIGUOUS → no register change, ambiguity recorded.
4. Certification status update as a NEW artifact citing `dc110ba` (frozen, unamended).
5. Full suite + expected_red re-run; commit; push to `testinggrut`.

*Protocol frozen 2026-09-13, before any advocate ran. The builder's evidence inventory
above was assembled before freezing; the advocates are charged to attack its selection.*
