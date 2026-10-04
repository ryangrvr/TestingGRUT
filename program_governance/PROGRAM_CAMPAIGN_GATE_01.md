# GRUT PROGRAM CAMPAIGN GATE 01

**Status:** program law, **awaiting owner review / freeze**.
**Home:** `grut-program-governance-1`, branched from `grut-selector-screen-1-frozen @ a1194546c08d81c110a6ad8cc7583fddeedbebd1`. This branch is the authoritative home for program-wide campaign governance.
**Companion:** `program_governance/CONJECTURE_MODE_CHARTER_01.md`, which holds the full Door C rules.
**Supersedes:** the non-authoritative pointer recorded in `verification_0/VER0_LEDGER.md` (branch `grut-independent-verification-0`, owner pointer at `3770ce7`). That pointer named this file as the rule's future home. VER0 is not amended by this file.

---

## 1. Scope

This gate governs the **opening** of every new scientific campaign in the GRUT program. It does not reopen, amend or re-grade any frozen or closed branch, and it creates no physics.

A proposed campaign must enter through exactly one of two legal doors. A proposal that fits neither door is refused: **NO CAMPAIGN**.

## 2. Door D — Derivation / Selector campaigns

> **No new DERIVATION/SELECTOR campaign may open unless it names the frozen residual entry it intends to eliminate or render downstream, and a named mechanism plausibly capable of distinguishing at least two frozen-admissible members.**

A Door D proposal must state, before any calculation:

| item | requirement |
|---|---|
| D-1 Named target | at least one exact entry in the frozen residual ledger (`residual_synthesis/MASTER_RESIDUAL_LEDGER.md` at the frozen base) that the campaign intends to eliminate or render downstream |
| D-2 Named mechanism | a specific mechanism able to distinguish at least two currently admissible members of that entry |
| D-3 S1 plausibility | plausibly passes the selector-screen **S1 information-elimination** test (`selector_screen_1/SCREEN_CHARTER.md`) |
| D-4 S2 plausibility | plausibly passes the selector-screen **S2** test ("selection, not another no-go": it chooses among ≥ 2 members that both survive every frozen admissibility constraint; referred to in owner rulings as witness discrimination) before calculation |

If any of D-1 to D-4 fails: **NO DERIVATION CAMPAIGN.**

The frozen selector-screen terminal (`selector_screen_1/GRUT_SELECTOR_SCREEN_01.md`) governs Door D.

## 3. Door C — Conjecture campaigns

> **No new CONJECTURE campaign may open unless it names the frozen residual entry it instantiates, an admissible GRUT origin (C0.2), an explicit postulate, its information price, a frozen observable projection, a dataset, and a preregistered kill condition whose numerical boundaries are fixed by a committed owner ruling before data access.**

A Door C campaign does **not** need to eliminate its residual entry. Every new law it introduces is labeled **POSTULATED** and can never be described as derived, selected, emergent or uniquely required.

The full requirements are in `CONJECTURE_MODE_CHARTER_01.md`:
- the seven C0 fields, including C0.2 Origin / GRUT rationale;
- the C1 firewall;
- the owner-lock gate (§4A), with its executor / reviewer / owner authority roles;
- the network preflight (§4B);
- the result states (§4C);
- the permanent Prior-Exposure Disclosure (§4D);
- the anti-curve-fitting rules C-F1–C-F10.

The selector-screen terminal governs Door D, not Door C. Door C does not permit unrestricted parameter fitting.

## 4. Separation of authority

- **Program law** (this file and the Conjecture Mode Charter) lives on `grut-program-governance-1`.
- **Experiments under that law** (card ledgers, equivalence registry, CONJECTURE-GENERAL register, card specs, owner threshold rulings, data provenance, code and results) live on a **separate descendant branch** (e.g. `grut-conjecture-mode-1`). That branch is created only after this governance state is reviewed and frozen, and only by separate authorization.
- An experiment branch may not amend program law. A change to program law requires a new versioned governance file and an owner ruling.

## 5. Standing constraints preserved

- No canonical GRUT modification.
- No modification of frozen or closed branches.
- No PR or merge unless separately authorized.
- Labels: "INTERNALLY PROVED / NOT EXTERNALLY REVIEWED"; numerics are "independent code path, not independent reviewer".
- No novelty claims.

## 6. State at commit

- Card #1 is **NOT EXECUTABLE**:
  - its exact response law is an empty C0 field;
  - thresholds S1–S8 are empty;
  - no DESI product has been accessed;
  - no network preflight has been run.
- The execution branch is not created.
- VER0 remains paused, with V0-5 / V0-6 outstanding and access-dependent.
- The Candidate #50 screen is not started.
- DESI DR3 / the complete-five-year release remains passive Lane 1.
