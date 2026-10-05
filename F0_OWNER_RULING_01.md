# F0 — OWNER RULING 01 (F0 DIRECTION — PASS · SPECKER AMENDMENT — REPAIR BEFORE FREEZE)

**Date:** 2026-10-05 · **Given from the owner's words** (ruling message, with the attached
reviewer ruling adopted verbatim). Any misstatement is corrected by owner edit, not defended.
**Precedent chain:** EA-0 (`EA0_OWNER_RULING_02.md`), the capability-contract rule
(`rrp/CAPABILITY_CONTRACT_01.md`), the hard stop (`STATE.md` boundary block).

## 1. Ruling

- **F0 DIRECTION — PASS.**
- **SPECKER AMENDMENT — REPAIR BEFORE FREEZE.**
- **Do NOT run the finite candidate-law search yet.**
- Revise the charter only; perform the targeted baseline audit; commit the repaired
  charter/baseline artifacts; **STOP for owner review.**
- Apply **R1–R10 from the reviewer handoff exactly**, plus the additional repair below.

## 2. R11 — ACCESS TRANSITION BEFORE CLOSURE

Do **not** make a `Gamma`-dependent closure operator the fundamental F0-B object. A closure
operator already assumes unearned structure: monotone access growth
(`C ⊆ Cl_Gamma(C)`) and saturation/idempotence (`Cl_Gamma(Cl_Gamma(C)) = Cl_Gamma(C)`).
Neither has been derived. A genuinely dynamical access structure may gain contexts, lose
contexts, exchange a compatibility relation for another, cycle, bifurcate, or depend on
history.

The more primitive F0-B object is a **lawful access-transition relation**:

- `R_Gamma(C, C') ∈ {0,1}`, or equivalently `A_Gamma(C) = { C' : R_Gamma(C, C') }` —
  given only legitimately available information in `(C, Gamma)`, which next access
  structures are allowed?

A closure operator `Cl_Gamma` may then be studied as **one special candidate** for
campaigns where monotone access enlargement is independently justified.

## 3. The revised hierarchy

- **F0-A — KINEMATICS.** `(C, Γ)` plus representation equivalence and composition.
- **F0-PHYS — PHYSICALITY OF ACCESS.** Establish whether `C` has an observer-independent
  physical meaning rather than merely encoding an experimenter's available measurement
  menu. No proto-geometry/locality claim before this gate.
- **F0-B — ACCESS TRANSITION / ADMISSIBILITY.** Construct `R_Gamma(C, C')` (or an
  equivalent rule identifying allowed/forbidden changes of access from currently
  legitimate information). This level may constrain or veto proposed changes; it need not
  positively choose one.
- **F0-C — POSITIVE ACCESS SELECTION.** Only if earned: `C' = S(C, Γ)`. Constraint is not
  selection.
- **F0-D — COUPLED DYNAMICS.** Only if both access and influence evolve under a specified
  law: `(C_{t+1}, Γ_{t+1}) = U(C_t, Γ_t)`. `Γ' ∈ Adm(C')` is not dynamics; an
  allowed-successor set is not dynamics; a veto rule is not dynamics.

## 4. Primary F0 question (prominent)

> **Can presently realized influence data constrain and positively generate future physical
> accessibility without secretly assuming response data for inaccessible contexts?**

## 5. Terminals

- `F0-ACCESS-LAW-PASS`
- `F0-CONSTRAINT-ONLY`
- `F0-LATENT-RELOCATION`
- `F0-OPERATIONAL-ONLY`
- `F0-RESTATED`
- `F0-OPEN`
- `F0-CLOSURE-ONLY` (added) — a valid monotone closure/admissibility rule exists but no
  positive successor-selection law is derived.

## 6. Bell/Specker reconnaissance status

The completed Bell access-enlargement calculation is recorded **only** as
**CALIBRATION / REPRODUCTION OF KNOWN LOCAL-TO-GLOBAL STRUCTURE** (`F0_BASELINE_AUDIT_01.md`).
Its earned statement:

> Given a proposed access enlargement and an independently specified closure rule, currently
> available contextual influence data can sometimes rule out the enlarged access structure
> without assigning a probability law to the nonexistent joint context.

It demonstrates **blocking/admissibility**. It does **NOT** demonstrate positive access
selection.

## 7. Hard stop

After the charter repair and the targeted literature baseline
(`F0_CHARTER_02_REPAIRED.md`, `F0_BASELINE_AUDIT_01.md`): commit and push **only** those
governance/baseline artifacts. Do **NOT** implement candidate `R_Gamma`, `Cl_Gamma`, `S`,
or `U`. Do **NOT** run fixed-point searches. Do **NOT** derive Born, geometry, or gravity.
**STOP for owner review.**
