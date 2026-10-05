# F0 — CHARTER 02 (REPAIRED)
**Status:** OWNER-DIRECTED REPAIR (`F0_OWNER_RULING_01.md`). Supersedes the prior F0 charter
draft wherever they conflict. R1–R10 from the reviewer handoff applied exactly; R11 added here.
**Date:** 2026-10-05 · **No computation. No candidate-law search. No implementation of
`R_Gamma`, `Cl_Gamma`, `S`, or `U` in this charter.**

## 1. The question (owner, verbatim, prominent)

> Can presently realized influence data constrain and positively generate future physical
> accessibility without secretly assuming response data for inaccessible contexts?

## 2. Hierarchy (F0-A / F0-PHYS / F0-B / F0-C / F0-D)

### F0-A — KINEMATICS
Define `(C, Γ)` plus representation equivalence and composition. `C` = an access structure
(what can be jointly co-instantiated / jointly accessed); `Γ` = the contextual influence
data legitimately available in the current structure.

### F0-PHYS — PHYSICALITY OF ACCESS (gate; may be the hardest gate)
Establish whether `C` has an **observer-independent physical meaning** — e.g. physical
co-instantiability, simultaneous intervention compatibility, or jointly realizable
transformations — rather than merely encoding the measurement menu some apparatus happens
to permit. If only the latter, `C` is operational metadata, not proto-geometry.
**No locality/geometry language precedes this gate.**

### F0-B — ACCESS TRANSITION / ADMISSIBILITY (R11)
Construct the lawful access-transition relation
`R_Gamma(C, C') ∈ {0,1}`, or the correspondence
`A_Gamma(C) = { C' : R_Gamma(C, C') }`:
using only information legitimately contained in the current `(C, Γ)`, which successor
access structures are allowed? F0-B may **constrain or veto**; it need not positively
choose. A `Gamma`-dependent closure operator `Cl_Gamma` is **one special candidate**,
admissible only where monotone access growth (`C ⊆ Cl_Gamma(C)`) and idempotence are
independently justified — neither is assumed here.

### F0-C — POSITIVE ACCESS SELECTION (only if earned)
`C' = S(C, Γ)` or equivalent. **Constraint ≠ selection.**

### F0-D — COUPLED DYNAMICS (only if both access and influence evolve under a specified law)
`(C_{t+1}, Γ_{t+1}) = U(C_t, Γ_t)`. `Γ' ∈ Adm(C')` is not dynamics; an allowed-successor
set is not dynamics; a veto rule is not dynamics.

## 3. Specker baseline repair (frozen distinctions)

### 3.1 Structural Specker principle
Use only at an **explicitly declared compatibility scope**. The standard relevant
principle concerns **sharp/projective measurements** (propositions): pairwise
compatibility/joint measurability of the sharp objects implies joint compatibility.
Do **NOT** assert it for arbitrary POVMs/unsharp measurements: in quantum theory,
**pairwise joint measurability of unsharp measurements does not imply global joint
measurability** (cf. arXiv:1712.01225). Unsharp pairwise-but-not-globally-compatible
quantum examples are **mandatory hostile controls**: Specker's triangle is an especially
good control precisely because arbitrary unsharp quantum measurements exhibit
pairwise-without-triplewise joint measurability. Use it to test whether F0 knows which
physical compatibility notion it is talking about; the sharp/projective sector and the
unsharp operational sector remain **distinct**.

### 3.2 Statistical exclusivity principles
Event-level **consistent exclusivity / compatible orthogonality / local orthogonality**
stays distinct from the full **structural** Specker principle. Almost-quantum correlations
satisfy consistent-exclusivity-type principles while lying beyond the quantum set
(arXiv:1403.4621) — a derivation of one is not automatically a derivation of the other.

### 3.3 Quantum correlation set
Both of the above stay distinct from **exact quantum realizability**. Almost-quantum
correlations are a **mandatory hostile control**: Gonda et al. show a theory yielding the
almost-quantum set cannot satisfy the stronger Specker principle, and they lie outside the
quantum set.

### 3.4 CHSH precision
Any statement that an exclusivity principle yields the Tsirelson CHSH bound must list the
**additional assumptions** of the cited result (Cabello's GPT result is exclusivity **plus
two further assumptions**, arXiv:1406.5656). Never write "exclusivity alone derives
Tsirelson." Price every additional assumption.

### 3.5 Existing derivations as null comparators
A new F0 derivation of Specker-type structure earns credit only if the whole-law linkage
and **information price** differ materially from known routes — e.g. Bacciagaluppi's
derivation of Specker's principle from maximal entanglement + non-maximal measurements +
no-signalling (arXiv:2305.07917) is a **null comparator**, not a precedent to be beaten by
restatement.

## 4. Bell/Specker reconnaissance status
The completed Bell access-enlargement calculation is recorded as
**CALIBRATION / REPRODUCTION OF KNOWN LOCAL-TO-GLOBAL STRUCTURE**. It demonstrates
blocking/admissibility; it does **not** demonstrate positive selection. It lands as an
instance of **F0-B blocking**: given a proposed enlargement and a closure assumption, the
current `Γ` can veto the enlargement without assigning outcomes to the nonexistent joint
context. It does **not** supply `S`.

## 5. Specker as benchmark, not target law
"Do not make 'derive Specker' the definition of F0 success." Sharp-sector Specker behavior
is one benchmark consequence a deeper access-transition law may derive. If Specker closure
is inserted directly into `R_Gamma`, `S`, `Adm`, initial data, or a hidden global response
object, classify the result **SUPPLIED/RELOCATED at that point**.

## 6. Quantum classification firewall
Every finite result labels quantum status as exactly one of:
`EXACT` · `EXPLICIT-QUANTUM-REALIZATION` · `GRAPH-THETA-AT-DECLARED-SCOPE` ·
`NPA-k-OUTER` · `UNRESOLVED`.
**Never upgrade an outer approximation to exact quantum membership.**

## 7. Terminals
`F0-ACCESS-LAW-PASS` · `F0-CONSTRAINT-ONLY` · `F0-LATENT-RELOCATION` ·
`F0-OPERATIONAL-ONLY` · `F0-RESTATED` · `F0-OPEN` · `F0-CLOSURE-ONLY`
(a valid monotone closure/admissibility rule exists but no positive successor-selection
law is derived).

## 8. Hard stop
This charter repair, plus `F0_BASELINE_AUDIT_01.md`, are the **only** artifacts of this
cycle. No candidate-law implementation, no fixed-point search, no Born/geometry/gravity.
**STOP for owner review.**
