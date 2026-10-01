# SCOUT-1 D3-SCOUT RESULT — representation / composition selection

**SCOUT-only premise experiment** (owner ruling after ZOOM_OUT_07). GRUT-RAI is untouched, SCOUT-0 is frozen, and
D2/D4 are closed.

**Files:**
- `D3_PREMISE_PRICE_TABLE.md`;
- `d3_ablation.py` / `.log`;
- the record bridge (read-only sweep; quotes in §5).

**Source grade.** Reconstruction theorems are **SECONDARY**: arXiv is blocked (curl and web fetch both refused).
The computable parts are checked ✓.

## 0. The two verdicts

> **PHYSICS VERDICT: C4-SELECTED-FROM-C5 (partial A).** Given an operational probabilistic framework and
> composition-level principles, the **representation layer** is uniquely selected:
> - complex scalars (target item 2);
> - the reversible group PU(d) (item 5);
> - the Tr(ρE) pairing (item 6, given the framework);
> - no superselection (item 9).
>
> The routes apply to a broader parent class (GPTs), and single-axiom ablation reopens named alternatives. **But by
> the frozen anti-circularity rule, the load-bearing axioms are themselves target items at the C5 level:**
> - local tomography = items 8/4;
> - purification = item 7;
> - the framework = items 1/10 and part of 6.
>
> So the selection is **real for C4** and **pays in C5 coin**: the wall moves up, it does not disappear.
>
> **Route C (Barnum–Wilce) is mostly CONTENT-REPACKAGED.** Jordan/HSD already reduces the field question to the
> finite Jordan list (Koecher–Vinberg), and self-duality already yields the Tr pairing.
>
> **GRUT VERDICT: SELECTOR-SUPPLIED / ABSENT.** None of the load-bearing principles is earned by GRUT:
> - composition / tensor factorization: SUPPLIED ("C-9 CRITERION-SMUGGLED");
> - local tomography: ABSENT (the RRP sibling program prices it as "a selection-type input");
> - purification / dilation: SUPPLIED (lift IRREDUCIBLE/SUPPLIED);
> - continuous reversibility: SUPPLIED / CONDITIONAL (the floor generator is dissipative);
> - the operational probabilistic framework: ABSENT;
> - Born: BORROWED.
>
> **No principle is INCOMPATIBLE.** GRUT's wall therefore remains at **C4/C5**.

## 1. Ablation — which single axioms carry the load (D3-3)

| Removed axiom | Alternative that reopens | Check |
|---|---|---|
| purification (A) | classical probability theory | ✓ a mixed bit has no pure joint extension |
| continuous reversibility (B) | classical theory (a discrete group) | standard |
| local tomography (A, B, C) | real QM, fermionic QT, quaternionic, spin factors | ✓ rebit `σ_y⊗σ_y` (W1-L); ✓ fermionic: local-even differences 0.0, global ±1; ✓ only complex matches `(d+1)²` (16 = 16; real 10 ≠ 9; quaternionic 28 ≠ 36) |
| Jordan / HSD (C) | general GPTs (boxworld …) | standard |
| framework affinity | non-Born outcome rules | ✓ `h_ε(p)` not affine (0.7023 vs 0.7500) |

**Load bearers:**
- **local tomography** for the field;
- **purification / continuous reversibility** for quantum vs classical;
- **the framework itself** for the Born pairing.

## 2. Purification hostile (D3-4): **PARTIAL SELECTOR**

- Purification holds in **real and complex** QM: ✓ real purification reproduces ρ, and two purifications are related
  by an orthogonal / unitary map on the purifier. It therefore **does not select the field**.
- It does select quantum-like over classical theories.
- It presupposes composites with entangled pure states. Its content (existence plus uniqueness up to reversible maps
  on the purifier) is essentially the Stinespring/Schmidt structure, which is **target item 7**.

**Classification:** PARTIAL SELECTOR. It is a genuine discriminator against classical/GPT alternatives, but it
repackages the dilation item.

## 3. Local tomography revisited (D3-5)

- **Within Jordan/HSD + a qubit + composites**, local tomography alone singles out complex QM (✓ counting).
- **Without Jordan/HSD**, local tomography needs CDP's or MM's other axioms to exclude non-Jordan GPTs.
- **Price comparison** (partial order, least to most upstream content):
  - **MM / Hardy:** finiteness + subspace equivalence + continuous reversibility + LT + all measurements. Mostly
    C/D/E-class; LT is the only G.
  - **CDP:** framework + five E-class axioms + **purification (G for item 7)** + LT.
  - **Barnum–Wilce:** **Jordan/HSD (G for items 1, 6, most of 2)** + LT + qubit.

  MM/Hardy is the cheapest in target-in-disguise content. Barnum–Wilce is the most expensive.

## 4. Born accounting (D3-6)

The Born pairing is **not** derived independently by any route. It is bought by the **framework's affine probability
premise** (plus self-duality or a Hilbert representation). That is the same premise as W1-G's Busch-type additivity
and SCOUT-0 P-15's decomposition independence.

**The "three primitives collapse at once" hypothesis** (complex field, Born, `p = |α|²`):
- It **survives only in the weak sense** that one package — framework + composition principles — yields all three.
- The Born / affine part is purchased by the framework, not derived.
- **Status: HIGH-VALUE CROSS-LAYER HYPOTHESIS → PARTIALLY CONFIRMED (physics), INAPPLICABLE (GRUT: package not
  earned).**

## 5. GRUT bridge (D3-7) — record evidence

| Principle | Status | Record evidence |
|---|---|---|
| composition / tensor | SUPPLIED | `EA0_OWNER_RULING_01.md:61-63` (verified): "Locality being an earned ingredient does not imply that the tensor-product structure … is uniquely or necessarily selected … C-9 = CRITERION-SMUGGLED" |
| local tomography | ABSENT | none in GRUT. `rrp/RRP_01_CONCEPTUAL_PASS_01.md:89` (verified): "local tomography (the composite rule — a selection-type input) … composite structure priced as axiom" |
| lift / complex structure | SUPPLIED | `L0_LIFT_SELECTION_OWNER_RULING_02.md:9` (verified): "L0 LIFT SELECTION = IRREDUCIBLE/SUPPLIED" |
| operational probability rule | ABSENT / consumed | `rrp/RRP_00_REQUIREMENTS_PICTURE_01.md:76` (verified): "The operational probability rule is consumed, not produced." |
| Born | SUPPLIED | `NO_GO_LEDGER.md:66,70` (agent-reported): "BORROWED" |
| no-signalling, perfect distinguishability, HSD / Jordan, qubit | ABSENT | agent-reported |
| reversibility | SUPPLIED / CONDITIONAL | agent-reported: the floor generator is dissipative and supplied (synthesis A-3); a conservative parent is CONDITIONAL (S5) |
| purification / dilation | SUPPLIED | agent-reported: the Sz.-Nagy dilation exists, but its physical selection is supplied |

**No circular use:** the supplied quantum lift was not used anywhere in this audit.

## 6. Information compression test (D3-8)

| Criterion | Holds? |
|---|---|
| applies to a broader parent class | ✓ (GPTs) |
| axioms independently motivated / operational | ✓ mostly (MM best) |
| **no axiom contains the target** | **✗**: LT (items 8/4), purification (item 7), framework (items 1/10) |
| conclusion unique | ✓ |
| ablation reopens alternatives | ✓ |

**Verdict:** a genuine C4 compression **relative to C5**. It is not an elimination of input information. The
irreducible layer is now:

> **the operational probabilistic framework + composition (tensor / local tomography) + one quantum-vs-classical
> principle (purification or continuous reversibility) + finiteness / dimension.**

## 7. TC-4′ consequence (D3-10)

- **Universal (physics) wall:** moves from C4 to **C5** (composition / operational framework / dimension) + basin.
- **GRUT's wall:** stays at **C4/C5**, because none of the C5 principles is earned. No primitive collapse is
  available to GRUT on the present record, so the **anti-circularity reproduction is not triggered**: the GRUT
  bridge is empty.

**Status: D3 COMPLETE — PHYSICS: C4-SELECTED-FROM-C5 (partial A; Barnum–Wilce route mostly repackaged). GRUT:
SELECTOR-SUPPLIED / ABSENT. Born bought by the framework.**
