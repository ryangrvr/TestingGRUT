# D3 PREMISE PRICE TABLE (SCOUT-only)

## How to read this table

**Source grade.** Theorem statements are SECONDARY: they come from training knowledge, because arXiv / primary texts
are blocked from this sandbox. Computational checks are in `probes/d3_ablation.py` / `.log`. The GRUT column comes
from a read-only record sweep, with key quotes verified (`probes/D3_RESULT.md` §5).

**Class codes:**
- A — operational / empirical;
- B — composition principle;
- C — symmetry;
- D — regularity / continuity;
- E — information principle;
- F — mathematical structural assumption;
- G — target-in-disguise, relative to the **frozen D3-0 target items 1–10**.

**GRUT status codes:**
- EARNED / CONDITIONAL / SUPPLIED / ABSENT / INCOMPATIBLE / UNFORMULABLE.
- The record has **no INCOMPATIBLE entries**. The floor is a classical deterministic substrate, and the lift is
  IRREDUCIBLE/SUPPLIED.

## Framework premise (shared by routes A and B)

| Premise | Class | Content | Excludes | Supplies target item | GRUT |
|---|---|---|---|---|---|
| Operational probabilistic framework: convex states, effects with an affine probability pairing, sequential/parallel composition | F + **G** | outcome probabilities are affine in the state; systems compose | every non-affine outcome rule. This includes the W1-G non-Born family: `h_ε` is not affine (0.7023 vs 0.7500, `d3_ablation.log` §4) | **1, 10, part of 6, part of 3** | **ABSENT** (GRUT's "convex" means a convex potential). RRP: "reconstructions consume the operational framework itself as primitive" |

## Route A — Chiribella–D'Ariano–Perinotti (SECONDARY)

| Axiom | Class | Content | Removing it admits (ablation) | Supplies target item | GRUT |
|---|---|---|---|---|---|
| Causality | A/B | no signalling from future choices; a unique deterministic effect | retro-causal GPTs | — | no-signalling **ABSENT**; universal cone **SUPPLIED** |
| Perfect distinguishability | E | every non-completely-mixed state is perfectly distinguishable from some state | theories with "fuzzy" boundaries (e.g. some restricted GPTs) | part of 1 | **ABSENT** |
| Ideal compression | E | lossless encoding into a minimal faithful subsystem | — | — | **ABSENT** |
| Local distinguishability (= local tomography) | B + **G** | composites determined by local statistics | **real QM, fermionic QT** (✓ fermionic: identical local-even statistics, global ±1) | **8, 4** | **ABSENT** (RRP prices it as "a selection-type input") |
| Pure conditioning | E | pure ⊗ pure-effect conditioning gives a pure state | — | — | **ABSENT** |
| **Purification** | B + **G** | every state has a purification, unique up to reversible maps on the purifier | **classical probability theory** (✓ no pure extension of a mixed bit); boxworld-type GPTs | **7** (literally the dilation item) | **SUPPLIED** (Sz.-Nagy dilation exists, but its physical selection is supplied; L0_LIFT_SELECTION_02) |

## Route B — Masanes–Müller / Hardy (SECONDARY)

| Axiom | Class | Content | Removing it admits | Supplies | GRUT |
|---|---|---|---|---|---|
| Finiteness / simplicity (`K = N^r`) | F/E | finite-dimensional state spaces; minimal r | infinite-dimensional or odd scalings | part of 1 | **ABSENT** |
| Local tomography / composite rule (`K_AB = K_A K_B`, `N_AB = N_A N_B`) | B + **G** | as above | real (10 ≠ 9) and quaternionic (28 ≠ 36) composites (✓ counting) | **4, 8** | **ABSENT** |
| Equivalence of subspaces | C/E | all two-level subsystems are equivalent | mixed or non-homogeneous theories | part of 5 | **ABSENT** |
| **Continuous reversibility** | C/D | continuous reversible transitivity on pure states | **classical theory** (a discrete permutation group) | **5** | floor is dissipative, **SUPPLIED**; a conservative parent is **CONDITIONAL** (S5) |
| All measurements allowed | A/F | no restriction on effects | restricted-effect GPTs | part of 6 | **ABSENT** |

## Route C — Barnum–Wilce (Jordan / HSD) (SECONDARY)

| Premise | Class | Content | Removing it admits | Supplies | GRUT |
|---|---|---|---|---|---|
| **Jordan / homogeneous self-dual state cones** | F + **G** | by Koecher–Vinberg, the state spaces are exactly the Euclidean Jordan algebras: real / complex / quaternionic / spin factors / octonionic 3×3 | general GPTs (boxworld …) | **1, 6 (self-duality ⇒ Tr pairing), most of 2** (it reduces the field question to a finite list) | **ABSENT** |
| Local tomography | B + **G** | as above | real, quaternionic, spin-factor composites (✓ only complex matches) | **8** | **ABSENT** |
| Existence of a qubit / suitable composite | A/F | a two-level system composes | classical-only | part of 3 | **ABSENT** |

## Load-bearing summary (ablation)

| Selection | Single load-bearing axiom | Evidence |
|---|---|---|
| quantum vs classical | **purification** (A) / **continuous reversibility** (B) | classical fails both (✓) |
| complex vs real / quaternionic / fermionic | **local tomography**, in all three routes | ✓ rebit (W1-L), ✓ fermionic, ✓ ball counting |
| Born / Tr pairing | the **framework's affine probability rule** (+ self-duality) | ✓ `h_ε` non-affine |
| no superselection | local tomography (fermionic parity fails it) | ✓ |
