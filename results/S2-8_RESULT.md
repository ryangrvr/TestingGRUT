# S2-8 RESULT — quantum Darwinism: does redundancy select the pointer AND the split?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-8. Pre-registered firewall: redundancy defined only after the split is
specified → **A-PRICED**.

**Files:** `probes/S2-8/s2_8_darwinism.py`, with log `probes/S2-8/s2_8_darwinism.log`.

**Wall attacked:** A.

**Model.** A system qubit S plus 8 environment qubits, with H = Σ g_k Z_S Y_k. S starts in |+⟩ and E in |0…0⟩.

**Labels:**
- **SELECTOR (pointer, given split)**;
- **SPLIT NOT SELECTED → A-PRICED** (and D-priced via S2-1b);
- KNOWN RESULT IMPORT: Zurek's quantum Darwinism; Page-curve behaviour of scrambled states — SECONDARY, with
  numerics ✓.

## 0. Verdict

**1. Given the split S | E₁ … E₈, redundancy works.**
- I(S:F) by fragment size is 0.92, 0.99, 1.00 … with a plateau at H(S) = 1, rising to 2 at the full environment.
- R₀.₁ = 8.

**2. Pointer selection given the split.** Holevo χ(basis : E₁):

| basis | χ(S:E₁) | χ(S:E₁E₂E₃) |
|---|---|---|
| Z | 0.871 | 0.9998 |
| X | 0.000 | 0.000 |
| (Z+X)/√2 | 0.350 | 0.399 |

The pointer is Z, the observable the interaction couples to. **SELECTOR (pointer, given split)**, priced by the
interaction H, which is D.

**3. Can redundancy pick which slot is the system?** No. With the same state and the same qubit TPS, R₀.₁ = 8.0 for
*every* choice of system slot j = 0…8, for both partial and perfect records. The record state is GHZ-like:
- each qubit is redundantly recorded by all the others;
- redundancy is symmetric under S ↔ E_k;
- **the split is not selected**.

**4. Redundancy across TPSs of the same global state.** The same GHZ state gives:
- frame 1: H(S) = 1, R = 8;
- frame 2 (CNOT-ladder frame): H(S) = 0, I = 0, R = 0.

Any two pure states of equal dimension are unitarily related. Maximizing redundancy over all TPSs therefore cannot
favour a state. It can only favour a TPS, and the TPS must come from elsewhere: H-locality (S2-1b, D) or access (A).

**5. Fragment individuation.**

| grouping | fragments | R₀.₁ |
|---|---|---|
| 1-qubit fragments | 8 | **8** |
| 2-qubit fragments | 4 | **4** |
| 4-qubit fragments | 2 | **2** |

R counts fragments, so its value depends on how the observer's accessible fragments are individuated. → **A**.

**6. Environment self-interaction.**
- With all-to-all scrambling in E, during or after record formation, I(S:F) by size becomes 0.06, 0.18, 0.42, 0.98,
  1.58 … The plateau disappears.
- The residual "R = 2" is the Page-curve floor (half the environment knows everything). It is not a record.
- Darwinism therefore needs a **non-scrambling environment of separately accessible fragments**: D (no env-env
  interaction at the record time scale) + A (fragment access).

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| pointer observable, given split and H | **forced** |
| which slot is the system | **not selected** (symmetric) |
| the TPS in which records are redundant | **not selected by redundancy**; frame-relative (GHZ ↔ product) |
| fragment individuation (the value of R) | **inserted: A** |
| a non-scrambling environment | **inserted: D** |

## 2. Accounting against T2-1′

**A is not eliminated by Darwinism.** Redundancy is defined only after:
- the split is given;
- the fragments are given.

Both are A-items, the split possibly D via S2-1b. Darwinism converts "which observable is classical" (a C5 readout
item) into "which observable the interaction couples to" (D). That conversion is real **only after A supplies the
split and the fragmentation**.

**Firewall outcome: A-PRICED.** The observer/access loop closes. The observer, as a reader of fragments, is
presupposed by the very redundancy that was meant to define objectivity for observers.

**Status: S2-8 COMPLETE — SELECTOR (pointer, given split); split not selected (slot symmetry; TPS-relative: R = 8 vs 0
for the same state); R depends on fragment individuation (8 / 4 / 2); needs a non-scrambling environment. A-PRICED +
D-priced.**
