> **AUDIT REPAIR 01 (cross-reference):** the owner option recorded in W1-L (three-primitive collapse) is
> downgraded to a HIGH-VALUE CROSS-LAYER HYPOTHESIS. W1-G's own premise-priced selection is unchanged.

# SCOUT-1 W1-G RESULT — quantum-foundations selectors for the outcome weight law (C10, C11)

**Charter:** `PROBE_CHARTERS.md` §W1-G. Question: do Gleason (d ≥ 3), Busch (POVMs, d ≥ 2) or composition +
no-signalling force `h(p) = p`?

**Files:** `w1g_born_selectors.py`, with log `w1g_born_selectors.log`.

**Test family:** `h_ε(p) = p + ε·sin(2πp)/2π`. It is monotone for |ε| ≤ 1 and satisfies `h(1−p) = 1 − h(p)`.

**Labels:**
- SELECTION PRINCIPLE (premise-priced);
- KNOWN RESULT IMPORT: Gleason 1957; Busch 2003; Gisin/Polchinski non-linearity ⇒ signalling — STANDARD-TEXTBOOK ✓;
- d = 2 projective counterexample family (REDISCOVERED-KNOWN);
- NON-DISTINCTIVE.

## 0. Verdict

> **SELECTION PRINCIPLE (premise-priced): `h(p) = p` is forced uniquely by each of three routes:**
> - noncontextual frame functions on a supplied Hilbert space with d ≥ 3 (Gleason);
> - additive probabilities on the full supplied effect algebra, already at d = 2 (Busch);
> - supplied tensor-product composition + no-signalling (applied to non-maximally entangled states).
>
> **Without those premises it fails.** In d = 2 with projective measurements only, **every** `h_ε` is a valid
> normalized frame function (deviation ≤ 1.6·10⁻¹⁵ over 2000 random bases).
>
> **Premise prices:**
> - the Hilbert space (supplied lift / quantum class, A-6);
> - noncontextuality, or admission of all POVMs;
> - the tensor-product composite.
>
> None is earned in the record (Born BORROWED, NO_GO 7).

## 1. Results

| Test | ε = 0 | ε = 0.05 | ε = 0.2 | ε = 0.5 |
|---|---|---|---|---|
| d = 2 projective, max `|Σh − 1|` | 1.1e−15 | 1.6e−15 | 1.6e−15 | 1.3e−15 |
| d = 3 projective (Gleason), max / median | 1.6e−15 | 2.1e−2 / 6.3e−3 | 8.3e−2 / 2.5e−2 | 2.1e−1 / 5.9e−2 |
| d = 2 trine POVM (Busch) | 6.7e−16 | — | 8.2e−2 | — |
| d = 2 five-star POVM | 6.7e−16 | — | 1.5e−1 | — |
| d = 2 tetrahedral POVM | 6.7e−16 | — | 8.3e−2 | — |
| no-signalling, Bell pair: spread of Bob's marginal over Alice's bases | 1e−16 | 2e−16 | 1e−16 | 0 |
| no-signalling, `cos 0.4 |00⟩ + sin 0.4 |11⟩` | 0 | **2.9e−3** | **1.25e−2** | **3.45e−2** |

**Functional-equation route (d = 3).** `h(x) + h(y) + h(1−x−y) = 1` on the simplex, with `h(0) = 0` from the
d = 2 embedding, gives `h(x+y) = h(x) + h(y)`, so a monotone h is `h(p) = p`. The residuals for `h_ε` are
non-zero, about ε × 0.2.

**Nuance found.** The maximally entangled Bell pair **does not** expose non-linearity: its symmetric overlap
pattern keeps Bob's marginal at ½. Signalling appears only with non-maximal entanglement. A no-signalling test
restricted to maximally entangled resources would **not** select Born.

## 2. Ten-point hostile test of "Gleason/Busch/no-signalling ⇒ `h(p) = p`"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | No, under each full premise set. Yes in d = 2 projective-only. |
| 2 | P silently contains Q? | **Partly.** Busch's additivity over effects is linearity in E, so `v(E) = Tr ρE` is close to a restatement. Gleason's content is real (d ≥ 3 needed). |
| 3 | Representation-dependent? | It needs the supplied Hilbert-space lattice of projectors/effects. |
| 4 | Physical or gauge? | Physical (frequencies). |
| 5 | Standard? | Yes. |
| 6 | Parent variation? | d = 2, 3; three POVMs; two entangled states. Dimension and entanglement matter as stated. |
| 7 | Composition? | Composition is the third route's premise. Under composition a nonlinear h signals (`cos/sin 0.4`). |
| 8 | Coarse-graining? | Coarse-graining outcomes (merging projectors) is additivity, which nonlinear h violates (d = 3 deviations). |
| 9 | Unique or stationary? | Unique. |
| 10 | Boundary condition selecting? | No. But **noncontextuality is the hidden selector**: contextual assignments evade Gleason (Kochen–Specker), and Bohmian mechanics gets Born only from the quantum-equilibrium **initial distribution** (a supplied preparation). |

## 3. What this means

- **Q5 (`h(p)`) can be fixed**, but only by premises about **supplied structure**: Hilbert space,
  measurement/effect algebra, composition. This matches SCOUT-0 P-15 (Born ⟺ martingale ⟺
  decomposition-independence) from the other side: decomposition-independence is the outcome-level
  noncontextuality/additivity premise.
- **TP-1 instance 6.** Selection = (supplied Hilbert/composition structure) × (noncontextuality or no-signalling
  principle) × (imported theorem). Unlike W1-A, W1-S and W1-R, **the principle here is not a symmetry**. It is
  a consistency condition (no-signalling) that is physically compelling **once composition is supplied**.
- **Nuance for the record:** the maximally-entangled case is blind to non-linearity. Any GRUT test of the
  outcome law via composite systems must use non-maximal entanglement.
- **Wave-2 spawn (W2-QE): Born as a dynamical attractor.** Valentini relaxation in de Broglie–Bohm would be
  a *fixed-point/attractor* selector (the W1-C corollary): it drives `h → p` dynamically from a supplied
  non-equilibrium start. It would test whether a dynamical principle can replace the supplied preparation.

**Status: W1-G COMPLETE — SELECTION PRINCIPLE (premise-priced: Hilbert space, d ≥ 3 or POVMs, noncontextuality /
composition). d = 2 projective non-Born family exhibited. Maximal entanglement is blind to non-linearity.**
