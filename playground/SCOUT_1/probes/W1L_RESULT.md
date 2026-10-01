> **AUDIT REPAIR 01:** the headline becomes **REAL-vs-COMPLEX DISCRIMINATOR / SELECTION WITHIN THE DECLARED FIELD
> CLASS**.
> - Kept: standard real QM fails local tomography; standard complex QM satisfies it (rebit witness).
> - Not claimed: that local tomography + tensor composition alone reconstruct complex QM from all probabilistic
>   theories. Full reconstruction needs further axioms (e.g. a Jordan / homogeneous-self-dual class, a qubit, and
>   composition / non-signalling assumptions, depending on the theorem).
> - The quaternionic case is **open** (not checked).
> - The owner option "composition + tomography + no-signalling collapses three primitives" is downgraded to a
>   **HIGH-VALUE CROSS-LAYER HYPOTHESIS**. It still has to show that the complex field, the Hilbert/effect
>   structure, `p = |α|²` and the Born weights all follow from one common premise set without circularly supplying
>   any of them.

# SCOUT-1 W1-L RESULT — do reconstruction axioms select the lift's number field? (C32)

**Charter:** `PROBE_CHARTERS.md` §W1-L. Preregistered outcomes: SELECTION PRINCIPLE (premise-priced) / NOT FIXED.

**Files:** `w1l_local_tomography.py`, with log `w1l_local_tomography.log`.

**Labels:**
- SELECTION PRINCIPLE (premise-priced);
- KNOWN RESULT IMPORT:
  - Hardy 2001 / Chiribella–D'Ariano–Perinotti 2011 / Masanes–Müller 2011 — STANDARD-TEXTBOOK ✓ for the real/complex
    counting;
  - quaternionic composites — SECONDARY (not checked);
  - Renou et al. 2021 network experiment excluding real QM — SECONDARY, not readable here;
- NON-DISTINCTIVE.

## 0. Verdict

> **SELECTION PRINCIPLE (premise-priced) — a supplied primitive is EXCHANGED, not removed.**
>
> **What is selected.** Local tomography (`K_AB = K_A·K_B`) together with a tensor-product composition rule
> selects **complex** QM over real QM, so the number field is selected:
> - complex counting is exactly multiplicative: `d²·d² = (d·d)²` for (2,2), (2,3), (3,3), (2,4);
> - real QM has global-only parameters, with deficit `K(d²) − K(d)² = 1, 9, 36, 100` for d = 2…5;
> - explicit witness: `ρ± = (I ± σ_y⊗σ_y)/4` are both valid **real** states. They agree on **every** local real
>   product observable (max difference 0.0) and are separated by the real global observable `σ_y⊗σ_y` (±1).
>
> **The price.** This removes SCOUT-0's "complex structure" lift price (P-08) **only by supplying a composition
> rule and the local-tomography principle**. It is the same pattern as W1-G: a principle that is physically
> natural once composition is supplied.
>
> **The empirical side.** Network tests (Renou et al. 2021, secondary) reportedly exclude real QM under an
> independent-sources assumption. If so, the field is fixed **empirically**, not by GRUT. That is an imported
> datum.

## 1. Ten-point hostile test of "local tomography ⇒ complex field"

| # | Question | Answer |
|---|---|---|
| 1 | Q varies while P holds? | No: real fails (witness) and complex passes (rank 16 = 4², 64 = 4³). |
| 2 | P silently contains Q? | Partly. Multiplicativity of `K = d^r` with r = 2 is close to the complex field. The reconstruction theorems need further axioms (continuous reversibility, purification) to exclude all non-quantum theories. |
| 3 | Representation-dependent? | No. |
| 4 | Physical or gauge? | Physical: the rebit witness is operationally distinguishable only by a global measurement. |
| 5 | Standard? | Yes. |
| 6 | Parent variation? | Several dimension pairs. The quaternionic case is not checked (SECONDARY). |
| 7 | Composition? | It *is* the premise. |
| 8 | Coarse-graining? | Not applicable. |
| 9 | Unique or stationary? | Unique among {ℝ, ℂ, ℍ} under the counting. |
| 10 | Boundary condition selecting? | No. |

## 2. What this means

- **Q1 lift price "complex structure" → "local tomography + tensor composition".** This is an exchange of
  supplied primitives. The new primitive is arguably more operational, and is shared by all physics, so the
  selection is **NON-DISTINCTIVE**.
- **TP-1 instance 7.** Again (supplied composition structure) × (consistency principle) × (imported theorem).
  W1-G and W1-L now form a **composition-principle subfamily**: both selectors become available the moment
  a tensor-product composite is supplied. That is the strongest *shared* premise of the quantum layer.
- **Owner-relevant (logged, not asked):** if the owner ever admits "composition + local tomography + no-signalling"
  as earned, then **three** supplied primitives would collapse to that one premise:
  - the number field (W1-L);
  - the Born weight `h(p)` (W1-G);
  - the `p = |α|²` identification (P-15 chain).

  That would be the largest reduction of supplied primitives seen in either campaign. It is a premise change
  and is recorded in `FRONTIER_QUEUE.md` as an owner option.

**Status: W1-L COMPLETE — SELECTION PRINCIPLE (premise-priced: tensor composition + local tomography). Rebit
witness exhibited. Supplied primitive exchanged, not removed. Composition-principle subfamily (W1-G + W1-L)
identified.**
