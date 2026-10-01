# S2-5 RESULT — what removes the hidden global degrees of freedom that break local tomography?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-5 (pre-registered).

**Files:** `probes/S2-5/s2_5_tomography.py`, with log `probes/S2-5/s2_5_tomography.log`.

**Labels:**
- **LT REDUCED TO A SHARED-REFERENCE STATE (C5-F → C5-H)**;
- STATE/BASIN-PRICED;
- KNOWN RESULT IMPORT: Aleksandrova–Borish–Wootters (real QM + a universal rebit); Bartlett–Rudolph–Spekkens
  (reference frames lift superselection) — SECONDARY, with the core computation ✓.

## 0. Verdict

> **The local-tomography failure of real QM is a missing shared reference frame.** The hidden global parameter is
> a correlation in the "imaginary unit" `J`. Two real states `ρ± = (I ± J⊗J)/4` are indistinguishable by local
> statistics (difference 0.0).
>
> **Give each party a local reference rebit, and LT is restored iff the references are in a shared, correlated
> state:**
>
> | References | ⟨J⊗J⟩_ref | Max local-statistics difference between ρ±, from party-local measurements on A R_A \| B R_B |
> |---|---|---|
> | shared real "singlet" | +1.000 | **2.000** (maximal: distinguishable) |
> | product \|00⟩ | 0.000 | **0.000** |
> | maximally mixed | — | **0.000** |
>
> **Complex QM = real QM + one globally shared reference rebit carrying J.** The ABW-type encoding reproduces a random
> complex expectation value exactly (−3.454168 both ways) with a real symmetric observable.
>
> **Consequence (structural):** local tomography, and with it the ℝ-vs-ℂ choice that SCOUT-1 D3 traced to LT, is
> **not a kinematic primitive**. It is a property of **which reference state the universe contains and shares**: a
> C5-H (basin / state) datum. **C5-F reduces to C5-H.**

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| real single-system kinematics | **inserted** (the weaker starting point) |
| local reference systems for each party | **inserted** (access: who holds which subsystem) |
| the shared correlated reference **state** | **inserted**: a state / basin datum |
| distinguishability of the global parameter by local statistics, given that state | **forced** |
| complex-QM statistics from real QM + global reference | **forced** (encoding identity) |

## 2. Accounting

- **The ℝ→ℂ selection is converted from an axiom (LT) into a state datum:** the universe contains a shared
  J-reference correlated across all parties.
- This is a genuine **move of the boundary**. The field choice joins basin / preparation data (C5-H) instead of the
  operational framework.
- It does **not** eliminate information. "Why is there a universally shared reference?" is now a C5-H question:
  cosmological / initial-state, S2-6.

**Scope:** shown for one hidden parameter on two rebits. Full LT (every global parameter, n parties) needs the
standard ABW construction (SECONDARY). The superselection / fermionic analogue (a parity reference) is cited, not
computed.

**Status: S2-5 COMPLETE — LT failure = a missing shared reference; LT restored by a shared correlated reference state.
C5-F (and with it ℝ vs ℂ) reduces to a C5-H state / basin datum.**
