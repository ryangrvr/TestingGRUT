> **REPAIR 01 (owner audit):** "LT iff ordered phase" is withdrawn. In the tested Ising reference-field model:
> - asymptotically shared reference coherence exists in the ordered phase and disappears at long distance in the
>   disordered phase;
> - this activates the S2-5 J-reference witness.
>
> Finite-distance correlations are not a binary theorem about LT: above T_c they decay but do not vanish at finite
> distance. At T = T_c correlations decay algebraically, and there is no non-zero long-range order parameter in the
> thermodynamic limit.

# S2-6 RESULT — can dynamics select basin / preparation / measure data?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-6 (pre-registered; the measure problem is central).

**Files:** `probes/S2-6/s2_6_basin.py`, with log `probes/S2-6/s2_6_basin.log`. A first run was killed by a worker
restart and rerun with smaller arrays, with no change to the method.

**Labels:**
- **BASIN SELECTED ONLY MEASURE-PRICED** (attractors, MaxEnt, physical measures);
- **SHARED-FRAME EXISTENCE = ORDERED-PHASE BASIN** (links S2-5 → C5-H);
- KNOWN RESULT IMPORT: SRB measures; Jaynes' MaxEnt coordinate problem; 2D Ising long-range order — STANDARD ✓.

## 0. Verdict

> **1. Coexisting attractors.** The dynamics `x' = x − x³ + 0.3` fixes the attractors, but not their "probabilities".
> P(basin of +1.13) is 0.587 under Lebesgue on [−2, 2] and 0.603 under an asymmetric reference measure. Basin weights
> are reference-measure volumes. **MEASURE-PRICED.**
>
> **2. Maximum entropy.** Unconstrained MaxEnt is uniform *in a chosen coordinate*:
> - uniform in x gives P(x < 0.5) = 0.499;
> - uniform in y = x² gives 0.249.
>
> MaxEnt needs a reference measure. **MEASURE-PRICED.**
>
> **3. Physical measure (logistic r = 4).** Two different a.c. preparations converge to the arcsine statistics
> (0.5001, 0.5002). The period-2 Dirac preparation stays at 1.000. Forgetting holds only within the a.c. class
> (the same as S2-3). **MEASURE-PRICED.**
>
> **4. Shared reference frames from long-range order** (the follow-up to S2-5). The parties' frames are taken as spins
> of a 2D Ising "frame field", and the restored local-tomography signal `2⟨f_A f_B⟩` is measured at distance 32:
>
> | T | 1.8 | 2.0 | 2.2 | 2.6 | 3.2 |
> |---|---|---|---|---|---|
> | `2⟨f_A f_B⟩` | 1.83 | 1.66 | 1.24 | −0.02 | 0.00 |
>
> T_c = 2.269. **A shared frame, and hence LT, exists iff the frame field is in its ordered basin.**
> - The *sign* of the order is irrelevant: J → −J is complex conjugation.
> - Its *existence* is a basin datum: T < T_c, plus d ≥ 2 (1D has no order at T > 0).

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| attractors, ordered vs disordered phases | **forced** by the dynamics |
| which basin / phase the universe is in | **inserted** (initial condition, temperature) |
| probabilities over basins | **inserted** (a reference measure) |
| the vacuum sign / frame orientation | **gauge** (irrelevant) |
| the frame field's dimension (d ≥ 2 for order) | **inserted** (C5-G) |

## 2. Consequence for the campaign

**C5-H is not selected by any mechanism tested.** Every basin / measure selection either:
- uses a supplied reference measure (attractor volumes, MaxEnt, SRB); or
- is earned only in non-forgetting uniquely ergodic dynamics (S2-3).

**S2-5 + S2-6 together:** local tomography, hence ℂ over ℝ, is equivalent to the universe being in the
**ordered phase of a shared reference field**. That is a basin datum.

The only freedom removed is the frame orientation, which is gauge.

**Status: S2-6 COMPLETE — basin / measure data not dynamically selected without a supplied measure. LT ⇔ ordered phase
of a shared frame field (basin datum; orientation gauge).**
