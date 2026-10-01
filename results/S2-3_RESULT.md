> **REPAIR 01 (owner audit):** "no measure needed" means **no independently chosen invariant measure**. Unique
> ergodicity *supplies* a unique invariant measure and uniform time averages for the relevant observable class. It
> does not eliminate the supplied topology/dynamics or the coarse partition (access, A).

# S2-3 RESULT — probability from deterministic microdynamics + coarse access

**Charter:** `probes/PROBE_CHARTERS.md` §S2-3 (pre-registered firewall: an invariant measure existing ≠ a probability
rule).

**Files:** `probes/S2-3/s2_3_probability.py`, with log `probes/S2-3/s2_3_probability.log`.

**Labels:**
- **MEASURE-PRICED** (mixing route);
- **TRUE DERIVATION, scoped** (uniquely-ergodic time averages);
- **DEFINITIONAL** (affinity);
- KNOWN RESULT IMPORT: Birkhoff; mixing; unique ergodicity (Weyl equidistribution); the residuality of non-normal
  numbers (SECONDARY) — numerics ✓.

## 0. Verdict

> **Three separate notions of "probability", each paying a different price.**
>
> **1. Ensemble relaxation (mixing chaos: the doubling map).** Absolutely continuous preparations forget their shape:
> - density 2x: frequency of A goes 0.248 → 0.375 → 0.439 → 0.484 → 0.500 by t = 8.
>
> Singular *invariant* preparations never move:
> - Bernoulli(0.3) stays at 0.70;
> - the fixed point stays at 1.00.
>
> **The limiting statistics are therefore fixed only once the admissible preparations are restricted to measures
> absolutely continuous w.r.t. Lebesgue.** → **MEASURE-PRICED.** The SRB / physical measure is "earned" only relative
> to a supplied reference measure on preparations.
>
> **2. Single-trajectory frequencies.**
> - **Mixing / positive-entropy case:** a Lebesgue-random point gives frequency ≈ 1/2 (0.547 at n = 190, still
>   converging). A constructed point with growing digit blocks oscillates forever (0.333 ↔ 0.67–0.71), with no limit.
>   Points without limiting frequencies form a topologically generic (residual) set (SECONDARY). **"Generic" is
>   ambiguous: Lebesgue-typical and Baire-typical points disagree.** The choice is supplied.
> - **Uniquely ergodic case (irrational rotation):** time averages converge to 0.5000 from **every** initial point
>   tested (0, 0.123456, 0.5, 1/3). No measure is needed. **This is a genuine measure-free frequency rule:**
>   **TRUE DERIVATION, scoped to uniquely ergodic dynamics.**
>   - **But ensembles never relax.** From an a.c. preparation on [0, 0.1), the frequency jumps 1, 0, 1, 1, 1, 0.017, 0,
>     1 forever. Preparation information is never forgotten.
>
> **3. Affinity.** The mixture's frequency equals the mixture of frequencies: 0.6401 vs 0.6401; 0.6417 vs 0.6414;
> 0.6404 vs 0.6415 (sampling noise). This is automatic because preparations are measures and dynamics acts by linear
> pushforward. **DEFINITIONAL:** affinity is inherited from treating preparations as measures, not derived from the
> dynamics.

## 1. Structural finding — candidate NO-GO (scoped), with its hostile

**Tension:**
- *forgetting preparations* (mixing) and *earning the measure* (unique ergodicity) pull in opposite directions in
  the tested systems;
- hyperbolic / positive-entropy chaos has a dense set of periodic orbits, hence infinitely many invariant measures,
  so it is **never uniquely ergodic**;
- the earned-measure system (rotation) does not mix.

**Candidate NO-GO (positive-entropy scope):** a positive-entropy (hyperbolic) system cannot both forget
preparations and earn its measure. The selection of the physical measure requires a supplied reference class of
preparations.

**Hostile (counterexample to any unscoped version):** horocycle flows on compact hyperbolic surfaces are **uniquely
ergodic and mixing** (Furstenberg; Marcus — SECONDARY). They have **zero entropy**. So "earned measure + forgetting"
is possible in parabolic, zero-entropy dynamics. → Spawn **S2-3b**: does a zero-entropy, uniquely ergodic, mixing
primitive dynamics deliver a measure-free, preparation-forgetting probability rule, and what does it cost (e.g. no
chaos, special geometry)?

## 2. Accounting

| Notion | Earned? | Price |
|---|---|---|
| ensemble limit (mixing) | no | a.c. preparation class (Lebesgue) — supplied |
| time-average frequency, positive entropy | no | "Lebesgue-almost every" vs "Baire-generic" — supplied |
| time-average frequency, uniquely ergodic | **yes** (all points) | uniquely ergodic dynamics; coarse partition (access); no forgetting |
| affine rule | definitional | preparations-as-measures |

**Access:** the coarse partition is supplied in every case.

**Observer loop:** the "frequency" notion presupposes a record of repeated coarse readings — a supplied readout.

**Status: S2-3 COMPLETE — probability is earned measure-free only for uniquely ergodic dynamics (which does not
forget preparations). Mixing chaos forgets preparations but needs a supplied reference measure. Affinity is
definitional. S2-3b spawned (zero-entropy uniquely-ergodic mixing).**
