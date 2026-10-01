# S2-4 RESULT — convexity / mixtures without supplied randomness

**Charter:** `probes/PROBE_CHARTERS.md` §S2-4 (pre-registered).

**Files:** `probes/S2-4/s2_4_convexity.py`, with log `probes/S2-4/s2_4_convexity.log`.

**Labels:**
- **CONVEXITY DERIVED (scoped), ACCESS-PRICED + INDEPENDENCE-PRICED**;
- KNOWN RESULT IMPORT: Weyl joint equidistribution; Sturmian sequences — STANDARD ✓.

## 0. Verdict

> **Convex mixtures emerge from fully deterministic ingredients, with no supplied randomness.**
> - Ingredients: a hidden uniquely ergodic coin `c_n = [nφ mod 1 < λ]`; deterministic preparations; deterministic
>   outcomes.
> - The coin's mixing weight is **earned** (S2-3: no measure needed): frequency 0.300001 for λ = 0.3, and 0.309017
>   for the irrational λ = φ/2.
>
> **The mixture is observed iff the tester's settings are statistically independent of the coin:**
>
> | Tester | Joint P(c=1, s=1) vs product | Setting-0 statistic vs convex mixture |
> |---|---|---|
> | (i) `s = [n√2 mod 1 < μ]`, rationally independent of φ | 0.1200 = 0.1200 | **0.3000 = 0.3000** ✓ (also for λ = φ/2: 0.3090 = 0.3090) |
> | (ii) `s = [2nφ mod 1 < μ]`, rationally dependent | 0.2000 ≠ 0.1200 | 0.1667 ≠ 0.3000 ✗ |
> | (iii) tester reads the coin | 0.3000 ≠ 0.0900 | 0.0000 ≠ 0.3000 ✗ |
>
> **Prices:**
> - **access:** the coin must be hidden from the tester;
> - **independence:** the tester's choices must be uncorrelated with the coin.
>
> Independence is **earned** when the rotation numbers are rationally independent (Weyl joint equidistribution). That
> rational independence is itself a **condition** on the dynamics. Calling it "generic" is Lebesgue-priced: almost
> every pair is independent, but the dependent pairs are dense.

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| mixture weights | **forced** (uniquely ergodic frequencies) |
| convex combination of statistics | **forced given independence** |
| independence of settings from preparations ("free choice") | **inserted as a condition** (rational independence) — or Lebesgue-priced if called generic |
| inaccessibility of the coin | **inserted** (access) |
| "state = equivalence class over a test set" | **definitional**: the test set (access) is supplied |

## 2. Accounting

The supplied *randomness* of the operational framework is replaced by three supplied structures:
1. a hidden uniquely ergodic subsystem;
2. an inaccessibility (access) rule;
3. a dynamical independence condition between preparer and tester.

This is a genuine re-expression in **dynamical** terms (no probability primitive). It is not an elimination:
**ACCESS-PRICED + INDEPENDENCE-PRICED.**

**Observer loop.** The "tester" is a second deterministic subsystem (here, another rotation). The operational notion
of mixture needs **two subsystems with incommensurate dynamics**. That ties C5-D/E (states, probability) back to
C5-A (individuation of a preparer and a tester): **the observer is a subsystem whose independence must be supplied
or earned dynamically.**

**Status: S2-4 COMPLETE — convexity emerges deterministically with earned weights, iff the tester is dynamically
independent of a hidden coin (rational independence; access rule). Randomness is replaced by
access + independence, not eliminated.**
