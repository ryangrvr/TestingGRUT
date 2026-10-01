# SCOUT-1 W2-FP RESULT — criticality without tuning? (BTW sandpile)

**Charter:** `PROBE_CHARTERS.md` §W2-FP. Preregistered outcomes: UNTUNED FIXED POINT / TUNING RELOCATED / NO.

**Files:** `w2fp_soc.py`, with log `w2fp_soc.log`.

**Setup:**
- Model: 2D BTW, open boundaries, parallel abelian toppling.
- Drive: 30 000 drives per size after a 3L² warm-up.

**Labels:**
- TUNING RELOCATED (to G-fixed points);
- KNOWN RESULT IMPORT: Bak–Tang–Wiesenfeld 1987; Dhar abelian sandpile; the "SOC needs conservation + slow drive" critique
  (Grinstein; Dickman–Muñoz–Vespignani–Zapperi) — STANDARD-TEXTBOOK ✓ (simulation);
- NON-DISTINCTIVE.

## 0. Verdict

> **TUNING RELOCATED — to exact zeros, which are exactly the G-fixed values the W1-C corollary allows a
> scale-free principle to select.**
>
> **The scale-free state.** With conservative bulk and infinitely slow drive, it arises with no dimensionful
> parameter tuned:
>
> | L | 16 | 32 | 64 |
> |---|---|---|---|
> | cutoff `s_c` | 92 | 464 | 2808 |
>
> So `s_c ∝ L^{2.46}`, with `⟨s⟩ ∝ L^{1.84}` (exact BTW: L²; finite-size).
>
> **But two control parameters must sit at exactly zero.**
> - **Bulk dissipation ε = 0.** At ε = 0.01 / 0.03 / 0.1, `s_c` = 183 / 58 / 18. At ε = 0.03 the cutoff is the
>   same for L = 32 and 64 (49.5 vs 57.5): **set by ε, not by L**. Criticality is lost.
> - **Drive rate r = 0** (infinite separation of time scales). At r = 10⁻⁴, `⟨s⟩` inflates from 96 to 240. At
>   r ≥ 10⁻³, activity never stops (supercritical, merged avalanches).
>
> **The prices.**
> - ε = 0 is an **exact conservation law** (a supplied symmetry/sector datum, S-7 class).
> - r = 0 is an **exact time-scale separation** (a supplied limit).
>
> Both are G-fixed values (zero rates), which is why they can be imposed by scale-free principles, whereas
> "h = h_c" (W1-R) cannot.

## 1. What this means

- **The tuning price of W1-R is not removed; it is moved.** A dimensionful tuning (`h − h_c = 0`, a non-fixed
  point of G) becomes **two exact zeros**, each a symmetry/limit statement. This sharpens the W1-C fixed-point
  corollary into a usable rule:

  > **Untuned criticality = criticality whose control parameters are G-fixed values (0 or ∞) imposed by
  > supplied exact symmetries or limits.**

- **GRUT relevance.** GRUT's parents have exact conservation where a sector is supplied (S-7), and the
  quadratic floor has no avalanche dynamics. A GRUT route to "earned criticality" would need its own exact
  conservation law **and** an infinite scale separation, and both are supplied in the current record. Note
  that the record's pin = 0 fork (E-2 Outcome A, held outside the floor) is precisely such a G-fixed value. It
  is the natural place an SOC-type principle would act. That is an **owner-level** fork and is not entered.
- **Hostile notes:**
  - The exponent fits (−0.99…−1.13 in short windows) are finite-size. The BTW τ is known to be ill-defined /
    multiscaling, and the claim does not depend on it.
  - The cutoff scaling is the robust signal.

**Status: W2-FP COMPLETE — TUNING RELOCATED to exact zeros (conservation ε = 0, separation r = 0): supplied
G-fixed-point premises. No untuned selection.**
