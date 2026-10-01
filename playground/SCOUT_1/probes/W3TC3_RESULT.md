# SCOUT-1 W3-TC3 RESULT — what Gaussian coarse-graining does to the quotient data

**Charter:** ZOOM_OUT_03 §4–5 (TC-3 conjecture).

**Files:** `w3tc3_gaussian_rg.py`, with log `w3tc3_gaussian_rg.log`.

**Labels:**
- TC-3 **REFINED** (the original conjecture was too strong);
- FAMILY THEOREM (parts 1 and 3; standard);
- KNOWN RESULT IMPORT: Schur complement / Feshbach; Gaussian fixed points; long-range hopping — STANDARD-TEXTBOOK ✓.

## 0. Verdict

> **TC-3 as conjectured ("every Gaussian quotient datum is exactly marginal; no Gaussian mechanism reduces data")
> is too strong, and is refined.**
>
> 1. **Exact decimation preserves the retained quotient.** Schur-complement decimation of every other site
>    reproduces `G_rr(z)` to ≤ 3·10⁻¹⁵. `μ_r` does not flow. This is E-3's "exact reduction".
> 2. **Gaussian IR flow *does* forget data.** Analytic higher-order dispersion terms are irrelevant. NNN hopping
>    with a quadratic minimum keeps the bulk γ = −0.499 / −0.498 / −0.498 for t₂ = 0 / 0.1 / 0.2. (t₂ = −0.2 sits
>    near the quartic point and shows a crossover at −0.546.) This is a genuine reduction of *irrelevant* data.
> 3. **But Gaussian fixed points form a continuum.**
>    - Fine-tuning t₂ = −1/4 gives a quartic minimum, γ = −0.739 (→ −3/4).
>    - Long-range hopping `J ∝ r^{−(1+s)}` gives dispersion exponents 0.51 / 0.99 / 1.44 / 1.65 for
>      s = 0.5 / 1 / 1.5 / 1.8, so `γ = 1/s − 1` varies continuously. The truncation at r_max = 2·10⁴ and the
>      log correction at s = 1 explain the deviations.
>    - Together with W1-S (radial dimension D, boundary class), the labels s, D and k are **fixed-point labels,
>      not flowing couplings**. Decimation cannot change a non-analytic tail (standard; asserted, not separately
>      computed).
>
> **Refined TC-3 (FAMILY THEOREM grade for parts 1 and 3, within the tested families):**
> - In the Gaussian class, coarse-graining forgets irrelevant analytic data but **never selects among fixed
>   points**: the IR labels form a continuum and nothing flows between them.
> - Selecting an **isolated** fixed point, and with it fixing a weight-0 IR datum without supplying its label,
>   requires interactions: W1-R (unitarity rigidity), W3-NL (KPZ), W2-DT (transmutation).

## 1. Consequence for GRUT

- GRUT's earned floor is Gaussian. So **any weight-0 IR datum it carries is a fixed-point label that must be
  supplied**: graph dimension, dispersion order, long-range tail, boundary class. This explains EDA-01's
  READ / SUP pattern from the RG side.
- The E-3 wall is now stated precisely:

  > **Linearity ⇒ continuum of fixed points ⇒ no dynamical selection of IR labels.**

**Status: W3-TC3 COMPLETE — TC-3 REFINED: the Gaussian class forgets irrelevant data but has a continuum of fixed
points; isolated fixed points (hence label selection) need interactions.**
