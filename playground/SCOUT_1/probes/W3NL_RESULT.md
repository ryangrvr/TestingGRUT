# SCOUT-1 W3-NL RESULT — does a nonzero nonlinearity collapse the IR quotient? (EW → KPZ)

**Charter:** ZOOM_OUT_03 §5.

**Files:** `w3nl_kpz.py`, with log `w3nl_kpz.log`.

**Setup:**
- Model: exclusion process on a ring ↔ single-step interface. The symmetric case (ε = 0) has a linear current; the
  asymmetric case (ε > 0) has a nonlinear current `j ∝ ε ρ(1−ρ)`.
- Flat start.
- Main runs: L = 8192, t ≤ 4000, 6 runs. Extended runs: L = 16384, t ≤ 40000, 3 runs.

**Labels:**
- CLASS COLLAPSE / **classical transmutation (NR, scoped)**;
- KNOWN RESULT IMPORT: KPZ 1986; ASEP ↔ KPZ; EW — STANDARD-TEXTBOOK ✓ (simulation; finite-time / statistical
  caveats);
- NON-DISTINCTIVE.

## 0. Verdict

> **Any nonzero nonlinearity drives the IR to a single fixed point (KPZ, β = 1/3). Its strength survives only as a
> crossover scale and an amplitude, which are G-moved (TC-1-unfixable) data. This is a classical analogue of
> dimensional transmutation.**

| ε | β[10–100] | β[100–1000] | β[1000–4000] | extended, latest window |
|---|---|---|---|---|
| 0 (SSEP) | 0.226 | 0.233 | 0.255 | — (EW 1/4) |
| 0.25 | 0.238 | 0.261 | 0.265 | β[15000–40000] = 0.388 (noisy; still crossing) |
| 0.5 | 0.249 | 0.296 | 0.351 | β[15000–40000] = 0.321 |
| 1.0 (TASEP) | 0.294 | 0.337 | 0.328 | — |

- The crossover moves to later times as ε decreases, as expected from `t_x ∝ ε^{−4}`.
- The late-time exponent approaches 1/3 for every ε > 0. The ε = 0.25 extended windows are noisy (3 runs); the
  table is evidence of crossover, not a precision exponent fit.
- The amplitude `W³/t` at t = 4000 is 0.051, 0.097, 0.172 for ε = 0.25, 0.5, 1. So ε sets the *scale*, while the
  exponent is ε-independent.

## 1. Information accounting

- **Linear class (ε = 0).** The exponent is also universal here (1/4), but only within short-range linear models.
  The linear / Gaussian class as a whole carries a **continuum** of IR data: any kernel or `μ_r` (W1-S radial
  chains: γ = D/2 − 1 for any real D; Lévy kernels: continuous z).
- **Nonlinear class.** The relevant nonlinearity removes the dimensionless microscopic ratio ε from the IR
  (**NR: ε consumed**). The basin is wide: short-range interactions, 1D, a conserved field with a nonlinear
  current.
- **Prices (supplied):**
  - the conservation law;
  - spatial dimension;
  - short-range noise and couplings (long-range-correlated noise gives continuously varying KPZ exponents, which
    is standard);
  - the existence of the nonlinearity itself.

## 2. GRUT relevance and honesty fence

- The record's **declared** nonlinear class (S2 C-B nonlinear drift, E-15) is a 0D (single retained coordinate)
  stochastic law. **The field extension used here is a SCOUT construction**, not a record object. So this shows
  what a nonlinear extension *would* do, not what the record does.
- Within the record, the S2 discriminator works at O(β), i.e. linear response in the nonlinearity (`Δc₂ = −24βT₁a`).
  W3-NL shows that in an extended system the nonlinearity would be **relevant**: its strength would reduce to a
  crossover time, and its IR consequence would become a universal exponent independent of β.
  **Consequence:** an earned field-level nonlinearity would trade the S2 β-dependence (a free coupling) for a
  universal IR class. This would be a GRUT-internal NR, **but only after admitting a spatially extended nonlinear
  conserved field**, which is a premise change (owner level).
- **Third instance of the E-3 wall.** The reduction appears precisely when linearity is abandoned.

**Status: W3-NL COMPLETE — CLASS COLLAPSE / classical transmutation: any ε > 0 flows to KPZ (β → 1/3); ε survives
only as a scale. NR scoped to the nonlinear extended class, which the record does not contain (premise change).**
