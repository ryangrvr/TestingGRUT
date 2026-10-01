# SCOUT_0 K1 OCCUPANCY — CORRECTION 01 (external audit with primary text; the original stays visible)

**Corrects:** `K1_OCCUPANCY_01.md` (commit `cbe300b`). That file is kept unchanged below its new banner,
so the error stays in the audit trail.

**Source grade change:** Linder, "Limited Modified Gravity" (JCAP 10 (2020) 042, arXiv:2003.10453),
formulas for Only Run (`G_matter = R(1+α_M)`, `G_light = R(1+α_M/2)`, `R = m_p²/M_*²`) and for general
luminal Horndeski (`μ, Σ` with `A = α_B + 2α_M`, `D = (2−α_B)A + 2α_B′`, `R′ = −α_M R`):
**IMPORTED-SECONDARY → PRIMARY-TEXT-VERIFIED** (by the external audit; the sandbox itself still cannot
reach arXiv). The small-α Horndeski ratio `(1+r)/(2+r)` remains secondary. It is now shown to be
exactly the `R = 1`, `α_B′ = 0`, leading-order limit of Linder's formulas (`k1h_constraint.log`).

## Error 1 — Only Run "on the line at z = 0 for every amplitude": WITHDRAWN

- **Exact relation:** `2Σ − μ = R`, so Only Run lies on K1 (`2Σ − μ = 1`) **iff `R = 1`**, i.e.
  `M_*² = m_p²` (re-derived: `k1h_constraint.log`).
- **What I did wrong:** I identified `R = 1` with `z = 0` without checking the model's normalization.
  Linder's benchmark restores GR in the **early** universe: `M_*²/m_p²` runs from 1 in the past to
  `e^{4c_M/τ}` in the future.
- **Correct reading:** the standard Only Run history is on K1 only near the GR-restoration epoch and
  **off** it after the Planck mass has run.
- **Withdrawn claim:** "low-redshift amplitude-only data cannot distinguish GRUT from Only Run". A
  present-day normalization `R₀ = 1` would be a separate assumption, and it would have to be checked
  for consistency with early-time GR restoration. Not done.

## Error 2 — "holding K1 at all `a` forces GR within leading-order Horndeski": WITHDRAWN

- **What I did wrong:** the claim used the small-α ratio, which drops `α_B′` and `R − 1`.
- **Exact K1 condition** (re-derived independently; the numerator equals minus the audit's form, so
  the zero set is the same):

  `2(R−1)α_B′ + [2(R−1) + α_B](α_B + 2α_M) = 0`.

  This is a first-order ODE for `α_B` given `α_M(a)`. At `R = 1` it reduces to
  `α_B(α_B + 2α_M) = 0` (Only Run and No Slip as instantaneous roots). It does **not** force GR. The
  question is now answered by K1-H (`K1H_RESULT.md`).

## What survives from K1_OCCUPANCY_01

- The convention reconstruction (§1): K1 ⟺ `μη = 1` ⟺ `2Σ − μ = 1`; inherited-conditional slip;
  constant-x cut.
- The classes with `Σ₀ = 0` (f(R), BD, DGP, coupled quintessence) and `Σ₀ = μ₀` (No Slip, No Run, cubic,
  clustering DE), which meet the line only at the origin.
- The non-evaluated list.

**Counterexample ledger:** CE-07.
