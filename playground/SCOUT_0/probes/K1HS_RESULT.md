# SCOUT_0 K1-HS RESULT — the strong-gravity K1 branch, solved as a separatrix/boundary problem

**Question (external audit):** does the exact early-GR-restoring `μ > 1` K1 branch admit a globally
regular, EFT-stable continuation to `a = 1`? Π₀ and `p_tt` were not entered. K1-Z stays on hold. No
ruling on the sign of `x` is made.

**Scripts:**
- `k1hs_separatrix.py` — high-order series + forward integration, and backward shooting;
- `k1hs_readings.py` — exact K1 identities, two stability readings, filter-corrected shooting.

Logs alongside.

**Inputs and grades:**
- **K1 ODE:** from Linder 2020 `μ, Σ` (PRIMARY-TEXT-VERIFIED via the external audit).
- **EFT stability:**
  - ghost-free: `α = α_K + 3α_B²/2 > 0` (`α_K` free);
  - gradient: `N = αc_s² = D/2 − (2−α_B)H′/H − ρ_m/(H²M_*²) > 0`.
  - **VERIFIED against published Horndeski equations:** Peirone et al., PRD 97, 043519 (2018), App.
    Eq. A18, Bellini–Sawicki `α_B` convention, `α_T = 0` (via the external audit).
- Background: ΛCDM, `Ω_m0 = 0.3`.

**Terminology (repaired):**
- `EFT-STABLE` = ghost- and gradient-stable.
- `ODE-ATTRACTING` / `ODE-REPELLING` = sensitivity inside the K1 constraint ODE.

Neither is called "stability" without qualification.

## 1. Structure of the K1 ODE near `a = 0` (why forward IVPs mislead)

With `u = α_B/aⁿ` and `α_M ≈ c aⁿ`, the K1 ODE becomes, at leading order, the autonomous equation
`u̇ = (n/2c)(u − b₀)(u − b₁)` in `ln a`, with `b₀ < 0 < b₁`.

- **Strong branch `b₁`:** a **source**. It is ODE-REPELLING forward in `a`, but every GR-restoring
  history in `(b₀, ∞)` **emanates from it** as `a → 0`.
- **Weak branch `b₀`:** ODE-ATTRACTING forward. As a past asymptote it is the exceptional trajectory.
- **Linearized departure rate from `b₁`:** `√(4n+1)` in `ln a`. In `t = aⁿ` the exponent is
  non-integer (`√13/3` for `n = 3`). So there is a unique analytic (series) strong solution, plus a
  one-parameter family `K·a^{√(4n+1)}`, all GR-restoring.
- **Consequence:**
  - forward IVPs from a leading-order start are dominated by truncation error amplified at that rate
    (the K1-H blow-ups);
  - **backward** integration from `a = 1` is ODE-attracting toward `b₁` and well-conditioned;
  - K1-H's statement "branch 0 generic attractor, branch 1 repelling" is true only for forward-time
    sensitivity, **not** for which histories restore GR.

## 2. Exact identities on K1 (new)

From Linder's formulas (`k1hs_readings.log`): `μ = R(1 + A²/D)` and `Σ = R(1 + (α_M+α_B)A/D)`. On K1:

> **`μ − 1 = 2α_M(1−R)/α_B`, `Σ − 1 = α_M(1−R)/α_B`, `D = α_B A R/(1−R)`.**

For running histories with `α_M > 0` (so `R < 1`), **`sign(μ − 1) = sign(α_B)`**, and `μ` has a pole
wherever `α_B` crosses zero.

## 3. Construction (I): high-order series + forward integration (`k1hs_separatrix.log`)

- **Method:** strong-branch series `α_B = Σ b_k t^k`, `t = aⁿ`, to order 5, solved recursively from
  the exact ODE with exact `R(t)`. Then forward integration from `a₀ ∈ {10⁻¹, 3·10⁻², 10⁻²}` at
  orders 1, 3, 5.
- **Leading order (order 1) always blows up** (`a ≈ 0.18–0.81`), confirming the audit's diagnosis.
- **Orders 3 and 5 converge in both `a₀` and order:**

| `α_M` history | c | `α_B(1)` (converged) | `μ(1) − 1` | `min N` | outcome |
|---|---|---|---|---|---|
| `a¹` | 0.05 | +0.15358 | +0.0318 | > 0 | regular, EFT-stable, `μ > 1` |
| `a¹` | 0.2 | +0.52955 | +0.1369 | > 0 | regular, EFT-stable, `μ > 1` |
| `a^1.5` | 0.05 / 0.2 | +0.11579 / +0.40046 | +0.0283 / +0.1247 | > 0 | regular, EFT-stable, `μ > 1` |
| `a²/(1+a²)` (non-power-law) | 0.05 / 0.2 | +0.03983 / +0.13781 | +0.0216 / +0.0972 | > 0 | regular, EFT-stable, `μ > 1` |
| Ω_DE-tracking | 0.1 | +0.02158 (a₀ = 0.1, 0.03); +0.02147 (a₀ = 0.01) | +0.516 / +0.519 | > 0 | regular, `μ > 1`; convergence ~5·10⁻³ in `μ`, because `α_B(1)` is near 0 (pole proximity) |
| Ω_DE-tracking | 0.3 | −0.01493 (converged) | −6.35 | > 0 (BS) | the analytic separatrix **crosses `α_B = 0`**, so `μ` has a pole at finite `a` (§2). This one member is singular in `μ` |

## 4. Construction (II): backward shooting from `a = 1` (`k1hs_readings.log`, filter-corrected)

- **Method:** family labelled by `α_B(1) ∈ [0.05, 1.5]`, integrated back to `a = 10⁻⁴`. Accepted
  members have:
  - a strong-branch origin (early ratio within 2 % of `b₁`);
  - `α_B > 0` throughout, so no `μ` pole;
  - `R > 0`, `N > 0`, and K1 residual ≤ `10⁻¹⁴`.
- **Results:**

| history | c | EFT-stable strong-origin members with `μ > 1` | `μ(1) − 1` range |
|---|---|---|---|
| Ω_DE-tracking | 0.1 / 0.3 | 30/30 / 30/30 | +0.007 … +0.22 / +0.06 … +1.90 |
| `a¹` | 0.05 / 0.2 | 30/30 / 30/30 | +0.003 … +0.10 / +0.05 … +1.45 |
| `a^1.5` | 0.05 / 0.2 | 30/30 / 30/30 | +0.002 … +0.07 / +0.03 … +1.00 |
| `a²/(1+a²)` | 0.05 / 0.2 | 30/30 / 30/30 | +0.0006 … +0.017 / +0.009 … +0.27 |

- Values of `α_B(1)` below the weak trajectory are singular in the past (not GR-restoring).
- *(The first run of `k1hs_separatrix.py` reported "no member" for Ω_DE-tracking. That was a
  normalization bug in the branch filter — it compared to `b/c` instead of `b/(c/Ω_m0)`. It is fixed in
  `k1hs_readings.py`.)*

## 5. Verdict

> **A. STRONG-K1-VIABLE.** Exact, global `μ > 1` K1 histories exist that are early-GR-restoring,
> regular to `a = 1`, EFT-stable (ghost-free with `α_K ≥ 0`; gradient-stable by the verified `N`) and
> non-GR.
>
> - This holds in **every** tested history: Ω_DE-tracking, `a¹`, `a^1.5`, and the smooth non-power-law
>   `a²/(1+a²)`. In each there is a whole one-parameter family; in four of them the analytic
>   separatrix itself also qualifies.
> - The one exception is the analytic separatrix of the Ω_DE-tracking `c = 0.3` history: it crosses
>   `α_B = 0` (a `μ` pole), but the same history's shooting family still has 30/30 viable members.
> - **The `μ > 1` half of K1 is Horndeski-occupied.**
>
> **K1-H's "μ > 1: GR-ONLY-AFTER-STABILITY" is WITHDRAWN.** (It was already demoted to UNRESOLVED by
> the audit; K1-HS resolves it to A.)

## 6. A denominator-consistency question (owed; affects only the weak side)

Linder's `μ = R(1 + A²/D)` has the standard quasi-static structure `μ = R(1 + A²/(2αc_s²))`, with `D`
in the place of `2αc_s²`. The verified Horndeski expression gives
`αc_s² = D/2 − (2−α_B)H′/H − ρ_m/(H²M_*²)`, which differs from `D/2` by background/matter terms
(`3Ω_m(1 − R − α_B/2)` on a ΛCDM background normalized with `m_p`).

| Reading | Strong side (`α_B > 0`) | Weak side (K1-H, `−2α_M < α_B < 0`) |
|---|---|---|
| **BS** (stability = verified `N`) | EFT-stable (30/30) | EFT-stable (`min N > 0`) |
| **L** (`αc_s² = D/2`, as Linder's `μ` structure implies) | EFT-stable (30/30): on K1, `αc_s² = α_B A R/(2(1−R)) > 0` | **gradient-UNSTABLE** (`min D/2` from −0.08 to −0.23) |

- **Either way, the strong half is occupied.**
- The weak half's EFT status depends on which reading Linder's equations assume. That is an owed
  primary-text check: whether his Eqs. (8)–(9) are exact quasi-static expressions on a general
  background, or assume a background in which these terms cancel.

## 7. Consequence (as pre-registered by the audit)

> **K1 is Horndeski-occupied on the `μ > 1` side by EFT-stable, early-GR-restoring luminal histories.**
> **`sign(x)` does not rescue distinctiveness.** Even a derived `p_tt` with `x > 0` would land K1 on a
> line that stable Horndeski already populates, with `μ − 1 = 2α_M(1−R)/α_B`.
> **The `p_tt` route loses much of its empirical uniqueness.** K1 joins the campaign's other results
> as a universality/degeneracy relation.

What could still distinguish GRUT's K1 from the Horndeski family is **not** the line itself. It is
GRUT's additional structure: a time-constant ratio with `μ − 1 = xα` constant on the constant-`x`
cut, against Horndeski's history-dependent `μ − 1 = 2α_M(1−R)/α_B`. But the record itself calls
constant `x` a cut ("x is a kernel, not a constant"), so that would need its own derivation before it
could count. **Not pursued here.**

**Status: K1-HS COMPLETE. A. STRONG-K1-VIABLE across all tested histories (series-converged and
shooting families). K1-H's μ > 1 classification withdrawn. The weak side's EFT status is
reading-dependent (owed primary-text check). K1-Z remains on hold; no sign(x) ruling; Π₀/p_tt not
entered.**
