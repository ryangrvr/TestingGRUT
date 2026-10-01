> **REPAIRED (external audit + `K1HS_RESULT.md`):**
> (1) **Stability source grade:** RECONSTRUCTED → **VERIFIED against published Horndeski equations** (Peirone et al.,
> PRD 97, 043519 (2018), App. Eq. A18, BS convention, α_T = 0).
> (2) **"μ > 1: GR-ONLY-AFTER-STABILITY" — WITHDRAWN.** "ODE-repelling" was conflated with "physically unstable", and the
> forward IVPs started from a leading-order coefficient, so truncation error was amplified at rate √(4n+1). K1-HS solves the
> branch properly (high-order series + backward shooting): **A. STRONG-K1-VIABLE** — EFT-stable, early-GR-restoring,
> non-GR μ > 1 K1 histories exist in every tested history.
> (3) **Terminology:** "attractor/repelling" below means ODE-ATTRACTING/ODE-REPELLING (forward in a), not EFT stability. As a
> past asymptote the strong branch is the generic GR-restoring origin and the weak branch the exceptional one.
> (4) Final audit ruling: there is only one stability criterion (`N`; Linder's own c_s² reduces to it exactly,
> `k1hs_denominator.log`), so "HORNDESKI-OCCUPIED for μ < 1" holds in the tested EFT-stable histories, with no reading caveat. The text below is kept for the audit trail.

# SCOUT_0 K1-H RESULT — is `2Σ − μ = 1` a viable non-GR luminal Horndeski submanifold?

**Question (external audit):** does the exact K1 constraint admit a globally viable, stable, non-GR
luminal Horndeski history with early-time GR restoration? This was run **before** any K1-Z forecast.
Π₀ and `p_tt` were not entered.

**Scripts:**
- `k1h_constraint.py` — symbolic re-derivation;
- `k1h_solutions.py` — asymptotics plus integration;
- `k1h_scan.py` — scan over histories.

Logs alongside.

**Inputs and grades:**
- Linder's `μ, Σ` are PRIMARY-TEXT-VERIFIED (external audit).
- The stability expression is **RECONSTRUCTED**, not primary-text-verified here. It is the
  Bellini–Sawicki form with `α_T = 0`:
  - `Dkin = α_K + (3/2)α_B² > 0`; `α_K` is free, so this is always satisfiable;
  - `N ≡ Dkin·c_s² = ½D − (2−α_B)H′/H − ρ_m/(H²M_*²)`.
  - Checks: in the quintessence limit it gives `c_s² = 1`, and its `α_B′` combination is Linder's
    denominator `D`.
- Background: ΛCDM, `Ω_m0 = 0.3`. Quasi-static, `α_T = 0`.

## 1. The constraint (step 1)

`2Σ − μ = R(2A + 2α_B′)/D`. So K1 ⟺ `2(R−1)α_B′ + [2(R−1) + α_B](α_B + 2α_M) = 0`, with `R′ = −α_M R`.

**Checks:**
- Only Run (`α_B = 0`): `2Σ − μ = R`.
- No Run (`α_M = 0`): `Σ = μ`.
- At `R = 1`: `α_B(α_B + 2α_M) = 0`.

## 2. Boundary conditions (step 2)

Early-time GR restoration: `α_M, α_B → 0` and `R → 1` as `a → 0` (Linder's benchmark convention).
Stability: `Dkin > 0` (choose `α_K ≥ 0`) and `N > 0` for all `a`. No present-day normalization of `R`
is imposed.

## 3. Solutions (steps 3–4)

**Early-time asymptotics (exact at leading order).** Take `α_M = c aⁿ`, `α_B = b aⁿ`,
`R − 1 = −(c/n)aⁿ`. The leading-order condition is `n b² − 2cb − 4c² = 0`, so
**`b = c(1 ∓ √(4n+1))/n`** — two non-GR branches that reach GR as `a → 0`.

| | branch 0 (`b/c = (1−√(4n+1))/n`) | branch 1 (`b/c = (1+√(4n+1))/n`) |
|---|---|---|
| `n = 3` ratio `α_B/α_M` | −0.869 | +1.535 |
| sign of `μ − 1` (for `c > 0`) | **negative** (`−0.768 c a³` at `n = 3`) | **positive** (`+0.434 c a³`) |
| early stability `N` (for `c > 0`) | positive for `n ≲ 4` (`+0.83 c a³` at `n = 3`; `+4.38 c a` at `n = 1`) | positive |
| dynamics | **attractor**: ±5 % perturbations converge onto it | **repelling separatrix**: `α_B′ ≈ α_B²/(2|R−1|)` gives a finite-time Riccati blow-up for any start above it; starts below fall to branch 0 |
| for `c < 0` (`α_M < 0`) | `μ > 1`, but `N < 0`: gradient-unstable | unstable |

Both branches differ from the named Limited-Modified-Gravity subclasses: No Slip has `α_B/α_M = −2`,
Only Run `0`, No Run `α_M = 0`.

**Numerical integration to `a = 1`** (`k1h_solutions.log`, `k1h_scan.log`). K1 residual
`|2Σ−μ−1| ≤ 10⁻¹³` on every completed run:

| `α_M` history | c | branch 0 at `z = 0` | stability (`min N`) | branch 1 |
|---|---|---|---|---|
| `c·Ω_DE(a)/Ω_Λ` | 0.1 | `μ−1 = −0.092`, `Σ−1 = −0.046` | **> 0** throughout | blows up |
| `c·Ω_DE(a)/Ω_Λ` | 0.3 | `μ−1 = −0.259`, `Σ−1 = −0.129` | **> 0** throughout | blows up |
| `c·a` | 0.05 / 0.2 | `μ−1 = −0.079 / −0.290` | **> 0** | blows up (−5 % start → branch 0) |
| `c·a^1.5` | 0.05 / 0.2 | `μ−1 = −0.060 / −0.225` | **> 0** | blows up |
| `c·a³` (pure power law) | 0.05 / 0.2 | `μ−1 = −0.038 / −0.147` | **< 0** late (Λ era): unstable | blows up |
| any of the above | −0.05 | `μ−1 > 0` | **< 0**: unstable | unstable or blows up |

## 4. Verdict (steps 5–6)

> **Sign-split.**
>
> - **`μ < 1` half-line (weak gravity): HORNDESKI-OCCUPIED.** For `α_M > 0` the K1 constraint
>   defines a generic (attracting), early-GR-restoring, gradient-stable, non-GR family of luminal
>   Horndeski histories. Examples: Ω_DE-tracking, `a¹` and `a^1.5` running. It satisfies
>   `Σ − 1 = (μ − 1)/2` identically in time and is distinct from GR and from every named LMG subclass.
>   Stability depends on the history (pure `a³` fails late).
> - **`μ > 1` half-line (strong gravity) — where GRUT's interior family sits for `x > 0`:
>   GR-ONLY-AFTER-STABILITY within the scanned class.** Every `μ > 1` K1 history found is either
>   gradient-unstable (`α_M < 0`) or the repelling Riccati separatrix (branch 1). The separatrix blows
>   up numerically. Whether its exact trajectory survives to today is **not established**; it is
>   non-generic in any case.

**Grade:**
- leading-order analytic (exact) plus a numerical scan (3 power laws, the Ω_DE-tracking shape,
  `c = ±0.05, ±0.2, 0.1, 0.3`);
- reconstructed stability formula; ΛCDM background; quasi-static; `α_T = 0`.
- **Not a theorem.** Not covered: histories with sign-changing `α_M`, non-ΛCDM backgrounds,
  beyond-Horndeski, DHOST, nonlocal, and an exact treatment of the separatrix.

## 5. What this does to K1 and to the `p_tt` question

- **The sign of `x` is now load-bearing.** The record's family is "between the banked endpoints"
  (`x ∈ [0, 1]`, so `μ ≥ 1`). But the `x ≥ 0` orientation lemma was **RULED UNBANKED**
  (`X_FLOOR_MAP.md`).
  - With `x > 0`, K1 lies on the half-line that stable, early-GR-restoring luminal Horndeski did not
    occupy in this scan. That makes K1 **a sharper conditional discriminator than previously
    thought**: on the line *and* `μ > 1`.
  - With `x < 0` allowed, K1 is Horndeski-degenerate. It is another universality relation.
- **Observable reading (conditional):** a detection on the line `Σ − 1 = (μ − 1)/2` with `μ < 1`
  would be consistent with stable luminal Horndeski. A detection on the line with `μ > 1` would not
  be produced by any stable early-GR-restoring luminal Horndeski history in the scanned class.
- **For the owner:** a `p_tt` derivation now needs to settle **three** things for K1 to become a
  distinctive prediction:
  1. the sector assignment (`p_tt`);
  2. the inherited-conditional slip `η = 1/μ` (R2);
  3. **the sign of `x`**.

  Without (3), deriving `p_tt` could land K1 on the Horndeski-occupied half.

**Before K1-Z (owed):**
- primary-text verification of the stability expression;
- a scan widening (sign-changing `α_M`, other backgrounds);
- the exact separatrix fate;
- the beyond-Horndeski/DHOST classes.

K1-Z stays on hold.

**Status: K1-H COMPLETE (one pass): HORNDESKI-OCCUPIED for `μ < 1`; GR-ONLY-AFTER-STABILITY for `μ > 1`
within the scanned class; sign of `x` exposed as load-bearing.**
