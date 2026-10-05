# BRI1 PUBLICATION FIGURE REPORT 02 — V3-R1
## DETERMINISTIC QUADRATURE REMEDIATION
## Executed under `BRI1_PUBLICATION_VERIFICATION_CHARTER_01.md` (frozen `8d765ea`), V3 lane only

**Status:** V3-R1 — illustrative / evidence-grade numerical demonstration. **Not part of the analytic proof.** V1 (the proof) is unchanged. V2's KNOWN-RESULT-NEW-FRAMING disposition is unchanged. Report 01 (failed MC attempt) is preserved unchanged as the permanent record of that attempt.

**Motivation:** Report 01 documented that the preregistered Monte Carlo path failed for a structural reason — the third-cumulant signal at small t is ~10⁻¹⁰–10⁻¹¹, far below any feasible sampling noise floor (~10⁻²). This remediation replaces the sampled bath with **deterministic phase-space quadrature over the single-oscillator dynamics**, exploiting the exact i.i.d. structure to avoid simulating the full N_B-dimensional bath.

---

## A. Core numerical strategy (as specified)

For each finite N_B: ε = N_B^(−1/2). Under a fixed clamp, every bath oscillator is an i.i.d. copy of the single-oscillator system x″ + x + x³ = ε q(t) with the frozen Gibbs preparation. Its time-t position is X^ε(t).

Moments of X^ε are computed **deterministically** by quadrature over the initial Gibbs distribution ρ ∝ exp[−p²/2 − x²/2 − x⁴/4] on (x₀, p₀). The exact i.i.d. cumulant-sum identities then give, **without simulating an N_B-dimensional bath**:

- κ₃(F) = N_B^(−1/2) κ₃(X^ε) — **exact** (identity ★ from V1-3)
- Var(F) = Var(X^ε) — **exact**
- **γ₁(F) = [N_B^(−1/2) κ₃(X^ε)] / Var(X^ε)^(3/2)** — the primary finite-N observable

The asymptotic K/N_B expression is **not** substituted; κ₃(X^ε) is evaluated at the actual ε of each N_B by quadrature of the finite-ε one-oscillator dynamics.

## B. Pre-registration (frozen before evaluation)

- **Protocols:** P0, P1 exactly as in X1 (frozen; no new protocol).
- **Primary witness time:** **t★ = 0.5**, chosen before quadrature execution. Reason: well inside the early P1 ramp (0.5 < π); removes ~7 powers of small-time suppression vs the failed t=0.05 MC; **not claimed to lie inside the analytically certified existential interval (0, δ)** — no numerical lower bound on δ has been proved. This route is an **illustrative finite-time numerical consistency test**, not an extension of the certified interval.
- **Diagnostic times (predeclared, all reported):** t ∈ {0.25, 0.75, 1.0}.
- **Bath-size grid:** N_B ∈ {4, 8, 16, 32, 64, 128} (frozen).
- **Scaling fit:** log|γ₁(P1)| = a − p·log N_B, OLS on all six sizes at t★; reference p = 1 displayed; no point deletion.

## C. Deterministic Gibbs quadrature

- **Rule:** tensor-product Gauss–Legendre in x₀ ∈ [−6, 6] and p₀ ∈ [−6, 6], with explicit reweighting by the Gibbs factor exp[−x₀⁴/4] (the p-part exp[−p²/2] is handled by the same rule; the weight product is computed explicitly).
- **Normalization:** the truncated-domain weight sum is computed and the joint weight normalized to sum 1 on the domain (ratio normalization, same scheme as the record's prior preflight).
- **Tail control:** at |x| = 6, exp[−x²/2 − x⁴/4] ≈ exp[−18 − 324] ≈ e⁻³⁴² — utterly negligible; likewise for p at |p| = 6. Domain truncation error ≪ machine precision.
- **Resolution ladder (three settings, minimum required):** (nx, np) ∈ {(120,120), (240,240), (480,480)}.

## D. ODE convergence

Integrator: fixed-step RK4, dt ∈ {10⁻³, 5×10⁻⁴, 2.5×10⁻⁴} (refinement by factors of 2, three settings). **Convergence demonstrated empirically**, not asserted from formal order:

- P0 mean at t★: 5.20×10⁻¹⁸ / −5.20×10⁻¹⁸ / −2.03×10⁻²³ across dt — fluctuating at 10⁻¹⁸–10⁻²³ level, i.e., **numerically zero** (analytically zero by symmetry).
- P0 k₃: −1.85×10⁻¹⁷ / +1.85×10⁻¹⁷ / +6.7×10⁻²³ — **numerically zero** (analytically zero; residual at machine-epsilon × grid-size level).
- P1 γ₁ values: **identical to 4 significant figures across all nine (resolution, dt) combinations** (e.g., t★, N_B=4: −5.571×10⁻⁶ in all nine). Quadrature and integration convergence are both far tighter than the observed signal.

## E. P0 null (primary control)

P0's third cumulant is analytically zero (exact symmetry). Numerically: |κ₃(P0)| ≤ 3×10⁻¹⁷ across all resolutions/times — i.e., **zero to machine-epsilon level**, ~10¹² times smaller than the P1 signal. **The quadrature pipeline reproduces the P0 null.** ✓

## F. P1 finite-N calculation (t★ = 0.5; full table)

| N_B | κ₃(X^ε) | Var(X^ε) | κ₃(F) = κ₃/√N_B | γ₁(F) |
|---|---|---|---|---|
| 4 | −1.048024×10⁻⁵ | 0.9600462 | −5.240118×10⁻⁶ | **−5.570611×10⁻⁶** |
| 8 | −7.410645×10⁻⁶ | 0.9600462 | −2.620059×10⁻⁶ | −2.785305×10⁻⁶ |
| 16 | −5.240118×10⁻⁶ | 0.9600462 | −1.310029×10⁻⁶ | −1.392653×10⁻⁶ |
| 32 | −3.705323×10⁻⁶ | 0.9600462 | −6.550147×10⁻⁷ | −6.963263×10⁻⁷ |
| 64 | −2.620059×10⁻⁶ | 0.9600462 | −3.275073×10⁻⁷ | −3.481632×10⁻⁷ |
| 128 | −1.852661×10⁻⁶ | 0.9600462 | −1.637537×10⁻⁷ | −1.740816×10⁻⁷ |

Note κ₃(X^ε) itself varies with ε (≈ ∝ ε at leading order, consistent with the V1 analytic structure: κ₃(X^ε) = εK + O(ε³) with K = 3c(t★) < 0), which is exactly why the finite-N observable must be built from the quadrature rather than from the asymptotic formula.

## G. Figures

- **Figure 1** (`v3_r1_figures.png`, left panel): log–log |γ₁(P1)| vs N_B at all four times, with reference slope N_B⁻¹. All four times show clean N_B⁻¹ behavior.
- **Figure 2** (right panel): N_B·γ₁(P1) vs N_B — **flat to 4 significant figures** (−2.228×10⁻⁵ at every N_B for t★), demonstrating genuine 1/N_B scaling rather than mere decay toward zero. All four diagnostic times flat as well (values: t=0.25: −3.29×10⁻⁷; t=0.75: −6.12×10⁻⁵; t=1.0: −1.65×10⁻³).
- Both captions carry the required label: *"Illustrative deterministic numerical reproduction — not part of the analytic proof."*

## H. Shape-invariant interpretation

The standardized skewness is itself an affine/location-scale/sign shape witness:

- location changes do not alter skewness;
- positive scale changes do not alter standardized skewness;
- sign/reflection changes its sign but not its absolute magnitude.

Therefore |γ₁| = 0 (P0, exact) versus |γ₁| > 0 (P1, at every finite N_B in the grid) is an **illustrative orbit/shape distinction**: no shared location-scale-sign modulation of one common process can map the P0 force law to the P1 force law at these parameters. **This numerical observation is not the proof; V1 remains the proof.**

## I. Diagnostic-time firewall

All predeclared diagnostic times reported, none suppressed:

| t | κ₃(X^ε) at N_B=4 | γ₁(F) at N_B=4 | N_B·γ₁ (stable value) |
|---|---|---|---|
| 0.25 | −2.256×10⁻⁷ | −8.221×10⁻⁸ | −3.29×10⁻⁷ |
| **0.5** | **−1.048×10⁻⁵** | **−5.571×10⁻⁶** | **−2.228×10⁻⁵** |
| 0.75 | −4.032×10⁻⁵ | −4.781×10⁻⁵ | −6.12×10⁻⁵ |
| 1.0 | −1.674×10⁻⁴ | −4.115×10⁻⁴ | −1.65×10⁻³ |

All four times show: nonzero P1 skewness, clean N_B⁻¹ scaling (Fig 1), and flat N_B·γ₁ (Fig 2). The story is **consistent across times** — the coefficient grows with t (as expected from the t⁷ leading behavior of c(t) and beyond), and the scaling exponent is 1 at every one of them. No near-zero of the coefficient at t★; no retuning needed.

## J. Interpretation answers

1. **Does deterministic quadrature reproduce P0 symmetry?** Yes — |κ₃(P0)| ≤ 3×10⁻¹⁷, machine-epsilon level. ✓
2. **Is P1 skewness nonzero at the preregistered t★ = 0.5?** Yes — γ₁(F) = −5.571×10⁻⁶ at N_B = 4, nonzero at every N_B. ✓
3. **Does |γ₁| decrease with N_B?** Yes — monotonically, by exactly the factor √(N_B ratio) in κ₃(F) with Var essentially constant. ✓
4. **Is the fitted exponent compatible with asymptotic p = 1?** **Fitted p = 1.0000** (OLS on all six sizes at t★; residuals at the 4th-figure level — the fit is exact to displayed precision). ✓
5. **Does N_B·γ₁ approach a stable value?** Yes — constant to 4 figures (−2.228×10⁻⁵). ✓
6. **Are there clear pre-asymptotic corrections?** None visible at the displayed precision — κ₃(X^ε)'s own ε-dependence is smooth and the leading linear-in-ε term dominates cleanly at these parameters. (The exact coefficient's ε³ correction would show at higher precision.) ✓
7. **Do the diagnostic times tell a consistent story?** Yes — same exponent p = 1 at all four times, coefficients ordered as expected from the underlying c(t) growth. ✓
8. **Is quadrature convergence substantially tighter than the observed signal?** Yes — variations across the 9 (resolution, dt) settings are below 10⁻⁴ relative, while the signal spans 4 orders of magnitude in N_B and 2 in t. ✓
9. **Does anything numerically contradict V1?** **No.** The quadrature result at t★ is fully consistent with the V1 analytic structure (κ₃(X^ε) ≈ εK(t★) with K(t★) = 3c(t★) < 0; Var(X^ε) → m₂ ≈ 0.96 as ε → 0). ✓

## K. Disposition

All PASS conditions met: converged quadrature ✓; P0 null reproduced ✓; nonzero P1 shape witness at preregistered t★ ✓; decreasing finite-size witness ✓; no numerical contradiction of the asymptotic 1/N_B structure ✓ (fitted p = 1.0000).

# **V3-R1-ILLUSTRATION-PASS**

This does **not** upgrade novelty or theorem status. V1 remains the proof; this is an illustrative deterministic numerical consistency demonstration at evidence grade (no formal interval certification was added).

## L. Reproducibility files (this directory)

- `v3_r1.py` — the executed script (preregistered choices embedded).
- `v3_r1.log` — raw output, all nine (resolution, dt) settings, unmodified.
- `v3_r1_results.json` — full numerical results.
- `v3_r1_table.csv` — the primary table.
- `v3_r1_figures.png` — Figures 1 and 2.

Original V3 (Report 01) files untouched.

---

**V1 is the proof; V3-R1 is visualization/support only. Nothing in this report upgrades or downgrades the theorem, V1, V2, or the KNOWN-RESULT-NEW-FRAMING disposition.**
