# BRI1 PUBLICATION FIGURE REPORT 01 — V3
## Executed under `BRI1_PUBLICATION_VERIFICATION_CHARTER_01.md` (frozen `8d765ea`), V3 lane only
## Illustrative numerical reproduction — supporting evidence only; V1 is the proof

---

## §1. Preregistration (recorded BEFORE any run)

Per charter §1, the following were fixed before execution:

- **Protocols:** frozen X1 family {P0, P1} only (P0: q ≡ 0; P1: quintic ramp
  s(t/π) on [0, π], then q = 1). No new protocol invented.
- **Witness time:** t_num = 0.05 (Path A: small-time regime, justified
  below).
- **Bath-size grid:** N_B ∈ {4, 8, 16, 32, 64}.
- **Sample count:** 200,000 Gibbs draws per (N_B, protocol) (original
  preregistration); reduced to 20,000 mid-run for compute feasibility —
  documented in §4 as a deviation with cause, not as a silent change.
- **Seed policy:** single fixed seed 20260101 (numpy PCG64 Generator).
- **Integrator:** RK4 fixed-step, dt = 1/512 (26 steps to t_num = 0.05).
  Formal local truncation O(dt⁵); formal global order O(dt⁴). **No actual
  step-refinement/convergence check was performed in this run, so no
  numerical error bound follows from the formal order alone**; the
  integration error is expected to be small relative to sampling noise
  but this was not verified.
- **Statistic:** third central moment κ₃(F_q(t_num)) and standardized
  skewness γ₁(F_q(t_num)).
- **Slope fit:** OLS of log|γ₁(P1)| on log(N_B). **Note: any slope
  produced from the noise-dominated output is not an estimate of the
  physical exponent and must not be interpreted as such.**
- **Exclusion criteria:** none — all results reported.
- **Expected behavior:** |γ₁(P1)| ∝ 1/N_B asymptotically; γ₁(P0) ≈ 0.

### Path A justification for t_num = 0.05 (as preregistered)

The V1-6 result c(t) = C₇t⁷ + O(t⁸) with C₇ = −Var(x₀²)/(28π³) < 0
establishes the sign for t ∈ (0, δ) where δ is existential. t_num = 0.05
is chosen on the structural expectation that the t⁷ term dominates t⁸ at
very small t (ratio t·C₈/C₇), making it a natural small-time probe. **This
is a heuristic smallness choice, not a certified lower bound on δ.**

---

## §4. Execution and result — the honest record

**The preregistered Monte Carlo design failed to detect the signal.** The
run completed (see `v3_output.log`) but the result is uninformative:

| N_B | κ₃(P0) | γ₁(P0) | κ₃(P1) | γ₁(P1) | Var(P0) | Var(P1) |
|---|---|---|---|---|---|---|
| 4 | 4.08e-03 | 1.29e-02 | −7.13e-03 | −2.24e-02 | 0.4640 | 0.4661 |
| 8 | −4.75e-03 | −1.46e-02 | 8.59e-03 | 2.72e-02 | 0.4727 | 0.4642 |
| 16 | −5.12e-03 | −1.61e-02 | −1.78e-03 | −5.49e-03 | 0.4665 | 0.4723 |
| 32 | −3.86e-03 | −1.22e-02 | 3.01e-03 | 9.41e-03 | 0.4636 | 0.4674 |
| 64 | 7.96e-03 | 2.51e-02 | 1.25e-02 | 3.86e-02 | 0.4644 | 0.4727 |

**Diagnosis — signal buried in noise, by construction:**

- Sampling noise floor on κ₃ with N_samples = 20,000 is ≈ 1.6 × 10⁻²
  (measured: 30 repeated draws of a standard Gaussian at n = 20,000 give
  κ₃ std ≈ 0.0165; the same order applies here because Var(F) ≈ 0.46 for
  all cases).
- The **true** signal at t_num = 0.05 is
  κ₃(F_{P1,N}(0.05)) ≈ 3c(0.05)/N_B ≈ −Var(x₀²)·(0.05)⁷/(28π³)·3/N_B ≈
  **10⁻¹⁰–10⁻¹¹** for N_B ∈ {4,…,64}.
- Signal-to-noise ratio ≈ 10⁻⁸–10⁻⁹. **Reaching SNR of order unity at
  t_num = 0.05 via Monte Carlo would require roughly 10²⁰–10²² samples,
  depending on where the true signal lies in the quoted 10⁻¹⁰–10⁻¹¹ range
  (order-of-magnitude estimate from the measured noise floor and signal
  range; not a precise power calculation). Not feasible.**

**This is a design error in the preregistered path, not a failure of the
theorem.** The theorem's small-time regime (t → 0) is precisely where the
third-cumulant signal is smallest — t⁷ suppression. A direct-MC
illustration of the theorem's *quantified interval* (0, δ) is
computationally out of reach for any reasonable t in that interval.

**The sample-count reduction from 200,000 to 20,000** (mid-run, documented
here) is immaterial to the conclusion: even at 200,000 the noise floor
(≈ 5 × 10⁻³) exceeds the signal by ≥ 7 orders of magnitude.

---

## §5. What the numerical run does show (null controls)

Per charter §9:

- **Null 1 (P0 symmetric):** κ₃(P0) values (|γ₁(P0)| ∈ [1.2, 2.5] × 10⁻²)
  are consistent with the measured Gaussian noise floor (≈ 1.6 × 10⁻²).
  **P0 is numerically symmetric within sampling uncertainty.** ✓
- **Null 2 (P1 shape difference decreases with N_B):** not assessable —
  the signal is invisible at every N_B in the grid. ✗ (blocked by §4)
- **Null 3 (variance positive/stable):** Var(P0) ≈ Var(P1) ≈ 0.46–0.47 for
  all N_B; standardization well-defined throughout. ✓

No evidence of numerical pathology in the integration itself (RK4 at
dt = 1/512; variance stable; no NaN/inf).

---

## §6. Charter §5–§7 figure/visualization requirements

**Not deliverable at this boundary.** The primary statistic (|γ₁(P1)| vs
1/N_B) is below the Monte Carlo noise floor for every N_B in the
preregistered grid; no meaningful figure can be produced from the
preregistered design. The quotient-visualization (§6) and affine-fit check
(§7) are likewise blocked: they require the standardized P1 distribution's
shape to be resolvable above noise, which it is not at t_num = 0.05 with
feasible sample counts.

**Alternative path (documented, not executed):** a deterministic Gauss–
Hermite (or validated trapezoid) quadrature over (x₀, p₀), following the
original preflight's approach, could evaluate c(t*) at t* chosen large
enough for the signal to be detectable (e.g., t* ~ 0.5–1, still small
compared to π), and directly visualize the 1/N_B scaling of κ₃(F)/N_B.
This is the correct numerical route but requires either (i) the
certified-quadrature machinery the original record documented as not yet
built, or (ii) an owner decision to accept quadrature evidence at
evidence-grade (not certified) for the figure.

---

## §12. V3 disposition

# **V3-ILLUSTRATION-OPEN**

The preregistered Monte Carlo design is computationally inadequate for the
theorem's small-time regime (signal ≈ 10⁻¹¹ below noise floor ≈ 10⁻²).
No figure is produced. The theorem itself is untouched — V1 (the proof)
stands; this is a visualization-infrastructure failure, not a scientific
one.

Per charter §10: failure of the numerical illustration triggers
investigation (this report is that investigation) but does not rewrite the
theorem. Per charter §12: V3-ILLUSTRATION-OPEN is recorded; V3 does not
upgrade novelty or theorem status.

**Path forward (owner decision required):** either (a) authorize the
quadrature-based route (evidence-grade, consistent with the original
preflight's PF4Q-I grade discipline), or (b) accept that the figure lane
is not deliverable at this boundary and record V3 as open with no figure,
letting the publication case rest on V1 + V2 alone (the theorem and
novelty positioning survive without a numerical illustration).

---

## Files

- `v3_fast.py` — the executed script (preregistered t_num, grid, seed;
  reduced sample count documented).
- `v3_output.log` — the raw output, unmodified.
- No figure is produced; no results CSV (no signal to record).

**V1 is the proof; V3 is visualization/support only. Nothing in this
report upgrades or downgrades the theorem, V1, V2, or the
KNOWN-RESULT-NEW-FRAMING disposition.**
