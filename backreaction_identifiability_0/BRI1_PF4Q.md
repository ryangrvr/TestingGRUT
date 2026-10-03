# BRI1-PF4Q — FINITE-TIME COEFFICIENT EVALUATION AT THE FROZEN τ

**Status:** owner-authorised continuation of the Candidate-1 preflight. **Not** the Candidate-1 numerical campaign.
- No Monte Carlo or quasi-Monte Carlo.
- No finite-N_B simulation.
- No numerical D_orb.
- No change to P0 / P1 / P2, τ = (π, 3π/2, 2π) or the N_B grid.

**Base:** `grut-backreaction-identifiability-0 @ 6220a9e`. At that boundary the owner accepted Scope Repair 02, the
Candidate-1 preflight and the terminal **X1-PF-INDETERMINATE**.

**Sections 0 – 3 are written and committed BEFORE any frozen-τ coefficient is evaluated.**

## 0. Owner ruling recorded

- **Accepted:** Scope Repair 02, the Candidate-1 analytic preflight, and X1-PF-INDETERMINATE at `6220a9e`.
- **PF4Q-01:** grade the symmetric third-cumulant tensor
  K_q[a,b,c] = c_q(t_a,t_b;t_c) + c_q(t_a,t_c;t_b) + c_q(t_b,t_c;t_a),
  where c_q(t_a,t_b;t_c) = Cov_Gibbs(x₀(t_a)x₀(t_b), y₁,q(t_c)).
  - All **10** components (1 ≤ a ≤ b ≤ c ≤ 3) are evaluated for **both** P1 and P2: **20 K-values**.
  - One-time cases: K_q[a,a,a] = 3c_q(t_a,t_a;t_a).
  - The individual c-terms are not graded.
- **PF4Q-02 … 08:** as stated by the owner, and implemented below. In summary:
  - the minimal lemma BRI1-R1;
  - the small-time result is kept but is never a witness;
  - certified quadrature with a six-item error budget, **or** evidence only;
  - validations V1 – V3;
  - the grading rule PF4Q-A / I / Z / R;
  - no claim inflation;
  - "D_orb ~ O(N_B⁻²)" is a heuristic only.

## 1. LEMMA BRI1-R1 — first-order Gibbs response justified (PROVED HERE; INTERNALLY PROVED / NOT EXTERNALLY REVIEWED)

**Setting.** A single clamped oscillator z = (x, p):
- ẋ = p, ṗ = −x − x³ + ε q(t);
- q ∈ {q₁, q₂}, with |q| ≤ Q := 2 on [0, 2π];
- initial data z₀ ~ ρ(z₀) ∝ exp(−H₀(z₀)), where H₀ = p²/2 + x²/2 + x⁴/4;
- |ε| ≤ 1.

Write X^ε(t) := x(t; z₀, ε). Let τ = (t₁, t₂, t₃) = (π, 3π/2, 2π).

**LEMMA BRI1-R1.**
- **(i)** For each z₀ the solution exists on [0, 2π] and is C^∞ in (z₀, ε).
- **(ii)** y₁ := ∂_ε X^ε|_{ε=0} solves ÿ₁ + (1 + 3x₀²)y₁ = q(t), with y₁(0) = ẏ₁(0) = 0.
- **(iii)** Every moment E[Π_{i≤3} ∂_ε^{k_i} X^ε(t_{j_i})] with Σk_i ≤ 3 is finite, uniformly in |ε| ≤ 1, and
  differentiation in ε commutes with E up to third order.
- **(iv)** κ₃(X^ε(t_a), X^ε(t_b), X^ε(t_c)) is an **odd** C³ function of ε, with derivative K_q[a,b,c] at ε = 0. Hence,
  for every a, b, c ∈ {1, 2, 3}:
  **κ₃(F_q(t_a), F_q(t_b), F_q(t_c)) = K_q[a,b,c]/N_B + R, with |R| ≤ C_q·N_B⁻²**,
  where C_q < ∞ is uniform on the frozen tuple.

**Proof.**

*(i) Energy bound and smoothness.*
1. Along the forced flow, dH₀/dt = ε q p, so d√H₀/dt = ε q p / (2√H₀) ≤ |ε| Q/√2.
2. Hence √H₀(t) ≤ √H₀(0) + √2·Q·π =: √E₀ + c_Q for t ≤ 2π.
3. The orbit stays in a compact set, so the solution is global on [0, 2π].
4. The vector field is polynomial in (z, ε), so the flow is C^∞ in (z₀, ε) on the compact time interval (standard
   smooth dependence on initial data and parameters).

*(ii) Variational equation.* Differentiate the equation in ε at ε = 0. x₀ := X⁰ is the unforced Gibbs trajectory.

*(iii) Moment bounds.*
1. Let w_k := ∂_ε^k z. Then:
   - w₁ obeys ẇ₁ = A(t)w₁ + (0, q), with A = [[0, 1], [−1 − 3x², 0]];
   - w₂ obeys ẇ₂ = A w₂ + (0, −6x(w₁ˣ)²);
   - w₃ obeys ẇ₃ = A w₃ + (0, −18x w₁ˣ w₂ˣ − 6(w₁ˣ)³).
2. Each starts from 0.
3. Since x⁴/4 ≤ H₀, we have x² ≤ 2√H₀ ≤ 2(√E₀ + c_Q). So ‖A(t)‖ ≤ 2 + 6(√E₀ + c_Q) =: a(E₀), and |x| ≤ a(E₀).
4. Grönwall on [0, 2π] gives:
   - |w₁| ≤ 2πQ·e^{2πa};
   - |w₂| ≤ 2π·6a·|w₁|²_max·e^{2πa} ≤ C·a·e^{6πa};
   - |w₃| ≤ C′·a²·e^{10πa};
   - and |X^ε| ≤ a.
5. Every integrand in (iii) is therefore bounded by P(E₀)·exp(30π·6√E₀) for a polynomial P, uniformly in |ε| ≤ 1.
6. Under ρ, E₀ = H₀(z₀). Since ∫ e^{−E₀ + β√E₀} P(E₀) dx₀ dp₀ < ∞ for every β (because e^{−E₀ + β√E₀} ≤ e^{β²/2}·e^{−E₀/2},
   and the phase-space area of {H₀ ≤ e} grows polynomially in e), these bounds are ρ-integrable.
7. Dominated convergence then justifies differentiating under E up to order 3, with continuous derivatives.

*(iv) Parity and the cumulant expansion.*
1. The map (x, p, ε) ↦ (−x, −p, −ε) maps solutions to solutions, and ρ is invariant under (x, p) ↦ (−x, −p). Hence
   (X^{−ε}(t))_t has the same law as (−X^ε(t))_t, jointly in t.
2. So κ₃(X^{−ε}) = −κ₃(X^ε): κ₃ is odd in ε, and in particular κ₃(x₀) = 0 and E x₀(t) = 0.
3. Expand κ₃ = E[X_aX_bX_c] − E X_a·E[X_bX_c] − E X_b·E[X_aX_c] − E X_c·E[X_aX_b] + 2E X_a·E X_b·E X_c, and
   differentiate at ε = 0 using (iii):
   - d/dε E[X_aX_bX_c] = E[y_a x_b x_c] + E[x_a y_b x_c] + E[x_a x_b y_c];
   - in each of the three subtracted products, the factor E x₀ = 0 kills one term, leaving E y_c·E[x_a x_b] and its
     permutations;
   - the triple product's derivative vanishes, since every term keeps a factor E x₀ = 0.
4. So d/dε κ₃|₀ = Σ_{3 perms}(E[x x y] − E[x x]·E[y]) = **K_q[a,b,c]**.
5. κ₃(X^ε) is C³ and odd, so Taylor's theorem gives κ₃(X^ε) = εK + R₃(ε), with |R₃| ≤ (|ε|³/6)·sup|∂_ε³κ₃| < ∞ by (iii).
6. Under the clamp the N_B oscillators are i.i.d. copies of X^ε with ε = N_B^{−1/2}. F_q = εΣ_jX_j, so the exact
   identity (★) of the Candidate charter gives κ₃(F_q) = N_B·ε³·κ₃(X^ε) = N_B^{−1/2}κ₃(X^ε).
7. Substituting: κ₃(F_q) = K/N_B + N_B^{−1/2}R₃(N_B^{−1/2}), with |N_B^{−1/2}R₃(N_B^{−1/2})| ≤ C_q·N_B^{−2}.

∎

**Scope.**
- C_q is finite but astronomically crude (Grönwall). It certifies the **order** of the remainder, not its size at
  small N_B.
- BRI1-R1 replaces the part of ASSUMPTION R needed for this gate (the third cumulant at τ, to first order with
  remainder). The all-orders ASSUMPTION R remains unproved and is not needed here.

## 2. Certification machinery: availability statement (made before any τ evaluation)

**Available here:**
- arb ball arithmetic (python-flint 0.9.0), including rigorous one-dimensional integration of analytic integrands
  (`acb.integral`).

**Not available here:**
- a validated ODE integrator: an interval / Taylor-model / Lohner-type method with wrapping control for box initial
  data;
- an error-certified two-dimensional quadrature for ODE-defined integrands.

A grade-bearing enclosure of K_q needs budget items 3 – 5 certified (ODE errors for x₀ and y₁, and quadrature error).
That would require building such a pipeline. **One viable route**, recorded for the owner and not built:
1. a complex-strip trapezoid rule in (x₀, p₀), using the analyticity of the Duffing flow for complex initial data within
   a uniform strip;
2. crude interval bounds on |integrand| along the strip boundary, where wrapping is tolerable because the error factor
   is e^{−2πa/h};
3. validated point integration at the nodes;
4. analytic Gibbs-tail bounds.

This is a substantial engineering task with many hours of compute, so it is outside this turn.

**Therefore, per PF4Q-04's fallback, every τ value computed below is NUMERICAL EVIDENCE ONLY. It is NOT certified and
CANNOT upgrade the terminal to X1-PF-A.**

Budget item status:

| item | status |
|---|---|
| 1. Gibbs normalisation | ratio estimator on the same nodes. The exact normaliser is certified with arb (V1) |
| 2. Tail truncation | Gauss–Hermite: none needed. Trapezoid: the domain is chosen so the omitted Gibbs mass is < 10⁻¹⁵ (analytic bound) |
| 3. ODE error (x₀) | estimated only (two tolerances) |
| 4. Variational / finite-difference error | estimated only |
| 5. Quadrature error | estimated only (resolution ladders) |
| 6. Rounding | float64, estimated |

## 3. Pre-declared evidence method (fixed before the τ run; not tuned on τ results)

**Two independent code paths:**

| | Method A — variational | Method B — finite difference |
|---|---|---|
| ODE system | (x, p, y₁⁽ᴾ¹⁾, v⁽ᴾ¹⁾, y₁⁽ᴾ²⁾, v⁽ᴾ²⁾) | the nonlinear clamped oscillator at ε = ±h_ε |
| how K is formed | from the c-terms | K ≈ κ₃(X^{h})/h (the O(h²) error follows from oddness), with Richardson extrapolation between h and h/2 |
| quadrature | tensor probabilists' Gauss–Hermite in (x₀, p₀), with factor e^{−x₀⁴/4} in the integrand; ratio-normalised weights | uniform trapezoid on [−4.5, 4.5] × [−8.5, 8.5] with factor e^{−H₀} |
| ODE solver | DOP853, rtol = atol = 10⁻¹³ | DOP853, rtol = atol = 10⁻¹² |
| restart | at t = π, so q's third-derivative jump is never stepped across | same |
| resolution ladder | n ∈ {48, 64, 96} nodes per dimension | grid spacing h ∈ {0.12, 0.08}; ε ∈ {10⁻³, 5·10⁻⁴} |

**Evidence criterion (declared now).** A component is recorded as "numerically non-zero (evidence)" only if both hold:
- A(n = 96) and B(finest, Richardson) agree to within 10⁻⁶ absolute;
- |K| > 100 × max(|A₉₆ − A₆₄|, |A₉₆ − B|).

**Validations** (run before the τ evaluation; they print no τ coefficient):
- **V1:** certified arb enclosures of Z_x, m₂ and m₄. Checks: m₂ + m₄ = 1 (enclosure contains 1) and
  Var(x₀²) = 1 − m₂ − m₂² > 0 (enclosure excludes 0). The quadrature ratio estimates of m₂ and m₄ are compared against
  the certified values.
- **V2:** P0 one-time variance E[x₀(t)²] at t ∈ {π, 3π/2, 2π}, compared with m₂.
- **V3:** at the **pre-declared validation time t_v = π/8** (code validation only, never a witness),
  c_P1(t_v, t_v; t_v) from A and B, compared with the exact series c = Σ_{k=7}^{14} c_k t^k from
  `bri1_preflight_symbolic.log`, evaluated at the certified m₂. The sign is expected negative.

**Grading under §2:** whatever the τ values, **no certified enclosure exists**. The terminal therefore cannot become
X1-PF-A in this turn (PF4Q-04 fallback). The applicable rule is **PF4Q-I**: BRI1-R1 is proved and finite-time
non-vanishing at τ is **not certified**. In this turn the reason is that certification machinery is absent, **not**
that the enclosures contain zero.

### 3.1 Validation-driven amendment (recorded BEFORE any τ evaluation; triggered only by V2 on P0)

**V2 failed for Method A as declared.**
- With tensor Gauss–Hermite at n = 64, the P0 one-time variance E[x₀(t)²] drifts from the certified m₂ by up to
  5·10⁻⁴ at t = 2π (`preflight/pf4q_validate.log`).
- The diagnostic (`preflight/pf4q_v2_diagnostic.log`, P0 only) shows phase mixing: x₀(t) becomes increasingly
  oscillatory in the initial data, so Gauss–Hermite converges slowly. At t = 2π the error is
  1.6·10⁻³ (n = 96), 2.0·10⁻⁴ (n = 128) and 6·10⁻⁸ (n = 192).
- The uniform trapezoid converges exponentially: 7·10⁻¹⁰ (h = 0.12), 1·10⁻¹³ (h = 0.08), ≤ 5·10⁻¹⁶ (h = 0.06).
- No τ coefficient was computed or inspected in this diagnosis.

**Amendment (supersedes §3 where they differ):**
- **Method A:** variational equations on the **uniform trapezoid** grid, h_A ∈ {0.08, 0.06}, same domain. Gauss–Hermite
  at n = 192 is reported as a **tertiary** check only.
- **Method B:** unchanged. Finite difference on the trapezoid with h_B ∈ {0.12, 0.08} and ε ∈ {10⁻³, 5·10⁻⁴},
  Richardson-extrapolated.
- **Independence:** A and B differ in formulation (variational vs finite difference in ε). The decisive comparison uses
  **different grids**: A at h = 0.06 against B at h = 0.08.
- **Evidence criterion (replaces §3):** a component is "numerically non-zero (evidence)" only if both hold:
  - |A₀.₀₆ − B_R(0.08)| < 10⁻⁶;
  - |K_A(0.06)| > 100 × max(|A₀.₀₆ − A₀.₀₈|, |A₀.₀₆ − B_R(0.08)|, |B_R(0.08) − B_R(0.12)|).

**Status of the validations:**
- V1 is certified and PASSES.
- V3 PASSES: series −4.4647·10⁻⁷; A −4.4647·10⁻⁷; B −4.4647·10⁻⁷; all negative.
- V2 PASSES for the amended Method A (trapezoid, h ≤ 0.08: |error| ≤ 10⁻¹³).

**Grading is unchanged:** evidence only; PF4Q-I.

## 4. Results — the 20 frozen K_q[a,b,c] (NUMERICAL EVIDENCE ONLY; NOT CERTIFIED)

**Run record** (`preflight/pf4q_run.py`, log `pf4q_run.log`, values `pf4q_results.json`):
- The first execution printed all 20 values, then crashed while writing JSON (a numpy-bool serialisation error), so
  `pf4q_run_attempt1_json_crash.log` is kept as emitted.
- After the one-line serialisation fix, the deterministic rerun reproduced **every printed value identically** (diff
  empty).

**Columns:**
- K_A: Method A (variational) on the trapezoid at h = 0.06;
- |A − B|: difference from Method B (finite difference, Richardson, trapezoid h = 0.08);
- max spread: the largest of |A₀.₀₆ − A₀.₀₈|, |A₀.₀₆ − B₀.₀₈| and |B₀.₀₈ − B₀.₁₂|.

| τ component | P1: K_A | P1: \|A − B\| | P1: max spread | P2: K_A | P2: \|A − B\| | P2: max spread |
|---|---|---|---|---|---|---|
| (π, π, π) | −0.14878490 | < 10⁻¹¹ | < 10⁻¹⁰ | −0.14878490 | < 10⁻¹¹ | < 10⁻¹⁰ |
| (π, π, 3π/2) | −0.04432639 | < 10⁻¹¹ | 9·10⁻¹² | −0.04446964 | 10⁻¹² | 8·10⁻¹² |
| (π, π, 2π) | −0.02861432 | < 10⁻¹¹ | 4·10⁻¹¹ | −0.11718637 | < 10⁻¹¹ | 1.5·10⁻⁹ |
| (π, 3π/2, 3π/2) | +0.06894372 | < 10⁻¹¹ | 9·10⁻¹¹ | +0.07089725 | < 10⁻¹¹ | 8·10⁻¹¹ |
| (π, 3π/2, 2π) | +0.08596046 | < 10⁻¹¹ | 5·10⁻¹⁰ | +0.14047054 | < 10⁻¹¹ | 2·10⁻⁹ |
| (π, 2π, 2π) | −0.12240520 | < 10⁻¹¹ | 1.2·10⁻⁹ | −0.12335677 | < 10⁻¹¹ | 1.2·10⁻⁸ |
| (3π/2, 3π/2, 3π/2) | −0.24181326 | < 10⁻¹¹ | 5·10⁻¹⁰ | −0.24953226 | < 10⁻¹¹ | 5·10⁻¹⁰ |
| (3π/2, 3π/2, 2π) | +0.12598512 | < 10⁻¹¹ | 2·10⁻⁹ | +0.09826085 | 10⁻¹² | 2.6·10⁻⁹ |
| (3π/2, 2π, 2π) | −0.09578872 | 3·10⁻¹² | 2.2·10⁻⁹ | −0.05959596 | 4·10⁻¹² | 1.5·10⁻⁸ |
| (2π, 2π, 2π) | +0.00565026 | 1.3·10⁻¹¹ | 8.7·10⁻⁸ | −0.14313464 | 2·10⁻¹¹ | 3.8·10⁻⁸ |

**Evidence criterion (§3.1): met by 20 / 20 components.**
- Methods A and B agree to ≤ 2·10⁻¹¹.
- The largest spread is 8.7·10⁻⁸, against |K| ≥ 5.65·10⁻³.
- The tertiary Gauss–Hermite n = 192 values agree to its own (lower) accuracy.

**Internal consistency:**
- **K(π, π, π) is identical for P1 and P2**, as causality requires, since q₁ = q₂ on [0, π].
- The P1 and P2 components involving later times differ.

**The smallest component**, P1 (2π, 2π, 2π) = +5.65·10⁻³:
- Its spread is the largest of all (8.7·10⁻⁸), and it still meets the criterion by a factor of about 65 000.
- The Gauss–Hermite n = 192 value (4.99·10⁻³) is visibly under-resolved there, consistent with the V2 diagnosis at
  t = 2π.

**Sign pattern** (all 20 components):

| components | P1 | P2 |
|---|---|---|
| (π,π,π), (π,π,3π/2), (π,π,2π), (π,2π,2π), (3π/2)³, (3π/2,2π,2π) | − | − |
| (π,3π/2,3π/2), (π,3π/2,2π), (3π/2,3π/2,2π) | + | + |
| (2π)³ | + | − |

The one-time components (K = 3c) are:
- K(π)³ = −0.1488 (both protocols);
- K(3π/2)³ = −0.2418 (P1) and −0.2495 (P2);
- K(2π)³ = +0.0057 (P1) and −0.1431 (P2).

**Leading-order reading (no claim inflation).** By BRI1-R1, κ₃(F_q(t_a), F_q(t_b), F_q(t_c)) = K_q[a,b,c]/N_B +
O(N_B⁻²), and every third cumulant of P0 vanishes exactly. If these values were certified, each would be a
leading-order, reflection-safe orbit violation (BRI-E2+O type) at the preregistered tuple. They do **not** establish a
finite-N_B detection threshold, do not validate D_orb, and do not require the N_B grid. They are not GRUT physics.

## 5. Grading (PF4Q-06)

**BRI1-R1:** PROVED (§1).

**Certified enclosures:** NONE. The machinery is unavailable (§2), so budget items 3 – 5 are only estimated.

**Per the owner's PF4Q-04 fallback, high-precision deterministic values are evidence only and may not upgrade PF-A.**

**Rule applied: PF4Q-I.** Retain **X1-PF-INDETERMINATE**, with the narrower statement: **"finite-time non-vanishing at τ
not certified."**
- **Precise reason:** no certified enclosure was computed (machinery absent). It is **not** that enclosures contain zero.
- **Strength of the evidence:** 20 / 20 components are non-zero, with two independent formulations and quadratures
  agreeing to about 10⁻¹¹ and magnitudes 5.7·10⁻³ – 0.25.
- The coefficients are **not** called zero. **PF-B is not assigned. PF4Q-Z does not arise.**

**Path to PF-A (owner decision):** build the certified pipeline sketched in §2 and certify at least one component. The
best candidate is P2 (3π/2)³ ≈ −0.250, the largest magnitude, which leaves the widest margin.

**D_orb (PF4Q-08):** "D_orb ~ O(N_B⁻²)" remains a feasibility heuristic. It is not derived here and not used.

## 6. Owner ruling on PF4Q and provenance wording repair (additive)

**Ruling.** BRI1-PF4Q is accepted as reported at `d5a0bdb`.
- **PF4Q-I, X1-PF-INDETERMINATE AT THE FROZEN τ**, for the stated reason: finite-time non-vanishing at τ was not
  certified.
- The 20 K-values remain **strong numerical evidence only, not certified**. They are not called zero, they do not upgrade
  PF-A, and the certification rule is not weakened.
- **The certified-τ pipeline is not to be built now.** It remains documented (§2) as a possible future verification
  project.

**Wording repair (provenance only).** Wherever §3 / §3.1 / §4 describe Methods A and B as "independent", read: **two
independent formulations / deterministic cross-checks with shared numerical infrastructure**.
- They are distinct in formulation: variational equations vs a finite difference in ε.
- They follow different declared resolution paths: trapezoid h = 0.06 vs 0.08 / 0.12.
- They share the DOP853 ODE solver, the trapezoid quadrature family and the core code (`pf4q_core.py`).
- They are not fully independent numerical reproductions.
