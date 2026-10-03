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
