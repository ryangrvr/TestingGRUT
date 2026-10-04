# VER0-B — Blind reproduction of BRI1-X1 (clamped Duffing bath, E₂± non-membership)

**Input:** `verification_0/specs/VER0_B_BRI1_X1_TARGET.md` only (see §11 for the full file-access log).
**Verdict:** the target theorem is **established**, and in fact holds in a **stronger form** (no large-N threshold is
needed: for P1 and all sufficiently small t the family is outside E₂± for **every** N_B ≥ 1). The asymptotic structure
(κ₃(F) = K(t)/N_B + O(N_B⁻²), first non-zero small-time order t⁷, negative sign, Gibbs combination Var_μ(x²)) is
reproduced. Two wording corrections are recorded in §10.

---

## 0. Setting, notation, standing hypotheses

* State z = (x, p) ∈ ℝ². Gibbs measure μ(dz) = Z⁻¹ exp[−(p²/2 + x²/2 + x⁴/4)] dz. Under μ, p ~ N(0,1) independent of x;
  x has density ∝ exp(−x²/2 − x⁴/4). Write m_k := E_μ[x^k] (m_odd = 0).
* **Moment identity (integration by parts).** For polynomial f, E_μ[f′(x)] = E_μ[f(x)(x + x³)] (boundary terms vanish).
  With f = x^{n−1}: (n−1) m_{n−2} = m_n + m_{n+2}. Hence m₄ = 1 − m₂, m₆ = 4m₂ − 1, m₈ = 6 − 9m₂, …
* Clamp q ∈ 𝒳, Q := sup|q|, T = 2π. For ε ∈ ℝ, X^ε_q(t)(z) is the x-component of the solution of
  ẋ = v, v̇ = −x − x³ + ε q(t), (x, v)(0) = z. Single-oscillator moments M_k(ε) := E_μ[X^ε_q(t)^k], cumulants κ_n(ε).
* Bath: X_j^ε(t), j = 1…N_B, are i.i.d. copies (same q, independent μ-distributed initial data; under the clamp the bath
  oscillators do not interact). F_q(t) = ε Σ_j X_j^ε(t) with ε = N_B^{−1/2}. Standardised skewness
  γ(Y) := κ₃(Y)/κ₂(Y)^{3/2} (for Var Y > 0, E|Y|³ < ∞).
* **Dominating class 𝒟.** 𝒟 := measurable functions D(z) with |D(z)| ≤ C (1 + |z|)^m exp[c(|x| + |p| + x²)] for some
  constants C, m, c. Every D ∈ 𝒟 is μ-integrable (the x-weight e^{−x⁴/4} beats e^{c x²}, the p-weight e^{−p²/2} beats
  e^{c|p|}), and 𝒟 is closed under sums and products. All constants below are uniform in |ε| ≤ 1 and s ∈ [0, T] unless
  stated, and depend on (q, T) only through Q.

KNOWN results used: smooth dependence of ODE solutions on initial data and parameters for a vector field continuous in
t and C^∞ in (state, parameter) (derivatives obey the variational equations); Liouville's theorem; Gronwall's
inequality; differentiation under the integral sign with an integrable dominating function (mean-value theorem +
dominated convergence); Taylor's theorem with Lagrange remainder; Lindeberg–Feller CLT for triangular arrays (via
Lyapunov's condition); Cramér–Wold; additivity and homogeneity of cumulants.

---

## 1. Basic lemmas (regularity and integrable bounds)

**Lemma A (global flow, a-priori bound).** Let E(t) = v²/2 + x²/2 + x⁴/4. Then Ė = ε q v, |Ė| ≤ |ε| Q √(2E), so
d/dt √(E+η) ≤ |ε|Q/√2 for every η > 0; letting η ↓ 0, √E(t) ≤ √E(0) + |ε| Q t/√2 (same backwards in time). Solutions
exist globally, the time-t map Φ^ε_t is a homeomorphism (indeed C^∞ diffeomorphism) of ℝ², and for |ε| ≤ 1, s ≤ T:

  |x(s)|, |v(s)| ≤ B(z) := √(2E₀) + Q T,  √(2E₀) = √(p² + x² + x⁴/2) ≤ |p| + |x| + x²/√2,

so B and every polynomial in B lie in 𝒟.

**Lemma C (energy estimate for the linearised operator).** Let a ∈ C¹, a ≥ 1, and Z″ + a Z = f, Z(0) = Z′(0) = 0. Put
W = Z′²/2 + a Z²/2. Then W′ = Z′ f + a′Z²/2 ≤ |f|√(2W) + (|a′|/a) W, hence (via √(W+η), η ↓ 0)

  √W(t) ≤ e^{Λ(t)} ∫₀ᵗ |f(s)| ds /√2,  Λ(t) := ∫₀ᵗ |a′|/(2a) ds,  and |Z(t)|, |Z′(t)| ≤ √(2W(t)).

For a = 1 + 3x(s)² along any trajectory of Lemma A: |a′|/a = 6|x||v|/(1+3x²) ≤ √3 |v| (since |x|/(1+3x²) ≤ 1/(2√3)),
so Λ(t) ≤ (√3/2) T B(z) and **e^{Λ} ∈ 𝒟**. (A naive Gronwall bound with sup a ~ B² would give e^{c B² t}, which is *not*
μ-integrable for large t; the energy weighting is what makes the interchange legitimate for every finite t.)

**Lemma B (ε-derivatives).** The vector field (v, −x − x³ + εq(t)) is continuous in t and polynomial in (x, v, ε), so
X^ε_q(t)(z) is C^∞ in (z, ε) (KNOWN). With a := 1 + 3X², Y_k := ∂_ε^k X^ε (all zero initial data):

  Y₁″ + a Y₁ = q,  Y₂″ + a Y₂ = −6 X Y₁²,  Y₃″ + a Y₃ = −6 Y₁³ − 18 X Y₁ Y₂.

(Differentiate X″ + X + X³ − εq = 0 three times.) By Lemma C, for |ε| ≤ 1, s ≤ T:
|Y₁(s)| ≤ e^{Λ} Q s, then |Y₂| ≤ e^{Λ}·6B·(e^{Λ}QT)²·T, etc.; every |Y_k|, |Y_k′| (k ≤ 3) is bounded by an element
of 𝒟. **Small-time refinement for P1** (used in §4–§5): on [0, π], 0 ≤ q(s) ≤ 10 (s/π)³ (because
s(u) = u³(10 − 15u + 6u²) and 0 < 10 − 15u + 6u² ≤ 10 on [0,1]); hence, for s ≤ π and |ε| ≤ 1,

  |Y₁(s)| ≤ D s⁴, |Y₂(s)| ≤ D s⁹, |Y₃(s)| ≤ D s¹³ (one D ∈ 𝒟),

since ∫₀ˢ|q| ≤ (5/(2π³)) s⁴, the Y₂-forcing is ≤ D s⁸ and the Y₃-forcing ≤ D (s¹² + s¹³).

**Corollary (moments are C³ in ε).** d^j/dε^j [X^k] (j ≤ 3, k ≤ 3) is a polynomial in X, Y₁, Y₂, Y₃ dominated, uniformly
in |ε| ≤ 1, by an element of 𝒟. By differentiation under the integral sign, M_k ∈ C³([−1, 1]) with
M_k^{(j)}(ε) = E_μ[∂_ε^j X^k]. Cumulants κ_n (n ≤ 3) are polynomials in M₁, M₂, M₃, hence C³ on [−1, 1].

**Lemma G (Gibbs invariance).** For ε = 0 (or q ≡ 0) the flow is Hamiltonian: it preserves Lebesgue measure (Liouville)
and the energy, hence μ. So X⁰(t) has the law of x under μ for every t; in particular E[X⁰(t)²] = m₂ for all t.

---

## 2. R1 — first-order response and parity structure

**Variation equation at ε = 0.** Y := ∂_ε X^ε|_{ε=0} solves Y″ + (1 + 3 X⁰(t)²) Y = q(t), Y(0) = Y′(0) = 0, i.e.
Y_q(t) = ∫₀ᵗ G(t, s; z) q(s) ds with G the Green function of the linearisation along the free trajectory X⁰.
Regularity and interchange with E_μ: Lemma B + Corollary.

**Parity (Z₂) structure.** Let S(x, v) = (−x, −v). If (x, v) solves the ε-equation then (−x, −v) solves the (−ε)-equation:
(−x)″ + (−x) + (−x)³ = −εq = (−ε)q. By uniqueness, X^{−ε}_q(t)(Sz) = −X^ε_q(t)(z) for all t (also jointly in t), and
μ∘S = μ. Hence, **as processes**, X^{−ε}_q(·) =ᵈ −X^ε_q(·), equivalently X^ε_{−q} =ᵈ −X^ε_q and F_{−q} =ᵈ −F_q. Therefore
κ_n(−ε) = (−1)ⁿ κ_n(ε): **κ₁, κ₃ are odd in ε, κ₂ is even.**

*What parity kills:* κ₃(0) = 0 and κ₃″(0) = 0 (and all even ε-orders of κ₃); κ₂′(0) = 0; M₁(0) = M₁″(0) = 0.

**Third-cumulant expansion.** κ₃ = M₃ − 3M₂M₁ + 2M₁³, so using M₁(0) = 0, M₃(0) = 0 and M₂(0) = m₂ (Lemma G):

  K_q(t) := κ₃′(0) = M₃′(0) − 3 M₂(0) M₁′(0) = **3 Cov_μ( X⁰(t)², Y_q(t) )**,

which is linear in q. Taylor with Lagrange remainder on the odd C³ function κ₃:

  **κ₃(X^ε_q(t)) = ε K_q(t) + R_q(ε, t),  |R_q| ≤ (|ε|³/6) sup_{|η|≤1} |κ₃‴(η)| =: C₃(t) |ε|³,  |ε| ≤ 1.**

**Uniform small-time form (P1, t ≤ π).** κ₃‴ = M₃‴ − 3Σ_{i+j=3} C(3,i) M₂^{(i)} M₁^{(j)} + 2(M₁³)‴, with
M₃‴ = E[6Y₁³ + 18XY₁Y₂ + 3X²Y₃], M₂‴ = E[6Y₁Y₂ + 2XY₃], M₂″ = E[2Y₁² + 2XY₂], M₂′ = E[2XY₁], M₁^{(j)} = E[Y_j] and
|M₁(ε)| ≤ |ε| sup|M₁′|. With the s⁴/s⁹/s¹³ bounds of Lemma B, every term is O(t¹²) uniformly in |ε| ≤ 1:
e.g. M₃‴ = O(t¹²), M₂″M₁′ = O(t⁸·t⁴), M₂′M₁″ = O(t⁴·t⁹), M₂M₁‴ = O(t¹³), (M₁³)‴ = O(t¹²). Hence

  **|κ₃(X^ε_{P1}(t)) − ε K_{P1}(t)| ≤ C* |ε|³ t¹²  for all |ε| ≤ 1, 0 < t ≤ π.**  (2.1)

---

## 3. R2 — cumulant scaling (exponent derived, not assumed)

X_j^ε i.i.d., so by additivity and homogeneity, κ_n(F_q(t)) = N_B εⁿ κ_n(X^ε_q(t)) = N_B^{1 − n/2} κ_n(ε). For n = 3:

  κ₃(F_q(t)) = N_B^{−1/2} κ₃(ε) = N_B^{−1/2}[0 + εK_q(t) + 0·ε² + O(ε³)] = **K_q(t)·N_B^{−1} + O(N_B^{−2})**,

with |remainder| ≤ C₃(t) N_B⁻². Parity removes the a-priori N_B^{−1/2} term (κ₃(0) = 0) and the N_B^{−3/2} term
(κ₃″(0) = 0). The exponent of the leading term is exactly −1 whenever K_q(t) ≠ 0 (shown for P1 below); the remainder
exponent is −2 (generically attained via κ₃‴(0)/6 ≠ 0; we only claim the bound). For n = 2:
Var F_q(t) = κ₂(ε) = m₂ + O(ε²) = m₂ + O(N_B⁻¹) (κ₂ even, C², κ₂(0) = m₂ by Lemma G). Hence for fixed t

  γ(F_q(t)) = K_q(t) / (m₂^{3/2} N_B) + O(N_B⁻²).

---

## 4. R3 — small-time coefficient for P1 (analytic proof)

Fix ε = 0, 0 < t ≤ π, where q(s) = 10(s/π)³ − 15(s/π)⁴ + 6(s/π)⁵. Write B₀ := √(2E₀) ∈ 𝒟, a(s) = 1 + 3X⁰(s)²,
and D₁, D₂, … ∈ 𝒟 for the dominating functions produced (all independent of t ∈ (0, π]).

1. Variation of constants: Y(t) = ∫₀ᵗ (t−s)[q(s) − a(s)Y(s)] ds. Exactly,
   Y_free(t) := ∫₀ᵗ (t−s) q(s) ds = t⁵/(2π³) − t⁶/(2π⁴) + t⁷/(7π⁵).
2. Lemma C: |Y(s)| ≤ D₁ s⁴. Then |Y − Y_free| ≤ (1+3B₀²) D₁ t⁶/30, so |Y(s)| ≤ D₂ s⁵ and
   Y(s) = s⁵/(2π³) + σ(s), |σ(s)| ≤ D₃ s⁶.
3. |X⁰(s)² − x²| ≤ 2 B₀² s (mean value; |d(X²)/ds| = 2|X v| ≤ 2B₀²), so a(s) = 1 + 3x² + ρ(s), |ρ(s)| ≤ 6B₀² s. Then
   ∫₀ᵗ (t−s) a(s) Y(s) ds = (1+3x²)(2π³)⁻¹ ∫₀ᵗ (t−s) s⁵ ds + O(D₄ t⁸) = (1+3x²) t⁷/(84π³) + O(D₄ t⁸).
4. Hence **Y(t) = α(t) − x² t⁷/(28π³) + R_Y(t)**, α(t) = t⁵/(2π³) − t⁶/(2π⁴) + t⁷/(7π⁵) − t⁷/(84π³) deterministic,
   |R_Y(t)| ≤ D₄ t⁸.
5. K(t) = 3 Cov_μ(X⁰(t)², Y(t)). Deterministic parts have zero covariance, so
   K(t) = −(3t⁷/(28π³)) Cov(X⁰(t)², x²) + 3 Cov(X⁰(t)², R_Y(t)), with
   |Cov(X⁰(t)², R_Y)| ≤ t⁸ E_μ[(B₀² + m₂)D₄] and Cov(X⁰(t)², x²) = Var_μ(x²) + Cov(X⁰(t)² − x², x²),
   |Cov(X⁰(t)² − x², x²)| ≤ 2t E_μ[B₀²(x² + m₂)]. All expectations finite (𝒟).

**Result.** For 0 < t ≤ π:

  **K_{P1}(t) = c₇ t⁷ + r(t),  c₇ = −(3/(28π³)) Var_μ(x²) = −(3/(28π³))(m₄ − m₂²) = (3/(28π³))(m₂² + m₂ − 1),
  |r(t)| ≤ C_K t⁸.**

* **Vanishing orders:** the coefficients of t⁰,…,t⁶ are zero. Y = O(t⁵) kills t⁰…t⁴; the t⁵ and t⁶ coefficients of Y
  are deterministic (they come only from ∫(t−s)q with the s³, s⁴ terms of q), and a deterministic summand has zero
  covariance with X⁰(t)². The first random part of Y (−x²t⁷/(28π³), from the anharmonic stiffness 3X⁰² acting on the
  leading response s⁵/(2π³)) produces the first non-zero order **t⁷**.
* **Sign: negative.** **Gibbs-moment combination:** Var_μ(x²) = m₄ − m₂² = 1 − m₂ − m₂² (using m₄ = 1 − m₂).
* **Strictly non-zero:** x² is a non-constant continuous function and μ has an everywhere-positive density, so x² is not
  μ-a.s. constant; hence Var_μ(x²) > 0 and c₇ < 0. (By-product: m₂ < (√5 − 1)/2.)

**Fresh symbolic cross-check (not the proof).** `code/ver0_b_small_time_series.py` (log `.log`) computes the exact
Taylor coefficients of X⁰, Y (truncated Picard iteration on polynomials — reproduces Taylor coefficients exactly up to
the truncation order) and reduces Gibbs moments with the IBP recursion. Output: E[X⁰(t)²] ≡ m₂, E[X⁰(t)³] ≡ 0;
Y-coefficients t⁵: 1/(2π³), t⁶: −1/(2π⁴), t⁷: −(3π²x² − 12 + π²)/(84π⁵) (agrees with step 4);
K: zero through t⁶, t⁷: 3(m₂² + m₂ − 1)/(28π³), t⁸: −9(m₂² + m₂ − 1)/(112π⁴), plus t⁹…t¹¹; and
c₇ − [−3/(28π³)(m₄ − m₂²)] = 0. The analytic remainder bound is steps 1–5 above (an alternative analytic route: on [0, π]
the augmented system (X⁰, V⁰, Y, Y′, τ) is autonomous polynomial, so the 8th t-derivative of X⁰²Y is a polynomial in the
state, dominated by 𝒟 via Lemmas A, C; Lagrange's remainder then gives O(t⁸) after E_μ).

**Labelled deterministic numerical sanity check (not proof; one pre-chosen time t* = 0.5, no scan).**
`code/ver0_b_numeric_sanity.py` (log `.log`), deterministic quadrature: m₂ ≈ 0.4679199170, Var_μ(x²) ≈ 0.3131310343,
c₇ ≈ −1.0820×10⁻³; K(0.5) from the variational ODE ≈ −6.8347×10⁻⁶ vs series through t¹¹ ≈ −6.8412×10⁻⁶;
κ₃(X^ε(0.5))/ε equals K(0.5) to ~10⁻¹⁴ for ε = 0.4, 0.2, 0.1 (consistent with (2.1): ε²t¹² is tiny); P0: κ₃ = 0.

---

## 5. R4 — sign interval (no scan)

From §4, for 0 < t ≤ π: K(t) ≤ c₇t⁷ + C_K t⁸. Set **δ := min(π, |c₇|/(2C_K)) > 0**. For t ∈ (0, δ):
K(t) ≤ c₇t⁷/2 < 0. Fixed sign: negative. No numerical δ is claimed.

Strengthened interval (used in §8): with C* from (2.1), set **δ* := min(δ, (|c₇|/(4C*))^{1/5})**. For t ∈ (0, δ*) and
every 0 < ε ≤ 1: κ₃(X^ε_{P1}(t)) ≤ ε[c₇t⁷/2 + C* t¹²] ≤ ε c₇ t⁷/4 < 0.

---

## 6. R5 — P0 control

For q ≡ 0 each X_j(t) is a free trajectory; by Lemma G, X_j(t) ~ law of x under μ, which is symmetric. The X_j are
independent, so F_P0(t) = ε Σ_j X_j(t) is symmetric; all odd cumulants vanish: **κ₃(F_P0(t)) = 0 and γ(F_P0(t)) = 0 for
every finite N_B and every t** (moments finite by Lemma A). Symmetry: the Z₂ parity (x, p) ↦ (−x, −p) of the bath
Hamiltonian H_B = Σ[p²/2 + x²/2 + x⁴/4] and of μ, which under the clamp is a symmetry of the dynamics only jointly with
q ↦ −q (ε ↦ −ε); P0 is its fixed point. (Also: F_{−q} =ᵈ −F_q for every q, so γ(F_{−q}) = −γ(F_q) — which is why the
witness must be |γ|, not γ, against a class that allows sign changes of G.)

---

## 7. R6 — non-degenerate variance

Var F_q(t) = N_B ε² κ₂(ε) = Var(X^ε_q(t)). The map z ↦ X^ε_q(t)(z) = π₁∘Φ^ε_t(z) is continuous, and surjective onto ℝ
because Φ^ε_t is a bijection of ℝ² (Lemma A, forward and backward global existence). If Var = 0 then X = c μ-a.s.;
but {z : X(z) ≠ c} is open and non-empty, so it has positive Lebesgue measure, hence positive μ-measure (μ ~ Lebesgue).
Contradiction. Therefore **Var F_q(t) > 0 for every q ∈ 𝒳, every t > 0 (indeed t ≥ 0) and every N_B ≥ 1**; moreover
Var F_q(t) = m₂ + O(N_B⁻¹) → m₂ > 0. All moments of F_q(t) are finite (Lemma A: |F| ≤ ε Σ_j B(z_j)).

---

## 8. R7 — E₂± orbit escape

**Signed-affine invariance.** If Y has finite third moment, Var Y > 0, α ∈ ℝ, λ ∈ ℝ∖{0}, then κ₂(α+λY) = λ²κ₂(Y),
κ₃(α+λY) = λ³κ₃(Y), so γ(α+λY) = sign(λ) γ(Y). Hence |γ| is invariant and in particular γ = 0 ⇔ γ ≠ 0 is preserved.

**Reduction.** Suppose, for a fixed N_B, the family {F_q : q ∈ 𝒳} were in E₂±: F_q(t) =ᵈ M_t[q] + G_t[q] ξ(t) with M, G
deterministic causal, G_t[q] ≠ 0, ξ of one q-independent law. (Equality of process laws implies equality of the time-t
marginals; a pathwise representation implies equality in law. Only the time-t marginal is used, so the argument covers
every reading.) Apply at P0, P1 ∈ 𝒳 (P1 is C² on [0, 2π]: s(1) = 1, s′(1) = s″(1) = 0; q(0) = q̇(0) = 0; bounded) and a
fixed t: with a_i = M_t[P_i], b_i = G_t[P_i] ≠ 0,
ξ(t) =ᵈ (F_P0(t) − a₀)/b₀ (so ξ(t) has all moments and positive variance by §7), and
F_P1(t) =ᵈ a₁ + (b₁/b₀)(F_P0(t) − a₀), a signed-affine image with slope b₁/b₀ ≠ 0. Then
|γ(F_P1(t))| = |γ(F_P0(t))| = 0 by §6.

**Contradiction.**
* Target form: fix t ∈ (0, δ). By §3, γ(F_P1(t)) = K(t)/(m₂^{3/2}N_B) + O(N_B⁻²) with K(t) < 0, so γ(F_P1(t)) < 0 for
  N_B ≥ N₀(t) := ⌊C₃(t)/|K(t)|⌋ + 1 (then |remainder of κ₃(F)| ≤ C₃/N_B² < |K|/N_B, and κ₂(F) > 0 by §7).
* **Stronger form:** fix t ∈ (0, δ*). By §5, κ₃(F_P1(t)) = N_B^{−1/2} κ₃(X^ε(t)) < 0 for **every** N_B ≥ 1
  (ε = N_B^{−1/2} ∈ (0, 1]), and κ₂ > 0, so γ(F_P1(t)) < 0 for every N_B ≥ 1.

Either way, the P1 time-t law is **not in the signed-affine orbit** of the P0 time-t law, so the family ∉ E₂±. Since
E₁ ⊂ E₂± (take G ≡ 1), the family is also ∉ E₁. Remark: the argument does not even use G ≠ 0 at P0 (b₀ = 0 would force
Var F_P0 = 0, contradicting §7) and b₁ = 0 would force Var F_P1 = 0; so the family is also outside the relaxed class
with zero allowed.

---

## 9. R8 — reservoir limit (triangular array), scoped

Fix q ∈ 𝒳 (in particular P0 or P1) and times t₁,…,t_k ∈ [0, T]; λ ∈ ℝᵏ. For N = N_B, ε_N = N^{−1/2}, put
S_ε := Σ_a λ_a X^ε_q(t_a) and row variables ζ_{N,j} := N^{−1/2}(S^{(j)}_{ε_N} − E S_{ε_N}), j ≤ N (i.i.d. within a row;
the law changes with N — a genuine triangular array).

* **Covariance convergence.** Σ_j Var ζ_{N,j} = Var S_{ε_N} → λᵀC₀λ, C₀(t_a, t_b) := Cov_μ(X⁰(t_a), X⁰(t_b)): X^ε → X⁰
  pointwise as ε → 0 (continuous dependence) with |X^ε X^ε'| ≤ B² ∈ 𝒟, so dominated convergence applies. C₀ is the
  equilibrium autocovariance of the free anharmonic oscillator, **independent of q**; C₀(t,t) = m₂ > 0.
* **Uniform moment / Lindeberg.** E S_ε⁴ ≤ (Σ|λ_a|)⁴ E_μ B⁴ =: L < ∞ uniformly in |ε| ≤ 1, so
  Σ_j E ζ⁴_{N,j} = N⁻¹ E(S − ES)⁴ ≤ 16L/N → 0 (Lyapunov ⇒ Lindeberg).
* **fdd Gaussian convergence.** Lindeberg–Feller gives Σ_j ζ_{N,j} ⇒ N(0, λᵀC₀λ) (if λᵀC₀λ = 0, Chebyshev gives ⇒ 0).
  Cramér–Wold: the centred fdd of F_q converge to those of a centred Gaussian process ξ with covariance C₀ — **the same
  law for every q** (here F_q − E F_q = Σ_j ζ_{N,j} with the corresponding λ).
* **Mean.** E F_q(t) = N^{1/2} M₁(ε_N) = M₁′(0) + O(ε_N²) = E_μ[Y_q(t)] + O(N⁻¹) (M₁ odd, C³), and
  m_q(t) := E_μ Y_q(t) = ∫₀ᵗ χ(t,s) q(s) ds, χ(t,s) = E_μ G(t,s;z): linear and causal in q.

Hence for each fixed q the fdd of F_q converge to those of m_q + ξ: an **E₁-type limit** (M_t[q] = m_q(t), shared
Gaussian ξ). The non-affine witness vanishes: γ(F_P1(t)) ~ K(t)/(m₂^{3/2}N_B) → 0, exactly at rate N_B⁻¹ for t ∈ (0, δ).
**Scope:** all constants depend on (Q, T); convergence is pointwise in q (trivially uniform on the finite frozen family,
and uniform on {q : sup|q| ≤ Q}), but not uniform over 𝒳, which contains clamps with arbitrarily large amplitude
(e.g. q_N = N^{1/2}φ makes εq_N = φ an O(1) forcing whose skewness does not vanish). Note that the pointwise limits do
assemble into an E₁ family {m_q + ξ}; what is not claimed is a uniform-in-q approximation of the finite-N family by it.

---

## 10. R9 — exact claim; corrections

**Confirmed (with edits):** "A finite reciprocal anharmonic environment, under the clamp protocol, generates an
interventional reduced-force **family** outside the shared causal affine exogenous class E₂± (and hence outside E₁);
the distinguishing witness — the one-time standardised skewness, P1 vs P0 — is of exact order N_B⁻¹ for fixed small t
and vanishes in the reservoir limit, where the per-protocol fdd limits are of E₁ form."

Corrections / precisions:
1. **Understatement (strengthening):** the large-N threshold is unnecessary. There is δ* > 0 such that for every
   t ∈ (0, δ*) and **every N_B ≥ 1**, γ(F_P1(t)) < 0 = γ(F_P0(t)); so the family is outside E₂± for **all** N_B ≥ 1.
   The spec's statement (with N₀(t)) is a true weaker corollary.
2. "an interventional reduced-force **law**" should read "**family**": any single law F_q is trivially in E₂±
   (ξ := F_q), and §3 of the spec itself makes membership a family property.
3. "reciprocal" should be read as: F_q = −∂H_int/∂q is the back-reaction of a bath driven by the prescribed q (the q
   dynamics, V and M play no role).
4. Not claimed, and correctly so: no escape from E_univ (F_q = 𝔉_t[q, U] with U = bath initial state), no primitive
   randomness, no ontology/GRUT content. The reservoir companion is pointwise in q (see §9 scope).
5. The witness only uses one-time marginals; no Gaussianity, Markov or stationarity assumption on ξ is needed or used.

### Precise theorem proved

Let H, μ, 𝒳, P0, P1, F_q be as in the spec; m₂ = E_μ[x²]. Then:
(i) for every q ∈ 𝒳, t ∈ (0, T], N_B ≥ 1: Var F_q(t) > 0, F_{−q} =ᵈ −F_q, κ₃(F_q(t)) = K_q(t)/N_B + O_{q,t}(N_B⁻²),
K_q(t) = 3Cov_μ(X⁰(t)², Y_q(t)), Var F_q(t) = m₂ + O(N_B⁻¹); κ₃(F_P0(t)) = 0;
(ii) K_{P1}(t) = −(3/(28π³)) Var_μ(x²) t⁷ + O(t⁸) as t ↓ 0, with Var_μ(x²) = 1 − m₂ − m₂² > 0;
(iii) there is δ* > 0 such that for every t ∈ (0, δ*) and every N_B ≥ 1, γ(F_P1(t)) < 0, and for fixed such t,
γ(F_P1(t)) = K_{P1}(t)/(m₂^{3/2}N_B) + O(N_B⁻²);
(iv) consequently, for every N_B ≥ 1, {F_q : q ∈ 𝒳} ∉ E₂± (and ∉ E₁);
(v) for each fixed q ∈ 𝒳, the fdd of F_q converge as N_B → ∞ to those of m_q + ξ, ξ a centred Gaussian process with
the q-independent covariance C₀(s,t) = Cov_μ(X⁰(s), X⁰(t)) and m_q linear causal in q.

---

## 11. File-access log

Read/opened: `/home/user/VER0/verification_0/specs/VER0_B_BRI1_X1_TARGET.md` only, plus my own files below.
Disclosure: a `mkdir -p … && ls` on `/home/user/VER0/verification_0/code/` printed the *names* of pre-existing files in
that directory (`v0_1_*`, `v0_2_*`, `v0_3_*`, `v0_4_*`); none was opened or read, and none bears on this model.

Created: `verification_0/VER0_B_BRI1_X1_REPRODUCTION.md`, `verification_0/code/ver0_b_small_time_series.py` (+ `.log`),
`verification_0/code/ver0_b_numeric_sanity.py` (+ `.log`).
