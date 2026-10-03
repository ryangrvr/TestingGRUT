# BRI1 — ANALYTIC ESCAPE THEOREM FOR X1 (analytic only; no new computation)

**Status:** INTERNALLY PROVED / NOT EXTERNALLY REVIEWED.
- No numerical calculation was made for this document. No witness time is chosen numerically, no new protocol and no
  new candidate are introduced.
- The frozen-τ PF4Q values are **not** used anywhere in the proof.

**Base:** `grut-backreaction-identifiability-0 @ d5a0bdb` (BRI1-PF4Q accepted as PF4Q-I).

**Inputs used:**
- Candidate-1 definitions (`BRI1_CANDIDATE_CHARTER.md` §1);
- the exact identity (★): κ_n(F) = N_B^{1−n/2}·κ_n(X^ε);
- the parity facts (§4.2 there);
- **LEMMA BRI1-R1** (`BRI1_PF4Q.md` §1);
- the **exact symbolic Taylor coefficients** of `preflight/bri1_preflight_symbolic.log` (coefficients of t⁰ … t¹⁴ of
  c(t, t; t) for P1);
- the reflection-safe invariants of E₂± (`BRI0_CHARTER.md` §R BRI-O5, §R2).

**Notation:**
- H₀(x, p) = p²/2 + x²/2 + x⁴/4 and E₀ := H₀(z₀). The initial law is ρ ∝ e^{−H₀}.
- m₂ = E x₀² and Var(x₀²) = m₄ − m₂² = 1 − m₂ − m₂², using the exact identity m₂ + m₄ = 1.
- For P1 on [0, π]: q(t) = s(t/π), a polynomial.
- c(t) := c_{P1}(t, t; t) = Cov(x₀(t)², y₁(t)) = E[(x₀(t)² − m₂)·y₁(t)]. The last form uses E x₀(t)² = m₂, by stationarity
  of the Gibbs law under the unforced flow.

## T1 — LEMMA BRI1-T1 (finite-order small-time expansion)

**LEMMA BRI1-T1.** There exists δ ∈ (0, 1] such that **c(t) < 0 for every 0 < t < δ**.

**Proof.**
1. **Pathwise smoothness.** For fixed z₀, the unforced orbit x₀ solves a polynomial ODE, and y₁ solves the linear ODE
   ÿ₁ = −(1 + 3x₀²)y₁ + q(t) with polynomial q on [0, 1]. Both are C^∞ in t. Hence g(t; z₀) := (x₀(t)² − m₂)·y₁(t) is
   C^∞ on [0, 1].
2. **Derivative bounds on [0, 1].**
   - Energy is conserved along x₀, so x₀² ≤ 2√E₀ and p₀² ≤ 2E₀.
   - With k(t) = 1 + 3x₀² ≤ k̄ := 1 + 6√E₀, the first-order system for (y₁, ẏ₁) has matrix norm ≤ 1 + k̄. Grönwall gives
     |y₁|, |ẏ₁| ≤ e^{1+k̄} on [0, 1], because |q| ≤ 1 there.
   - Every time derivative of order ≤ 8 of x₀, ẋ₀, y₁ or ẏ₁ is, by repeated use of the ODEs, a polynomial in
     (x₀, p₀, y₁, ẏ₁) and in q, …, q⁽⁸⁾, the latter bounded on [0, 1]. So
     sup_{t∈[0,1]} |∂_t⁸ g(t; z₀)| ≤ P(E₀)·e^{1+k̄} =: D(z₀) for a polynomial P.
3. **Integrability.** Since −E₀ + 6√E₀ ≤ −E₀/2 + 18, and phase-space volume grows polynomially in E₀, D is ρ-integrable.
   Likewise |∂_t^k g(0; z₀)| is a polynomial in z₀ for k ≤ 7, hence integrable.
4. **Taylor with remainder, then expectation.** For each z₀,
   g(t) = Σ_{k=0}^{7} g⁽ᵏ⁾(0)t^k/k! + r(t; z₀), with |r| ≤ D(z₀)·t⁸/8!.
   Taking E (no derivative/expectation interchange is needed, because the identity holds pointwise and every term is
   integrable):
   **c(t) = Σ_{k=0}^{7} C_k t^k + R(t), with |R(t)| ≤ C_R·t⁸, where C_k := E[g⁽ᵏ⁾(0)]/k! and C_R := E D/8! < ∞.**
5. **The coefficients are exactly the symbolic ones.**
   - g⁽ᵏ⁾(0)/k! is the k-th Taylor coefficient of g at 0. It is the polynomial in (a, b) = z₀ produced by the exact
     recursion in `bri1_preflight_symbolic.py`, which builds the Taylor coefficients of x₀ and y₁ order by order from the
     ODEs. It is exact through t¹⁴, because each coefficient of order k − 2 of the right-hand side depends only on
     coefficients of order ≤ k − 2.
   - The script evaluates E of these polynomials exactly, using Gaussian moments of b and the a-moments reduced by
     m_{k+1} + m_{k+3} = k·m_{k−1}. That recursion is an integration by parts against e^{−a²/2 − a⁴/4}, whose boundary
     terms vanish.
   - Its output: **C₀ = … = C₆ = 0 and C₇ = (m₂² + m₂ − 1)/(28π³) = −Var(x₀²)/(28π³).**
   - The script computes E[x₀(t)²y₁] − E[x₀(t)²]·E[y₁]. This equals E[g] because E[x₀(t)²] = m₂ for all t, which is
     exact by invariance.
6. **Sign.** x₀² is not a.s. constant under ρ, since its law has a density, so **Var(x₀²) > 0** and C₇ < 0. (It is also
   certified numerically by arb in V1, Var(x₀²) ∈ 0.3131310343256930878… ± 4·10⁻²¹, but the proof does not need that.)
7. **Conclusion.** Choose δ := min(1, |C₇|/(2C_R)). For 0 < t < δ, c(t) ≤ C₇t⁷ + C_R t⁸ ≤ (C₇/2)·t⁷ < 0. ∎

**Remarks.**
- δ is existential; no value is claimed.
- No complex-time radius and no all-orders expansion is used.

## T2 — Non-zero third cumulant for all large finite N_B

Fix any t* ∈ (0, δ). By T1, **K_{P1}(t*, t*, t*) = 3c(t*) < 0.**

**BRI1-R1 at t*.** The proof of BRI1-R1 (`BRI1_PF4Q.md` §1) uses only:
- t ≤ 2π;
- the polynomial vector field;
- the Grönwall bounds;
- Gibbs integrability;
- parity.

It therefore applies verbatim to the one-time tuple (t*, t*, t*), giving

  κ₃(F_{P1,N}(t*)) = K_{P1}(t*, t*, t*)/N_B + R_N, with |R_N| ≤ C(t*)·N_B⁻², and C(t*) < ∞.

**Claim.** Let N₀(t*) := ⌊C(t*)/|K(t*)|⌋ + 1. Then for every integer N_B ≥ N₀:
- |R_N| ≤ C·N_B⁻² < |K|·N_B⁻¹;
- hence **κ₃(F_{P1,N}(t*)) ≠ 0**, with the sign of K, i.e. negative.

∎

N₀ is existential. **No value is claimed, and it is not claimed to lie inside the frozen grid {1, …, 64}.**

## T3 — Conversion to E₂± escape: THEOREM BRI1-X1

**P0 third cumulant.** F_{P0,N}(t) = N_B^{−1/2}·Σ_j x₀,j(t), with i.i.d. symmetric summands, since the Gibbs law and the
unforced flow are invariant under (x, p) ↦ (−x, −p). **So κ₃(F_{P0,N}(t)) = 0 exactly for every t and every N_B.**

**Positive variances.**
- Var F_{P0,N}(t) = m₂ > 0.
- Var F_{P1,N}(t*) = Var X^ε(t*), where ε = N_B^{−1/2}.
- The forced flow z₀ ↦ z(t*; z₀, ε) is a C¹ diffeomorphism of ℝ² (a time-t* map of a smooth global flow), and ρ has a
  density. So (x, p)(t*) has a density, x(t*) is non-degenerate, and **Var F_{P1,N}(t*) > 0**, which is finite by the
  BRI1-R1 moment bounds.
- So t* lies outside the degeneracy set for both protocols, and the standardised forces are defined.

**Skewness.**
- **P0:** γ(F_{P0,N}(t*)) = 0.
- **P1, for every N_B ≥ N₀:** γ(F_{P1,N}(t*)) = κ₃ / Var^{3/2} ≠ 0.

**Reflection safety.** Suppose the X1 family were in E₂±. By BRI-C2± (§R2 form), there is a shared causal sign functional
S with s_q·Z_q sharing one law across 𝒳. At the single time t*, this gives Z_{P1}(t*) ~ ±Z_{P0}(t*). Then
|γ(Z_{P1}(t*))| = |γ(Z_{P0}(t*))| = 0, contradicting γ(F_{P1,N}(t*)) ≠ 0. This is a one-time **reflection-orbit
violation** (BRI-O4 (a)), with no degeneracy involved, so it is not BRI-DEG.

**THEOREM BRI1-X1.** There exist:
- the frozen Candidate-1 protocol **P1**;
- a time **t* ∈ (0, δ)**;
- a finite integer **N₀**;

such that **for every finite bath size N_B ≥ N₀, the X1 interventional force family lies outside E₂±**. Equivalently,
X1_{N_B} ∉ E₂± for all sufficiently large finite N_B. The witness is a **BRI-E2+O** (orbit / shape) violation against the
reference protocol P0. ∎

## T4 — Reservoir limit (kept separate)

**Third cumulant.** For fixed t*, κ₃(F_{P1,N}(t*)) = K/N_B + O(N_B⁻²) → 0 and Var → Var x₀(t*) = m₂, so the standardised
skewness → 0 at rate N_B⁻¹.

**Limit law (proved here without ASSUMPTION R).** Fix q ∈ {P0, P1, P2} and a finite tuple of times in [0, 2π]. Write
F̊_{q,N} = N_B^{−1/2}·Σ_{j=1}^{N_B} (X_j^ε − E X^ε), with ε = N_B^{−1/2}: a triangular array of i.i.d. rows.
- By BRI1-R1 (iii), Cov(X^ε(t_a), X^ε(t_b)) is continuous in ε, so it tends to C₀(t_a, t_b) = ⟨x₀(t_a)x₀(t_b)⟩_Gibbs.
- Fourth moments are uniformly bounded for |ε| ≤ 1.
- The Lyapunov condition therefore holds: N_B·E|N_B^{−1/2}(X − E X)|⁴ = N_B⁻¹·O(1) → 0.
- The multivariate Lindeberg–Feller CLT (KNOWN) gives fdd convergence of F̊_{q,N} to the **centred Gaussian with
  covariance C₀, the same for every protocol**.
- So the limiting family has a q-invariant centred law: it is in **E₁** at the fdd level (BRI-C1). The deterministic mean
  E F_{q,N} = √N_B·E X^ε → E y₁ is the linear (FDT) response, absorbed by M.

**Resulting structure (stated carefully):**
- X1 is **outside E₂± at every sufficiently large finite bath size**, by THEOREM BRI1-X1.
- Its non-affine witness (standardised skewness at t*) is O(N_B⁻¹) and **vanishes** as N_B → ∞.
- The limiting centred force is the protocol-independent Gaussian of the E₁ / linear-response class.

**Not claimed:**
- any topological discontinuity of E₂± membership (no topology or closure statement is proved);
- that the finite-N_B effect stays experimentally observable as N_B grows.

## T5 — Relation to the frozen τ (firewall)

**The PF4Q grade is unchanged.** It stays **PF4Q-I, X1-PF-INDETERMINATE AT THE FROZEN τ**: finite-time non-vanishing at
τ = (π, 3π/2, 2π) is not certified.

THEOREM BRI1-X1 is a **distinct** class-level result. It is **not data fishing**:
- the non-zero leading coefficient C₇ = −Var(x₀²)/(28π³) was derived analytically in the preflight (`6220a9e`), before
  any numerical coefficient existed;
- **no numerical scan was used to select t***: the statement is existential over the interval (0, δ), which is forced by
  the sign of the first non-zero Taylor coefficient;
- the protocol P1 was frozen before any computation;
- the PF4Q numerical values play no role in the proof.

The frozen-τ witness and its certification route stay documented as a possible future verification project, and are not
pursued.

## T6 — Grade

**BRI1-X1-THEOREM — ANALYTIC E₂± ORBIT ESCAPE AT FINITE BATH SIZE: PROVED** (INTERNALLY PROVED / NOT EXTERNALLY
REVIEWED).

**BRI-E2+O — EARNED AT CLASS-THEOREM LEVEL FOR X1**, with the precise qualifier: *existence for all sufficiently large
finite N_B; no explicit threshold N₀ (nor δ); frozen-τ certification remains indeterminate.*

**The earned statement.** A finite reciprocal anharmonic environment (X1, finite Duffing bath, 1/√N_B coupling) produces
an interventional force law that **cannot be represented as one shared causal affine modulation of a
protocol-independent exogenous process**. At the same time, the distinguishing non-affine contribution **decays to zero
in the reservoir limit**, where the centred force returns to the protocol-independent Gaussian (E₁) class.

**This is an identifiability theorem relative to E₂±. It is not:**
- primitive randomness or a unique ontology;
- a GRUT-specific law or a GRUT empirical prediction;
- consciousness or TRUE COMPRESSION;
- an escape from E_univ (impossible: BRI-UPPER; X1 ∈ E_univ, PF-2).

**Relation to SCOUT-0 (P-17 / P-18).** The P-17 / P-18 non-identifiability (class 𝓗 ⊂ E₁; one-way drivers) **does not
extend through finite reciprocal anharmonic baths** under the affine competitor class E₂±. The escape is O(N_B⁻¹) and
disappears macroscopically. This is consistent with, and quantifies, why the harmonic / linear-response description is
hard to escape for large reservoirs. The SCOUT-0 saturation verdict is not reopened: BRI0 is a new-premise campaign.

**X1 numerical evidence:** the 20 frozen-τ PF4Q values are banked unchanged as supporting evidence only. They were not
rerun and not used. No certified pipeline, no N_B grid run, and no D_orb.
