# BRI0 · CANDIDATE 1 — CHARTER AND ANALYTIC PREFLIGHT (NOT NUMERICALLY RUN)

**Status:** CHARTER + ANALYTIC PREFLIGHT ONLY. **No Monte Carlo, no quadrature, no numerical D_orb, no sampling.**

**Labels:**
- symbolic checks are **exact power-series identities** (`preflight/bri1_preflight_symbolic.py`, log `.log`);
- analytic statements are **PROVED HERE (INTERNALLY PROVED / NOT EXTERNALLY REVIEWED)**;
- stated hypotheses are marked **ASSUMPTION**.

**Governing documents:** `BRI0_CHARTER.md` (Charter 0, plus §R Scope Repair 01, plus §R2 Scope Repair 02).

**Base:** `grut-backreaction-identifiability-0 @ 8df4b88`. Frozen parent: `scout-0 @ ab2da47`.

## 1. Frozen candidate definitions

### C1 — exact harmonic nonlinear-coupling control

**Hamiltonian:**
H = H_S(q, p_q) + Σ_j [p_j²/2 + (ω_j²/2)(x_j − c_j A(q)/ω_j²)²], with **A(q) = q + q³/3**.
- A′(q) = 1 + q² ≥ 1 > 0 for all q, so there are no sign or zero artefacts.
- **Preparation:** at q(0) = 0 (so A(0) = 0), the bath is thermal at temperature T₀:
  x_j(0) ~ N(0, T₀/ω_j²) and p_j(0) ~ N(0, T₀), independently.
- **Force convention:** F = −∂H_B/∂q.

### X1 — reciprocal finite Duffing environment

**Hamiltonian:**
H = p_q²/(2M) + V(q) + Σ_{j=1}^{N_B} [p_j²/2 + x_j²/2 + x_j⁴/4] − (g/√N_B)·q·Σ_j x_j, with **g = 1** and
**V(q) = q²/2 + q⁴/4**.
- **Sign convention (fixed here):** the coupling carries a minus sign, so that the clamp equations take the owner's
  form, x_j″ + x_j + x_j³ = q(t)/√N_B, and the environment force on q is
  **F_q(t) = −∂H_int/∂q = (1/√N_B)·Σ_j x_j(t)**.
- V is irrelevant to every clamp calculation; it only defines the autonomous parent.
- **Preparation:** at q(0) = 0, the pairs (x_j, p_j) are i.i.d. with ρ ∝ exp[−(p²/2 + x²/2 + x⁴/4)]. There is no
  protocol-dependent re-preparation.
- **Notation:** ε := N_B^{−1/2}.

### Protocols (frozen; T = 2π)

With s(u) = 10u³ − 15u⁴ + 6u⁵:
- **P0:** q ≡ 0.
- **P1:** q = s(t/π) on [0, π]; q = 1 on (π, 2π].
- **P2:** q = s(t/π) on [0, π]; q = 1 + s((t − π)/π) on (π, 2π].

q1 and q2 agree exactly on [0, π]. All three lie in 𝒳: they are C², with q(0) = q̇(0) = 0. **No protocol is added or
changed.**

### Witness tuple, primary witness and bath-size axis (frozen)

**Witness tuple:** τ = (π, 3π/2, 2π). **Carrier reference:** P0. It is non-degenerate on τ for X1, since
Var F₀(t) = m₂ > 0 (§3), and for C1, since Var ξ(t) = T₀Σ_j c_j²/ω_j² > 0.

**Primary population witness:**
D_orb(q, 0) = min_{s∈{±1}³} ED(Law Z_q(τ), Law diag(s)Z₀(τ)), with ED the population energy distance.
- D_orb ≥ 0.
- **D_orb = 0 iff the two 3-dimensional laws agree under some coordinatewise reflection**, given finite first moments.
  This uses the characterisation property of the energy distance (KNOWN: Székely–Rizzo).
- There is no bandwidth or fitted scale.
- **D_orb > 0 is a direct BRI-E2+O witness.**
- It is **not computed** in this preflight.

**Secondary analytic witnesses:** |standardised skewness|, excess kurtosis, and |mixed standardised cumulants|. They
cannot replace D_orb.

**Bath-size axis (frozen for any later numerical stage):** N_B ∈ {1, 2, 4, 8, 16, 32, 64}. No sample count is chosen.

## 2. PF-1 — C1 is exactly in E₂± (PASS)

**PROP BRI1-C1.** Under any clamp q ∈ 𝒳,
F_q(t) = M_t[q] + A′(q(t))·ξ(t), where:
- ξ(t) := Σ_j c_j [x_j(0) cos ω_j t + (p_j(0)/ω_j) sin ω_j t];
- M_t[q] := −A′(q(t))·∫₀ᵗ γ(t − s)·A′(q(s))·q̇(s) ds;
- γ(t) := Σ_j (c_j²/ω_j²) cos ω_j t.

**Proof.**
1. The bath obeys ẍ_j = −ω_j² x_j + c_j A(q(t)), which is linear and driven by the deterministic path A(q(·)). This is
   P-17's proof with q replaced by A(q).
2. Variation of constants plus one integration by parts give
   Σ_j c_j (x_j − c_j A(q)/ω_j²) = ξ(t) − ∫₀ᵗ γ(t − s) (d/ds)A(q(s)) ds, using A(q(0)) = 0.
3. The force is F = −∂H_B/∂q = A′(q)·Σ_j c_j (x_j − c_j A(q)/ω_j²).
4. ξ depends only on the bath's initial data, whose law is protocol-independent. M is deterministic and causal.
   G_t[q] := A′(q(t)) = 1 + q(t)² is deterministic, causal and **≥ 1, never zero**. ∎

**Exact symbolic check.** For one mode and a generic clamped q = q₃t³ + q₄t⁴, the difference
force − A′(q)[ξ − ∫γ dA] vanishes identically through t¹¹ (`bri1_preflight_symbolic.log`, line 1).

**Further consequences:**
- D_q = ∅ (thermal Var ξ(t) > 0 is constant).
- **C1 ∉ E₁** for any protocol with q ≢ 0, because Var F_q(t) = (1 + q(t)²)²·Var ξ(t) changes with the protocol.
- So **C1 is structurally BRI-E1** (E₁ escaped, E₂± not escaped). This is analytic, and the control needs no
  simulation.

**PF-1 terminal: C1 IN E₂± EXACTLY.** The preflight continues to X1.

## 3. PF-2 — X1 well-posedness and reciprocity (PASS)

1. **Bounded below.** For each j, x_j²/2 − ε q x_j ≥ −ε²q²/2, so Σ_j[x_j²/2 + x_j⁴/4 − ε q x_j] ≥ −N_B ε² q²/2 = −q²/2.
   Hence H ≥ p_q²/2M + q⁴/4 + Σ_j p_j²/2 ≥ 0 for every N_B.
2. **Reciprocal.** The single term −ε q Σ_j x_j produces both the force on the bath (+ε q) and the force on the system
   (+ε Σ_j x_j).
3. **Clamp equations.** These are ẍ_j = −x_j − x_j³ + ε q(t). With E_j = p_j²/2 + x_j²/2 + x_j⁴/4,
   dE_j/dt = ε q(t) p_j ≤ ε|q|√(2E_j), so √E_j grows at most linearly. Solutions are global on [0, 2π] and smooth in the
   initial data. X1 is in the BRI-UPPER parent class, hence **X1 ∈ E_univ**.
4. **Preparation.** ∫exp(−x²/2 − x⁴/4)dx < ∞ and p is Gaussian, so it is normalisable. At q(0) = 0 the coupling
   vanishes, so the canonical q = 0 law is the product of single-Duffing Gibbs laws at temperature 1. It is the same for
   every protocol.

**PF-2 terminal: PASS.**

## 4. PF-3 — Bath-size perturbation structure (PROVED HERE, under ASSUMPTION R)

### 4.1 Exact independence structure (no approximation)

Under a clamp, the N_B oscillators have i.i.d. initial data and are each driven by the **same deterministic** ε q(t).
Their trajectories are therefore i.i.d. copies of one single-oscillator process X^ε_q(·), and
**F_q = ε·Σ_{j=1}^{N_B} X_j exactly.** Hence, for every finite time tuple and every order n,

  **κ_n(F_q(t₁), …, F_q(t_n)) = N_B^{1 − n/2}·κ_n(X^ε_q(t₁), …, X^ε_q(t_n)).**   (★)

**Corollaries:**
- The equilibrium force fluctuation is O(1): Var F₀(t) = N_B ε²·Var x(t) = m₂ := E_Gibbs[x²], with m₂ ∈ (0, 1) (§6).
- The perturbation felt by each oscillator is O(ε) = O(N_B^{−1/2}).

This **derives** the scaling claimed in the candidate selection.

### 4.2 Parity (exact)

The single-oscillator dynamics and the Gibbs law are invariant under (x, p, q) ↦ (−x, −p, −q). Write
X^ε_q = x₀ + ε y₁ + ε² y₂ + ε³ y₃ + …, where:
- x₀ is the unforced Gibbs trajectory;
- ÿ₁ + (1 + 3x₀²) y₁ = q(t);
- ÿ₂ + (1 + 3x₀²) y₂ = −3 x₀ y₁²;
- …, all with zero initial data.

By induction, **y_k is odd in x₀ for even k and even in x₀ for odd k**:
- y₁ depends on x₀ only through x₀², so it is even;
- y₂ is driven by x₀y₁², which is odd;
- each later y_k is driven by products whose x₀-parity alternates.

y_k is homogeneous of degree k in q. Expectations of x₀-odd functionals vanish.

### 4.3 ASSUMPTION R

The ε-expansion of X^ε_q(t) holds in L^p for every p on the finite horizon, with remainder O(ε^{K+1}).
- *Plausibility:* the flow is smooth in parameters on a finite horizon, and the Gibbs tails exp(−x⁴/4) give all moments.
- A proof would need moment bounds on the linearised Duffing flow uniformly over high-energy initial data. **Not proved
  here.**

### 4.4 Order of each effect (from (★) and parity)

**Deterministic mean response:**
- E F_q(t) = √N_B·E[X^ε_q] = E y₁(t) + N_B^{−1}·E y₃(t) + O(N_B^{−2}).
- The O(1) term is the equilibrium linear response, E y₁(t) = ∫₀ᵗ χ(t − s) q(s) ds. By the classical FDT at temperature
  1, χ(t) = −(d/dt)⟨x(t)x(0)⟩_Gibbs (KNOWN).
- Deterministic and causal, so absorbed by M.

**Cumulant departures from P0:** by (★), κ_n(F_q) − κ_n(F₀) = N_B^{1−n/2}·[κ_n(X^ε) − κ_n(x₀)]. The first non-zero
ε-order of the bracket is **ε¹ for odd n** (n·κ(x₀^{n−1}, y₁), which is x₀-even) and **ε² for even n** (the ε¹ term
κ(x₀^{n−1}, y₁) is x₀-odd and vanishes). Therefore:

| cumulant order n | change relative to P0 | role |
|---|---|---|
| 2 (covariances) | **O(N_B^{−1})** | per-time variances are amplitude (absorbed by G); a change of the **correlation coefficients** between distinct times would be non-affine |
| 3 (third cumulants, incl. mixed) | **O(N_B^{−1})** | zero at P0 exactly (parity), so any non-zero value is non-affine |
| 4 | O(N_B^{−2}) | |
| 5 | O(N_B^{−2}) | |
| n | N_B^{(1−n)/2} for odd n; N_B^{−n/2} for even n | |

**There is no O(N_B^{−1/2}) change of the force law.** The pathwise correction ε²Σ_j(y₁ − E y₁) = N_B^{−1/2}·W₁ exists,
but it enters every cumulant only through x₀-odd cross terms, which vanish.

**Reservoir limit.** As N_B → ∞ with q fixed, the multivariate CLT applies to the i.i.d. sums. The centred force
converges to a Gaussian with q-independent covariance C₀(t, t′) = ⟨x(t)x(t′)⟩_Gibbs. **X1 → E₁ in the reservoir limit**,
i.e. Caldeira–Leggett-type linear-response universality (KNOWN as a mechanism).

**The registered mesoscopic prior is therefore derived:** any E₂± escape in X1 is an **O(N_B^{−1})** effect and vanishes
in the reservoir limit.

**The three-way separation PF-3 asks for:**
- **deterministic mean response:** O(1) linear plus O(N_B^{−1}) cubic, absorbed by M;
- **pure amplitude:** O(N_B^{−1}) per-time variance change, absorbed by G;
- **genuinely non-affine:** O(N_B^{−1}) third (odd) standardised cumulants, and any O(N_B^{−1}) change of inter-time
  correlation coefficients.

## 5. PF-4 — First reflection-safe departure at the frozen τ

**Lower orders vanish identically (PROVED):**
- every odd standardised (mixed) cumulant of Z_q(τ) is **exactly zero at P0 for all N_B** (parity), and its departure for
  P1 / P2 is O(N_B^{−1}) or smaller;
- the even-cumulant shape changes (excess kurtosis) are O(N_B^{−2});
- there is **no reflection-safe departure at O(N_B^{−1/2})**, for any protocol or time.

**The first possible order is O(N_B^{−1}).** Its coefficients are explicit:
- third mixed cumulants: κ₃(F_q(t_a), F_q(t_b), F_q(t_c)) = N_B^{−1}·[c(t_a, t_b; t_c) + c(t_a, t_c; t_b) +
  c(t_b, t_c; t_a)] + O(N_B^{−2}), with **c(t_a, t_b; t_c) := Cov(x₀(t_a)x₀(t_b), y₁(t_c))**;
- standardised skewness: γ₁[F_q(t)] = 3c(t, t; t) / (N_B·m₂^{3/2}) + O(N_B^{−2}).

Any non-zero value is a **reflection-orbit violation**: P0's odd cumulants vanish, a reflection changes only signs, and
so |·| > 0 cannot be mapped to 0. That makes it a **BRI-E2+O** witness at O(N_B^{−1}).

**The coefficient is not identically zero (PROVED, exact symbolic series).** For P1 on [0, π], with Gibbs moments reduced
to m₂ by the exact identity m_{k+1} + m_{k+3} = k·m_{k−1}:

  **c(t, t; t) = −Var(x₀²)·t⁷/(28π³) + 3Var(x₀²)·t⁸/(112π⁴) + O(t⁹)**, with Var(x₀²) = m₄ − m₂² = 1 − m₂ − m₂² > 0.

- The log prints the t⁷ coefficient as (m₂² + m₂ − 1)/(28π³), which is exactly −Var(x₀²)/(28π³). It is the same
  expression the hand derivation gives: −3·Var(x₀²)·(10/π³)·t⁷/840.
- Var(x₀²) > 0 because x₀² is non-degenerate under Gibbs.
- So **c(t, t; t) < 0 for all sufficiently small t > 0: X1 has a genuine non-affine O(N_B^{−1}) response.** Physically,
  a stiffer, higher-energy oscillator responds less, so response amplitude anticorrelates with x₀².

**Obstruction at the frozen tuple.** τ begins at t = π. There, c(·) and its mixed versions are Gibbs averages of the
linearised Duffing flow over finite times. They have no closed form, and the small-t series cannot certify their value
at t ≥ π:
- the Duffing solutions are elliptic, and their t-series radius shrinks with energy;
- the Gibbs ensemble includes arbitrarily high energies.

So **the analytic preflight cannot certify that the O(N_B^{−1}) coefficients are non-zero at τ = (π, 3π/2, 2π).**
- They are generically non-zero: they are non-identically-zero functionals, analytic in t on (0, π) and (π, 2π).
- An accidental zero at all ten mixed third cumulants on τ is not excluded analytically.
- **No time or protocol has been changed to rescue the candidate.**

## 6. PF-5 — CSI possibility at the frozen tuple (PROVED: NOT POSSIBLE via BRI-O11)

**Reference symmetry.** The reference law Z₀(τ) is invariant under the global flip −I, exactly for every N_B (parity,
§4.2).

**Consequence.** By the §R2 structural lemma, each non-empty A_q(τ) projects onto **both** signs at t* = π. Condition 4
of BRI-O11 (disjoint projections) can therefore never hold. **P1 / P2 cannot yield a BRI-E2+C certificate at
τ = (π, 3π/2, 2π)**, because their shared prefix contains only one time of τ.

- This is a structural property of X1 with carrier P0, not a sampling issue.
- CSI is **not** manufactured.
- The primary target remains BRI-E2+O.

**Moment fact used in §4 and §5:** m₂ ∈ (0, 1) follows from m₂ + m₄ = 1 with m₄ > 0. Var(x₀²) = 1 − m₂ − m₂² > 0
follows from non-degeneracy, which also bounds m₂ < (√5 − 1)/2.

## 7. PF-6 — Preflight terminal

**X1-PF-INDETERMINATE**, with the following precise analytic obstruction.

**Proved:**
- C1 lies in E₂± exactly (PF-1).
- X1 is well posed, reciprocal, and has an identical preparation (PF-2).
- Exact structure (★), parity, the **orders** of every effect, and the reservoir limit X1 → E₁ (PF-3, under ASSUMPTION R).
- **No reflection-safe departure below O(N_B^{−1})**.
- The O(N_B^{−1}) non-affine coefficient is **non-zero as a functional**: c(t, t; t) ≈ −Var(x₀²)t⁷/(28π³) < 0 for small t
  under P1.
- CSI is impossible at the frozen τ (PF-5).

**Not proved:** non-vanishing of the explicit O(N_B^{−1}) coefficients **at the frozen tuple τ**. They are finite-time
Gibbs averages of the linearised Duffing flow, which have no closed form and lie beyond the convergence of the exact
small-t series.

**Why not the other terminals:**
- **Not X1-PF-A:** the departure is not certified at τ.
- **Not X1-PF-B:** a genuinely non-affine response is proved to exist, so "affine to computed order" would be false.
- **Not X1-PF-C:** no E₂± representation exists generically, and none is constructed.
- **Not X1-PF-D:** no structural failure.

**Proposed resolution (owner decision; not done here).** A **deterministic, certified-error evaluation** of the ten
coefficients c_abc[q](τ) for P1 and P2 would decide between X1-PF-A (some coefficient ≠ 0 at τ) and an accidental
O(N_B^{−1}) zero at τ.
- Each coefficient is a two-dimensional Gibbs integral over (x₀, p₀) of ODE solutions, so a sampling-free quadrature is
  possible.
- This is numerical, so it awaits owner authorisation.

**Feasibility notes for any later numerical charter** (no numbers computed):
- By (★), an N_B scan needs only **single-oscillator** ensembles at each ε, not N_B-body simulations.
- The witness magnitudes scale as N_B^{−1} (|skewness|) and D_orb ~ O(N_B^{−2}) (quadratic in a small law
  perturbation). Detecting them at N_B = 64 requires very large sample sizes, so it should be planned explicitly.

## 8. Claim firewall

Even an eventual X1-PF-A, or a BRI-E2+, establishes only this: **finite nonlinear reciprocal environments can have
interventional response laws outside a shared causal affine exogenous representation.** It establishes no GRUT physics,
GRUT prediction, primitive randomness, unique ontology, consciousness or TRUE COMPRESSION, and no escape from E_univ.

**Literature context (metadata only; no novelty claimed):** nonlinear generalised Langevin descriptions and general
nonlinear system–bath couplings are established territory, which is why C1 is a control.

**Hard stop.** No Monte Carlo, no quadrature, no numerical D_orb, no protocol / τ / N_B changes. Stopping for owner
review.
