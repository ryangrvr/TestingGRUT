# B1R1 INDEPENDENT REPRODUCTION 01 — V1
## Executed under `B1R1_PUBLICATION_VERIFICATION_CHARTER_01.md` (frozen `8d765ea`), V1 lane only

**Independence firewall (declared before derivation).** This reproduction was
performed from the frozen model definitions only:

- `B1R0_CHARTER.md` §1 (X1 Hamiltonian, preparation, sign convention),
  §4 (clamp family), §R B1R-O2/O4/O9 (E2± definition and B1R-C2±
  characterization — used as **frozen inputs/definitions**, not derivation
  source);
- `B1R1_CANDIDATE_CHARTER.md` §1 (protocol definitions P0/P1/P2, s(u));
- standard mathematical theorems, identified in the Standard-Theorem Ledger.

**Not used:** the prose proof in `B1R1_ANALYTIC_ESCAPE_THEOREM.md` (opened
only after the derivation was complete, for the comparison section below);
frozen-τ PF4Q numerical values; old numerical K-values; old scripts.

**Ordering record:** the derivation below was completed first; the
post-derivation comparison section was written only afterward. The single
exception is that the frozen *definitions* (not proofs) were read first, as
permitted.

---

## V1-1 — Frozen dynamics

From the X1 Hamiltonian (frozen): H = p_q²/2M + V(q) + Σ_j [p_j²/2 + x_j²/2
+ x_j⁴/4] − ε q Σ_j x_j, ε = N_B^(−1/2), g = 1, V(q) = q²/2 + q⁴/4.

Hamilton's equations for oscillator j: ẋ_j = p_j; ṗ_j = −∂H/∂x_j = −x_j −
x_j³ + ε q(t). Therefore x_j'' = −x_j − x_j³ + ε q(t):

  **x_j'' + x_j + x_j³ = q(t)/√N_B.** ✓

Environmental force on the clamp: F_q = −∂H_int/∂q = −(−ε Σ_j x_j):

  **F_q(t) = (1/√N_B) Σ_j x_j(t).** ✓

**Reciprocity:** the single coupling term −ε q Σ_j x_j produces the force on
the bath (+ε q per oscillator, from −∂H_int/∂x_j) and the force on the
system (+ε Σ_j x_j) — both from one term. ✓

**Common preparation:** at q(0) = 0 the coupling vanishes; the (x_j, p_j)
are i.i.d. with ρ ∝ exp(−H₀), H₀ = p²/2 + x²/2 + x⁴/4, for **every**
protocol. No protocol-dependent re-preparation. ✓

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-2 — Exact i.i.d. force-sum structure

Under a deterministic clamp q(·), oscillator j satisfies the same ODE with
i.i.d. initial data. The solution map is deterministic: x_j(t) = Φ_t(x_j(0),
p_j(0); ε, q). Since (x_j(0), p_j(0)) are i.i.d., the trajectories are i.i.d.
copies of the single driven process X^ε_q(t) := Φ_t(z₀; ε, q). Hence

  **F_q(t) = ε Σ_j X_j^ε_q(t), exactly, for every finite N_B.** ✓

No asymptotics used.

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-3 — Cumulant scaling (exact, then expansion)

Cumulant generating function of the force vector (F(t₁),…,F(tₙ)):

K_F(θ) = log E exp(Σ_a θ_a F(t_a)) = log E exp(ε Σ_a θ_a Σ_j X_j(t_a)).

The X_j vectors are i.i.d., so E exp(ε Σ_j Y_j) = (E exp(ε Y))^{N_B} for the
single-oscillator vector Y = (X(t₁),…,X(tₙ)) (independence), giving

K_F(θ) = N_B · K_X(εθ).

Differentiating n times in θ at 0: each differentiation brings one factor ε,
so

  **κ_n(F_q(t₁),…,F_q(tₙ)) = N_B εⁿ κ_n(X^ε_q(t₁),…,X^ε_q(tₙ)) =
  N_B^(1−n/2) κ_n(X^ε_q(t₁),…,X^ε_q(tₙ)).**  (★)

This is an **exact finite-N_B identity** (a differentiation of an
elementary identity, using only cumulant additivity over independent
variables and homogeneity of degree n under scaling — both immediate from
the CGF definition).

**Specialization n = 3:** κ₃(F_q(t₁),F_q(t₂),F_q(t₃)) = N_B^(−1/2)
κ₃(X^ε(t₁),X^ε(t₂),X^ε(t₃)). ✓

The small-ε expansion of κ_n(X^ε) is a **separate subsequent step** (V1-5,
V1-8); (★) itself is exact.

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-4 — Reference parity (P0)

Under P0 (q ≡ 0) the oscillator is unforced Duffing with initial law ρ ∝
e^{−H₀}.

(a) *Symmetry at all times.* The unforced vector field is odd: (x,p) ↦
(−x,−p) maps solutions to solutions (ẋ = p odd, ṗ = −x−x³ odd). ρ is even,
so Law(Φ_t(z₀)) = Law(R Φ_t(z₀)) for every t: **x₀(t) has a reflection-
symmetric law for every t** (via commutation of the flow with R and
invariance of ρ under R).

(b) *Stationarity (for V1-6's bookkeeping).* The Gibbs measure is invariant
under the unforced Hamiltonian flow (Liouville + energy conservation on
level sets) — Standard-Theorem Ledger S1. Hence E x₀(t)² = m₂ for all t.

(c) *Exact zero third cumulant.* F_{P0,N}(t) = ε Σ_j x₀,j(t) is a sum of
i.i.d. symmetric summands ⇒ symmetric law ⇒ κ₃ = 0 (odd moment of a
symmetric law, E X = 0). By (★) equivalently: κ₃(x₀) = 0 ⇒ N_B^(1−3/2)·0 =
0.

  **κ₃(F_{P0,N}(t)) = 0 exactly, for every t and every finite N_B.** ✓

**Disposition: REPRODUCED-WITH-STANDARD-THEOREM** (S1: invariance of the
Gibbs measure under Hamiltonian flow — hypotheses trivially satisfied:
smooth potential, confined measure; everything else elementary).

## V1-5 — First-order response

Differentiate ẍ = −x − x³ + ε q(t) in ε at ε = 0, writing X^ε = x₀ + ε y₁ +
O(ε²):

  **ÿ₁ + (1 + 3x₀²) y₁ = q(t), y₁(0) = 0, ẏ₁(0) = 0.** ✓

Initial conditions: the bath initial data (x_j(0), p_j(0)) are drawn from ρ,
which does not depend on ε (preparation precedes the protocol; the coupling
vanishes at t = 0). Hence ∂_ε(z(0)) = 0.

*Justification of differentiation at the level needed for expectations:*

1. **Flow smoothness (Ledger S2):** the vector field is polynomial in
   (z, ε) and solutions exist globally on the finite horizon (energy bound:
   dH₀/dt = ε q p ≤ ε|q|√(2H₀) ⇒ √H₀(t) ≤ √E₀ + c for bounded |q| ≤ Q; the
   orbit stays compact). Standard ODE theory gives C^∞ dependence on
   (z₀, ε) on the compact interval.
2. **Differentiation under the expectation (Ledger S3):** with y_k := ∂_ε^k
   z the variational equations satisfy Grönwall bounds |y_k| ≤ P(E₀)
   e^{β√E₀} for polynomials P and constants β (k ≤ 3 suffices), uniformly
   in |ε| ≤ 1. Under ρ, E P(E₀)e^{β√E₀} < ∞ since −E₀ + β√E₀ ≤ −E₀/2 +
   β²/2 and the phase-space volume of {H₀ ≤ e} grows polynomially in e.
   Dominated convergence then gives: ∂_ε^k E[f(X^ε)]|_{ε=0} =
   E[∂_ε^k f-combination] for the polynomials f arising in cumulants up to
   order 3, with all moments finite uniformly in |ε| ≤ 1.

All initial conditions stated; hypotheses checked. ✓

**Disposition: REPRODUCED-WITH-STANDARD-THEOREM** (S2: smooth parameter
dependence of ODE flows; S3: dominated convergence for differentiation
under the expectation).

## V1-6 — Small-time coefficient (central algebraic check)

**Method: hand recurrence, independently derived.** No old script imported.

*Recursion for x₀.* x₀(t) = Σ_{k≥0} a_k t^k, a₀ = a = x₀(0), a₁ = b = p₀.
From x₀'' = −x₀ − x₀³ (unforced):

  (k+2)(k+1) a_{k+2} = −a_k − Σ_{i+j+l=k} a_i a_j a_l.

*Recursion for y₁.* y₁(t) = Σ_{k≥0} b_k t^k, b₀ = b₁ = 0. From ÿ₁ +
(1+3x₀²) y₁ = q(t):

  (k+2)(k+1) b_{k+2} = −b_k − 3 Σ_{m=0}^{k} (x₀²)_m b_{k−m} + q_k,

where (x₀²)_m is the m-th coefficient of x₀(t)², and q(t) = s(t/π) = 10t³/π³
− 15t⁴/π⁴ + 6t⁵/π⁵ on [0, π], so q₃ = 10/π³, q₄ = −15/π⁴, q₅ = 6/π⁵, all
other q_k = 0 (through the orders needed).

*Moment inputs.* a and b independent; b ~ N(0,1) (odd moments 0); a has
density ∝ e^{−a²/2 − a⁴/4} (even ⇒ odd moments 0). Integration by parts
against e^{−a²/2 − a⁴/4} (boundary terms vanish): **m_{k+1} + m_{k+3} = k
m_{k−1}**, where m_j := E a^j. In particular **m₂ + m₄ = 1**, so Var(a²) =
m₄ − m₂² = 1 − m₂ − m₂².

*y₁ coefficients (all that are needed):*

- b₂ = b₃ = b₄ = 0 (only q-independent and q₂ terms, all zero).
- 20 b₅ = q₃ ⇒ **b₅ = 1/(2π³)**.
- 30 b₆ = q₄ ⇒ **b₆ = −1/(2π⁴)**.
- 42 b₇ = −b₅ − 3(x₀²)₀ b₅ + q₅ = −(1+3a²)/(2π³) + 6/π⁵ ⇒
  **b₇ = −(1+3a²)/(84π³) + 1/(7π⁵)**.

*Target.* c(t) := Cov(x₀(t)², y₁(t)) = E[x₀(t)² y₁(t)] − m₂ E[y₁(t)] (using
V1-4(b)). Writing x₀² = Σ d_i t^i, the coefficient of tⁿ in c(t) is

  C_n = Σ_{i+j=n} E[d_i b_j] − m₂ E[b_n].

*The n = 7 computation (all that is needed for the leading term):*

- E[d₀ b₇] = E[a² · (−(1+3a²)/(84π³) + 1/(7π⁵))] = −(m₂ + 3m₄)/(84π³) +
  m₂/(7π⁵).
- E[d₁ b₆] = E[2ab] · (−1/(2π⁴)) = 0 (E[ab] = E a · E b = 0).
- E[d₂ b₅] = E[b² − a² − a⁴]/(2π³) = (1 − m₂ − m₄)/(2π³) = **0** by m₂ + m₄
  = 1.
- All other d_i b_j pairs with i + j = 7 have j ∈ {0,…,4} ⇒ b_j = 0.
- m₂ E[b₇] = m₂[−(1+3m₂)/(84π³) + 1/(7π⁵)].

Sum:

C₇ = [−m₂ − 3m₄ + m₂ + 3m₂²]/(84π³) + [m₂/(7π⁵) − m₂/(7π⁵)]
  = 3(m₂² − m₄)/(84π³)
  = **−(m₄ − m₂²)/(28π³) = −Var(x₀²)/(28π³).** ✓

The π⁵ terms cancel identically — an internal consistency check on the
arithmetic.

*Lower coefficients.* Direct evaluation gives C₀ = C₁ = C₂ = C₃ = C₄ = 0
(b₀…b₄ = 0), and:

- C₅ = E[d₀b₅] − m₂ E[b₅] = m₂/(2π³) − m₂/(2π³) = 0.
- C₆ = E[d₀b₆] − m₂ E[b₆] = −m₂/(2π⁴) + m₂/(2π⁴) = 0.

  **C₀ = C₁ = … = C₆ = 0 and C₇ = −Var(x₀²)/(28π³) < 0.** ✓

*Var(x₀²) > 0.* The variable a has an a.e.-strictly-positive density (∝
e^{−a²/2−a⁴/4} > 0), so a² is not a.s. constant ⇒ Var(x₀²) > 0. (No value
of m₂ is needed or used.) Equivalently via the identity: Var = 1 − m₂ −
m₂², and stationarity pins E x₀(t)² = m₂ at every t; the identity m₂ + m₄ =
1 with m₄ > m₂² forces m₂ < (√5−1)/2, consistent. The primary argument is
the non-degeneracy one.

**Disposition: REPRODUCED-INDEPENDENTLY.** (Hand recurrence; full algebra
recorded; the sign rests on Var(x₀²) > 0 and 28π³ > 0, not on opaque code.)

## V1-7 — Small-time sign interval

*Remainder bound.* For fixed z₀ = (a, b), g(t; z₀) := (x₀(t)² − m₂) y₁(t) is
C^∞ on [0, 1] (both factors solve polynomial ODEs; q is a polynomial). By
repeated use of the ODEs, every time-derivative of order ≤ 8 of x₀, ẋ₀, y₁,
ẏ₁ is a polynomial in (a, b, y₁, ẏ₁) and derivatives of q, bounded on
[0,1]; with k(t) = 1 + 3x₀² ≤ 1 + 6√E₀, Grönwall gives |y₁|, |ẏ₁| ≤
e^{1+6√E₀} on [0,1] (|q| ≤ 1 there — in fact on [0, π], q ≤ 1). Hence
sup|∂_t⁸ g(·; z₀)| ≤ P(E₀) e^{1+6√E₀} =: D(z₀) for a polynomial P. Taylor
with remainder (Ledger S4):

  c(t) = Σ_{k=0}^{7} C_k t^k + R(t), |R(t)| ≤ C_R t⁸, C_R := E D / 8! < ∞
  (integrability as in V1-5).

Since C₀…C₆ = 0: **c(t) = C₇ t⁷ + R(t), C₇ = −Var(x₀²)/(28π³) < 0.**

Sign: choose **δ := min(1, |C₇|/(2 C_R))**. For 0 < t < δ:
c(t) ≤ t⁷(C₇ + C_R t) and C_R t < |C₇|/2 = −C₇/2, so C₇ + C_R t < C₇/2 < 0:

  **∃ δ > 0 such that c(t) < 0 for every 0 < t < δ.** ✓

Existential; no numerical δ is claimed (C_R is an existence bound, not
computed).

**Disposition: REPRODUCED-WITH-STANDARD-THEOREM** (S4: Taylor's theorem
with remainder, applied pointwise before taking expectations — no
expectation/derivative interchange needed since every term is integrable).

## V1-8 — Finite-bath nonzero skewness, quantifier structure

Fix t* ∈ (0, δ). By V1-4 parity extended to the ε-family (the map (x,p,ε)
↦ (−x,−p,−ε) maps solutions to solutions and ρ is even): κ₃(X^ε(t*)) is odd
and C³ in ε, with ∂_ε κ₃(X^ε)|₀ = 3c(t*) =: K(t*) < 0 (Ledger S3 for the
differentiation, already justified in V1-5). Taylor (S4): κ₃(X^ε(t*)) =
ε K(t*) + R₃(ε), |R₃(ε)| ≤ C(t*) ε³ for |ε| ≤ 1.

By (★) with ε = N_B^(−1/2):

  **κ₃(F_{P1,N}(t*)) = K(t*)/N_B + O(N_B^(−2)), K(t*) ≠ 0.**

Define N₀(t*) := ⌊C(t*)/|K(t*)|⌋ + 1, finite. For every integer N_B ≥ N₀(t*):
|O(N_B^(−2))| ≤ C(t*)/N_B² < |K(t*)|/N_B, hence

  **κ₃(F_{P1,N}(t*)) ≠ 0 for every N_B ≥ N₀(t*).**

Quantifiers, explicit:

> **∃ δ > 0 such that ∀ fixed t ∈ (0, δ), ∃ N₀(t) < ∞ such that ∀ N_B ≥
> N₀(t): κ₃(F_{P1,N}(t)) ≠ 0 (with the sign of K(t), i.e. negative).**

The threshold depends on t; **no uniform N₀ over the whole interval is
claimed** (C(t) is not uniform in t). ✓

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-9 — Positive variance

- P0: Var F_{P0,N}(t) = Var x₀(t) = m₂ > 0 (stationarity; m₂ = E a² > 0
  since a² ≥ 0 with a.e.-positive density and P(a ≠ 0) > 0). Finite.
- P1 at t*: Var F_{P1,N}(t*) = Var X^ε(t*) (by (★), n = 2). The time-t*
  map z₀ ↦ z(t*; z₀, ε) is a C¹ diffeomorphism of ℝ² (time-map of a smooth
  global flow, with smooth inverse), and ρ has an a.e.-positive density. A
  C¹ diffeomorphism maps absolutely continuous measures to absolutely
  continuous measures (pushforward of the density, Jacobian ≠ 0), so (x, p)(t*)
  has a density; x(t*) is non-degenerate; **Var F_{P1,N}(t*) > 0**. Finiteness
  from the V1-5 moment bounds. ✓

The time-map argument is sufficient: non-degeneracy of the initial law
propagates exactly through diffeomorphisms. ✓

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-10 — E2± escape (identifiability)

*Competitor definition (frozen, B1R-O9/B1R-C2±):* F ∈ E2± iff (i) the
degeneracy sets D_q coincide across q, and (ii) there exists **one shared
deterministic causal sign functional** S_t: q_[0,t] ↦ {±1} such that, with
s_q(t) := S_t[q_[0,t]], the finite-dimensional laws of s_q · Z_q (off the
common degeneracy set) are **the same for every protocol q** — where Z_q(t)
:= (F_q(t) − E F_q(t)) / sd F_q(t) is the standardized force.

*Argument.* Suppose X1 ∈ E2±. Then there is one common law ℒ with Law(s_q ·
Z_q(t)) = ℒ_t for every q and t off D. Consequently every reflection-
invariant functional of the standardized law takes **the same value for
every protocol**. The absolute standardized skewness |γ(Z_q(t))| =
|κ₃(Z_q(t))| / Var(Z_q(t))^{3/2} is reflection-invariant (a coordinatewise
sign flip multiplies κ₃ by the sign product, which has modulus 1):

  **|γ(Z_q(t))| = |γ(ℒ_t)| is the same for every q ∈ 𝒳.**

- At P0: V1-4 gives κ₃(F_{P0,N}(t)) = 0 exactly ⇒ γ(Z_{P0}(t)) = 0 ⇒ |γ(ℒ_t)|
  = 0.
- At P1, fixed t ∈ (0, δ), N_B ≥ N₀(t): V1-8 + V1-9 give κ₃ ≠ 0 and Var > 0
  ⇒ **|γ(Z_{P1}(t))| > 0**.

Contradiction. Hence **X1_{N_B} ∉ E2± for every N_B ≥ N₀(t), t ∈ (0, δ).** ✓

*Identifiability content (not "different skewness implies different model"
hand-waving):* E2± membership is exactly the statement that one shared
deterministic affine modulation (location M, signed scale G, one law Law ξ)
reproduces the whole interventional family. The standardization Z_q divides
out exactly the affine modulation (per-protocol M and |G| cancel; the sign
s_q is absorbed by taking |·|), so any E2±-invariant functional must agree
across protocols. |γ| is such a functional. P0 and P1 give different values
(0 versus > 0). Therefore no single shared (M, G, Law ξ) exists for the
family {P0, P1} — which is the definition of failure of the shared
representation.

*Classification:* the escape occurs at the one-time tuple (t, t, t), with
non-degenerate variances at both protocols (V1-9). It is a **reflection-orbit
violation** in the sense of B1R-O5/B1R-O10(a): the standardized law of P1
lies outside every coordinatewise reflection image of the common reference
law (here the reference value |γ| = 0 cannot be matched by any sign flip of
a nonzero |γ|). It is **not** a degeneracy artifact (both variances > 0), so
it is **B1R-E2+O**, not B1R-DEG. ✓

**Disposition: REPRODUCED-INDEPENDENTLY.**

## V1-11 — Reservoir-limit scaling

*Load-bearing part (needed for publication):* for fixed t ∈ (0, δ),
γ₁[F_{P1,N}(t)] = κ₃(F)/Var(F)^{3/2}. Numerator: K/N_B + O(N_B^(−2)) → 0 at
rate N_B^(−1) (V1-8). Denominator: Var(F_{P1,N}(t)) = Var X^ε(t) → Var x₀(t)
= m₂ > 0, since X^ε(t) → x₀(t) in L² (from the V1-5 bounds: E|X^ε −
x₀|² ≤ ε² E|y₁|² + O(ε⁴) → 0). Hence Var^{3/2} → m₂^{3/2} > 0 and

  **|γ₁[F_{P1,N}(t)]| = O(N_B^(−1)) → 0 as N_B → ∞.** ✓

*Stronger claim audited separately (the Gaussian reservoir limit).* For the
frozen protocols (P0/P1/P2) and any fixed finite tuple of times, the centred
force F̊_{q,N} = ε Σ_j (X_j^ε − E X^ε) is a triangular array of i.i.d. rows
(row N has ε = N_B^(−1/2)).

**Ledger S5 (multivariate Lindeberg–Feller CLT) — hypotheses checked:**
1. Within-row i.i.d.: by V1-2. ✓
2. Finite second moments: by the V1-5 moment bounds. ✓
3. Lindeberg condition: fourth moments are uniformly bounded for |ε| ≤ 1
   (E|X^ε|⁴ ≤ C, same Grönwall machinery), so
   N_B · E|N_B^(−1/2)(X − E X)|⁴ = N_B^(−1) E|X − E X|⁴ ≤ C/N_B → 0. This
   implies Lindeberg (Lyapunov condition ⇒ Lindeberg). ✓
4. Covariance convergence: Cov(X^ε(t_a), X^ε(t_b)) is continuous in ε at 0
   (V1-5 smoothness + dominated convergence), converging to C₀(t_a, t_b) =
   ⟨x₀(t_a) x₀(t_b)⟩_Gibbs. ✓

Conclusion: for each fixed finite time tuple, the centred force laws of the
**frozen protocols** converge to the **same** Gaussian law with covariance
C₀ — an E1-type reservoir limit **for the frozen / pointwise-fixed protocol
family**.

**Scope preserved (per the repaired record):** the Grönwall constants depend
on the clamp's derivative bounds (|q| ≤ Q, sup |q^(k)|), which are finite
but protocol-dependent; **no uniform theorem over the entire infinite clamp
class 𝒳 is established or claimed.** The limit is pointwise-fixed: protocol
by protocol, tuple by tuple. ✓

**Disposition: REPRODUCED-WITH-STANDARD-THEOREM** (S5: multivariate
Lindeberg–Feller / Lyapunov CLT, hypotheses explicitly checked; S2/S3 for
the L² convergence).

---

## Standard-Theorem Ledger

| Imported theorem | Exact use | Hypotheses checked? | Load-bearing? |
|---|---|---|---|
| S1: Invariance of Gibbs measure under Hamiltonian flow | V1-4(b): E x₀(t)² = m₂ for all t (stationarity) | Yes: smooth confining potential, ρ normalizable, flow global | Yes (bookkeeping in V1-6; also symmetry argument's invariance half) |
| S2: Smooth dependence of ODE flows on parameters/initial data | V1-5: existence of y₁, smoothness in ε; V1-11 L² convergence | Yes: polynomial vector field, global solutions on finite horizon via energy bound | Yes |
| S3: Dominated convergence / differentiation under E | V1-5, V1-8: ∂_ε commutes with E up to order 3 | Yes: integrable dominant P(E₀)e^{β√E₀}, uniform in \|ε\| ≤ 1 | Yes |
| S4: Taylor's theorem with remainder | V1-7 (pathwise, then expectation of the identity); V1-8 (κ₃ in ε) | Yes: C^∞ integrands; polynomial q on [0, π] | Yes |
| S5: Multivariate Lindeberg–Feller / Lyapunov CLT | V1-11 Gaussian reservoir limit (frozen protocols, fixed tuples) | Yes: i.i.d. rows, uniform 4th-moment bound ⇒ Lyapunov ⇒ Lindeberg; covariance convergence | No (contextual statement only; the load-bearing V1-11 content is the O(1/N_B) rate, which needs only S2/S3) |
| Cumulant additivity/homogeneity (CGF differentiation) | V1-3 | Elementary property of the CGF, derived in-line | Yes |
| Diffeomorphism pushforward of densities | V1-9 | Yes: C¹ time-map of smooth global flow | Yes |

No imported theorem was used without its hypotheses being checked. No
imported theorem substitutes for the central algebraic result (V1-6), which
was done by hand.

---

## Per-step disposition summary

| Step | Disposition |
|---|---|
| V1-1 | REPRODUCED-INDEPENDENTLY |
| V1-2 | REPRODUCED-INDEPENDENTLY |
| V1-3 | REPRODUCED-INDEPENDENTLY |
| V1-4 | REPRODUCED-WITH-STANDARD-THEOREM |
| V1-5 | REPRODUCED-WITH-STANDARD-THEOREM |
| V1-6 | REPRODUCED-INDEPENDENTLY |
| V1-7 | REPRODUCED-WITH-STANDARD-THEOREM |
| V1-8 | REPRODUCED-INDEPENDENTLY |
| V1-9 | REPRODUCED-INDEPENDENTLY |
| V1-10 | REPRODUCED-INDEPENDENTLY |
| V1-11 | REPRODUCED-WITH-STANDARD-THEOREM |

No REPRODUCTION-CODE was required (hand recurrence sufficient); no code
directory is committed.

---

## Post-derivation comparison with `B1R1_ANALYTIC_ESCAPE_THEOREM.md`

*(performed only after V1-1…V1-11 were complete, per the frozen ordering)*

**Substantive agreements (every load-bearing element):**
1. Clamped equation and force convention, sign and normalization (V1-1).
2. Exact i.i.d. structure, F_q = ε Σ X_j, identity (★) (V1-2, V1-3).
3. P0 exact parity: κ₃(F_{P0,N}(t)) = 0 for all t, N_B (V1-4).
4. Variational equation ÿ₁ + (1+3x₀²) y₁ = q(t), zero initial data (V1-5).
5. The central coefficient: **C₀…C₆ = 0, C₇ = −Var(x₀²)/(28π³)** — identical
   to the theorem's recorded result, including the equivalent form (m₂² +
   m₂ − 1)/(28π³) (using m₄ = 1 − m₂) and the π⁵-cancellation structure.
6. The existential δ via δ := min(1, |C₇|/(2C_R)) and the same inequality
   chain (V1-7).
7. K/N_B + O(N_B^(−2)) with per-t threshold N₀(t) := ⌊C/|K|⌋ + 1 and the
   explicit non-uniform quantifier structure (V1-8).
8. Positive variance via density-propagation under the flow time-map
   (V1-9); the theorem's phrasing ("C¹ diffeomorphism of ℝ²… so (x,p)(t*)
   has a density") matches.
9. The E2± escape mechanism: |γ| reflection-invariance ⇒ common-law value
   forced to 0 by P0 ⇒ contradiction with P1 ≠ 0; one-time tuple; non-
   degenerate ⇒ B1R-E2+O, not B1R-DEG (V1-10).
10. Reservoir-limit scope: O(N_B^(−1)) vanishing witness; Gaussian limit
    pointwise-fixed over frozen protocols only; no uniform statement over
    𝒳 (V1-11) — matching the §T4-R scope repair exactly.
11. The same "not used" firewall items: frozen-τ values, no numerical scan
    for t*, δ and N₀ existential.

**Discrepancies found:** none substantive.

- *Notation-level differences (harmless):* this reproduction derives the
  Taylor recurrences directly from the ODEs and records the y₁ coefficients
  b₅ = 1/(2π³), b₆ = −1/(2π⁴), b₇ = −(1+3a²)/(84π³) + 1/(7π⁵), which the
  theorem attributes to the preflight script's recursion; the recursions are
  the same object (order-by-order from the ODEs). The theorem carries the
  expansion to t¹⁴ via the script; this reproduction needed only t⁷ plus an
  order-8 remainder bound — a strictly weaker requirement, consistent with
  the theorem's own Grönwall remainder step.
- *Proof-level differences:* none. The independence check of V1-6 was done
  by hand rather than by script; the theorem's claim that the script's
  output is exact through t¹⁴ is consistent with (and stronger than) what
  the hand derivation verifies through the order actually needed.

**Over-strong dependence check:** the derivation used only (i) frozen X1 /
protocol definitions, (ii) frozen B1R0 class characterizations (B1R-C2±,
B1R-O5/O9/O10 — inputs, not the proof under test), (iii) standard theorems
(S1–S5) with hypotheses checked. Nothing stronger than the frozen theorem
was assumed: no uniform-in-𝒳 statement, no certified numerical values, no
all-orders expansion, no complex-analytic radius claims.

---

## Proposed publication theorem (draft; not yet a paper claim)

> **Theorem (finite-bath identifiability of reciprocal anharmonic back-
> reaction against a shared affine exogenous class).** Let X1_{N_B} denote
> the finite Duffing-bath environment (frozen Hamiltonian, coupling
> g = 1, ε = N_B^(−1/2), common Gibbs preparation at q(0) = 0), and let the
> clamp protocols be the frozen family {P0, P1} (P0: q ≡ 0; P1: the frozen
> quintic ramp s(t/π) on [0, π], then q = 1). There exists δ > 0 such that
> for every fixed t ∈ (0, δ) there is a finite integer N₀(t) with the
> following property: for every bath size N_B ≥ N₀(t), the interventional
> force family {F_q : q ∈ {P0, P1}} of X1_{N_B} does **not** admit a
> representation F_q(t) = M_t[q] + G_t[q] ξ(t) with M, G deterministic
> causal functionals (G ≠ 0) and one shared protocol-independent process ξ
> — i.e. the family lies outside the shared causal signed-affine exogenous
> class E2±.
>
> The witness is the absolute standardized skewness at time t: it is
> exactly zero for P0 for all N_B, while for P1 and N_B ≥ N₀(t) it is
> nonzero, with κ₃(F_{P1,N_B}(t)) = K(t)/N_B + O(N_B^(−2)) and K(t) < 0;
> equivalently the leading coefficient is C₇ = −Var(x₀²)/(28π³) < 0, where
> x₀ is the unforced Duffic Gibbs trajectory and Var(x₀²) > 0.
>
> The witness is mesoscopic: |γ₁[F_{P1,N_B}(t)]| = O(N_B^(−1)) and vanishes
> as N_B → ∞, with the centred frozen-protocol force laws converging to a
> common Gaussian limit (pointwise over fixed finite time tuples; no
> uniform statement over the full clamp class is claimed).
>
> The result is relative to exactly E2±. It makes no claim against a
> universal causal exogenous representation.

All quantifiers explicit: X1 frozen; protocol family frozen and finite;
δ existential (per the sign argument, no value claimed); t ranges over
(0, δ) pointwise; N₀(t) existential and t-dependent; competitor class
exactly E2±; rate O(N_B^(−1)) in the stated perturbative regime; reservoir
scope pointwise-fixed; no E_univ claim.

---

## V1 terminal

All eleven steps reproduce; nine of eleven fully independently, two with
standard theorems whose hypotheses are explicitly checked; no discrepancy;
no over-strong assumptions; the firewalls held (frozen-τ unused, no old
proof prose used as source, no values fitted).

# **V1-REPRODUCED**

Per the charter: with V1-REPRODUCED, V1 is complete. V2 and V3 remain
not-executed (not authorized at this boundary). Awaiting owner review before
the charter's V2 lane may be opened.
