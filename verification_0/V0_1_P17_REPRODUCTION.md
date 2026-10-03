# V0-1 — Blind reproduction of P-17 (harmonic-bath elimination ⇒ equality of reduced path laws)

**Inputs used:** only `verification_0/specs/V0_1_P17_TARGET.md` (the spec), plus general textbook knowledge
(marked KNOWN). No other project file was opened. Everything below was derived and checked fresh.

**Code / logs (all new):**
- `code/v0_1_taylor_checks.py` — exact rational Taylor checks T1–T8. Logs: `code/v0_1_taylor_checks_2modes.log`
  (two modes, n ≤ 8) and `code/v0_1_taylor_checks_1mode.log` (one mode, n ≤ 12).
- `code/v0_1_elimination_checks.py` — symbolic checks of steps (a)–(d), (g). Log: `code/v0_1_elimination_checks.log`.

**Verdict in brief:** every acceptance test T1–T8 PASSES, and every expected expression in the spec is reproduced
exactly. The theorem holds as stated. Section 4 lists some clarifications and limits on scope. None of them is a
numerical error in the spec, but (P) depends on how F is split off, the identifiability claim needs side conditions,
and the B → A direction needs 𝒫 to be supported on the bath's trigonometric span.

---

## 1. Derivation

Notation: q ∈ ℝⁿ, M is symmetric positive definite, and V is smooth and bounded below (confining). Bath modes are
j = 1..m with ω_j > 0 and c_j ∈ ℝⁿ.

  H = ½pᵀM⁻¹p + V(q) + Σ_j [ p_j²/2 + ½ω_j²(x_j − c_jᵀq/ω_j²)² ].

Hamilton's equations are:

  q̇ = M⁻¹p,  ṗ = −∇V(q) + Σ_j c_j (x_j − c_jᵀq/ω_j²),
  ẋ_j = p_j,  ṗ_j = −ω_j² x_j + c_jᵀq.

(The counterterm contributes −Σ c_j c_jᵀ q/ω_j² to ṗ, and the bath cross term contributes +Σ c_j x_j.)

**Global well-posedness of A (KNOWN).** The vector field is smooth, hence locally Lipschitz. H is conserved and
bounded below. The bath part is ≥ 0 and V is bounded below. So every coordinate stays bounded on bounded time
intervals, and the solution is global and unique.

### (a) Exact elimination of the bath

Fix the system path q(·). Each bath mode is a driven harmonic oscillator, ẍ_j + ω_j² x_j = c_jᵀq(t). Variation of
constants gives (KNOWN):

  x_j(t) = x_j(0) cos ω_j t + (p_j(0)/ω_j) sin ω_j t + (1/ω_j) ∫₀ᵗ sin ω_j(t−s) c_jᵀq(s) ds.

`v0_1_elimination_checks.log`, item (a), confirms that this satisfies the ODE and both initial conditions for an
arbitrary function q.

**Integration by parts.** Use d/ds[cos ω(t−s)] = ω sin ω(t−s). Then

  (1/ω)∫₀ᵗ sin ω(t−s) q(s) ds = (1/ω²)[q(t) − q(0) cos ωt] − (1/ω²)∫₀ᵗ cos ω(t−s) q̇(s) ds.

Item (b) of the elimination log checks this exactly for three concrete paths. It follows that

  x_j(t) − c_jᵀq(t)/ω_j² = ξ_j cos ω_j t + (η_j/ω_j) sin ω_j t − (1/ω_j²)∫₀ᵗ cos ω_j(t−s) c_jᵀq̇(s) ds,

  with **ξ_j := x_j(0) − c_jᵀq₀/ω_j²** (the bath displacement measured from its q₀-shifted equilibrium) and
  **η_j := p_j(0)**.

The counterterm −c_jᵀq(t)/ω_j² exactly cancels the instantaneous term (1/ω_j²)c_jᵀq(t) produced by the
integration by parts. That is its whole role: with the counterterm, no renormalisation of the potential
(no −Σ c_j c_jᵀ q/ω_j² shift) survives. Without it, the system would feel V_eff = V − ½Σ (c_jᵀq)²/ω_j².

Substituting into ṗ gives the closed equation, exactly and pathwise:

  **M q̈ = −∇V(q) − ∫₀ᵗ γ(t−s) q̇(s) ds + F(t),**

### (b) Memory kernel

  **γ(t) = Σ_j (c_j c_jᵀ/ω_j²) cos ω_j t.**

This matches the spec's γ. It is even in t, and γ(0) = Σ c_j c_jᵀ/ω_j² is exactly the counterterm matrix. The kernel
is a property of the Hamiltonian alone and does not depend on the preparation.

### (c) Free force

  **F(t) = Σ_j c_j [ ξ_j cos ω_j t + (η_j/ω_j) sin ω_j t ],  with ξ_j = x_j(0) − c_jᵀq₀/ω_j² and η_j = p_j(0).**

Write this as F = Φ(ξ, η). Φ is linear. It depends on the bath initial data, and on q₀ only through the shift inside
ξ. It never depends on p₀. When the ω_j are distinct and the c_j ≠ 0, the functions {c_j cos ω_j t, c_j sin ω_j t}
are linearly independent, so Φ is injective. (KNOWN: distinct-frequency trigonometric functions are linearly
independent.)

### (d) The preparation condition and the initial slip

From (c), the law of F given (q₀, p₀) is Φ_#(law of (ξ, η) given (q₀, p₀)). So:

**(P) holds ⇔ Φ_#Law(ξ, η | q₀, p₀) does not depend on (q₀, p₀).** When Φ is injective (distinct frequencies), this
is the same as: the conditional law of the *shifted* bath data (x_j(0) − c_jᵀq₀/ω_j², p_j(0)) does not depend on
(q₀, p₀). With degenerate frequencies, only the law of the relevant linear combinations has to be invariant.

Natural example: the joint Gibbs measure exp(−H/T) conditioned on (q₀, p₀) is exactly the bath Gibbs measure
centred at x_j = c_jᵀq₀/ω_j². This holds because the bath part of H is ½Σ[p_j² + ω_j²ξ_j²]. This is the
"shifted-Gibbs" preparation, and it satisfies (P).

**Product preparation.** Take x_j(0) ~ N(0, T/ω_j²) and p_j(0) ~ N(0, T), independent of q₀. Then
ξ_j = x_j(0) − c_jᵀq₀/ω_j² ~ N(−c_jᵀq₀/ω_j², T/ω_j²). Write ξ_j = ξ_j^th − c_jᵀq₀/ω_j² with ξ_j^th ~ N(0, T/ω_j²).
Then

  F(t) = F^th(t) − Σ_j c_j (c_jᵀq₀/ω_j²) cos ω_j t = **F^th(t) − γ(t) q₀**,

where F^th has the thermal law (mean 0, covariance Tγ). The law of F depends on q₀ through its mean, so **(P) fails**
whenever γ(·)q₀ ≢ 0. The spec's slip term **−γ(t)q₀ is CONFIRMED, exactly and pathwise** (elimination log, item (d)).
Under the product preparation, A therefore equals the GLE driven by an exogenous thermal F plus the deterministic
forcing −γ(t)q₀.

*Remark (adversarial).* Integrating by parts once more,

  −∫₀ᵗγ(t−s)q̇ ds − γ(t)q₀ = −γ(0)q(t) − ∫₀ᵗ γ̇(t−s) q(s) ds.

So under the product preparation, A is exactly a GLE written in the "q-kernel" form, with a renormalised potential
term −γ(0)q, driven by F′ = Σ c_j[x_j(0) cos + (p_j(0)/ω_j) sin], whose law *is* q₀-independent. "(P) fails" is
therefore a statement relative to the split into a q̇-memory term plus the free force F. It is not intrinsic to the
dynamics. The spec's wording ("A equals B plus a deterministic initial slip") is consistent with this.

### (e) Equality of reduced path laws (A ≡ B)

**Regularity needed (R).** Fix t_max. For every (q₀, p₀) and every forcing path f in a set S ⊂ C([0, t_max]; ℝⁿ)
carrying 𝒫, the GLE M q̈ = −∇V(q) − ∫γ(t−s)q̇ ds + f(t), q(0) = q₀, q̇(0) = M⁻¹p₀ has a unique C² solution
q = G_{q₀,p₀}(f). The map f ↦ G_{q₀,p₀}(f) must also be Borel measurable from S to C([0, t_max]).

(R) holds here, for all continuous f:
- *Local existence and uniqueness:* this is a Volterra integro-differential equation with a smooth bounded kernel and
  a locally Lipschitz drift (standard Picard argument, KNOWN).
- *No blow-up:* the energy E = ½q̇ᵀMq̇ + V(q) satisfies dE/dt = q̇ᵀ(f − ∫γq̇). With ‖γ‖_∞ < ∞ and f bounded on
  [0, t_max], Gronwall bounds |q̇|, and V bounded below then bounds E.
- *Measurability:* Gronwall stability makes G continuous in f under the sup norm, hence Borel.

For the bath-generated f = Φ(ξ, η), existence is immediate in any case, because the Hamiltonian trajectory supplies a
solution.

**Proof.** (1) By (a)–(c), for every bath initial datum, the system component of the (unique) Hamiltonian solution
solves the GLE with f = Φ(ξ, η). By uniqueness in (R), q_A = G_{q₀,p₀}(F_A) pathwise, with F_A = Φ(ξ, η).
(2) In B, q_B = G_{q₀,p₀}(F_B) with F_B ~ 𝒫 independent of (q₀, p₀). (3) By (P), Law(F_A | q₀, p₀) = 𝒫 = Law(F_B).
The same measurable map G_{q₀,p₀} is applied to equal laws, so the push-forward laws on C([0, t_max]) coincide:
Law(q_A(·) | q₀, p₀) = (G_{q₀,p₀})_#𝒫 = Law(q_B(·) | q₀, p₀).

If (q₀, p₀) is itself random, (P) in its conditional form means F_A is independent of (q₀, p₀). The joint laws of
(q₀, p₀, q(·)) then also agree. Every reduced statistic (any measurable functional of the system path, at any number
of times) therefore agrees. The p-path agrees too, since p = Mq̇. ∎

(iii) **Deterministic bath microstate.** If ξ, η are deterministic and q₀-independent (that is, x_j(0) = c_jᵀq₀/ω_j²
+ ξ_j* with ξ* fixed), then 𝒫 = δ_{Φ(ξ*, η*)} and q_A = q_B pathwise. CONFIRMED.

### (f) Ontology / identifiability conclusion and its scope

**Statement.** Within 𝓗, assume (P), with M and V known. Suppose the reduced data consist of the laws of the system
path q(·) on [0, t_max] for all (q₀, p₀). Then the data determine the pair (γ, 𝒫) and nothing more about where the
forcing comes from. Both models below give identical reduced data for every (q₀, p₀):
- a Hamiltonian bath whose randomness is only ignorance of its initial microstate, distributed so that Φ_#Law = 𝒫;
- an exogenous, "primitively random" process with law 𝒫.

The same holds for any two bath microstate laws with the same image 𝒫, for example Gaussian versus non-Gaussian
microstate laws, provided their Φ-images agree. For a Dirac 𝒫, a hidden deterministic microstate and a deterministic
external drive are likewise indistinguishable.

**Why (γ, 𝒫) is identified, and what this requires.**
- *𝒫 given γ:* on each path, F(t) = Mq̈ + ∇V(q) + ∫γq̇, so F is recovered path by path and its law is determined.
- *γ:* suppose (γ, 𝒫) and (γ′, 𝒫′) produce the same path laws. Then F′ = F + ∫₀ᵗ Δγ(t−s) q̇(s) ds, with
  Δγ = γ′ − γ, must have a (q₀, p₀)-independent law. Differentiate its mean in p₀ (this needs E|q̇| < ∞ and
  differentiation under E, both fine here). The result is ∫₀ᵗ Δγ(t−s) K(s) ds ≡ 0, where K(s) = E[∂q̇(s)/∂p₀] is
  continuous and differentiable with K(0) = M⁻¹. This is a first-kind Volterra equation with an invertible kernel
  at 0. Differentiating gives a second-kind equation whose only solution is Δγ ≡ 0 (KNOWN).

  The spec's assertion that the data identify "the memory kernel" therefore needs the reduced data to include
  variation over initial momenta (or something equivalent) and M, V to be known. It is not true for a single fixed
  initial condition. **This is a qualification the spec leaves implicit.**

**Scope.**
1. Class 𝓗 only: finitely many harmonic modes, coupling linear in x_j and in q, with the counterterm; classical.
2. Preparation (P) is required. The product preparation is covered only after the deterministic slip is added.
3. Converse direction B → A ("primitive randomness could have been a bath"): within 𝓗 with kernel γ, F_A always lies
   in the 2m-dimensional span V_γ = span{c_j cos ω_j t, c_j sin ω_j t}. So an exogenous 𝒫 is bath-realisable with
   that γ **iff 𝒫 is supported on V_γ**. The law of (ξ, η) is then Φ⁻¹_#𝒫, using Φ injective (otherwise pick any
   preimage law).
   - The Gaussian law with covariance Tγ *is* supported on V_γ: it has rank 2m.
   - White noise, or any F not in V_γ, is **not** realisable by an 𝓗-bath with that γ.

   The non-identifiability is two-sided only for 𝒫 on V_γ, or in a many-mode / continuum limit, which lies outside the
   finite-mode 𝓗.
4. Reduced data means system observables only. Any access to bath variables, or any intervention on the bath, breaks
   the equivalence.
5. Claim (v) (the C-B discriminator) relies on a model definition that is not in the spec. It **cannot be verified
   here and is not reproduced.** My derivation supports only the generic part: reduced statistics depend on the
   forcing through its law, i.e. its cumulant hierarchy, and not on its ontological origin.

### (g) Thermal shifted-Gibbs law of F

ξ_j ~ N(0, T/ω_j²) and η_j ~ N(0, T), all independent. Φ is linear, so F is a centred Gaussian process with

  E[F(t)F(s)ᵀ] = Σ_j c_j c_jᵀ [ (T/ω_j²) cos ω_j t cos ω_j s + (T/ω_j²) sin ω_j t sin ω_j s ] = T Σ_j (c_j c_jᵀ/ω_j²) cos ω_j(t−s) = **T γ(t−s)**.

So mean = 0 and covariance = Tγ(t−s), which is stationary. This is the classical fluctuation–dissipation relation.
CONFIRMED (elimination log, item (g)).

"Matched covariance suffices" is true **only among Gaussian exogenous forcings**: a Gaussian law is fixed by its mean
and covariance (KNOWN).

The random-phase preparation gives F_j(t) = c_j √(2T)/ω_j cos(ω_j t + φ_j). This has the same mean (0) and the same
covariance (Tγ), but it is not Gaussian. Its fourth cumulant at equal times, per mode, is −(3/2)T²c_j⁴/ω_j⁴ (log,
item (g)). Matched covariance therefore does **not** suffice for a non-Gaussian 𝒫. T4/T6 show this explicitly:
one mode, n = 12.

### (h) Hostile scope checks (analytic, short)

1. **Non-harmonic bath** (for example ½ω²y² + λy⁴ with y = x − cq/ω²). The bath equation
   ÿ + ω²y + 4λy³ = −(c/ω²)q̈ is nonlinear, so y is no longer an affine functional "free part + linear response to q".
   The linear response of a mode about its own orbit has an amplitude-dependent frequency. The "memory kernel" then
   becomes random and is correlated with the "noise", and the force does not split into a q-independent additive
   F plus a deterministic linear memory term.
   - **The theorem in form B fails** in general.
   - What survives is the trivial abstract statement: q_A = 𝒢_{q₀,p₀}(bath data), with bath data independent of
     (q₀, p₀). So *some* exogenous-randomness model reproduces A (the bath data themselves act as the exogenous
     input), and the determinism-versus-randomness non-identifiability persists in that weak sense.
   - The specific GLE structure and the "(γ, 𝒫) only" conclusion do not survive.

2. **Coupling nonlinear in q** (x_j coupled to A_j(q), with counterterm ½ω_j²(x_j − c_jA_j(q)/ω_j²)²). The bath is
   still linear and the elimination goes through. With F_j(t) = c_j[(x_j(0) − c_jA_j(q₀)/ω_j²) cos ω_j t + (p_j(0)/ω_j) sin ω_j t],
   the reduced equation is

   M q̈ = −∇V − Σ_j ∇A_j(q(t)) [ ∫₀ᵗ (c_j²/ω_j²) cos ω_j(t−s) (d/ds)A_j(q(s)) ds − F_j(t) ].

   - The noise becomes **multiplicative**, the memory becomes state-dependent, and the forcing has to be described as
     a vector of channel processes (F_j), not their sum.
   - The theorem **survives with B replaced by this multiplicative GLE**. (P) becomes: the law of
     (x_j(0) − c_jA_j(q₀)/ω_j², p_j(0)) is (q₀, p₀)-independent.
   - The spec's additive B with kernel γ **fails** for this class. The spec excludes it by defining 𝓗 as linear in q.

3. **q₀-dependent preparation.** If Law(ξ, η | q₀) depends on q₀, then A = GLE driven by F ~ 𝒫_{q₀}, and the exact
   identity holds with B_{q₀}. The "exogenous, independent of (q₀, p₀)" form of B fails.
   - The violation is absorbable into deterministic memory only when the dependence is a deterministic, q₀-linear
     shift of the mean. Explicitly, ξ = ξ⁰ + Lq₀ with Law(ξ⁰) fixed gives F = F⁰ + [Σ c_j (L_j q₀) cos ω_j t + …],
     a deterministic forcing linear in q₀. For L_j = −c_jᵀ/ω_j² this is exactly the product slip −γ(t)q₀.
   - A q₀-dependent temperature, covariance or higher cumulant is **not** absorbable. The forcing law then encodes
     q₀, and no q₀-independent exogenous 𝒫 reproduces A.
   - p₀-dependence of the preparation behaves the same way, even though p₀ never enters Φ.

4. **Extra restrictions on the exogenous force.** Exactness requires Law(F_B) = 𝒫 *exactly* (all cumulants).
   - *Gaussian F_B* while 𝒫 is non-Gaussian (random phase, Dirac): this fails. Demonstrated by T4/T6, one mode,
     n = 12: E_A − E_R = +138510 T²aβ².
   - *Markov / white F_B:* Tγ is a finite cosine sum, not δ, so white forcing fails. F_A is Markov only as the output
     of the 2m-dimensional rotating state (ξ, η). It is not Markov as an n-dimensional process.
   - *Fixed covariance C ≠ Tγ* (Gaussian, mean 0): the leading discrepancy is already at n = 6. By the same
     computation as in T3, the coefficient changes to E_A[q⁽⁶⁾] − E_{B′}[q⁽⁶⁾] = 6!·[−βa(Tγ(0) − C(0,0))/(10M³)],
     which is non-zero whenever C(0,0) ≠ Tγ(0) and βa ≠ 0. If C(0,0) = Tγ(0), the mismatch shows up at higher n,
     through derivatives of C at 0.

**T3 derivation (analytic).** Write q = a + δ. Then δ = F(0)t²/(2M) + O(t³), so E[δ²] = Tγ(0)t⁴/(4M²) + O(t⁵).
The term −4βq³ contributes −12βaE[δ²] to ME[q̈] − ME[q̈_D] at leading order, which is −3βaTγ(0)t⁴/M². Integrating
twice and dividing by M gives E[q − q_D] = −βaTγ(0)t⁶/(10M³). Hence the leading coefficient is
6!·[−βaTγ(0)/(10M³)], CONFIRMED. With γ(0) = 25/64 (two modes) this gives −225Taβ/8; with γ(0) = 1/4 (one mode)
it gives −18Taβ.

---

## 2. Computational method (code designed independently)

Model as in the spec: n = 1, M = 1, V = (23/10)q²/2 + βq⁴, so V′ = (23/10)q + 4βq³. Modes are (1, 1/2) and
(2, 3/4). Initial data: q(0) = a, p(0) = 0. Statistic: E[q⁽ⁿ⁾(0)] = n!·E[q_n], where q(t) = Σq_k t^k.

All arithmetic is exact, in sparse multivariate polynomial rings over ℚ (`sympy.polys.rings`, FLINT-backed `QQ`).
The recursions are linear in the order, and only the cubic q³ needs a Cauchy product, (q³)_m = Σ_{i+j+k=m} q_i q_j q_k.
There is no symbolic differentiation and no expression swell; runs take < 1 s.

- **A, R, product (full Hamiltonian, all 2 + 2m phase-space variables).** Power-series recursion
  z_{m+1} = (rhs)_m/(m+1) for (q, p, x_j, p_j). Coefficients lie in ℚ[a, β, ξ_j, η_j]:
  - shifted preparations use x_j(0) = c_j a/ω_j² + ξ_j;
  - the product preparation uses x_j(0) = ξ_j.

  E is then taken monomial by monomial:
  - Thermal / product: E[ξ^i] = (i−1)!!(T/ω²)^{i/2} and E[η^k] = (k−1)!!T^{k/2} (zero if i or k is odd).
  - Random phase: E[ξ^i η^k] = (2T)^{(i+k)/2} ω^{−i} (−1)^k E[cos^iφ sin^kφ], with
    E[cos^iφ sin^kφ] = (i−1)!!(k−1)!!/(i+k)!! for i, k both even, and 0 otherwise.
- **B (Gaussian exogenous F), D, B + slip: the GLE solved directly**, with no bath variables. Recursion:

  q_{m+2} = [−(23/10)q_m − 4β(q³)_m − mem_m + F_m + S_m]/((m+1)(m+2)),
  mem_m = Σ_{k+l+1=m} g_k (l+1)q_{l+1} · k!l!/m!,

  where g_k are the Taylor coefficients of γ and the factor comes from ∫₀ᵗ(t−s)^k s^l ds = k!l!t^{k+l+1}/(k+l+1)!.
  - F is kept symbolic through f_k = F⁽ᵏ⁾(0), with F_m = f_m/m!.
  - E is computed by **Wick/Isserlis pairing** (memoised on sorted index multisets), using
    Cov(f_k, f_l) = T(−1)^l γ⁽ᵏ⁺ˡ⁾(0). This follows from the stationary covariance Tγ(t−s) alone, with no reference
    to the bath.
  - D sets F = 0. "B + slip" adds S(t) = −γ(t)a.
- **B_R (exogenous random-phase F):** the same GLE recursion with
  F(t) = Σ c_j/ω_j (u_j cos ω_j t − v_j sin ω_j t), u_j = √(2T) cos φ_j, v_j = √(2T) sin φ_j, φ_j i.i.d. uniform,
  and exact phase moments.

A (Hamiltonian, bath moments) and B (GLE, Wick from γ) share no intermediate objects. T1 is therefore a genuine
independent check of the elimination, the counterterm cancellation and the FDT covariance.

**Interpretation of T7 "product − B (shift-matched)":** product preparation = Hamiltonian with x_j(0) ~ N(0, T/ω_j²),
p_j(0) ~ N(0, T), q(0) = a. It is compared with B, the same GLE with thermal exogenous F, and with B + slip (B plus the
deterministic forcing −γ(t)a).

---

## 3. Acceptance-test results (verbatim from the logs)

### Two modes, n ≤ 8 (`v0_1_taylor_checks_2modes.log`; γ(0) = 25/64)

| test | obtained | spec | status |
|---|---|---|---|
| T1 A − B | 0 for n = 0..8 | 0 | **PASS** |
| T2 A − D | 0 for n ≤ 5 and n = 7; n=6: `-225*T*a*beta/8`; n=8: `15*T*a*beta*(63360*a**2*beta + 25973)/256` | identical | **PASS** |
| T3 | leading n = 6 coefficient `-225*T*a*beta/8` = 6!·(−βaTγ(0)/10) | identical | **PASS** |
| T4 A − R | 0 for n = 0..8 | 0 | **PASS** |
| T5 R − B_R | 0 for n = 0..8 | 0 | **PASS** |
| T6 B_R − B | 0 for n = 0..8 | 0 | **PASS** |
| T7 product − B | non-zero at n = 2, 4, 6, 8 (n=2: `-25*a/64` = −γ(0)a); 0 at odd n | non-zero at {2,4,6,8} | **PASS** |
| T7 product − (B + slip) | 0 for n = 0..8 | 0 | **PASS** |
| T8 | β=0 ⇒ A−D ≡ 0: True; T=0 ⇒ A−D ≡ 0: True | as stated | **PASS** |

T7 product − B, all non-zero entries (two modes):
- n=4: `a*(19200*a**2*beta + 7633)/4096`
- n=6: `-a*(368640000*a**4*beta**2 + 229816320*a**2*beta + 10264813)/1310720`
- n=8: `3*a*(23040000000*T*beta + 3303014400000*a**6*beta**3 + 3341613465600*a**4*beta**2 + 821672551680*a**2*beta + 4587466931)/419430400`

### One mode, n ≤ 12 (`v0_1_taylor_checks_1mode.log`; γ(0) = 1/4)

| test | obtained | spec | status |
|---|---|---|---|
| T1 A − B | 0 for n = 0..12 | 0 | **PASS** |
| T2 A − D | 0 for n ≤ 5 and all odd n; n=6: `-18*T*a*beta`; n=8: `24*T*a*beta*(495*a**2*beta + 193)/5`; n=10: `-3*T*a*beta*(27302400*a**4*beta**2 + 18202640*a**2*beta + 1976871)/200`; n=12: `3*T*a*beta*(46170000*T*beta + 15956352000*a**6*beta**3 + 15791424000*a**4*beta**2 + 4083013860*a**2*beta + 136763237)/500` | identical | **PASS** |
| T3 | leading n = 6: `-18*T*a*beta` = 6!·(−βaT/40) | identical | **PASS** |
| T4 A − R | 0 for n ≤ 11; n=12: `138510*T**2*a*beta**2` | identical | **PASS** |
| T5 R − B_R | 0 for n = 0..12 | 0 | **PASS** |
| T6 B_R − B | 0 for n ≤ 11; n=12: `-138510*T**2*a*beta**2` | "non-zero only at n = 12" | **PASS** (value obtained) |
| T7 product − B | non-zero at n = 2,4,…,12 (n=2: `-a/4` = −γ(0)a); 0 at odd n | non-zero at {2,…,12} | **PASS** |
| T7 product − (B + slip) | 0 for n = 0..12 | 0 | **PASS** |
| T8 | β=0 ⇒ 0: True; T=0 ⇒ 0: True | as stated | **PASS** |

Consistency: T4 + T5 + T6 close, since (A−R) + (R−B_R) + (B_R−B) = A − B = 0, i.e. 138510 + 0 − 138510 = 0. The
sign of A − R > 0 matches the negative fourth cumulant of the random-phase force. The fourth cumulant first enters at
order β²T²a t¹², via the chain δ₁ ~ Ft², δ₂ ~ βδ₁³t² ~ t⁸, and −12βa·2E[δ₁δ₂] ~ t¹⁰ in q̈.

---

## 4. Disagreements, corrections and clarifications relative to the spec

- **No numerical disagreement.** Every expected expression in T1–T8 was reproduced exactly.
- **(P) depends on the split** (§1(d) remark). Under the product preparation, (P) fails for the q̇-memory/free-force
  split. A q-kernel split (−γ(0)q − ∫γ̇(t−s)q(s)ds plus the unshifted free force) restores a q₀-independent forcing
  law exactly. The spec's (ii) is correct as stated, but its "(P) fails" refers to that split only.
- **"Reduced data identify the memory kernel"** (iv) needs side conditions: M and V known, and data across initial
  momenta (or equivalent). §1(f) gives a sufficient argument. For a single initial condition it is not established.
- **The non-identifiability is two-sided only on V_γ.** An exogenous 𝒫 is realisable by an 𝓗-bath with the same γ
  iff 𝒫 is supported on the 2m-dimensional trigonometric span. Gaussian Tγ, random phase and Dirac laws are; white
  noise or a generic external drive is not.
- **"Matched covariance suffices"** (i) is true only inside the Gaussian family. T4/T6 at one mode, n = 12, give an
  explicit counterexample for non-Gaussian 𝒫 with matched covariance.
- **Claim (v) (C-B)** cannot be checked from the spec alone and is not reproduced. The spec itself defers the C-B
  identities.
