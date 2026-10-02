# S2-D-arrow D0 — theorem-first: no every-state arrow in finite closed unitary dynamics

**Status:**
- **TRUE DERIVATION** (elementary; proof below);
- **KNOWN RESULT IMPORT** for the ingredients: quantum recurrence (Bocchieri–Loinger 1957), Dirichlet simultaneous
  approximation, Fannes continuity — STANDARD.

What SCOUT adds is the scoping and the consequences for the arrow wall.

## Setting (the scope; nothing outside it is claimed)

- ℋ is finite-dimensional, dim N.
- 𝒟 is the set of density matrices, with the trace norm ‖·‖₁.
- H is a Hermitian Hamiltonian, U_t = e^{−iHt}, and ρ(t) = U_t ρ U_t†.
- An **arrow functional** is any map A: 𝒟 → ℝ that is continuous in ‖·‖₁.
- No boundedness assumption is needed: continuity on the compact set 𝒟 already makes A bounded.

**Covered examples:**
- global and **subsystem** von Neumann entropies S(Tr_E ρ), which are continuous by Fannes;
- Rényi entropies;
- mutual informations;
- trace distance to a fixed state;
- **coarse-grained / observational / Boltzmann-type entropies**, i.e. the Shannon entropy of the outcome distribution
  of a fixed POVM (continuous in ρ);
- the S2-8 redundancy diagnostics built from these.

## Lemma 1 (recurrence)

For every ρ ∈ 𝒟, ε > 0 and T > 0 there is τ > T with ‖ρ(τ) − ρ‖₁ < ε.

*Proof.*
- Write H = Σ_j E_j P_j and fix a step s > T.
- By Dirichlet's simultaneous approximation theorem, for every integer Q ≥ 1 there is an integer 1 ≤ q ≤ Q^N such
  that for all j the distance of q·E_j s/2π from the nearest integer is ≤ 1/Q.
- Put τ = q s ≥ s > T. Then |e^{−iE_jτ} − 1| ≤ 2π/Q for all j, so ‖U_τ − 𝟙‖_op ≤ 2π/Q.
- Hence ‖ρ(τ) − ρ‖₁ ≤ 2‖U_τ − 𝟙‖_op ≤ 4π/Q. Choose Q > 4π/ε. ∎

## Theorem 1 (no orbit-wise eternal arrow)

Let A be continuous. If t ↦ A(ρ(t)) is monotone (non-decreasing or non-increasing) on [0, ∞) for some ρ, then it is
constant on [0, ∞).

*Proof.*
- Suppose A is non-decreasing and A(ρ(t₁)) = A(ρ) + δ with δ > 0.
- By continuity there is η > 0 with |A(σ) − A(ρ)| < δ/2 whenever ‖σ − ρ‖₁ < η.
- By Lemma 1 there is τ > t₁ with ‖ρ(τ) − ρ‖₁ < η, so A(ρ(τ)) < A(ρ) + δ/2 < A(ρ(t₁)).
- This contradicts monotonicity. The non-increasing case is symmetric. ∎

## Corollary 1 (every-state hostile)

If A is monotone along **every** orbit, then A is conserved: A(ρ(t)) = A(ρ) for all ρ and t.

Equivalently: for any A that is not conserved there is an **admissible state whose forward evolution strictly
decreases A**. To find it:
1. take any orbit on which A is nonconstant;
2. by Theorem 1 there are t₁ < t₂ with A(ρ(t₁)) > A(ρ(t₂));
3. the admissible initial state ρ′ = ρ(t₁) decreases over the window [0, t₂ − t₁].

→ **EVERY-STATE ARROW = IMPOSSIBLE** for every continuous functional: global, subsystem, coarse-grained or
record-based.

## Corollary 1′ (every state, any fixed nonzero horizon) — REPAIR 04 clarification

*Owner-requested clarification; not a new result.*

Suppose a continuous A were non-decreasing on one fixed interval [0, T], T > 0, for **every** admissible state.

1. Apply this to the admissible states U_{kT} ρ U_{kT}† for k = 0, 1, 2, …. Then A(ρ(t)) is non-decreasing on each
   [kT, (k+1)T].
2. The intervals share endpoints, so A(ρ(t)) is non-decreasing on [0, ∞).
3. By Theorem 1, A(ρ(t)) is constant for every ρ, i.e. A is conserved.

The same argument works for non-increasing. Hence, in finite closed unitary dynamics:

> **EVERY-STATE ARROW OVER ANY FIXED NONZERO HORIZON = IMPOSSIBLE** for any non-conserved continuous functional.

## Corollary 2 (window mirror under a Σ-respecting time reversal)

Suppose there is an antiunitary Θ with ΘHΘ⁻¹ = H and A(ΘρΘ⁻¹) = A(ρ). Then for every ρ and t, the state
ρ′ = Θρ(t)Θ⁻¹ satisfies ρ′(s) = Θρ(t − s)Θ⁻¹, so

    A(ρ′(s)) = A(ρ(t − s))   for 0 ≤ s ≤ t.

So **every increase of A over a window [0, t] is matched by an equal decrease, over the same window, from another
admissible state.**

**Σ-dependence.** For subsystem functionals, A(ΘρΘ⁻¹) = A(ρ) holds when Θ is a **product of local antiunitaries in the
TPS that defines A**, e.g. complex conjugation in a product basis times local unitaries. The mirror therefore depends on
Θ being compatible with Σ. A time-reversal symmetry that is non-local in Σ gives no subsystem mirror. This is an
explicit **Σ-coupling of the arrow question**.

## Proposition 2 (a state fixed by H alone cannot carry an arrow)

Let F map Hamiltonians to states **unitary-covariantly**: F(VHV†) = V F(H) V† for all unitaries V. Then F(H) is
stationary under H.

*Proof.*
- Take V_s = e^{−iHs}. Then V_s H V_s† = H, so F(H) = V_s F(H) V_s† for all s.
- Differentiate at s = 0: [H, F(H)] = 0. ∎

**Consequence.** No arrow can start from a state fixed by the law H alone:
- ground states;
- Gibbs states;
- any functional-calculus state.

A *structurally* fixed special state **with dynamics** needs extra covariant structure. The minimal case is Σ. For
example, "the product state (relative to Σ) of minimal energy" is fixed by (H, Σ) and is generically non-stationary.
This sharpens D10 as follows:

- a boundary condition fixed by **(H, Σ)** is possible;
- one fixed by **H alone** is impossible.

This is consistent with REPAIR 02: H-only covariant selectors see at most commutant classes.

## Consequences (D14 theorem target, in the proved scope)

> In a finite closed unitary universe, with a fixed TPS:
> 1. exact global information is preserved (‖ρ_a(t) − ρ_b(t)‖₁ constant);
> 2. no nonconstant continuous arrow functional — subsystem entropy, coarse-grained entropy, mutual information or a
>    record diagnostic — is monotone along any single orbit forever (Theorem 1). None that is non-conserved is monotone on a fixed nonzero
>    horizon [0, T] for every admissible state (Corollary 1′, by concatenating U_{kT}-shifted windows);
> 3. if the dynamics has a Σ-local time-reversal symmetry, every arrow window has an exactly mirrored anti-arrow window
>    (Corollary 2).
>
> Effective subsystem arrows therefore require **restricting the initial global states** (and/or correlations), the
> **observables** (coarse-graining, which still does not escape Theorem 1), the **time window** (finite horizon,
> shorter than the recurrence times), or the **factorization** used to define them.

**Outside scope (not claimed):**
- infinite-dimensional systems with continuous spectrum, where Lemma 1 fails (scattering, infinite baths);
- idealized Markov semigroups (limits of the above);
- non-unitary or non-compact invertible flows;
- discontinuous functionals.

These are exactly the places where effective arrows are usually manufactured. The numerics below check which supplied
item does the work in finite models.

**Verdict for D13:** outcome **A (CLOSED-UNITARY ARROW SELECTED) is excluded in scope** by Corollary 1. The question
becomes which of B / C / D / E describes the effective arrows that do occur.
