# SCOUT_0 COUNTEREXAMPLE LEDGER

One decisive counterexample per broad claim attempt. Smallest first.

## CE-01 (P-06): "CM-ness of a kernel is isolated in its deformation space"

- **Claim attacked:** that the S5-1 parent's kernel might be CM-isolated (the charter's null).
- **Smallest counterexample:** the family S_p(x) ∝ (x²+κ²)^{-p} at p=2: kernel f(t)=e^{-κt}(t+1/κ), f''(t)=e^{-κt}(κ²t−κ)<0 on (0,1/κ). Exact, artifact-free, CM fails on an open interval.
- **Nearest counterexample:** p=1.05 — non-CM detected at t=0.08 (order 6), nothing at t≥0.5. Boundary location between 1.0 and 2.0 is bracketed, not proven.
- **Same assumptions, opposite result:** p<1 members are all CM at tested resolution — the opposite side of the same family.
- **Same result, fewer assumptions:** none needed beyond the recorded κ and GR transform table.

## CE-02 (P-06b audit): "if g is CM and α ∈ (0,1], then z^α g(z) is CM" (lemma used in `81ad661`)

- **Claim attacked:** the lemma carrying the `1/2 < p ≤ 1` leg of the banked P-06b proof.
- **Smallest counterexample:** g ≡ 1 (CM, Bernstein measure δ₀); z^α is increasing, not CM.
- **Nearest counterexample:** g = e^{-z}, α = 1: (z e^{-z})′ = e^{-z}(1−z) > 0 on (0,1). The script's integral form of the lemma is also false (RHS/LHS = 1.5 at α=1/2, s=1, z=1; exact RHS = e^{-sz}[s z^{α−1} + (1−α) z^{α−2}]).
- **Same assumptions, opposite result:** none — the lemma is simply false; the specific conclusion it was used for (z^{2ν}·ℒ[(s²−1)^{ν−1/2}] is CM for 0<ν≤1/2) is nevertheless true, because that function equals z^νK_ν(z), which is CM by the unified Bernstein representation.
- **Same result, fewer assumptions:** the unified representation z^νK_ν(z) = √π2^ν/Γ(1−p)·∫_1^∞e^{-zs}(s²−1)^{-p}ds proves CM on all of 0<p<1 with no lemma at all (DLMF 10.32.8 at index 1/2−p). See `probes/P06B_CORRECTION_01.md`.
- **Lesson for the scout:** every numerical check in `81ad661` tested the (true) conclusion, never the lemma. Lemmas introduced to bridge a gap need their own smallest-counterexample pass before banking.

## CE-03 (P-08 audit): "the Wick map is not multiplicative ⇒ the bosonic lift carries no faithful representation of ℝ[x]" and "Sz.-Nagy represents ℝ[x] faithfully by the inclusion H ⊂ K" (`cf0f7fe`)

- **Claims attacked:** the two Reading-2 verdicts of the banked P-08 that fixed its survivor set {cotangent, Sz.-Nagy}.
- **Smallest counterexample (bosonic):** one real site vector under the recorded doubling complex structure: `[Φ(e₁),Φ(e₂)] = i·Im⟨e₁,e₂⟩ = 0`, so `f ↦ f(Φ(e₁),…,Φ(e_N))` (= multiplication by `f` on `L²(ℝᴺ)`) is an injective unital homomorphism. Non-multiplicativity of the Wick map (a linear quantization map) says nothing about the existence of a homomorphism.
- **Smallest counterexample (Sz.-Nagy):** the recorded readout is the amplitude `⟨e₁,U(τ)e₁⟩`, a linear functional, not an operator; and `⟨x⊕0, A(x⊕0)⟩` is even in `x` for every operator `A`, so no operator image reproduces `x₁`. The inclusion `H ⊂ K` maps vectors, not the algebra.
- **Same assumptions, opposite result:** under the map-level reading the bosonic lift *survives*, and under the admissible classical function/Poisson reading inherited from the Λ-H row Sz.-Nagy survives by the pullback `f∘P_H` — so `cf0f7fe` reached a defensible verdict for Sz.-Nagy by an argument that establishes nothing.
- **Same result, fewer assumptions:** the fermionic failure of the conditional ℝ[x] criterion survives with a map- and image-independent proof (finite-dimensional observable algebras; Cayley–Hamilton), replacing the map-specific `x_i ↦ a_i`.
- **Lesson for the scout:** test faithfulness against the *recorded* source-observable map AND the admissible algebra-type reading for the lift (the Λ-H row's "Poisson" entry is a row-level reading, not a per-construction certificate); distinguish "this particular identification is not multiplicative" from "no faithful representation exists"; and apply the criterion uniformly across lifts (granting the cotangent lift `f∘π` while denying Sz.-Nagy `f∘P_H` was non-uniform). See `probes/P08_CORRECTION_01.md`.

## CE-04 (P-08 correction self-audit): "multiplicative and expectation-reproducing are satisfied by different maps and by no single map" (this memo's own first draft)

- **Claim attacked:** the bosonic "dichotomy" claim written in the first draft of `P08_CORRECTION_01.md`, before its adversarial pass.
- **Smallest counterexample:** `π_a(x_i) = √2·a_i`. The `a_i` commute, so `π_a` is an injective unital homomorphism of `ℝ[x]`; coherent states are *joint eigenvectors* of the `a_i`, so `⟨x|π_a(f)|x⟩ = f(x)` for **every** polynomial (checked: `⟨π_a(x)⟩ = 0.9`, `⟨π_a(x²)⟩ = 0.81`, `⟨π_a(x³)⟩ = 0.729` at `x = 0.9`, exactly, with no vacuum variance); and the quasi-free Heisenberg map acts by substitution `f(a) ↦ f(Ta)`, hence multiplicatively and Koopman-intertwined. One map, all three sub-readings.
- **Same assumptions, opposite result:** add the hypothesis that `π` is a `*`-representation (self-adjoint coordinate images) and the dichotomy is restored — then *no* symmetric choice works, by `⟨x|x'⟩ = e^{−|x−x'|²/4} ≠ 0`.
- **Same result, fewer assumptions:** the price is simply recorded instead: `π_a(x_i)* ≠ π_a(x_i)`, the bosonic analogue of the recorded Λ-F odd-image price LS-8. The record's own V-3 intertwiner (`W M_x W⁻¹ = a†`, "not the Segal field") points at this image class.
- **Lesson for the scout:** an impossibility claim must state every hypothesis it uses. The draft's proof silently assumed self-adjointness (it used `⟨π²⟩ = ‖π|x⟩‖²`) while the claim was stated for all maps. Two of fourteen refuters found this independently; it would have been banked otherwise.

## CE-05 (P-17): "the C-B retained-mean discriminator distinguishes ongoing randomness from a deterministic hidden environment"

- **Claim attacked:** a tempting physical reading of S2-1 (not S2-1 as worded, which says "same deterministic state space").
- **Smallest counterexample:** one system coordinate with quartic potential, one harmonic bath mode, linear coupling with counterterm, bath shifted-Gibbs. The retained-mean Taylor coefficients `E q⁽ⁿ⁾(0)` of the autonomous Hamiltonian equal those of the GLE driven by exogenous Gaussian forcing with covariance `Tγ`, exactly, for every `n ≤ 12` — while `S − D ≠ 0` (first at `t⁶`, `−18Taβ`). The discriminator fires identically for a deterministic bath and for primitive noise.
- **Nearest counterexample:** a random-phase (fixed-energy, non-Gaussian) bath: differs from the Gaussian bath at `t¹²` (`138510T²aβ²`) but equals the exogenous random-phase process at every order. What differs is the law, never the ontology.
- **Same assumptions, opposite result:** none inside class 𝓗 (theorem). Outside it (anharmonic bath, coupling nonlinear in bath coordinates) — open; P-18.
- **Same result, fewer assumptions:** the Gaussian/thermal hypothesis is not needed — only that the free-force law is independent of the system preparation (otherwise the deterministic slip `−γ(t)q₀` appears, and is absorbed into the memory term).
- **Lesson for the scout:** "ongoing noise vs initial uncertainty" is a statement relative to a declared state space. Before reading a reduced-data discriminator as ontology, ask whether an enlarged deterministic parent with the same forcing law exists. See `probes/P17_RESULT.md`.

## CE-06 (P-15): "a nonlinear self-consistent drift on branch weights yields Born-like outcome weights"

- **Claim attacked:** the P-15 charter's hope that nonlinear branching dynamics generates the weight law.
- **Smallest counterexample:** `dp = 2p(1−p)(2p−1)dt + √(2p(1−p)) dW` (frozen λ = 2): exact `h(0.25) = 0.21955 ≠ 0.25`, `h(0.75) = 0.78045` (MC-verified within 0.2–0.4 SE).
- **Nearest counterexample:** any `b ≢ 0` — the backward equation forces `h(p) = p ⇒ b ≡ 0`.
- **Same assumptions, opposite result:** `b = 0` with any admissible σ gives `h = p` exactly (martingale convergence), including asymptotic-only absorption (QSD).
- **Same result, fewer assumptions:** decomposition independence alone (outcome frequencies a function of `ρ₁₁` for all mixtures) forces `h` affine, hence Born.
- **Lesson for the scout:** nonlinearity in the branch coordinate is the enemy of Born, not its source; Born fixation is a martingale (linearity) property of a supplied stochastic law. See `probes/P15_RESULT.md`.
