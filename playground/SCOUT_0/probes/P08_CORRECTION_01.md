# SCOUT_0 W1 P-08 CORRECTION 01 — the conditional `ℝ[x]` branch, re-audited against the RECORDED lift maps

**Scope:** audit of the banked P-08 (`cf0f7fe`, `probes/P08_RESULT.md`) per the auditor's hold.
**Robust result unchanged:** at the earned reading the faithfulness test is **not formulable as an
operator-representation test at all** (the "full source observable algebra" is not an earned
object; EA-0 `n4`; bridge-ruling output-map warning), so it supplies **no earned selector** and
excludes nothing. **Withdrawn:** the conditional claim "Reading 2 prunes 4 → 2, survivors
{cotangent, Sz.-Nagy}" *as proved* — both of its load-bearing arguments are defective.
**Replaced by:** one theorem-grade, map-independent failure of the conditional criterion (the fermionic lift) plus an explicit
demonstration that **every finer discrimination is reading- and image-dependent**. Terminal
(IRREDUCIBLE/SUPPLIED) untouched; no new selector; faithfulness failure remains a constraint, not
a disqualification (charter fences).

**Self-audit disclosure.** This memo's own first draft was adversarially checked before banking
(14 independent refuters, 3 lenses per claim: pure mathematics, record fidelity, overclaim/grade).
**Three of its six claims were refuted and are corrected here** — see §5. The sub-readings below
are **post hoc relative to `cf0f7fe`** and are *not* pre-registered; they are used only to demote
exclusions, never to create one.

## 0. One-line verdict

> **Earned scope: the test is not formulable as stated; it reduces to D-1 at the declared
> coordinate class and excludes nothing. All four lifts survive.**
> **Conditional `ℝ[x₁,…,x_N]` branch (NEW ASSUMPTION): exactly one failure holds under every
> admissible sub-reading and every admissible image class — the complex fermionic lift, by
> finite-dimensionality of its observable algebra. For the other three, survival is
> reading-dependent: no further pruning is reachable without a further unpriced choice.**

The second sentence is the scientifically load-bearing one: **granted its own new assumption, the
faithfulness constraint still cannot select.** That reinforces the lift terminal rather than
denting it.

## 1. What was wrong in `cf0f7fe`

| Lift | `cf0f7fe` verdict (Reading 2) | Defect |
|---|---|---|
| Λ-B complex bosonic | FAILS — "the observable map is the Wick map; Wick is not multiplicative" | **Category error.** Non-multiplicativity of one *linear quantization map* does not show that no faithful *homomorphism* exists. Under the recorded canonical complex structure (doubling, V-7/LS-9) the coordinate images `Φ(e_i)` **commute** and functional calculus gives an injective unital homomorphism (§2.2). Its premise `[Φ(f),Φ(g)] ≠ 0` is false for the real site vectors. Its "Wick" object was neither normal ordering nor the record's `W = e^{−Δ/2}` but the Gaussian (Isserlis) pairing expansion, hard-coded rather than computed; and it cites LS-5 against LS-5's content (LS-5 records NO MERGE and that *no* algebra isomorphism intertwines `U` and `P`). |
| Λ-F complex fermionic | FAILS — "`x_i ↦ a_i`, `a_i² = 0`" | **Conclusion correct, proof map-specific.** The record's natural coordinate image is the **odd** generator (LS-8: `γ_i`/`B(e_i)`), not necessarily `a_i`. A map-independent proof is given in §2.1. |
| Λ-H Sz.-Nagy | survives — "the source algebra acts on the dilation space by the inclusion `H₀ ⊂ K`: injective and multiplicative" | **Names no map.** The inclusion `H ⊂ K` maps **vectors**; the record represents the coordinates by **amplitudes** `⟨e_i,·⟩` (functionals). A faithful representation does exist under the recorded reading — but it is the **pullback `f ↦ f∘P_H`**, which `cf0f7fe` did not identify (§2.3). |
| Λ-H cotangent | survives — `f ↦ f∘π` | Correct (§2.4). |

Also wrong: `cf0f7fe`'s H4 fallback ("if only injectivity of a linear observable map is meant,
Reading 2 prunes nothing — all four injective on `𝒜_src`"). Λ-F fails even that (§2.1).

## 2. Lift-by-lift, using the recorded maps

**Reading 2** = "the lift restricted to `𝒜_src = ℝ[x₁,…,x_N]` is an injective unital algebra
homomorphism into the lift's observables". Three sub-readings, applied **uniformly** to all four
lifts:

- **(2-map)** such a homomorphism exists, extending the lift's *recorded* coordinate image;
- **(2-state)** (2-map) and the lift's *recorded* state identification reproduces polynomial
  expectations;
- **(2-dyn)** (2-map) and the Heisenberg/transport dynamics on the image intertwines the source
  flow.

Two further choices cut across these, and the record fixes **neither** — this is the heart of the
corrected result:

- **(Image class)** whether the coordinate images must be **self-adjoint** (a `*`-representation).
  The record fixes no coordinate image for the quasi-free lifts: verification §1 lists three
  readouts side by side ("the amplitude; the coherent-state field mean; the Mehler/OU conditional
  mean"), evaluation §1 marks it "declared per V-1" (a price), and V-3 records the canonical
  intertwiner as `W M_x W⁻¹ = a†`, "**not** the Segal field".
- **(Algebra type)** for the Λ-H dilations, whether the lift is read as a **classical function
  (Poisson) algebra** or a **quantum operator algebra**. `L0_LIFT_SELECTION_01.md` §1 places the
  Sz.-Nagy dilation in the Λ-H row whose Algebra column is "**Poisson** (quantizable with priced
  choices)" — the same row as the cotangent lift — and verification §4 says classical and quantum
  dilations "differ in algebra type" without fixing which.

### 2.1 Λ-F complex fermionic — fails under every sub-reading AND every image class (map-independent)

Recorded structure: `Γ_F(e^{−Kt})` over the recorded one-particle space, `N = 23` sites — `ℝ²³`
(real Clifford variant) or `ℂ²³` after the priced doubling (complex Fock variant; V-7, LS-9).
**Every recorded fermionic observable algebra is finite-dimensional** (distinguishing Hilbert-space
from algebra dimension): `Cl(ℝ²³)` and `Λ(ℝ²³)` have real dimension `2²³ = 8 388 608`; the Fock
space `Λ(ℂ²³)` has Hilbert dimension `2²³`, so the CAR algebra `B(Λℂ²³) = M_{2²³}(ℂ)` has complex
dimension `2⁴⁶` (real `2⁴⁷`); the parity-superselected **even** subalgebras — the record's physical
net (V-6/LS-8) — have dimension `2²²` resp. `2⁴⁵`. `ℝ[x₁,…,x₂₃]` is infinite-dimensional over `ℝ`,
so **no injective `ℝ`-linear map** — a fortiori no injective homomorphism — into any of them
exists, whatever the coordinates are sent to (odd `γ_i`, `B(e_i)`, `a_i`, or even bilinears).

Two explicit forms:
- *Dimension count.* Polynomials of degree ≤ 9 already exceed `dim Cl(ℝ²³)`
  (`C(32,9) = 28 048 800 > 2²³`); degree ≤ 27 exceeds `dim_ℂ M_{2²³}(ℂ)`
  (`C(50,23) > 2⁴⁶`); degree ≤ 28 exceeds it real-linearly (script, C1).
- *Cayley–Hamilton (sharper, no bookkeeping).* For **any** image `π(x₁) = A` in `M_{2²³}(ℂ)`,
  `χ_A(A) = 0` with `deg χ_A ≤ 2²³`, so a nonzero **univariate** polynomial already lies in
  `ker π`. For the natural images the annihilating polynomial has degree 2
  (`γ_i² = 1`, `B(e_i)² = 1`, `a_i² = 0`, `n_i² = n_i`).

Secondary (redundant given the above, recorded because it survives to infinite dimension), for the
recorded **odd** images (LS-8): if `π(x_i) = u`, `π(x_j) = v` anticommute, multiplicativity on the
commutative source forces `uv = vu = −uv`, so `uv = 0`, contradicting `γ_iγ_j` invertible (real
Clifford) or `a_ia_j ≠ 0` (complex Fock); and `u² = c` scalar gives `π(x_i² − c) = 0`. Checked on
explicit 2-mode CAR/Majorana matrices (script, C1).

Because (2-state) and (2-dyn) each contain (2-map) as a conjunct, failure of (2-map) propagates to
all three. **This is the one invariant failure in the whole audit** (conditional criterion; not an R-1 exclusion).

*Scope notes.* (i) **Map-independent is not lift-variant-independent:** the dimension count rests
on the recorded *finite* one-particle space. A fermionic-noise dilation of the CP semigroup, or a
continuum/thermodynamic limit, has an infinite-dimensional `CAR(h)` which *does* contain faithful
copies of `ℝ[x]` as an abstract algebra; there only the odd-image obstruction remains, and that one
*is* map-dependent. The recorded Λ-F is the finite-`N` semigroup, so the exclusion is correct as
scoped — carry the hypothesis, never drop it in summary. (ii) **Not transferable to Λ-B:** bosonic
Fock space over `ℂ²³` is infinite-dimensional (§2.2). (iii) Terminology: "fails the constraint" is
a *pruning under the un-earned Reading 2*; faithfulness is not among D-1…D-7 (R-0 closed), so this
is **not** an R-1 exclusion (ruling 02 §2), and the terminal is untouched.

### 2.2 Λ-B complex bosonic — (2-map) survives; the rest is **image-dependent**, not failing

Recorded structure: `ℝ²³` carries no complex structure (V-7); the canonical route is **doubling**
`ℝᴺ → ℂᴺ` (or `T*ℝᴺ`), added directions being `π`-fields conjugate to the site coordinates
(V-7, LS-9). The record writes no literal operator image of `x_i`; the coherent-state field-mean
readout (verification §1) is what *identifies* it as `Φ(e_i)` — an **inference step**, recorded as
such. Recorded quadratic image: `V = ½xᵀPx ↦ dΓ(P)`, "ordering choice" flagged (verification §5,
evaluation §3).

**(2-map): a faithful representation exists, and is unique given the image.**
`[Φ(f),Φ(g)] = i·Im⟨f,g⟩` vanishes for real site vectors (`ℝᴺ ⊂ ℂᴺ` is Lagrangian under doubling).
In the Schrödinger representation `Fock(ℂᴺ) ≅ L²(ℝᴺ)`, `Φ(e_i)` is multiplication by `x_i`, and
`f ↦ f(Φ(e_1),…,Φ(e_N)) = M_f` is an injective (a nonzero polynomial is nonzero a.e.; equivalently
the vacuum law is a nondegenerate Gaussian) unital homomorphism into operators affiliated with the
Weyl algebra, on the common invariant domain of finite-particle/Schwartz vectors. Since `ℝ[x]` is
free, it is the **unique** homomorphism extending `x_i ↦ Φ(e_i)` — so `cf0f7fe`'s criterion and its
Wick premise were inconsistent with each other. Numerically `‖[Φ(e₁),Φ(e₂)]‖ = 0` under doubling,
while under a non-doubling `J` on `ℝ²` (`e₂ = Je₁`, even `N` only) the commutator is `i` and the
images do **not** commute (script, C2).

**(2-state) and (2-dyn): the verdict depends on the image class — `cf0f7fe`-style "FAILS" is not
available.**

- *Self-adjoint (`*`-representation) images.* No `*`-representation reproduces polynomial
  expectations in the recorded coherent states: `π(x_i²) = π(x_i)²` with
  `⟨x|π(x_i)²|x⟩ = x_i² = ⟨x|π(x_i)|x⟩²` would force `|x⟩` to be a joint eigenvector for every `x`,
  hence `⟨x|x'⟩ = 0` for `x ≠ x'`, contradicting `⟨x|x'⟩ = e^{−|x−x'|²/4} ≠ 0` (script, C2b: `0.9139`
  at `x = 0.9, x' = 0.3`). This rules out **every** symmetric generator choice at once — the right
  proof; the earlier "`Φ(e_i)` has continuous spectrum" argument covered only the pinned choice.
  Concretely `⟨x|Φ(e_i)²|x⟩ = x_i² + ½`. The Wick map `:f(Φ):` and the recorded `dΓ(P)` are
  expectation-reproducing but not multiplicative (`dΓ(P) ≠ ½ΣP_ijΦ_iΦ_j` as operators, though
  `⟨x|dΓ(P)|x⟩ = ½xᵀPx`; script, C3). And (2-dyn) fails: the quasi-free Heisenberg map restricted
  to `{f(Φ)}` is the Mehler/OU semigroup, non-multiplicative.
- *Non-self-adjoint images.* `π_a(x_i) = √2·a_i` is an injective unital homomorphism (the `a_i`
  commute; injective because `f(√2a) = 0` forces the symmetric tensor `f` to vanish), coherent
  states are **joint eigenvectors**, so it reproduces **every** polynomial expectation exactly
  (script, C2b: `⟨π_a(x)⟩ = 0.9`, `⟨π_a(x²)⟩ = 0.81`, `⟨π_a(x³)⟩ = 0.729` — no vacuum variance),
  and the quasi-free Heisenberg map acts by substitution `f(a) ↦ f(Ta)`, i.e. multiplicatively and
  **Koopman-intertwined**. So a *single* map satisfies all three sub-readings. Its **price**:
  `π_a(x_i)* ≠ π_a(x_i)`, the exact bosonic analogue of the recorded Λ-F odd-image price LS-8; and
  V-3's own intertwiner `W M_x W⁻¹ = a†` points at this image class, so it is record-grounded, not
  invented. (Residual: an R-3 question — a Koopman-intertwined faithful image sits close to the
  Carleman/Doi–Peliti form that LS-7 calls representational. Not adjudicated here.)

*Normalisation.* In the `Φ = (a+a†)/√2` convention `Var_vac = ½` and the Mehler defect on the
`Φ`-image is `(1−T²)/2`; LS-5's `1−T²` is the unit-variance `L²(γ)` convention. Both appear in
`cf0f7fe`-era text; they cannot be read in one convention (script, C2b).

*Vacuity caveat.* A faithful representation of `ℝ[x]` exists on *every* infinite-dimensional
separable Hilbert space, so map-level survival has content only because it is tied to a recorded
readout image, as above.

### 2.3 Λ-H Sz.-Nagy — **reading-dependent** (not "unsupported": a faithful representation does exist under the recorded algebra type)

Recorded structure: `H` = the one-particle space with the classical state as a **vector** `x ∈ H`;
`K ⊃ H`, unitary `U_t`, `P_H U_t|_H = e^{−Kt}` (needs `K` accretive, R-4); recorded readout = the
**amplitude** `⟨e₁,U(τ)e₁⟩` with zero bath component (verification §1); recorded source-descended
observable = "compression to `H`, `⟨v,Pv⟩`" (verification §5); recorded `(R-b)` image of a quadratic
`V` = `P⊕0`, value `⟨U(t)v,(P⊕0)U(t)v⟩` (evaluation §3). The Algebra column of its row is
**Poisson** (`L0_LIFT_SELECTION_01.md` §1).

- **Under the admissible classical function/Poisson reading inherited from the Λ-H row: (2-map) and (2-state) SURVIVE.**
  The dual of the inclusion `H ⊂ K` is the projection `P_H: K → H`, and
  **`π(f) = f∘P_H`** is an injective unital algebra homomorphism (pullback along a surjection) with
  `δ_{x⊕0}(π(f)) = f(x)` for **all** polynomials — odd ones included (script, C4a: `1`, `x₁`,
  `x₁x₂`, `x₁²x₃ − 2x₂` all reproduced exactly). This is the **same mechanism** the record and this
  memo credit to the cotangent lift (`f∘π`), so denying it to Sz.-Nagy while granting it to the
  cotangent lift would apply the criterion non-uniformly. `cf0f7fe`'s sentence still stands
  withdrawn — it named the inclusion of *vectors*, not this pullback — but the **verdict** it
  reached (survives) is recoverable under the recorded reading, by an argument it did not give.
  **(2-dyn)** fails as an identity on `K`: `P_H U_t` is a compression, not a flow on `K` projecting
  to `φ_t`.
- **Under a quantum operator / vector-state reading: FAILS.** Relative to the recorded operator
  class (compressions `B(H)⊕0`), `dim B(H) = 23² = 529 < dim ℝ[x]_{≤3} = C(26,3) = 2600`: no
  injective linear map (script, C4b). The recorded quadratic image is not multiplicative either —
  commuting source quadratics get non-commuting images (`‖[P_{x₁²}, P_{x₁x₂}]‖ = 1 ≠ 0`) and
  `image(x₁²)·image(x₂²) = 0` while `x₁²x₂² ≠ 0`. And no operator `A` reproduces an odd polynomial
  in the recorded vector states, since `⟨x⊕0, A(x⊕0)⟩` is homogeneous of degree 2 in `x` (degree 0
  normalised), hence even (script, C4c); unnormalised, even `f = 1` fails.
- **No degree-1 operator image is recorded anywhere**, so relative to an operator reading the
  (2-map) status is strictly **NOT-ESTABLISHED** rather than refuted.

*Withdrawn from the first draft:* "a faithful representation can be obtained **only** by composing
with a Koopman or bosonic Fock layer." `K` is infinite-dimensional, so `B(K)` carries many
non-canonical faithful copies of `ℝ[x]`; and `f∘P_H` supplies one directly under the Poisson
reading. The honest statement is "the record-descended constructions that supply one are the
pullback, or a Koopman/Fock layer over `K`".

### 2.4 Λ-H cotangent — survives every sub-reading

Recorded map: `f ↦ f∘π` (verification §5, evaluation §3). Pullback along the surjection
`π: T*ℝᴺ → ℝᴺ` is an injective unital homomorphism and, on `ℝ[x]`, the unique one extending the
recorded coordinate image; the recorded point states `δ_{(x,p)}` reproduce `f(x)`; `H = pᵀf(x)` gives
`ẋ = ∂H/∂p = f(x)`, so the flow projects to the source forward semiflow and
`(f∘π)∘Φ_t = (f∘φ_t)∘π` (script, C5: the `x`-projection is independent of `p₀`). After Weyl/KvN
quantization its degree-0-in-`p` part is `x_i ↦ M_{x_i}` on `L²(ℝᴺ)` (common domain: Schwartz),
still injective and multiplicative — note `ℝ[x]×1` is the associative subalgebra inside LS-7's
`≤1`-in-`p` **Lie** subalgebra.

*Recorded asymmetries, disclosed.* (i) (2-state) survives because the recorded states are classical
**point evaluations** (characters), not Hilbert-space vector states; under a normalisable
vector-state reading `M_{x_i}` has continuous spectrum and the quantized cotangent image would fail
exactly as Λ-B's `*`-image does. (ii) The faithful representation lives entirely in the commutative
pullback/Koopman layer — by R-3 logic this is "survival by containing a representational layer",
the same objection §2.3 weighs for Sz.-Nagy. (iii) (2-dyn) on `ℝ[x]` is strictly ill-posed at
`β > 0` for *every* lift, since `ℝ[x]` is not flow-invariant (the exact Bernoulli solution is not
polynomial); it is well-defined only after extending `π` to `C^∞(ℝᴺ)`, which the pullback does
canonically. (iv) If a *genuinely noncommutative* image of the coordinates were demanded, no lift
could comply (the source algebra is commutative) — that reading is empty, not selective.

## 3. Survivor sets by reading

| Reading / image class | cotangent | Sz.-Nagy | Λ-B complex | Λ-F complex |
|---|---|---|---|---|
| 1 (earned) — no lift excluded | ✓ | ✓ | ✓ | ✓ |
| 2-map, recorded Poisson reading (Λ-H) / any image (Λ-B) | ✓ | ✓ (`f∘P_H`) | ✓ | ✗ |
| 2-map, quantum operator reading | ✓ | ✗ (dim) / NOT-ESTABLISHED (no degree-1 image) | ✓ | ✗ |
| 2-state, recorded Poisson + non-`*` images | ✓ | ✓ | ✓ (`√2 a_i`) | ✗ |
| 2-state, `*`-images / normalisable vector states | ✓ (point states) | ✗ | ✗ | ✗ |
| 2-dyn | ✓ | ✗ (on `K`) | ✓ (`√2 a_i`) / ✗ (`Φ`) | ✗ |

**Exactly one cell pattern is invariant: Λ-F is excluded everywhere; no other lift is.** The
`{cotangent, Sz.-Nagy}` set of `cf0f7fe` *is* recoverable — in the row combining the recorded
Poisson reading for Sz.-Nagy with `*`-images for Λ-B — so the first draft's claim that it "holds in
no column" is **withdrawn**. What remains withdrawn is `cf0f7fe`'s *reasoning* for it, and any
suggestion that it is the unique maximal admissible reading.

Applying the gate's own R-1 logic (an exclusion holding under one admissible reading only is not
identity/theorem grade), and the gate's precedent of listing NOT-ESTABLISHED lifts separately from
the survivor count (evaluation §4, FKM):

> **Λ-F fails the conditional ℝ[x] faithfulness criterion at theorem grade. Established survivors under the recorded readings:
> {cotangent, Sz.-Nagy, Λ-B}; which of the latter two survives a *stronger* reading depends on a
> choice the record does not make.**

The count "3" therefore records *non-exclusion*, not three verified representations; "theorem
grade" attaches to the Λ-F failure alone, never to the count.

## 4. Corrected status

> **P-08: CONDITIONAL-ON-NEW-ASSUMPTION. EARNED SCOPE: the faithful-operator-representation test is
> not formulable (the "full source observable algebra" is not an earned object — the charter's own
> BLOCKED/ambiguity criterion); it reduces to D-1 at the declared coordinate class and excludes
> nothing. CONDITIONAL `ℝ[x]` BRANCH: formulable; exactly one failure of the criterion is theorem-grade and
> image-independent — the complex fermionic lift over the recorded `N = 23` one-particle space, by
> finite-dimensionality. Remaining survivor set READING- AND IMAGE-DEPENDENT / PROVISIONAL.
> Terminal untouched.**

**Earned-scope label.** `cf0f7fe` called Reading 1 "CONSTRAINT-EMPTY". That is a *positive null*
for a test which, at earned scope, is not an operator-representation test at all — the charter's
failure criterion (ambiguity in "the full observable algebra" → BLOCKED) is the one that actually
fires. Recorded here as **NOT-FORMULABLE-AT-EARNED-SCOPE (constraint empty a fortiori)**. Also
corrected: `cf0f7fe`'s "all four survivors are faithful on the earned declared coordinate class"
overstates — Sz.-Nagy's recorded coordinate image is a functional, not an operator, and Λ-F's is
odd with vanishing one-time expectation under superselection (LS-8). What holds is that **no lift
is excluded**, at the readout-identification (D-1/LS-1) level, the identification itself being a
price.

**P-09 inheritance (binding):** run the earned reading on **all four** lifts. Any `ℝ[x]` branch is
run separately, with Λ-F failing the ℝ[x] criterion and Λ-B / Sz.-Nagy flagged provisional, and must declare its
image class (`*` or not) and algebra type (Poisson or operator) **in advance and uniformly**, per
verification §5 and ruling 02 §2. Every added notion of quantum locality is priced as NEW
ASSUMPTION.

## 5. What this memo's own audit changed (adversarial pass, 14 refuters)

| Claim | Verdict | Change |
|---|---|---|
| Λ-F fails map-independently | survived 3/3 | kept; added Cayley–Hamilton form; fixed dimension targets (`2⁴⁶` for CAR, not `2²³`); carried the finite-`N` caveat into the status line |
| Λ-B faithful rep exists (doubling) | survived 3/3 | kept; added uniqueness, domain, vacuity caveat, record-fidelity notes |
| Λ-B "multiplicative and expectation-reproducing by no single map" | **REFUTED ×2** | **false without a `*`-hypothesis**: `π_a(x_i) = √2a_i` satisfies all three sub-readings. Restated as image-dependent; `*`-case reproved via coherent-state non-orthogonality; normalisation inconsistency fixed |
| Sz.-Nagy survival unsupported ⇒ NOT-ESTABLISHED / FAILS | **REFUTED** | **the record's Λ-H algebra type is Poisson**; `f∘P_H` is faithful and reproduces all polynomials at `δ_{x⊕0}`. Restated as reading-dependent; the non-uniform treatment vs the cotangent lift is corrected; the unproved "only" is dropped |
| cotangent survives | survived | kept; disclosed the point-state vs vector-state asymmetry, the R-3 "representational layer" objection, and the `β>0` ill-posedness of (2-dyn) on `ℝ[x]` |
| synthesis 4→3 | survived (grade lens) | count-vs-grade separated; "holds in no column" **withdrawn**; earned label corrected to NOT-FORMULABLE; post-hoc status of sub-readings disclosed |

Two agents (one C6 math refuter, the completeness critic) died on a usage limit; the C6 **overclaim/grade**
lens did run. So the synthesis has grade-level adversarial cover but **no completeness critic** —
recorded as a residual gap. The strongest known open question it would have been asked:
whether any admissible reading of "faithful operator representation of the full source observable
algebra" escapes the three sub-readings above.

## 6. Gate check (theorem gate, for the one theorem-grade item)

- Assumptions: Reading 2 (`𝒜_src = ℝ[x₁,…,x_N]`, NEW ASSUMPTION, not earned); the recorded Λ-F
  structure over the `N = 23` one-particle space.
- Class: the recorded fermionic quasi-free lift (complex Fock or real Clifford variant), including
  its even/superselected subalgebra.
- Conclusion: no faithful — indeed no injective `ℝ`-linear — representation of `𝒜_src` exists.
- Proof: finite-dimensionality (dimension count; Cayley–Hamilton), §2.1; independent of the
  coordinate map and of the image class.
- Assumption removal: an infinite-dimensional fermionic construction (bath dilation, continuum
  limit) escapes the count, leaving only the map-dependent odd-image obstruction.
- Known? Elementary (REDISCOVERED-KNOWN as mathematics); KNOWN-BUT-NEW-IN-GRUT as a constraint
  audit — it is LS-4's Sym²-vs-Λ²/nilpotency invariant in map-independent form.
- Math vs GRUT interpretation: a constraint under a priced assumption. Not a disqualification of
  the fermionic route, not a selector (R-1, charter), terminal untouched.

## 7. Files changed

- `probes/P08_RESULT.md` — correction header; Reading-2 table and verdict rewritten; "scout map"
  section rewritten (P-09 inherits all four at earned scope); correction-script log added.
- `probes/p08b_recorded_map_checks.py` — new; checks C1, C2, C2b, C3, C4a/b/c, C5.
  `p08_faithfulness_checks.py` kept unchanged as the audited artifact.
- `FRONTIER_QUEUE.md` — P-08 and P-09 rows.
- `COUNTEREXAMPLE_LEDGER.md` — CE-03 (the two defective `cf0f7fe` arguments) and CE-04 (this
  memo's own refuted "no single map" claim).

## 8. Wording note 01 (auditor, post-banking; mathematics unchanged)

1. **"Excluded" is reserved for an EARNED R-1 selector** (`L0_LIFT_SELECTION_01.md` §5 R-1). `ℝ[x]` is a
   NEW ASSUMPTION, so the Λ-F result is stated as: **"Λ-F fails the conditional `ℝ[x]` faithfulness
   criterion at theorem grade."** Wherever this memo, `P08_RESULT.md`, the frontier queue or the
   counterexample ledger said "excluded"/"exclusion" for Λ-F, read "fails the conditional criterion".
   The summary-grade lines were reworded in place; §§1–5 keep their audit-trail wording under this note.
2. **Sz.-Nagy algebra type.** The canonical table puts Sz.-Nagy in the Λ-H row whose Algebra column
   says "Poisson", but that row bundles several heterogeneous conservative realizations (cotangent,
   Bateman, FKM baths, Sz.-Nagy). The successful branch is therefore stated as: **"under the
   admissible classical function/Poisson reading inherited from the Λ-H row, `f ↦ f∘P_H` is
   faithful"** — not as a Poisson observable algebra independently established for the Sz.-Nagy
   construction. Phrases above such as "the record's Λ-H algebra type is Poisson" (§2, §5) are to be
   read this way. The conclusion (reading-dependent survival) is unchanged.
