# V0-2 P-15 — Blind reproduction

**Firewall.**
- The only project file read was `verification_0/specs/V0_2_P15_TARGET.md`.
- No P-15 source, script or log was opened. No other project directory was listed or searched, and git was not used.
- Everything below was derived independently. Textbook results are cited as **KNOWN**.

**Code.**
- Script: `verification_0/code/v0_2_controls.py`. It uses sympy (exact) and mpmath at 40 digits; no Monte Carlo.
- Output: `verification_0/code/v0_2_controls.log`.

---

## 0. Standing notation, hypotheses and known results

The process is dp = b(p)dt + σ(p)dW on [0,1], and ζ := inf{t : p_t ∉ (0,1)} (with ζ = ∞ allowed).

**Outcome event.** The event is O₁ := {lim_{t↑ζ} p_t = 1}, and h(p) := P_p(O₁).
- If ζ < ∞, O₁ is "hit 1".
- If ζ = ∞, O₁ is a *tail* event (asymptotic fixation). It is not realised at any finite time: only P_p(p_t > 1−ε) → h(p) holds.

**Hypotheses** (each is used only where it is named):

- **(H0) Interior non-degeneracy and integrability.** σ > 0 on (0,1), σ(0) = σ(1) = 0, and (1+|b|)/σ² ∈ L¹_loc(0,1).
  - The original class implies (H0). That class is: σ² bounded away from 0 on compacts of (0,1), and b locally bounded.
  - (H0) is strictly weaker.
- **(A) Absorbing endpoints.** The process is stopped at ζ, so h(0) = 0 and h(1) = 1.
  - Writing b(0) = b(1) = 0 is **not** automatically enough for non-Lipschitz coefficients: the SDE can leave a zero of both coefficients, e.g. dp = 3p^{2/3}dt.
  - So (A) is stated as a convention. If b ≡ 0 on [0,1], (A) is automatic, because p and 1−p are nonnegative continuous local martingales (hence supermartingales), and those stay at 0 once they hit it.
- **(SF) Both endpoints attracting.** s(0+) > −∞ and s(1−) < +∞. Here s is the scale function: s′(x) = exp(−∫_{1/2}^x 2b/σ²).
- **(B) Drift as an L¹_loc class.** b is identified up to Lebesgue-null subsets of (0,1).

**Regularity of s.** Under (H0), 2b/σ² ∈ L¹_loc. So s′ is locally absolutely continuous and strictly positive, s ∈ C¹, and s″ = −(2b/σ²)s′ a.e.
- If b is locally bounded and σ² is locally bounded below (the original class), then s′ is locally Lipschitz, i.e. s ∈ W^{2,∞}_loc.
- s ∈ C² **iff** b/σ² has a continuous version; b and σ continuous is sufficient.

**KNOWN results used:**
- **K1 (Engelbert–Schmidt; Karatzas–Shreve Thm 5.5.15).** Under (H0), from every x ∈ (0,1) there is a weak solution up to ζ, and it is unique in law.
  - Pathwise uniqueness is not available in general.
  - h depends only on the law.
- **K2 (Karatzas–Shreve Prop 5.5.22).** Under (H0):
  - (a) If s(0+) = −∞ and s(1−) = +∞, then ζ = ∞ a.s., the process is recurrent, and liminf = 0 and limsup = 1, so there is no limit.
  - (b) If s(0+) > −∞ and s(1−) = +∞, then p_t → 0 a.s. as t↑ζ. The mirror case gives p_t → 1 a.s.
  - (c) Under (SF), p_t → {0,1} a.s. as t↑ζ, and P_x(→1) = (s(x)−s(0))/(s(1)−s(0)).
- **K3 (Feller boundary classification).** Let S = ∫₀ s′, Σ = ∫₀ S(0,x] dM(x) and N = ∫₀ M(0,x] dS(x), where m = 2/(s′σ²).
  - Regular: Σ < ∞ and N < ∞.
  - Exit: Σ < ∞ and N = ∞.
  - Entrance: Σ = ∞ and N < ∞.
  - Natural: both infinite. A natural boundary is *attracting* iff S < ∞.
  - Accessible (reached in finite time) ⇔ Σ < ∞.
- **K4 (occupation-time formula).** ∫₀^t f(p_s) d⟨p⟩_s = ∫ f(a) L_t^a da. Also, the local time L_t^a is a.s. > 0 for every a strictly inside the range of p on [0,t] (regular 1-D diffusion).
- **K5 (Dambis–Dubins–Schwarz).** A continuous local martingale is a time-changed Brownian motion.
- **K6 (Jensen's functional equation; Kuczma, ch. 13).** A solution of f((x+y)/2) = (f(x)+f(y))/2 on an interval is a + A(x), with A additive. If f is bounded on a set of positive measure (in particular, bounded), A is linear.
- **K7.** A bounded continuous local martingale is a martingale, and it converges a.s. and in L².

---

## T1 — Backward equation

**Statement (established with qualification).** Assume (H0), (A) and (SF). Then:

  h(p) = (s(p) − s(0)) / (s(1) − s(0)),   p ∈ [0,1],

and h has the following properties:
- (i) h ∈ C¹[0,1) ∩ C¹(0,1], h′ > 0, and h′ is locally absolutely continuous on (0,1). Under the original class, h ∈ W^{2,∞}_loc.
- (ii) ½σ²h″ + b h′ = 0 **Lebesgue-a.e.** on (0,1). Equivalently, (h′/s′)′ = 0 in the sense of distributions, which is the Feller form (d/dm)(d/ds)h = 0 and holds everywhere.
- (iii) h(0+) = 0 and h(1−) = 1, with h(0) = 0 and h(1) = 1 by (A).
- (iv) h is the **unique** W^{2,1}_loc solution of the boundary-value problem with boundary limits (iii), because every a.e. solution is α + βs.
- The ODE holds **pointwise / classically (h ∈ C²)** only where b/σ² is continuous. With b merely locally bounded, it holds a.e. only.

**Proof.** The formula for h is K2(c). For (ii), differentiate h′ = s′/(s(1)−s(0)) and use s″ = −(2b/σ²)s′ a.e.

**Which boundary assumptions are needed:**
- **(SF) is necessary and sufficient** for the BVP formulation.
  - If (SF) fails, K2(a)/(b) give h ≡ 0 or h ≡ 1 on (0,1).
  - The ODE's bounded solutions are then constants, so no solution has h(0+) = 0 and h(1−) = 1.
- **Accessibility is NOT needed.** Each endpoint may be:
  - exit (finite-time absorption);
  - regular (with absorbing convention (A));
  - natural *attracting* (S < ∞, Σ = ∞: never reached, only approached as t → ∞; see C1-variant).
- Entrance boundaries and non-attracting natural boundaries violate (SF).
- (A) is needed for the endpoint values h(0) = 0 and h(1) = 1. For example, if b(0) > 0 at an accessible 0, the process re-enters and "outcome" is ill-defined.

## T2 — Born ⇒ zero drift

**Statement (established with qualification).** Assume (H0) and (A), and let h(p) = p for every p ∈ (0,1). Then **b = 0 Lebesgue-a.e. on (0,1)**.
- No C² regularity and no a-priori (SF) is used.
- b(0) and b(1) are not determined, beyond what (A) requires.
- If b has a continuous version on (0,1), then that version is ≡ 0 pointwise.

**Proof.**
1. By K2, if (SF) fails then h is constant on (0,1). That contradicts h = p, so (SF) holds.
2. By T1, (s − s(0))/(s(1) − s(0)) = p. So s is affine, s′ is constant, and ∫_{1/2}^x 2b/σ² = 0 for all x.
3. Hence b/σ² = 0 a.e., and since σ > 0, b = 0 a.e. ∎

The equivalent "backward-equation" route also works a.e.: h′ = 1 and h″ = 0 a.e., so b·1 = 0 a.e.

**Why only a.e. (sharp).** Let N ⊂ (0,1) be Lebesgue-null, and put b̃ = b + c·1_N (c ≠ 0 bounded).
- By K4, ∫₀^t 1_N(p_s) ds = ∫ 1_N(a) L_t^a σ^{-2}(a) da = 0.
- So the drift integrals agree, the laws agree (K1 uniqueness), and h̃ = h.
- Hence h determines b only as an L¹_loc class. Example: b = 1 on the rationals, 0 elsewhere, gives h = p while b ≢ 0.

**Accessible states.** Under (H0) every interior point is reached with positive probability from every interior start (regular diffusion). So there is no additional "inaccessible interior set" qualification: the only exceptional sets are Lebesgue-null sets, plus the endpoint values.

**Local and pointwise versions:**
- h = p on a subinterval J implies b = 0 a.e. on J only.
- h(p₀) = p₀ at isolated points implies nothing. C2 has h(½) = ½ with b ≢ 0.

## T3 — Zero drift ⇒ Born

**Statement (established).** Assume σ > 0 on (0,1) with 1/σ² ∈ L¹_loc(0,1) (implied by (H0)), b = 0 a.e. on (0,1), and (A). Then for every p ∈ (0,1):
- p_{t∧ζ} is a [0,1]-valued, hence bounded, martingale;
- p_∞ := lim_{t↑ζ} p_t exists a.s. and lies in {0,1} a.s.;
- h(p) = P_p(p_∞ = 1) = E_p p_∞ = p.

No finiteness of ζ is assumed. Both exit endpoints (ζ < ∞) and natural-attracting endpoints (ζ = ∞) are covered.

**Proof.**
1. b = 0 a.e. and K4 give ∫ b(p_s) ds = 0, so p_{t∧ζ} = p + ∫₀^{t∧ζ} σ dW.
2. This is a continuous local martingale in [0,1], hence a martingale (K7), converging a.s. and in L¹. So E p_∞ = p.
3. Limit in {0,1}. By K5, p_t = β_{⟨p⟩_t} for t < ζ, with β a Brownian motion from p. Let A := ⟨p⟩_ζ ≤ H := inf{u : β_u ∉ (0,1)}.
4. Time change: ζ = ∫₀^A σ^{-2}(β_u) du = ∫ L^a_A(β) σ^{-2}(a) da.
5. Suppose A < H on some event. Then β[0,A] is a compact subset of (0,1), L^a_A(β) is bounded with compact support in (0,1), and 1/σ² ∈ L¹_loc. So ζ < ∞.
6. But ζ < ∞ means p_ζ = β_A ∈ {0,1}, which contradicts A < H. Hence A = H and p_∞ = β_H ∈ {0,1}. ∎

**Simpler variant under the original class.** E⟨p⟩_ζ = E(p_∞ − p)² ≤ 1, so ⟨p⟩_ζ < ∞ a.s. On {ζ = ∞, p_∞ = c ∈ (0,1)}, σ² ≥ ε > 0 near c would make ⟨p⟩_ζ = ∞, a contradiction.

**Condition that forces the limit into {0,1}.** There is no interior trap, i.e. 1/σ² is locally integrable on (0,1). This is sharp:
- If σ(c) = 0 with 1/σ² ∉ L¹ near c, then c is an attracting natural trap, and the martingale converges to c with positive probability.
- Control (H1) in the log: σ = 4√κ p(1−p)|2p−1| and b = 0 give h = 0 on (0, ½] and h = 2p − 1 on [½, 1). This is a **martingale but not Born**.

## T4 — Meaning of "Born ⟺ martingale"

**Strongest justified statement (established with qualification).** Under (H0) and (A), the following are equivalent:
- (a) b = 0 Lebesgue-a.e. on (0,1), i.e. b = 0 as an L¹_loc class;
- (b) for **every** x ∈ (0,1), p_{t∧ζ} is a P_x-martingale;
- (b′) for **some** x ∈ (0,1), p_{t∧ζ} is a P_x-martingale;
- (c) h(x) = x for all x ∈ (0,1).

**Proof:**
- (a ⇒ b) is T3. (b ⇒ b′) is trivial.
- (b′ ⇒ a):
  1. The finite-variation part ∫₀^t b(p_s) ds of a continuous martingale vanishes. So b(p_s) = 0 for a.e. s, and ∫₀^t 1_{b≠0}(p_s) ds = 0.
  2. By K4, L_t^a = 0 for a.e. a ∈ {b ≠ 0}.
  3. From x₀ the process covers any [c,d] ⊂ (0,1) with positive probability, with L^a > 0 on the interior of the range. So {b ≠ 0} ∩ (c,d) is null.
- (a ⇔ c) is T2 + T3.

**Distinctions the spec asked for:**
- **Pointwise b ≡ 0 as a generator statement is NOT equivalent** to (b)/(c). It is sufficient but not necessary, because of Lebesgue-null modifications.
- "Martingale from every interior initial condition" is equivalent to the a.e. statement, and so is "from one interior initial condition".
- Exceptional sets:
  - Lebesgue-null interior sets are invisible.
  - There are no inaccessible interior points under (H0).
  - Endpoint drift values matter only through (A).
- Outside (H0), the equivalence (b) ⇒ (c) fails (interior zero of σ, H1).
- A *measure-valued* drift concentrated on a null set (e.g. skew Brownian motion, a drift ∝ dL^c) is not a function. It is excluded by local boundedness or (H0), and it *does* change h.
- "Bounded" in "bounded martingale" adds nothing, because the state space is [0,1].

## T5 — Decomposition independence

Here h : [0,1] → [0,1] is an arbitrary function (a probability, hence bounded). The ensemble frequency Σ wᵢ h(pᵢ) is the SUPPLIED classical mixing rule.

**(i) Real weights, points in [0,1] (endpoints allowed).** Σ wᵢ h(pᵢ) depends only on p̄ for all finite decompositions **⟺ h is affine on [0,1]. No regularity at all is needed.**
- Proof (⇒): {(1, x)} and {(1−x, 0), (x, 1)} have the same p̄ = x, so h(x) = (1−x)h(0) + x h(1).
- Proof (⇐): linearity.
- With (A): h(0) = 0 and h(1) = 1, so h(p) = p.

**(ii) Only 50/50 decompositions (or rational weights).** The condition gives Jensen's equation.
- **Midpoint-affine + bounded ⇒ affine** (K6). Boundedness is automatic, since h ∈ [0,1].
- Without boundedness or measurability, Hamel-type solutions exist, but they cannot occur here.
- The original's "midpoint-affine plus bounded" is therefore correct, though more than needed when arbitrary real weights are allowed.

**(iii) Decompositions restricted to interior points pᵢ ∈ (0,1).** One gets only "h is affine on (0,1)". **This does NOT yield h = p**, even with h(0) = 0 and h(1) = 1.
- Counterexample (H2, in the log): σ² = 2p(1−p) and b = 3(½ − p).
  - Both endpoints are entrance, s(0+) = −∞ and s(1−) = +∞, and the process is recurrent.
  - So h ≡ 0 on (0,1) and h(1) = 1. h is constant, hence affine, in the interior, with decomposition-independent value 0. Yet h ≠ p and b ≠ 0.
- The repair is one of:
  - include endpoints in the decomposition class, as in (i); or
  - assume (SF), giving h(0+) = 0 and h(1−) = 1. Under (H0) + (A), K2 shows h on (0,1) is either constant 0, constant 1, or the normalised scale. So "affine and non-constant on (0,1)" already forces (SF) and h = p.

**Continuity.** No continuity of h is required in (i)–(ii). Under (H0) + (A), h is in fact continuous on (0,1), and can be discontinuous only at an endpoint (H2, or one-sided K2(b)).

**Density-matrix version** (the original corollary). The condition "frequency is a function of ρ = Σ wᵢ|ψᵢ⟩⟨ψᵢ|" is *weaker* than "function of p̄ = ρ₁₁". It still forces affinity, with no regularity:
- ρ = diag(x, 1−x) equals both {(x, |1⟩), (1−x, |0⟩)} and {½, √x|1⟩ ± √(1−x)|0⟩}.
- So h(x) = x h(1) + (1−x) h(0).

The corollary "Born is the unique hitting law making statistics a function of ρ" is therefore **true** under (A), given the supplied facts:
- the pure-state outcome law depends only on p;
- the classical mixing rule.

**But C3 as posed is not a same-ρ comparison** (see C3).

## T6 — Complete chain

**Verdict: one common hypothesis class does work, after scope repair (established with qualification).**

**Hypothesis class H\*:**
- (H0) — implied by the original class;
- (A) — absorbing endpoints;
- (B) — drift as an L¹_loc class;
- the decomposition clause ranges over all finite decompositions with points in **[0,1]**, endpoints allowed. Alternatively, keep interior points and add (SF).

**Under H\*, the following are equivalent:**
1. **Born:** h(p) = p for all p ∈ (0,1);
2. **zero drift:** b = 0 Lebesgue-a.e. on (0,1);
3. **martingale:** p_{t∧ζ} is a (necessarily bounded, a.s. {0,1}-convergent) P_x-martingale for every x ∈ (0,1), equivalently for one x ∈ (0,1);
4. **decomposition independence:** Σ wᵢ h(pᵢ) depends only on p̄ for every finite decomposition. Equivalently (4′): it depends only on ρ.

**Proofs:**
- 1⇔2⇔3 is T2–T4.
- 1⇒4 and 1⇒4′ follow by linearity.
- 4⇒1 is T5(i) with (A). 4′⇒1 is the ρ argument.

(SF) is a *consequence* of each item, not a hypothesis. It is needed only to state T1's BVP for general b.

**What the original class does and does not deliver.** "σ² bounded away from 0 on compacts, b locally bounded" implies (H0), but on its own it is **insufficient / overstated**:
- (a) "zero drift" holds only a.e. (T2). Pointwise b ≡ 0 requires b continuous.
- (b) Endpoint absorption (A) is not stated; b at the endpoints is unconstrained.
- (c) With interior-only decompositions, the decomposition leg fails (H2).
- (d) The "via backward equation" route needs (SF) and gives h ∈ W^{2,∞}_loc, not C², with the ODE a.e.
- (e) The original's own non-degeneracy hypothesis is what makes Born ⇐ martingale work. Dropping it (an interior zero of σ) breaks 3⇒1.

If one insists on *pointwise* "b = 0", the correct statement is the conditional set:
- (b ≡ 0) ⇒ 1, 3, 4;
- 1 ⇒ (b = 0 a.e.), and pointwise if b is continuous.

---

## Controls (see `code/v0_2_controls.log`)

### C1 — Driftless Wright–Fisher: σ² = 2Dp(1−p), b = 0

**Hitting law.** s(x) = x, so h = p. This holds by T3, and the exact ODE residual is 0.

**Boundary type.** Feller numerics at ε → 0 give S → ¼, Σ → 0.287682 (finite) and M, N ~ ln(1/ε) → ∞.
- So 0 is an **exit** boundary, and by symmetry so is 1. Absorption happens in finite time.

**Expected absorption time.** Solve ½σ²T″ = −1 with T(0) = T(1) = 0:
- T″ = −(1/D)(1/p + 1/(1−p)), which gives **E_p τ = −[p ln p + (1−p) ln(1−p)]/D**.
- The symbolic residual is 0, and T(0+) = T(1−) = 0.
- Justification:
  1. T(p_{t∧τ}) + t∧τ is a martingale (Itô; T ∈ C² inside and bounded on [0,1]).
  2. So E[t∧τ] ≤ T(p). Monotone convergence then gives E τ ≤ T(p) < ∞.
  3. Bounded convergence gives E τ = T(p).
- Independent Green-function quadrature, E_p τ = ∫ min(p,y)(1−max(p,y)) · 2/σ²(y) dy with D = 1, agrees to about 1e-41:

| p | E_p τ |
|---|---|
| 0.1, 0.9 | 0.32508297339144823951 |
| 0.25 | 0.56233514461880835029 |
| 0.5 | 0.69314718055994530942 = ln 2 |

- **Agrees with the spec.**

**Variant σ = 4√κ p(1−p), b = 0.**
- Feller: s(x) = x, so S = ¼ < ∞ and the boundary is attracting.
- M, Σ and N are all infinite: Σ grows like (1/(8κ)) ln(1/ε).
- So both endpoints are **natural, attracting**:
  - never hit, so ζ = ∞ a.s. and E τ = ∞;
  - the equation ½σ²T″ = −1 has no bounded solution (T → −∞ at 0+).
- h = p holds by T3 (the time-change argument covers ζ = ∞).
- **Independent check (logit).** Let y = ln(p/(1−p)). Then dy = 8κ tanh(y/2) dt + 4√κ dW.
  - s_y′ = sech²(y/2), so s_y = 2 tanh(y/2) = 2(2p−1), with finite limits ±2.
  - P(y → +∞) = e^y/(1+e^y) = p, exactly (sympy).
- The spec claims ("natural boundaries, asymptotic-only absorption", h = p) are **confirmed**. Precision: the boundaries are natural *attracting*, and "outcome 1" here is the tail event p_t → 1.

### C2 — b = λp(1−p)(2p−1), σ² = 2Dp(1−p)

**Scale function.** 2b/σ² = (λ/D)(2p−1), so s′(x) = exp(−(λ/D)(x−½)²). This gives

  h(p) = [erf(√(λ/D)(p−½)) + erf(½√(λ/D))] / [2 erf(½√(λ/D))].

For λ = 2 and D = 1: h(p) = [erf(√2(p−½)) + erf(1/√2)] / [2 erf(1/√2)].
- The exact ODE residual is 0.
- Both endpoints are **exit** (S → 0.187850, Σ → 0.260918, N → ∞), and b = σ = 0 at the endpoints. So T1's hypotheses hold, with finite-time absorption.
- h ∈ C^∞, since the coefficients are smooth.

| p | h(p) (erf, 20 d; quadrature agrees to 1e-41) | spec | h − p |
|---|---|---|---|
| 0.1 | 0.077927293835155296592 | 0.07793 ✓ | −0.022072706 |
| 0.25 | 0.21954678740599844567 | 0.21955 ✓ | −0.030453213 |
| 0.4 | 0.38390079186514348908 | 0.38390 ✓ | −0.016099208 |
| 0.5 | 0.5 | 0.50000 ✓ | 0 |
| 0.6 | 0.61609920813485651092 | 0.61610 ✓ | +0.016099208 |
| 0.75 | 0.78045321259400155433 | 0.78045 ✓ | +0.030453213 |
| 0.9 | 0.92207270616484470341 | 0.92207 ✓ | +0.022072706 |

- **All 7 spec values agree to 5 dp.**
- h(p) + h(1−p) = 1 exactly, so h = p only at p ∈ {0, ½, 1}.

### C3 — Pure versus mixture under C2

| case | outcome-1 frequency | spec |
|---|---|---|
| pure p = 0.4 | 0.383900791865143 | 0.38390 ✓ |
| 50/50 mixture of p = 0.2 (h = 0.169327146748297) and p = 0.6 (h = 0.616099208134857) | **0.392713177441577** | 0.39271 ✓ |
| difference | 0.008812385576 | |

- Under C1 both frequencies are 0.4 exactly. **Agrees with the spec.**

**Correction / precision.** The pure state with p = 0.4 and the 50/50 mixture of pure p = 0.2 and p = 0.6 states have **different density matrices**:
- the pure state has |ρ₀₁| = √0.24 ≈ 0.49;
- the mixture has |ρ₀₁| ≤ 0.445.

So C3 refutes p̄-dependence, not ρ-dependence. A genuine same-ρ refutation, computed in the log for ρ = diag(0.4, 0.6) under C2:

| decomposition of ρ = diag(0.4, 0.6) | outcome-1 frequency |
|---|---|
| {\|1⟩, \|0⟩}, weights 0.4 / 0.6 | 0.4 |
| p = 0.4, phases 0 and π, weights ½ / ½ | 0.3839007919 |
| p = 0.2 and p = 0.6, phases 0 and π each, weights ¼ | 0.3927131774 |

So the C3 numbers do certify ρ-dependence, but only once phases are chosen to equalise ρ.

---

## Hostile assumption audit (spec §4)

1. **σ has an interior zero.**
   - If 1/σ² ∉ L¹ near the zero c, then c is an attracting trap. With b = 0 the process is still a bounded martingale, but h ≠ p. In H1, h(0.4) = 0 and h(p) = 2p − 1 for p > ½.
   - If 1/σ² ∈ L¹ near c, Engelbert–Schmidt uniqueness fails (sticky solutions exist), so h is not well defined.
   - Either way the chain breaks at martingale ⇒ Born. Non-degeneracy (H0) is essential.
2. **A boundary is inaccessible.**
   - **Inaccessible-but-attracting** (natural, S < ∞) is harmless: T1–T6 all hold, with "outcome" as an asymptotic tail event (C1 variant).
   - **Non-attracting** (entrance or natural with S = ∞) breaks T1. Then h ≡ 0 or 1 on (0,1), Born fails, and b ≢ 0, because b = 0 forces s = x and hence attraction.
   - The equivalences survive (both sides are false). The interior-only decomposition leg fails, though (H2).
3. **Absorption/fixation is not almost sure.**
   - Under (H0) this occurs only in K2(a), the recurrent case where both scale ends are infinite. Then h ≡ 0 on (0,1).
   - Under (H0) with b = 0, fixation is automatically a.s. (T3).
   - Outside (H0), interior traps give convergence to an interior point (item 1).
   - So "a.s. fixation" need not be assumed for the chain; it is a consequence of b = 0 + (H0). It *is* needed (in the form SF) for T1's BVP.
4. **b ≠ 0 only on a polar/inaccessible set.**
   - In (0,1) no point is polar under (H0). But any Lebesgue-null set is *occupation-null* (K4), so b may be ≠ 0 there with an identical law, h = p, and the martingale property intact.
   - Hence only b = 0 a.e. is characterised.
   - Endpoint values of b matter only through (A).
   - A *measure* drift on a null set (skew BM) is visible, but it is excluded by local boundedness.
5. **Is boundedness alone enough for the decomposition argument?**
   - With arbitrary real weights, nothing is needed: not even boundedness.
   - With only midpoint (50/50) decompositions, midpoint-affinity + boundedness suffices (K6), and boundedness is automatic since h ∈ [0,1].
   - Boundedness does **not** rescue interior-only decompositions: endpoints or (SF) are required (H2 is bounded and interior-affine, but not Born).
6. **Does h = p characterise zero drift pointwise under the stated regularity?**
   - **No.** It gives b = 0 Lebesgue-a.e. on (0,1) only, since b locally bounded is not continuous.
   - Pointwise holds iff b is taken continuous, or as its canonical L¹_loc class.
   - b(0) and b(1) are undetermined.
7. **Is p = |α|² derived?** No. It is used only as the name of the coordinate, and it remains SUPPLIED. The ρ-version of T5 additionally uses the supplied assumptions that the dynamics, and hence h, depend on the pure state only through p.
8. **Does the chain derive probability?** No.
   - The law P_p (Wiener measure driving dW), the identification of h with an outcome frequency, and the mixing rule Σ wᵢ h(pᵢ) are all SUPPLIED.
   - The chain shows only consistency: *given* that law, frequencies equal p iff the drift vanishes a.e.

## Corrections to the spec / original

1. **T2 is overstated** if read pointwise. The correct conclusion is b = 0 a.e.; pointwise only for continuous b.
2. **The backward equation holds a.e., not classically.** With b only locally bounded, h ∈ W^{2,∞}_loc, not C². The BVP needs (SF), and accessibility is not required.
3. **The decomposition leg needs endpoints (or SF).** Midpoint-affine + bounded is correct, but unnecessary with real weights.
4. **The hypothesis class must add** absorbing endpoints (A) and the a.e. reading of "zero drift". (H0) is a strictly weaker sufficient replacement for "σ² bounded below on compacts + b locally bounded".
5. **C3 compares states with different ρ**, so as stated it refutes p̄-dependence, not "function of the density matrix". Same-ρ numbers are given above. The density-matrix corollary itself is true under (A).
6. **All numerical values (C1, C2 table, C3) are confirmed.** The C1 variant has natural *attracting* boundaries with E τ = ∞.
