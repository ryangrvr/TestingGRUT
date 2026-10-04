# VER0-B — COMPARISON PHASE (BRI1-X1; the originals were unsealed only after the reproduction was committed at `046a945`)

**Original:** closed branch `grut-backreaction-identifiability-0 @ f2e6999`:
- `BRI1_ANALYTIC_ESCAPE_THEOREM.md` (T1 – T6, §T4-R);
- `BRI1_PF4Q.md` §1 (LEMMA BRI1-R1);
- `BRI1_CANDIDATE_CHARTER.md` §4 (cumulant structure).

**Reproduction:** `VER0_B_BRI1_X1_REPRODUCTION.md` and `code/ver0_b_*`. Produced by a context-isolated sub-agent from
the statement-only spec. The orchestrator re-ran both scripts, and both logs are byte-identical. The reproduction is
**not** edited below.

## 1. Item-by-item comparison

| item | original | independent reproduction | classification |
|---|---|---|---|
| Variational equation | ÿ₁ + (1 + 3x₀²)y₁ = q, zero initial data | identical | **exact agreement** |
| Gibbs domination | Grönwall with ‖A‖ ≤ 2 + 6(√E₀ + c_Q), using x² ≤ 2√H₀; bounds ~P(E₀)e^{O(√E₀)}, integrable since −E₀ + β√E₀ ≤ −E₀/2 + β²/2 | energy-estimate class 𝒟 of dominating functions; notes that a *naive* Grönwall bound exp(cB²t) with B ~ \|z\| would not be integrable | **equivalent derivation.** The original's √E₀ route is valid; the naive-bound warning does not apply to it |
| Parity | (x, p, ε) ↦ (−x, −p, −ε) plus Gibbs symmetry, so κ₃ is odd in ε; E x₀ = 0 | identical (X^{−ε} ~ −X^ε as processes; F_{−q} ~ −F_q); also kills κ₃″(0), κ₂′(0) | **exact agreement** |
| K_q | K = Σ c-terms; one-time K = 3c, with c = Cov(x₀(t)², y₁(t)) | K_q(t) = 3Cov_μ(X⁰(t)², Y_q(t)) | **exact agreement** |
| Small-time coefficient | c(t) = −Var(x₀²)t⁷/(28π³) + O(t⁸); C₀ … C₆ = 0; Var(x₀²) = 1 − m₂ − m₂² > 0 | first non-zero order **t⁷**; c₇[K] = −(3/(28π³))·Var_μ(x²), i.e. **3 × the original's c-coefficient**; orders t⁰ … t⁶ vanish; same Gibbs combination | **exact agreement**, derived independently |
| Finite-order remainder | pathwise Taylor through t⁷ with an integrable C⁸ bound, then E (no interchange needed) | variation-of-constants bound with integrable dominating functions | **equivalent derivation** |
| Cumulant scaling | κ_n(F) = N_B^{1−n/2}κ_n(X^ε); κ₃(F) = K/N_B + O(N_B⁻²) | identical, including parity removal of the N_B^{−1/2} and N_B^{−3/2} terms | **exact agreement** |
| **Threshold logic** | for fixed t\* ∈ (0, δ), N₀(t\*) = ⌊C/\|K\|⌋ + 1 exists and the escape holds for **N_B ≥ N₀** (existential) | **uniform remainder** \|κ₃(X^ε_{P1}(t)) − εK\| ≤ C\*\|ε\|³t¹², for all \|ε\| ≤ 1 and t ≤ π. With δ\* = min(δ, (\|c₇\|/(4C\*))^{1/5}), the escape holds for **every N_B ≥ 1**: **no threshold is needed** | **STRENGTHENING (S-1).** The original statement is true but weaker. The orchestrator checked the t¹² power counting term by term: each ε-derivative of X carries a factor of q = O(s³), so all products in κ₃‴ are O(t¹²) uniformly in \|η\| ≤ 1 |
| Variance positivity | the time-t\* map is a C¹ diffeomorphism and ρ has a density, so x(t\*) is non-degenerate | the flow is a bijection, X^ε(t) is continuous and onto ℝ, and μ has full support | **equivalent derivation** |
| P0 control | κ₃(F_P0) = 0 exactly (symmetric i.i.d. summands) | identical; the Z₂ symmetry of the clamp holds together with q → −q, and P0 is its fixed point | **exact agreement** |
| E₂± orbit escape | \|γ\| is reflection-safe: 0 for P0 vs non-zero for P1, so the family ∉ E₂± (BRI-E2+O) | identical (γ(α + λY) = sgn(λ)γ(Y)); also notes the family ∉ E₁ (immediate from E₁ ⊂ E₂±) and that the conclusion survives even if G = 0 were admitted | **exact agreement** |
| Reservoir limit | Lindeberg–Feller with the R1 bounds; E₁-type for the frozen / pointwise-fixed family; no uniform whole-𝒳 claim (§T4-R) | identical scope; adds that convergence is uniform on {sup\|q\| ≤ Q}, and not uniform over 𝒳, which contains arbitrarily large clamp amplitudes | **agreement**, with a sharper scope statement |
| Earned-claim wording | "…produces an interventional force **law** that cannot be represented…"; the BRI0 synthesis uses "reduced-force laws" | should read "**family**": a single law is trivially in E₂±, and membership is a property of the family over 𝒳 | **precision (P-1)** |
| Non-claims | no primitive randomness, ontology, escape from E_univ, GRUT physics or TRUE COMPRESSION | identical; escape from E_univ fails trivially (take U = the bath initial state) | **exact agreement** |

## 2. Grade

**VER0-B-A — BRI1-X1 THEOREM INDEPENDENTLY REPRODUCED — VER-I1 (orchestrator-exposed / context-isolated reproducer).**

**Every element survives exactly or by an equivalent derivation:**
- the theorem and its assumptions;
- the parity structure;
- K = 3Cov(x₀², y₁);
- the first non-zero order t⁷ with negative coefficient −(3/(28π³))Var(x₀²);
- κ₃(F) = K/N_B + O(N_B⁻²);
- P0 = 0;
- positive variance;
- the E₂± orbit escape (BRI-E2+O);
- the scoped reservoir limit.

**No original claim fails.**

**Recorded alongside the grade (not corrections of error):**
- **S-1 (strengthening, reproducer-proved, orchestrator-checked at proof level).** The finite-N_B threshold is
  unnecessary. There is δ\* > 0 such that for every t ∈ (0, δ\*) and **every N_B ≥ 1**, γ(F_P1(t)) < 0, so the X1 family
  is outside E₂± at **every** finite bath size. The original existential-N₀ statement is implied.
- **P-1 (precision).** "Force law" should read "interventional force **family**" in the earned-claim sentence.

**Owner options:**
- grade VER0-B-C if P-1 is treated as a scope correction;
- decide whether S-1 should be adopted into the BRI record (on a new commit to the closed BRI branch, or only in a
  future paper / synthesis). It is not adopted here.

**Frozen-τ firewall:** unchanged, **PF4Q-I / X1-PF-INDETERMINATE at the frozen τ**. No PF4Q values were used.

**Independence:**
- **VER-I1:** the derivation is independent, but the reproducer is the same model running in a context-isolated
  sub-agent, so this is not external review.
- **Orchestrator-exposed:** I authored BRI1.
- **Not VER-I2.**

## 3. Owner ruling (additive; review of `ae56ed0`)

**VER0-B ACCEPTED.** **Final grade: VER0-B-A — BRI1-X1 THEOREM INDEPENDENTLY REPRODUCED — VER-I1 (orchestrator-exposed /
context-isolated reproducer).**

**P-1:** "interventional force law" → "interventional force **family**" is an **editorial / mathematical precision, not
grade-bearing**. E₂± was already defined at family level, and the proof uses the P0 / P1 family comparison correctly. The
grade remains A.

**S-1 accepted as a VER0 independently reproduced strengthening** (INTERNALLY REPRODUCED / NOT EXTERNALLY REVIEWED):
- There exists δ\* > 0 such that for every 0 < t < δ\* and **every finite N_B ≥ 1**, the X1 interventional force family
  lies outside E₂±.
- **Reason:** K_P1(t) = −c t⁷ + O(t⁸), with c > 0, and |κ₃(X^ε_P1(t)) − εK_P1(t)| ≤ C\*ε³t¹², uniformly for 0 < ε ≤ 1.
  The five-power small-time separation permits one δ\* independent of ε (equivalently, of N_B).

**The closed BRI0 branch is not amended or reopened.** `grut-backreaction-identifiability-0 @ f2e6999` remains the
historical BRI0 theorem boundary. S-1 is banked in VER0 for use in a later standalone paper or synthesis.

**Scientific statement now supported.** A finite reciprocal anharmonic environment can generate an interventional
reduced-force **family** outside the shared causal affine exogenous class E₂± at **every finite bath size**, on a common
non-empty small-time interval. The distinguishing standardised-skewness witness vanishes as O(N_B⁻¹) in the reservoir
limit.

**Firewalls preserved:**
- no primitive randomness;
- no unique ontology;
- no escape from E_univ;
- no GRUT-specific prediction;
- no TRUE COMPRESSION;
- frozen τ remains **PF4Q-I / X1-PF-INDETERMINATE**.
