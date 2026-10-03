# V0-2 — COMPARISON PHASE (P-15; the original was unsealed only after the reproduction was committed at `bde523e`)

**Original:** `playground/SCOUT_0/probes/P15_RESULT.md` (blob `0e0bb089…`).

**Reproduction:** `V0_2_P15_REPRODUCTION.md` and `code/v0_2_controls.*`. These were produced by a context-isolated
sub-agent from the spec alone. The reproduction is **not** edited below.

## 1. Scientific comparison

| item | original | independent reproduction | classification |
|---|---|---|---|
| Hypothesis class | σ² bounded away from 0 on compacts; b locally bounded; σ(0) = σ(1) = 0; "absorbing branch diffusion"; h = s/s(1) "whenever s(0+), s(1−) finite and p_t reaches {0, 1} a.s." | H0: σ > 0 and (1 + \|b\|)/σ² ∈ L¹_loc (strictly weaker than the original's class). It requires **A** (absorbing endpoints, stated explicitly, since σ(0) = 0 alone does not guarantee absorption for non-Lipschitz coefficients) and **SF** (finite scale at both ends) for T1 | **equivalent** in substance. The original puts SF and absorption in the weight definition and the class name rather than in the theorem's hypothesis line. The reproduction makes them explicit |
| Backward equation (T1) | b h′ + ½σ²h″ = 0, h(0) = 0, h(1) = 1; scale function | the same, but valid **a.e.** (h ∈ C¹ with absolutely continuous h′; classical C² only where b/σ² is continuous). SF is necessary and sufficient; accessibility is not needed | **agreement plus precision** |
| **Born ⇒ zero drift (T2)** | "h = p ⇒ h″ = 0, h′ = 1, so the backward equation forces **b ≡ 0**" | under the original's own hypotheses (b only locally bounded), h = p gives **b = 0 a.e.** only. Counterexample: b = 1 on the rationals gives the same law and h = p. Pointwise b ≡ 0 holds only for continuous b | **CORRECTION (CR-1)**: the original overstates the conclusion as pointwise |
| Zero drift ⇒ Born (T3) | bounded martingale; a.s. convergence; finite quadratic variation forces ∫σ²(p_t)dt < ∞, so the limit lies in {σ = 0} = {0, 1}; h = p; finite-time or asymptotic | the same mechanism; the "no interior trap" condition is identified as 1/σ² ∈ L¹_loc; finite-time and asymptotic both handled | **equivalent derivation** |
| "Born ⟺ martingale" (T4) | b ≡ 0 ⟺ p is a (bounded) martingale | under H0 + A: b = 0 a.e. ⟺ martingale from every (or one) interior start ⟺ h = p. Pointwise b ≡ 0 is sufficient, not necessary | **follows CR-1** |
| Decomposition independence (T5) | iff h is affine (midpoint-affine plus bounded); with the boundary conditions, iff h = p | real weights: affine with **no** regularity needed; 50 / 50 mixtures: midpoint-affine plus bounded (bounded automatically). **Must include the endpoints (or SF)**: an interior-only counterexample has h ≡ 0 | **agreement plus scope precision (CR-2).** The original's use of h(0) = 0, h(1) = 1 implicitly takes the closed interval, and its weight definition assumes SF, so this states an implicit hypothesis rather than overturning a claim |
| "Function of the density matrix" | the corollary says Born is the unique hitting law with outcome statistics a function of ρ; control F is phrased with "same ρ₁₁" | control F's pure-vs-mixture pair has **different ρ** (off-diagonals ≈ 0.49 vs ≤ 0.445), so it refutes dependence on p̄ (= ρ₁₁), not the ρ-dependence corollary. A **same-ρ** counterexample was supplied: 0.4 / 0.38390 / 0.39271. The corollary itself is true | **precision (CR-3)**: control F tests ρ₁₁-dependence, as its own wording says. The class-result sentence "depend on the state only through ρ₁₁" and the corollary's "function of the density matrix" should be kept distinct |
| Common hypothesis class (T6) | presented as a single biconditional chain | a single class exists: H0 + A + B (b up to null sets), with decompositions including the endpoints. With that class the full chain holds; the original's stated class alone does not give it | **follows CR-1, CR-2** |
| SUPPLIED items | p = \|α\|² supplied; the stochastic law supplied; "never probability derived from determinism" | identical | **exact agreement** |
| C1 | h = p; E[τ] = −(p ln p + (1 − p)ln(1 − p))/D; the A′ variant is asymptotic, with Born holding | identical (symbolic, plus a Green-function quadrature agreeing to ~10⁻⁴¹); A′ endpoints are natural and attracting, with E τ = ∞ | **exact agreement** |
| C2 table | 0.07793, 0.21955, 0.38390, 0.5, 0.61610, 0.78045, 0.92207 | closed form h = [erf(√2(p − ½)) + erf(1/√2)] / [2 erf(1/√2)]; all 7 agree to 5 decimals (and to 10 digits) | **exact agreement** |
| C3 | 0.38390 vs 0.39271 | 0.383900792 vs 0.392713177 | **exact agreement** |

**Out of scope, not compared:** Monte Carlo counts; control C / D / E detail; the Lindblad unravelling; the multibranch
extension.

## 2. Grade

**V0-2-C — P-15 CHAIN REPRODUCED WITH CORRECTION — VER-I1 (orchestrator-exposed).**

**Core chain survives.** Under the explicit common class H0 + A + B, with decompositions including the endpoints:
Born ⟺ b = 0 a.e. ⟺ bounded martingale ⟺ decomposition-independent outcome statistics. All controls C1 – C3 are
reproduced exactly.

**Corrections:**
- **CR-1 (grade-bearing).** "Born ⇒ b ≡ 0" must read **"Born ⇒ b = 0 Lebesgue-a.e."** under the stated hypotheses. It is
  pointwise only for continuous b.
- **CR-2 (scope).** Absorbing endpoints and finite scale (SF) must be **theorem hypotheses**, not only parts of the weight
  definition. The decomposition leg must range over [0, 1], including the endpoints.
- **CR-3 (precision).** Keep dependence on ρ₁₁ distinct from dependence on ρ. Control F refutes ρ₁₁-only dependence. The
  "function of ρ" corollary holds, with a same-ρ counterexample available under control B.

**Independence:** VER-I1, orchestrator-exposed. The orchestrator had seen the original's proof sketch during target
extraction. The reproduction was done by a context-isolated sub-agent from the spec only, with a file-access report.
**Not VER-I2.**

**Frozen `scout-0` is not modified.**

## 3. Owner ruling (additive; review of `55a9aa8`)

**Final grade: V0-2-C — P-15 CHAIN REPRODUCED WITH CORRECTION — VER-I1 (ORCHESTRATOR-EXPOSED).** Criterion-2 item 3:
**COMPLETE AT VER-I1 WITH CORRECTION.**

**Accepted corrected theorem.** The common hypothesis class is:
- an interior diffusion coefficient non-degenerate on (0, 1);
- Engelbert–Schmidt / local-integrability conditions;
- the absorbing-endpoint convention;
- drift taken as an L¹_loc / Lebesgue-a.e. equivalence class;
- decomposition tests over [0, 1], endpoints included.

Under it, **h(p) = p ⟺ b = 0 Lebesgue-a.e. ⟺ p_t is the bounded martingale branch coordinate ⟺ outcome statistics are
decomposition-independent.** Pointwise b ≡ 0 needs stronger regularity, e.g. continuous b.

**Preserved:**
- p = |α|² is SUPPLIED;
- the stochastic law is SUPPLIED;
- probability is not derived;
- absorption supplies selection;
- the martingale structure supplies the Born weights, conditional on those supplied ingredients.

**CR-2 wording repair (supersedes §2 where they differ).** SF is **not** an extra hypothesis of the final biconditional
chain.
- SF (finite scale, i.e. attracting endpoints) is required for the **general** backward boundary-value formula
  𝓛h = 0, h(0) = 0, h(1) = 1, for arbitrary drift.
- **Within the repaired equivalence theorem, SF follows** when the Born / zero-drift / martingale conditions hold.
- Absorbing endpoints remain part of the common class.

**CR-3 accepted.** Dependence only on ρ₁₁ and dependence on the full ρ are kept distinct. The original control-F pair
tests the former. The reproduction's phase-balanced same-ρ decompositions are banked as the correct witness for the ρ
version.

Frozen `scout-0` is not modified.
