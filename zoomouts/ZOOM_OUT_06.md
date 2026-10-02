# SCOUT-2 ZOOM-OUT 06 (after REPAIR 03 and S2-D-arrow)

## The owner's six questions

**1. Is H_env really the residual item, or should it be replaced by H_corr / freshness?**
Replace it.
- A **maximally mixed** fresh bath gives a Markovian arrow, an entropy rise (0.633 → 0.997) and records (0.225).
- A **6.9-bit** (near-maximal-entropy) bath **with S–E correlations** gives an exact **anti-arrow** (S_S 0.955 → 0).
- Bath-internal correlations with maximally mixed marginals break Markovianity.
- A reused (non-fresh) bath kills the semigroup (M = 1 repurifies the system).

The residual boundary datum is

    H_arrow = H_corr = independence / low inter-subsystem correlation (incl. freshness of
              not-yet-interacted degrees of freedom), relative to Σ, at a chosen time.

Low entropy is **neither necessary** (maximally mixed bath) **nor sufficient** (correlated low-entropy-marginal
counterexamples are possible by the same construction).

**2. Is the arrow definable independently of Σ?** **No.**
- The same global state is product in one frame (the Householder frame, S_S(0) = 0, then rising) and at equilibrium in
  another (Clifford frame, flat at ≈ 1.9 – 2.0).
- "Low correlation" is a Σ-relative property.
- The exact mirror theorem (D0 Corollary 2) itself needs a Σ-local time reversal.

→ **The arrow price is jointly H + Σ.**

**3. Did any unitary mechanism work for every state?** **No, and it cannot.** D0 Theorem 1 / Corollary 1: in finite
closed unitary dynamics no continuous arrow functional is monotone along any orbit forever, nor for every state over
any window. Outcome A is excluded in scope.

**4. Did typicality contribute anything beyond equilibrium?** **No.** Typical states (Haar, or Haar on an energy shell)
are already at equilibrium: ⟨S_S⟩ = 1.98 – 1.99 of 2, and P(increase) = 0.499 ± 0.013. Typicality gives the
equilibrium *appearance*, relative to a supplied measure and constraint. It gives no direction.

**5. Is A still required for every irreversibility claim?** **Yes.**
- Global distinguishability is constant (0.0158 at every t; global TD 0.992… in H2).
- Every arrow observed lives in reduced, coarse or fragment observables (A_res).
- A_res is necessary, not sufficient: coarse entropies are themselves covered by Theorem 1, and they mirror exactly.

**6. What genuinely irreducible specification remains?**

    D_dyn          — dynamics (+, for compressions, a non-injective reduced law obtained by tracing)
    H_corr|Σ       — the independence / freshness boundary condition relative to Σ, at an epoch
                     (replaces H_env; H_measure → D and H_state → D mod gauge were shown in S2-3b / H2
                      as scoped compressions)
    Σ              — (d, #factors, preference order, scale) / Sym(H)
    A_res          — reduced / coarse / fragment access (needed for every observed irreversibility)
    orientation    — NOT selected: arrows are two-sided about the special condition (Janus, D10)

**New coupling found:** the boundary condition and Σ are **one object**. A special state can be fixed structurally by
(H, Σ) + a selector functional, e.g. the mean-field product state, unique modulo the reflection symmetry. It can
**never** be fixed by H alone (Proposition 2: H-only covariant states are stationary). And different selector
functionals pick different states (overlap 0.000), the same preference pricing as S2-Σ.

## The ten questions (condensed)

| # | question | answer |
|---|---|---|
| 1 | derived vs renamed | D0 is a **TRUE DERIVATION** (scoped no-go). H_env → H_corr\|Σ is a sharpening, not an elimination |
| 2 | measure | typicality (Haar / shell): equilibrium only |
| 3 | access | needed for every arrow (A_res) |
| 4 | unique decomposition | unchanged; and now **coupled to the boundary condition** |
| 5 | composition | unchanged |
| 6 | basin | H_corr\|Σ at an epoch |
| 7 | observer | records need H_corr (they un-form from the reversed state, 6 → 0) |
| 8 | information decreased? | **yes, in a scoped sense**: "low-entropy universe" is replaced by the more specific "independence relative to Σ", and the middle condition can be compressed into (H, Σ) + selector |
| 9 | boundary moved | yes: the C5-H remainder is now H_corr\|Σ, plus an unselected orientation |
| 10 | hostile that would overturn this | (i) an *infinite* or continuous-spectrum model where an arrow emerges from every state of a physically motivated class without an independence condition (outside D0's scope; e.g. scattering); or (ii) a selector functional on (H, Σ) singled out by a principle with no preference choice |

## Candidate T2-5 (classification + scoped theorem; **CANONICAL-CANDIDATE: NO**)

> **(Theorem part, finite closed unitary scope.)** No every-state arrow exists for continuous functionals (D0). States
> fixed covariantly by H alone are stationary (Proposition 2).
>
> **(Classification part, tested models.)** Every effective subsystem or record arrow found requires an
> independence / freshness boundary condition relative to a factorization Σ (H_corr|Σ), together with reduced or
> coarse access (A_res). It does not require low environment entropy. The special condition can be fixed by (H, Σ)
> plus a selector functional, never by H alone. Arrows are two-sided about it, so the time orientation is not
> selected.

T2-4 (OPEN CANDIDATE) is **resolved in favour of the owner's suspicion**: the price is correlation / freshness relative
to Σ, not "low-entropy environment".

## Next-phase decision

| option | assessment |
|---|---|
| **S2-G2 / G3 / G4** (dimension notions) | **recommended next.** Σ's local-dimension component (G1) is a primitive. G2 / G3 / G4 test whether larger-scale dimension notions are equally supplied or emerge from D + H_corr\|Σ |
| S2-Σb (epoch-integrated criteria) | still held. The epoch now appears **inside H_corr\|Σ**; Σb could be merged into a later "boundary selector" probe |
| new D-hostile: infinite / continuous-spectrum arrow (outside D0's scope) | valuable later; it is the main loophole in T2-5 |
| saturation | **not reached** |
