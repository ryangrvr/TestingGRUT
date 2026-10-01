# SCOUT_0 W3 P-15 RESULT — branch selection, branch probabilities and Born weights (one analytic pass)

**Charter:** `PROBE_CHARTERS.md` P-15, sharpened by the auditor into a theorem hunt on the absorbing
branch diffusion. **Three questions kept separate throughout** (as Experiment P required):

1. decoherence;
2. a definite single outcome (selection = absorption);
3. the weight law over outcomes (hitting probability `h`).

**Script:** `p15_branch_hitting.py` (log `p15_branch_hitting.log`); leak check `p15_leak_check.py`.
**Constants frozen in the script header before any evaluation:** `D = 1`, control B `λ = 2`, control D
`γ = 4, p₀ = 1/2`, control E `κ = 1/4`. Monte Carlo is used **only** to verify exact hitting
probabilities. NEW HYPOTHESIS CLASS — `EXPERIMENT_P_FRONTIER = CLOSED_AT_TESTED_MODEL_CLASS` is
untouched. P-15 deliberately changes the premise (nonlinear/stochastic, non-unitary selection).

## 0. Verdict (layered labels)

| Layer | Label |
|---|---|
| S2 C-B variable → branch population | **E. UNFORMULABLE-FROM-S2** — the toy is an analogy only |
| driftless absorbing diffusion (A) | **C. BORN-CONDITIONAL-MARTINGALE** |
| frozen nonlinear drift (B) | **B. SELECTION-NON-BORN** — a NEGATIVE result for that model |
| noise off (C) | **D. DETERMINISTIC-BASINS / MEASURE-REQUIRED** |
| noise without absorption (D); dephasing / random-unitary (E) | **A. NO-SELECTION** |

> **Class result (REDISCOVERED-KNOWN mathematics, banked as such):** for the absorbing branch
> diffusion `dp = b dt + σ dW`, **Born hitting probabilities `h(p) = p` for all `p` ⟺ `b ≡ 0` ⟺ `p` is
> a (bounded) martingale.** Equivalently, Born ⟺ the outcome frequencies depend on the state only
> through `ρ₁₁` (decomposition independence).
>
> - Selection comes from **absorption**.
> - Trajectory probabilities come from the **supplied stochastic law**.
> - Born weights come from the **martingale property**.
> - `p = |α|²` is **supplied** quantum structure.
>
> A nonlinear drift in the branch coordinate generically **destroys** Born weights. It does not
> produce them.

## 1. Pre-flight: does S2 define a branch weight?

**No.** The C-B potential `V = ½xᵀK_bx + βΣx⁴` is strictly convex (`min eig K_b = 0.3045 > 0`,
`β ≥ 0`). It has one equilibrium, no two attractors, no bounded `[0, 1]` coordinate, and no canonical
map to `|α|²`. **UNFORMULABLE-FROM-S2.** Everything below is a NEW TOY HYPOTHESIS CLASS, priced:

- `p ∈ [0, 1]` is a **branch coordinate**, with `p = 0, 1` the two definite branches;
- reading `p = |α|²` is SUPPLIED;
- `p` is not called a probability until the hitting theorem has been applied.

## 2. Exact core and theorem

**Generator:** `𝓛 = b(p)∂_p + ½σ²(p)∂_p²` on `(0, 1)`, with `σ > 0` inside and `σ(0) = σ(1) = 0`.
**Weight:** `h(p) = P_p[p_t → 1]` solves `b h′ + ½σ² h″ = 0`, `h(0) = 0`, `h(1) = 1`. With scale
function `s′(p) = exp(−∫2b/σ²)`, `h = s/s(1)` whenever `s(0+)`, `s(1−)` are finite and `p_t` reaches
`{0, 1}` a.s.

**Theorem (two-branch Born criterion).** Assume `σ²` is bounded away from 0 on compact subsets of
`(0, 1)` and `b` is locally bounded.

- **(⇐)** If `b ≡ 0`, `p_t` is a bounded martingale. It converges a.s. Finite quadratic variation
  forces `∫σ²(p_t)dt < ∞`, so the limit lies in `{σ = 0} = {0, 1}`. Bounded convergence gives
  `h(p) = E_p[p_∞] = p`. This holds whether absorption is in finite time or only asymptotic.
- **(⇒)** If `h(p) = p` on `(0, 1)`, then `h″ = 0` and `h′ = 1`, so the backward equation forces
  `b ≡ 0`.
- **Corollary (decomposition independence).** A mixture of branch coordinates `{(w_i, p_i)}` with
  the same `ρ₁₁ = Σw_i p_i` gives outcome frequency `Σw_i h(p_i)`. That equals `h(ρ₁₁)` for every
  mixture iff `h` is affine (midpoint-affine plus bounded). With the boundary conditions, iff
  `h(p) = p`. **Born is the unique hitting law under which outcome statistics are a function of the
  density matrix.** This is the known no-signalling/linearity link (Pearle 1976; Gisin 1984/1989;
  Bassi–Ghirardi, Phys. Rep. 379, 2003).

**Multibranch (noted, not run):** on the simplex with absorbing vertices, Born ⟺ every `p_i` is a
bounded martingale (all drifts zero); multi-allele Wright–Fisher and diffusive QSD are standard
examples.

## 3. Controls

- **A. Driftless, `σ² = 2Dp(1−p)`:** `s′ = 1`, so `h(p) = p` exactly.
  - Absorption occurs in **finite** time: `E[τ] = −(p ln p + (1−p) ln(1−p))/D` (exit boundaries;
    verified symbolically).
  - **A′ (`σ = 4√κ p(1−p)`, the diffusive QSD unraveling):** natural boundaries, so absorption is
    only **asymptotic**, and still `h = p` by martingale convergence.
  - **Label:** BORN-CONDITIONAL-MARTINGALE. Probabilities enter through `dW`, and `|α|²` enters
    through the supplied identification.
- **B. Frozen `b = λp(1−p)(2p−1)`, `λ = 2`, same σ:** exact `s′ = e^{c·y(1−y)}`, `c = λ/D = 2`.

  | `p` | 0.1 | 0.25 | 0.4 | 0.5 | 0.6 | 0.75 | 0.9 |
  |---|---|---|---|---|---|---|---|
  | `h(p)` | 0.07793 | 0.21955 | 0.38390 | 0.50000 | 0.61610 | 0.78045 | 0.92207 |
  | `h − p` | −0.022 | −0.030 | −0.016 | 0 | +0.016 | +0.030 | +0.022 |

  The majority branch is amplified. Monte Carlo (20 000 paths, all absorbed) agrees within
  0.2–0.4 standard errors (`p₀ = 0.25`: 0.2201 vs 0.2195; `p₀ = 0.6`: 0.6149 vs 0.6161). Control A
  agrees within 0.8 SE. **SELECTION-NON-BORN.** Not retuned.
- **C. σ = 0:** fixed points at 0, ½ (unstable) and 1. Initial conditions with `p₀ > ½` go to 1, the
  rest to 0, exponentially and never in finite time. A "branch-1 frequency" `μ({p₀ > ½})` needs a
  supplied measure `μ` on initial conditions; basin volume is not a probability unless that measure
  is given. **DETERMINISTIC-BASINS / MEASURE-REQUIRED.**
- **D. `b = γ(p₀ − p)`, σ as in A, `γ = 4, p₀ = ½`:** stationary Beta(2, 2), both boundaries
  unattainable (Feller exponents ≥ 1). The process is stochastic and fluctuating but never selects.
  The Monte Carlo "absorbed" fraction is Euler leakage, shrinking with step size: 0.0070 → 0.0040 →
  0.0000 at dt = 2e-4, 5e-5, 1.25e-5. **NO-SELECTION.**
- **E. One dephasing Lindblad (`L = √κ σ_z`), two unravelings:**
  - Ensemble coherence at `t = 4`: random-unitary 0.0715, QSD 0.0635, exact `√0.21·e^{−2κt}` =
    0.0620 (random-unitary within ~2σ of its MC noise).
  - Populations: the random-unitary unraveling keeps `p_t = 0.3` on every trajectory. The QSD
    unraveling (`dp = 4√κ p(1−p)dW`, driftless, derived here from the QSD equation) has mean
    0.3013 with **90.7 %** of trajectories within 0.05 of a branch.
  - **Identical decoherence, selection only in one unraveling.** The ensemble master equation does
    not determine selection.
- **F. Decomposition test:** under B, pure `p = 0.4` gives `h = 0.38390`, while a 50/50 mixture of
  `p = 0.2, 0.6` (same `ρ₁₁`) gives 0.39271. The frequencies depend on the decomposition. Under A
  both give 0.4.

## 4. Measure relocation (charter wording interrogated)

"Unequal weights without postulating a measure" splits in two:

1. **no separate postulate of outcome weights** — achievable: the hitting problem turns the initial
   branch coordinate into outcome frequencies;
2. **no probability measure anywhere** — not achieved: `dW` is a measure on histories (and in C a
   measure on initial conditions is required).

So the strongest positive statement is **BORN WEIGHTS GENERATED AS HITTING PROBABILITIES CONDITIONAL
ON A SUPPLIED STOCHASTIC LAW (MARTINGALE BRANCH COORDINATE)**, never "probability derived from
determinism".

**Circularity audit:**
- Noise is not weighted by `|α|²`.
- Initial states are not Born-sampled (each run starts at a fixed `p₀`).
- Drift and σ were frozen before evaluation.
- The absorbing boundaries are the bare endpoints and encode no `p`.
- `p` was called a coordinate until §2.

## 5. Hostile attack

| Question | Answer |
|---|---|
| Did the probability enter through the noise law? | **Yes**, always (A, A′, E-QSD). |
| Did `|α|²` enter through a supplied identification? | **Yes** (`p = |α|²`). |
| Change `b` at fixed σ: is Born destroyed? | **Yes** (B; theorem ⇒). |
| Change σ at fixed `b = 0`: is Born preserved? | **Yes** (A vs A′; theorem ⇐ for any admissible σ). |
| Absorption guaranteed or asymptotic? | Depends on σ near the boundaries: finite time (A) vs asymptotic (A′). Born holds in both. |
| Many microscopic parents, one hitting law? | **Yes.** `h` depends only on the scale function `s` — i.e. on `2b/σ²` — so every parent with the same `b/σ²` has the same weights, and every `b = 0` parent is Born. **Quotient variable: the scale function of the branch generator.** |

## 6. Experiment-P comparison and minimal sufficiency

| Ingredient | Supplies | Experiment P's unitary class | Earned in GRUT? |
|---|---|---|---|
| decoherence (dephasing) | loss of coherence | yes | no (borrowed, NO_GO 7) |
| absorption (non-unitary, `σ → 0` at the branches) | **selection** | absent (hence no outcome) | no |
| stochastic law `dW` | trajectory probabilities | absent | no |
| martingale branch coordinate (`b = 0`) | **Born fixation** (= decomposition independence) | n/a | no |
| `p = |α|²` | the number fed in | inherited from initial amplitudes | no (supplied) |

**Minimal sufficiency, as asked:**
- **Definite outcomes:** absorption plus a stochastic or deterministic driver.
- **Born frequencies:** in addition, the branch coordinate must be a martingale. That is the same as
  requiring outcome statistics to be linear in `ρ`.

**Classification:** REDISCOVERED-KNOWN (scale functions and hitting probabilities; Wright–Fisher
fixation; collapse-model martingale and no-signalling arguments). KNOWN-BUT-NEW-IN-GRUT: the
three-layer decomposition placed on the Experiment-P ledger, and the quotient (scale function) that
fixes the outcome law.

**Status: P-15 COMPLETE (one analytic pass). S2: UNFORMULABLE. Born ⟺ martingale ⟺
decomposition-independent (theorem, class-scoped, known). Nonlinear drift ⇒ SELECTION-NON-BORN.
σ = 0 ⇒ MEASURE-REQUIRED. Dephasing/stationary noise ⇒ NO-SELECTION. Same decoherence, different
selection (unraveling non-uniqueness). Every ingredient supplied.**
