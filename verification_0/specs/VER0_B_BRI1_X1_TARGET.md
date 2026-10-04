# VER0-B TARGET SPEC — BRI1-X1 (statement only)

**Contents:** definitions, model, preparation, clamp protocols, competitor-class definitions, the theorem statement,
assumptions and acceptance tests. **NO PROOF, NO COEFFICIENTS, NO NUMERICAL VALUES.**

**Source** (closed branch `grut-backreaction-identifiability-0 @ f2e6999c71b8b621868a019299085f64b80dce46`):

| file | blob |
|---|---|
| `backreaction_identifiability_0/BRI0_CHARTER.md` | `726aee36…` (class definitions as repaired: §4 clamp family, §R E₂±, §R2) |
| `backreaction_identifiability_0/BRI1_CANDIDATE_CHARTER.md` | `3a8bfb82…` (§1 model, preparation and protocols) |
| `backreaction_identifiability_0/BRI1_ANALYTIC_ESCAPE_THEOREM.md` | `ee50f05e…` (theorem statement as accepted, with §T4-R scope) |

**The orchestrator authored all of these: it is FULLY EXPOSED.** From the commit of this spec, all BRI0 / BRI1 files are
**SEALED** until the independent reproduction is committed.

## 1. Model

Retained coordinate q with momentum p_q, coupled reciprocally to N_B independent identical anharmonic oscillators
(x_j, p_j), j = 1 … N_B:

  H = p_q²/(2M) + V(q) + Σ_j [ p_j²/2 + x_j²/2 + x_j⁴/4 ] − N_B^{−1/2}·q·Σ_j x_j,  with V(q) = q²/2 + q⁴/4.

**Clamp.** The system trajectory q(·) is externally prescribed (V and M are irrelevant under the clamp). With
ε := N_B^{−1/2}, each bath oscillator obeys

  ẍ_j + x_j + x_j³ = ε·q(t),

and the environment force returned to q is **F_q(t) := ε·Σ_j x_j(t)**, i.e. −∂H_int/∂q.

**Preparation.** At q(0) = 0, the pairs (x_j(0), p_j(0)) are i.i.d. with density ∝ exp[−(p²/2 + x²/2 + x⁴/4)] (canonical
Gibbs at temperature 1). The preparation is identical for every protocol, with no re-preparation.

## 2. Protocol family and frozen protocols

**Clamp family 𝒳.** All C² trajectories q: [0, T] → ℝ with q(0) = 0, q̇(0) = 0, and q, q̇, q̈ bounded; T = 2π. One shared
environment preparation.

**Frozen protocols:**
- **P0 (rest):** q ≡ 0.
- **P1:** q(t) = s(t/π) for 0 ≤ t ≤ π, and q(t) = 1 for π < t ≤ 2π, where s(u) = 10u³ − 15u⁴ + 6u⁵.

## 3. Competitor classes

These are shared across all of 𝒳: they may depend on the prescribed path q_[0,t], but not on a protocol label.

| class | form | conditions |
|---|---|---|
| **E₁** | F_q(t) = M_t[q] + ξ(t) | M is an arbitrary deterministic causal functional; ξ is one exogenous process with one fixed, protocol-independent law |
| **E₂± (primary)** | F_q(t) = M_t[q] + G_t[q]·ξ(t) | M is deterministic causal; **G is deterministic causal with G_t[q] ∈ ℝ \ {0}** (sign changes allowed, zero not allowed); ξ is one shared, protocol-independent law. No nonlinear transformation of ξ |
| **E_univ** | F_q(t) = 𝔉_t[q, U] | U is one exogenous random object with a protocol-independent law; 𝔉 is causal in q |

**No restriction anywhere** to finite memory, Gaussianity, Markovity or stationarity. **Membership of E₂± is a property
of the whole family over 𝒳.**

## 4. Theorem to establish or refute

> **THEOREM (target).** There exist a non-empty small-time interval (0, δ), and for every fixed t in that interval a finite
> threshold N₀(t), such that for every integer N_B ≥ N₀(t) the interventional force family {F_q : q ∈ 𝒳} of this model
> lies **outside E₂±**.

**The intended witness** is the one-time standardised skewness, comparing P1 with P0. The reproducer must derive this
independently, including the order in N_B and the sign.

**Companion claim (reservoir limit), scoped.** For the frozen protocol family (and pointwise for any separately fixed
admissible bounded clamp, under the same estimates), the centred finite-dimensional force laws converge as N_B → ∞ to a
common Gaussian law, i.e. an **E₁-type** limit. **A uniform E₁ representation over all of 𝒳 is NOT claimed.**

**Earned claim if proved:** a finite reciprocal anharmonic environment can generate an interventional reduced-force law
outside the shared causal affine exogenous class E₂±, while the distinguishing non-affine witness vanishes in the
reservoir limit. **Not claimed:**
- primitive randomness;
- unique microscopic ontology;
- an escape from E_univ;
- GRUT-specific physics or predictions;
- TRUE COMPRESSION.

## 5. Acceptance tests (owner R1 – R9)

| test | required |
|---|---|
| **R1 — first-order response** | Derive the ε-variation equation at ε = 0. Establish the regularity to differentiate the finite-time flow in ε and to interchange with the Gibbs expectation. Obtain a single-oscillator third-cumulant expansion (κ₃(X^ε_q(t)) = ε K_q(t) + remainder) with a remainder strong enough for the aggregate result. **Derive the parity structure.** |
| **R2 — cumulant scaling** | From i.i.d. bath copies (F_q = ε Σ_j X_j^ε), derive the exact cumulant scaling, and the N_B-dependence of κ₃(F_q(t)). **Do not assume the exponent.** |
| **R3 — small-time sign** | For P1, derive the first non-zero small-time coefficient of the quantity controlling the third cumulant: which orders vanish, the first non-zero order, its sign, and the Gibbs-moment combination. **Prove the coefficient is strictly non-zero** from properties of the Gibbs measure. Symbolic algebra is allowed if written fresh, but the finite-order Taylor remainder must be justified analytically. |
| **R4 — sign interval** | Prove there exists δ > 0 such that the coefficient has a fixed non-zero sign on (0, δ). No numerical δ is required; **no scan over t**. |
| **R5 — P0 control** | κ₃(F_P0(t)) = 0 for every finite N_B and t. Identify the symmetry. |
| **R6 — non-degenerate variance** | Strict positivity of Var F at the chosen t, by a measure-theoretic / flow argument (not numerics). |
| **R7 — E₂± orbit escape** | From the §3 definitions, show signed affine modulation preserves zero vs non-zero \|standardised skewness\|. Then decide whether P0 = 0 and P1 ≠ 0 for all large finite N_B. If so, conclude the P1 law is not in the P0 signed-affine orbit, hence the family ∉ E₂±. |
| **R8 — large-N limit** | Handle the triangular array (ε = N_B^{−1/2}): covariance convergence; a uniform moment / Lindeberg condition; fdd Gaussian convergence. Scope as in §4. |
| **R9 — exact claim** | Confirm, correct or reject the §4 earned-claim wording. |

## 6. Firewall for the reproducer

**Allowed:** this spec; its own new files; standard mathematical knowledge and local libraries (sympy, numpy, scipy,
mpmath).

**Forbidden:**
- anything under `/home/user/BRI0` (all BRI0 / BRI1 files, scripts, logs and PF4Q material);
- any other branch or working tree;
- anything in `verification_0/` other than this spec and its own outputs;
- git;
- searching the filesystem.

**A file-access report is mandatory.** No frozen-τ or PF4Q values are used: the frozen-τ grade (PF4Q-I) is unaffected by
this reproduction.
