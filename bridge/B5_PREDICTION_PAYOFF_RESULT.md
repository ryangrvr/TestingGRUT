# B5 — PREDICTION / PAYOFF FIREWALL (result)

**Question (owner).** After every input Bridge-1 still calls supplied is allowed to vary, does any non-trivial
observable relation remain forced, and is it distinctive to GRUT rather than ordinary harmonic-network / open-system
mathematics?

**Script:** `bridge/b5/b5_payoff.py` (+ `b5_payoff.log`), seed 20261005. Label: **independent code path, not independent
reviewer.** Canonical GRUT was read-only.

**Variation ledger:** `B5_VARIATION_LEDGER.md`.

**Baseline before B5:** ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS.

## B5-0 Preregistered success bar (all seven required)

1. **Forced** (from earned / Bridge-reviewed structure, not inserted).
2. **Residual-input invariant.**
3. **Convention invariant**: access labelling, closure rule, threshold, sampling, sign convention, normalization, epoch
   coordinate.
4. **Non-trivial**: not energy conservation, covariance propagation, a matrix identity, locality consistency, Gibbs
   stationarity or standard harmonic-network mathematics.
5. **Distinctive**: not generic in the ordinary comparison class.
6. **Observable**: an operationally specified quantity.
7. **Falsifiable.**

The bar was fixed before any candidate was evaluated (owner text). Candidates were classified mechanically, without
ranking (B5-14).

## Candidates

### P1 — R1 frame overdetermination: **CONSISTENCY-ONLY / STANDARD-STRUCTURE**

- **Test A.** When the net, the on-site drift and the noise are declared in one site frame, the reconstructed frames
  agree exactly (overlap 1.000000). The drift axes are read from the Hessian of the quartic force at a probe state; the
  noise axes are the eigenframe of Q.
- **Test B.** Writing the drift on-site in an **independently rotated** frame gives an admissible model: V is still
  convex, and the flow relaxes from |x| = 4.15 to 7e-4. Its drift axes then match the rotated frame (overlap 1.000000) and
  not the net frame (0.349). **Nothing earned forbids the independent rotation.**
- **Test C.** The ordinary nonlinear Langevin network with site-local V₄ and noise is the same equation, so the agreement
  is generic.
- **B5-11 firewall:** neither side is derived independently of the other. Agreement is a consistency condition between
  two supplied declarations, **not a prediction.**

Fails criteria 1, 2 and 5.

### P2 — generalized integrated transfer: **STANDARD-STRUCTURE (identity) / INPUT-DEPENDENT (value)**

**Provenance.** The equations of motion give d(E₁ + E_int)/dt = −g·q₁p₂. So X(T) = E₁(0) − E₁(T) + g(Q12(T) − Q12(0)) for
**any** quadratic chain and any state. This is **energy accounting**: X(T) is exactly the bath self-energy gained. The
T → ∞ form then adds local return to equilibrium (LS-1 / B2-P1).

**Variation.**
- Across 5 chains (κ, g, c varied) and 4 preparations (T_s, T_b, α varied), the closed form matches to ≤ 1.5e-4, and the
  **value moves with every supplied input.**
- Outside the absolutely continuous class (κ = 1.0, with a bound state at λ = 0.231), X(T) has no limit (it wanders over
  [1.35, 2.00]).

**Elimination attempt.** As an identity among measured quantities, the relation is already parameter-free, because it
*is* energy conservation (fails criterion 4). As a number, it cannot be freed of H_marginals and H_cross (fails 2;
B5-12: CONDITIONAL-RELATION at best).

**Comparator.** A random 30-node harmonic network with a 3-site system satisfies ∫ bath flux = ΔE_B to 1e-6.

### P3 — equal-temperature offset ½T_b r²: **INPUT-DEPENDENT + STANDARD-STRUCTURE**

**General closed form** (equal T, product preparation, q₁–q_B cross α × Gibbs cross):

    X(∞) = g·Q12_G (½ − α) = −⟨E_int⟩_G (½ − α)

This was confirmed for κ ∈ {2.3, 3.0, 2.6}, g ∈ {1, 0.6, 1.2} and c ∈ {2.3, 2.5, 2.8} (ratio ±0.4998 … ±0.5000).

| part | verdict |
|---|---|
| **the coefficient ½** | the static quadratic-network identity E_S(prod) − E_S(Gibbs) = ½⟨E_int⟩_G. It is exact in random 25-node networks with a 4-site system (−1.112204 = −1.112204, etc.): equipartition for quadratic Hamiltonians. **Standard.** It survives only at α = 0 |
| **r²** | (K⁻¹)₁₂ = 0.339, 0.241, 0.167, 0.319, 0.224 across the chains. **Not universal** |
| **sign** | reverses for α > ½ |

### P4 — A_closure compression: **NO-OBSERVABLE**

The closure dimension depends on the seed and on R_closure (B3-1: 16 vs 4; classical 24 / 12 / 2). No operational
measurement turns an algebra dimension into a number without a further supplied readout. Bookkeeping only.

### P5 — access-relative geometry: **INPUT-DEPENDENT (access-priced)**

Holding D fixed and varying the seed, partition and resolution changes what geometry is recoverable:
- unique under full access;
- a non-isometric family under single-site access;
- order-relative topology (GS1 canonical, B3-4, B3-7).

No geometry invariant survives every access choice.

### P6 — correlation-boundary sign: **BOUNDARY-STATE-DEPENDENT** (strong negative)

A full sweep of the PSD-admissible q₁–q₂ correlations at the S6 marginals uses extremal rank-1 crosses saturating
|Q12(0)| ≤ √(Q11·Q22), with min eigenvalues checked.

| (T_s, T_b) | X(∞) range | sign |
|---|---|---|
| (1, 1) | **[−0.334, +0.673]** | both signs |
| (½, 1) | [−0.686, +0.025] | both signs |
| (2, 1) | [+0.458, +1.881] | fixed |
| (10, 1) | [+7.58, +10.76] | fixed |

- **No parameter-free arrow sign exists.**
- Where the sign is fixed, it is fixed only by a Cauchy–Schwarz / PSD bound with supplied marginals. That is a generic
  positivity bound and is excluded by the bar.
- This matches the known open-system result that initial correlations can reverse heat flow, up to bounds set by the
  correlations. **KNOWN-RESULT IMPORT** (Partovi 2008; Jennings & Rudolph 2010; Micadei et al. 2019), cited from memory
  and not re-fetched in this environment.

### B5-9 Dimensionless combinations

| ratio | values | verdict |
|---|---|---|
| X/(g·Q12_G) | ½ only at (α = 0, T_s = T_b); −0.200 at α = 0.7; 3.45 / 4.65 / 10.48 at (2, 1) | killed |
| X/(T_s − T_b) | 1.169 / 1.120 / 1.050 | K-dependent |

Only combinations suggested by the exact equations were tested (no fishing).

## Mechanical scorecard (criteria 1 – 7)

| cand. | 1 forced | 2 input-inv. | 3 conv.-inv. | 4 non-trivial | 5 distinctive | 6 observable | 7 falsifiable | terminal |
|---|---|---|---|---|---|---|---|---|
| P1 | ✗ (declared) | ✗ | ✗ (frame declaration) | ✗ | ✗ | ✓ | ✓ | CONSISTENCY-ONLY |
| P2 | ✓ | ✗ (value) | ✓ | ✗ (energy conservation) | ✗ | ✓ | ✓ | STANDARD-STRUCTURE |
| P3 | ✓ | ✗ | ✓ | ✗ (equipartition) | ✗ | ✓ | ✓ | INPUT-DEPENDENT |
| P4 | ✓ | ✗ | ✗ (R_closure) | — | — | ✗ | ✗ | NO-OBSERVABLE |
| P5 | ✓ (given access) | ✗ | ✗ (access labelling) | ✗ | ✗ | ✓ | ✓ | INPUT-DEPENDENT |
| P6 | ✓ | ✗ (both signs) | ✓ | ✗ (PSD bound) | ✗ | ✓ | ✓ | BOUNDARY-STATE-DEPENDENT |

**No candidate passes all seven.**

**B5-13 (empirical mapping)** is therefore not triggered. No platform was invented.

## Terminal

> **NO DISTINCTIVE PAYOFF FOUND AT CURRENT BRIDGE-1 PREMISE ENVELOPE.**
> **ZERO CONFIRMED DISTINCTIVE GRUT QUANTITATIVE PREDICTIONS** (preserved).

## Scope

- classical Gaussian harmonic chains and networks, plus one nonlinear on-site class (P1);
- finite-N cross-checks of infinite-N closed forms;
- no lift / ħ / outcome / gravity consumed.

A distinctive payoff outside this envelope (quantum lift, non-Gaussian bath, gravity sector) is **not excluded by B5**.
It is also not suggested by any Bridge-1 relation.
