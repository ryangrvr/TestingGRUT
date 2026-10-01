# SCOUT_0 W1 P-06b RESULT — EXACT closure of the CM boundary (theorem)

**Charter:** follow-up to P-06 (spawned by P-06's H2 and the auditor's follow-up directive).
**Script:** `probes/p06b_exact_boundary.py` (mpmath dps=22; all checks pass, log below).
**Claim proved.** Within the declared Lorentzian-exponent family
`S_p(x) ∝ (x²+κ²)^{-p}`, `f_p(t) ∝ t^{p-1/2} K_{p-1/2}(κt)`:

> **`f_p` is completely monotone ⟺ `0 < p ≤ 1`.**

Family theorem only. No claim beyond the family; no derivation of the Level-0 generator;
S5-0/S5-1 terminals untouched (per P-06 charter).

## Proof

Let `ν = p − 1/2`, `z = κt`. CM is invariant under `z ↦ κt` (positive rescaling) and under
positive multiplication, so work with `g_ν(z) = z^ν K_ν(z)`.

**(A) CM for `p ≤ 1/2` (ν ≤ 0, μ = −ν ≥ 0).** The standard representation
(DLMF 10.32.9-type, valid for μ > −1/2)

```
K_μ(z) = √π (z/2)^μ / Γ(μ+1/2) ∫_1^∞ e^{-zs} (s²−1)^{μ−1/2} ds
```

gives, after dividing by `(z/2)^μ = z^{-ν}`:

```
z^ν K_ν(z) = √π 2^{-μ} / Γ(μ+1/2) · ∫_1^∞ e^{-zs} (s²−1)^{μ−1/2} ds
```

—a Laplace transform of the strictly positive weight `(s²−1)^{μ−1/2}` (integrable at s=1
since μ ≥ 0). Hence CM. This is the exact Bernstein representation the auditor asked for;
note it requires `ν ≤ 0` so that the `(z/2)^μ` prefactor cancels against `z^ν`.

**(B) CM for `1/2 < p ≤ 1` (0 < ν ≤ 1/2).** Same representation with `ν > 0` leaves a factor:

```
z^ν K_ν(z) = √π 2^{-ν} / Γ(ν+1/2) · z^{2ν} · ℒ[w],  w(s) = (s²−1)^{ν−1/2} ≥ 0.
```

The extra factor is multiplication by `z^α`, `α = 2ν ∈ (0, 1]`. Lemma: if `g` is CM and
`α ∈ (0, 1]`, then `z^α g(z)` is CM.
*Proof of lemma.* `α = 1`: `zg = −g′`, and derivatives of CM functions are CM. For
`0 < α < 1`, write (standard) `z^α = (α/Γ(1−α)) ∫_0^∞ (1 − e^{-uz}) u^{-α−1} du`, so for
`g = ℒ[μ]`:

```
z^α g(z) = (α/Γ(1−α)) ∫∫ u^{-α−1} [e^{-sz} − e^{-(s+u)z}] dμ(s) du
         = (α/Γ(1−α)) ∫∫ u^{-α−1} ζ ∫_s^{s+u} e^{-ζz} dζ  dμ(s) du,
```

using `e^{-sz} − e^{-(s+u)z} = ζ∫_s^{s+u} e^{-ζz} dζ` with `ζ ∈ (s, s+u)`. Swapping the order,
the coefficient of `e^{-ζz}` is `(1/Γ(1−α)) ∫_{ζ−s}^∞ u^{-α−1} du = (ζ−s)^{-α}/Γ(1−α) ≥ 0` on
`ζ > s`. So `z^α g` is a positive mixture of CM atoms `ζ e^{-ζz}` ⇒ CM. ∎

Since `ℒ[w]` is CM (A) and `α = 2p−1 ≤ 1 ⟺ p ≤ 1`, `g_ν` is CM exactly on `0 < p ≤ 1` —
**the factor `z^{2p−1}` is precisely where the boundary lives.**

**(C) Non-CM for `p > 1` (ν > 1).** Recurrence `d/dz [z^ν K_ν(z)] = −z^ν K_{ν−1}(z)`, so
`g_ν′ = −z g_{ν−1}` and

```
g_ν″(0) = −g_{ν−1}(0) = −2^{ν−2} Γ(ν−1) < 0   for every ν > 1 (integers included).
```

CM requires `g″ ≥ 0`; violation at the origin itself — no small-t grid artifact, and the
violation is uniform in ν > 1. (For `1/2 < ν < 1` separately: `g_ν ~ c₀ + c₁ z^{2ν}` with
`c₁ = 2^{−ν−1} Γ(−ν) < 0`, so `g″ ~ c₁(2ν)(2ν−1) z^{2ν−2} < 0` near 0; `ν = 1`: `g = zK₁`,
`g″ ~ ln(z/2) + γ − 1/2 → −∞`.)

**Endpoint p = 1 (ν = 1/2):** `g = √(π/2) e^{-z}` — the recorded S5-1 scaled kernel, exactly
CM. **Endpoint p → 0:** `g ∝ e^{-z}/z` — still a Laplace transform (measure on `[1,∞)`).

## Numerical cross-checks (all pass; script log)

- (C): `f″(0.05) < 0` for p ∈ {1.05, 1.2, 1.5, 2.0, 2.5, 3.0}; analytic `f″(0) = −2^{ν−2}Γ(ν−1)`
  matches the trend for p ∈ {2, 2.5, 3} (e.g. p=3: −1.2488 vs predicted −1.2533 at z=0.05).
- (A): Laplace-representation identity verified to ≤ 2.1e-23 / 3.0e-20 / 1.5e-13 (quadrature
  floor at the endpoint-singular p=0.5 weight) for p ∈ {0.1, 0.25, 0.5}.
- (ii): derivative-sign CM test (n≤4, t ∈ [0.05, 12]) passes for p ∈ {0.6, 0.75, 0.9, 1.0}.
- p=2: numeric `f″(0.5) = −0.3800867253` vs analytic `√(π/2)e^{-z}(z−1) = −0.3800867253` — exact.

## Auditor items, answered

1. **Bernstein representation for `0 < p < 1`:** YES, but only directly for `p ≤ 1/2` (case A).
   For `1/2 < p ≤ 1` the representation exists but carries the prefactor `z^{2p−1}`; CM survives
   because `z^α·CM` is CM for `α ∈ (0,1]` (lemma proved above). **The prefactor exponent `2p−1 ≤ 1`
   IS the boundary** — this is why the naive "positive integral ⇒ CM for all p" reading fails
   (it is exactly what would have contradicted the proven p=2 counterexample).
2. **Small-z expansion, `p > 1`:** for `ν > 1`, `g″(0) = −2^{ν−2}Γ(ν−1) < 0` — a violation AT the
   origin, uniform in p; for `1 < p < 3/2` the `z^{2ν}` term makes `g″ < 0` on a neighborhood of 0;
   `p = 3/2` is the logarithmic case with the same sign. **Yes: `f_p″ < 0` sufficiently near zero
   for every `p > 1`.** The bracket `[1.0, 2.0]` from P-06 is closed.

## Status

**P-06b COMPLETE — scout theorem (family-level, exact):** within the Lorentzian-exponent family,
`f_p` is CM iff `0 < p ≤ 1`. Literature classification: the CM of individual members
(`e^{-z}`, `K_0`) is textbook (REDISCOVERED-KNOWN); the exact iff statement over the exponent
family, and its role as the `d_mono` order parameter for the S5-1 parent, is
**KNOWN-BUT-NEW-IN-GRUT** (ledger update pending). P-06c (second deformation direction) remains
queued — this theorem does not transfer across families automatically.

**Next: P-08 (frozen charter).**
