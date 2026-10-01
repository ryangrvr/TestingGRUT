# SCOUT_0 W1 P-06b RESULT — EXACT closure of the CM boundary (theorem)

**Charter:** follow-up to P-06 (spawned by P-06's H2 and the auditor's follow-up directive).
**Script:** `probes/p06b_exact_boundary.py` (mpmath; all checks pass, log below).
**Correction history:** `81ad661` banked this theorem with an invalid proof of the `1/2 < p ≤ 1` leg
(false lemma). Repaired in `P06B_CORRECTION_01.md`; this file is the corrected statement.
**Claim proved.** Within the declared Lorentzian-exponent family
`S_p(x) ∝ (x²+κ²)^{-p}`, `f_p(t) ∝ t^{p-1/2} K_{p-1/2}(κt)`:

> **`f_p` is completely monotone ⟺ `0 < p ≤ 1`.**

Family theorem only. No claim beyond the family; no derivation of the Level-0 generator;
S5-0/S5-1 terminals untouched (per P-06 charter).

## Proof

Let `ν = p − 1/2`, `z = κt`. CM is invariant under `z ↦ κt` (positive rescaling) and under
positive multiplication, so work with `g_ν(z) = z^ν K_ν(z)`.

**(A) CM for all `0 < p < 1` — one Bernstein representation.** Put `μ = −ν = 1/2 − p ∈ (−1/2, 1/2)`.
DLMF 10.32.8, valid for `Re μ > −1/2`:

```
K_μ(z) = √π (z/2)^μ / Γ(μ+1/2) ∫_1^∞ e^{-zs} (s²−1)^{μ−1/2} ds
```

With `K_ν = K_{−ν} = K_μ`, the prefactor `(z/2)^μ = (z/2)^{−ν}` cancels against `z^ν` for either
sign of `ν`, and `μ − 1/2 = −p`, `μ + 1/2 = 1 − p`:

```
z^ν K_ν(z) = √π 2^ν / Γ(1−p) · ∫_1^∞ e^{-zs} (s²−1)^{-p} ds ,      0 < p < 1 .
```

The weight `(s²−1)^{-p}` is strictly positive, integrable at `s = 1` iff `p < 1`, and exponentially
damped at `∞` — a Laplace transform of a positive locally-finite measure on `[1,∞)`. Hence CM
(Bernstein). This is exactly the auditor's formula; for `p ≤ 1/2` it coincides with the banked case
(A) of `81ad661`, and for `1/2 < p < 1` it is the same formula with `μ ∈ (−1/2, 0)`, which DLMF
permits.

**Endpoint p = 1 (ν = 1/2):** `g = √(π/2) e^{-z}` — the recorded S5-1 scaled kernel, exactly CM.
As `p → 1⁻` the normalised weight `(s²−1)^{-p}/Γ(1−p)` tends to `δ_1`, the Bernstein measure of
`e^{-z}` — the endpoint is the continuous limit of (A). **Endpoint p → 0:** `g ∝ e^{-z}/z` — still a
Laplace transform (measure on `[1,∞)`).

**(B) — retracted.** `81ad661` proved `1/2 < p ≤ 1` via a prefactor `z^{2p−1}` and the lemma
"`z^α·CM` is CM for `α ∈ (0,1]`". That lemma is false (`z e^{-z}` is not CM) and the leg is
superseded by (A). See `P06B_CORRECTION_01.md` §1.

**(C) Non-CM for `p > 1` (ν > 1/2).** Recurrence `d/dz [z^ν K_ν(z)] = −z^ν K_{ν−1}(z)`
(DLMF 10.29.4), so `g_ν′ = −z g_{ν−1}` and

```
g_ν″(0) = −g_{ν−1}(0) = −2^{ν−2} Γ(ν−1) < 0   for every ν > 1 (integers included).
```

CM requires `g″ ≥ 0`; violation at the origin itself — no small-t grid artifact, and the
violation is uniform in ν > 1. (For `1/2 < ν < 1` separately: `g_ν ~ c₀ + c₁ z^{2ν}` with
`c₁ = 2^{−ν−1} Γ(−ν) < 0`, so `g″ ~ c₁(2ν)(2ν−1) z^{2ν−2} < 0` near 0; `ν = 1`: `g = zK₁`,
`g″ ~ ln(z/2) + γ − 1/2 → −∞`.)

**Where the boundary lives.** The Bernstein (relaxation-rate) measure of `f_p` has density
`∝ (s²−1)^{-p}` on `[1,∞)` in the scaled rate `s = rate/κ`: the UV spectral exponent `p` is the
endpoint exponent of that measure at the gap edge `s = 1`. CM is lost exactly when the measure stops
being locally finite there (`p ≥ 1`), and (C) independently certifies `f″ < 0` beyond. (The earlier
"prefactor `z^{2p−1}`" reading of the boundary is withdrawn.)

## Numerical cross-checks (all pass; script log)

- (A): unified Laplace (Bernstein) representation verified to ≤ 3.8e-23 relative (dps=22) over
  `z ∈ {0.2, 1, 4}` for `p ∈ {0.1, 0.25, 0.5, 0.6, 0.75, 0.9, 0.95, 0.99}` (endpoint singularity
  removed by `s = 1 + t^{1/(1−p)}`).
- Lemma counterexample recorded: `(z e^{-z})′(0.5) = +0.3033 > 0`; the script-line-8 identity of
  `81ad661` has RHS/LHS = 1.5 at `(α,s,z) = (1/2, 1, 1)`.
- (C): `f″(0.05) < 0` for p ∈ {1.05, 1.2, 1.5, 2.0, 2.5, 3.0}; analytic `f″(0) = −2^{ν−2}Γ(ν−1)`
  matches the trend for p ∈ {2, 2.5, 3} (e.g. p=3: −1.2488 vs predicted −1.2533 at z=0.05).
- Corroboration: derivative-sign CM test (n≤4, t ∈ [0.05, 12]) passes for p ∈ {0.6, 0.75, 0.9, 1.0}.
- p=2: numeric `f″(0.5) = −0.3800867253` vs analytic `√(π/2)e^{-z}(z−1) = −0.3800867253` — exact.

## Auditor items, answered (corrected)

1. **Bernstein representation for `0 < p < 1`:** YES, directly and uniformly — the auditor's formula
   `z^ν K_ν(z) ∝ ∫_1^∞ e^{-zs}(s²−1)^{-p} ds` holds on the whole interval (DLMF 10.32.8 with index
   `1/2 − p ∈ (−1/2, 1/2)`). The `81ad661` claim that it works "only for `p ≤ 1/2`" was wrong. The
   representation fails at `p = 1` because the weight stops being integrable at `s = 1`; `p = 1` is
   the explicit exponential endpoint.
2. **Small-z expansion, `p > 1`:** for `ν > 1`, `g″(0) = −2^{ν−2}Γ(ν−1) < 0` — a violation AT the
   origin, uniform in p; for `1 < p < 3/2` the `z^{2ν}` term makes `g″ < 0` on a neighborhood of 0;
   `p = 3/2` is the logarithmic case with the same sign. **Yes: `f_p″ < 0` sufficiently near zero
   for every `p > 1`.** The bracket `[1.0, 2.0]` from P-06 is closed.

## Status

**P-06b COMPLETE — scout theorem (family-level, exact), proof repaired in CORRECTION 01:** within
the Lorentzian-exponent family, `f_p` is CM iff `0 < p ≤ 1`. Literature classification: the CM of
individual members (`e^{-z}`, `K_0`) is textbook (REDISCOVERED-KNOWN); the exact iff statement over
the exponent family is standard Matérn-type material and its role as the `d_mono` order parameter
for the S5-1 parent is **KNOWN-BUT-NEW-IN-GRUT** (explicit citation still owed; only DLMF
10.32.8 / 10.29.4 / 10.30 are load-bearing here). P-06c (second deformation direction) remains
queued — this theorem does not transfer across families automatically.

**Next: P-08 (frozen charter) — already executed (`cf0f7fe`).**
