# SCOUT_0 W1 P-06b CORRECTION 01 — proof repair (theorem stands; one leg retracted and replaced)

**Scope:** audit of the banked P-06b (`81ad661`, `probes/P06B_RESULT.md`) against the auditor's
analytic sketch. **Theorem statement unchanged. Case (C) unchanged. Case (B) retracted and replaced.**
S5-0/S5-1 terminals untouched; no claim beyond the declared Lorentzian-exponent family.

## Verdict in one line

> `f_p` is CM ⟺ `0 < p ≤ 1` is **true and now correctly proved**; the banked proof of the
> `1/2 < p ≤ 1` leg was **invalid** (it rested on a false lemma), and the auditor's original one-line
> route — which the banked text claimed only works for `p ≤ 1/2` — in fact covers all of `0 < p < 1`.

## 1. What was wrong in `81ad661`

**(B-lemma) "if `g` is CM and `α ∈ (0,1]`, then `z^α g(z)` is CM" — FALSE.**

- Smallest counterexample: `g ≡ 1` (CM; Bernstein measure `δ₀`). `z^α` is increasing, not CM.
- Nearest counterexample to the use made of it: `g = e^{-z}`, `α = 1`:
  `(z e^{-z})′ = e^{-z}(1 − z) > 0` on `(0, 1)`. Numerically `(z e^{-z})′(0.5) = +0.3033`.
- Where the written proof breaks: the step `α = 1: zg = −g′` is false
  (`−g′ = ℒ[s·dμ(s)]`, whereas `zg = ∫ z e^{-sz} dμ(s)`; `z e^{-sz}` is not CM in `z`).
  For `0 < α < 1` the step `e^{-sz} − e^{-(s+u)z} = ζ∫_s^{s+u} e^{-ζz} dζ` should read
  `z∫_s^{s+u} e^{-ζz} dζ` — the stray factor is `z`, not `ζ`, and it is exactly what destroys the
  "positive mixture of CM atoms" conclusion.
- The script's version of the lemma (`p06b_exact_boundary.py` line 8, `81ad661`):
  `z^α e^{-sz} = (1/Γ(1−α)) ∫_s^∞ ζ e^{-zζ} (ζ−s)^{-α} dζ` — also false. Exact RHS is
  `e^{-sz}[ s z^{α−1} + (1−α) z^{α−2} ]`. Checked: RHS/LHS = 1.5 at `(α,s,z)=(1/2,1,1)` and
  3.2653 at `(0.8, 2, 0.7)`.

Consequences for the banked text:
- "Auditor items, answered — item 1" (claiming the direct Bernstein representation exists "only
  for `p ≤ 1/2`" and that for `1/2 < p ≤ 1` it "carries the prefactor `z^{2p−1}`") is **wrong**.
- The interpretive claim "**the factor `z^{2p−1}` is precisely where the boundary lives**" and the
  FRONTIER_QUEUE line "mechanism = prefactor `z^{2p−1}` with `α = 2p−1 ≤ 1 ⟺ p ≤ 1`" are
  **unsupported** and are withdrawn. (The function `z^{2ν}·ℒ[w]` *is* CM for `0<ν≤1/2` — but only
  because it equals `g_ν`, which is CM by the correct argument below; nothing about `α ≤ 1` proves it.)

Why the error was not caught by the numerics: every numerical check in `81ad661` tested the
*conclusion* (which is true), never the lemma.

## 2. Replacement proof of the CM side — one formula for all `0 < p < 1`

Let `ν = p − 1/2`, `g_ν(z) = z^ν K_ν(z)`, and set `μ = −ν = 1/2 − p`. For `0 < p < 1`,
`μ ∈ (−1/2, 1/2)`. DLMF 10.32.8 (valid for `Re μ > −1/2`):

```
K_μ(z) = √π (z/2)^μ / Γ(μ + 1/2) · ∫_1^∞ e^{-zs} (s² − 1)^{μ − 1/2} ds .
```

Using `K_ν = K_{−ν} = K_μ`, the prefactor `(z/2)^μ = (z/2)^{−ν}` cancels against `z^ν` **for every
sign of `ν`**, and `μ − 1/2 = −p`, `μ + 1/2 = 1 − p`:

```
┌────────────────────────────────────────────────────────────────────────┐
│  z^ν K_ν(z) = √π 2^ν / Γ(1−p) · ∫_1^∞ e^{-zs} (s² − 1)^{-p} ds,   0<p<1  │
└────────────────────────────────────────────────────────────────────────┘
```

The weight `(s²−1)^{-p}` is strictly positive on `(1,∞)`, integrable at `s = 1` **iff `p < 1`**, and
exponentially damped at `∞`. So `g_ν` is the Laplace transform of a positive locally-finite measure
on `[1,∞)` ⇒ CM by Bernstein's theorem. This is exactly the auditor's formula. For `p ≤ 1/2` it is
literally the banked case (A); for `1/2 < p < 1` it is the *same formula* with `μ ∈ (−1/2, 0)` —
the banked restriction "requires `ν ≤ 0`" was spurious (DLMF only needs `μ > −1/2`).

**Endpoint `p = 1`** (`ν = 1/2`): `g = √(π/2) e^{-z}` — the recorded S5-1 scaled kernel; CM.
Consistently, as `p → 1⁻` the normalised weight `(s²−1)^{-p}/Γ(1−p)` concentrates at `s = 1`
(`1/Γ(1−p) → 0` while `∫` diverges at the endpoint), i.e. the Bernstein measure tends to `δ_1`,
which is the Bernstein measure of `e^{-z}`.

**Numerical confirmation of the formula** (`p06b_exact_boundary.py`, corrected; mpmath dps=22,
endpoint singularity removed by `s = 1 + t^{1/(1−p)}`): max relative mismatch over
`z ∈ {0.2, 1, 4}` is ≤ 3.8e-23 for `p ∈ {0.1, 0.25, 0.5, 0.6, 0.75, 0.9, 0.95, 0.99}`.

## 3. Non-CM side (`p > 1`) — unchanged

Case (C) of `81ad661` is sound and is kept verbatim: `d/dz[z^ν K_ν] = −z^ν K_{ν−1}` (DLMF 10.29.4)
gives `g_ν″(0) = −2^{ν−2}Γ(ν−1) < 0` for all `ν > 1` (`p > 3/2`); for `1/2 < ν < 1` the
`Γ(−ν) 2^{−ν−1} z^{2ν}` term (`Γ(−ν) < 0`) gives `g″ → −∞` as `z → 0`; `ν = 1` is the logarithmic
case with the same sign. Hence `f_p″ < 0` near `0` for every `p > 1`.

## 4. The boundary mechanism, correctly stated

Within the family, in the scaled rate `s = rate/κ`:

- `0 < p < 1`: positive, absolutely continuous Bernstein measure on `[1,∞)` with density
  `∝ (s²−1)^{-p}` — **the UV spectral exponent `p` is the endpoint exponent of that measure at
  the gap edge `s = 1`**;
- `p = 1`: the density becomes non-integrable, but the normalised measure has the singular endpoint
  limit `δ_1`, which is the Bernstein measure of the exponential — **still CM**;
- `p > 1`: CM fails (§3).

The boundary is a property of the Bernstein measure, not of a "prefactor exponent". (Auditor's
wording note, applied: do not say CM is "lost when the measure ceases to be locally finite at
`p ≥ 1`" — `p = 1` is CM.)

This is a family theorem. It does not derive the Level-0 generator, does not rescue S5-0, and
does not establish `p` as a universal GRUT order parameter (P-06c remains the test of
family-specificity).

## 5. Gate check (theorem gate)

- Assumptions: `S_p(x) ∝ (x²+κ²)^{-p}`, `κ > 0`, `p > 0`; GR transform table as in P-06.
- Quantified class: the one-parameter Lorentzian-exponent family (declared; no wider class).
- Conclusion: `f_p` CM ⟺ `0 < p ≤ 1`.
- Proof rather than examples: §2 (Bernstein rep + DLMF 10.32.8), endpoint, §3 (recurrence/small-z).
- Counterexamples when an assumption is removed: the `z^α`-lemma counterexample (§1) shows
  CM is **not** preserved under the prefactor manipulation the banked proof relied on.
- Already known? The CM range of the Matérn-type function `t^ν K_ν(t)` is standard material in the
  positive-definite-function literature; classification stays **KNOWN-BUT-NEW-IN-GRUT**, with an
  explicit literature citation still owed (only DLMF 10.32.8 / 10.29.4 / 10.30 are load-bearing here).
- Math vs GRUT interpretation: separated (§4 last paragraph).

## 6. Files changed by this correction

- `probes/P06B_RESULT.md` — case (B) replaced by the unified representation; auditor item 1 and
  the status paragraph rewritten; pointer to this memo added.
- `probes/p06b_exact_boundary.py` — false lemma removed from the proof comment; check (ii) now
  verifies the unified representation across `0 < p < 1` (including `p = 0.9, 0.95, 0.99`) and
  records the lemma counterexample; all other checks unchanged.
- `FRONTIER_QUEUE.md` — P-06b row and Wave-2 line corrected.
- `COUNTEREXAMPLE_LEDGER.md` — CE-02 added (the false lemma).
