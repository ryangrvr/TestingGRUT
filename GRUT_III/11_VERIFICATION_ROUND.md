# 11 — Sixth-round verification: independent re-derivation + adversarial audit of `v8`, `v11`, `v12`, `v13`, `v14`

**Method.** A workflow of agents, two per script: lens A reads only the module docstring, re-derives the claims in its own
code, then runs the script and compares; lens B reads the code and the citing documents and tries to refute them
(mutation tests, convention and branch checks, scope of statements). A second workflow put three skeptics (mathematics,
physics/scope, numerics) on `v14`. Results are structured (verdict, checks, issues with evidence, overclaims) and were
acted on by the author; this file is the record.

**Status (2026-10-12).** First launch: 3 of 10 round-6 agents and 0 of 3 `v14` agents ran; the rest were refused at
launch by a session usage limit. Both workflows were relaunched after the limit reset (completed agents replay from
cache). This file records the three completed audits and the fixes made; it is updated when the relaunch reports.

## Completed audits

### `v8_residual_framework.py` — lens A: CONFIRMED; lens B: PARTIAL (computation right, documentation overstated)

Independently reproduced (lens A, own sympy/mpmath, before reading the code): the FDT identity; the pointwise improved
stress-tensor correlator ⟨T₀₀T₀₀⟩_c = 1/(4π⁴(t−iε)⁸) in Minkowski by Wick contraction; its Weyl map to
H⁸/(1024π⁴ sinh⁸(H(t−iε)/2)); the recursion for the Fourier transforms of sinh^{−2n} and hence
G⁺_T(ω) = ω(ω²+H²)(ω²+4H²)(ω²+9H²)/(10080π³(1−e^{−2πω/H})) exactly (direct quadrature agrees to 15 digits); the odd polynomial
ρ_T, N_T(0) = H⁷/560π⁴, the local exponents 1 → 7 with the script's intermediate slopes; the kms_defect helper's linear
detection of planted non-KMS perturbations; every entry of the §C tables to all printed digits. Lens B reproduced ρ_T from
the short-distance Laurent structure alone and G⁺_T by independent quadrature to 30 digits, and confirmed the §C thresholds
are live (perturbed models fail).

Issues found and **fixed** (commits 271319a, and this one):
- `07` §5 misquoted `v8` §C: "low-frequency exponents 3 and 0–1 against 7; a log-log slope residual of −4 for ω³". The
  medium's low-frequency exponent is 1 (7 is the ultraviolet one) and −4 was only the helper's self-test on pure power
  laws. Rewritten with the numbers the script prints.
- `07` §4 described the `v8` §D template as taking "ρ_f^X (or state, or κ_n)"; it takes only a model ρ and forces N by
  KMS. Restated as the ρ-level template; state-level (T2 type (b)) and cumulant claims need their own N or κ_n.
- `v8`'s printed conclusion said the old kernels differ from the medium "of the same (observable, state)" while its own
  template labels both "unspecified". Rewritten: they are read as ρ-shapes; C1's ω_c is set to H by fiat.
- The (D) check "implied noise KMS-consistent by construction" is a tautology; relabelled as helper self-consistency.
- Docstring: "NEVER ω³" (the local exponent passes through 3 at ω = 1.21H) and the shape_residual "iff" softened; the
  EMPTY verdict is now stated as grid-limited; the "nothing is left to supply" line is scoped (two-point level, this
  observable and state); the two-point/Gaussian scope added to the thesis.
- `07` T5: sign of the 1/(2 sinh²) term (1 − Z = −2 sinh²(HΔt/2) for timelike separation).
- `05` O5 "The rule": KMS fixes N given ρ; conventional operator content is what fixes ρ — both clauses now stated.
- `02` D2 attributed to `v8` B checks that only `v1` performs; split.

Open (minor, recorded, not changed): the Drude reading's shape test is at ω = 1H, 10H rather than 0.1H, 10H (disclosed in
the docstring; a per-model criterion); the (A) mutation test's threshold is at a single ω = 0.5H; C2's odd extension gives
N ~ 1/ω (a logarithmically divergent equal-time variance), which is one more reason it is "a reading", not the record's
object.

### `v11_quasi_ds_kms_defect.py` Part B (general a(t)) and Part D (ΛCDM) — lens A: CONFIRMED with a scope correction

Independently derived (own functional for ln a(t̄+δ) = Hδ + f(δ)): the Ḣ term alone gives ln R = −2ε[(x/2)coth(x/2) − 1],
exactly the shift H → H(1−ε) = ä/ȧ, remainder identically zero (symbolic); the Ḧ term gives a non-shift residual, so the
theorem needs Ḧ/H³ = ε(2ε − η) second order, i.e. **η first-order small too**. Numerics confirm f₂(x) (power law
p = 20…320), the O(ε) residual if H̄ is used, the εη coefficient to all digits, Part C's D₂(ω) to 3–5e-5, and Part D's
ΛCDM numbers (ε₀ = 0.4725, H_eff = 0.5275H₀, deviations 0.12%, 0.56%, 3.64%, 19.2%) exactly.

Issues found and **fixed** (commit 958bdb7):
- "Theorem (first order, any a(t))" dropped the hypothesis on η. Counterexamples: a = exp(t + At³) (ε = 0, first-order
  non-thermal residual ∝ A); a = exp(t + A sin t) at t̄ = 1; and ΛCDM today (η₀ = −2.06). The heading of D11, `09` §1,
  `05` O5 and the `v11` docstring now state the hypothesis.
- The `v11` docstring said Part D computes a Wigner spectrum at t₀; it computes only the time-domain comparison (the
  worldline exists only for |Δ| < 2t₀). Corrected.
- D11's "factors 4–18 above the second-order prediction" used the η = 0 form; with the εη term the factors are 4–22;
  qualified (conclusion unchanged).
- R11 kept the superseded revision-2 numbers without an in-place pointer to R15; marked.
- "Exactly KMS to first order" is a statement at fixed x = HΔ, not uniform for ω ≲ εH; the clause was added to D11.
Open (minor): Part B's first-order identity is checked at random parameter points rather than simplified symbolically.

### `v11` Part D lens B, `v12` (both lenses), `v13` rev. 2 (both lenses), `v14` (three lenses)

Not yet run (refused at launch; relaunched). Until they report, `v12`, `v13` rev. 2 and `v14` remain author-level, as
`08` states.

## What the round changes

Nothing in the physics of the branch was refuted: `v8`'s identities and spectra and `v11`'s first-order theorem were each
re-derived independently and agree to all digits. What the audit corrected is the prose around them — four misquotations
or overstatements in `07`/`05`/`02`/`09`, one missing hypothesis in a theorem statement, one stale docstring — all now
fixed. The pattern is the same as in rounds two to five: the computations hold; the claims made about them had to be cut
back to what was computed.
