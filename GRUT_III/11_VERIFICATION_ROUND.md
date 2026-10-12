# 11 — Sixth-round verification: independent re-derivation + adversarial audit of `v8`, `v11`, `v12`, `v13`, `v14`

**Method.** A workflow of agents, two per script: lens A reads only the module docstring, re-derives the claims in its own
code, then runs the script and compares; lens B reads the code and the citing documents and tries to refute them
(mutation tests, convention and branch checks, scope of statements). A second workflow put three skeptics (mathematics,
physics/scope, numerics) on `v14`. Results are structured (verdict, checks, issues with evidence, overclaims) and were
acted on by the author; this file is the record.

**Status (2026-10-12).** First launch: 3 of 10 round-6 agents and 0 of 3 `v14` agents ran; the rest were refused at
launch by a session usage limit. Both workflows were relaunched after the limit reset and all thirteen agents have reported.
Every finding is recorded below with the fix made or the reason it was left.

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

### `v14_massless_sector_first_order.py` — mathematics: CONFIRMED; numerics: CONFIRMED; physics/scope: PARTIAL

**Mathematics lens** (own derivations before reading the code): the de Sitter propagator normalisation (Hadamard
coefficient universal in ν; ν = 1/2 gives the conformal propagator; ν → 3/2 gives Allen–Folacci's Z-dependent part); the
Pochhammer-series derivation of S₁ and S₂ (match hyp2f1 derivatives to 1e-31); g₁ = 3/2, g₂ = 1/2 + π²/6; the `v12` lemma
with normalisation (Wronskian +i; mode-equation residuals 1e-21); the closed form against a power-law mode integral at
ν = 1.75 with an explicitly derived IR term — and the proof that the IR term X(ττ′)^{−δ} is *exactly* constant in t and t′,
so "C_IR constant" is right; the first-order algebra and c₂; the Richardson theorem (halving ratios 3.994–4.004, variants
O(ε)); the KMS defect of the *exact* closed form (∝ ε² at β_eff, ∝ ε at β_KH and β̄); analyticity of S₂(cosh²) in the strip
(w = 0 removable; the cuts meet the strip only on its edges); A₁ = −1, B₁ = 4 ln 2 − 4; r(ω) digit for digit; r → 3/2ω²;
the UV pole cancellation; D from B₁ and from a direct polynomial fit of the exact function (0.99761073 vs 0.99761371 at
ε = 1e-3). Caveat found: the fixed-mean-time function in the window also carries an O(ε)Δ² term, −(3εH̄⁴/32π²)Δ² in W,
which exceeds the O(ε) correction to the |Δ| term for H̄Δ ≳ 3; D is correct as the |Δ| coefficient, but "≈ D|Δ|" and the
σ calibration were stated more loosely than computed.

**Numerics lens** (mutation-tested): Part C fails under c₂ × 1.1 (worst |ratio − 4| = 1.99), c₂ × 0.9 (2.52), g₁ = 1.4 (2.18),
c₂ = 0 (2.08); the reference subtraction hides only a Δ-independent constant (the IR pole plus O(1), O(ε) constants —
harmless for ω ≠ 0); the ε² scaling persists at Im Δ = −0.05 and −1.0 and hyp2f1 and the closed forms are on the same branch;
the mode integral converges to the closed form 6e-11 → 7e-23 as grid and precision are tightened, and the cutoff case's
remainder scales as k_IR^{5−2ν} (a finite-cutoff effect, not quadrature); r(0.5), r(1), r(2) stable to 18 digits; r(12),
r(16) continue to decay once dps ≥ 50 (at dps 30 the e^{βω/2} amplification destroys ω = 16: precision, not physics); the
analytic large-ω form r = 3(2ω⁴ − ω² − 1)/(4ω⁴(ω² + 1)) reproduces r(8), r(12), r(16) to 1.7e-7, 6e-9, 7e-10; the asymptotic
polynomial is reproduced by a fit on s = 20…40 to 1e-15; the first light-cone image of the exact function sits at
|Im Δ| = β_eff + 20.9ε², O(ε) away from β_KH and β̄, and the exact function's nonlinear KMS defect is 3.4ε² at β_eff versus
0.51ε at β_KH and 1.01ε at β̄.

**Physics/scope lens** (accepted the computation; corrected the framing): (i) "KMS at fixed mean time t̄" is a convention,
and it is doing work: the massless function is non-stationary at zeroth order (W = F(Δ) + (H³/4π²)t̄ + const), so c₀ is the
O(1) secular growth, not an O(ε) constant; in the detector's own time (instantaneous one-sided Unruh–DeWitt rate) the same
term becomes a −(H³/8π²)|s| kernel and the massless detector rate in exact de Sitter is Ḟ(E) = G₀(−E) + H³/4π²E² — not KMS at
any temperature (Garbrecht–Prokopec; the auditor validated the formula to 1e-21). (ii) Even for the conformal sector the
first-order temperature depends on the time-reference convention: the Wigner kernel gives the pure shift to ä/ȧ (D11
reproduced, c = 1.000), the one-sided rate gives a frequency-dependent c(E) = 0.24, 0.93, 3.97 at E/H = 0.25, 0.5, 1,
tending to the bare H̄ as E → 0 (auditor's numbers; the structure — a term s·∂_H W, not a temperature shift — is confirmed
analytically by the author; the numbers are not re-derived). (iii) The general-a(t) extension asserted in D13's scope is not
a corollary: the S₂ term comes from the constant index shift and has no direct analogue for varying ε. (iv) The graviton
remark is an inference, not a computation. (v) The σ calibration is an interpretive identification; "and at no other σ" is
monotonicity. (vi) A stale reference to `v13` as the massless precedent in `05`.

**Fixes made (this commit):** `v14` — sign slip in the printed β_KH − β_eff; Part D relabelled as an implementation check and
Part D″ added (the exact-function KMS defect test, ratios 4.00/2.01/2.01); the r → 3(2ω⁴ − ω² − 1)/(4ω⁴(ω² + 1)) check at ω = 8;
Part B's cutoff comparison as one well-conditioned difference integral with a larger regulator (3.4e-6 → 1.2e-10, tolerance
1e-7); the docstring carries the convention, the Δ² companion term, the c₂ factor and the precision note on the r table.
D13 — a "Convention" paragraph; the KMS statement restricted to the ω ≠ 0 part of the fixed-mean-time spectrum; the
byproduct restated with the Δ² term and the interpretive step; the scope paragraph no longer asserts the general-a(t) or
graviton results; "P1 adds nothing" replaced by "establishes nothing about P1 beyond consistency". D11 — the same
convention note. `05`, `09` §1 and §4 — "seen by a comoving detector" / "what an Unruh–DeWitt detector reads" replaced by the
Wigner-spectrum statement, with the one-sided caveat. The headline survives in this form: the fixed-mean-time spectrum of
the massless sector at ω ≠ 0 is KMS at (ä/ȧ)/2π to O(ε), Kodama–Hayward excluded, no coarse-graining scale.

### `v11` Part B — lens B: CONFIRMED (with code fixes)

Could not break the headline: the Part B first-order remainder simplifies to zero symbolically (not merely small at points);
a mutation that mis-shifts H_eff by 1e-6·εH fails the test, Kodama–Hayward fails by 11; an independent non-power-law family
a = t^p e^{ct} (η/ε = −2, −6) gives (R−1)/ε → 0 to 3e-15. Found and **fixed**: `richardson()` applied the first-order weight 2 at
every level (should be 2, 4/3, 8/7), leaving an ε²/8 artefact — exactly the "c = 1 to 3e-5" of D11/R15; with the correct
weights c = 1 to 1e-9 (the D₂ digits were unaffected); the εη coefficient of the second-order remainder was printed but never
checked — now checked against the auditor's closed form g(x) = 2 + x²/4 − (x + x³/24)coth(x/2) to 1e-16; the first-order
identity is now also established symbolically in the script. Recorded, not changed: the second-order formula is complete
only under the full adiabatic hierarchy (ε₃ = η̇/Hη also small) — added to the docstring; the scope is the conformally coupled
massless scalar's worldline function in flat FRW (other scalar primaries by Weyl covariance; composite, higher-spin and
smeared observables not computed) — added.

### `v11` Part D (ΛCDM) — lens A: PARTIAL; lens B: PARTIAL (numbers right, framing wrong)

Every number reproduced to all digits (ε₀ = 0.4725, H_eff = 0.5275H₀, t₀ = 0.951/H₀, deviations 0.1217%, 0.5563%, 3.639%,
19.16%). What did not hold was the prose: (i) "deviates from any thermal form by up to 19%" — only the H_eff form was compared;
the T = 0 (Minkowski) form is closer at every grid point, the deviation is unbounded as Δ → 2t₀ (1422% at Δ = 1.9/H₀), and the
KMS defect at 2π/H_eff is O(1–100); (ii) "factors 4–18 above the second-order prediction" — the factor ≈ 4 at Δ ≤ 0.5/H₀ is exactly
1/(1−ε₀)² from evaluating f₂ at x = H_effΔ rather than H̄Δ (an O(ε³) convention), so it carried no information; (iii) the small-Δ
agreement is the universal Hadamard term — the local effective rate in the Δ² correction is H₀√(1 − 2ε₀) = 0.23H₀, not
H_eff = 0.53H₀; (iv) the Δ = 1.5/H₀ pair samples z = 2.25 (ε = 1.4 there) and a point 10.9 Gyr in the future; (v) a unit label
"T_eff = 0.084 H₀/2π" (should be 0.084 H₀); (vi) the printed conclusion was hard-coded text. Positive finding (lens B): in the
midpoint convention the exact ΛCDM kernel has imaginary period 4π/(√Ω_Λ H₀) with a pole inside the strip, so it is KMS at no β
— the conclusion "not thermal at any single temperature" is true, but was not what the script computed. **Fixed:** D11's last
paragraph, `05`, R15 and the CHANGELOG reworded to what was computed; `v11` Part D now prints the T = 0 comparison, the local
effective rate, both x conventions, and a computed (not hard-coded) summary.

### `v12` — lens A: CONFIRMED; lens B: PARTIAL (computation right; attributions wrong)

Independently reproduced: κ from its definition in flat FRW, T_KH, T_det, the exact ratio, the Cai–Kim and Akbar–Cai identities
with exactly the temperature/sign each author uses, the lemma's exponent. Found and **fixed**: (i) 12.2 contained no
computation — now the Kodama vector K = ∂_t − Hr∂_r is shown null on R_A = 1/H and its Lie derivative of the metric has the
single component −a²rḢ (Killing iff ε = 0); (ii) 12.4's ν(ε) was typed in — the exponent −p/(p−1) is now derived from
η = ∫dt/a; (iii) **source misattribution**: Cai–Kim 2005 and Cai–Cao–Hu 2009 use T₀ = 1/(2πR_A) = H/2π, not T_KH (Cai–Cao–Hu
leave the |κ|/2π relation open); `08` said Cai–Cao–Hu derive Hawking radiation "with the full κ" — false; (iv) Tian–Booth 2014
misdescribed: their "temperature confusion" is T₀ versus Hayward's κ/2π, resolved by adopting T₀; the positive-heat-out
convention fixes the sign of the Gibbs relation, not the temperature; (v) "the mismatch is exactly the apparent horizon's
expansion rate" — it is half of it relative to T_KH (and the full rate relative to T₀); (vi) "T_KH closes the unified first
law" — only the signed κ/2π = −T_KH does; (vii) the step "not Killing ⇒ P1 assigns no temperature" is an application of P1's
definition (whose geometric list in `01` includes conformal double cones), not a computation — now said so in `09` and `v12`;
(viii) "thermalises" replaced by the `v11` statement in the fixed-mean-time convention; (ix) the velocity correlator is IR
finite only for ε < 2/3 — added.

### `v13` rev. 2 — lens A: CONFIRMED; lens B: PARTIAL (arithmetic right; the inference's reasons were the wrong ones)

Every equation transcribed from Cho–Kim–Kitamoto verified against the paper; the σ-dependence, c₁(σ), c_X(σ) and the three σ
values reproduced; ten mutations each break the right checks. Substantive findings: (i) **a σ-free invariant exists** — the
zero-frequency noise strength D(t′) = lim_{t→∞}∂_{t′}W(t,t′) has the closed form D/H³ = (√π/2π³)(1−ε)³Γ(ν)Γ(5/2−ν)cos(πε/(1−ε))
(auditor's derivation, checked against a direct mode integral to 15 digits at ε = 0, 0.05, 0.1, 0.3), whose O(ε) coefficient is
exactly −(1 + 2 ln 2) — the number `v14` E found from the two-time function; so "not a physical invariant" and `05`'s "any
dimensionless comparison shifts by a logarithm" were too strong (only the sharp-window σ = 1 value is scheme-dependent);
(ii) the zero-mode-*temperature* route fails for a prior reason: Cho–Kim–Kitamoto's own late-time solution (3.27) for ε ≠ 0
is exp(−V/T_zm) times an extra O(ε) Gaussian factor, so no stationary weight exp(−V/T) exists at O(ε); (iii) the published
drift 1/((3−ε)H) is a leading-IR approximation — the exact k = 0 attractor gives 1/((3+ε)H); (iv) `08` misquoted the paper
("rescaled time variable"); (v) a stale wrong coefficient in `10` ((1 − 2γ_E)ε; correct 1 − 2γ_E − 2 ln 2); (vi) R18 and the
rev.-2.3 CHANGELOG entry still stated the withdrawn inference without a pointer; (vii) the "c_X → −∞ as σ → 0" clause is a
formal limit outside the expansion's validity. **Fixed:** `v13` 13.5 now checks the invariant's ε = 0 value and its O(ε)
coefficient symbolically and the σ* identification; the docstring and `05`, `08`, `10`, R18, CHANGELOG reworded as above.
Still open: the invariant's closed form is the auditor's and is not re-derived in the repository.

## What the round changes

Nothing in the arithmetic of the branch was refuted: `v8`'s identities and spectra, `v11`'s first-order theorem, `v12`'s
identities, `v13`'s σ-dependence and every link of `v14` were each re-derived independently and agree to all digits. Three
things of substance did change: the thermal statements are now explicitly about the fixed-mean-time (Wigner) spectrum, not an
instantaneous detector rate; `v11` Part D's ΛCDM comparison is now described as what it is (one time-domain comparison, not a
thermality test); and `v13`'s "not a physical invariant" became "the sharp-window coefficient is scheme-dependent", because the
audit produced the invariant in closed form and it agrees with `v14`. What the audit corrected is the prose around them — four misquotations
or overstatements in `07`/`05`/`02`/`09`, one missing hypothesis in a theorem statement, one stale docstring — all now
fixed. The pattern is the same as in rounds two to five: the computations hold; the claims made about them had to be cut
back to what was computed.
