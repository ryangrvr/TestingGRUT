# GRUT III — changelog

## Revision 2.1 (2026-10-11, after the second round of reviews)

- `v11` rev. 3: revision-2 headline (first-order c ≈ 1.01 plus 1%-of-ε frequency dependence) **retracted** — it was
  O(ε²) contamination. Now proved: first-order exact thermality at T_eff = (ä/ȧ)/2π for any a(t); second-order
  coefficient f₂(x) = x²/4 − (x³/12)coth(x/2) (reviewers') confirmed; Richardson D/ε → 2πω/H̄ (c = 1 to 3e-5); the true
  O(ε²) frequency-dependent defect D₂(ω); present-day (ε₀ = 0.47) exact ΛCDM deviation 0.1–19% from any thermal form.
- T5/T6 scopes tightened (asymptotic Brownian covariance; representation ≠ physical field theory).
- `08`: the independent review (branch `codex/grut-causal-state-review`) recorded; RESPONSE R15.

## Revision 2 (2026-10-11, after three external reviews and an owner's note)

Corrections (all accepted; see `RESPONSE_TO_REVIEWS.md`):
- η = c³/16πG relabelled as the horizon's surface viscosity, units kg s⁻¹ (was Pa·s).
- "Memory time 1/(4H)" relabelled as the asymptotic decay rate of the pointwise ⟨T̂T̂⟩; memory now defined only for
  smeared observables (`07`, `v7`).
- a₀ Jacobson-route argument retracted; `06` is a coincidence appendix crediting Milgrom (1999).
- w ≡ −1 reframed as a falsifier of the minimal framework shared with ΛCDM; "DESI DR3" removed; DR2 figures cited.
- P1 restated in modular (Tomita–Takesaki) form; no factorisation or influence functional claimed from it.
- Ledger relabelled: 5 postulates, imports tagged IMPORTED (zero credit), data; "compression" withdrawn.
- P5 derivation via HJW withdrawn; reconstruction literature cited.
- v4 wording: √s vs frequency; staircase marked illustrative; on-shell s = 0.
- Priority: dS/conformally-flat noise kernels credited to Pérez-Nadal–Roura–Verdaguer and Cho–Hu; "not on the record"
  downgraded to "corollary; placement".
- Record correction D1 now proved (the record's S(w) is G⁺(w) by its own definition); gap quoted as 0.55%.

Additions:
- `07_ACCESS_RESOLVED_RESPONSE.md`: the Δρ residual framework (FDT identity, trichotomy, memory dictionary).
- `calc/v7`: memory from access geometry (reviewer's smeared-field formula; smeared stress-tensor kernel closed form;
  KMS for smeared observables).
- `calc/v8`: residual framework; the record's kernels as Δρ ≠ 0 by fiat.
- `calc/v9`: induced-gravity / entanglement-entropy matching per species; SM species sums (induced 1/G tiny or negative).
- `calc/v10`: gravitational sector = stationary thermal process + Brownian zero mode, D = H³/4π².
- `calc/v11`: conventional O(ε) KMS defect under slow roll (the first non-KMS baseline).
- `calc/requirements.txt`; `run_all.py` extended to v1–v11.
- `08_PRIORITY_AND_SOURCES.md`: nine papers checked at equation level (PNRV is ξ = 0 only; conformal kernel is Cho–Hu's;
  Milgrom's coefficients are definite and overshoot; HLZ widths contain D5's laws; Jacobson confirms R1; NPT/Deser–Levin
  confirm the dS Unruh temperature); Stage-3 charter located (DEF-16, SCOPE-Q).
- Verification status: the adversarial three-lens pass on v7–v10 did not run (usage limit); author-level only; owed.

## Revision 1 (2026-10-11, commit 862e666)

Initial principle-first reconstruction: P0–P5, D1–D8, v1–v6.
