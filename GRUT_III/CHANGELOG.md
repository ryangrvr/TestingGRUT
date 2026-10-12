# GRUT III — changelog

## Revision 2.4 (2026-10-11, sixth round) — final; baseline banked 2026-10-12

- `12_BASELINE_AND_HANDOFF.md`: verification rounds stopped on the reviewer's recommendation; the honest ledger of what is
  new (two or three candidate conventional results) versus rediscovery; the seven constraints a Stage-3 Card 1 must meet,
  derived from the baseline; the three-way candidate space with each route priced; the handoff. README status → BANKED
  BASELINE; `05` O7 closed as a list.
- `v13` rev. 2: the reviewer's coarse-graining-scale objection accepted; the σ-dependence A₀(σ) = A₀(1)σ^{−2ε/(1−ε)} derived
  from the mode function; c₁(σ) = 1 − 2γ_E − 2 ln(2σ), c_X(σ) = 0.4685 + ½ ln σ; rev. 1's inference ("no horizon temperature
  matches") **withdrawn** — the O(ε) zero-mode coefficient is not a physical invariant. Items 1 and 3 still close on P1 and
  conventional QFT. `05` O5, `10` §4, `08`, RESPONSE R19.
- `05` O5: option (b) of the fork identified with a Stage-3 𝒦 card (Card 1 under the frozen charter on `grut2-stage3`).
- Sixth-round verification workflow over `v8`, `v11` (general a(t), ΛCDM), `v12`, `v13` rev. 2: `11_VERIFICATION_ROUND.md`.
- `calc/v14_massless_sector_first_order.py` (PASS) and D13: the σ-free computation. Closed form for the massless minimally
  coupled comoving two-point function in power-law FRW (dS function of index 3/2 + ε/(1−ε) × scale-factor ratio; checked
  against a mode integral); first-order theorem W = (1 + 3ε/2)W_mmc(Δ; ä/ȧ) + ε(H̄²/8π²)S₂(cosh²(H̄Δ/2)) + ω = 0 terms
  (Richardson); midline KMS symmetry ⇒ **thermal at (ä/ȧ)/2π at first order, Kodama–Hayward excluded, no σ**; spectral
  shape change r(ω) (r(0) = 2 − 2 ln 2, r(∞) = 0); σ-free diffusion coefficient D = (H̄³/4π²)[1 − (1 + 2 ln 2)ε]
  (= stochastic value at σ* = e^{1−γ_E}). `05` O5 item 1 computed; `09` §4; `08`; RESPONSE R20; `run_all` runs v1–v14.
- `v14` audited (three lenses, `11`): mathematics and numerics confirmed; framing corrected — the thermal statement is about
  the fixed-mean-time (Wigner) spectrum at ω ≠ 0, not an instantaneous one-sided detector rate (which differs at O(ε) in both
  sectors and is non-thermal for the massless sector even in exact dS unless derivative-coupled); general-a(t) and graviton
  extensions marked not computed; D's Δ² companion term stated; sign slip fixed; exact-function KMS test and r asymptote
  check added to the script. D11, D13, `05`, `09` reworded.
- Round-6 audit of `v8`, `v11`, `v12`, `v13` (ten agents, `11`): no arithmetic refuted; fixes — `v11` Richardson weights (c = 1 now
  to 1e-9), symbolic first-order identity, εη closed-form check, Part D reworded and recomputed (T = 0 comparison, local rate,
  both x conventions); `v12` 12.2 now a computation (Kodama vector), 12.4 exponent derived, attributions corrected (Cai–Kim and
  Cai–Cao–Hu use T₀ = 1/2πR_A; Tian–Booth adopt it), "half the expansion rate"; `v13` 13.5 adds the σ-free invariant
  D/H³ = (√π/2π³)(1−ε)³Γ(ν)Γ(5/2−ν)cos(πε/(1−ε)) with O(ε) coefficient −(1 + 2 ln 2) = `v14` E, and the no-stationary-weight
  point; D11, `05`, `08`, `09`, `10`, R15, R18 reworded accordingly.

## Revision 2.3 (2026-10-11, fourth round)

- `10_SECULAR_IR_SURVEY.md`: literature check for O5 item 3 (five agents; 30+ papers, most at full text). Scalars: all
  resummations agree at leading log and NLO — closed. Gravitons: live dispute over which observables carry secular logs
  (Tsamis–Woodard school vs Garriga–Tanaka/Tanaka–Urakawa/Fröb relational observables); a quantum-gravity question, not a
  GRUT opening. Item 1 largely published (Abbott–Wise, Lucchin–Matarrese, Garbrecht–Prokopec, Cho–Kim–Kitamoto 2015 [arXiv:1508.07877]; the 2018 derivative-interaction paper 1811.01830 is Kitamoto alone).
- `09` §4 and `05` O5: lemma scope (conformal-time mode equations only), item-1 priority, the decisive `v13` computation.
- `calc/v13_zero_mode_temperature.py` (PASS): the published O(ε) Fokker–Planck coefficients of the massless zero mode in
  power-law FRW (Cho–Kim–Kitamoto 2015, eqs. 3.14–3.16; confirmed at equation level) define a zero-mode temperature
  T_zm = (3−ε)A₀H⁴/2 = 3H⁴/8π² + O(ε) whose O(ε) coefficient (2/3 − 2γ_E − 2 ln 2) is the Euclidean weight of **no**
  horizon temperature (detector c = 1, Kodama–Hayward c = 1/2, naive c = 0; required c_X = 0.4685). P1 assigns none and
  adds nothing; conventional theory fixes the coefficients. **O5 items 1 and 3 close.** `05` O5, `08`, RESPONSE R18.
  [Inference withdrawn in revision 2.4 (σ-dependence; no stationary weight at O(ε)); entry kept as the record.]
- `run_all.py` now runs v1–v13; `08` gains the Cho–Kim–Kitamoto source entry and the survey's source list pointer;
  verification status updated (v12, v13 author-level).

## Revision 2.2 (2026-10-11, third round)

- `09_TWO_TEMPERATURES.md` + `calc/v12`: the reviewer's O(ε) question — detector temperature (ä/ȧ)/2π vs
  apparent-horizon Kodama–Hayward temperature (H/2π)(1 − ε/2) — answered from P1: only the state's KMS parameter is a
  bath temperature; the apparent horizon is non-Killing for ε ≠ 0; Cai–Kim and Akbar–Cai horizon first laws verified
  to be identities of the Einstein equations (symbolic); no heat-flow law; no O(ε) non-thermal window for the conformal
  sector. Lemma for the massless/graviton sector: power-law massless mode equation ≡ dS massive with ν = 3/2 + ε/(1−ε).

## Revision 2.1 (2026-10-11, after the second round of reviews)

- `v11` rev. 3: revision-2 headline (first-order c ≈ 1.01 plus 1%-of-ε frequency dependence) **retracted** — it was
  O(ε²) contamination. Now proved: first-order exact thermality at T_eff = (ä/ȧ)/2π for any a(t); second-order
  coefficient f₂(x) = x²/4 − (x³/12)coth(x/2) (reviewers') confirmed; Richardson D/ε → 2πω/H̄ (c = 1 to 3e-5); the true
  O(ε²) frequency-dependent defect D₂(ω); present-day (ε₀ = 0.47) exact ΛCDM deviation 0.1–19% from the first-order thermal form at H_eff (time domain; what this
  does and does not show is in `11`).
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
