# m1_ledger — adversarial skeptic pass (key `m1_ledger_skeptic`)

1 November 2026. Claude, issue #13, brief [BRIEF_M1.md](../../BRIEF_M1.md). Target: the thirteen claims and the new item of
[NOTES.md](NOTES.md) (track `m1_ledger`). Codex packets read only through `git show` (d4747cb, 1fcc6b7, d71cc01, ed73791,
bf8ab3b, c481530, e160356, 8986d34, b5174e5). **Own code only, algebra and quadrature, no engine runs, no Codex code run**:
`skeptic/code/skeptic_algebra.py` (mpmath, 50 digits, every number rebuilt from the raw packet inputs), `skeptic/code/skeptic_born.py`
(the Born pair integral with a third numerical method), `skeptic/code/make_summary.py`. Results: `skeptic/results/*.json`.
One thread, foreground, about 4 s. Nothing here validates, excludes or predicts; same-session AI work on one supplied model.

**Outcome: no claim of m1_ledger is refuted.** Every key number reproduces by a different route; four small corrections
and seven new issues are listed, the first two of which bear on how much the "+0.7009%" calibration means.

## Methods that differ from m1_ledger's

- Frozen test, old target, λ chain, A, Γt₀: rebuilt from `CANDIDATE_M1.json`, `DEFORMATION_RESULTS.json`, `ABUNDANCE_LEDGER.json`,
  `FIXED_YIELD_RESULTS.json`, `COUPLING_RELATIONS.json`, `COLD_LOSS_RESULTS.json`, `PORTAL_RESULTS.json` in 50-digit arithmetic.
- Yield convention without any engine: Codex's rows (N_X, a, T_vis, H) satisfy Y = |N_X| H³/(a³ s), s = 2π²g*T³/45, g* = 106.75,
  to 4×10⁻¹⁶ on all seven rows (`codex_yield_formula_from_rows`), so the "+0.7009% by actual ODE rerun" is the 2/3 root of the
  baseline yield converted with the scaled H_i. The exponent 3/2 is taken from Codex's own **scale_0.5 / scale_2.0** control rows
  (a factor 4 in H_i, which m1_ledger did not use), not from our engine.
- A from the M1 heavy parameters through hw_gates' α₆ = κ²μ⟨Σ⟩/(M_Y⁴ − 4μ²⟨Σ⟩²), ⟨Σ⟩ = h/√2 (not only the λ⁻² scaling).
- Born integral: the profile's **native grid** (no spline; M1 dr = 0.05 to 160, old dr = 0.075 to 120), Simpson transform,
  adaptive `quad` with the algebraic endpoint weight (q_max − q)^{1/2}; a second spline/trapezoid variant; the Parseval identity
  ∫q²|F̃₄|²dq = 8π³∫r²F⁸dr as an absolute normalisation check (closes to 2×10⁻¹⁶); a Gaussian control of the transform (1×10⁻¹⁶).
  The prefactor 900/(32π³) was re-derived from scratch (C(6,2) = 15, contraction 2, two-body factor (π/2)β, identical-particle ½).

## Verdict table (my number | m1_ledger | Codex)

| # | Claim | Label | Mine (this pass) | m1_ledger | Codex | Verdict |
|---|---|---|---|---|---|---|
| 1 | frozen coefficients, endpoints | [D] | 0.0667804674, 2.133560935; η ≥ 0.8502556152, p ≤ 0.0046869999526 | 0.06678047, 2.13356093; 0.8502556, 0.0046870 | 0.0667805, 2.13356; 0.85026, 0.004687 | stands (2e-15) |
| 2 | Y_Q = D/ε_Q | [D] | 5.079414232725842e-9 | 5.079414233e-9 | 5.07941423e-9 | stands |
| 3 | old target 4.8294e-9 ↔ 0.936494 D | [D] | D/Y_Q,old = 0.0911086092 GeV = 0.985928 m_X = d71cc01 late mass/charge (identity 3e-16) = **object 0.923316 + 2×0.031306 released pairs**; 0.9364942; ×1.0678122; kin f 2.028545 | 0.98592 m_X "v-calibrated supplied"; same numbers | 0.936494; 1.0678122 | stands, origin sharpened |
| 4 | +0.7009% from the old point | [C] | λ = 1.0070090495735 (2/3 root of Codex's baseline row); H = 1.342842735e11 GeV | 1.0070090496 (own engine) | 1.00700905 | stands; **clock band 1.00484–1.00834** (N2), D band ±0.67% (N1) |
| 5 | old-object repairs; (5,1) M1 | [C] | 1.0172235041614; **1.0162835956** (m1_ledger: 1.016277, slip); 1.006078579112 | 1.0172235042, 1.016277, 1.0060785791 | 1.0172235042, 1.0162835956, 1.0060785791 | stands; one digit slip (N4) |
| 6 | Y ∝ H_i^{3/2} | [D] | Codex scale_0.5/2.0 rows: exponent 1.5000000003 (Y), 0.4999999996 (T_vis) over ×4; N_X, a unchanged to 5e-10 | 1.5000000002 (own engine, λ ≤ 1.05) | "verified by ODE rerun" | stands (N5) |
| 7 | baseline yields | [C] | from Codex row numbers: 5.026475516056737e-9, 5.03345022205823e-9 (4e-16) | 5.02647551601e-9, 5.03345022182e-9 (engine) | 5.026475516e-9, 5.033450222e-9 | stands (engine not re-run) |
| 8 | our 1.0454 vs Codex 1.0408 | [C] | `rp_reheating/results/extra.json` yr_b3_g0.5: end 80, window N_X −8156.7346; ratio 8156.7346/8120.7401 = 1.0044324 | 1.0454 (τ 80) → 1.0408 (τ 800); 1.00443 | 1.040807 | stands; Codex's instantaneous clocks give 1.00382 (different estimator) |
| 9 | populations vs the test | [E,A] | pulse rows with η = 1 + (N̄−N₊)/Q₀: 0.00432/0.01749/0.07124/0.3015; kin rows and crossing 5.824 → 5.768 reproduced | 0.00445/0.01802/0.07354/0.3113 | η, p "remain to compute" | stands; pulse rows 3% high (N3); excitation energy ignored (N6) |
| 10 | cold completions 0.0275/0.1127/0.4903/2.482% | [D] | (E(Q_f) + N_free m_X)/M_ref − 1 = 2.751558e-4/1.127453e-3/4.903148e-3/2.481765e-2 (≤ 2e-13); frozen/actual 0.529–0.551 | same; ×1.9 | same | stands |
| 11 | A → 0.98613 A | [D] | 2.079605512546138e-31 from κ, μ, h, M_Y (1e-16); A_old,full/λ² = Codex (4e-17); r² = 0.2025 | 2.0796057e-31 (from the 6-digit A_old) | 2.0796055e-31 | stands |
| 12 | Born pair channel | [C] | I_M1 = 2868466.5179 (native grid, +2.4e-9), 2868466.5114 (spline/trapz, +1.5e-10); I_old = 3149626.631 (+1.3e-8); c₁ rates 9.47134e-47, 1.03997e-46; Parseval 2e-16 | 2868466.51, 3149626.59 | 2868466.5109 | stands |
| 13 | Γt₀ 6.27e-17 (from 7.08e-17) | [C] | new 6.2716294e-17 (+2.4e-9); old 7.0836787e-17 = e160356's 7.08367857848652e-17; ratio 0.885363 = 0.9104478 × 0.9724483 | 6.27163e-17; old 7.0837e-17 (6-digit A, 2e-7 off) | 6.2716e-17; 7.08368e-17 | stands |

Units checked: ε_Q, D, m_X in GeV; A in GeV⁻²; c6 = A M_red² dimensionless (1250390.42 old, 1233044.90 new); Born rate (900/32π³)A²(v/√2)⁸I/m_X³
in GeV, /ħ (6.582119569e-25 GeV s) → s⁻¹; fractional core rate 4 events/Q (core loses 4, −2 escapes); age 13.8 Gyr Julian = 4.3549488e17 s.
Conventions checked: the three-reservoir ledger (free X = (1 − η + p)Q, X̄ = pQ, charge constraint 1 − η + p ≥ 0); Codex's late window τ ≥ 0.8·end
(trapezoid), the same as our engine's; η² = (λH₀/M_red)² scaled in Codex's `match.py`; r = 2|A|/B = 0.45 invariant.

## New issues (none refutes a claim; N1–N2 bear on the meaning of the calibration)

- **N1 [E, R] The target D = 4.4e-10 GeV is a rounded value.** Ω_c h² = 0.1200 ± 0.0012 (Planck 2018 VI), ρ_c/h² = 1.05368e-5 GeV cm⁻³
  and s₀ = 2891.2 cm⁻³ (PDG; [R] recalled, not fetched) give D = 4.3733e-10 GeV: Codex's value is 0.61% high with a ±1% band. At the
  Planck value λ moves by −0.40%, and the band is ±0.67%. **The +0.7009% calibration step is inside the target's own rounding and
  uncertainty.** Neither Codex nor m1_ledger states this; it does not touch the 1% program benchmark, which Codex declares not observational.
- **N2 [C] λ is quoted to 1e-8 but is clock-dependent at the 10⁻³ level.** Codex's own `charge_clock` entries for the baseline
  (instantaneous N_X at τ = 20/50/80/200/800: −8117.9/−8104.6/−8147.1/−8133.9/−8116.1) give λ = 1.00484 (τ 80) to 1.00834 (τ 50);
  the window average gives 1.00701. m1_ledger reported the 0.44% yield drift; in λ it is −0.22/+0.13% absolute, a third of the step.
- **N3 [D] Pulse rows are 3% high.** Reading emitted X̄ as p with η = 1 ignores that the core charge rises (η = 1 + (N̄ − N₊)/Q₀, the
  ed73791 remark m1_ledger itself confirmed); the frozen term is (c_p − c_η)N̄/Q₀ + c_η N₊/Q₀ = 2.0668 N̄/Q₀ + 0.0668 N₊/Q₀:
  10% pulse 0.01749 (×1.75), not 0.01802 (×1.8). Verdicts unchanged.
- **N4 [D]** m1_ledger's "(5,1): 1.016277" is 1.0162836 (Codex row `b5_g1_matched`); +1.628% either way; superseded number.
- **N5 [C]** Codex's scale_0.5 / scale_2.0 rows already fix the exponent over a factor 4 (1.50000000026, T 0.4999999996, N_X and a
  unchanged to 5e-10); the own-engine measurement in m1_ledger confirms but adds nothing, and the "ODE rerun" cannot fail by construction.
- **N6 [A]** The frozen formula assigns ε_Q to a pulsed core; c481530 reports 2.305% excess energy in the full 10% pulse and a 0.2237%
  unprocessed remainder, so the pulse-row readings (not only the kinetics row) are lower bounds on the mass excess.
- **N7 [A]** T_vis(τ = 800) = 7.7e12 GeV exceeds M_Y = 4.03e12 GeV; g* = 106.75 is the SM count and holds only because the conservative
  candidate has zero SM portals; the "separately tested portal extension" would need g*, the entropy and the yield recomputed.

## What this pass does not do

No engine run (so claim 7's 1e-11 engine agreement is accepted as a same-formula two-engine comparison, which is what it claims); the
coherent channel (6.4e-5 of the included rate) is Codex's number; the population rows of claim 9 are translations of other tracks'
numbers, checked for arithmetic and bookkeeping only; nothing about whether η or p can be computed is decided here.

Reproduce: `cd skeptic/code; OMP_NUM_THREADS=1 python3 -B skeptic_algebra.py; python3 -B skeptic_born.py; python3 -B make_summary.py`
(profiles are the `git show` copies, sha256 c08176e9… and 7a737d67…, kept in the scratchpad under review/m1_ledger_skeptic/data/).
