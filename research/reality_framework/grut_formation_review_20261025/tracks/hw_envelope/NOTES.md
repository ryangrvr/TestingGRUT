# Track hw_envelope: Codex 2908e91's envelope bound against the ls_scaling heavy window

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication.

Packet: 2908e91 (`origin/codex/grut-concept-checkpoint-20261025`, `research/reality_framework/grut_concept_checkpoint_20261025/`),
read only through `git archive` into the scratchpad. No Codex code was run; Codex's CHECKS.json numbers are quoted [Q]
for comparison only. Our own ls_scaling results (`tracks/ls_scaling/results/scan_B.json`) are quoted [Q, own track].
Nothing was committed.

**Scope.** One supplied tree-level working model; homogeneous classical release, instant at a = 1 with instant radiation;
supplied CP phase; the late-particle decay ledger and the case-L proxy with the 1% allowance, both quoted from earlier
tracks. Nothing here is evidence that X, Y, Σ or the vector exist. It is not a validation, an exclusion or a viable point.

Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued.

## Question

Codex 2908e91 (i) shows algebraically that a nonzero charge source does not by itself imply residual heavy quanta, and
(ii) bounds the charge of the quiet (forced-response, EFT) regime at the fixed published couplings by
|a³n_X| ≤ 12|α6|R⁶/H_i, R = r_i(1+2/c_H)^{1/4}, giving ≤ 1.74e-4 Y_Q even at r_i = M_Y. It asks whether producing
enough lasting charge necessarily generates excessive processed late debris. Our ls_scaling track claims a heavy-mediator
window (M_Y/H_i ≈ 13–43) where the spot point (H = 1.334e11, M_Y = 4.00e12, h = 5.33e12, μ = 9.55e11 GeV, κ = 3.46e-3,
r = 0.45, Y-relaxed preparation) gives ε = −0.146, yield 1.00 Y_Q, R = 1.14 against R_L = 88. Are the two consistent?

## Answer

**Yes. The spot point satisfies Codex's bound (margin 2.65 with Codex's leading R, 4.07 with R built from the actual
prepared amplitude), sits inside Codex's own matching regime, and the two statements evaluate the same bound on two
different lines in coupling space.** Codex moves H_i at fixed (M_Y, h, μ, κ): along that line the bound never reaches
Y_Q while the matching diagnostics hold. ls_scaling moves (κ, r, M_Y/H_i, h/H_i) at fixed η = α6X0²/λ_X: along that
line α6 in H_i units is 528× larger and H_i is 9.4× larger, which is the whole factor 15,279 between 1.74e-4 Y_Q and the
bound's 2.65 Y_Q at the spot point. Within ls_scaling's declared domain the implication "enough charge ⇒ excessive late
debris" fails: our own release gives 1.00 Y_Q with 0.35 heavy quanta per comoving H_i³ (1.3e-4 per unit charge) and
R = 1.14 against R_L = 88 [Q].

## 1. Codex's four control groups, reproduced by other routes [D, C]

All in `code/controls.py` (2 s). No Codex code imported.

| Group | Codex | Ours | Verdict |
|---|---|---|---|
| Exact sources (−6 Im K, 2 Im K − 4 Im L, 2 Im L); gauge sum 0; n_D sourced by −12 Im L | symbolic, phase coordinates | symbolic in complex variables with arbitrary radial g(\|z\|²) added (they drop out); plus a dynamical check: a random three-field release with H = 1/(2t) integrated, d(a³n_i)/dt against the quadrature of the source: rel. error 6e-7 (trapezoid), gauge drift 1e-11 | CONFIRMED |
| 48 complex samples: force source, angular finite difference, gauge sum | 1.5e-14, 7.9e-11, 8e-14 | different seed, plus a numerical dV/dz* (3e-11): 1.9e-14, 1.4e-10, 1.1e-13 | CONFIRMED |
| Forced response Y_* = (−AJ + BJ*)/D; ΔV = −A\|κ\|²\|X\|⁶/D + [μκ²Σ*X⁶/D + h.c.]; zero displacement action; four amplitudes | numeric | sympy proof of Y_*, of V_Y(Y_*) = ΔV, and of −6 Im(κY_**X³) = −12 Im(α6X⁶) = −12 Im(μΣ*Y_*²) identically; numerics by BFGS minimum and fsolve with λ_Y = 1: Codex's four source values to all digits (dev 0) | CONFIRMED; **added [D]**: on the forced branch the n_Y source 2 Im K − 4 Im L vanishes identically, so the forced Y transfers charge without accumulating any |
| ρ_phase = n_X²/(2f²) + n_X²/(72h²), n_X = 6Iδ̇, n_Σ = −n_X/6 | derivation | sympy elimination of A_0 from the gauged phase Lagrangian | CONFIRMED |
| Yield bound table 1.2201e-8 / 5.4890e-8 / 3.8583e-5 / 1.7358e-4 Y_Q; ∝ r_i^{7/2}; extrapolation r_i/M_Y = 18.24 / 11.87 with Y-quartic 477 / 9.16 and Σ-shift 78,344 / 1505 | computed | own evaluation, same prompt-radiation entropy (identical to our scal.py convention): all values to the printed digits (dev < 1e-12); α6 = 3.5533e-32 GeV⁻² = rκ²/(2M_Y²(1−r²)) | CONFIRMED |

## 2. The bound at the ls_scaling spot point and ladder [C]

`code/spot_bound.py` evaluates the closed form on our scan_B rows [Q, own track]; `code/release_spot.py` is an
independent integrator (Python DOP853, cosmic time, physical fields, own multi-start preparation) that re-runs the
release and tests each premise on the solution.

**Own release at the spot point (m = M_Y/H_i = 30, r = 0.45, h = 40H_i, η = 0.5 η_edge, δ = π/2, Y-relaxed, a_end = 30;
H_i = 1.3335e11 GeV for the yield):**

| Quantity | ls_scaling | Ours (own integrator) |
|---|---:|---:|
| Prepared \|X_i\|, Y, \|Σ\| | 76.3998, (−2.11151, −0.93933), 28.29314 | 76.39978, (−2.11151, −0.93933), 28.29314 |
| ε (late window) | −0.14572 (a = 30–100) | −0.14576 (a = 22–30); −0.14567 (a = 75–100 control) |
| a³n_X | −2788.18 | −2790.0 (a = 22–30); −2788.19 (a = 75–100 control) |
| E_X | 38,386 | 38,407 |
| N_YH, N_YL, N_σ | 0.01725, 0.3345, 3.2e-9 | 0.017253, 0.33447, 3.2e-9 |
| Yield / Y_Q | 1.000 | 1.0007 |
| R_Y, R_casc, R_central (closed-form ledger [D]) | 0.806, 0.339, 1.144 | 0.805, 0.339, 1.144 |
| R_L (case L, normalised yield) | 88 [Q, sp_critic] | not rechecked |

**Codex's bound there.** α6 = rκ²/(2M_Y²(1−r²)) = 3.750e-9 H_i⁻² = 2.11e-31 GeV⁻² (c6 = 1.25e6), η = α6X0²/λ_X = 0.01875,
same entropy convention as Codex (own yield from the scan's a³n_X reproduces every ls_scaling yield to 4 digits).

| R used | R (H_i) | bound on \|a³n_X\| | simulated | margin | bound yield / Y_Q | ε_max |
|---|---:|---:|---:|---:|---:|---:|
| leading relation r_i = √(c_H/2λ_X) H_i = 70.71 → R = 74.01 | 74.01 | 7394 | 2790 | **2.65** | 2.65 | 0.465 |
| actual prepared r_i = 76.40 → R = (r_i⁴ + r_i²/λ_X)^{1/4} = 79.48 | 79.48 | 11,344 | 2790 | **4.07** | 4.07 | — |

**The spot point satisfies the bound.** Premise by premise on the full solution:

| Premise | Status at the spot point |
|---|---|
| Leading radial relation H_i = √(2λ_X/c_H) r_i | Off by 8.0%: the prepared \|X_i\| is 1.0805 × √(c_H/2λ_X) H_i because the companion β6 = −M_Y²κ²/D deepens the Hubble-held minimum. Same 5–8% along the whole ladder. |
| Envelope \|X\| ≤ R/a | With R_lead: violated by 8.0% in amplitude (max a\|X\| = 79.97 at a = 1.09; late envelope from E_X 78.7, i.e. R⁶ is low by 1.59). With R from the actual r_i: holds to 0.6% (79.97 vs 79.48; late 78.7 < 79.48). |
| Fixed broken Σ | Holds: Codex's Σ forced-shift diagnostic 6.4e-4, N_σ = 3e-9, \|Σ\| moves by 3e-4 of h/√2 at preparation. |
| Zero initial physical velocity | Holds by construction (both tracks). |
| EFT matching (forced response) | Holds: the α6 source −12 Im(α6X⁶) integrated along the actual trajectory gives 0.980 of the actual charge; the one-field matched EFT release gives ε = −0.1387 vs −0.1458 full (5% low); Codex's diagnostics κR²/M² = 0.021, ω_X/gap = 0.149, ω_X/m_σ = 0.059, Y-quartic force 0.032, λ_Y\|Y\|²/M² = 0.006; \|Y − Y_*\| = 0.047 at preparation against \|Y_*\| = 2.36. |
| Phase coherence (the triangle inequality) | Nearly saturated: with the actual \|X(t)\| the integral 12∫a³\|α6\|\|X\|⁶dt is 3841 and the charge 2783 (0.72; 0.985 by a = 3). The slack of the R-bound is the oscillation average ⟨\|X\|⁶⟩ ≈ 0.4 R⁶ of a nearly radial orbit, not phase cancellation. |

So the bound is a factor ≈ 2.6 estimator of the real charge in this regime, and ε ≈ 0.38 ε_max ≈ 7.8 η (ls_scaling: 7–10 η).

**Ladder (own closed form on scan_B rows; m = 10, 30, 40 also re-run with our integrator):**

| m | κ | α6 (H_i⁻²) | bound (R_lead) | sim \|a³n_X\| | margin | sim ε | own ε | own N_Y | own R_central | R_L [Q] | Codex Y-quartic diag. | EFT/actual charge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 1.153e-3 | 3.75e-9 | 7394 | 1998 | 3.70 | −0.1130 | −0.1126 | 22.42 | 144.5 (fails) | 60 | 2.63 | 1.24 |
| 15 | 1.729e-3 | 3.75e-9 | 7394 | 2468 | 3.00 | −0.1337 | — | — | — | 77 | 0.52 | — |
| 20 | 2.306e-3 | 3.75e-9 | 7394 | 2687 | 2.75 | −0.1424 | — | — | — | 85 | 0.16 | — |
| 30 | 3.459e-3 | 3.75e-9 | 7394 | 2788 | 2.65 | −0.1457 | −0.1458 | 0.352 | 1.144 | 88 | 0.032 | 0.98 |
| 40 | 4.612e-3 | 3.75e-9 | 7394 | 2791 | 2.65 | −0.1455 | −0.1455 | 0.0888 | 0.251 | 88 | 0.010 | 0.98 |

The bound is the same along the ladder (η fixed ⇒ α6 in H_i units fixed), and every row satisfies it (margins 1.8–3.7
over all 13 scan_B rows, 3.1–5.1 with the actual r_i). The lower edge of the window coincides with the end of Codex's
matching regime: at m = 10 the Y-quartic diagnostic is 2.6, the EFT overestimates the charge by 24%, the prepared
\|Y − Y_*\| is 2.4, and N_Y = 22.

**Preparation control (Codex prescription at the spot couplings, own run):** ε = −0.1588, N_YH = 76.8, N_σ = 0.039,
R_central = 1415 (ls_scaling 1414.9). The prepared free mode is \|Y − Y_*\| = 1.47 (relaxed: 0.047). This is the
dynamical instance of Codex's point (i): at a heavy point the residual heavy quanta are the prepared displacement from
the forced response, not the charge source. Codex's excitation form a³ρ_exc/M_Y is 0.27 (relaxed) against 93 (Codex
preparation) for the same charge.

**Why the fixed-coupling case is 15,279× lower (bound 1.74e-4 Y_Q at r_i = M_Y, c_H = 10 versus 2.65 Y_Q):**

| Factor | Value | Origin |
|---|---:|---|
| α6 in H_i units | ×527.7 | r/(1−r²): ×7.94 (0.0707 → 0.45); κ²: ×11.96 (1e-3 → 3.46e-3); (M_Y/H_i)⁻²: ×5.56 (70.7 → 30) |
| R⁶ | ×1.000 | same c_H = 10 and λ_X ⇒ same r_i/H_i = 70.7 |
| H_i^{3/2} (a³n_X H_i³ over s ∝ H_i^{3/2}) | ×28.95 | 1.414e10 → 1.3335e11 GeV |
| Entropy convention | ×1 | identical (prompt radiation, g* = 106.75, M_red = 2.435e18) |
| Product | 15,279 | equals the ratio of the two bounds to 1e-12 |

GeV view of the same ratio: α6 ×5.93 (κ² ×12 · r-factor ×7.94 / M_Y² ×16), r_i/M_Y = 2.36 instead of 1 (×20.1 via
r_i^{7/2}), overall scale M_Y ×4 at fixed r_i/M_Y (×128 via M^{7/2}).

**Fixed couplings along m (Codex's own line, c_H = 10):** at M_Y/H_i = 70.7 / 43 / 30 / 13 / 5.96 / 2.30 the bound gives
ε_max = 0.0009 / 0.0024 / 0.0049 / 0.026 / 0.124 / 0.83 and yield ≤ 1.7e-4 / 9.9e-4 / 3.5e-3 / 0.065 / 1.0 / 28 Y_Q.
The matching diagnostics pass only down to m ≈ 30–43 (at fixed h = 1e11 GeV, Σ goes light: ω_X/m_σ = 1.8 and the
Σ-shift 13.9 already at m = 13). So at the fixed couplings the quiet regime gives ε ≲ 0.005 and ≲ 3.5e-3 Y_Q, and Y_Q is
reached only at m ≈ 6 (Y-quartic 9.2, Σ-shift 1505), on the way to the resonant stated candidates at m = 2.3. **Codex's
statement is confirmed and sharpened.**

## 3. Relation between the two statements [A, on the numbers above]

- Codex: "the quiet regime at the fixed couplings is unproductive" = the bound's value on the line
  (M_Y, h, μ, κ) fixed, H_i varied. True; the bound and the matching diagnostics both say so, and ls_scaling's own
  one-at-a-time scan at the stated point agrees (no single move passes).
- ls_scaling: "a heavy window exists after a coordinated coupling change" = the same bound (and the full release) on the
  line η = α6X0²/λ_X = 0.5 r/12 fixed, with r 0.0707 → 0.45, κ ×3.46, M_Y/H_i 70.7 → 30 (and h/H_i raised to 40 so Σ
  stays heavy). There the bound allows 2.65–4.07 Y_Q, the full release gives 1.00 Y_Q, every Codex diagnostic is < 0.05,
  and the processed debris is 1.3% of R_L.
- They are consistent; they differ only in domain. Codex's packet already states that its failure "does not extend to
  changed couplings". The window is exactly the case Codex's point (i) anticipates: the forced Y carries the X⁶ torque,
  the residual free Y is the sudden-release kick (N_Y ∝ m^{−4.6}, 0.35 at m = 30) plus whatever the preparation puts in.
- Codex's decisive question, within ls_scaling's declared domain (homogeneous, classical, instant release and reheating,
  supplied phase, Y-relaxed preparation, case-L proxy, no vector/Goldstone production): **the implication "enough charge
  ⇒ excessive processed late debris" is refuted** by the m = 15–40 points. Outside that domain it is not checked.

## 4. Verdicts

| Claim | Verdict |
|---|---|
| Codex: exact sources, gauge sum, n_D source (−12 Im L) | CONFIRMED [D, C] |
| Codex: 48-sample force/phase/gauge checks | CONFIRMED [C] |
| Codex: forced response, ΔV, nonzero torque at zero displacement action | CONFIRMED [D, C]; added: n_Y source vanishes on the forced branch |
| Codex: ρ_phase = n_X²/(2f²) + n_X²/(72h²) | CONFIRMED [D] |
| Codex: yield-bound table, r_i^{7/2}, extrapolations 18.24 / 11.87 | CONFIRMED [C] |
| Codex: quiet regime unproductive at the fixed couplings | CONFIRMED and sharpened: ε ≲ 0.005, ≲ 3.5e-3 Y_Q for M_Y/H_i ≥ 30 [C] |
| Codex: a nonzero source does not imply residual heavy quanta | CONFIRMED with a dynamical instance: 1.3e-4 quanta per unit charge at the spot point; 77 with Codex's preparation [C] |
| Codex's bound applied to the ls_scaling spot point | SATISFIED, margin 2.65 (R_lead) / 4.07 (actual r_i) [C]; the leading radial relation is the only strained premise (8%) |
| Implication "enough lasting charge ⇒ excessive surviving debris", within ls_scaling's domain | REFUTED [C] (own release: 1.00 Y_Q, ε = −0.146, R = 1.14 vs R_L = 88 [Q]) |
| ls_scaling spot point and ladder (ε, N_Y, R, yield) | CONFIRMED [C] with an independent integrator at m = 10, 30, 40 and the Codex-preparation control |
| ls_scaling: ε ≈ (7–10) η | CONFIRMED: 7.8 (m = 30, 40), 6.0 (m = 10) [C] |
| Case-L thresholds R_L; ledger widths; vector/Goldstone production; phase holding; strings; running; non-instant release | NOT CHECKED (quoted from sp_critic / ls_scaling or open there) |

**Corrections for Codex.** (1) Write R with the actual prepared r_i; the leading relation is 5–8% low whenever the β6
companion matters, and the envelope then fails by 8% in amplitude. (2) The bound is a factor ≈ 2.6 estimator, not only
an inequality: the triangle inequality is 98% saturated in the first e-fold; the slack is ⟨|X|⁶⟩/R⁶ ≈ 0.4. (3) State the
fixed-coupling conclusion with its domain: raising α6 in H_i units (r, κ, smaller M_Y/H_i at fixed heaviness) and h/H_i
reaches Y_Q inside the quiet regime, which is ls_scaling's window.

## 5. Limits

- Same benchmark assumptions as ls_scaling; the window's open gates (release-induced vectors/Goldstones, phase holding,
  strings, running, instant release and reheating) are unchanged and not checked here.
- R_L and the decay-averaged ledger are quoted from our own earlier tracks, not re-derived.
- The m = 15, 20 rows are our closed-form evaluation on ls_scaling's data, not re-run.
- Our release uses a = 22–30 window means (a_end = 30); ls_scaling used a = 30–40 (a_end = 40). Agreement to 0.1%; the
  a_end = 100 control (window a = 75–100) gives a³n_X = −2788.19, yield 1.00002 Y_Q, R_central = 1.1440.

## Reproduction (about 4 min on 4 cores)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/hw_envelope
L=research/reality_framework/grut_formation_review_20261025/tracks/ls_scaling
python3 -B $T/code/controls.py $T/results/controls.json                       # 2 s
python3 -B $T/code/spot_bound.py $L/results/scan_B.json $T/results/spot_bound.json
python3 -B $T/code/release_spot.py $T/results/release_m30_relaxed.json m=30 a_end=30       # 40 s
python3 -B $T/code/release_spot.py $T/results/release_m30_codexprep.json m=30 a_end=30 prep=codex
python3 -B $T/code/release_spot.py $T/results/release_m10_relaxed.json m=10 a_end=30
python3 -B $T/code/release_spot.py $T/results/release_m40_relaxed.json m=40 a_end=30
python3 -B $T/code/release_spot.py $T/results/release_m30_relaxed_a100.json m=30 a_end=100  # control
```
