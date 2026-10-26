# Track rp_reheating: independent check of Codex e160356 handoff item 1 (physical curvature reheating)

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. Packet read only through
`git archive e160356` into the scratchpad (`origin/codex/grut-curvature-boundary-20261026`,
`research/reality_framework/grut_reheating_preparation_20261026/`). No Codex code was run; Codex's numbers are quoted
[Q] from its REHEATING_RESULTS.json / INDEPENDENT_REHEATING.json / REHEATING_VERIFICATION.json. Nothing was committed.
The kinetic tracks (kin_*) were not touched.

**Scope.** One supplied tree-level working model at the heavy-mediator spot point (H_i = 1.3334962e11 GeV, M_Y = 30 H_i,
h = 40 H_i, μ = 7.15945616 H_i, κ = 3.458684721e-3, λ_X = 1e-3, λ_Y = λ_Σ = 1, ξ_X = −5/6 ↔ c_H = 10; light sector and
Y_Q = 4.8294e-9 unchanged), homogeneous and classical, with a quadratic inflaton as a post-inflationary surrogate. Nothing
here is evidence that X, Y, Σ or the inflaton exist; it is not a validation, an exclusion or a viable point. Labels:
[D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued.

**Own code** (`code/`): `rpreheat.c` (cosmic time, physical fields, H *reconstructed from the Friedmann constraint* at
every evaluation and the Raychaudhuri equation integrated into an auxiliary H_aux as the geometry check; Dormand–Prince
5(4) with the coefficient table of our hw_verify/hwrel.c), `rp.py` (own preparation solve from an analytic starting
guess, driver, yield and dilution bookkeeping), `rp_scipy.py` (a second formulation with H evolved, as Codex does, under
SciPy DOP853 and Radau), `runs.py`, `analysis.py`. Reuse from our earlier tracks, stated: the DP5(4) tableau and the
Y free-mode diagnostic definition (hw_verify), the decay ledger and its (Rr/3)^{1/4} reheating correction (hw_verify,
hw_gates), the yield convention (ls_scaling). Codex's engine differs: C DOP853 with H evolved and Friedmann checked.

## 1. The equations, re-derived [D]

Action in signature (−+++): S = ∫√−g [½F R − ½(∂I)² − Σ|∂z|² − ½m_I²I² − V_h + L_r], F = M_r² − 2ξ|X|², so the
nonminimal term is −ξR|X|² in L, i.e. +ξR|X|² in the effective potential: −10H² in de Sitter (R = 12H²), 0 in ideal
radiation, −2.5H² in ideal matter (R = 3H²). Varying g gives F G_μν = T_μν + ∇_μ∇_νF − g_μν□F with T from everything
but the ½FR term. In flat FLRW (∇_0∇_0F = F̈, □F = −F̈ − 3HḞ, ∇_i∇_jF = −HḞ g_ij):

- 00: **3FH² + 3HḞ = ρ**, ρ = ½İ² + ½m_I²I² + Σ|ż|² + V_h + ρ_r. The 3HḞ term is the only ξ-dependent correction to
  the constraint; it is present in Codex's `3fh_c² + 3h_cf' = ρ`.
- ij: −F(2Ḣ + 3H²) = p + F̈ + 2HḞ; with the constraint, **Ḣ = −(ρ + p + F̈ − HḞ)/(2F)**, ρ + p = İ² + 2Σ|ż|² + 4ρ_r/3.
  Codex's `h_c' = −(j² + 2η²K + 4r/3 + f'' − h_cf')/(2f)` is this equation.
- trace: −FR = T − 3□F → FR = ρ − 3p − 3F̈ − 9HḞ. With Ḟ = −4ξRe(X*Ẋ), F̈ = −4ξ[|Ẋ|² + Re(X*Ẍ)] and the X equation
  Ẍ = −3HẊ − ∂V/∂X* − ξRX, the H·Re(X*Ẋ) terms cancel exactly and
  **(F + 12ξ²|X|²) R = ρ − 3p + 12ξ|Ẋ|² − 12ξ Re(X* ∂V/∂X*)**, ρ − 3p = (m_I²I² − İ²) + (4V_h − 2Σ|ż|²).
  This is Codex's R formula term by term (its `4U − j²` is m_I²I² − İ² in M_r² H_i² units).
- fields: z̈ + 3Hż + ∂V/∂z* + ξRX δ_zX = 0; Ï + (3H + Γ)İ + m_I²I = 0; ρ̇_r = Γİ² − 4Hρ_r (the friction term transfers
  exactly Γİ²).
- energy and Bianchi: d(K + V_h)/dt = −6HK − 2ξR Re(X*Ẋ) (work done by the curvature mass); ρ̇ = −3H(ρ + p) + ½ḞR;
  d/dt(3FH² + 3HḞ) = −3H(ρ + p) + 3Ḟ(2H² + Ḣ) = −3H(ρ + p) + ½ḞR = ρ̇, so the constraint is preserved by the
  Raychaudhuri + matter equations: the Bianchi identity closes with no extra term.

Nothing is missing from MATHEMATICS.md and no sign differs. The one remark to add: during the oscillation era R/H² =
3 − 9w_I runs over [−6, +12] (w_I from +1 to −1), so **ξR/H² runs from −10 to +5** with matter average −2.5; the realised
ranges in the runs are −10 … +4.3 (m_I = 2) and −10 … +4.7 (m_I = 10) [C]. The time average of R/H² over τ = 5–15 is
1.5–1.6 (less than 3 because radiation is already 20–50% of ρ).

Numerical verification of the derivation: my C engine (constraint-solved H) and my scipy formulation (H evolved) both
reproduce Codex's recorded 17 comparable state variables at τ = 10 to ≤ 3.6e-12 and ≤ 9.9e-12 relative
(`results/scipy_short.json`); over the full b2_g0.3 run (τ = 83.3) the two formulations agree to 4.6e-10 with Friedmann
residual 4e-11 (`results/scipy_b2_g0.3.json`). Own preparation solve (Newton from x² = c_H/2λ_X, Y = −κx³/M², s = h/√2,
with the actual R(0) = 12.000) lands on Codex's initial state to the last printed digit (76.39978266076902,
−2.111509854283515 − 0.9393290404896636 i, −28.293142125868524 i; residual 1.2e-12).

## 2. Reproduction of the six runs and the controls [C]

Own engine, rtol 1e-10 / atol 1e-13, late means τ-weighted over the last 20% as in Codex. `results/reproduce.json`.

| run | ε ours / Codex | N_X ours / Codex | N_YH | N_YL ours / Codex | yield/Y_Q ours / Codex | \|N_X\|/2788.18 | D = (ρ_r a⁴/3)^{3/4} | a_f |
|---|---|---|---|---|---|---|---|---|
| m_I = 2, Γ = 0.03 | −0.19128 / −0.19128 | −12425.5 / −12425.5 | 2.016e-05 | 3.5043e-03 / 3.5043e-03 | 0.3327 / 0.3327 | 4.456 | 13.40 | 96.30 |
| 2, 0.1 | −0.19877 / −0.19877 | −12263.8 / −12263.8 | 2.016e-05 | 3.5065e-03 / 3.5065e-03 | 0.5596 / 0.5596 | 4.399 | 7.86 | 44.19 |
| 2, 0.3 | −0.21868 / −0.21868 | −12041.9 / −12041.9 | 2.015e-05 | 3.5116e-03 / 3.5116e-03 | 0.7843 / 0.7843 | 4.319 | 5.51 | 22.70 |
| 10, 0.03 | −0.18148 / −0.18148 | −3309.4 / −3309.4 | 3.245e-03 | 1.7398 / 1.7398 | 0.1847 / 0.1847 | 1.187 | 6.43 | 75.39 |
| 10, 0.1 | −0.18153 / −0.18153 | −3314.6 / −3314.6 | 3.197e-03 | 1.7071 / 1.7071 | 0.3250 / 0.3250 | 1.189 | 3.66 | 34.25 |
| 10, 0.3 | −0.17971 / −0.17971 | −3339.2 / −3339.2 | 3.068e-03 | 1.6142 / 1.6142 | 0.5182 / 0.5182 | 1.198 | 2.31 | 17.01 |

Worst relative difference over all compared quantities (ε, N_X, N_YH, N_YL, yield, a_f, H_f, eight runs): 1.4e-4 (the
N_YH of the ξ = 0 control); everything else ≤ 3e-5. Controls: rtol 1e-12 changes N_X, E_X, ε by 2e-12 (Codex 3e-10);
rtol 1e-8 by 1e-11; phase ±1e-8: response 3.9e-10 (Codex 4.4e-10), d ln N_X/dδ = 0.0386 (Codex 0.0386); ξ = 0 control
ε −0.15982, 0.1184 Y_Q (Codex same); long run τ = 500: ΔN_X/N_X = 1.40e-4, Δε/ε = 2.21e-3 (Codex 1.40e-4, 2.21e-3);
τ = 1500: 2.4e-4 and 3.4e-4 with ε = −0.19884, so ε wanders in [−0.1983, −0.1988] and "about −0.199" is the right
quotation. Geometry: |H_aux − H|/H ≤ 5e-13 (m_I = 2) and ≤ 1.9e-8 (m_I = 10, τ = 833; the fast oscillation accumulates
a secular error, the mirror of Codex's 1.27e-9 Friedmann residual); Gauss ≤ 1.4e-12; hidden work identity ≤ 1e-12.
RD reference (de Sitter-prepared state released into pure radiation, η² = 0, a-window 30–40): ε −0.1457174,
N_X −2788.1768, N_YH 0.0172529, N_YL 0.334465 against the parent's late values to ≤ 1e-5 (and hw_verify's instant
release). All of Codex's reported numbers for item 1 reproduce.

## 3. The yield deficit: dilution of a charge that is fixed early [C, D]

Decomposition yield/Y_Q = (|N_X|/2788.18) / D with D = (ρ_r a⁴/3)^{3/4} the comoving-entropy ratio to the instant release
(hw_verify's bookkeeping [D]). Profiles (`results/reproduce.json`, scratch .npy): the charge is generated during the first
roll-down of X while R falls from 12H² through zero — 90% of the final |N_X| by τ = 2.6 (m_I = 2) or τ = 0.5
(m_I = 10), when radiation is 15% / 5% of ρ and |X| is still 32–38 — and 99% by τ ≈ 3, long before reheating (radiation
fraction ½ at τ ≈ 9.5 for Γ = 0.1). The entropy is produced afterwards: D = 2.3–13.4 in the six runs. Two effects set
the charge ratio: the time-averaged curvature mass (−2.5H²) alone gives the instant-release charge (fluid inflaton with
ξ⟨R⟩, de Sitter-prepared: |N_X|/2788 = 1.02, ε −0.17 at every Γ); the oscillation of R adds 19% (m_I = 10); the slow
first roll of a m_I = 2 inflaton (w ≈ −1 for τ ≲ 0.8, R crossing zero at τ ≈ 0.9) keeps |X| large longer and gives
4.3–4.5× the charge but also ≈ 0.26 extra e-folds of quasi-de Sitter expansion, which raises D by 2.2 over the pure
matter-era value 0.926 (1.5/Γ)^{1/2} (analytic, checked at m_I = 10: 3.59 vs 3.66). Net: m_I = 2 loses to dilution
(4.4/7.9 = 0.56); m_I = 10 has little extra charge and loses to dilution (1.19/3.66 = 0.33).

**Does a different m_I or Γ_I restore the yield without retuning?** Yes, inside the model's own two parameters and inside
perturbative decay (Γ_I ≤ 0.3 m_I). `results/gamma_scan.json`, `results/extra.json` (yr_*):

| m_I/H_i | Γ_I/H_i (Γ/m_I) | |N_X|/2788 | D | yield/Y_Q | ε | N_YL | R_central [E] |
|---|---|---|---|---|---|---|---|
| 3 | 0.5 (0.17) | 2.93 | 2.80 | **1.045** | −0.199 | 0.0068 | 0.008 |
| 5 | 1.0 (0.20) | 1.88 | 1.79 | **1.047** | −0.185 | 0.018 | 0.03 |
| 10 | 2.0 / 3.0 (0.2 / 0.3) | 1.28 / 1.32 | 1.37 / 1.30 | 0.937 / **1.019** | −0.174 / −0.173 | 0.99 / 0.76 | — / 1.7 |
| 15 | 4.5 (0.3) | 1.20 | 1.20 | 1.004 | −0.164 | 0.74 | 2.2 |
| 20 | 6 (0.3) | 1.15 | 1.15 | 1.007 | −0.160 | 0.42 | 1.9 |
| 30 | 6 / 9 (0.2 / 0.3) | 1.10 / 1.11 | 1.13 / 1.10 | 0.980 / 1.014 | −0.157 / −0.156 | 0.27 / 0.29 | — |
| 50 | 10 / 15 (0.2 / 0.3) | 1.08 / 1.08 | 1.08 / 1.06 | 1.003 / 1.018 | −0.154 / −0.152 | 0.24 / 0.27 | 0.8 / — |

The yield-restoring line is Γ_I ≈ 0.2–0.3 m_I for every m_I from 3 to 50 H_i (a matter era of at most about one
Hubble time), with ε = −0.15 to −0.20 and the light sector, scale and phase untouched. Codex's grid had Γ_I/m_I ≤ 0.15
(m_I = 2) and ≤ 0.03 (m_I = 10). Rows with Γ_I > 0.3 m_I (overdamped inflaton, e.g. m_I = 2 with Γ ≥ 1) are kept in the
JSON but flagged; they are outside perturbative decay and at Γ ≫ m_I the field creeps with w ≈ −1 (m_I = 2, Γ = 30:
a_f = 1.3e4), which is not a reheating history. Reheating temperature on the line [E]: T_RH ≈ 0.54 √(Γ_I M_r) =
2.2e14 GeV at Γ_I = 0.5 H_i, 5.3e14 GeV at 3 H_i (Codex's T_vis,final is the temperature at the end of the run, not T_RH).

**Relation to hw_verify's finite release** (Γ_φ/H_i = 1: ε −0.279, R 0.032, diluted 1.16 Y_Q [Q own]). Its model is a
fluid inflaton with the X mass −c_H ρ_φ/3 = −10H² *throughout* the matter era. My mode 2 reproduces it: Γ = 0.1 / 0.3 /
1 / 3 / 10 → ε −0.3809 / −0.3465 / −0.2782 / −0.2137 / −0.1731 and yields 0.545 / 0.836 / 1.157 / 1.304 / 1.252
(hw_verify −0.3804 / −0.3456 / −0.2791 / −0.2135 / −0.1720; 0.545 / 0.835 / 1.16 / 1.30 / 1.24) — no numerical
discrepancy. The physical difference is the Hubble mass during the matter era: the covariant ξ_X R|X|² gives −2.5H² on
average (oscillating −10 … +5), not −10H², so the slow-release charge enhancement is 1.0–1.3× (m_I ≥ 5) or 4.4×
(m_I = 2, from the slow first roll) instead of 1.75–2×, and |ε| stays at 0.15–0.22 instead of rising to 0.38. At the
same Γ the yields nearly coincide (Γ = 0.1: 0.545 hw_verify, 0.560 Codex m_I = 2, 0.325 m_I = 10; Γ = 0.3: 0.835 /
0.784 / 0.518); the "1.16 Y_Q" is at Γ = H_i, which Codex did not run, and in the covariant model Γ = H_i gives 0.78
(m_I = 10). hw_verify's finite-release row is therefore CORRECTED in physical content (its −c_H H² matter-era mass is
not what ξ_X R gives) while its numbers stand in its stated model.

## 4. Preparation, the oscillating R and the Y excitation [C]

The initial state is ls_scaling's Y-relaxed conditional minimum at the held phases (arg X = 0, arg Σ = −π/2) with the
de Sitter curvature mass ξR(0) = −10H² (R(0) = 12.000 since the inflaton starts at rest), zero velocities, a = 1, H = H_i,
and the inflaton amplitude fixed by the constraint (I₀/M_r = 1.2247 for m_I = 2). No holder term is in the dynamics:
the phases are held by fiat up to τ = 0 and released together with the inflaton, so the preparation assumes a holder that
switches off exactly when the de Sitter stage ends. The phase value remains supplied (hw_gates/hw_verify: unheld, it
drifts 0.25 rad per e-fold during inflation). The state carries 0.048 light and 0.0017 heavy free Y quanta at τ = 0
(the quartic and Σ back-reaction displace it from the forced solution; hw_verify's 0.047), i.e. it is the quiet
preparation; nothing new about the holder is established by the packet, as it says.

**Does the oscillating R excite Y?** Yes, through X (ξR couples only to X; κX³ drives Y) and resonantly. Scan at
Γ = 0.1 H_i (`results/b_scan.json`): N_YL = 0.0026, 0.0035, 0.0071, 0.0022, 0.028, 0.16, 0.41, 0.76, 1.27, 1.71, 2.03,
2.07, 2.15, **2.22**, 1.91, 1.52, 1.15, 0.87, 0.68, 0.55, 0.47, 0.41, 0.34, 0.28 for m_I/H_i = 1, 1.5, 2, 3, 4, 5, 6, 7, 8,
9, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25; N_YH peaks at 0.16 for m_I = 19. The light-Y
frequency is ω_L = (M_Y² − 2μ|Σ|)^{1/2} = 22.25 H_i and the heavy one ω_H = 36.1 H_i; the driving frequency of R is
2m_I, so the peaks at 2m_I = 24 and 38 are the two Y modes (the shift from 22.2/36.1 is within the Hubble-broadened
width). Codex's m_I = 10 sits on the light-mode resonance shoulder, m_I = 2 is far off it: that is why "m_I/H_i = 10
produces more Y excitation". The occupation is set by τ ≲ 2 while |X| ≳ 30, then constant to 1e-4 over the last 20%
("undamped" = conserved free-mode number; no Y decay is in the equations, and Γ_L ≈ 1e-7 H_i would not matter on these
time scales). N_YL depends weakly on Γ at fixed m_I (1.74 / 1.71 / 1.61) because the inflaton decays earlier.

**Debris implication [E, reused hw_verify ledger with (Rr/3)^{1/4} on ⟨a⟩]**: per light/heavy Y quantum the deposit is
8.79e4 / 7.54e5 E_X-units at a = 1 (m_L = 2.97e12 GeV, Γ_L = 1.34e4 GeV; m_A = 3.2e6 GeV, Γ_A = 2.12e-10 GeV). m_I = 2
runs: R_central = 0.004 (quiet, 300× below the instant 1.14). m_I = 10 runs: R_central = 8.0 / 6.5 / 5.2 for Γ = 0.03 /
0.1 / 0.3 — passes the isentropic proxy (R_L ≈ 82, larger at ε = 0.18) by ≥ 10 but sits at the energy-conserving proxy
(7–8 at ε = 0.146): the resonant Y excitation turns a quiet release into a marginal one under that proxy. On the
yield-restoring line: R = 0.008 (m_I = 3), 0.03 (5), 1.7 (10), 2.2 (15), 1.9 (20), 0.8 (50). The quiet-and-sufficient
corner of this model is therefore m_I/H_i ≲ 5 with Γ_I ≈ 0.2–0.3 m_I, or m_I ≳ 30 with Γ_I ≳ 0.2 m_I. These are
late free-mode occupations converted with the decay ledger, not a particle-production or kinetic calculation.

## 5. Verdicts

| # | Codex claim (e160356 item 1) | Codex value [Q] | This track (own code) | Verdict |
|---|---|---|---|---|
| 1 | Jordan-frame equations: Friedmann with 3HḞ, Raychaudhuri, trace-eliminated R, field/radiation transfer, work identities, Bianchi | MATHEMATICS.md | re-derived line by line; no missing ξ term, no sign difference; two own formulations agree with Codex's recorded states to ≤ 1e-11 (τ = 10) and with each other to 5e-10 (τ = 83) [D, C] | CONFIRMED |
| 2 | Curvature mass ξR with ξ = −5/6: −10H² dS, 0 RD, −2.5H² matter; oscillates in the inflaton era | stated qualitatively | ξR/H² ∈ [−10, +5] (w_I ∈ [−1, 1]); realised −10 … +4.3/+4.7; τ-average 1.5–1.6 R/H² at τ = 5–15 [D, C] | CONFIRMED (range made explicit) |
| 3 | Six runs: ε, N_YL, yields (table) | −0.1913…−0.2187 / −0.1815…−0.1797; 0.00352–0.00353 / 1.61–1.74; 0.333–0.784 / 0.185–0.518 | identical to all printed digits; worst relative difference 1.4e-4 over eight runs × seven quantities [C] | CONFIRMED |
| 4 | Radau vs DOP853 3.21e-12; tighter tolerance 3e-10 | as stated | own DOP853/Radau/DP5 vs Codex states ≤ 9.9e-12; tighter 2e-12 [C] | CONFIRMED |
| 5 | Friedmann residual ≤ 1.27e-9; Gauss 8e-14; work 2e-14 | as stated | mirrored check \|H_aux − H\|/H ≤ 5e-13 (m_I = 2), 1.9e-8 (m_I = 10, τ = 833); Gauss 1.4e-12; work 1e-12 [C] | CONFIRMED (secular growth with m_I noted) |
| 6 | Long run: charge 0.014%, ε 0.22%; quote ε ≈ −0.199 | as stated | 1.40e-4 and 2.21e-3 at τ = 500; τ = 1500: ε = −0.1988; band [−0.1983, −0.1988] [C] | CONFIRMED |
| 7 | Phase ±1e-8: 4.36e-10 relative | as stated | 3.9e-10; slope 0.0386 [C] | CONFIRMED |
| 8 | ξ = 0 control, RD reference reproduces old late values to 4e-7 | −0.1598 / 0.118; errors ≤ 6e-7 | −0.15982 / 0.1184; RD errors ≤ 1e-5 (window) [C] | CONFIRMED |
| 9 | Undamped Y occupation 0.00353 (m_I = 2); "much larger" at m_I = 10 | 0.00352–0.00353; 1.62–1.74 | same; explained as the 2m_I ≈ ω_YL (22.25 H) resonance, peak 2.22 at m_I = 12; heavy mode peak 0.16 at m_I = 19 [C] | CONFIRMED and explained |
| 10 | None of the six reaches Y_Q; no retuning performed | 0.18–0.78 Y_Q | true for the six; the same model gives 0.94–1.15 Y_Q along Γ_I ≈ 0.2–0.3 m_I for m_I = 3–50 H_i without retuning (m_I = 3, Γ = 0.5: 1.045, ε −0.199, N_YL 0.0068) [C] | CORRECTED in scope: deficit is the grid, not the model |
| 11 | Cause of the deficit (handoff asks) | not analysed by Codex | entropy dilution D = 2.3–13.4 of a charge fixed by τ ≈ 0.5–3; charge ×4.4 (m_I = 2, slow first roll) or ×1.19 (m_I = 10); matter-averaged ξ⟨R⟩ alone gives ×1.02 [C] | new |
| 12 | hw_verify finite release Γ_φ/H_i = 1: ε −0.279, R 0.032, 1.16 Y_Q [own] | as stated | reproduced in its model (−0.2782, 1.157); covariant ξR gives −2.5H² not −10H² in the matter era: Γ = H_i → 0.78 Y_Q, ε −0.18 (m_I = 10) [C] | CORRECTED (physical content of our own earlier row) |
| 13 | Preparation: conditional minimum with actual R(0); phases supplied; no holder | as stated | identical state from own solve; R(0) = 12.000; holder switch-off implicitly at τ = 0; 0.048 light quanta at start [C] | CONFIRMED as conditional |
| 14 | Debris of the m_I = 10 runs (not in packet) | — | R_central 5.2–8.0 vs 1.14 instant (ledger reuse) [E] | new: marginal under the energy-conserving proxy |
| 15 | Particle production, finite-k, decays, fragmentation, thermal bath | not computed by Codex | not computed here | NOT CHECKED |

## 6. Limits

- Homogeneous, classical, tree level; no finite-k modes, mediator decays in the dynamics, backreaction of daughters,
  thermal corrections or fragmentation; the inflaton is a quadratic surrogate with a constant Γ (perturbative decay,
  no preheating). Runs with Γ_I > 0.3 m_I are flagged, not used for conclusions.
- The de Sitter stage is not evolved: the state at τ = 0 is the held-phase conditional minimum released with zero
  velocities; the holder and the phase value are supplied.
- The yield uses instantaneous thermalisation of the decay products and constant g* = 106.75; D is the comoving
  entropy relative to the instant release in the same convention. T_RH values are the standard estimate [E].
- The debris numbers reuse our hw_verify ledger (decay-averaged ⟨a⟩ in radiation, (Rr/3)^{1/4} correction) and the two
  equilibrium proxies of ls_critic/hw_critic; neither is a kinetic threshold (kin_* tracks, not touched here).
- Two formulations and three integrators, all written by the same session; this is not outside replication.

## Reproduction (≈ 1 min on 4 cores; the C engine is compiled into the scratchpad build directory)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/rp_reheating/code
cd $T
python3 -B runs.py ../results reproduce gamma_scan b_scan models extra     # 85 C runs, ~15 s
python3 -B rp_scipy.py ../results/scipy_short.json 10 2 0.1                # second formulation, tau = 10 (~10 s)
python3 -B rp_scipy.py ../results/scipy_b2_g0.3.json 83.333333333333333 2 0.3   # full run, Radau ~70 s
python3 -B analysis.py ../results                                           # tables above -> results/summary.json
```
Credit: Dine–Randall–Thomas Hubble-induced terms and Affleck–Dine [R]; scalar–tensor field equations in the Jordan
frame (standard, e.g. Faraoni) [R]; Dormand–Prince 5(4); Codex e160356 [Q]; our tracks ls_scaling, hw_verify, hw_gates,
hw_critic [Q where cited].
