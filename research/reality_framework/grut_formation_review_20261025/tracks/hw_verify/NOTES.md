# Track hw_verify: independent verification of the ls_scaling heavy-mediator window (different integrator), and the four uncertified items

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication. Resumed run (code/ and the first
four result files were produced by the earlier attempt of this track; prep, release_*, analytic, rl_actual_fine and
summary were produced in this run; everything was re-read and checked for consistency).

Packets read through `git archive` into the scratchpad only: Codex 2908e91 (`origin/codex/grut-concept-checkpoint-20261025`).
No Codex code was run. ls_scaling's numbers are quoted [Q, own track] for comparison; none of its code was imported.
The case-L threshold uses our own sp_yledger/gp_leftovers benchmark code (reuse, labelled). Nothing was committed.

**Scope.** One supplied tree-level working model: homogeneous classical release into radiation, a supplied CP phase,
the late-particle decay ledger, the case-L proxy with the 1% allowance. Nothing here is evidence that X, Y, Σ or the
vector exist; it is not a validation, an exclusion or a viable point. Labels: [D] derived, [C] computed, [E] estimate,
[Q] quoted, [A] argued.

## Method (own code, `code/`)

- `hwrel.c`: homogeneous X, Y, Σ release in **conformal time with conformal fields** χ = a z (χ'' = (a''/a)χ − a³F(χ/a),
  own derivation), **Yoshida 8th-order symmetric composition** of drift–kick–drift in extended phase space, step
  c/Ω with Ω a Gershgorin bound; a Dormand–Prince 5(4) path on the same system as an internal control. ls_scaling
  used cosmic time, physical fields and DP5(4) (sp_yledger's yrelease.c with runtime couplings); hw_envelope used
  Python DOP853; hw_gates a DP5 C code. Two backgrounds: (0) radiation from a = 1 (instant release); (1) inflaton
  matter decaying at width Γ_φ into radiation, with the Hubble-induced X mass −c_H(ρ_φ/3)|X|² (finite release).
  Gauss's law N_X + 3N_Y + 6N_Σ is exact to rounding (residual ≤ 5e-10 in every run).
- Preparations (own minimisers, Newton-polished, Hessian-checked): `codex` = aligned inflationary minimum, Σ rotated
  by δ with Y left at its δ = 0 value; `relaxed` = conditional minimum over (|X|, Y, |Σ|) at the held phases.
- Ledger [D, re-derived from Codex b448f25 / ls_scaling]: R_Y = (N_YH + N_YL)(m_L/H)⟨a⟩_L/E_X, R_Σ = 2N_σ(m_A/H)⟨a⟩_A/E_X,
  R_casc = N_YH(m_A/H)⟨a⟩_A/E_X, ⟨a⟩ = Γ(3/2)√(2H_i/Γ), Γ_L = 3κ²m_L/256π³, Γ_A = e²m_A(1 − 4m_X²/m_A²)^{3/2}/48π.
- Thresholds: (a) R_L(ε) from our sp_critic case-L table (log-log interpolation; the value every hw track quotes);
  (b) **direct** case-L runs at the actual ε = 0.1457 and at both yields (`rl_actual.py`, our thermal_relic code).
- Yield: |a³n_X| H_i³/s(T_R), prompt radiation, g* = 106.75, M_red = 2.435e18 GeV (the convention of ls_scaling and
  Codex); Y_Q = 4.8294e-9. c6 = α6 M_red², α6 = μκ²⟨Σ⟩/(M_Y⁴ − 4μ²⟨Σ⟩²).

## 1. Reproduction of the ladder with the different integrator [C]

c_H = 10, r = 0.45, h = 40H, η = 0.5 η_edge (κ = 1.1529e-4 m), Y-relaxed, δ = π/2, H_i = 4.3487e11 GeV, a_end = 100,
window a = 30–100 (ls_scaling: a_end = 40, window 30–40). `results/ladder.json`.

| m = M_Y/H | ε ours / ls | N_YH | N_YL | N_σ | R_Y | R_Σ | R_casc | R_central ours / ls | R_L table / ls | yield/Y_Q ours / ls |
|---:|---|---:|---:|---:|---:|---:|---:|---|---|---|
| 10 | −0.11290 / −0.113 | 2.098 | 20.32 | 3e-10 | 98.8 | 1e-8 | 45.8 | **144.5 / 145** (fails) | 59.8 / 60 | 4.22 / 4.2 |
| 15 | −0.13375 / −0.134 | 0.402 | 7.234 | 3e-9 | 25.9 | 1e-7 | 8.29 | 34.2 / 34 | 77.2 / 77 | 5.21 / 5.2 |
| 20 | −0.14233 / −0.142 | 0.1257 | 2.124 | 9e-9 | 6.43 | 4e-7 | 2.51 | 8.94 / 8.9 | 84.8 / 85 | 5.67 / 5.7 |
| **30** | **−0.14567 / −0.146** | 0.01725 | 0.3345 | 3e-9 | 0.805 | 1e-7 | 0.339 | **1.144 / 1.14** | 87.8 / 88 | 5.89 / 5.9 |
| 40 | −0.14541 / −0.145 | 0.00387 | 0.0849 | 6e-7 | 0.175 | 2e-5 | 0.0757 | 0.251 / 0.25 | 87.6 / 88 | 5.89 / 5.9 |

c6 = 1.1757e5 at H_i (ls 1.18e5). Spot point (m = 30): N_X = −2788.1, E_X = 38,401. **Rescaled by the symmetry to
the yield-normalised scale** (s = 5.889^{−2/3} = 0.30664, run directly in GeV: `spot_Hnorm`): H = 1.3335e11,
M_Y = 4.0005e12, h = 5.3340e12, μ = 9.5471e11 GeV, κ = 3.4587e-3 (0.701 κ_max), ε = −0.14567, yield 1.0000 Y_Q,
R = 1.144, c6 = 1.2504e6. All ls_scaling spot numbers reproduced.

Controls at the spot point: Yoshida c = 0.1 / 0.05 / 0.025, DP5(4) rtol 1e-11 / 1e-13, a_end = 300 (window
100–300): ε = −0.145669 and R = 1.1440 in all six; Codex preparation: ε = −0.15882, N_YH = 76.8, N_σ = 0.039,
**R = 1415** (ls 1414, hw_envelope 1415). The δ = π/4 rungs are in ladder.json (R 124 → 0.36, same pattern).

Heavy-limit slopes between rungs (m ≥ 15): N_Y ∝ m^{−4.2, −4.6, −4.8}, R ∝ m^{−4.7, −5.1, −5.3} (ls: 4.2→4.8 and
4.9→5.3). Window edges: debris R = R_L crossing at m = 12.3 (12.1 with the corrected R_L below); boundedness
κ = κ_max at m = 42.8 (ls: ≈ 13 to ≲ 43).

**Codex 2908e91's envelope bound, own implementation of the quoted formula** (`analytic.py`): the four fixed-coupling
values 1.2201e-8 / 5.4890e-8 / 3.8583e-5 / 1.7358e-4 Y_Q and the extrapolation r_i/M_Y = 11.87 (c_H = 10) reproduced
to the printed digits; that extrapolation sits at M_Y/H_i = 5.96, the resonant regime. On the ladder (R_env from the
actual prepared |X_i|): |N_X|/bound = 0.20, 0.23, 0.24, 0.24, 0.24 for m = 10–40 (yield bound 21–25 Y_Q at H_i against
4.2–5.9 simulated); Codex's matching diagnostics at m = 30: κR²/M² = 0.025, ω_X/gap = 0.16, ω_X/m_σ = 0.063,
Y-quartic 0.05, Σ shift 1e-3 (at m = 10 the Y-quartic diagnostic is 3.5: outside matching). Agrees with hw_envelope.

## 2. R_L at the actual ε and yield [C, direct runs of our case-L code]

`results/rl_actual.json`, `rl_actual_fine.json`. ε = 0.1457, case L, f(R):

| yield | R = 70 | 80 | 82 | 83 | 85 | 90 | 91 | 93 | **R_L (f = 0.01)** |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 Y_Q | 1.2e-10 | 0.00303 | 0.0150 | 0.0293 | 0.0916 | 0.638 | | | **81.5** (bracket 80–82) |
| 5.889 Y_Q (raw at H_i) | 4e-34 | 3.1e-10 | | | | 0.00717 | 0.0167 | 0.0694 | **90.4** (bracket 90–91) |

The table value 87.8 quoted by every hw track is **7.8% high**: the sp_critic table brackets f = 0.01 between R = 75
and 100 at ε = 0.15 and ln f(R) is strongly concave (slopes 1.7, 0.68, 0.39 per unit R over 70–80, 80–85, 85–90), so
the log-linear crossing lands high. The spot-point margin is therefore R_L/R = **71** (not 77). At the raw yield
the threshold is 90.4 (more hidden charge per entropy tolerates more debris). Codex's calibrated values were ≈ 9%
below our table [Q ls_scaling]; the direct value nearly removes that difference.

## 3. Item (1): is the Y-relaxed state the physical preparation? [C; `prep.py`, `results/prep.json`]

Spot point, δ = π/2.

**Relaxation rate.** Linearised modes about the relaxed conditional minimum with the two phases held (canonical
4×4 Hessian): ω/H = 4.09 (X radial), **22.96** (Y light Majorana mode; M_Y√(1−r) = 22.2), **36.33** (Y heavy;
M_Y√(1+r) = 36.1), 56.6 (Σ radial). All are underdamped (ω > 3H/2), so every amplitude decays as e^{−3N/2} and
every free-mode number as e^{−3N}: **the Y relaxation rate during inflation is Hubble-limited, 3H/2, independent of
M_Y**; e-folding time of the amplitude 0.67 e-folds.

**Residual displacement after N e-folds.** The Codex-prescription state differs from the relaxed one by
|ΔX| = 3.80, |ΔY| = 1.63, |ΔΣ| = 0.018 H; its displacement from the forced response is |Y − Y_*| = 1.47 H (relaxed:
0.047), carrying N_YH = 77.8 free heavy quanta at a = 1, which is what the instant release emits (76.8). De Sitter
evolution from the Codex state with the phases held (DOP853, then instant release carrying the acquired velocities):

| N (e-folds) | 0 | 0.1 | 0.3 | 0.5 | 0.75 | 1 | 1.5 | 2 | 3 | 4 | 5 | relaxed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free N_Y at release | 77.9 | 61.8 | 33.3 | 16.9 | 8.3 | 3.9 | 1.0 | 0.29 | 0.065 | 0.049 | 0.050 | 0.050 |
| ε | −0.159 | −0.152 | −0.141 | −0.138 | −0.143 | −0.142 | −0.146 | −0.145 | −0.146 | −0.146 | −0.146 | −0.146 |
| R_central | 1415 | 1163 | 684 | 385 | 204 | 88.3 | 22.8 | 7.87 | 1.65 | 1.06 | 1.12 | 1.144 |

R drops below R_L after **N = 1.02 e-folds** of held-phase relaxation, is within 50% of the relaxed value by N = 3 and
indistinguishable by N ≥ 4 (the ±5% at N = 4, 5 is the X-radial mode, ω = 4.1 H, still ringing). Fitted
d ln N_Y/dN = −2.4 over N = 0.2–2 (−3 for a pure mode; the slower rate is the X-radial admixture). Controls: the relaxed
state held for 0.5 and 2 e-folds stays at R = 1.144; de Sitter quantum excitation of a 30 H field is e^{−2πm/H} = 1e-82
[E]; even the massless cap (H/2π)² would give 0.76 quanta per H_i³ (R ≤ 17 if all heavy).

**Intermediate preparations** (zero velocity, linear interpolation Codex → relaxed, fraction f): N_YH = 77.8(1 − f)²
to 20%; R = 1415, 1286, 1163, 824, 376, 93, 13.3, 1.14 at f = 0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1. **The pass needs
f ≥ 0.75**, i.e. a residual |Y − Y_*| ≤ 0.37 H or N_YH ≤ 4. Allowance at the spot point with the direct R_L: 3.7 heavy
quanta (each costs 21.9 R: 2.29 own decay + 19.6 cascade vector) or 35 light quanta (2.29 each) per comoving H_i³,
i.e. 1.3e-3 / 1.3e-2 per unit charge.

**What the holding assumption hides.** With the phases free the relaxed state is not stationary: dV/dθ_X = −8693 H⁴
and dV/dθ_Σ = +1449 H⁴ (slow-roll drift 0.25 and −0.30 rad per e-fold; phase-sector mass² 1.0 and 5.0 H²). Released
after free evolution for N = 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2 e-folds: ε = −0.153, −0.163, −0.178, −0.204, −0.239,
−0.185, **+0.039** (R 1.16–1.49). The supplied phase is not kept at the 10% level in ε for even 0.1 e-fold without a
holder (ls_scaling §5.2; hw_gates gate 2 quantifies the holder).

**Verdict.** The Y-relaxed state is the inflationary attractor of a heavy Y **given** held phases and a holder that
does not enter the Y equation: a displaced Y relaxes at 3H/2 and the debris verdict flips after 1 e-fold. It is not
"what a heavy Y does" unconditionally: nothing in the stated action holds the phase, and a holder that enters the Y
equation (the DRT-type A-term, hw_gates D1) makes the Codex state the attractor instead; then the switch-off rate,
not the preparation, decides (hw_gates). CORRECTED in scope, as hw_gates found by a different route; the new numbers
are the rate, the e-fold count and the tolerance |Y − Y_*| ≤ 0.37 H.

## 4. Item (2): robustness across the supplied phase and the edge fraction [C; `phases.py`, `results/phases.json`]

**Phase, Y-relaxed** (δ = 5°…175° in 5° steps): R_central = 0.59–1.47 at every phase; R/R_L ≤ 0.040 wherever R_L is
defined (|ε| ≥ 0.068: δ = 20°–145°). |ε| ≥ 0.10 for δ ∈ [25°, 120°], ≥ 0.13 for [40°, 100°]; best |ε| = 0.1546 at
δ = 70° (R = 1.30, R_L = 96.5). CP mirror δ = 270°: ε = +0.14567 exactly. ls_scaling ("15°–165°: R ≤ 1.5; best
ε −0.154, R_L 96"): CONFIRMED. **Codex preparation** across phases: R = 51–2700, R/R_L ≥ 1.39 at every phase where
R_L is defined — the Codex-prepared heavy point fails at every phase, not only at π/2.

**Edge fraction** (δ = π/2, relaxed): frac 0.1 / 0.25 / 0.5 / 0.75 / 0.9 / 1.0 / 1.2 → ε = −0.031 / −0.074 / −0.146 /
−0.229 / −0.294 / −0.352 / −0.626 and R = 0.55 / 0.86 / 1.14 / 1.27 / 1.35 / 1.49 / 2.08 (R_L 33 at 0.25, 164 at 0.75;
beyond the table for |ε| > 0.24, where R_L(0.24) = 172 and rising). κ/κ_max = 0.99 at frac 1.0 and 1.09 at 1.2
(marginal minimum, smallest Hessian eigenvalue 0.02 H²); at frac ≥ 1.3 (κ ≥ 1.13 κ_max) **no aligned inflationary
minimum exists** (only saddles and the unbounded direction with |Σ| → 0), so the preparation fails. ls_scaling
("0.25–0.9: ε 0.074–0.29, R ≤ 1.4"): CONFIRMED and extended; at m = 30 the usable edge fraction is bounded by
κ = κ_max near frac ≈ 1 (ε ≈ −0.35), not by the companion cap.

## 5. Item (3): the exact scaling symmetry [D, C; `scaling.py`, `results/scaling.json`]

[D] In H_i units the release depends only on λ's, κ, c_H, M_Y/H, h/H, μ/H; the ledger deposits on m_L/H, m_A/H and
H/Γ with Γ_L ∝ m_L, Γ_A ∝ m_A; so ε, N_i, E_X and every R are invariant under (M_Y, h, μ, H_i) → s(…). Yield ∝ H^{3/2}
through s(T_R) ∝ (H M_red)^{3/2}; c6 ∝ α6 ∝ H^{−2}. Two conditions: constant g*, and (1 − 4m_X²/m_A²)^{3/2} ≈ 1
(1 − 1e-14 at s = 0.1, where m_A = 1.04e6 GeV).

[C] Re-run in GeV at s = 0.1, 0.30664, 3, 10 (couplings scaled and re-divided by H, a non-trivial floating-point path):
|Δε| ≤ 4.2e-13, |ΔN_X/N_X| ≤ 1.6e-12, |ΔR/R| ≤ 1.3e-10, yield ratio/s^{1.5} − 1 ≤ 1.6e-12, c6 ratio·s² − 1 ≤ 3e-16.
**CONFIRMED.** One consequence the symmetry does not cover: R_L at the *actual* yield moves with s (81.5 at Y_Q,
90.4 at 5.89 Y_Q), so "R vs R_L" is s-dependent only through the threshold, by ≈ 10%.

## 6. Item (4): does the instant-release assumption matter more or less in the heavy regime? [C; `release.py`, `results/release_*.json`]

Background 1: inflaton matter at a = 1 decaying at Γ_φ into radiation, X Hubble mass −c_H(ρ_φ/3)|X|² decaying with it,
zero initial velocities, a_end = 100, window 30–100. Late bookkeeping [D]: the comoving entropy is (Rr/3)^{3/4}
times the instant value (yield_diluted) and the ledger's ⟨a⟩ is (Rr/3)^{1/4} larger (R_reheat), with Rr = ρ_r a⁴ at the
end (3 for the instant case).

**Heavy spot (m = 30, relaxed, δ = π/2):**

| Γ_φ/H_i | ∞ (instant) | 100 | 30 | 10 | 3 | 1 | 0.3 | 0.1 |
|---|---|---|---|---|---|---|---|---|
| ε | −0.1457 | −0.1492 | −0.1563 | −0.1720 | −0.2135 | −0.2791 | −0.3456 | −0.3804 |
| free N_Y | 0.352 | 0.346 | 0.256 | 0.087 | 0.020 | 0.011 | 0.008 | 0.007 |
| R_central (R_reheat) | 1.144 | 1.07 | 0.69 | 0.20 | 0.051 (0.054) | 0.032 (0.037) | 0.029 (0.039) | 0.029 (0.045) |
| R_L (table) | 87.8 | 91.1 | 98.3 | 115 | 153 | > 172 | > 172 | > 172 |
| yield/Y_Q at H_i (diluted) | 5.89 | 6.14 | 6.60 | 7.32 | 7.66 | 6.83 | 4.92 | 3.21 |

A slower release gives **more charge and less debris**: |X| tracks the shrinking Hubble-held minimum (|X| = 46, 12 at
a = 1.6, 2.5 for Γ = 0.1 H against 7.0, 14.3 for the instant release), so the X⁶ torque acts longer at large |X|, while
the forced Y follows adiabatically (ω_Y/Γ_φ ≥ 230) and is hardly kicked. The window passes for every Γ_φ tried, with
the margin growing from 71 to > 5000; the normalising H_i moves from 1.33e11 to 2.0e11 GeV at Γ_φ = 0.1 H_i (yield
dilution 3.56). hw_gates saw the same sign when removing the X Hubble mass slowly from its D3 state (ε −0.20…−0.33).

**Stated point (Codex candidate 10 couplings, c_H = 10, δ = π/6), Codex preparation:** ε = 0.157, 0.144, 0.118, 0.106,
0.105, 0.064, −0.037, −0.057; free N_Y = 39, 50, 50, 140, 210, 164, 76, 22; R = 242, 335, 655, 623, 939, 1001, 704,
507 (R/R_L = 2.4, 3.9, 10, 11, 18, then R_L below the table's 29.8: > 17–40); diluted yield 5.7 → 0.34 Y_Q. Y-relaxed
preparation: the same pattern (ε 0.165 → −0.061, R 192 → 520–1200). (Instant reference against ls_scaling: ε 0.15749
vs 0.15741; R 242 vs 228 because N_σ at this resonant point is 29.4 in the 30–100 window against 25.0 at 100–300; the
6% window effect is immaterial for the Γ comparison.)

**Verdict.** The instant-release assumption matters **less** for the heavy window's pass/fail (it passes for all
Γ_φ/H_i from 0.1 to ∞, and the sensitivity is favourable: slower ⇒ |ε| up to 2.6×, R down to 0.026×) and **more** for
the published candidates (slower ⇒ charge down to 0.4× by Γ_φ = H_i and sign-reversed at ≤ 0.3 H_i, debris up 2–6×).
In this model of the release the instant assumption is conservative for the heavy window's debris verdict and
load-bearing for the stated point's already-failing charge.

## 7. Codex's question, within this track's evidence

"Does producing sufficient lasting charge necessarily generate excessive surviving (processed, late) debris?" Within
ls_scaling's declared domain (homogeneous, classical, supplied held phase, Y-relaxed preparation, decay ledger,
case-L proxy, no vector production): **no** — ε = −0.146 and 1.00 Y_Q with R = 1.14 against a direct R_L = 81.5, and
the finite-release scan gives ε = −0.38 with R = 0.03. The implication **does** hold in two places this track measured:
the resonant regime (m ≤ 12, and the stated couplings at every release duration), and any preparation that leaves
|Y − Y_*| > 0.37 H at the heavy point. Codex's point (i) — the forced Y carries the X⁶ torque with no excitation — is
realised dynamically (N_Y = 0.007 per comoving H_i³ for 5416 units of charge at Γ_φ = 0.1 H_i).

## 8. Verdicts

| # | Claim | Source value | This track (own code unless stated) | Verdict |
|---|---|---|---|---|
| 1 | ls_scaling ladder m = 10–40: ε, R, R_L, yield | −0.113…−0.145; 145, 34, 8.9, 1.14, 0.25; 60–88; 4.2–5.9 | −0.1129…−0.1454; 144.5, 34.2, 8.94, 1.144, 0.251; 59.8–87.8 (table); 4.22–5.89 [C, conformal Yoshida-8 vs their cosmic DP5] | CONFIRMED (≤ 0.3% in ε, ≤ 1% in R) |
| 2 | Spot point at the normalised scale: H 1.334e11, M_Y 4.00e12, h 5.33e12, μ 9.55e11 GeV, κ 3.46e-3, ε −0.146, yield 1.00, R 1.14, c6 1.25e6 | as stated | H 1.3335e11, M_Y 4.0005e12, h 5.3340e12, μ 9.5471e11, κ 3.4587e-3, ε −0.14567, yield 1.0000, R 1.144, c6 1.2504e6 [C] | CONFIRMED |
| 3 | Ledger components, decay-averaged ⟨a⟩ | R_Y 0.81, R_casc 0.34 (hw_envelope) | R_Y 0.805, R_Σ 1.2e-7, R_casc 0.339; N_YH 0.01725, N_YL 0.3345, N_σ 3.2e-9 [D, C] | CONFIRMED |
| 4 | R_L = 88 at the spot point (sp_critic table) | 88 [Q] | **81.5** at Y_Q (direct, bracket 80–82); 90.4 at the raw 5.89 Y_Q [C, our case-L code]; table 7.8% high (coarse bracket) | CORRECTED; margin 71 not 77 |
| 5 | Codex preparation at the spot: R ≈ 1414, N_YH 77 | 1414, 77 | 1415, 76.8 [C]; the 77.8 prepared free quanta are the emitted ones | CONFIRMED |
| 6 | Y-relaxed is "what a heavy Y does during inflation" (ls_scaling) | [A] | relaxation at 3H/2 (ω_Y = 23, 36 H underdamped); R < R_L after 1.02 held e-folds, converged by 3–4; pass needs |Y − Y_*| ≤ 0.37 H (f ≥ 0.75); but the phase drifts 0.25 rad/e-fold unheld and a Y-coupled holder changes the attractor [C] | CORRECTED (scope; conditional on a Y-blind holder) |
| 7 | Phases 15°–165°: R ≤ 1.5, best ε −0.154 (R_L 96) | as stated | 5°–175°: R 0.59–1.47, R/R_L ≤ 0.04; best −0.1546 at 70° (R_L 96.5); mirror exact [C]; Codex prep fails at every phase (R/R_L ≥ 1.4) | CONFIRMED (+ new) |
| 8 | Edge fraction 0.25–0.9: ε 0.074–0.29, R ≤ 1.4 | as stated | 0.074–0.294, R 0.86–1.35; extended: 1.0 → −0.35 (κ = 0.99 κ_max), 1.2 → −0.63 (κ > κ_max), ≥ 1.3 no inflationary minimum [C] | CONFIRMED (+ edge = κ_max at m = 30) |
| 9 | Exact scaling: ε, N, R invariant; yield ∝ s^{1.5}; c6 ∝ s^{−2} | machine precision at s = 0.3066 | s = 0.1, 0.3066, 3, 10: Δε ≤ 4e-13, ΔR/R ≤ 1e-10, yield 1.6e-12, c6 3e-16; conditions: constant g*, m_A ≫ m_X [D, C] | CONFIRMED |
| 10 | Heavy-limit slopes N_Y ∝ m^{−(4.2→4.8)}, R ∝ m^{−(4.9→5.3)}; window m ≈ 13–43 | as stated | 4.2, 4.6, 4.8; 4.7, 5.1, 5.3; debris edge 12.1–12.3, stability 42.8 [C] | CONFIRMED |
| 11 | Instant release (ls_scaling §5.4, open) | not computed | heavy: passes for Γ_φ/H_i = 0.1–∞, R 1.14 → 0.03, |ε| 0.146 → 0.38; stated: ε 0.157 → 0.06 (Γ = H), sign flip ≤ 0.3 H, R ×2–6 [C, one release model] | new: less sensitive (and favourably) in the heavy regime |
| 12 | Codex 2908e91 (ii): envelope table, r_i^{7/2}, extrapolation 11.87 | as stated | reproduced to the printed digits; extrapolation at M_Y/H_i = 5.96; ladder sits at 0.20–0.24 of the bound [C] | CONFIRMED |
| 13 | Codex 2908e91 (i): forced Y carries the torque with zero excitation | [A] | N_Y = 0.007 per H_i³ at 5416 units of charge (Γ_φ = 0.1 H_i) [C] | CONFIRMED dynamically |
| 14 | "Enough lasting charge ⇒ excessive processed late debris" | Codex's open question | fails inside ls_scaling's domain and for finite release durations; holds in the resonant regime and for preparations with |Y − Y_*| > 0.37 H [C] | REFUTED within the domain; NOT CHECKED outside it |
| 15 | Vector production, strings, running, light sector, holder designs | hw_gates | not re-done here | NOT CHECKED (see hw_gates) |

## 9. Limits

- Homogeneous, classical, tree level; no lattice, fragmentation, finite-k mediator dynamics or backreaction.
- The finite-release model is one choice: perturbative inflaton decay with the X Hubble mass tied to ρ_φ, instant
  thermalisation of the decay products, no inflaton oscillation or preheating effects; the reheating corrections to the
  yield and to ⟨a⟩ are applied as the factors (Rr/3)^{∓3/4, 1/4}, not re-derived through the decay epoch.
- The de Sitter stage holds the two phases by fiat (no holder term); hw_gates models the holders.
- The ledger and the case-L proxy are reused from our earlier tracks (the direct R_L is a reuse of our own benchmark
  code, not an independent threshold); R_L is undefined outside |ε| ∈ [0.068, 0.24] in our table.
- The stated-point window (30–100) leaves N_σ 6% unconverged; the heavy rows are converged (a_end = 300 control).
- Codex's envelope numbers are our implementation of Codex's quoted formula (reproduction of a formula, not an
  independent derivation of the bound).

## Credit

Yoshida (1990) symmetric compositions; Dormand–Prince 5(4); Dine–Randall–Thomas Hubble-induced terms and Affleck–Dine
[R]; Q-ball formation and fragmentation (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki) [R]; Codex 2908e91
and b448f25 [Q]; our tracks ls_scaling, sp_critic, sp_yledger, gp_leftovers, hw_envelope, hw_gates [Q where cited].

## Reproduction (~25 min on 4 cores; `hwrel.c` is compiled with gcc -O3 into the scratchpad build directory)

```
T=research/reality_framework/grut_formation_review_20261025/tracks/hw_verify
cd $T/code
python3 -B ladder.py  ../results/ladder.json            # ~4 min (18 releases)
python3 -B phases.py  ../results/phases.json            # ~3 min (55 releases)
python3 -B scaling.py ../results/scaling.json           # ~1 min
python3 -B prep.py    ../results/prep.json              # ~1 min (33 releases + de Sitter stages)
for g in heavy stated_codex stated_relaxed; do python3 -B release.py ../results/release_$g.json $g; done   # ~7 min
python3 -B rl_actual.py ../results/rl_actual.json 0.1457 1.0:70,80,85,90 5.889:50,60,70,80,90      # ~5 min
python3 -B rl_actual.py ../results/rl_actual_fine.json 0.1457 1.0:82,83 5.889:91,93                # ~4 min
python3 -B analytic.py ../results                        # 1 s  -> analytic.json
python3 -B summary.py  ../results                        # 1 s  -> summary.json (every number quoted above)
```
