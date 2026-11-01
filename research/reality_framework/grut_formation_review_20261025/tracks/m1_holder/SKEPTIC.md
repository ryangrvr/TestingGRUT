# Skeptic's review of track m1_holder (Codex f8186fa Option 1A holder; 94f660d sextic)

1 November 2026 · Claude (adversarial reviewer, key m1_holder_skeptic), issue #13, BRIEF_M1.md · same-session AI work, not
outside replication. Codex packets read only through `git show f8186fa:…`, `94f660d:…`, `d4747cb:…`; nothing checked out,
nothing committed. Remote re-checked on 2026-11-01 (`git ls-remote origin`, `git log --all --since=2026-10-29`): the newest
`codex/grut-m1-*` head is f8186fa; the brief's "newer isocurvature packet" (4e-12 / 3.02 / 3.435e-70 / 38 checks) does not
exist on the remote — the track's NOT CHECKED stands.

**Independence.** Every number marked "skeptic" below comes from code in `skeptic/code/` written for this review, using
methods different from the track's and from Codex's: `skeptic_point.py` (40-digit mpmath Newton on the three amplitude
equations, Hessian by `mp.diff` and `mp.eigsy`, symbolic re-derivation of every inequality of the bound with sympy,
numerical rephasing check of the CP phase, UV-completion and effective-scale estimates), `skeptic_evolve.py` (third
implementation of the homogeneous reheating system: scipy `solve_ivp` DOP853 on Codex's *prescribed* background written from
DERIVATIONS.md, own right-hand side, dt = 0.05 trapezoid windows, running 80-wide windows, charges of all three fields),
`skeptic_power.py` (closed-form small-x expansion of the de Sitter modes and an ODE in the conformal variable x with WKB start;
neither mpmath Hankel functions nor the N-variable integration), `skeptic_extra.py` (own re-solve of the 94f660d sextic point,
effective scales, exchange summary). Codex's and the track's numbers are quoted [Q] for comparison; none of their code was run.
Labels [D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued.

## 1. Verdicts on the track's 20 claims

| # | Track claim (verdict) | Codex [Q] | Track | Skeptic (method) | Skeptic verdict |
|---|---|---|---|---|---|
| 1 | (ξ, α) ⇔ μ_eff(R), dimension 5 (CONFIRMED) | 30.746728738599167, 2.192144074179184 | identical | **30.74672873859916466, 2.19214407417918401** (40-digit mpmath from the GeV pins; lock = (−4.4430e-12 + 6.2057e-12 i) GeV⁻¹); with the brief's μ = 7.2 one would get 30.806, 2.1948 | stands |
| 2 | ξ ≈ 31 not loop-natural, needs tree input (NEW [E]) | "dimension five, not UV-complete" | loop 1.5e-3·c, required/loop 2e4 | loop estimate confirmed (with the heat-kernel 1/6: 2.5e-4, ratio 1.2e5). **Sharper: the operator's own scale is Λ_eff = M_Y/ξ = 0.976 H_i = 1.31e11 GeV**, below every amplitude (|Y|/Λ_eff = 5.9, |Σ|/Λ_eff = 29, |X|/Λ_eff = 89) and below √R (R/Λ_eff² = 12.6); relative to the declared cutoff 1e15 GeV the coefficient is c_Λ = 7632. The abandoned sextic had Λ/ξ₆^{1/5} = 142 H_i | stands, strengthened |
| 3 | EFT control at the prepared point, "CONFIRMED for Y (margin 5), caveat X" | not stated | |Y|/M_Y = 0.19 | **overstated.** The margin is measured against M_Y, but the only known coefficient in the tower is ξ = 31, i.e. Λ_eff ≈ H_i. Any tree-level singlet completion with g₁g₂ = ξM_S/M_Y induces −(2Re g₂Σ*Y²)²/(2M_S²), which at the point is −1693/g₁² × λ_Y|Y|⁴ (independent of M_S); even at g₁ = 60 it is 47 % of λ_Y|Y|⁴ (−516 H_i⁴) | **CORRECTED → "not established" (no EFT margin relative to the operator's own scale)** |
| 4 | α is the supplied rephasing-invariant CP phase; no selection law (CONFIRMED) | stated | argument | numerical rephasing check: after making κ, μ real by U(1)³ rephasings, arg(lock) − α = 0 for five random rephasings | stands |
| 5 | Stationary point, δ = −π/2 exactly (CONFIRMED) | (87.19556213287578, 5.75424796328559, 28.387181944361622) | same to 3e-15 | **87.195562132875782802, 5.7542479632855904690, 28.387181944361621657** (mpmath findroot, 40 digits, residuals 1e-37; rel. diff to Codex 5e-17); V = −32480.585311406205 (rel. diff 9e-18) | stands |
| 6 | Hessian (13.0914, 14.5892, 538.147, 1546.66, 3239.89) + gauge zero (CONFIRMED) | as stated | sympy exact, 5e-14 | mp.diff + eigsy: **13.091408829453044, 14.589249978667072, 538.146643427554, 1546.659721973108, 3239.885251158698** (rel. diff to Codex ≤ 1.7e-14); gauge null residual 3e-39; all positive | stands |
| 7 | V ≥ −285370 H_i⁴ at R = 12 (CONFIRMED, loose 8.8×) | −285369.6842985777, B = 0.5283522987 | chain re-derived | every inequality re-derived symbolically (min_x[c x⁴ − 2κyx³] = −27κ⁴y⁴/(16c³) etc.); bound −285369.6842985776 (3.5e-16); of the reserves 0.05/0.1/0.2/0.3 the chosen 0.2 is the tightest (0.1: −375217, 0.3: −382175; ≥ 0.4 gives B < 0). **Fragility [D]: the proof needs (1 − reserve)λ_X > 6.23e-4, a margin of 1.6 in λ_X (1.42 in κ)**; with the completion-induced δλ_Y of claim 3, B > 0 needs g₁ > 56.6 | stands, with fragility noted |
| 8 | Effective Planck ratio 1 + 4.81e-11 (CONFIRMED) | 1.0000000000480667 | same | dV/dR = −7902.4525029, ratio 1.0000000000480667 | stands |
| 9 | Sextic holder 94f660d: δ = −2.09993, lightest 10.669 (CONFIRMED) | (86.657, 3.823, 28.338) | (86.6572, 3.8228, 28.3376), −2.099927, 10.66925 | own BFGS+LM on an own V with FD gradient: **(86.65719, 3.82277, 28.33761), δ = −2.0999272, lightest 10.66924** (FD Hessian) | stands |
| 10 | Late charge −8267.437 vs −8120.740, +1.80645 % (CONFIRMED) | −8267.43703208 / −8120.73969565, 1.0180645288 | −8267.43781046 / −8120.74008245, 1.0180645762 | third code: **−8267.43780895 / −8120.74008174, ratio 1.0180645761** (agrees with the track to 1.8e-10 / 8.7e-11; Codex's values are 9.4e-8 / 4.8e-8 off — Codex's rtol 2e-9 and 0.49-spaced trapezoid; harmless) | stands; Codex's digits beyond the 7th are not reproducible |
| 11 | Window dependence 1.5e-3, caused by 3.08 light-Y quanta exchanging charge with X (CORRECTED) | 3.7e-5 (holder N_X alone) | ratio 1.018065 / 1.018293 / 1.019612 (640–800 / 720–800 / 800–1000) | ratios **1.0180645761 / 1.0182929463 / 1.0196123997** (= track to 1e-9), 900–1000: 1.0195239; running 80-wide windows 400–1000: **1.01749 … 1.02082** (±0.17 %); N_X peak-to-peak 0.229 % (holder), 0.104 % (original). **Cause refuted:** corr(N_X, 6N_Σ) = −0.9999998 (holder), −0.999999999999 (original); corr(N_X, 3N_Y) = −0.01 / +0.23; N_X + 6N_Σ constant to 5e-3 (holder) and 9e-6 (original) across all windows while 3N_Y moves by the same tiny amounts. The partner is **Σ** (charge 6), not Y; the original run has the same modulation at half amplitude with 450× fewer light-Y quanta | **magnitude confirmed (larger); mechanism REFUTED → X ↔ Σ exchange** |
| 12 | H_i refit 0.98813541 → 1.32691046e11 GeV, ±0.1 % (CONFIRMED) | 0.9881354099844972 | 0.9881353793 | 0.9881353794 from the 640–800 ratio; the 80-window spread gives **0.98636 … 0.98850, i.e. H_i(new) = 1.3245e11 … 1.3274e11 GeV (±0.11 %)**; the fully-cooled 0.98263384 reproduced | stands (the ±0.1 % is slightly understated: ±0.11–0.17 %) |
| 13 | Σ-phase slope 0.554384 (CONFIRMED, 4e-6) | 0.5543843049 tangent; FD 0.5543843711 / 0.5543843215 | 0.5543819906 / 0.5543819413 | **0.5543819908 / 0.5543819413** (steps 1e-3 / 5e-4): identical to the track to 2e-10; Codex is 4.3e-6 high in both FD and tangent | stands; the 4e-6 is Codex's error, so "0.554384305" is CORRECTED to 0.554382 |
| 14 | Response row agrees with Codex to 3.9e-6 (CONFIRMED) | 12-row | own FD row | not recomputed (20 pairs × 2 runs); consistency only: both rows reproduce the power table under skeptic mode functions (claim 17) | NOT CHECKED here |
| 15 | finite-k 7.33e-5 / 7.33e-7 | as stated | not run | not run | NOT CHECKED |
| 16 | Time shift = adiabatic mode (CONFIRMED [D]) | 4e-12 / 3.02 unpushed | theorem + −9e-7 | argument re-read: δy = εẏ solves the tangent system including −V_rR δR; agreed; no recomputation. Note the slow X ↔ Σ exchange (claim 11) is what makes d ln N_X/dt ≠ 0 at τ ≈ 760 (−9e-7 is of the right size for the running-window slopes) | stands |
| 17 | Power 1.23575e-10 (4), 7.59235e-16 (8), 1.76119e-31 (20), 1.44301e-70 (50) (CONFIRMED) | as stated | mpmath Hankel 1.4e-9 | **x-ODE (conformal variable, WKB start x₀ = 2000): Codex table reproduced to 1.4e-9 at N = 0, 2, 4, 8**, track's own-row table to 1.6e-11; **two-branch closed form reproduces N = 8, 20, 50, 60 to 3e-9** (Codex row) and 1e-14 (track row); leading law P = e^{−3N}/(4π²) Σ[(L_pos − 1.5L_vel)² + μ²L_vel²]/(μ(1 − e^{−2πμ})) is within 3.7e-5 of the exact values. The relative oscillation of |S|² is 2e^{−πμ} = 6.4e-5 / 3.2e-5 (the track wrote e^{−πμ}). 82 % of the late power sits in the 14.59 mode, 39 % comes from velocity response | stands (factor-2 wording fix) |
| 18 | 3.435e-70 (unpushed) NOT CHECKED | 3.435e-70 | 2.38× f8186fa | 2.3804× the closed form; not derivable from f8186fa data | stands |
| 19 | "No resolved signal" is the generic heavy-spectator expectation (CONFIRMED as expectation [A]) | stated | items (i)–(vi) | agreed; **add (vii): the X ↔ Σ (longitudinal-vector) reservoir split is itself history dependent at the 0.2 % level at τ = 800 (claim 11), so even the homogeneous "charge in X" input to any transfer is not frozen** | stands, one item added |
| 20 | Holder excites Y 450×; R_central ≈ 10 under the ledger scaling (NEW [E]) | not reported | N_YL 0.0068 → 3.08 | excitation confirmed independently: late fast N_Y oscillation std 0.0015 → 0.83 (×540); **its role as the channel of the slow charge drift is refuted** (claim 11); the R_central ≈ 10 [E] is not checked here | excitation stands; exchange statement REFUTED; debris estimate NOT CHECKED |

## 2. The slow charge exchange (new mechanism) [C, E]

Skeptic runs to H_i t = 1000 (rtol 1e-10, atol 1e-12, 0.68e6 and 1.12e6 rhs calls). Running means over 80-wide windows:

| window start | original N_X | original 6N_Σ | holder N_X | holder 6N_Σ |
|---|---|---|---|---|
| 400 | −8122.57 | 8122.57 | −8264.67 | 8264.67 |
| 520 | −8115.76 | 8115.76 | −8283.59 | 8283.59 |
| 640 | −8122.87 | 8122.87 | −8267.75 | 8267.75 |
| 760 | −8116.41 | 8116.41 | −8270.98 | 8270.98 |
| 880 | −8119.41 | 8119.41 | −8279.97 | 8279.97 |

N_X + 6N_Σ is constant to 5e-3 (holder) and 9e-6 (original); the Y charge takes no part (3N_Y moves by the same tiny
amounts; its fast oscillation of amplitude 0.83 in the holder run averages out). Extrema of the holder's N_X at window centres
560, 720, 920 (original 520, 660, 840): modulation period ≈ 300–360 H_i t, i.e. 5–6 late X orbits (quartic-oscillator period
≈ 63 H_i t at |X|_max = 1.85 H_i; original 71 at 1.65). The channel is structural: Σ couples to the rest only through
u Σ*Y², and integrating out the heavy Y gives the tree-level sextic (κ²u/M_Y⁴) Σ*X⁶ + h.c. with coefficient 1.06e-10 H_i⁻⁵ at
R = 0 — half the c₆ = 2.10e-10 of the 94f660d sextic holder, and with the phase of u(R). Its torque on Σ,
2a³ Im(Σ*uY²) ∝ |X|⁶ sin(6θ_X − θ_Σ), is sharply peaked at apoapsis of the eccentric X orbit (ε = −0.15), where the X phase is
nearly stationary, so each orbit delivers a coherent kick whose sign rotates with 6θ_X; the kick per orbit falls as a⁻² and
the cumulative remaining drift as τ^{−1/2}, so the exchange has not frozen at τ = 800 and was not followed to freeze-out.

Consequences: (a) Codex's "+1.80645 %" and the track's "+1.806458 %" are the 640–800 window values of a quantity that varies by
±0.17 % between windows at this epoch; the H_i refit is 1.3269e11 ± 0.0015e11 GeV (Codex's nine digits and the track's eight
are not meaningful beyond four). (b) The reheating "late charge in X" is not a frozen quantity in this model until the X ↔ Σ
(longitudinal-vector) exchange decouples; the formation input N_X/s inherits this. (c) The track's attribution to light-Y
quanta must be removed from NOTES §3 and §5 (claims 11, 20); the Y excitation itself and the debris [E] are separate statements.

## 3. Other corrections and new issues

1. **Effective scale of Option 1A [E].** (ξ/M_Y) R Σ*Y² = R Σ*Y²/Λ_eff with Λ_eff = 0.976 H_i. Every amplitude and √R exceed
   Λ_eff; a dimension-6 term R Σ*Y²|Y|²/Λ_eff² or R²Σ*Y²/Λ_eff³ with an O(1) coefficient would be 35× or 12.6× the kept term.
   Codex's "dimension five, not UV-complete" and the track's "not loop-natural" both understate this; the right statement is that
   the holder is a Hubble-scale operator whose tower is not under control at the prepared point.
2. **Tree completion versus the bound [E].** L ⊃ −½M_S²S² + S[g₁M_S R + (g₂Σ*Y² + h.c.)] gives ξe^{iα}/M_Y = −g₁g₂/M_S and the
   dimension-6 term −(2Re g₂Σ*Y²)²/(2M_S²) = −(1693/g₁²) λ_Y|Y|⁴ at the point for any M_S. B = λ_Y − 0.4717 > 0 requires
   g₁ > 56.6 (a large linear nonminimal coupling; its de Sitter tadpole ⟨S⟩ = 12g₁/M_S ≈ 24 H_i for M_S = M_Y); with
   g₁ ≤ 4π the induced quartic is −10.7 λ_Y and the fixed-R boundedness proof of DERIVATIONS.md no longer applies. A completion
   that avoids this is not excluded, but none is given.
3. **Bound-proof margin [D].** 27κ⁴/(16c³) < λ_Y needs c = (1 − reserve)λ_X > 6.23e-4: a factor 1.6 in λ_X (1.42 in κ).
4. **Codex's numerical precision [C].** Two independent codes (track rpholder.c, skeptic solve_ivp) agree on N_X to 1e-10 and
   on the phase slope to 2e-10; Codex's N_X values are 5–9e-8 off and its slope 4.3e-6 off (rtol 2e-9 background, 0.494-spaced
   trapezoid). Codex's "0.554384305" is 0.554382; immaterial for physics, but the quoted digits are not reproducible.
5. **R_min.** The track's "−1.0636 vs Codex −1.0633" is not a trajectory difference: Codex's `Ricci_range` is the minimum over
   dt = 0.05 samples (skeptic sampled value −1.06327166535, identical to 13 digits); the true minimum is −1.0636058 at τ = 0.70357.
6. **Wording.** The residual oscillation of |S|² for the two light modes is 2e^{−πμ} = 6.4e-5 and 3.2e-5, not e^{−πμ}.
7. **Brief.** μ = 9.614e11 GeV / 1.342842735e11 GeV = 7.159456 H_i; "7.2" would shift (ξ, α) to (30.806, 2.1948). Track's correction stands.

## 4. What this review did not do

No lattice or finite-k runs; the response row (claim 14) and the time-shift restarts (claim 16) were not recomputed; the
debris-ledger estimate R_central ≈ 10 (claim 20) was not re-run; the exchange was not followed beyond τ = 1000, so the
freeze-out value of the X ↔ Σ split is unknown; the completion of item 2 is one example, not a theorem; nothing here is evidence
that X, Y, Σ exist, and nothing is a validation or exclusion.

## Reproduction (one core, ≈ 3 min)

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd research/reality_framework/grut_formation_review_20261025/tracks/m1_holder/skeptic/code
python3 -B skeptic_point.py                      # -> ../results/skeptic_point.json   (1 s)
python3 -B skeptic_evolve.py --end 1000 --phase  # -> ../results/skeptic_evolve.json  (150 s)
python3 -B skeptic_power.py                      # -> ../results/skeptic_power.json   (11 s)
python3 -B skeptic_extra.py                      # -> ../results/skeptic_extra.json   (1 s)
```
Credit: Bunch–Davies/Hankel small-argument asymptotics (standard; e.g. Birrell–Davies); the adiabatic-mode theorem as cited by
the track (Weinberg 2003; Wands–Malik–Lyth–Liddle 2000); Codex's bound chain (DERIVATIONS.md) re-derived, not reused.
