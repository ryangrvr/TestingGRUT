# Skeptic's review of track fx_holder (inflaton-trilinear phase holder, λ_Y = λ_Σ = 0.5 repair)

1 November 2026 · Claude (adversarial reviewer, key fx_holder_skeptic), issue #13 working-theory repair pass (BRIEF.md of
this directory) · same-session AI work on one supplied model, not outside replication. Nothing here is evidence that X, Y, Σ
or an inflaton coupling exist; nothing is a validation or an exclusion. **The repair under review is a model choice**, and so
is every alternative boundary condition considered below. Labels: [D] derived, [C] computed with our own code, [E] estimate,
[Q] quoted, [A] argued.

**Independence.** Every number marked "skeptic" comes from code written for this review in `skeptic/code/`. None of the
track's `fx_common.py`, `rptri.c`, `tri_*.py` and none of Codex's code was imported or run; the track's JSON results are only
read for comparison (`sk_summary.py`). Method differences:
- `sk_common.py`, `sk_point.py`: own Wirtinger derivatives; stationary points by 40-digit `mpmath.findroot` on the five real
  gradient equations (no amplitude ansatz); Hessians by `mp.diff` + `mp.eigsy` (the track used sympy + BFGS/LM); multi-start BFGS
  with an own analytic gradient; joint (φ_I, X, Y, Σ) quartic positivity by rescaled direct minimisation.
- `sk_evolve.py`: scipy `solve_ivp` DOP853 (the track: hand-written DP5(4)) on the **prescribed** Einstein background
  H² = (J²/2 + b²I²/2 + ρ_r)/3, R = 2b²I² − J² (the track: constraint-solved Jordan engine with exact back-reaction). The holder
  is written **without g_I or η²**: Gφ̂(t) = (10i − μ)·I(t)/I(0), so u(t) = μ + (10i − μ) I(t)/I_ref. Own preparations, own
  window means on dt = 0.05 samples, own free-Y diagnostic.
- Reproduction: one core, about 13 minutes in total (33 runs to H_i t = 1000, two at a time).

## 1. Verdicts on the track's claims

| # | Track claim | Skeptic (method) | Verdict |
|---|---|---|---|
| FX1 | g_I = 8.3067×10⁻⁷ (b = 3), 1.3845×10⁻⁶ (b = 5); α_I = 2.192144074 = α(1A); g_Iφ_I(0) = 12.2987 H_i | 40-digit GeV computation, φ_I(0) = √6 H_i M_r/m_I: **8.306742069498×10⁻⁷**, 1.384457011583×10⁻⁶ (rel. diff 3×10⁻¹², η²-order terms), α_I identical, 12.298691 H_i | stands |
| FX2 | inflaton-mass corrections ≤ 2×10⁻⁹ H_i² | not recomputed beyond orders of magnitude; with an M_P (not 10¹⁵ GeV) cutoff the sunset is 8.8×10⁻³ H_i², still ≪ 9 [E] | stands |
| FX3 | β(λ_IΣ) = g_I²/4π², β(λ_IY) = 2×; at φ_I(0) δm_Σ² ≈ +1.2, **δm_Y² ≈ 7.7** H_i² per unit log; portals ≲ 4.6×10⁻¹⁵ "natural at one loop"; "phase-blind and below the mass gaps" | Tr M⁴ from numerical Hessians: coefficients exactly 8\|u\|²\|Σ\|² + 16\|u\|²\|Y\|², so the induced terms are (\|u\|² − μ²)(\|Σ\|² + 2\|Y\|²)/4π² per unit ln Q. **δm_Y² = 2.47 H_i²** per unit log (the track kept the φ-linear A_IY piece for Σ but dropped it for Y: 7.66 − 5.19). The yield is **not** insensitive to them (section 2: −1.7 % in N_X per unit ln Q at λ = 0.5; a Planck boundary condition moves H_i by −17 %). The bound 4.6×10⁻¹⁵ is exceeded by the track's own one-loop running after 0.26 units of ln Q (2.3×10⁻¹³ from M_P) | **REFUTED in part** (δm_Y² number; "natural"/harmless reading) |
| FX4 | φ_I → YY, ΣYY closed; 6X ~ 10⁻³⁸ H_i | thresholds re-derived: 44.50, 101.07 (λ = 1); **84.50 and σσ 80.0 at λ = 0.5**, all closed for b = 3, 5; 6X rate not recomputed | stands |
| FX5 | conjugate −N_X; opposite branch (125.94, 16.62, 29.69), δ = 2.532, N_X = −3.93× | own prep identical to 2×10⁻¹⁶; opposite run **N_X = 26118.6073** (track 26118.6073, 1.1×10⁻⁹), ratio −3.930682, N_YL 41.763. But it is **not a different minimum**: the extended scan (FX6) shows one smooth branch through s = 0 into s < 0 | stands; wording (see N7) |
| FX6 | prep (87.19556, 5.75425, 28.38718), Hessian ≥ 13.09; 7×7 shifts ≤ 10⁻⁹; unique and smooth for s ∈ [0, 2] | mpmath: (87.19556213287578, 5.75424796328559, 28.38718194436162), Hessian (13.091408829453, 14.58924997867, 538.14664342755, 1546.6597219731, 3239.8852511587), ≤ 2×10⁻¹⁴ from the track; 7×7 inflaton shift −7.96×10⁻¹⁰ (λ = 1), **−4.6×10⁻⁹ at λ = 0.5**; scan extended to **s ∈ [−1.5, 4]**: always the unique minimum (80–150 starts), lightest mass² ≥ 4.47 H_i² | stands, strengthened |
| FX7 | trilinear/original = 0.81840, 1A/original = 1.01907, H_i(tri) = 1.5350×10¹¹ | **0.8184043101** (track 0.8184043102), 1A **1.0190698**; N_X(tri) −6644.802238 vs −6644.802239 (1.5×10⁻¹⁰); Y/Y_Q 0.818251, H_i 1.53497×10¹¹; charge fraction at H_i t = 0.8: 0.490; min Im u = −2.63 | stands |
| FX8 | (5, 1): tri 1.06817, 1A 1.35428; λ = 0.5 tri 1.27145, 1A 1.60639 | 1.0681700, 1.3542848, 1.2714477, 1.6063936 (≤ 1.2×10⁻⁹) | stands |
| FX9 | holders make the yield history-dependent (15–20 % in H_i) | reproduced at fixed prepared state; **at fixed Lagrangian the trilinear moves the other way** (N3) | stands; incomplete |
| FX10 | κ/κ_max 0.834, λ_Y,eff 0.2585, margin 1058, chain −9.0×10⁵, unique prep | arithmetic re-derived; 150-start search finds only the prep at λ = 0.5 (V = −33293.717302) | stands (its "quartic ≥ 3.8×10⁻²⁴" is FX15) |
| FX11 | repaired λ = 0.5: prep (90.79987, 6.88773, 28.57624), lightest 10.19; Y/Y_Q 0.99073; **H_i = 1.3512×10¹¹ ± 0.13 %** | prep identical (2×10⁻¹⁶), lightest 10.191481317; N_X −8045.485138 (2.6×10⁻¹⁰); **Y/Y_Q = 0.990733, H_i = 1.35120×10¹¹**. The ±0.13 % is windowing only; the unstated renormalization condition of the induced inflaton portals moves H_i by **±1.1 % per unit ln Q**, and to **1.124×10¹¹** if they vanish at M_P (N1) | number stands; **error bar not established** |
| FX12 | debris N_YL 0.567 (λ = 0.5), 1A 12–24× more | own diagnostic: N_YL 0.56664, 0.25583, 3.0757, 41.763 (≤ 10⁻⁶ rel.); 1A/tri = 12.0, 12.1, 23.8, 16.2 (late N_Y std ratio 12.9–22.7); ledger not re-run | stands |
| FX13 | band: Y/Y_Q 1.556 … 0.818 | λ = 0.25: 1.556048 (H_i 1.00002×10¹¹), 0.75: 0.879656, unheld 0.25: 1.021627 (≤ 7×10⁻¹⁰) | stands |
| FX14 | d ln\|N_X\|/d ln g_I = 1.71 (tri), 1.34 (1A) | **1.71241** reproduced at λ = 1; but at the proposed λ = 0.5 it is **3.342 (tri)** and **2.648 (1A)** (N2) | stands at λ = 1; does not apply to the working point |
| FX15 | λ_I ≥ g_I⁴/(λ_Σλ_Y²) = 4.8×10⁻²⁵ / 3.8×10⁻²⁴; runaway \|Σ\| > 2.55×10⁶ H_i | the X sector lowers the effective Y quartic: the exact threshold (λ_I/4 φ⁴ normalisation, worst-case phases, numeric = analytic) is g_I⁴/(λ_Σλ_Y,eff²) = **8.28×10⁻²⁵ (λ = 1), 1.42×10⁻²³ (λ = 0.5)**; runaway at **2.22×10⁶ (λ = 1), 1.30×10⁶ H_i (λ = 0.5)**; the one-loop g_I⁴/16π² = 3×10⁻²⁷ cannot supply it | **REFUTED** (factor 1.7–3.7; physically inconsequential) |
| FX16 | improvement in control; costs 1–5; "1A worse on" history; working point H_i = 1.351×10¹¹ | costs 4 and 5 understated (N1, N2); at λ = 0.5 1A and the trilinear have comparable history dependence (Δln H = −0.207 vs −0.176 at fixed prep; −0.207 vs +0.193 at fixed Lagrangian) and the trilinear is the **more** clock- and portal-sensitive one; no failure found, so no fallback forced | stands as a model choice, with corrections |
| FX17 | finite-k Y production and exchange freeze-out OPEN | not computed here | OPEN |

Two codes that share no line agree on every N_X to ≤ 1.5×10⁻⁹ (median 2×10⁻¹⁰). The residual is the size of the omitted
O(η²V) ≈ 10⁻¹⁰ hidden energy in the prescribed background. The track's rtol control (≤ 1.2×10⁻¹²) therefore understates the
cross-code spread but is immaterial. Our convention reproduces the M1 calibration: the unheld λ = 1 run gives
Y/Y_Q = 1.0000000004 in the 640–800 window. Option 1A's long-window +1.907 % is confirmed.

## 2. New issues

**N1. The repaired H_i depends on an unstated renormalization condition of the inflaton portals [C].** The trilinear's own
one-loop running, which the track computed (FX3), induces exactly the phase-blind operators
c(Q)(|u|² − μ²)(|Σ|² + 2|Y|²), with dc/d ln Q = 1/4π². These are λ_IΣφ²|Σ|², A_IΣφ|Σ|², λ_IY, A_IY in fixed ratios. Their value
at the hidden scale is a new input. The track set it to zero at Q = M_Y without saying so. Re-solved preparation plus full
evolution (`sk_evolve.py`, cport):

| run | c(M_Y) in units of 1/4π² | N_X / N_X(c = 0) | Y/Y_Q | H_i (10¹¹ GeV) |
|---|---|---|---|---|
| trilinear λ = 0.5 (repaired point) | +1 / −1 | 0.98341 / 1.01740 | 0.97430 / 1.00797 | 1.3664 / 1.3358 |
| | **−13.31 (zero at M_P)** | **1.31827** | **1.30606** | **1.1239** (N_YL 1.41) |
| | +13.31 | 0.82995 | 0.82227 | 1.5300 |
| trilinear λ = 1 | +1 / −1 | 0.99190 / 1.00828 | — | ±0.55 % per unit |
| | −13.31 | 1.12671 | 0.92193 | 1.4176 |
| Option 1A λ = 0.5 / 1 (same structure in u(R)) | +1 | 0.98825 / 0.99468 | — | — |

The sensitivity is d ln N_X/d(unit ln Q) = −1.70 % at λ = 0.5, so d ln H_i = +1.13 % per unit. With the minimal ultraviolet
condition (portals vanish at M_P and run to M_Y with the track's own β functions) the repaired point moves to
**H_i = 1.124×10¹¹ GeV, −16.8 %**. The finite one-loop CW part at Q = M_Y has the same form with O(1) logarithms, so even the
most favourable choice leaves about ±1 % [E]. The honest statement is that the trilinear brings four renormalizable portal
couplings whose value at M_Y is a supplied input. It must be fixed to about 0.9 units of the one-loop running (≈ 1.5×10⁻¹⁴ in λ_IΣ) for
1 % in H_i. That is a tuning against the coupling's natural size, not "natural at one loop". 1A has the same structure with
2/3 of the sensitivity, so this is not a trilinear-only cost. The track's FX3 "phase-blind, below the mass gaps" (read as
harmless) and FX11's ±0.13 % do not hold.

**N2. The clock sensitivity at the working point is twice the quoted value [C].** Re-prepared ±1 % runs give:

| | d ln\|N_X\|/d ln(holder coupling) |
|---|---|
| trilinear λ = 1 | 1.7124 (reproduces the track) |
| **trilinear λ = 0.5** | **3.342** |
| Option 1A λ = 0.5 | 2.648 |
| Option 1A λ = 1 | 1.337 (track) |

At the proposed working point H_i ∝ φ_I(slice)^(−2.23): a 1 % uncertainty in the inflaton's field value at the end of
inflation is a 2.2 % uncertainty in H_i. The trilinear is the more clock-sensitive holder at both couplings.

**N3. History dependence at fixed Lagrangian [C].** The track compares (3, 0.5) and (5, 1) at a fixed prepared state
u(0) = 10i, re-tuning g_I to the new φ_I(0) = √6 H_i M_r/m_I. At a fixed Lagrangian (g_I set once by the (3, 0.5) history), the
(5, 1) inflaton prepares u(0) = 2.864 + 6.000i (δ = −2.016, not −π/2):

| λ | H_i: (3, 0.5) | (5, 1) re-tuned g_I | (5, 1) fixed g_I | 1A (5, 1) at fixed ξ |
|---|---|---|---|---|
| 1 | 1.5350 | 1.2842 (−16.3 %) | **1.6621 (+8.3 %)** | 1.0963 vs 1.3262 (−17.3 %) |
| 0.5 | 1.3512 | 1.1329 (−16.2 %) | **1.6384 (+21.3 %)** | 0.9694 vs 1.1924 (−18.7 %) |

For the trilinear, g_I must know the inflaton mass (g_I ∝ m_I at fixed u(0)), just as 1A's ξ must know H_i. The Lagrangian-level
history dependence of the two holders is comparable at λ = 0.5. The unheld run differs by only 0.08 % across histories.

**N4. δm_Y² per unit log is 2.47 H_i², not 7.7 [C].** Covered in FX3; the track overstated this cost.

**N5. Joint boundedness threshold [D, C].** Covered in FX15: use λ_Y,eff = λ_Y − 27κ⁴/(16λ_X³), giving 1.42×10⁻²³ at the
repaired point. It must be a tree-level inflaton quartic; any real one suffices; nothing changes dynamically.

**N6. Forced-Y diagnostic in rpholder.c and rptri.c [D, C].** Both engines use yp = (−A J + B J*)/D with B = 2μ_eff Σ and the
rotation e^{−i arg B/2}. The Y equation, ∂V/∂Y* = M²Y + κX³ + 2μ_eff* Σ Y*, needs c = 2 conj(μ_eff) Σ. They agree only for real
μ_eff. At the prepared point (μ_eff = 10i) the diagnostic's forced Y is −1.562 against the correct −6.901 (actual −5.754), and
the light and heavy directions are swapped. Early-time N_YL/N_YH in both tracks are therefore wrong while Im μ_eff ≠ 0
(H_i t ≲ 20–50). No reported number is affected: all are 400–1000 means with μ_eff real to 10⁻⁶. Our diagnostic, built with the
correct c, reproduces every late N_YL to ≤ 10⁻⁶. It is a latent bug for any future early-time use.

**N7. The opposite branch is the same branch [C].** The prepared family is one smooth, unique minimum for s = φ_I/φ_I(0) from
−1.5 to 4 at fixed H. |u| runs from 5.8 to 45.4 H_i and the lightest mass² stays ≥ 4.47 H_i². The forced light-Y mass² at the vev
crosses zero near s ≈ 1.5 and s ≈ −0.85: a crossover into the Y-heavy regime, not a separate phase. The "Y-condensed opposite
minimum" (FX5) is s = −1 of this family. There is no first-order trap if φ_I was larger earlier in inflation (as in any
large-field history), which strengthens the track's adiabatic-tracking argument [A].

**N8. The reflection-even inflaton sector is a real assumption [A].** The CP-sign statement "α_I ≃ α_I + π" needs the
inflaton's decay to respect φ_I → −φ_I. The common scalar channel A φ_I|H|² is odd, and Γ_I = 0.5 H_i needs A = (8π m_I Γ_I)^{1/2}
≈ 6.1 H_i [D]. With it the two branches are not degenerate: radiative odd terms appear, and at φ_I ≈ 0.8 M_r the Higgs mass² is
±Aφ_I ≈ ±9×10⁷ H_i², so one branch has a tachyonic Higgs. A fermionic decay with a chiral rotation can keep the reflection. Either
way, the trilinear makes the inflaton's decay sector part of the CP and stability bookkeeping, which it was not for 1A.

## 3. What this review did not do

- No CW-corrected preparation or evolution: N1 uses the running direction only, and the finite part is an [E].
- The debris ledger (R_central) was not re-run; only its inputs were reproduced.
- No 1000–2000 windows (our runs end at H_i t = 1000).
- No finite-k Y production (FX17 stays OPEN).
- The 6X width and the two-loop sunset were not recomputed.
- Nothing here is evidence for or against the real world.

## Reproduction (one core per process, two processes, about 13 min)

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_holder/skeptic/code
python3 -B sk_point.py                       # static items                 -> ../results/sk_point.json   (25 s)
# runs: TAG MODE LAM B GAM SCALE CPORT END [BREF]; CP = 1/(4 pi^2) = 0.0253303, CPP = ln(M_P/M_Y)/(4 pi^2) = 0.3371983
python3 -B sk_evolve.py tri_l05_b3 tri 0.5 3 0.5 1 0 1000          # repaired point            -> ../results/runs/*.json (~40 s each)
python3 -B sk_evolve.py tri_l05_b3_cpPlanck tri 0.5 3 0.5 1 -0.3371983231 1000
python3 -B sk_evolve.py tri_l05_b5_fixg tri 0.5 5 1 1 0 1000 3      # (5, 1) at fixed g_I
xargs -P 2 -L 1 python3 -B sk_evolve.py < jobs.txt   # all 33 runs (orig/1A/tri; lam 1, 0.5, 0.25, 0.75; b 3, 5; cport 0, +-1, +-13.31 units; scale 0.99/1.01)
python3 -B sk_summary.py                     # comparison with the track's JSON -> ../results/sk_summary.json
python3 -B sk_diag_check.py                  # N6                            -> ../results/sk_diag_check.json
```
