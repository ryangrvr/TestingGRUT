# Track fx_holder: a renormalizable phase holder (inflaton trilinear), perturbative heavy couplings, recalibrated yield

1 November 2026 · Claude, issue #13 working-theory repair pass (BRIEF.md of this directory, repair 1) · same-session AI
work on one supplied model, not outside replication. Nothing here is evidence that X, Y, Σ or the inflaton coupling exist,
and nothing is a validation or an exclusion. **The repair is a model choice.** It does not follow from closure or from
anything else in GRUT.

Labels: [D] derived, [C] computed with our own code, [E] estimate, [Q] quoted, [A] argued, [R] recalled.

**Own code** (`code/`):
- `rptri.c` is a copy of our `m1_holder/code/rpholder.c` (itself our `rp_reheating/code/rpreheat.c`). It has
  μ_eff = μ + lock·R + G·φ̂_I with φ̂_I = φ_I/H_i = I/√η², and the exact inflaton back-reaction
  √η² ∂V_h/∂φ̂ = √η² · 2Re(GΣ*Y²) in the inflaton equation. The interaction energy enters ρ, the trace equation and the
  constraint exactly. Other changes: work ledger term (∂V_h/∂φ̂)·φ̂̇, a self-consistent constraint start, 20-wide
  time bins, and columns for u, φ̂ and the forced light-Y mass².
- `fx_common.py` holds the model and drivers. It copies the held-phase preparation of `rp_reheating/code/rp.py::prepare`
  and the exact-Hessian stationary-point method of `m1_holder/code/holder_point.py`, both generalised to any
  λ_Y = λ_Σ. The yield convention and decay ledger are those of rp_reheating.
- `tri_point.py` covers items 1, 2 and 5 plus the static parts of item 4. `tri_runs.py` covers items 3 and 4.
  `lam_scan.py` runs the coupling band. `controls.py` runs the rtol and independent-binary checks, and `summary.py`
  collects the headline numbers.

None of Codex's code was run.

The candidate is L ⊃ −(g_I e^{iα_I} φ_I Σ*Y² + h.c.), which gives μ_eff(t) = μ + G φ̂_I(t) with G = g_I e^{iα_I}. It is
fixed by μ_eff = 10i H_i at the preparation slice, which is the same requirement as Option 1A. The rest of M1 is
unchanged (H_i = 1.342842735×10¹¹ GeV, M_Y = 30, h = 40, μ = 7.159456 H_i, κ = 3.4587×10⁻³, λ_X = 10⁻³, ξ_X = −5/6).

## 1. Coupling, dimension, naturalness, decays, CP [D, C, E]

**Coupling [D, C].** The engine fixes the inflaton at the slice from the Friedmann constraint (a = 1, H = H_i, at rest),
giving I₀ = √6/b M_r:
- (3, 0.5): φ_I(0) = 0.81650 M_r = 1.988×10¹⁸ GeV, so g_I = |10i − μ|/φ̂_I(0) = **8.3067×10⁻⁷**;
- (5, 1): φ_I(0) = 0.48990 M_r, so g_I = 1.3845×10⁻⁶.

**α_I = arg(10i − μ) = 2.192144074 rad.** This equals Option 1A's α, because R(0) = 12H² > 0 and φ_I(0) > 0 enter the
same way. g_I φ_I(0) = 12.2987 H_i, which is |lock|·12 of 1A. The operator φ_I Σ*Y² has dimension 4 and g_I is
dimensionless, so there is no EFT cutoff. (1A's own scale was M_Y/ξ = 0.976 H_i.)

**Inflaton potential [C].** All corrections are far below m_I² = 9 H_i²:

| Correction | Value (H_i²) |
|---|---|
| Hidden-sector one-loop CW (MS-bar, Q = M_Y) at the late vacuum: δm_I² | −8.7×10⁻¹³ (analytic Y-loop check agrees to 10⁻¹⁰) |
| The same, running | 2.8×10⁻¹¹ per unit ln Q |
| Brief's g_I²M_Y²/16π² | 3.9×10⁻¹² |
| Two-loop sunset with a 10¹⁵ GeV cutoff | 1.5×10⁻⁹ |
| Tree level, hidden fields integrated out at their φ_I-dependent minimum | −4×10⁻¹⁰ |

The tadpole shifts the inflaton by δφ̂ = −2.8×10⁻⁷ H_i. g_I is the only coupling that breaks the inflaton reflection
φ_I → −φ_I, so its small value is 't Hooft natural.

**Hidden potential [C].** The induced portal operators run per unit ln Q as follows. The Tr M⁴ coefficients
reproduce 8g_I² and 16 Re(μG*) exactly.
- λ_IΣ φ_I²|Σ|²: g_I²/4π² = 1.75×10⁻¹⁴.
- λ_IY φ_I²|Y|²: 3.5×10⁻¹⁴.
- A_IΣ φ_I|Σ|²: Re(μG*)/2π² = −1.75×10⁻⁷ H_i.

At φ_I(0) these give δm_Σ² ≈ +3.8 − 2.6 = +1.2 H_i² and δm_Y² ≈ 7.7 H_i² per unit log. Over ln(M_P/M_Y) = 13.3 that is
at most about 16 and 100 H_i². These are phase-blind (U(1) invariant) and absent from X at one loop. They are small
against the Y and Σ gaps (M_Y² = 900, m_Σ² = 1600–3200), but comparable to m_I² and to 10–60% of the holder's own
tree-level effect on the light Y (495 → 334 H_i²).

The price is **sequestering**. Every tree-level inflaton–hidden portal must satisfy λ_IΣ, λ_IY, λ_IX ≲ 1/φ̂_I(0)² =
4.6×10⁻¹⁵ per H_i² of effect. This is natural at one loop, but it is a supplied zero, like M1's "zero portals".
Planck-suppressed c φ_I²Σ*Y²/M_P must have c ≲ 10⁻⁸ for 1% holder accuracy. That needs an approximate inflaton shift
symmetry [A]. The near-Planckian φ_I(0) is the price of renormalizability.

**Boundedness with the inflaton [D, C].** With no inflaton quartic, minimising over φ_I gives an effective sextic
−2g_I²|Σ|²|Y|⁴/m_I². The ½m²φ² surrogate alone then runs away for |Σ| > b√(λ_Y/2)/g_I = 2.55×10⁶ H_i
(3.4×10¹⁷ GeV) at λ = 1. Any inflaton quartic λ_I ≥ g_I⁴/(λ_Σλ_Y²) = 4.8×10⁻²⁵ (λ = 1) or 3.8×10⁻²⁴ (λ = 0.5) makes
the joint quartic form positive. Any real inflaton potential has that.

**Decays [D, E].** In the Higgs phase the trilinear is a φ_I Y Y coupling with A_I = g_I h/√2 = 2.35×10⁻⁵ H_i.
- φ_I → YY needs m_I > 2m_YL = 44.5 H_i, and φ_I → ΣYY needs 101 H_i. Both are closed for m_I = 3 and 5 H_i.
- The open channel φ_I → 6X through two virtual Y (c₆ = A_Iκ²/m_YL⁴ = 1.1×10⁻¹⁵ H_i⁻³) has Γ ≈ 2×10⁻³⁸ H_i at b = 3
  and 8×10⁻³⁷ H_i at b = 5 [E].
- Inflaton–σ mixing through the induced A_IΣ is of order 10⁻⁹ and negligible.

The inflaton's visible decay Γ_I is untouched.

**CP phase [D, C].**
- With κ and μ made real by rephasing Y and Σ, **α_I = arg(G/μ) is the rephasing-invariant CP phase**. The numerical
  check over random U(1)³ rephasings gives residuals ≤ 4×10⁻¹⁶. φ_I is real and cannot be rephased.
- Unlike R = 12H² > 0, φ_I has a sign. If the inflaton sector is reflection-even (½m²φ², Γ_I), then α_I ≃ α_I + π. The
  physical CP-odd input is then **α_I together with the sign of φ_I during inflation**: a Lagrangian phase plus a
  horizon-coherent initial condition.
- If φ_I is CP-odd, the CP violation is α_I − π/2 = 0.621 rad.
- Computed: the conjugate, u(0) = −10i, gives exactly −N_X (ratio −1.000000000).
- The opposite branch, φ_I(0) < 0 with u(0) = 2μ − 10i and |u| = 17.47, prepares a different, Y-condensed minimum:
  (125.94, 16.62, 29.69) H_i, δ = 2.532 rad, lightest mass 22.5 H_i², the only minimum. It yields
  **N_X = +26118.6 = −3.93×** the reference, with N_YL = 41.8 (R_Y = 9.7).
- So the asymmetry's sign is not fixed by the Lagrangian alone unless the inflaton potential is asymmetric.

## 2. The preparation point [C]

At μ_eff = 10i (λ = 1) the conditional minimum at fixed φ_I(0) is **(87.19556213287578, 5.754247963285591,
28.38718194436162) H_i**:
- phase −π/2 exactly; V = −32480.5853114062;
- equal to Option 1A's point to 2.4×10⁻¹⁵, V to 1×10⁻¹⁵;
- Hessian (13.09141, 14.58925, 538.147, 1546.66, 3239.89) H_i² plus a gauge zero, equal to 1A to 4×10⁻¹⁴;
- **every physical mass² is ≥ 13 H_i²**;
- a 300-start global search finds only this minimum.

**With the inflaton as a seventh direction:** the mixing vector has norm 3.9×10⁻⁴ H_i². It shifts the inflaton mass² by
−8×10⁻¹⁰ and the hidden eigenvalues by ≤ 6×10⁻¹⁰. The hidden force on the inflaton is 9.5×10⁻¹² of its tree force.

**At fixed H, the prepared family is smooth in s = φ_I/φ_I(0) over [0, 2].**
- It is always the global minimum, with δ = arg u + π.
- The lightest mass² is ≥ 4.47 H_i² > 9/4 (minimum at s = 0.35) and ≥ 13 H_i² for s ≥ 1.
- There is no Y instability up to |u| = 21 H_i, so the holder tracks the inflaton field adiabatically.
- dδ_eq/d ln φ_I = μ/10 = 0.716, the same as 1A's dδ_eq/d ln R.

## 3. Evolution through reheating [C]

The runs are the unheld original, Option 1A and the trilinear, each run to H_i t = 2000 ((3, 0.5)) or 1000 ((5, 1)) with
rtol 10⁻¹⁰. "Late N_X" is the 400–1000 mean, with the spread of running 80-wide windows over 400–1000.

| history, λ | run | N_X (400–1000) | ratio to original | 80-window ratio range | 1000–2000 ratio | N_YL | R_Y | R_central (e = 10⁻⁷) | Y/Y_Q at M1 H_i | H_i (10¹¹ GeV) |
|---|---|---|---|---|---|---|---|---|---|---|
| (3, 0.5), 1 | original | −8119.22 | 1 | — | 1 | 0.0068 | 0.008 | 0.008 | 0.99981 (640–800: 1.000000) | 1.3430 |
| | Option 1A | −8274.05 | **1.01907** (640–800: 1.018065) | 1.0175–1.0208 | 1.01896 | 3.08 | 2.54 | 3.23 | 1.01888 | 1.3262 |
| | **trilinear** | −6644.80 | **0.81840** (640–800: 0.81711) | 0.8161–0.8201 | 0.81791 | 0.256 | 0.199 | 0.245 | 0.81825 | **1.5350** |
| (5, 1), 1 | original | −5208.34 | 1 | — | — | 0.018 | 0.029 | 0.029 | 1.00101 | 1.3419 |
| | Option 1A | −7053.57 | **1.35428** | 1.3527–1.3557 | — | 39.1 | 42.5 | 69.4 | 1.35565 | 1.0963 |
| | trilinear | −5563.39 | **1.06817** | 1.0654–1.0707 | — | 1.65 | 1.76 | 2.98 | 1.06925 | 1.2842 |

**Checks.**

| Check | trilinear runs | Option 1A runs |
|---|---|---|
| Gauss | ≤ 1.9×10⁻¹¹ | ≤ 2.3×10⁻¹⁰ |
| constraint vs Raychaudhuri H | ≤ 6×10⁻¹² | ≤ 1.8×10⁻⁸ (its holder F-term is dropped, as in m1_holder) |
| hidden work ledger | ≤ 1.2×10⁻¹² | ≤ 8.6×10⁻¹¹ |

Other controls:
- rtol 10⁻¹² changes N_X by ≤ 1.2×10⁻¹² in all five rtol controls ((3, 0.5) trilinear; (5, 1) 1A and trilinear at λ = 1 and 0.5).
- 1A in `rptri.c` reproduces the unmodified `rpholder.c` binary bit for bit at (3, 0.5) and (5, 1), λ = 1 and 0.5.
- The original and 1A reproduce m1_holder's −8120.740082 and −8267.437810 (W80) to 5×10⁻¹³.

**Exchange [C].** In every run corr(N_X, 6N_Σ) = −1.000000 (≥ 0.9986 for the strongly excited 1A (5, 1)). In the
trilinear run std(N_X + 6N_Σ) = 7×10⁻⁴ against std(N_X) = 10.1: the partner is Σ, as the M1 skeptic found. The
±0.2–0.25% windowing spread averages down: the 400–1000 and 1000–2000 ratios agree to 6×10⁻⁴ (trilinear) and 1×10⁻⁴
(1A). **Option 1A's long-window value is +1.90% (1.01907 / 1.01896)**, not the 640–800 window's +1.806%. That gives
H_i = 1.3262×10¹¹ GeV, inside the review's 1.3269 ± 0.0015.

**Mechanism [C].** The yield is set by how non-adiabatically the CP-odd phase is switched off while X makes its first
roll. The charge fraction reached at H_i t = 0.8 is 0.87 (original), 0.71 (1A) and 0.49 (trilinear).
- Near the slice φ_I/φ_I(0) ≈ 1 − (bt)²/2 while R/12 ≈ 1 − 3(bt)²/2, so the trilinear holds u near 10i about √3 times
  longer during the first X half-orbit. The condensate partly follows adiabatically, and less charge is kicked.
- φ_I then changes sign. Im u reaches −2.62 (−3.15 at b = 5), against −0.71 for 1A, which briefly reverses the
  torque.
- A faster inflaton (b = 5) gives a less adiabatic release and more charge for both holders.
- Switch-off envelope |u − μ| at H_i t = 2, 5, 10, 20, 50: 3.2, 0.57, 0.083, 4.8×10⁻³, 1.1×10⁻⁶ (trilinear, ∝ inflaton
  amplitude) against 0.84, 0.026, 5.5×10⁻⁴, 1.8×10⁻⁶, 2×10⁻¹³ (1A, ∝ amplitude²). **The switch-off is built in:**
  u → μ (real) as the inflaton decays.

**Debris [C, E with our ledger].** The trilinear excites the light Y 12× less than 1A at (3, 0.5) (0.256 against 3.08)
and 24× less at (5, 1) (1.65 against 39.1). The ledger gives R_central = 3.2 for 1A at (3, 0.5), not the ≈ 10 that
m1_holder estimated by scaling [E]. 1A at (5, 1) reaches 69 (e = 10⁻⁷) and 51 (e = 10⁻⁶), which fails the
energy-conserving proxy of 7–8.

**History dependence (NEW) [C].** With no holder, (3, 0.5) and (5, 1) agree to 0.1% (1.0000 and 1.0010). With a holder
that switches off with the inflaton, the dimensionless yield depends on the inflaton history:
- Option 1A: +1.9% → +35%;
- trilinear: −18% → +7%.

So the M1 calibration's H_i is tied to the supplied (3, 0.5) history at the 15–20% level in H_i. This holds for both
holders.

## 4. Perturbative couplings λ_Y = λ_Σ = 0.5 [D, C]

**Bounds.**
- κ/κ_max = κ/(16λ_X³λ_Y/27)^{1/4} = **0.834** at tree level. cf_scale gives 0.83 including the running [Q].
- λ_Y,eff = 0.2585 > 0, so the hidden quartic form is positive.
- Late-vacuum large-Σ margin 1058 > 0 (2445 at λ = 1). A 200-start search confirms the vacuum (0, 0, h/√2) is global.
- At the slice (|u| = 10) the sufficient large-Σ test fails (−16.3). Codex's chain still gives V ≥ −9.0×10⁵ H_i⁴
  (reserve 0.1). B > 0 needs reserve < 0.215, so the λ_X margin is 1.27 (1.6 at λ = 1). The 300-start search finds only
  the prepared point.
- No 4π below M_P for λ ≤ 0.566 [Q cf_scale]. g_I enters the hidden β-functions only at O(g_I²) ~ 10⁻¹².

**Re-solved preparation.**

| Preparation | Point (H_i) | Phase | Mass² (H_i²) |
|---|---|---|---|
| Holder | (90.79987, 6.88773, 28.57624) | −π/2 | 10.191, 18.757, 471.16, 1536.28, 1665.74 |
| Held-phase original | (76.4675, 2.3399, 28.3025) | — | — |

For the holder preparation, the lightest mass² is **10.19 H_i² < 13** but > 9/4, so it is still heavy and
de Sitter-suppressed. g_I is unchanged (8.3067×10⁻⁷).

**Yield at λ = 0.5** (same windows and checks as above; Gauss ≤ 1.9×10⁻¹¹, H ≤ 6×10⁻¹², work ≤ 1.2×10⁻¹² for the
trilinear):

| history | run | ratio to original | Y/Y_Q at M1 H_i (80-window range) | **H_i (10¹¹ GeV)** | N_YL | R_Y | R_central e = 10⁻⁷ / 10⁻⁶ |
|---|---|---|---|---|---|---|---|
| (3, 0.5) | original | 1 | 1.01363 | 1.3308 | 0.0067 | 0.008 | 0.008 |
| | Option 1A | 1.17906 | 1.19513 | 1.1924 | 6.84 | 4.99 | 6.5 / 5.5 |
| | **trilinear (repaired point)** | **0.97741** | **0.99073 (0.98890–0.99258)** | **1.3512 (1.3495–1.3529)** | **0.567** | **0.379** | **0.44 / 0.40** |
| (5, 1) | original | 1 | 1.01492 | 1.3296 | 0.031 | 0.051 | 0.05 |
| | Option 1A | 1.60639 | 1.63037 | 0.9694 | 69.5 | 64.9 | 111 / 80 |
| | trilinear | 1.27145 | 1.29042 | 1.1329 | 4.30 | 3.96 | 7.1 / 4.9 |

The recalibration uses Y_Q ∝ H_i^{3/2} at fixed dimensionless history: H_new = H_M1 (Y/Y_Q)^{−2/3}. The target is
Y_Q = 5.07941423×10⁻⁹ at E/Q = 0.9374 m_X. Our convention reproduces the M1 point exactly: the λ = 1 original gives
Y/Y_Q = 1.000000 in the 640–800 window. The 640–800 window gives 1.3498 and the 1000–2000 window 1.3508 for the repaired
point. **Repaired point: H_i = 1.351×10¹¹ GeV ± 0.13% (windowing)**, for the supplied (3, 0.5) history only. The
(5, 1) history gives 1.133×10¹¹.

**Across the band (3, 0.5)** (`lam_scan.json`):

| λ | Y/Y_Q | H_i (10¹¹ GeV) | lightest mass² (H_i²) | N_YL | κ/κ_max |
|---|---|---|---|---|---|
| 0.25 | 1.556 | 1.000 | 3.5 | 1.85 | 0.99 |
| 0.3 | 1.261 | 1.150 | 6.5 | 1.24 | 0.95 |
| 0.4 | 1.077 | 1.278 | 8.9 | 0.77 | 0.88 |
| 0.5 | 0.991 | 1.351 | 10.2 | 0.57 | 0.83 |
| 0.566 | 0.952 | 1.387 | 10.8 | 0.49 | 0.81 |
| 0.75 | 0.880 | 1.463 | 12.0 | 0.35 | 0.75 |
| 1 | 0.818 | 1.535 | 13.1 | 0.26 | 0.70 |

The unheld original moves only from 1.022 to 1.000 over the same range. Near λ = 0.5, d ln Y/d ln λ ≈ −0.35, so
d ln H_i/d ln λ ≈ +0.24. **The holder makes the yield about 20× more coupling-sensitive.**

**Clock sensitivity [C].** The re-prepared response to a 1% shift of the holder coupling is
d ln|N_X|/d ln g_I = 1.71 (trilinear) and d ln|N_X|/d ln ξ = 1.34 (1A). m1_holder's "0.4% per 1%" [E] counted only the
phase. The amplitude change (|u| moves 1% per 1%) makes the true response about 3× larger. For the trilinear, the
inflaton field value at the end of inflation is now a calibration input, in place of H.

## 5. Verdict

**The repaired holder is a genuine improvement in control, and it costs a recalibration and a relocated CP sign. It does
not fail, so no fallback is forced.**

Gains:
- It is renormalizable and perturbative (g_I = 8.3×10⁻⁷), with no Hubble-scale cutoff [D].
- It is technically natural. Inflaton mass corrections are ≤ 2×10⁻⁹ H_i² against m_I² = 9 [C].
- It reproduces 1A's prepared point and its ≥ 13 H_i² spectrum exactly [C].
- Its switch-off is built in: u → μ as the inflaton decays [C].
- It produces 12–24× less light-Y debris than 1A, so the release stays quiet [C, E].

Costs:
1. **Yield −18.2%** at λ = 1 (H_i → 1.535×10¹¹), −2.3% at λ = 0.5 (H_i = 1.351×10¹¹) [C].
2. **The inflaton must be sequestered and shift-symmetric.** All other portals ≲ 5×10⁻¹⁵, Planck-suppressed couplings
   ≲ 10⁻⁸, and an inflaton quartic ≥ 5×10⁻²⁵ is needed for boundedness. These are supplied assumptions, natural at one
   loop [D, A].
3. **The CP input is now α_I plus the sign of φ_I** if the inflaton sector is reflection-even. The other branch gives the
   opposite sign and 3.9× the charge [D, C].
4. **The clock is φ_I, not H**: d ln N_X/d ln φ_I(slice) = 1.71 [C].
5. **History and coupling sensitivity**: +29% between (3, 0.5) and (5, 1) at λ = 0.5, and −0.35 per ln λ [C].

Option 1A is worse on 5 and on debris. On its own terms 1A remains usable at (3, 0.5) with λ = 0.5: H_i = 1.192×10¹¹,
R_central 6.5, i.e. marginal under the energy-conserving proxy. Its caveat stands: an uncontrolled operator tower at
Λ_eff ≈ H_i.

Proposed working holder: the inflaton trilinear with λ_Y = λ_Σ = 0.5 and H_i = 1.351×10¹¹ GeV, on the supplied
(3, 0.5) history. Both the holder and the couplings are model choices.

## 6. Limits

- The calculation is homogeneous and classical, with the quadratic-inflaton preparation slice of our engine. The
  preceding inflation is represented by the static conditional minimum at φ_I(0). The adiabatic-tracking statement
  rests on the s-scan (masses > 9/4 H²) [A].
- There are no finite-k Y modes. The trilinear modulates the light-Y mass² by up to 160 H_i² at the slice (495 → 334),
  oscillating at m_I with the decaying inflaton amplitude. That would also drive finite-k Y production. It is adiabatic
  (ω_Y ≈ 22 ≫ m_I) and not computed.
- The X ↔ Σ exchange was not followed to freeze-out beyond H_i t = 2000.
- Debris uses our decay ledger, which is an estimate.
- CW corrections are hidden-sector loops only, at Q = M_Y. Portal and Planck-suppressed bounds are order-of-magnitude
  estimates.
- The opposite-branch run starts from a Y-condensed state where the forced-Y diagnostic is undefined at t = 0
  (min_mL2 < 0). Its N_YL is indicative only.

## Reproduction (one core, about 2.5 min)

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_holder/code
python3 -B tri_point.py      # items 1, 2, bounds, CP, s-scan        -> ../results/tri_point.json  (25 s)
python3 -B tri_runs.py       # items 3, 4: 19 engine runs              -> ../results/tri_runs.json   (50 s)
python3 -B lam_scan.py       # band scan, 14 runs                      -> ../results/lam_scan.json   (22 s)
gcc -O3 -std=c11 -o $FX_SCRATCH/build/rpholder_ref ../../../../grut_formation_review_20261025/tracks/m1_holder/code/rpholder.c -lm
python3 -B controls.py       # rtol 1e-12, rpholder.c cross-check      -> ../results/controls.json   (40 s)
python3 -B summary.py        #                                          -> ../results/summary.json
```
(FX_SCRATCH defaults to the session scratchpad `wt2/fx_holder`; the engine binary is built there.) Credit: the
Coleman–Weinberg potential and 't Hooft naturalness are standard [R]. Our earlier tracks are reused as stated above.
