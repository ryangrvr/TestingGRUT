# Track hw_gates: the open gates of the heavy-mediator window, tested at the ls_scaling spot point

26 October 2026 · Claude, issue #13 · same-session AI work, not outside replication.

Packet read: Codex 2908e91 (`origin/codex/grut-concept-checkpoint-20261025`), through `git archive` into the scratchpad
only. No Codex code was run. Our own earlier tracks are reused where the brief allows (sp_vector, sp_phaseprep,
gp_stability, gp_matching methods; ls_scaling numbers as [Q, own track]). Nothing was committed.

**Scope.** One supplied tree-level working model. Homogeneous classical release, instant at a = 1 into radiation; a
supplied CP phase; the late-particle decay ledger and the case-L threshold of our sp_critic tables. Nothing here is
evidence that X, Y, Σ or the vector exist. It is not a validation, an exclusion or a viable point. Comparisons with
PTA/LVK/Planck data are retrospective.

Labels: [D] derived, [C] computed, [E] estimate, [Q] quoted, [A] argued, [R] recalled.

## The spot point (ls_scaling B_base_Hscaled) [Q, own track]

H_i = 1.3335e11 GeV, M_Y = 30 H = 4.00e12, h = 40 H = 5.33e12, μ = 7.159 H = 9.55e11 GeV (r = √2μh/M_Y² = 0.45),
κ = 3.459e-3, λ_X = 1e-3, λ_Y = λ_Σ = 1, e = 1e-7, c_H = 10; Y-relaxed preparation at δ = π/2; yield 1.00 Y_Q.

**Own reproduction of the point** (`code/hwgrel.c`, DP5(4), cosmic time; window a = 30–60): ε = −0.14568,
N_X = −2788.2, E_X = 38398, N_YH = 0.01725, N_YL = 0.33446, N_σ = 3.2e-9, R_Y = 0.805, R_casc = 0.339, R_central = 1.144,
R_L = 87.85, yield 1.0000 Y_Q [C]. ls_scaling: ε = −0.14572, R = 1.144, R_L = 87.9. Codex's preparation at the same point:
R = 1415 (ls_scaling 1414.9). Gauss residual 6e-12, energy ledger 3e-10.

**Inflationary state** (own 6×6 canonical Hessian, gauge mode projected): aligned minimum |X|, |Y|, |Σ| = 80.20, 3.442,
28.311 H; physical masses² 6.41, 14.77, 573, 1334, 3210 H². Higgs norm F_inf = 266.1 H against F_vac = 6h = 240 H.

## Gate 1: release-induced longitudinal vector / Goldstone production

[PENDING: filled from results/gate1_vector.json and results/gate1b_pieceB.json]

## Gate 2: phase holding, and what sets the phase value (`code/gate2_phase.py`, `code/spotmodel.py`)

**The relative phase and its mass [C].** The charge-sourcing phase is δ = θ_Σ − 6θ_X (= ψ₂ + 2ψ₁). At the aligned
inflationary minimum the gauge-projected 6×6 Hessian gives the physical modes m² = 6.41 (light phase, pure angular:
dθ_X, dθ_Y, dθ_Σ = 0.0080, 0.0061, −0.0107 per canonical H), 14.8 (X radial), 573 (ψ₁-type phase), 1334 (Y radial),
3210 (Σ radial) H². So the held relative phase has **m = 2.53 H, f_eff = 1/|dδ/dφ_c| = 17.1 H**.
ls_scaling's 1.84 H [D] used the leading |X| = 70.7 (actual 80.2) and f² = |X|²/18 without the Σ kinetic weight:
**CORRECTED to 2.53 H** (both are above 3H/2, so the qualitative statement stands). Along the Y-relaxed path the
conditional potential V_rel(δ) has its minimum at δ = 0, half-barrier B = 1473 H⁴, d²V/dδ² = 1796 ≈ m²f² = 1873.

**Nothing in the action chooses δ = π/2 [C].** The de Sitter equilibrium P ∝ exp(−8π²V/3H⁴) gives
P(π/2)/P(0) = exp(−43,277): the phase sits at the torque-free minimum δ = 0, where ε = 0 exactly. At the three
supplied states the Hessian has a negative or small curvature along δ (Codex state −0.17, Y-relaxed −0.83, κ-rotated
+0.37 H², all with nonzero phase gradient): none is held. A term that moves the minimum, with a CP phase relative to
the real couplings, is required (as for 373d37e [Q sp_phaseprep]); the stated action has no physical CP phase at
dimension ≤ 4 with Z2_φ [Q gp_phase].

**A held phase is quiet [C].** If a holder fixes δ = π/2 with the same mass, the exact Bunch–Davies coarse-grained
spread is σ_δ = 3.7e-3 rad (x_c = 1), 1.2e-4 (x_c = 1/10, modes still outside the horizon when the charge forms at
a ≈ 2–10), 1.5e-5 (1/40). The band with |ε| ≥ 0.10 is δ ∈ [0.44, 2.15] rad (|ε| ≥ 0.13: [0.69, 1.83]) from our own
phase scan, so P(band) = 1 for any holder mass ≳ 1 H (σ_δ ≤ 3.6e-3 at m = 1 H). **Isocurvature:** with
T = d ln N_X/dδ = 0.31 per rad at π/2 (smooth map) and the e^{−3N/2} decay of a mode with m > 3H/2, S/S_max ≈ 1e-30
(m = 2.53 H, N_* = 50); even m = 1 H gives 1.3e-6. The isocurvature obstruction of the 373d37e point (×16–95) is
absent here.

**Classical drift if the holder stops before inflation ends [C]** (1D, f²δ̈ + 3f²δ̇ + V_rel′ = 0): δ = 1.55, 1.49,
1.40, 1.17, 0.48, −0.10 rad after 0.1, 0.2, 0.3, 0.5, 1, 2 e-folds, i.e. ε from the scan −0.146, −0.148, −0.151,
−0.151, −0.104, −0.03. The holder must stay on until ≲ 0.7 e-folds before the end of inflation (373d37e: 0.02 e-folds
[Q sp_phaseprep]).

**Three holders and what they leave behind [D, C].** Each makes a supplied state an exact (D1, D2) or
near-exact (D3, gradient 9e-4) inflationary minimum; with the holder on, the held phase mass is 2.53–2.67 H.

| Design | Term during inflation | Size needed for δ = π/2 | State after instant removal |
|---|---|---|---|
| D1 μ-rotation | μ_P Σ*Y² + h.c. (Hubble-induced A-term, DRT type) | \|μ_P\| = √2 μ = 10.1 H = 1.35e12 GeV (a_H ≈ 10; δ = arctan(a_H H/μ) for an orthogonal A-term, so \|ε\| ≥ 0.10 needs a_H ≳ 3) | Codex's preparation (Y real) |
| D2 κ-rotation | κ_P Y*X³ + h.c. (dimensionless, with a phase: an inflaton–dark coupling of size κ) | \|κ_P\| = 0.77 κ = 2.6e-3 | aligned moduli, Y and Σ rotated by δ/2, δ |
| D3 X⁶Σ*-rotation | c7 X⁶Σ* + h.c. (dimension 7; does not enter the Y equation) | \|c7\| = √2 μκ²/D = 1.9e-10 H⁻³, i.e. Λ₇ = 2.3e14 GeV = 58 M_Y | ls_scaling's Y-relaxed preparation |

**Releases from the three designs** (own code; window a = 30–60; R_L from our case-L table; "×" = R_L/R_central):

| Switch-off of the holder | D1 (μ_P) | D2 (κ_P) | D3 (c7) |
|---|---|---|---|
| instant (Codex's assumption) | ε −0.159, N_Y 77.9, **R 1415 vs 101: FAIL ×14** | ε −0.159, N_Y 50.6, **R 657 vs 101: FAIL ×6.5** | ε −0.159, R 1.70 vs 101: pass ×60 |
| exp, Δ = 0.1 (in a; ≈ 0.1 H⁻¹) | ε −0.117, R 161 vs 63: FAIL ×2.5 | ε −0.127, R 53 vs 72: pass ×1.4 | ε −0.125, R 1.33 vs 70: pass ×52 |
| exp, Δ = 0.3 | ε −0.092, R 11 vs 45: pass ×4, yield 0.73 | ε −0.102, R 7.8 vs 51: pass ×6.6, yield 0.80 | ε −0.099, R 1.33 vs 49: pass ×37 |
| exp, Δ = 1 | ε −0.065, R 2.4 (ε below table), yield 0.51 | ε −0.076, R 2.1 vs 35: pass, yield 0.59 | ε −0.074, R 1.34 vs 33, yield 0.58 |
| ∝ H (DRT-like persistence) | ε −0.075, R 4.6 vs 34: pass ×7, yield 0.59 | ε −0.088, R 3.7 vs 42: pass ×11, yield 0.69 | ε −0.084, R 1.33 vs 40: pass ×30, yield 0.66 |
| permanent | integration failed at a ≈ 20 (unphysical case, not pursued) | same | ε −0.028: too little charge |

- **The debris of the heavy window is a sudden-quench effect.** A holder that enters the Y equation (D1, D2) and is
  removed in much less than Y's response time leaves Y displaced from its forced position; the displacement is
  radiated as N_Y ≈ 50–78 quanta. Removed over ≳ 0.3 H⁻¹ (≫ 1/M_Y = 0.033 H⁻¹) the same holders leave N_Y ≲ 0.5–3 and
  pass by 4–11; the price is ε ≈ 0.08–0.10 instead of 0.15 and a yield of 0.6–0.8 Y_Q at this H (the exact scaling
  symmetry restores Y_Q with H × 1.2–1.4 and c6 × 0.5–0.7). This is Codex's point (i) realised dynamically: the forced
  Y carries the torque without excitation when the preparation changes adiabatically on the Y time scale.
- ls_scaling's "the preparation is decisive" is CONFIRMED, but its reading that the Y-relaxed preparation is "what a
  heavy Y does during inflation" holds only for a holder that does not enter the Y equation (D3). For the natural
  DRT-type A-term (D1) the Y-relaxed state is not the inflationary minimum; what rescues D1 and D2 is the switch-off
  rate, not the preparation. **CORRECTED (scope).**
- Removing the X Hubble mass slowly (Δ = 0.1, 0.3, ∝ H) from the D3 state gives ε = −0.20, −0.33, −0.31 with R ≈ 0.4:
  the window does not need the sudden end of inflation that the 373d37e point needed [Q sp_phaseprep].

**What still sets the phase value.** The CP phase of the holder (arg μ_P, arg κ_P or arg c7), supplied by hand; its
size (a_H ≈ 3–10, or κ_P ≈ 0.2–0.8 κ, or Λ₇ ≈ 2e14 GeV); and its switch-off on a time scale between ~0.3 H⁻¹ and the
Hubble-tracking limit. Once held, the phase is quiet (spread ≤ 1e-4 rad, no isocurvature, P(band) = 1).

**Verdict:** a DRT-type term is still required (CONFIRMED, now quantified: e^{−43,000} without it); mass 1.84 H
CORRECTED to 2.53 H; the Y-relaxed preparation's physical status CORRECTED (holder-dependent); the cost is a
supplied CP phase plus a_H ≈ 3–10 and a non-instant switch-off; no isocurvature or band-probability obstruction.

## Gate 3: Z6 strings at h = 5.33e12 GeV (`code/gate3_strings.py`)

**Tension [C].** β = m_s²/m_v² = 2λ_Σ/(36e²) = 5.56e12 does not depend on h. The large-β law of the winding-1
Nielsen–Olesen string, ε(β) = ½ ln β + c, was re-derived with our gp_stability solver: ε − ½ ln β = 0.2378, 0.15318,
0.150301, 0.150241 at β = 1e2, 1e4, 1e6, 1e8, so c = 0.15024(6). Then ε = 14.823, μ = π h² ε = 1.325e27 GeV² and

  **Gμ = 8.89e-12** (ls_scaling's scaled estimate 9e-12 [E]: CONFIRMED). Control at the stated point: 3.124e-15,
  the gp_stability value.

If e were raised (gate 1), β falls and Gμ with it: 3.4e-12 (e = 1e-3), 2.4e-12 (5e-3), 2.0e-12 (1e-2).

**Formation conditions [C, A].** A network needs the hidden U(1) restored and re-broken.
- During inflation Σ sits at its vacuum modulus: |Σ|/(h/√2) = 1.0009 (the μY² forcing is 6e-4 at this point), with
  radial mass 56.6 H. X (charge 1) is also condensed. The U(1) is broken throughout inflation, and the Σ fluctuations
  are frozen out (m ≫ H), so there is no Kibble mechanism from inflation.
- Hidden temperature after release: the whole release energy is ρ = 3.51e4 H⁴. Even if all of it thermalised at once
  in the hidden sector (g_h = 10) the temperature would be 1.36e12 GeV = **0.25 h**. The symmetry is never restored
  thermally in the hidden sector.
- The visible reheating temperature is 3.08e14 GeV = 58 h, but Σ has no visible coupling (the Higgs portal acts on φ,
  κ_portal ≲ 1e-9 [Q 7737487]); it reaches the visible bath only through X ↔ φ and Y. No restoration.
- Along the release |Σ| never leaves its vacuum value (min |Σ|/vev = 1.0000 to 1e-4) and the Σ phase winds by only
  −0.238 rad in total (0.037 windings). No phase slips, so no strings from the release either (contrast the 373d37e
  point, where |Σ| fell to 1e-5 of its vacuum value in 7–40% of ensemble members [Q sp_vector]).

**So at the spot point no string network forms, and Gμ is not a constraint.** The tension is recorded in case a
different history restores the U(1).

**Retrospective limits, with read levels.**
- LVK O3 (arXiv:2101.12248): abstract re-read this session [Q]: "for one loop distribution model, we set the most
  competitive constraints to date, Gμ ≲ 4e-15". Per-model ranges [Q via our gp_stability track]: model A (BOS)
  excluded above 9.6e-9–1e-6; model B (LRS) above 4.0–6.3e-15; C-1 above 2.1–4.5e-15. The spot-point Gμ would be
  2.2e3 (model B) and 4.2e3 (C-1) times over those edges, and 1e-3 of model A's.
- NANOGrav 15-yr new physics (arXiv:2306.16219, §V.4 read this session [Q]): stable field-theory strings give Bayes
  factors 0.28–0.38 (or 0.76–0.89) against the SMBHB interpretation; when fitted to the signal the posteriors peak at
  log10 Gμ ≈ −10.5…−10.0; treating the signal as a bound gives log10 Gμ < −9.67, −9.87, −10.10, −9.71 (and
  −9.88…−10.26) across loop models. The spot point (log10 Gμ = −11.05) is 6–24 times below those limits and 10–30
  times below the strings-as-signal range. PPTA 2016 Gμ < 1.5e-11 for one model [R]; Planck 2013 XXV 1.5e-7 (NG) [Q via
  gp_stability].

**Verdict:** Gμ = 8.9e-12 CONFIRMED [C]; formation: not triggered at this point [C, A] (NOT a constraint); if a network
of model-B/C type formed, the point would be retrospectively excluded by LVK O3; PTA limits would not exclude it.

## Gate 4: one-loop running of κ against the stability bound (`code/gate4_running.py`)

Own tensor one-loop RG for the 8 real fields (X, Y, Σ, φ, D), 16-operator closed basis, pure scalar sector (the
U(1) terms are O(e²) and negligible for e ≤ 1e-2); control β(λ) = 20λ²/16π² reproduced to 3e-16; projection residual
3e-16. Bound [D, ls_scaling]: along the aligned large-field direction λ_eff = λ_X − (3/2)κ^{4/3}(2λ_Y)^{-1/3} > 0, i.e.
κ < κ_max = (λ_X/1.1906)^{3/4} λ_Y^{1/4}.

| μ (GeV) | λ_X | λ_Y | κ | κ/κ_max | λ_eff | V4 min on sphere |
|---|---|---|---|---|---|---|
| M_Y = 4.0e12 (matching) | 1.0000e-3 | 1.000 | 3.4587e-3 | **0.701** | 3.77e-4 | +6.7e-5 |
| \|Σ\|_inf = 3.8e12 | 0.99984e-3 | 0.993 | 3.4587e-3 | 0.702 | 3.76e-4 | +6.7e-5 |
| \|X\|_inf = 1.07e13 | 1.00271e-3 | 1.142 | 3.4589e-3 | **0.677** | 4.07e-4 | +7.1e-5 |
| 3\|X\|_inf = 3.2e13 | 1.00575e-3 | 1.358 | 3.4592e-3 | 0.647 | 4.43e-4 | +7.6e-5 |
| 1e14 | 1.0089e-3 | 1.688 | 3.4595e-3 | 0.611 | 4.86e-4 | +8.1e-5 |
| 1e15 | 1.0153e-3 | 3.325 | 3.4602e-3 | 0.514 | 5.98e-4 | +9.3e-5 |

- κ has no λ_Y term at one loop and runs by 7.6e-5 per e-fold; λ_Y grows (β = 20λ²/16π²) and λ_X grows (κ² and g_X⁴
  terms), so κ/κ_max **falls** from 0.70 at M_Y to 0.68 at the inflationary |X| and 0.61 at 1e14 GeV. λ_eff stays
  positive everywhere and the numerical minimum of V4 on the unit sphere stays positive.
- If the couplings are instead stated at |X|_inf and run down to M_Y, κ/κ_max = 0.79 there. Running formally below M_Y
  (where Y decouples and the bound no longer applies) gives 0.77 at H_i and 0.85 at 1e9 GeV.
- λ_Y reaches 4π at 5.7e15 GeV (analytic 5.73e15), far above every field value of the point (1.07e13). The generated
  portal λ_XY is 2.8e-6 at |X|_inf and enters λ_eff positively. Young's sufficient condition is not met (ratio 0.26–0.30)
  but it is only sufficient; the aligned-direction criterion is the binding one.

**Verdict: CONFIRMED.** κ stays below its boundedness bound between M_Y and the inflationary field values at one loop,
with a margin of 0.68–0.79 depending on the (unstated) matching scale. Not checked: two loops, the Hubble-induced
mass in the effective potential during inflation, and thresholds at M_Y.

## Gate 5: consequences for the light sector (`code/gate5_light.py`)

**Unchanged.** The Candidate 1 object depends only on g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV and
Q = 200000: own fixed-charge solver (gp_stability code) gives Ω = 0.861195, E = 17064.537 GeV, binding margin
1417.255 GeV, ∫n² = 47241.83 GeV³, Q = 2 mod 6 — the gp_stability values [C]. The Z6 (charges 1, 3, 6; invariant
(Q_X + 3Q_Y) mod 6 on a closed slice) and the string type β are unchanged.

**Changed with h = 5.33e12 GeV (and e = 1e-7)** [C, D]:

| Quantity | stated h = 1e11 | spot h = 5.33e12 |
|---|---|---|
| m_A = 6eh | 6.0e4 GeV | **3.20e6 GeV** |
| screening length 1/m_A; m_A r_1/2 | 1.7e-5 GeV⁻¹; 2.5e6 | 3.1e-7 GeV⁻¹; 1.3e8 |
| local screened gauge energy ∫n²/(72h²) | 6.56e-20 GeV (4.6e-23 of margin) | **2.31e-23 GeV (1.6e-26 of margin)** |
| unscreened Coulomb bound (e only) | 3.40e-7 GeV (2.4e-10 of margin) | same |
| m_A / E_object | 3.5 | **188** (the "four objects exceed m_A" caveat disappears) |
| Γ(A → XX̄) = e²m_A/48π; τ; T_vis at H = Γ | 4.0e-12 GeV; 1.7e-13 s; 1.7 TeV | 2.1e-10 GeV; 3.1e-15 s; **12 TeV** |
| Gμ (if strings formed) | 3.1e-15 | 8.9e-12 (gate 3) |

If e were raised to 5e-3–1e-2 (gate 1 remedy), the unscreened Coulomb bound would grow to 0.6–2.4 of the binding
margin, so the Higgs screening (which is e-independent, m_A r_1/2 ≫ 1) becomes load-bearing for the object; m_A becomes
1.6–3.2e11 GeV and the vector decays promptly.

**Late survival c6 [D, C].** α6 = μκ²⟨Σ⟩/(M_Y⁴ − 4μ²⟨Σ⟩²) = 2.109e-31 GeV⁻², c6 = α6 M_red² = **1.2504e6** (ls_scaling
1.2504e6). With the unit-coefficient X⁶ core-loss rates (coherent 3.92e-50 /s, pair 1.040e-46 /s [Q xcharge /
Codex 8986d34]) the rate is 1.63e-34 /s and Γt₀ = 7.1e-17; the 1% ceiling is c6 = 1.486e13, a margin of 1.2e7.
Planck-suppressed, gauge-invariant operators shift the late α6 by: Σ*X⁶/M_red³ (dimension 7) 1.2e-12; Σ*YX³/M_red
(dimension 5, mixes with κ after Y is integrated out) **1.6e-3**; the dimension-5/6 radial operators at the core
field 1.7 GeV, < 1e-30. The U(1) gauge symmetry forbids every ΔQ_X ≠ 0 mod 3 operator (X⁵-type operators do not exist),
and with Y, Σ heavy the lowest X-only operator is X⁶ (Z6): c6 = 1.25e6 is unchanged to 0.2% by Planck-suppressed
operators with unit coefficients. The radial companion in the core is 2.7e-27 of λ_X|X|⁴.

**Verdict: CONFIRMED** (object unchanged; vector mass, screening energy, decay epoch and Gμ change as tabulated; c6
unchanged by Planck operators).

## Verdicts

[PENDING: table]

## Limits

[PENDING]

## Reproduction

```
T=research/reality_framework/grut_formation_review_20261025/tracks/hw_gates
cd $T/code
python3 -B gate1_vector.py ../results/gate1_vector.json        # ~15 min (Proca modes, e scan, Floquet)
python3 -B gate1b_pieceB.py ../results/gate1_vector.json ../results/gate1b_pieceB.json   # ~3 min, 4 cores
python3 -B gate2_phase.py ../results/gate2_phase.json          # ~15 min, 4 cores (24 releases)
python3 -B gate3_strings.py ../results/gate3_strings.json      # ~3 min (uses ../../gp_stability/code)
python3 -B gate4_running.py ../results/gate4_running.json      # 20 s
python3 -B gate5_light.py ../results/gate5_light.json          # 70 s (uses ../../gp_stability/code)
```
`hwgrel.c` is compiled automatically (gcc -O3) into the scratchpad build directory.
