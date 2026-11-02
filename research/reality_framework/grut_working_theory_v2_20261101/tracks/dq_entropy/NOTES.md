# Track dq_entropy: end-to-end S/Q of the M1r release, and the X ↔ Σ exchange to freeze-out

2 November 2026 · Claude, decisive pass (BRIEF_DECISIVE.md of this directory; issue #13 context). This is same-session AI
work on one supplied model, not outside replication. Nothing here is observational evidence, a validation, an exclusion
of the real world or a new theory. **Every input is a model choice, not a consequence of closure.** That covers M1r, the
trilinear holder, λ_Y = λ_Σ = 0.5, the (3, 0.5) inflaton history, e, ε and the portal renormalisation condition. Git
was read-only and nothing was committed. Labels: [D] derived, [C] computed with our own code, [E] estimate, [Q] quoted
from a reviewed track, [A] argued, [R] recalled.

## Code and reuse

| File | What | Reuse |
|---|---|---|
| `code/rpdq.c` | engine | copy of `tracks/fx_holder/code/rptri.c` with three changes. **DQ1:** the forced-Y diagnostic uses c = 2 conj(u)Σ; its time derivative includes u̇ and the portal mass; DQ1′ makes the forced point self-consistent with the Y quartic, so the free-Y number is exactly 0 at the prepared minimum. This fixes the known conj bug (fx_holder SKEPTIC N6). **DQ2:** optional induced-portal counterterm c(M_Y)(\|u\|² − μ²)(\|Σ\|² + 2\|Y\|²), with exact inflaton back-reaction (SKEPTIC N1). **DQ3:** conformal-time bins |
| `code/dq_common.py` | driver | imports fx_holder's `fx_common.py` read-only (preparation, G, yield convention, ledger) |
| `code/long_run.py`, `code/freeze_analysis.py` | item 1 | own |
| `code/hartree_dq.py` | Hartree solver | copy of `tracks/fx_chain/skeptic/code/hartree_indep.py`, identical to wt_critic's `hartree_copy.py`. Only change: `set_release()` points it at another (N_X, E_X, H_i) |
| `code/sq_scan.py` | item 2 | own; thresholds read from `fx_chain/skeptic/results/survival_threshold.json` |
| `code/entropy_increments.py` | item 3 | own; imports `fx_chain/code/heat_onset.py` read-only (m1_vacuum homogeneous functional) |
| `code/summary.py` | item 4 + ledger | own; p(S/Q) mapping read from `wt_critic/results/chain_consistency.json` [E] |

Results: `results/freeze.json`, `sq_scan.json`, `entropy_increments.json`, `summary.json` (the `for_critic` block holds
S/Q(R), R_crit, W_crit and the closed form). Long-run series are in the scratchpad
(`wt2/dq_entropy/runs/*.npz`).

**Reproduction of fx_holder [C].** With cport = 0, rpdq.c reproduces the fx_holder λ = 0.5 trilinear run bit for bit:
- N_X(400–1000) = −8045.485140302082;
- N_YL = 0.566640911258934, N_YH = 0.010377131381878653;
- E_X = 188706.29275778163.

The conj fix changes only the early diagnostic.
- Fixed free-Y numbers: N_YL at H_i t = 0 / 0.2 / 0.4 / 1 / 2 / 5 / ≥ 10 is 0 / 2.05 / 0.84 / 0.47 / 0.64 / 0.556 / 0.567.
- rptri gave 0 / 1.61 / – / 0.48 / 1.32 / 0.560, and N_YH = 1009 at t = 0, with the directions swapped.
- Late values agree to 10⁻¹⁵.

## Result in one paragraph

**The X ↔ Σ exchange freezes as a⁻²:** d ln std(N_X)/d ln a = −2.03 (benchmark) and −2.01 (Planck). Extrapolated
over H_i t ≤ 10⁵ (a = 630):
- N_X = −8047.84 ± 0.22, which is +0.029% on the 400–1000 window;
- Y/Y_Q = **0.991023 ± 0.000027** within the homogeneous model;
- E_X = 188737.7 ± 1.1 (conformal), ε = 0.127374, E/E_c = 15.6036;
- R_med(10⁻⁶) = 0.3970 [C].

**The relaxed Hartree plasma depends only on W = (E/E_c)(1 + R).** A direct solve on the M1r release and the
spot-normalised solver agree to 10⁻¹⁵ [C]. It obeys S/Q = 0.190255 (W − 0.97494)^0.75013 − 0.00428 (residual
≤ 7×10⁻⁶) [C].

**Critical R for this release** (R_crit = W_crit/(E/E_c) − 1 for any release):

| threshold | R_crit, Hartree | R_crit, with entropy production | W_crit (Hartree) |
|---|---:|---:|---:|
| bulk, S/Q = 3.4716 | 2.1442 | 2.136–2.139 | 49.06 |
| Q = 2×10⁵, S/Q = 2.9284 | 1.5195 | 1.511–1.514 | 39.31 |

**Post-relaxation entropy production is small [C/E]:**
- spinodal latent heat +0.004–0.0055;
- surface and ripening heat ≤ 0.0013;
- crossover and evaporation ≤ 10⁻⁴ (with vev-opened chemistry);
- X̄ burn-out: zero at equilibrium.

**Branch** with the consistent transfer of the provisional heat (M1r's own debris plus spot vectors in E_c units,
ν = 0.960):

| e | S/Q (band) | branch |
|---|---:|---|
| 10⁻⁶ | 4.685 | free X |
| 2×10⁻⁶ | 3.370 | **Q-matter** |
| 3×10⁻⁶ | 2.894 | Q-matter |
| 5×10⁻⁶ | 2.493 | Q-matter |

- e_crit = 1.86×10⁻⁶ (band) / 2.20×10⁻⁶ (broad).
- The literal provisional R(e) gives e_crit = 3.19×10⁻⁶ instead. It double-counts the spot's Y_L debris.
- With portals vanishing at M_P, the recomputed release has S/Q(10⁻⁶) = 4.01 and e_crit = 1.41×10⁻⁶.

## 1. X ↔ Σ charge exchange to freeze-out [C]

**Runs.** M1r prepared state (90.79987, 6.88773, 28.57624) H_i at u = 10i, (3, 0.5), rtol 10⁻¹⁰, to H_i t = 10⁵
(η = 317, a = 630, 3.9×10⁷ steps, 106 s).
- Checks: Gauss ≤ 5.6×10⁻¹⁰, constraint ≤ 3.3×10⁻¹⁰, hidden work ≤ 2.7×10⁻¹⁰.
- Control: rtol 10⁻¹² to 2×10⁴ changes every common conformal window mean of N_X by ≤ 2.0×10⁻⁵.

**Conformal-window means** (Δη = 10; benchmark, portals zero at M_Y). The X orbit has a conformal period of about 1.7.

| η | a | ⟨N_X⟩ | bin std of N_X | ⟨E_X⟩ |
|---|---:|---:|---:|---:|
| 10–20 | 30 | −8058.55 | 37.0 | 188761.6 |
| 40–50 | 90 | −8048.44 | 3.31 | 188734.4 |
| 70–80 | 149 | −8047.70 | 1.20 | 188741.0 |
| 100–110 | 209 | −8047.66 | 0.655 | 188735.2 |
| 160–170 | 328 | −8047.87 | 0.252 | 188736.0 |
| 220–230 | 447 | −8047.83 | 0.129 | 188736.5 |
| 300–310 | 606 | −8047.84 | 0.077 | 188736.3 |

The fluctuation obeys the a⁻² law: d ln std/d ln a = −2.03 for a > 60. The exchange partner is Σ: N_X + 6N_Σ + 3N_Y ≤ 10⁻⁵.

**Frozen values.** The adopted fit is N_inf + c₂a⁻² over a ≥ 210 (Δη = 20 windows). The uncertainty is the largest
deviation among the a⁻² fits over a ≥ 105 and ≥ 315, the a⁻¹ + a⁻² fit over a ≥ 105 and the constant over a ≥ 315,
combined with the fit's own error.

| quantity | benchmark (c(M_Y) = 0) | portals zero at M_P (c = −13.31/4π²) |
|---|---:|---:|
| N_X frozen | **−8047.84 ± 0.22** | −10603.762 ± 0.006 |
| frozen / 400–1000 window | 1.000293 | 0.999824 |
| N_Σ frozen | +1341.307 | +1767.294 |
| N_Y late | −5×10⁻⁷ (→ 0) | −1×10⁻⁶ |
| Y/Y_Q at the M1 H_i | **0.991023 ± 0.000027** (fx_holder 0.990733) | 1.305763 ± 0.000001 |
| E_X late (conformal, a = 1 at the slice) | **188737.7 ± 1.1** | 233801.8 ± 7.0 |
| E_X/N_X (release energy per charge, H_i) | 23.452 | 22.049 |
| ε = N_X/N_c(E_X) | 0.127374 | 0.142929 |
| E/E_c = ε^{−4/3} | **15.6036** | 13.3816 |
| N_YL, N_YH, N_σ (late) | 0.56642, 0.010304, 1.99×10⁻⁵ | 1.4073, 0.02334, 1.87×10⁻³ |
| R_med(e) = R_Y + (R_casc + R_Σ)(10⁻⁶/e)^{1/2} | 0.3786 + 0.0184 x^{1/2} | 0.7581 + 0.0389 x^{1/2} |

Notes on the table:
- The engine's E_X = a⁴(\|Ẋ\|² + λ\|X\|⁴) contains a non-conformal piece. Its leading cross term is
  2(aH)Re(χ̄χ′) ∝ 1/a, which explains the ±1.5% instantaneous oscillation and the −1.7×10⁻⁴ offset of the 400–1000
  mean. The late limit is the conformal energy.
- The constant c = N_c/E^{3/4} = 2√(2λ_X)(3λ_X)^{−3/4} = 6.977568 in every run. It is a definition.
- The physical content of "universal" is the scale invariance confirmed in item 2.

**Sensitivity per unit ln Q** of the induced-portal condition (±1/4π² runs to 2×10⁴) [C]:
- d ln N_X = −1.703% (skeptic: −1.70%);
- d(E/E_c) = +0.143;
- dR_med(10⁻⁶) = −0.0274;
- d ln H_i = +1.14%;
- d(S/Q at 10⁻⁶) = +0.045.

**Systematic beyond the homogeneous model [E].** The M1 condensate fragments at H_i t ≈ 60–70 [Q REVIEW_TIMING; not
recomputed for M1r]. Relative to the frozen value, N_X at that epoch is:
- instantaneous: +1.0% to −2.2% over H_i t ∈ [50, 200];
- 2-orbit running mean: −1.5% to +0.7%;
- 4-orbit running mean: −0.47% to +0.17%.

How this enters depends on the fragmentation:
- If the exchange stays adiabatic through fragmentation, its dressing is returned as a⁻² and the homogeneous frozen
  value holds.
- If fragmentation is non-adiabatic on the orbit time and dephases the exchange (X⁶Σ* summed over patches), N_X keeps
  the dressing it had then.

So the "0.2%" target of WORKING_THEORY §8 is met inside the homogeneous model (2.7×10⁻⁵ dynamical; +0.03% total
late drift). The fragmentation coupling of the exchange is the remaining uncertainty: plausibly −0.5% … +0.2% in N_X,
up to −1.5% … +0.7%, i.e. H_i within about −0.5% … +1%. Thermal wash-out of the X charge through the vev-induced
3X ↔ 3X̄ operator (coefficient κ²μh/M_Y⁴ ≈ 4×10⁻⁹ H_i⁻²) has Γ/H ~ 10⁻¹⁵ a⁻³ [E], so it is negligible.

## 2. S/Q(R) of the relaxed Hartree plasma for this release [C]

The setup:
- circular condensate plus Hartree X, X̄, φ gas at fixed N_X = 8047.84 and energy (1 + R) E_X, E_X = 188737.7;
- φ in equilibrium;
- onset where g²(r_c² + D) + 3λ_φP = λ_φv²a²/H_i².

| R | 0 | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 | 3.5 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S/Q | 1.4193 | 1.9575 | 2.4498 | 2.9109 | 3.3489 | 3.7685 | 4.1731 | 4.5650 | 4.9460 | 5.6805 | 6.3846 | 7.7215 |
| X̄ per charge | 0.115 | 0.161 | 0.203 | 0.242 | 0.280 | 0.316 | 0.351 | 0.385 | 0.418 | 0.482 | 0.543 | 0.660 |
| condensed fraction | 0.966 | 0.957 | 0.950 | 0.943 | 0.937 | 0.932 | 0.927 | 0.922 | 0.917 | 0.909 | 0.901 | 0.886 |
| T_h/m_X at onset | 5.569 | 6.177 | 6.638 | 7.013 | 7.332 | 7.610 | 7.858 | 8.082 | 8.286 | 8.650 | 8.966 | 9.502 |
| n_Q at onset (m_X³) | 159.24 | 157.76 | 156.50 | 155.38 | 154.37 | 153.45 | 152.58 | 151.77 | 151.00 | 149.58 | 148.27 | 145.92 |

The full 0.25 grid, with increments, is in `summary.json` under `benchmark.SQ_of_R_table`.
- a_on = 5.17×10¹² (R = 0) to 5.33×10¹² (R = 8), in rpdq's normalisation (a = 1 at the slice) at H_i = 1.293×10¹¹
  GeV. It scales ∝ H_i, so the Q-matter H_i of 1.51×10¹¹ gives ×1.17.
- A self-consistent onset, keeping −κa² inside the φ Hartree mass, moves a_on by +0.6–1.0%, T_h/m_X by −0.7–1.0%
  and S/Q by +0.1–0.3%.

**Universality [C].** The direct M1r solve and the spot-normalised solver (kin_source release, scaled with the same
E/E_c) agree at R = 0, 1, 2.5, 5, 8 to:
- |ΔS/S| ≤ 1.9×10⁻¹⁵;
- |ΔT/T| ≤ 9×10⁻¹⁶;
- |Δn/n| ≤ 2×10⁻¹⁵.

The relaxed state therefore depends only on W = (E/E_c)(1 + R). Closed form: S/Q = 0.190255 (W − 0.97494)^{0.75013}
− 0.00428. Its maximum residual is 2.3×10⁻⁶ on W ∈ [15.6, 140], and 6.9×10⁻⁶ on the independent Planck scan
W ∈ [13.4, 120]. That is the massless-gas law S ∝ (E − E_c)^{3/4} [C/D].

**Thresholds** (survival σ_min from fx_chain skeptic; R in units of this release's E_X):

| threshold | σ_min | R_crit (Hartree) | W_crit | T_h/m_X at onset | n_Q at onset | a_on (H_i = 1.293×10¹¹) |
|---|---:|---:|---:|---:|---:|---:|
| bulk Q-matter | 3.4716 | **2.1442** | 49.061 | 7.416 | 154.10 | 5.230×10¹² |
| Q = 3×10⁷ | 3.3840 | 2.0410 | 47.451 | 7.356 | 154.30 | 5.228×10¹² |
| **Q = 2×10⁵** | 2.9284 | **1.5195** | 39.313 | 7.026 | 155.34 | 5.216×10¹² |
| binding Q = 7.13×10⁴ | 2.5531 | 1.1095 | 32.916 | 6.726 | 156.24 | 5.206×10¹² |

At e = 10⁻⁶ (reading B, R = 3.648): T_h/m_X = 8.14, n_Q = 151.5 m_X³, X̄ 0.394 per charge. That reproduces wt_critic's
8.15 / 151.5 / 0.40.

**Induced-portal condition (vanish at M_P).** The release itself changes: E/E_c = 13.3816, R_med(10⁻⁶) = 0.797. R_crit
(bulk / Q = 2×10⁵) = 2.666 / 1.938, at the same W_crit.
- H_i (free X) = 1.076×10¹¹ and ν = 0.664.
- **S/Q(10⁻⁶) = 4.01 (band) / 4.34 (broad)**: free X, margin ×1.16.
- wt_critic's 3.89 assumed the M_Y release state.
- e_crit (bulk) = 1.41×10⁻⁶ (band) / 1.66×10⁻⁶ (broad).

## 3. Entropy that S/Q must still include [C/E]

| process | mechanism | Δ(S/Q) at R_crit (bulk) | at e = 10⁻⁶ (R = 3.65) | label |
|---|---|---:|---:|---|
| φ crossover, coherent | unique tilted root; along the homogeneous branch up to the spinodal m_eff² ≥ 2.26 m_X² (T = 0) / 2.51 (thermal); H/m_X = 1.43×10⁻¹³ at T_vis ≈ 99 GeV, so H/m_eff ≈ 9.5×10⁻¹⁴; coherent excitation ≲ (H/m_eff)² | ~10⁻²⁶ | ~10⁻²⁶ | C/D |
| φ crossover, chemical lag | as φ_eq ≥ 0.07v opens vev-induced 1↔2, 2↔3 channels (Γ/H ~ 10⁴–10⁷), Δs ~ (H/Γ) per affected quantum | ≲ 5×10⁻⁴ | ≲ 5×10⁻⁴ | E |
| **spinodal decomposition (supercooled condensate's latent heat; φ rolls from φ_eq(n_sp) ≈ 0.34v to the vacuum between lumps)** | T = 0: E/Q_hom(233.22) − Ω_b = 0.82727 − 0.77586 = **0.05140 m_X** per condensed charge, released at T_h = 8.7 (bulk crossing at 0.685 a_on, spinodal at 0.852 a_on). With thermal masses (vacuum-referenced bulk tangent, frozen at the crossing): 0.0361, n_sp = 205, n_b = 378 | **0.0055** (T = 0) / 0.0041 (thermal) | 0.0049 / 0.0033 | C (T = 0) / E |
| surface (Gibbs–Thomson) free energy of the first lumps | Q_frag = 4×10¹¹–1.2×10¹³: 3.6–11×10⁻⁴ m_X per charge, dissipated by ripening (T ~ 3–9) or at worst at T_ev ≈ 0.8 | 4×10⁻⁵ – 1.3×10⁻³ | same | E |
| evaporation | quasi-static two-phase evaporation is reversible (the vapour is held at saturation, μ_gas = ω); the irreversible lag per evaporated charge is −ln(z/z_sat) ~ H/Γ_ex ≤ 10⁻⁵ [Q Γ_ex/H ≥ 10⁵ fx_formation] | ≤ 10⁻⁵ | ≤ 10⁻⁵ | E |
| X̄ burn-out (2Ω or m_X + Ω per absorbed X̄) | two-phase era: absorption keeps μ_X̄ = −ω = −μ_X, which is reversible. Free-X branch: no objects below T_ev, so the gas-frozen p is never absorbed. Q-matter branch: capture keeps equilibrium to T_f ≈ 0.005–0.015 m_X, where n_X̄,eq ~ e^{−120} | **0** | **0** | D/A |
| (hypothetical burn-out of a gas-frozen remnant) | (m_X + Ω_b)p/T: p = 5.6×10⁻⁴ gives 0.010 / 0.020 / 0.099 at T = 0.1 / 0.05 / 0.01; p = 1.3×10⁻³ gives 0.023 / 0.046 / 0.23 | not realised in M1r | – | D |

**Total:** S/Q_final = S/Q_Hartree + 0.004–0.007, i.e. +0.1–0.2%.
- R_crit moves by −0.008 (T = 0 + surface upper) / −0.006 (thermal).
- e_crit moves by about +0.4%.
- The one uncomputed piece that could matter at the percent level is the chemical lag if only the pre-onset massless
  channels acted (Γ/H ≈ 10–10⁴ [Q REVIEW_TRANSPORT]). That would give Δ(S/Q) up to ≈ 0.1 at Γ/H = 10. The
  vev-opened channels make it ≲ 5×10⁻⁴ [E].

## 4. S/Q(e) and the branch [C on E inputs]

Two readings of the provisional heat (x = 10⁻⁶/e; S/Q includes the item-3 increments, T = 0 upper):
- **(A) literal:** R = 0.8056 + 4.0538x + 0.10716x^{1/2} plugged in as R of the M1r release. This keeps the spot's Y_L
  debris (0.8056, against M1r's own 0.3786) and its vector heat per spot E_X. It is an upper reading.
- **(B) transferred:** R = 0.3786 + 0.0184x^{1/2} + ν (13.0412/15.6036)(4.0538 or 4.7905)x. ν = (H_i/H_spot)²
  (Y_spot/Y)^{4/3} = 0.960 [E wt_critic]. ν is independent of the relic's mass per charge (both factors go as
  M^{∓4/3}), so it is the same in both branches.

| e | (A) R / S/Q / branch | (B band) R / S/Q / branch | (B broad) R / S/Q / branch |
|---|---|---|---|
| 1×10⁻⁶ | 4.967 / 5.66 / free X | 3.648 / **4.685** / **free X** (×1.35) | 4.239 / 5.13 / free X |
| 2×10⁻⁶ | 2.908 / 4.11 / free X | 2.017 / **3.370** / **Q-matter** (bulk; not Q = 2×10⁵) | 2.313 / 3.62 / free X |
| 3×10⁻⁶ | 2.219 / 3.54 / free X (2% margin) | 1.473 / **2.894** / **Q-matter** (Q ≥ 2×10⁵ survive) | 1.670 / 3.07 / Q-matter |
| 5×10⁻⁶ | 1.664 / 3.06 / Q-matter | 1.037 / **2.493** / **Q-matter** (binding objects too) | 1.155 / 2.60 / Q-matter |

| e_crit (S/Q crosses the threshold) | (A) | (B band) | (B broad) | (B band), portals 0 at M_P |
|---|---:|---:|---:|---:|
| bulk (3.4716) | 3.19×10⁻⁶ | **1.86×10⁻⁶** | 2.20×10⁻⁶ | 1.41×10⁻⁶ |
| Q = 2×10⁵ (2.9284) | 6.12×10⁻⁶ | 2.90×10⁻⁶ | 3.42×10⁻⁶ | 2.29×10⁻⁶ |

**For the critic.** Take dq_heat's R(e), in units of E_X = 188737.7 of this frozen release (a = 1 at the slice, H_i
units). Then:
- read `summary.json → benchmark.SQ_of_R_table`, or use the closed form with W = 15.6036(1 + R);
- add +0.004–0.007;
- compare with 3.4716.

Equivalently, free X iff R(e) > 2.136–2.139 (bulk), and Q = 2×10⁵ objects only if R < 1.511–1.514. If dq_heat's
vectors are given in spot units, convert with R_M1r = R_med,M1r + ν(13.0412/15.6036)R_vec,spot.

**H_i ledger [D on C].** H_i = H_M1 (0.9374/(M · Y/Y_Q))^{2/3} with Y/Y_Q = 0.991023:

| case | M per charge (m_X) | H_i (GeV) |
|---|---|---:|
| free X | 1.00114 (p = 5.7×10⁻⁴ [E]) | 1.2930×10¹¹ |
| Q-matter | 0.782–0.802 [Q fx_formation] | 1.499–1.524×10¹¹ |
| Q = 2×10⁵ convention | 0.9374 | 1.35094×10¹¹ (fx_holder: 1.35120) |
| portals zero at M_P, free X | – | 1.0764×10¹¹ (−16.8%) |
| portals zero at M_P, Q-matter | – | 1.247–1.268×10¹¹ |

## Limits

- The release is the homogeneous, classical rptri model. The frozen N_X is exact within it. Fragmentation-time
  dephasing of the exchange is the dominant systematic (§1, [E]). It is not computed: that needs a lattice with Y and Σ.
- S/Q is the relaxed Hartree value, which assumes thermalisation and chemical relaxation before onset [Q]. The vector
  heat is not computed here: reading B is wt_critic's transfer [E]; dq_heat's R(e) supersedes it.
- The thermal-mass latent heat freezes the onset fluctuations and scales them as a⁻² [E]. The condensed charge is
  taken ∝ a⁻³ through the crossing (as in fx_chain).
- p(S/Q) and the Q-matter mass per charge are quoted [E]. The g*s ≈ 96 used in the clock is approximate [R].
- The REVIEW_TIMING fragmentation epoch is M1's, not recomputed for M1r.

## Reproduction (OMP_NUM_THREADS=1; about 6 min, 2 processes at most)

```
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
C=research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_entropy/code; cd $C
python3 -B long_run.py base_1e5 0 100000 & python3 -B long_run.py planck_1e5 -13.31 100000 & wait      # 2 min (builds rpdq)
python3 -B long_run.py base_2e4_rtol12 0 20000 1e-12 & python3 -B long_run.py cp_plus1_2e4 1 20000 & wait  # 70 s
python3 -B long_run.py cp_minus1_2e4 -1 20000                                                           # 30 s
python3 -B freeze_analysis.py ../results/freeze.json base_1e5 planck_1e5 base_2e4_rtol12 cp_plus1_2e4 cp_minus1_2e4
python3 -B -W ignore sq_scan.py ../results/freeze.json ../results/sq_scan.json                          # 25 s
python3 -B -W ignore entropy_increments.py ../results/sq_scan.json ../results/entropy_increments.json    # 12 s
python3 -B -W ignore summary.py ../results                                                              # 5 s
```
Scratch: `$SP/wt2/dq_entropy` (binary in `build/`, series in `runs/`; DQ_SCRATCH overrides).

Credits: spinodal decomposition and Q-ball formation (Kusenko–Shaposhnikov; Kasuya–Kawasaki, hep-ph/9909509), Q-ball
evaporation thermodynamics (Laine–Shaposhnikov 1998), Bose thermodynamics (Kolb–Turner) [R].
