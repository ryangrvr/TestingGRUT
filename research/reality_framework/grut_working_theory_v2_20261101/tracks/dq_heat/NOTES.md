# Track dq_heat: the hidden heat budget R(e) of the actual M1r release

2 November 2026 · Claude, decisive pass (BRIEF_DECISIVE.md of this directory; issue #13 context) · same-session AI work
on one supplied model, not outside replication. Nothing here is evidence that X, Y, Σ, the hidden vector or the
inflaton coupling exist, and nothing is a validation or an exclusion. M1r and every input below (e, the inflaton
history, the portal renormalisation condition, the static preparation slice) are **model choices**; none follows
from closure. Labels: [D] derived, [C] computed with our own code, [E] estimate, [Q] quoted from a reviewed track,
[A] argued, [R] recalled.

**Scope.** This track computes R(e) for the release that M1r actually has:
- λ_Y = λ_Σ = 0.5, the inflaton-trilinear holder, the (3, 0.5) inflaton;
- prepared state (90.800, 6.888, 28.576) H_i.

Nothing is transferred from the old instant-release spot point (λ = 1, no holder). R is the deposited conformal energy
in units of the X release energy E_X of the same run.

## 0. Result

| e | R_vec | R_Y,hom | R_Y,k | R_cascade | **R_total** | S/Q (copied Hartree) |
|---|---:|---:|---:|---:|---:|---:|
| 1×10⁻⁶ | 3.608 | 0.379 | 4×10⁻⁴ | 0.0186 | **4.006** | 4.95 |
| 1.5×10⁻⁶ | 2.401 | 0.379 | 4×10⁻⁴ | 0.0152 | **2.796** | 4.01 |
| 2×10⁻⁶ | 1.799 | 0.379 | 4×10⁻⁴ | 0.0131 | **2.191** | 3.51 |
| 3×10⁻⁶ | 1.196 | 0.379 | 4×10⁻⁴ | 0.0107 | **1.586** | 2.99 |
| 5×10⁻⁶ | 0.715 | 0.379 | 4×10⁻⁴ | 0.0083 | **1.103** | 2.55 |
| 1×10⁻⁵ | 0.356 | 0.379 | 4×10⁻⁴ | 0.0059 | **0.741** | 2.20 |

All [C]; the S/Q column is [C] on these R. Benchmark: induced portals vanish at M_Y.

Reading the table:
- R_vec ∝ e^{−1.007} [C] and R_cascade ∝ e^{−1/2} [D].
- The objects' survival threshold S/Q = 3.4716 (bulk) is crossed at **e_crit = 2.06×10⁻⁶**. For Q = 2×10⁵ it is
  3.20×10⁻⁶; for binding objects 4.96×10⁻⁶.
- At e = 10⁻⁶ the release is in the free-X branch, with a margin of ×1.43.
- With the portals vanishing at M_P: R_total(10⁻⁶) = 3.548, S/Q = 4.09 and e_crit = 1.47×10⁻⁶.

Every source finishes before the φ onset: by 895–2.8×10⁴× in a (vectors, 1 − 10⁻⁶ decayed) and by 2.2×10⁸× (Y_L), at
H_i = 1.29×10¹¹ GeV. At 1.51×10¹¹ these factors are ×1.17 larger.

**The largest open sensitivity is a supplied input [C on the static family].** The vector abundance is set by the
hidden Higgs norm when the modes near k_* left the horizon, 2.7–3.8 e-folds before the slice. If the inflaton was
larger then (s = φ_I/φ_I(0) > 1 at fixed H), the X condensate is larger and N_A falls as (F(1)/F(s))². At e = 10⁻⁶
**the free-X verdict holds only for s < 1.50**. M1r as specified holds the slice state static (s = 1).

## Code (own; `code/`)

| File | What | Reuse |
|---|---|---|
| `rptri_dq.c` | the M1r engine | **copy of fx_holder/code/rptri.c**. Changes: (a) the forced-Y conj bug fixed: c = 2 conj(u)Σ, its derivative, the rotation and the (arg u)′ frame term (fx_holder SKEPTIC N6); (b) a dense "D" state line. The dynamics are untouched |
| `rpmodes.c` | `rptri_dq.c` plus linear finite-k fluctuation modes of (Y, Σ) or (X, Y, Σ), integrated in the same DP5(4) state | mode block new |
| `rptri_cp.c` | `rptri_dq.c` plus the induced portals c(\|u\|² − μ²)(\|Σ\|² + 2\|Y\|²); c from the environment | term as fx_holder SKEPTIC N1 |
| `dq_model.py` | model, preparation and driver | **copy of fx_holder/code/fx_common.py**; paths changed and D lines parsed |
| `background.py`, `planck_bc.py` | backgrounds (dense 0–200 and coarse 0–2000), validation | – |
| `vector_modes.py`, `vector_band.py` | longitudinal Proca modes on the recorded background | `evolve_mode` of **hw_gates/code/gate1_vector.py**, generalised |
| `sep_universe.py` | piece B, separate-universe Goldstone angle | method of **hw_gates/code/gate1b_pieceB.py**, with exact BD amplitudes |
| `ymodes.py` | finite-k mediator production (drives `rpmodes.c`) | – |
| `heat_table.py`, `finf_family.py`, `finf_sq.py`, `summary.py` | ledger, timing, S/Q, sensitivities | `hartree_copy.py` is an **unchanged copy of wt_critic/code/hartree_copy.py** (itself fx_chain skeptic hartree_indep.py) |

Results are in `results/*.json`; `summary.json` holds the headline numbers.

## 1. The M1r background [C]

The fixed engine reproduces fx_holder's λ = 0.5 trilinear run **bit for bit**: N_X(400–1000) = −8045.485140,
E_X = 188706.29, ε = −0.127359, N_YL = 0.566641, N_YH = 0.010377, N_σ = 2.0×10⁻⁵ (relative differences 0). This is
expected: the fix touches only the diagnostic, and u is real at late times.
- Checks: Gauss 6.8×10⁻¹², H 6.3×10⁻¹², work 3.4×10⁻¹³.
- The early N_YL with the fixed diagnostic (0.29, 0.47, 0.63 at τ = 0.5, 1, 2) differs from the old one. At τ = 0 it is
  still not meaningful (17.9), because the diagnostic's forced point omits the λ_Y\|Y\|² term while \|Y\| = 6.9.
- Prepared point: Hessian (0, 10.19, 18.76, 471.2, 1536.3, 1665.7) H_i². F_inf = (2Σq_i²\|z_i\|²)^{1/2} = **275.93 H_i**,
  of which X carries 21.7% of F², against F_vac = 6h = 240.
- F(t): 272.2 at τ = 0.25, 257.6 at 0.5, 240.04 at 1 (X near periapsis, \|X\| = 2.98), 245.6 at 1.5. The X part of F²
  then decays as a⁻². \|Σ\| = h/√2 to 10⁻⁶ after τ ≈ 3.
- **Radiation constant K = a²H = 1.98581** (R_r = 11.8303). The oscillating-inflaton era with Γ_I = 0.5 makes a²H grow
  from 1 to 1.986 before radiation domination. The spot point was an instant release, with K = 1.

## 2. Longitudinal hidden vectors (item 1) [C]

**Method.** These are the longitudinal Proca modes of gate 1, with m(η) = eF(η) on the recorded background:
- de Sitter during inflation with F_inf;
- the run's a(η) and a²F² as cubic splines up to τ = 2000;
- then exact radiation with the smoothed X part.

Settings:
- 61 modes over k = 10⁻⁵–30 H_i, from an adiabatic start at k/aH = 1000;
- read-out with the conformal number of v = zB, at 20 and 30 a_* (a_* = (K/m_vac)^{1/2}), or one unit of η after the run
  for k ≥ 0.3;
- the (B, Π) number oscillates by ±2–3% from a (z′/z)/ω cross term, while the v number is conserved to 10⁻⁵.

**Validation.** With the instant background and the spot point's F_inf = 266.06, the solver reproduces gate 1's n_k at
k = 1.2×10⁻³–1 to 1×10⁻⁶. Read-outs at 20 and 30 a_* differ by ≤ 1.7×10⁻³.

| e | N_A (H_i³) | N_A·(e/10⁻⁶)^{1/2} | m_A = 240e H_i | Γ_A/H_i | ⟨a_A⟩ = (πK/2Γ_A)^{1/2} | deposit/vector (E_X) | R_vec |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10⁻⁷ | 6.444 | 2.038 | 2.4×10⁻⁵ | 1.59×10⁻²¹ | 4.43×10¹⁰ | 5.630 | 36.29 |
| 10⁻⁶ | 2.026 | 2.026 | 2.4×10⁻⁴ | 1.59×10⁻¹⁸ | 1.40×10⁹ | 1.7805 | 3.608 |
| 1.5×10⁻⁶ | 1.652 | 2.023 | 3.6×10⁻⁴ | 5.37×10⁻¹⁸ | 7.62×10⁸ | 1.4538 | 2.401 |
| 2×10⁻⁶ | 1.429 | 2.020 | 4.8×10⁻⁴ | 1.27×10⁻¹⁷ | 4.95×10⁸ | 1.2590 | 1.799 |
| 3×10⁻⁶ | 1.163 | 2.015 | 7.2×10⁻⁴ | 4.30×10⁻¹⁷ | 2.69×10⁸ | 1.0280 | 1.196 |
| 5×10⁻⁶ | 0.898 | 2.008 | 1.2×10⁻³ | 1.99×10⁻¹⁶ | 1.25×10⁸ | 0.7963 | 0.715 |
| 10⁻⁵ | 0.631 | 1.996 | 2.4×10⁻³ | 1.59×10⁻¹⁵ | 4.43×10⁷ | 0.5630 | 0.356 |

Γ_A = e²m_Aβ³/48π with β = 1 to 10⁻³⁰. In GeV, at H_i = 1.29×10¹¹: m_A = 3.10×10⁷ GeV × (e/10⁻⁶) and Γ_A = 2.05×10⁻⁷ GeV
× (e/10⁻⁶)³. There are 2.5×10⁻⁴ vectors per net charge at 10⁻⁶.

**The scaling question, resolved [D, C].**
- Graham–Mardon–Rajendran give Ω_A ∝ m^{1/2}H_I² [R], so n_A/s = Ω_Aρ_c/(m s) ∝ m^{−1/2}. In comoving terms the relic is
  set where H = m, at a_* = (K/m)^{1/2}: N_A ~ a_*³ m (H/2πF_inf)² F_vac² ∝ K^{3/2} m^{−1/2}.
- The deposit per vector is m_A⟨a_A⟩, with ⟨a_A⟩ ∝ Γ_A^{−1/2} ∝ e^{−3/2}, so it scales as e^{−1/2}.
- Hence **N_A ∝ e^{−1/2} and R_vec ∝ e^{−1}**.
- Computed: d ln N_A/d ln e = −0.507, d ln(deposit)/d ln e = −0.500, d ln R_vec/d ln e = −1.007 over 10⁻⁶–10⁻⁵;
  R(10⁻⁷)/R(10⁻⁶) = 10.06.
- **gate1_vector.py's docstring ("N_A ∝ e^{1/2}", "R_A e-independent") is wrong.** Its own table (2.469, 0.771, 0.235 at
  10⁻⁷, 10⁻⁶, 10⁻⁵) shows e^{−1/2}, and its NOTES already corrected R_A to ∝ 1/e. kin_source and fx_chain use 1/e, which
  is right.

**Why M1r differs from the spot point [C].** Two factors.
- **K.** At e = 10⁻⁶ the actual reheating gives 2.724× the instant-release N_A (K^{3/2} = 2.798) and √K = 1.41× the
  decay scale factor. R_vec therefore carries ≈ K² = 3.9 relative to an instant release with the same E_X.
- **The frozen-angle factor (F_vac/F_inf)².** F_inf is 275.9 here against 266.1 at the spot, because the M1r condensate
  is larger (\|X\| = 90.8 against 80.2), which gives 0.93. The no-X control (X removed from F everywhere) gives a k ≤ 1
  number 1.283× larger, against the predicted 1.276.
- Against the GMR formula N_A is 0.70 (gate 1: 0.76 at the spot).
- In release-invariant units W_vec = |ε|^{−4/3}R_vec = 56.3 at 10⁻⁶. wt_critic's transfer gave 50.7 (band) and 60.0
  (broad). **The direct value is 0.94 of the broad transfer**; the transfer missed the X part of F_inf. The §8 failure
  condition (vector heat < 0.5 of the transfer) is not met at s = 1.

**Piece B, release-era production [C].**
1. *Sub-horizon modes (k ≥ 1)* are pumped by the oscillating inflaton (a″/a) and by the X-orbit modulation of F.
   - A resolved band (145 modes, k = 1–30) gives N_A(k ≥ 1) = 0.0512, independent of e (Goldstone regime).
   - It is dominated by a resonance at k = 4.45 (n_k = 0.042) near half the radial frequency of the eccentric X orbit.
   - The no-X profile gives 0.0236 (no resonance).
   - This is 2.5% of R_vec at 10⁻⁶ and 8% at 10⁻⁵, and it scales as e^{−1/2}. The full-versus-no-X difference
     (±0.05 in R at 10⁻⁶) is the model uncertainty of the single-Higgs-direction approximation [A].
2. *Separate-universe Goldstone angle* (gate1b method).
   - Responses dα_G/dq along the five heavy inflationary eigenmodes: −8.8×10⁻⁴, 3.35×10⁻³, −3.2×10⁻⁵, 3.6×10⁻⁶,
     −1.2×10⁻⁶ rad per H_i, linear between amplitudes 10⁻³ and 10⁻². The total Σ winding is −0.246 rad.
   - With exact de Sitter BD superhorizon powers: N_A,B = 3×10⁻⁴, i.e. **1.5–4.3×10⁻⁴ of piece A**. Negligible.

**Portals vanish at M_P [C].** Release from `rptri_cp.c`:
- N_X = −10606.160 (skeptic's independent code: −10606.160, 2.7×10⁻¹⁰), E_X = 233872, ε = 0.1429, N_YL = 1.408,
  N_YH = 0.0234, N_σ = 1.9×10⁻³, F_inf = 283.85.
- N_A(10⁻⁶) = 1.914, which is ×0.945 = (275.93/283.85)².
- R_vec(10⁻⁶) = 2.751. It is lower because E_X is larger, while R_Y doubles to 0.758.
- The c = 0 control reproduces rptri_dq exactly.

## 3. Mediator debris (item 2) [C]

**Homogeneous.** N_YL = 0.56664 and N_YH = 0.010377 (late, fixed diagnostic). The Y_H cascade first, so the Y_L parents
number N_YL + N_YH:
- m_L = 22.249 H_i, Γ_L = 3κ²m_L/256π³ = 1.006×10⁻⁷ H_i, ⟨a_L⟩ = 5569 with K;
- the deposit per Y_L is 0.6565 E_X;
- **R_Y,hom = 0.3788** (fx_holder's ledger: 0.37884).

**Finite-k Y production.** `rpmodes.c` integrates δr̈ + 3Hδṙ + (k²/a² + Hess V(t))δr = 0 for the canonical (Y, Σ)
components ("YS", X as background) and for all of (X, Y, Σ) ("XYS"). The Hessian is taken on the instantaneous
background and includes μ_eff(t) = μ + Gφ̂_I(t), the Y and Σ backgrounds, 3κX² and ξR.

Initial state:
- every inflationary eigenmode starts in its exact de Sitter Bunch–Davies state, using Hankel functions of imaginary
  order, with Wronskian error 2×10⁻¹⁶;
- 40 values of k over 0.03–80;
- quanta are read off with the rescaled variable v = a^{3/2}δr, Ω² = k²/a² + M² − 9H²/4 − 3Ḣ/2;
- rtol 10⁻¹¹ and 10⁻¹³ agree.

Findings:
- **For p = k/a ≳ 13 the "occupation" is the adiabatic ambiguity (H/ω)²** at the early read-outs. It decays as H falls;
  it is not particle production.
- **A band at physical momentum p = 11.2 at every time is the cascade Y_H → Y_L + Goldstone**, with p* = (m_H² −
  m_L²)/2m_H = 11.21. The linear theory feeds it from the homogeneous Y_H, which does not deplete. Its rate is the tree
  vertex √2μ G y_H y_L, Γ_H = μ²(m_H² − m_L²)/(8πm_H³) = 0.0350 H_i, which is kin_source's Goldstone-equivalence width,
  here exact. The k ≤ 80 plateau 5.1×10⁻³ is mostly these products born before a = 80/p* = 7.1 (4.5×10⁻³). They are
  already counted through N_YH and are not counted again.
- **Genuine release production of Y_L**: 0.94×10⁻⁴ (YS, k < 15, cascade-free); 2.6×10⁻⁴ (XYS: the δX coupling adds
  ×2.8 while δX self-resonates); 5.7×10⁻⁴ (subtraction estimate). The table uses the largest, giving **R_Y,k = 3.8×10⁻⁴**,
  i.e. 0.1% of R_Y,hom. The hard upper bound, all k ≤ 80 including the double-counted cascade, is 3.4×10⁻³.
- Finite-k Y_H and σ: N ≈ 5×10⁻⁶ and 3×10⁻⁹ at k < 15; above that it is noise at ±2×10⁻⁵.
- **So the holder's time-dependent Y mass (330 ↔ 534 H_i², at the inflaton frequency) leaves only ≤ 10⁻³ of the
  homogeneous Y_L in finite-k quanta** (ω_Y ≈ 20–23 ≫ m_I = 3). fx_holder's "adiabatic" guess holds.

**Cascade and σ vectors.**
- Y_H → Y_L + A_L at ⟨a_H⟩ = 9.4 (Γ_H = H at a = 7.5). The A_L carries E = 11.2 H_i, is non-relativistic well before it
  decays, and decays with the primordial vectors.
- σ → A_L A_L, with Γ_σ = λ_Σ²h²/8πm_σ = 0.40 H_i (σ → Y_LY_L is closed: 40 < 44.5), is prompt.
- R_cascade = (N_YH + 2N_σ)·m_A⟨a_A⟩/E_X = 0.0186 (10⁻⁶/e)^{1/2}. It scales as e^{−1/2} at fixed N_YH [D], as the
  fx_chain skeptic found.

## 4. Timing relative to the φ onset [C]

The onset comes from the copied relaxed-Hartree solver at the computed W = |ε|^{−4/3}(1 + R_total), converted to this
frame by the release-state invariant: a_on = a_on,spot (N_X/N_spot)^{1/3}(H_i/H_spot). Check: at e = 10⁻⁶ the onset
density n_Q = 151.0 m_X³ gives a_on = (H_i/m_X)(N_X/151)^{1/3} = 5.253×10¹² directly. **In this frame a_on ≈ 5.2×10¹²
(H_i = 1.29×10¹¹) and 6.1×10¹² (1.51×10¹¹).** The brief's 3.8×10¹² is the spot frame (N_X = 2788, K = 1, H_i = 1.3335×10¹¹).

| source | a at emission/decay | T_vis there (1.29 / 1.51×10¹¹) | a_on / a_source at H_i = 1.29×10¹¹ (a_source: 1 − 10⁻⁶ completion for Y_L and vectors; Γ = H for σ; ⟨a⟩ for Y_H) | same, 1.51×10¹¹ |
|---|---|---|---|---|
| σ → A_LA_L | ≈ 2.2 | – | 2.4×10¹² | 2.8×10¹² |
| Y_H → Y_L A_L | ⟨a⟩ = 9.4 | 4.5 / 4.9×10¹³ GeV | 5.6×10¹¹ (⟨a⟩) | 6.5×10¹¹ |
| Y_L → 3X | ⟨a⟩ = 5.57×10³, done 2.34×10⁴ | 7.7 / 8.3×10¹⁰ GeV | 2.2×10⁸ | 2.6×10⁸ |
| vectors, e = 10⁻⁶ | ⟨a⟩ = 1.40×10⁹, done 5.87×10⁹ | 3.0 / 3.3×10⁵ GeV | 895 | 1050 |
| 1.5×10⁻⁶ | 7.62×10⁸, 3.20×10⁹ | 5.6 / 6.1×10⁵ | 1640 | 1920 |
| 2×10⁻⁶ | 4.95×10⁸, 2.08×10⁹ | 8.6 / 9.3×10⁵ | 2510 | 2940 |
| 3×10⁻⁶ | 2.69×10⁸, 1.13×10⁹ | 1.6 / 1.7×10⁶ | 4610 | 5390 |
| 5×10⁻⁶ | 1.25×10⁸, 5.25×10⁸ | 3.4 / 3.7×10⁶ | 9890 | 1.16×10⁴ |
| 10⁻⁵ | 4.43×10⁷, 1.86×10⁸ | 9.6 / 10.4×10⁶ | 2.79×10⁴ | 3.27×10⁴ |
| φ onset | 5.25 / 6.15×10¹² | 81 / 75 GeV (g* = 106.75 held) | – | – |

Cascade and σ vectors decay with the primordial ones. Every source completes at least 895× before the onset. The
Hartree pre-onset state at e = 10⁻⁶ has T_h = 8.29 m_X, n_Q = 151 m_X³ and X̄ = 0.42 per charge. At 10⁻⁵ these are
6.41 m_X, 157 m_X³ and 0.18.

## 5. Sensitivities and the decision

1. **Inflaton value at the horizon exit of k_* [C on the family, A on its relevance].** On the prepared family at fixed
   H (`finf_family.py`), F_inf = 275.9, 321.7, 402.9, 578.2 at s = 1, 1.25, 1.5, 2. The light Y softens, so the X
   condensate grows. The k ≤ 1 vector number scales as (F(1)/F(s))² = 1, 0.736, 0.469, 0.228.

   | s | S/Q(10⁻⁶) | e_crit (bulk) |
   |---|---:|---:|
   | 1 | 4.95 | 2.06×10⁻⁶ |
   | 1.25 | 4.24 | 1.53×10⁻⁶ |
   | 1.5 | 3.47 (threshold) | 1.0×10⁻⁶ |
   | 2 | 2.72 | 6.1×10⁻⁷ |

   The engine represents the earlier inflation by the static slice state, so M1r as specified has s = 1. A real
   inflaton history, which M1r does not supply (m_I = 3H_i is not a slow-roll field), decides this. If H also grew with
   φ_I, as in large-field models, the factor would be ≈ s²(F(1)/F(s))² ≈ 1.06–1.15 [E]: nearly cancelling.
2. **Portal renormalisation condition [C].** Vanishing at M_P instead of M_Y gives R_total(10⁻⁶) 3.548 instead of 4.006,
   S/Q 4.09 instead of 4.95, and e_crit 1.47×10⁻⁶ instead of 2.06×10⁻⁶.
3. **Single-Higgs-direction model of the sub-horizon vectors [A]:** ±0.05 in R at 10⁻⁶ (full against no-X).
4. **Finite-k Y:** ≤ 3.4×10⁻³ in R even at the hard upper bound.
5. Not included: the X ↔ Σ exchange to freeze-out, and nonlinear (post-fragmentation) Y or vector production. XXX → Y_L
   is kinematically closed once the X quanta have p ≲ m_L/3 [A].

**Decision for the branch (item 1 of BRIEF_DECISIVE):**
- **At the declared e = 10⁻⁶ and s = 1, S/Q = 4.95 > 3.47: free X** [C].
- The Q-matter branch needs e > 2.06×10⁻⁶ (bulk) or > 3.2×10⁻⁶ (Q = 2×10⁵). With M_P portals the bulk bound is
  1.47×10⁻⁶. Alternatively it needs a horizon-exit inflaton ≥ 1.5 φ_I(0) at e = 10⁻⁶.

## 6. Limits

- The calculation is homogeneous and classical for the background, and linear for the fluctuations. Metric and inflaton
  fluctuations are dropped; their mixing is ≤ 4×10⁻⁴ H_i² and η² = 3×10⁻¹⁵.
- The longitudinal mode uses the single-Higgs-direction Proca action (as gate 1). The multi-field eaten-Goldstone mixing
  with δX is not followed.
- Inflation is de Sitter at H_i with the static slice state (see §5.1).
- g* = 106.75 is held through every decay (T_vis > 10⁴ GeV) and at the onset (81 GeV; the true g* there is about 5–15%
  lower) [E].
- The deposit ledger assumes A → XX̄ and Y_L → 3X with branching 1, decay at rest, and full conformal energy into the X
  plasma [Q kin_source].
- The S/Q column uses the copied relaxed-Hartree solver at W (the wt_critic convention). The end-to-end S/Q is the
  dq_entropy track's job.
- Process note: one vector-mode batch ran past the 10-minute call limit and was moved to the background by the harness
  (2 workers). It finished before any further work. One short tolerance check ran two processes concurrently. Neither
  exceeded 2 processes.

## Reproduce (one or two cores)

```bash
export OMP_NUM_THREADS=1; cd research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_heat/code
python3 -B background.py                       # 6 s, builds rptri_dq in $SP/wt2/dq_heat/build
for s in "full:1e-6 full:1e-5" "full:1e-7 full:1.5e-6 full:2e-6" "full:3e-6 full:5e-6 noX:1e-6 noX:1e-5" "inst:1e-7 inst:1e-6"; do
  python3 -B vector_modes.py run $s; done        # about 50 s per set (inst 1e-7: 200 s)
python3 -B vector_modes.py assemble && python3 -B vector_band.py      # 4.5 min
python3 -B planck_bc.py && for s in "full:1e-6 full:1.5e-6 full:2e-6" "full:3e-6 full:5e-6 full:1e-5"; do DQ_BG=planck python3 -B vector_modes.py run $s; done
DQ_BG=planck python3 -B vector_modes.py assemble
python3 -B ymodes.py YS 100 YS_r11 1e-11 80 && python3 -B ymodes.py XYS 40 XYS_r11 1e-11 80 0.03 40   # 20 s each
python3 -B ymodes.py YS 30 uv_r11 1e-11 400 20 12 && python3 -B ymodes.py YS 30 uv_r13 1e-13 400 20 12   # tolerance check
python3 -B sep_universe.py && python3 -B heat_table.py && python3 -B finf_family.py && python3 -B finf_sq.py && python3 -B summary.py
```

Credits: longitudinal-vector production during inflation (Graham–Mardon–Rajendran, arXiv:1504.02102) [R]; de Sitter
Bunch–Davies modes of massive fields (Hankel functions) [R]; Goldstone equivalence [R]. Our own tracks are reused as
stated above.
