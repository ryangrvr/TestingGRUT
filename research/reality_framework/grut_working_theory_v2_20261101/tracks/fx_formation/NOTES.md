# Track fx_formation: charge fraction in objects (η) and the object distribution under M1 (repair 3)

1 November 2026 · Claude, issue #13 (BRIEF.md of grut_working_theory_v2_20261101) · same-session AI work on one supplied
model; not outside replication, not a validation or an exclusion. Labels: [D] derived, [C] computed with our code,
[E] estimate, [A] argued, [R] recalled, [Q] quoted. **A repair is a model choice, not a consequence of closure.**
Git read-only; nothing committed.

## Result in one paragraph

Equilibrium holds (exchange rates exceed H by 10⁵–10¹¹ until T_h ≈ 0.005 m_X). That makes thermodynamics, not the
formation ensemble, set η, and one number decides it: the hidden entropy per net charge S/Q. M1's hidden sector has no
light neutral state (m_φ = √6 m_X; the vector and Y, Σ are heavy). So the latent heat of condensation can only be
carried by the X vapour. The vapour in coexistence with Q-matter holds at most n/s = **0.288** (at T_h = 0.336 m_X)
[C]. **If S/Q > 3.47 the condensed phase evaporates completely between T_ev and T_s and cannot re-form.** The e = 10⁻⁶
band, which gives p ≈ 0.005, has S/Q = 4.95 (heavy release) or 3.61 (pilot release) [E, from fx_chain's heat budget].
There η → 0, the relic is free X plus pairs at 1.0045 m_X per charge (+7.2% against 0.9374), and the frozen test reads
0.072 (×7). With less heat (R ≲ 2.8, e ≳ 2.2×10⁻⁶ on the heavy release) the objects survive, with η = 0.89–0.97 (set
by an entropy floor) and p ≈ 0. But the slow passage and ripening make **⟨Q⟩ ≈ 10¹²–10²⁰, not 2×10⁵**. Then
⟨E/Q⟩ = 0.7759–0.7763 m_X, and the relic's mass per charge is 0.782–0.802 m_X: **−14.5% to −16.6%** against the
calibration's 0.9374. Neither branch reproduces the reference object to 1%. η ≥ 0.92 is plausible only in the low-heat
branch, and there it does not rescue the calibration.

## 1. Object family and equilibrium solitosynthesis

### 1a. E(Q), Ω(Q) from Q_min to 10⁸ [C] (`code/family.py`, 125 s)
Solver: tracks/m1_vacuum `qball_m1.py`, copied unchanged to `code/qball_m1_copy.py`. Fixed-Q continuation (h = 0.025,
L = 60–140) plus a fixed-Ω branch (L = 100, h = 0.05) for the minimum.

| Q | Ω | E/(Q m_X) | r₉₀ m_X |
|---|---|---|---|
| **Q_min = 56 709** | 0.97281 | 1.01240 | ≈ 7.5 |
| 7.1263×10⁴ (binding) | 0.94017 | 0.99999933 | 5.89 |
| 10⁵ / 2×10⁵ / 10⁶ | 0.91610 / 0.88172 / 0.83430 | 0.978960 / **0.937400** / 0.867027 | 5.62 / 5.94 / 8.48 |
| 10⁷ / 3×10⁷ | 0.80197 / 0.79379 | 0.816078 / 0.803285 | 16.9 / 24.1 |
| 5×10⁷ / 7×10⁷ / **10⁸** (new) | 0.79093 / 0.78930 / 0.78777 | 0.798849 / 0.796336 / **0.793978** | 28.6 / 31.9 / 35.9 |

Every m1_vacuum value is reproduced to its quoted digits. Fit on Q ≥ 3×10⁵: **E/Q = 0.775865 + 8.1429 Q^{−1/3} +
125.79 Q^{−2/3} − 2847 Q^{−1}**, with maximum residual 3.2×10⁻⁵. A two-term fit gives 8.369 and 70.13 (4.6×10⁻⁴). The
derivative of the fit, Ω_b + (2/3)A Q^{−1/3} + (1/3)B Q^{−2/3}, matches the solver's Ω to ≤ 2×10⁻⁵ (first law). The
surface tension from A is σ = 98.8 m_X³. Gibbs–Thomson: μ(R) = Ω_b + 2σ/(n_B R).

### 1b. Gas, liquid and coexistence (`code/thermo.py`)
- **Vapour** [C]: ideal relativistic Bose gases with the vacuum masses that hold outside objects: X (m_X, +μ), X̄
  (m_X, −μ), φ (√6 m_X, μ = 0), integrated by quadrature. Checks: massless n, p and s to 10⁻⁹; net charge μT²/3 to
  5×10⁻⁶; the NR limit with its 15T/8m correction. Thermal masses are omitted (δm²/T² ~ 10⁻³) [E].
- **Liquid**: bulk Q-matter. The zero-T pressure P₀(μ) = 100 max_{S,F}[½μ²F² − U] gives P₀ ≈ n_B(μ − Ω_b) [C]. A
  one-loop thermal pressure is added from its three excitation branches: a phonon and two gaps, 1.740 and 2.134 m_X at
  k = 0, from the cold-fluid cubic of m1_vacuum (`code/coldfluid_m1_copy.py`, changed only to guard its main block) [E].
- **Coexistence** μ_c(T): P₀ + p_liq = p_vap. μ_c = 0.775865 for T ≤ 0.3, 0.775844 at 1, 0.7744 at 4 and 0.763 at 10
  (the thermal pressures nearly cancel). The shift is irrelevant where η is decided (T ≈ 0.1–1).
- **Saha objects** [C]: n(Q) = g_int (E T/2π)^{3/2} exp[−(E − μQ)/T] summed over Q ∈ [Q_min, Q_max], with g_int =
  10⁻¹⁰, 1, 10¹⁰ (baseline g_int = 1: ground state, collective-mode partition function omitted). **η is identical to
  four digits across 10^{±10}**: the shift ln g/Q ≤ 4×10⁻⁴ T. With ω(Q) decreasing, the weight has a minimum at the
  critical object Ω(Q*) = μ and is dominated by Q_max whenever η > 0 (mean Q = Q_max in every condensed row).
  Equilibrium therefore wants the largest objects, and Q is set by kinetics.

### 1c. Equilibrium history at fixed Y_h = n_Q/s_h (`code/eq_history.py`; results/eq_history.json, eq_saha.json)
Y_h is read from fx_chain heat_onset.json: onset and spinodal-crossing values, massless 3-dof entropy [E]. The visible
clock uses s_vis = n_Q/Y_Q and an SM g*(T) [R], with a matter term below 0.8 eV. Charge in the vapour at coexistence
divided by vapour entropy, n/s at μ_c [C]: 0.027 (T = 7), 0.084 (2.1), 0.146 (1.2), 0.235 (0.65), **0.288 (0.336,
maximum)**, 0.264 (0.2), 0.157 (0.06), 0.040 (0.01). This gives Y_crit = 0.288. Finite objects coexist at μ = Ω(Q)
and raise it: 0.293 (Q = 10⁸), 0.314 (10⁶), **0.341 (2×10⁵)**.

| Case (fx_chain label) | R | Y_h (onset / spinodal) | S/Q | phases as T_h falls | T_ev / T_s (m_X) | η_eq min |
|---|---:|---|---:|---|---|---:|
| heavy R0 | 0 | 0.810 / 0.784 | 1.24 | liquid → two-phase throughout | – | 0.64 (T ≈ 0.3) |
| heavy med_e1e6 | 0.91 | 0.484 / 0.462 | 2.07 | two-phase throughout | – | 0.40 |
| **heavy e1e6 band** | 4.97 | **0.202 / 0.185** | 4.95–5.4 | two-phase → **vapour only** → two-phase | **0.82–0.94 / 0.078–0.10** | **0** |
| **pilot e1e6 band** | 4.97 | **0.277** | 3.61 | same, narrow window | **0.47 / 0.24** | **0** |
| heavy e1e7 band | 41.7 | 0.046 / 0.035 | 22–29 | vapour only from T ≈ 4.3–5.7 | 4.3–5.7 / 0.009–0.011 | 0 |
| pilot R0 | 0 | 1.139 | 0.88 | two-phase | – | 0.75 |

Equilibrium η along the e1e6 heavy history: 0.92 (T = 8), 0.66 (2.7), 0.43 (1.5), 0.07 (0.88), 0 (0.5–0.1), then 0.26
(0.055), 0.49 (0.032) and 0.79 (0.011). Finite-Q Saha (Q_max = 10⁸) gives η = ½ at T = 1.66 and 0.029; with Q_max =
10⁶, 2.10 and 0.018. **Condensation temperature T_s in the Griest–Kolb sense** [C]: 0.08–0.10 m_X for the e1e6 band.
For S/Q < 3.47 (R0 and the medium case) there is no T_s: objects are favoured at every T. **Entropy floor** [D]: below
T ≈ 0.1 all entropy sits in the vapour (s/n ≈ 2.5 + Δ/T with Δ = m_X − μ_c = 0.224 m_X), so 1 − η_eq ≈
1/[Y_h(2.5 + Δ/T)]. η reaches 1 only as T → 0, and T_h falls only logarithmically in a, because latent heat keeps the
vapour warm (T_h/T_vis ≈ 0.0096 at the spinodal, but > 1 by T_vis ~ keV).

## 2. Kinetics against equilibrium [E] (`code/kinetics.py`; results/kinetics.json)
Geometric absorption πR² is used (r₉₀ up to 10⁸, bulk sphere above; tracks/m1_transport's 332/m_X² at Q = 2×10⁵ is
3× larger), with sticking s = 1 and 0.01.
- **Exchange per object charge** Γ_ex/H: 10¹⁰–4×10¹¹ (Q = 3×10⁵), 10⁵–4×10⁶ (Q = 10²⁰) and 1–400 (Q = 10³²) for
  T_h ≥ 0.07. Capture per vapour X is similar. Equilibrium is maintained to T_h ≈ 0.005–0.01 for every Q ≲ 10²⁶.
- **Evaporation** (e1e6 band): at T_ev every object with Q < 2.2×10³⁷ (heavy) or 3.3×10³⁶ (pilot) evaporates within
  1/H. That is 10⁻⁷ of the 1.5×10⁴⁴ charges per Hubble volume, and the largest grown object is 7×10¹⁹. Evaporation is
  complete. **Re-nucleation** below T_s from the supersaturated vapour (critical object Ω(Q*) = μ, or Q_min once
  μ > 0.9728): barrier/T ≥ 8.4×10⁴ for T = 0.002–0.08, i.e. never. Fugacity e^{(μ−m)/T} = 0.08–0.32 stays far from
  BEC (ζ(3/2)). There are no small bound states (Q_min = 5.67×10⁴; the φ-exchange Yukawa α ≈ 2×10⁻⁴). This is the
  Griest–Kolb/Postma bottleneck with no way round it [R for the literature, C for the numbers].
- **Merger**: Fuchs kernel with Epstein drag in the vapour; diffusive at T ≳ 0.1 (K ≈ 10⁻³–10⁻¹ m_X⁻²), ballistic
  below (K ≈ 2–8). Mean-field coagulation runs at rate/H ≈ 1 (self-similar) and gives Q ≈ 10¹¹ at T_h ≈ 6, rising to
  10¹³–8×10¹³ at freeze-out (s_m = 1; 10¹² for s_m = 0.1).
- **Ostwald ripening** (Wagner, interface-limited, ballistic vapour; σ = 98.8): Q ≈ 6×10¹⁶–2×10¹⁸ by T_h ≈ 6–3 and
  **2×10¹⁶–7×10¹⁹** at the end (s = 0.01–1). It dominates, and it is fastest exactly where the vapour holds a large
  share of the charge (T ≈ 0.3–3).
- **Freeze-out of vapour capture** (surviving cases; relaxation ODE along the equilibrium T(a)): T_f = 0.005–0.008
  m_X (T_vis ≈ 3 keV–0.3 MeV) for ⟨Q⟩ = 10¹⁴–10²⁰. 1 − η_f is 0.028–0.040 (R0, s = 1), up to 0.066 (s = 0.01,
  Q = 10²⁰); and 0.048–0.067, up to 0.113, for R = 0.91. p ≈ 0, because X̄ are absorbed at the same rate. The residual
  free X leave warm, v ≈ 0.13–0.16c at T_f. Their comoving speed is 15–1700× that of Codex's 6.6% gas, giving
  λ_fs ~ 0.5–40 kpc for a 3–11% fraction. This is a new transfer source that needs its own calculation [E].

## 3. The spinodal route [E] (kinetics.json: spinodal_passage, cavitation)
- The onset densities are below the spinodal, as the brief says. Own estimate: n = √λ_X r³ with r² = λ_φv²/g² =
  300 m_X² gives **164.3 m_X³** (untilted; fx_chain: 164.3, M1 curvature 76.4). Under M1, however, φ never vanishes
  and the homogeneous condensate crosses n_B = 449 at T_h ≈ 12.3 and the spinodal 233.2 at T_h ≈ 9.9 (e1e6;
  fx_chain [Q]). Decomposition therefore happens at the 233 crossing, not at onset.
- Thin-wall cavitation in the stretched liquid (P from −2 to −35 m_X⁴) has barriers ≥ 1.3×10⁴ m_X against 120 T ≈ 1200
  m_X, so there is no nucleation through most of the metastable band. Thin-wall fails near the spinodal; nucleation just
  above it would only make the structures larger.
- **Slow passage**: Γ_peak = 0.774 δ m_X and k_peak = 0.958 √δ m_X near the spinodal [C]. H/m_X = 1.31×10⁻¹³ at
  n_Q = 233 (T_vis = 95 GeV), and the amplification is A = c₁δ²/6H. For A* = 10–60: δ* = 3–8×10⁻⁶, k* =
  1.7–2.7×10⁻³ m_X and Γ*/H = 2–5×10⁷, so **Q_frag ≈ 4×10¹¹–1.2×10¹³** (half to full wavelength cube; k* ∝ H^{1/4}).
  The 2.96×10⁵ cube of the n = 88 m_X³ sudden quench is not the physical regime. Kasuya–Kawasaki (hep-ph/9909509,
  hep-ph/0002285) and Enqvist–McDonald (hep-ph/9711514, 9803380) fragment with H ~ m at the start of oscillation, which
  sets the charge at the fastest-growing scale [R]. Here the passage is quasi-static, and the volume fraction at the
  crossing is 233/449 = 0.52 (bicontinuous), so hydrodynamic coarsening (not computed) can only add to the
  ripening above.

## 4. Verdict [E unless marked]

| Branch | η | ⟨Q⟩ | charge-weighted ⟨E/Q⟩ | p | mass per charge (m_X) | vs 0.9374 | frozen test (reference form) |
|---|---|---|---|---|---|---|---|
| e = 10⁻⁶ band (heavy and pilot), S/Q = 3.6–5.4 | **≈ 0** (evaporated, no re-formation) | – | – | 0.0023 (kin_solve gas-only) | **1.0045** | **+7.2%** | **0.072 (fails ×7)** |
| low heat R ≲ 2.8 heavy / 4.7 pilot (e ≳ 2.2×10⁻⁶ / 1.1×10⁻⁶) | **0.89–0.97** | **10¹²–10²⁰** | **0.7759–0.7763** | ≈ 0 | **0.782–0.802** | **−14.5…−16.6%** | 0.002–0.008 (passes as written) |
| e = 10⁻⁷ band | 0 | – | – | 0.74 | 2.47 | +164% | fails |

- **Is η ≥ 0.92 plausible?** For the case the test was built on (p ≈ 0.005, e = 10⁻⁶): **no**; the condensed phase
  does not survive T_h ≈ 0.5–0.9 m_X. In the low-heat branch η ≥ 0.92 is plausible (0.93–0.97 at R = 0), but the
  objects are macroscopic Q-matter and the mass per charge is 16% low. The supplied Q = 2×10⁵ is never the formed
  population.
- **Repairs, each a model choice:** (i) recalibrate D on the computed relic, M/Q ≈ η Ω_b + (1 − η) m_X ≈ 0.78–0.80 m_X.
  Y_Q rises by 1.17–1.20 and H_i by about ×1.11–1.13 (Y_Q ∝ H_i^{3/2}), with the vapour fraction and the warm free-X
  tail kept in the ledger. (ii) Require S/Q < 3.47 at T_h ≈ 0.34 m_X (e ≳ 2.2×10⁻⁶ on the heavy release). (iii) Add a
  light hidden neutral that can carry the latent heat (η → 1). That is outside M1 and has its own ΔN_eff and coupling
  ledger. Fixing Q = 2×10⁵ would need a size-limiting mechanism M1 does not have: the hidden vector is heavy, so there
  is no Coulomb Q_max [A].
- **Which calculation decides it.** Survival is decided by S/Q at T_h ≈ 0.34 m_X against 3.47. Any entropy produced in
  the φ crossover, the decomposition or the evaporation cycles raises S/Q. No 3D formation-ensemble run can therefore
  rescue the e = 10⁻⁶ band, and the expanding run of REVIEW_M1 ends ~10¹¹ dynamical times before η is set. The one 3D
  calculation that would decide the remaining inputs (coexistence vapour density with its Gibbs–Thomson shift, and
  sticking) is a **48³ non-expanding box, dx = 0.5/m_X (L = 24/m_X)**. It holds one Q = 2×10⁵ M1 object (r₉₉ = 8.4/m_X)
  in an X/X̄/φ bath at T_h = 0.5 m_X, with total charge and energy fixed at S/Q = 5.0 (test) and 3.0 (control). It runs
  ≥ 10⁴/m_X, beyond the exchange time ≈ 2×10³–10⁴/m_X, and measures dQ/dt and the saturated n_X at μ = Ω(Q) =
  0.8817. Our numbers stand if n_X matches the ideal-gas value to 10% and sticking ≥ 10⁻⁴. The classical-statistical
  treatment is marginal there (occupations 0.5–4), so a quantum-kinetic vapour coupled to the classical object is the
  correct version.

## Limits
Y_h comes from fx_chain's heat budget with massless 3-dof entropy [E]. Q-matter is used at tree level with a one-loop
thermal estimate, and the vapour is ideal at vacuum masses. The rates are geometric, with sticking set by hand. Merger,
ripening and freeze-out are order-of-magnitude models, not a solved coagulation equation. The bicontinuous
hydrodynamic stage and gas–core momentum coupling are not computed. Literature content is recalled [R], with no
numbers taken from it. No 3D runs were made.

## Reproduction (scratch = $SP/wt2/fx_formation; OMP_NUM_THREADS=1)
```
T=research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_formation
python3 -I $T/code/family.py     $T/results/family.json            # 125 s
python3 -I $T/code/eq_history.py $T/results/eq_history.json hist   # 35 s
python3 -I $T/code/eq_history.py $T/results/eq_saha.json    saha   # 95 s
python3 -I $T/code/kinetics.py   $T/results/kinetics.json          # 60 s
python3 -I $T/code/summary.py                                      # writes results/summary.json
```
