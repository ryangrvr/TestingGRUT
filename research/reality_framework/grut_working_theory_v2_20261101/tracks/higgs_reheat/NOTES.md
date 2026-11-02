# Track higgs_reheat: Higgs vacuum survival through inflation and reheating in M1r with ξ_H = 0.1

2 November 2026 · Claude, decisive pass (BRIEF_DECISIVE.md item 3; GitHub issue #13 context). This is same-session AI
work on one supplied model, not outside replication. Nothing here is a validation, an exclusion of the real world or
a new theory. ξ_H = 0.1 is a **supplied model choice**; closure does not fix it. Data enter only through the
retrospective PDG 2026 inputs (M_h = 125.13, α_s = 0.1180, M_t = 172.60 ± 0.27 GeV) reused from cf_vacuum.

Labels: [D] derived, [C] computed with our own code (`code/`, output in `results/`), [E] estimate or criterion choice,
[A] argued, [Q] quoted (abstract level unless stated), [R] recalled.

**Convention.** V ⊃ ½ ξ R h² (conformal value 1/6; ξ > 0 gives a positive mass in de Sitter). "ξ_H = 0.1" is read as
ξ(M_t) = 0.1, the literature's ξ_EW. One-loop running gives ξ(H_i) = 0.0878 (§1). The alternative reading
ξ(H_i) = 0.1 is also reported where it matters.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer on one screen

| # | Result | Label |
|---|---|---|
| 1 | **SM potential** (reused 3-loop RGEs, SMDR matching, 2-loop λ_eff of cf_vacuum; reproduced to 4 digits): h_max = 2.69×10¹¹ GeV (M_t = 172.60), 8.0×10¹¹ (−1σ), 1.06×10¹¹ (+1σ); Λ_I = 3.49×10¹¹, 1.04×10¹², 1.36×10¹¹. At H_i = 1.29–1.51×10¹¹ GeV this gives h_max/H_i = 2.09–1.78 (central), 6.2–5.3 (−1σ), 0.82–0.70 (+1σ). **ξ(M_t) = 0.1 runs to ξ(H_i) = 0.0878** (16π²β_ξ = (ξ − 1/6)(12λ + 2Y₂ − 3g′²/2 − 9g²/2)). | [C] |
| 2 | **Inflation (de Sitter, R = 12H²).** Curvature mass m² = 0.937 H_i² (one loop; tree 12ξ(H_i) = 1.054). This is below 9H²/4, so the Higgs fluctuates: ⟨h_i²⟩ = 3H⁴/(8π²m²) = 0.0406 H_i² per real component (the exact Bunch–Davies super-horizon value is 0.0322). The de Sitter barrier sits at h_b = 15.5–24.6 H_i (2.0–3.6×10¹² GeV). The Hawking–Moss exponent is **B = 1564–4298** (one loop; 1087–2096 tree RG-improved), against survival of e^{3N} patches needing B > 3N + ln(4π/9) = 180 (N = 60). Expected AdS bubbles in our past light cone: 10⁻⁶⁰¹ to 10⁻¹⁷⁸⁸. **ξ_min(M_t) = 0.044 (one loop) / 0.050 (tree) at M_t = 172.60; 0.049 / 0.055 at +1σ; 0.037 / 0.046 at −1σ** (N = 60, either H_i; N = 50 lowers these by 0.002). So ξ_H = 0.1 is sufficient during inflation by a factor 1.8–2.7. | [C] |
| 3 | **But the end-of-inflation Hubble-patch tail is large.** With M_t ≥ 172.60, 10⁵⁷–10⁷⁶ of the e^{180} Hubble patches end inflation with \|h\| above the *flat-space* barrier h_max (Fokker–Planck tail). Whether ξ = 0.1 survives is therefore decided after inflation, not during it. | [C] |
| 4 | **Background after inflation** (our rpreheat.c engine, mode 0, (3, 0.5); own scipy integrator agrees to 1e-13). R(0) = 12; R < 0 first in τ = H_i t ∈ [0.54, 0.96] with **R_min = −1.0636 H_i² at τ = 0.705** (it reproduces the M1 review's −1.06), then −0.191 at τ = 1.81, then decaying. Radiation equality at τ = 1.79. T = 4.0–4.4×10¹³ GeV already at τ = 0.05; **T_max = 1.59×10¹⁴ (H_i = 1.29e11) to 1.72×10¹⁴ GeV (1.51e11) at τ = 0.81** (g* = 106.75). This is below the 2–5×10¹⁴ quoted in the brief. | [C] |
| 5 | **Linear k-resolved Higgs modes** (Bunch–Davies start, oscillating R, Hubble damping, Hartree back-reaction, 260 modes). At ξ(H_i) = 0.0878 the tachyonic window is shallow (m² = −0.093 H_i² for ≈ 0.4 H_i⁻¹) and nothing is amplified: the long-wave variance falls from 0.0287 to 0.0041 H_i² by τ = 7, because the blue spectrum (n − 1 ≈ 0.95) re-enters and redshifts. Resonance matters only for ξ(H_i) ≳ 10, where max σ² is 0.56 H_i² at ξ = 10 and 16 H_i² at ξ = 31. | [C] |
| 6 | **Decay criterion (ours).** A region with smoothed field h̄ > h_max can roll only once its tachyonic rate m_c = √(−V″ − ξ⟨R⟩ − m_p²) exceeds H, and only if it is larger than 1/m_c. At h_max, √(−V″) = 0.075 H_i. D = log₁₀[e^{3N} k_c³ P(\|h\| > h̄)] with the mode spectrum at that time; survival needs D < 0. | [E] |
| 7 | **Gravitational-only window** (no plasma mass ever; the literature's conservative setup). N = 60, ξ_EW: **[0.074, 26] at (172.60, 1.29e11); [0.079, 25] at (172.60, 1.51e11); [0.100, 20] at (172.87, 1.29e11); [0.104, 19] at (172.87, 1.51e11); [0.043, 32] at (172.33, 1.51e11).** The upper edge is 11 (central) or 2.2–3.0 (+1σ) if instead one uses the mean-field criterion 2σ_max = h_max. The benchmark gives **D = −229 / −168 at central M_t (safe), +1.1 / +18 at +1σ (marginal / fails).** | [C][E] |
| 8 | **With the plasma mass of the actual radiation** (m_T² = κ_T T², κ_T = 0.107–0.116): m_p² ≈ 1e4 H_i² by τ = 0.05 and 1.6e5 H_i² by τ = 0.7, far above ξ\|R_min\| and \|V″\|. There is no growth for any ξ(H_i) from 0.05 to 1000 (computed). The post-inflation lower edge disappears, and **the window becomes [ξ_min from inflation = 0.044–0.049, ≳ 10³]**. The estimated true upper edge is ~10⁴, where M_P/ξ falls to T_max, or ~10⁵, where ξ\|R_min\| reaches κ_T T². At +1σ and 1.51e11 the benchmark stays safe if the plasma mass appears by τ ≈ 15 (D = −20 at τ_th = 10, +16 at 20). | [C][A] |
| 9 | **Thermal vacuum decay** (O(3) bounce, full one-loop J_B/J_F, couplings at πT). The thermal barrier sits at h = 2.8–3.5 T. **S3/T = 581–896 for T = 10¹³–3×10¹⁴ GeV** (600 at T_max, M_t = 172.87). Expected thermal bubbles in our past light cone: 10⁻¹⁶⁶ to 10⁻²⁸⁰. Negligible. | [C][E] |
| 10 | **Verdict.** ξ_H = 0.1 is inside the allowed window. It is safe during inflation at every M_t and H_i considered (margin ×1.8–2.7), and safe through the oscillating-R era. Once the reheating radiation gives the Higgs its plasma mass, survival through reheating is unconditional. Without a plasma mass (sterile decay products, or thermalisation later than τ ≈ 15), it is safe at central M_t (lower edge 0.074–0.079) but on the edge at M_t + 1σ (needs 0.100–0.104). The upper edge (≥ 2–26 even without plasma) is irrelevant for 0.1. | [C][A] |

---------------------------------------------------------------------------------------------------------------------

## 1. SM effective potential with ξR — `code/smpot.py` [C]

**Reused, copied unchanged.**
- `sm_rge.py`: track vacua's 3-loop RGEs, the Higgs anomalous dimension and the Landau-gauge 2-loop λ_eff
  (Buttazzo et al. 1307.3536 App. A/B).
- `c1_fit_grid.json`: cf_vacuum's 11 SMDR fits at PDG 2026.
- `build_matching()`: copied from cf_vacuum's c2_sm_relation.py.

**Validation.** At M_t = 172.60 the code gives Λ_I(2-loop) = 3.4883×10¹¹ and h_max = 2.6915×10¹¹ GeV. cf_vacuum's
sm_relation.json has 3.4883×10¹¹ and 2.6925×10¹¹.

| M_t | Λ_I (GeV) | h_max (GeV) | h_max/H_i (1.29 / 1.51e11) |
|---|---|---|---|
| 172.33 | 1.04×10¹² | 8.04×10¹¹ | 6.23 / 5.32 |
| 172.60 | 3.49×10¹¹ | 2.69×10¹¹ | 2.09 / 1.78 |
| 172.87 | 1.36×10¹¹ | 1.06×10¹¹ | 0.82 / 0.70 |

**New pieces.**
- **One-loop running of ξ.** 16π²β_ξ = (ξ − 1/6)(12λ + 2Y₂ − 3g′²/2 − 9g²/2), with Y₂ = 3y_t² + 3y_b² + y_τ²
  [Q: 1804.02020 eq. (4.21), read in the PDF]. It is solved exactly along the 3-loop trajectory. ξ(M_t) = 0.1 gives
  ξ = 0.0877 (10¹⁰ GeV), 0.0878 (H_i), 0.0879 (√R = 4.5×10¹¹) and 0.0883 (10¹³). The map is
  (ξ_H − 1/6) = 1.1837 (ξ_EW − 1/6).
- **De Sitter one-loop potential.** We use Markkanen–Nurmi–Rajantie–Stopyra 1804.02020, eqs. (5.2)–(5.3) and
  Tables 1–2, read in the PDF; Landau gauge.
  - Field-dependent masses: W, Z, t, b, τ, h and the Goldstones, each with its H² shift and its n′H⁴ log term.
  - Couplings run at 3 loops, following Mantziris–Markkanen–Rajantie 2011.03763 [Q abstract].
  - Field-independent pieces cancel because ΔV = V(h) − V(0) is taken at the same μ (our prescription).
  - Scale: μ² = h² + R; the alternative μ² = h² changes B by < 0.2%.
  - Flat-space check: for H → 0 the code reproduces cf_vacuum's one-loop λ_eff to 1e-9.

## 2. During inflation — `code/inflation.py` → `results/inflation.json` [C]

H = H_i is constant over the last N e-folds; that is an assumption [A]. For H_i ≈ 1.3–1.5×10¹¹ GeV and near-instant
reheating, our present horizon left at N ≈ 61 [E].

| M_t, H_i | m²/H² (1-loop) | h_b/H | B (1-loop / tree) | ξ_EW,min N=60 (1-loop / tree) | ξ(H_i),min | Patches above flat h_max at end (log₁₀) |
|---|---|---|---|---|---|---|
| 172.33, 1.29e11 | 0.938 | 24.6 | 4298 / 2096 | 0.037 / 0.045 | 0.013 | −127 |
| 172.33, 1.51e11 | 0.938 | 23.7 | 3965 / 2016 | 0.038 / 0.046 | 0.014 | −71 |
| 172.60, 1.29e11 | 0.937 | 18.9 | 2402 / 1461 | **0.044 / 0.050** | 0.021 | +57 |
| 172.60, 1.51e11 | 0.937 | 18.4 | 2278 / 1418 | 0.044 / 0.051 | 0.022 | +63 |
| 172.87, 1.29e11 | 0.935 | 15.8 | 1627 / 1113 | 0.049 / 0.054 | 0.027 | +76 |
| 172.87, 1.51e11 | 0.936 | 15.5 | 1564 / 1087 | 0.049 / 0.055 | 0.028 | +76 |

**Variances.** Per real component:
- Stochastic, light-field: 3H⁴/(8π²m²) = 0.0406 H².
- Exact Bunch–Davies super-horizon (k < aH; Hankel modes, N = 50 or 60 give the same): 0.0322 H². It is 0.011 for
  k < 0.3 aH and 0.0048 for k < 0.1 aH.
- The spectrum is blue: P(k) ∝ k^{3−2ν} with 3 − 2ν ≈ 0.95. The Hubble-scale field is therefore carried by the last
  e-folds.
- For ξ above 3/16 the field is heavy and the super-horizon variance falls only as ≈ 1/m²: 0.0127 H² at m² = 2,
  0.0058 H² at m² = 4.

**Literature consistency [Q abstracts].**
- Herranen et al. 1407.3141: ξ_EW ≳ 0.06 stabilises; ≲ 0.02 destabilises.
- Mantziris et al. 2011.03763: ξ ≳ 0.051–0.066, for inflation models with higher H.
- Figueroa et al. 1709.00398: ξ ≳ 0.06.
- Our lower bound of 0.044–0.049 at H_i ≈ 1.3–1.5×10¹¹ is smaller, as expected. At the de Sitter barrier scale
  (~2×10¹² GeV), λ_eff is only about −0.002, against about −0.0035 to −0.006 at the 10¹³–10¹⁴ GeV scales of those
  models (2-loop flat λ_eff at M_t = 172.60).

**The answer to "is ξ_H = 0.1 sufficient at H_i ≈ 1.3–1.5×10¹¹ during inflation?" is yes.** B exceeds the
requirement by ×8.7–24, and ξ exceeds ξ_min by ×1.8–2.7 (×2.0–2.3 at central M_t). Under the reading ξ(H_i) = 0.1 the
margin is larger still (ξ(H_i),min = 0.013–0.028).

**Why this is not yet the whole answer.** For M_t ≥ 172.60 the flat-space barrier h_max is only 0.7–2.1 H_i. Many
patches therefore end inflation beyond it, while still well inside the de Sitter barrier h_b ≈ 15–25 H_i. Their fate
is settled after inflation (§3).

## 3. After inflation

### 3.1 Background — `code/background.py` → `results/background.json` [C]

We use our rpreheat.c, copied unchanged, in mode 0 with b = 3 and Γ = 0.5: an oscillating quadratic inflaton released
at rest with H = H_i, decaying into radiation.
- The hidden fields are put at a static point (X = Y = 0, Σ at its vev). Their contribution to R is η²V ~ 10⁻¹⁰ H_i²
  [D], so R(t) is that of the M1r release.
- An independent scipy integration of the same inflaton + radiation system agrees to 7×10⁻¹⁴ (R), 1×10⁻¹³ (H) and
  7×10⁻¹³ (ρ_r).

| Quantity | Value |
|---|---|
| R(0) | 12.000 H_i² (continuous with de Sitter) |
| R < 0 windows (τ) | [0.54, 0.96], [1.64, 2.04], [2.70, 3.10], [3.76, 4.16], ... |
| R minima | −1.0636 (τ = 0.705, H = 0.530), −0.1906 (1.805), −0.0568 (2.875), −0.0208 (3.94) H_i² |
| R maxima after τ = 0 | 0.831 (1.235), 0.204 (2.32), 0.068 (3.39) H_i² |
| Radiation equality | τ = 1.79 (a = 2.63) |
| T(τ = 0.05), T(0.3), T(1), T(10) | 4.0×10¹³, 1.21×10¹⁴, 1.55×10¹⁴, 6.8×10¹³ GeV (H_i = 1.29e11) |
| T_max | 1.59×10¹⁴ (H_i = 1.29e11), 1.72×10¹⁴ GeV (1.51e11), at τ = 0.81 |

### 3.2 Higgs modes — `code/modes.py` → `results/modes.json` [C]

**Equations.** For each of the 4 real components and comoving k (a(0) = 1):

  h_k″ + 3H h_k′ + [k²/a² + ξR + m_p² + 6λ_eff σ²] h_k = 0

- The background is integrated in the same system.
- The initial modes are exact Bunch–Davies de Sitter modes at τ = 0 with m² = 12ξ H². The Wronskian error is ≤ 7×10⁻¹⁵.
- σ² = σ_L² (k < 1, unsubtracted) + σ_P² (k ≥ 1, zeroth-order adiabatic subtraction).

**The benchmark ξ(H_i) = 0.0878** (ξ_EW = 0.1; M_t = 172.6, H_i = 1.29e11), from `series_*` in modes.json:

| τ | 0 | 0.70 | 2.0 | 7.0 | 20 | 40 |
|---|---|---|---|---|---|---|
| ξR + … (H_i²) | 1.055 | −0.093 | −0.004 | 0 | 0 | 0 |
| σ_L² (H_i²) | 0.0287 | 0.0149 | 0.0090 | 0.0041 | 0.0024 | 0.0016 |

- The tachyonic window is too short and too shallow: |m²|^{1/2} Δτ ≈ 0.12. Occupation numbers stay below 1 (max
  n_k = 0.91 at ξ = 0.088).
- The variance decreases monotonically, because modes near the end-of-inflation horizon re-enter and redshift.
- At ξ = 10: σ_P² rises to 0.55 H_i² at τ = 0.85 (2σ = 0.71 h_max) and n_k ≈ 3.4.
- At ξ = 31: σ² reaches 16 H_i² (2σ = 3.8 h_max). There is still no Hartree runaway, because the produced quanta have
  k/a ~ H_i ≫ √|V″| = 0.075 H_i.

**The decay criterion D [E]** (defined in the modes.py docstring).
- A region of comoving size 1/k_c and mean field h̄ ≥ h_max rolls only once m_c² = −V″(h̄) − ξ⟨R⟩ − m_p² exceeds H².
  ⟨R⟩ is the one-period running mean.
- At that first time τ_c: k_c = a(τ_c) m_c. σ_c² is the variance of the modes with k < k_c at τ_c; all of them are
  super-horizon then.
- D = log₁₀[e^{3N} k_c³ e^{−x}(1 + x)], with x = h̄²/(2σ_c²); the maximum is taken over h̄ ∈ [h_max, 20 h_max].
- Survival needs D < 0.
- This replaces the cruder "Hubble patch above h_max at the end of inflation" count (row 3). That count ignores two
  facts: those Hubble-scale fluctuations are gradient-dominated (k/a ≈ H_i ≫ m_c), and they redshift away. The
  replacement is an argued criterion, not a lattice result [E/A].
- Sensitivity: replacing m_c > H by m_c > H/2 or m_c > 2H moves the lower edge by < 3% and the upper edge by about
  ±25%.

**Numerical controls** (`code/controls.py` → `results/controls.json`, at both edges and the benchmark).
- At the lower edge, D moves by ≤ 0.8 under these changes: n = 400 modes, K_max ×1.5, rtol 1e-10, double sampling,
  Hartree off. This shifts the edge by ≪ 1% in ξ.
- At the upper edge, D moves by −4.2 → +1.4, a shift of ≈ 1–2% in ξ. Hartree off moves it to −28.8.

### 3.3 The window — `code/edges.py` → `results/edges.json` [C][E]

**Gravitational only** (no plasma mass ever; bisection on D_N = 0):

| M_t, H_i | Benchmark D₆₀ (τ_c) | ξ_EW lower (N=60 / 50) | ξ_EW upper, D (N=60) | ξ_EW upper, mean field 2σ = h_max |
|---|---|---|---|---|
| 172.33, 1.51e11 | −1087 (2.9) | 0.043 / < 0.043 (no edge above ξ_H = 0.02) | 32 | 20 |
| 172.60, 1.29e11 | **−229** (7.0) | **0.074 / 0.071** | 26 | 11 |
| 172.60, 1.51e11 | −168 (8.1) | 0.079 / 0.075 | 25 | 9.7 |
| 172.87, 1.29e11 | +1.1 (16.0) [−6.7 / −5.7 at m_c = H/2, 2H] | 0.100 / 0.098 | 20 | 3.0 |
| 172.87, 1.51e11 | **+18** (18.7) | 0.104 / 0.101 | 19 | 2.2 |

Upper-edge literature for other inflation models [Q abstracts]:
- Herranen et al. 1506.04065: ξ ≳ 1 is dangerous.
- Figueroa–Rajantie–Torrenti 1709.00398 (lattice): instability for ξ ≳ 4–5 at M_t = 173.3, and ≳ 11–12 at 172.1.
- Ema–Mukaida–Nakayama 1602.00483: upper bounds.
- Kohri–Matsui 1607.08133: H < Λ_I is needed. Ours: H_i/Λ_I = 0.37–0.43 at central M_t, but 0.95–1.11 at +1σ,
  exactly where we find marginality.

Our mean-field edges are of the same order as these. The higher D-edge reflects our much lower H_i/h_max.

**With the radiation's plasma mass.** m_p² = κ_T T(τ)², with T from the engine's ρ_r(τ); results in modes.json under
"plasma …".
- Prompt (τ_th ≤ 1):
  - σ_L² never grows and D = −∞, for ξ(H_i) = 0.05 … 100 at both (172.6, 1.29e11) and (172.87, 1.51e11).
  - For ξ = 100, 300 and 1000 to τ = 2.5, the total mass² stays ≥ 12ξ H_i² (its τ = 0 value) throughout.
- Delayed onset, for the benchmark at the worst corner (172.87, 1.51e11):
  - D₆₀ = −935 (τ_th = 0.5), −775 (2), −120 (5), −20 (10), +16 (20).
  - So the plasma mass must appear by **τ ≈ 15 H_i⁻¹ (≈ 7×10⁻³⁵ s)**.
- Large ξ with delayed plasma: ξ(H_i) = 100 survives τ_th = 2 at central M_t but not τ_th ≥ 5.

**Window statement.**
- Lower edge: with a prompt plasma mass, it is the inflationary one, ξ_EW = 0.044–0.049 (one loop) or 0.050–0.055
  (tree), for M_t = 172.60–172.87.
- Upper edge: ≥ 10³ with prompt plasma (computed). The estimate is ~10⁴, set by EFT control (M_P/ξ = T_max at
  ξ ≈ 1.5×10⁴) or by ξ|R_min| ≈ κ_T T² at ξ ~ 10⁵ [E].
- Without a plasma mass the window is the gravitational-only table above.

### 3.4 Thermal effects — `code/thermal.py` → `results/thermal.json` [C][E]

**Thermal mass.** κ_T = m_T²/T² = 3g²/16 + g′²/16 + y_t²/4 + λ/2 = 0.110 at T = 10¹⁴ GeV (couplings at πT).

**One-loop finite-T potential.** We use full J_B/J_F for W, Z and t; Higgs/Goldstone loops and daisy resummation are
omitted [E]; λ_eff is taken at max(h, T), which is conservative.
- The thermal barrier sits at h = 2.85–3.51 T.
- The bounce solver reproduces S3 = 18.897 m/|λ| for the pure quartic, and the Derrick identity to 1e-8.

| M_t | S3/T at 10¹³ | 6×10¹³ | 1.6×10¹⁴ | 3×10¹⁴ | log₁₀ thermal bubbles in past light cone (N = 60; H_i = 1.29 / 1.51e11) |
|---|---|---|---|---|---|
| 172.33 | — | — | 867 | 829 | −280 / −280 |
| 172.60 | 896 | 759 | 708 | 682 | −214 / −213 |
| 172.87 | 735 | 638 | 601 | 581 | −167 / −166 |

The danger threshold is S3/T ≈ 215. Thermal nucleation is therefore negligible even if T_RH were 3×10¹⁴ or higher;
S3/T falls only slowly with T. The plasma's real role is the stabilising mass of §3.3.

## 4. Verdict with numbers

1. **During inflation [C]: ξ_H = 0.1 is sufficient** at H_i = 1.29–1.51×10¹¹ GeV for M_t = 172.33–172.87.
   - B_HM = 1564–4298 against 180 needed.
   - ξ_EW,min = 0.044 (central M_t; tree 0.050) and 0.049 (+1σ; tree 0.055).
   - The Higgs is light (m² = 0.94 H_i²), with σ = 0.18–0.20 H_i per component.
2. **Through the oscillating-R era [C]:** nothing grows at ξ = 0.1.
   - R_min = −1.064 H_i² gives ξR = −0.093 H_i² for ~0.4 H_i⁻¹.
   - Curvature-driven resonance needs ξ(H_i) ≳ 10–30 before it threatens: mean-field 2σ = h_max at ξ_EW ≈ 11
     (central) or 2.2–3 (+1σ); D-edge at 19–26.
3. **Through reheating [C][A]:** the decisive input is whether the inflaton's decay products give the Higgs a plasma
   mass.
   - In the engine, radiation with T ≈ 4×10¹³ GeV exists by τ = 0.05, and m_T² ≈ 10⁴–10⁵ H_i² ≫ everything.
   - If those products are Higgs-coupled SM states, survival is unconditional for ξ_EW ∈ [0.044–0.049, ≳ 10³].
   - If they are not (or thermalise later than τ ≈ 15 H_i⁻¹), the gravitational-only window applies:
     - ξ = 0.1 is safe at M_t ≤ 172.60 (lower edge 0.074–0.079; D = −229/−168);
     - it is marginal at M_t = 172.87 with H_i = 1.29e11 (needs 0.100);
     - it fails at M_t = 172.87 with H_i = 1.51e11 (needs 0.104), the Q-matter calibration.
   - ξ_EW ≥ 0.11 would remove the marginal corner in every case we computed. That is an observation, not a proposed
     change.
4. **Thermal decay [C]:** S3/T ≥ 581 against ~215 needed. Negligible.
5. **WORKING_THEORY.md's "ξ_H ≈ 0.05–1" can be sharpened for this background:**
   - **[0.044–0.049, ≳ 10³]** with plasma;
   - **[0.074–0.104, 19–26]** without, or [.., 2.2–11] on the mean-field criterion.
   - ξ_H = 0.1 is inside in the first case, and inside except at the +1σ top mass in the second.

## 5. What is not established

- **The criterion D.** It is a coarse-grained classical argument: Gaussian statistics, a sharp onset when m_c = H,
  and regions counted as k_c³. It is not a lattice simulation. It accounts for the end-of-inflation Hubble-patch tail
  dissolving, which we argue rather than simulate.
  - The cruder Hubble-patch count says "fails" for M_t ≥ 172.60 without plasma.
  - The literature's mean-field criterion gives a lower upper edge.
  - Both readings are reported; neither changes the verdict for ξ = 0.1 with plasma.
- **Inflation side.**
  - H is constant through the last 60 e-folds.
  - The Hawking–Moss prefactor is set to H⁴.
  - The one-loop de Sitter potential is valid only in the UV-expansion sense of 1804.02020.
  - The tree-versus-loop spread (ξ_min 0.044 vs 0.050) is the stated scheme error.
- **Plasma side.**
  - Instant thermalisation is assumed, in the HTL thermal-mass form.
  - The decay channel of the inflaton is not specified by M1r. WORKING_THEORY requires a fermionic channel. Sterile
    intermediates with lifetimes > 15 H_i⁻¹ would put the model in the gravitational-only case.
  - Non-thermal hard decay products would give a larger induced mass (∝ ρ_r/E²), not a smaller one [E].
- **Not included.**
  - The hidden sector's own effect on the Higgs (zero portals): none at tree level [D].
  - Higgs–inflaton direct couplings: none in M1r by choice.
  - Gravitational-wave or isocurvature consequences of the Higgs fluctuations.
- **Our T_max of 1.6–1.7×10¹⁴ GeV** comes from the engine's rule ρ_r → T with g* = 106.75. The brief's 2–5×10¹⁴ is a
  different estimate. Thermal stability holds across both.

## 6. Files and reproduction

| File | Content |
|---|---|
| `code/smpot.py` | matching, 3-loop running, ξ running, flat λ_eff, de Sitter one-loop ΔV, barrier finder |
| `code/sm_rge.py`, `code/c1_fit_grid.json` | reused unchanged from cf_vacuum (track vacua's RGEs; SMDR fit grid) |
| `code/rpreheat.c` | our reheating engine, copied unchanged from rp_reheating |
| `code/inflation.py` → `results/inflation.json` | item 2 (B_HM, ξ_min, variances, tails) |
| `code/background.py` → `results/background.json` | R(τ), H, a, ρ_r, T; engine vs scipy |
| `code/modes.py` → `results/modes.json` | k-resolved modes; gravitational-only scan; plasma scan; benchmark series |
| `code/edges.py` → `results/edges.json` | window edges by bisection, benchmark D with m_c-criterion variants |
| `code/thermal.py` → `results/thermal.json` | κ_T, finite-T potential, O(3) bounces, thermal bubble counts |
| `code/controls.py` → `results/controls.json` | numerical controls of modes.py |
| `code/summary.py` → `results/summary.json` | headline numbers |

**Reproduction, in `code/`, one thread:**

```
python3 -B smpot.py                     # validation, ~20 s
python3 -B inflation.py                 # ~40 s
python3 -B background.py                # builds rpreheat in the scratchpad, ~3 s
python3 -B modes.py scan                # ~4 min
python3 -B modes.py plasma              # ~3 min
python3 -B modes.py series
python3 -B modes.py prompt_large
for a in "172.6 1.29e11" "172.6 1.51e11" "172.87 1.29e11" "172.87 1.51e11" "172.33 1.51e11"; do
    python3 -B edges.py $a              # ~2 min each
done
python3 -B thermal.py 172.87            # ~4 min
python3 -B thermal.py 172.6             # ~4 min
python3 -B thermal.py 172.33 1.6e14 3e14
python3 -B controls.py                  # ~1 min
python3 -B summary.py
```

The scratch build lives at `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/higgs_reheat/build/`. The text of
1804.02020 (Tables 1–2, eq. 4.21) was read from the arXiv PDF, saved in the scratchpad (`dl/mnrs/`), not in the
repository.

**Credit.**
- SM running and matching: Buttazzo et al. 1307.3536; SMDR, Martin–Robertson 1907.02500; cf_vacuum and vacua tracks.
- Curved-space potential: Markkanen–Nurmi–Rajantie–Stopyra 1804.02020.
- Inflationary stability: Herranen et al. 1407.3141; Kearney–Yoo–Zurek 1503.05193; Espinosa et al. 1505.04825;
  Mantziris–Markkanen–Rajantie 2011.03763; Markkanen–Rajantie–Stopyra 1809.06923 (review).
- Post-inflation: Herranen et al. 1506.04065; Ema–Mukaida–Nakayama 1602.00483; Kohri–Matsui 1607.08133;
  Figueroa–Rajantie–Torrenti 1709.00398.
- Hawking–Moss / stochastic: Starobinsky–Yokoyama [R].

Nothing here is GRUT-native: the calculation applies published physics to M1r's supplied background.
