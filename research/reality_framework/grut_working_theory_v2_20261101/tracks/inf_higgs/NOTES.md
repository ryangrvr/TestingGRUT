# Track inf_higgs: Higgs survival on the inflaton benchmarks with ξ_H = 0.1

3 November 2026 · Claude, inflaton-sector pass (BRIEF_INFLATON.md, item "Higgs"; GitHub issue #13 context). This is
same-session AI work on one supplied model, not outside replication. Nothing here is a validation, an exclusion of the
real world or a new theory. **Every choice below is a model choice, not a consequence of closure:** ξ_H = 0.1, the
inflaton benchmarks, the decay history, the reheating channel and its couplings (y, M_N), the R-even portal λ_IH, the
renormalisation of the induced A-term, and the survival criterion. Git was read-only and nothing was committed. Every
run was in the foreground with one worker process.

Labels: [D] derived, [C] computed with our own code (`code/`, output in `results/`), [E] estimate or criterion choice,
[A] argued, [Q] quoted from a reviewed track or an abstract, [R] recalled.

**Inputs.**
- M_t = 172.57 ± 0.29, M_h = 125.20 ± 0.11 and α_s = 0.1180 (PDG 2024 [R]). The M_h = 125.13 rows connect to higgs_reheat.
- ξ_H = 0.1 is read as ξ(M_t) = ξ_EW, as in higgs_reheat. It runs to ξ(H_i) = 0.0877.
- SM running and λ_eff: higgs_reheat/code/smpot.py, imported read-only.
- Higgs modes: the skeptic's normalisation-fixed modes_canon.py, generalised. The barrier and V″ are in the canonical
  field φ = e^{Γ(H_i)} h₀.
- Survival criterion: higgs_reheat's scale-resolved D, with survival for D < 0. It is an argued criterion [E], not a
  lattice result.

The benchmarks are inf_potential's HT6t, T1 and E1 at H_i = 1.29×10¹¹ GeV (and 1.50×10¹¹). H_i is H at the slice S1,
one e-fold before ε_H = 1. τ = H_i t is counted from S1, with a(S1) = 1. **ε_H = 1 is at τ = 1.00**, so subtract 1.0 to
measure from the end of inflation.

---------------------------------------------------------------------------------------------------------------------

## 0. Answer on one screen

| # | Result | Label |
|---|---|---|
| 1 | **Regression.** On the M1r surrogate the generalised code gives D₆₀ = +12.51 at (172.87, 1.29×10¹¹, M_h = 125.13), against the skeptic's +12.47, and −181.9 against −182.2 at 172.60. τ_c = 15.03 (skeptic 15.07); deadline τ = 9.95 (skeptic ≈ 10) | C |
| 2 | **Inflation with the benchmark's H(N).** R/H² falls from 12 to 6 at ε_H = 1, and H by 1.3% (T1, E1) or 4.3% (HT6t). At ε_H = 1 the Hawking–Moss exponent falls to about 0.27 of its constant-H value (tree: 430 against 1580 at 172.57; 310 against 1160 at 172.86). The bubble count of Mantziris–Markkanen–Rajantie 2011.03763 is then dominated by the end of inflation. **ξ_EW,min = 0.070 / 0.056 (tree / one loop) at 172.57 and 0.079 / 0.067 at 172.86**, against 0.048 / 0.041 and 0.053 / 0.047 for constant H. The three potentials agree to 2×10⁻⁴. ξ_H = 0.1 is still safe during inflation, with a margin of ×1.27–1.5 at +1σ. The non-equilibrium Gaussian count of patches beyond the instantaneous de Sitter barrier is 10⁻⁷⁰⁰–10⁻²⁰⁰⁰ | C [E criterion] |
| 3 | **Post-inflation R(t) of the real oscillation.** Every R < 0 window reaches R_min = −6.0 H² (M1r surrogate: −3.8 H²). The windows last only 0.028 / 0.011 / 0.005 H⁻¹ (HT6t / T1 / E1; surrogate 0.24), with a tachyonic exponent of 0.015 / 0.006 / 0.003 per window (surrogate 0.10). The anharmonic shape changes the fraction of time with R < 0 (0.37 HT6t, 0.35 T1/E1, 0.39 harmonic). At ω/H = 30–255 none of this compounds (Floquet q ≈ 9ξH²/4ω² ~ 10⁻⁴–10⁻⁶), so the Higgs sees ⟨ξR⟩ = 3ξH² > 0 in the matter era. **Shape and fragmentation move D by ≤ 2 between HT6t, T1 and E1, and by ≤ 1.2 for the w = 0 fragmentation switch** | C |
| 4 | **D at ξ_EW = 0.1, gravitational only (no Higgs-coupled plasma), H_i = 1.29×10¹¹, N = 60.** Slow decay, i.e. a matter era (inf_reheat's Γ_I/m_I ≲ 10⁻⁶): **D = −2023 / −365 / −23.8 at M_t = 172.28 / 172.57 / 172.86**, the same on all three benchmarks to ±2. At the benchmarks' own N_obs = 54.7 these are −2030 / −372 / −30.6. With Γ_I = 0.5 H_i: −1970 / −312 / −6.2. The M1r surrogate gives −1617 / −310 / −5.1. **ξ_H = 0.1 is safe at M_t + 1σ on every benchmark, but only by tens in D.** The criterion's own band (k_c ×½, ×2; onset at H/2, 2H) at 172.86 is [−71, −12] for the matter era and [−56, +14] for Γ = 0.5 | C [E criterion] |
| 5 | **Critical top mass** (D = 0 at ξ_EW = 0.1, gravitational only). **Matter era: 172.92 GeV (N = 60) / 172.94 (N_obs), i.e. +1.2σ**, with criterion band 172.89–173.00; M_h ± 1σ shifts it by ±0.06 (same sign). Γ_I = 0.5: 172.875 / 172.884 (+1.05σ), equal to the M1r surrogate (172.874). At the Q-matter calibration H_i = 1.50×10¹¹: 172.87–172.89 (matter era) and 172.83 for HT6t at Γ = 0.5, which fails at +1σ with D = +12 | C [E] |
| 6 | **Plasma deadline** (above the critical M_t; τ from S1, sharp onset). Matter era: τ = 16 (172.95), 11 (173.05), 8.6 (173.15), 5.7 (173.45). Γ = 0.5: 12 (172.90), 9.7 (172.95), 6.6 (173.15), 4.9 (173.45). The surrogate's 9.95 is at 172.87/125.13 | C |
| 7 | **The benchmark's actual plasma timing** (inf_reheat's S1 channel −(y/2)φ_I N̄N with seesaw N → LH, transcribed onto the real oscillation). Instant preheating at each zero of m_N (HT6t: yΦ > m_I/2 at y ≳ 3×10⁻⁵), or an open perturbative channel (T1/E1: m_I = 128–256 H_i against Φ ≈ 0.03 M_P, open for y ≲ 10⁻⁴), gives a hard-gas induced Higgs mass m_ind² ≥ 10⁻³ H_i² **by τ ≈ 0.95–1.7** for y ≥ 10⁻⁵ (3.4 for HT6t at y = 10⁻⁵, M_N = m_I/16). That is 0–2.4 H_i⁻¹ after the end and **before every deadline** of row 6. With m_ind²(τ) put into the criterion, **D < 0 at every M_t ≤ 174.6 GeV for y ≥ 3×10⁻⁵** on all three benchmarks; at y = 10⁻⁵ HT6t crosses at 173.1–173.4 and T1 (M_N = m_I/16) at 173.4. A fully thermal bath (inf_reheat's pessimistic LPM time) does not form by τ = 60 | E (rates C on the background) |
| 8 | **Catch: the induced R-odd A φ_I\|H\|² of S1.** During inflation it shifts the Higgs mass by A φ_I on one sign branch. Survival needs **\|A\| ≤ 12 TeV (HT6t), 34 TeV (T1), 35 TeV (E1) at +1σ** (tree; one loop +10–35%), i.e. \|Aφ_I\| ≲ 0.6–0.9 H_i² on the plateau. The natural size A_ind ≈ y y_ν² M_N/8π² × [1, ln(M_P²/M_N²)] **exceeds this at y_nat with M_N = m_I/4** (11–320 TeV HT6t; 150–4100 T1; 630–16000 E1). It passes for small y and/or M_N = m_I/16 (e.g. HT6t y ≤ 10⁻⁴; T1 y ≤ 3×10⁻⁵; E1 y = 10⁻⁵, all at m_I/16, including the logarithm). The plasma still comes by τ ≈ 1.0–1.4 there (3.4 at y = 10⁻⁵ for HT6t), and D stays < 0 to 173.1–174.6 GeV | C/E |
| 9 | **R-even portal (λ_IH/2)φ_I²\|H\|² (m_H² = λφ_I²/2).** λ_min = (20/3)(m_I/M_P)² = 3.4×10⁻¹¹ / 3.1×10⁻¹⁰ / 1.2×10⁻⁹ (HT6t / T1 / E1) keeps ⟨m_H²⟩ ≈ 10 H² in the matter era (checked) and makes the Higgs heavy during inflation (390–1750 H_i² at S1). q(λ_min) = 5×10⁻³–6×10⁻⁴ ≪ 1, and λ ≤ 7×10⁻⁹ / 3×10⁻⁷ / 2×10⁻⁶ for q ≤ 1. With it, **D ≤ −3300 at every M_t ≤ 174.6 (≤ −820 at 0.1 λ_min)**. Its own φφ → HH* annihilation gives m_ind² ≈ 0.007–0.012 H_i² in the first Hubble time | C/E |
| 10 | **ξ window at central M_t (+1σ).** **With a prompt Higgs-coupled plasma: [0.056–0.070 (0.067–0.079), ~25 (9)] in ξ_EW**: the lower edge is the inflationary one with H(N) (row 2); the upper edge is the mean-field 2σ_prod = h_max on HT6t. **Gravitational only: [0.067 (0.095), same]**: the lower edge is post-inflationary, 0.114 at +2σ. Upper edges: HT6t ξ(H_i) ≈ 29 (11), T1 ≈ 70 (37), M1r surrogate ≈ 13 (5). The scale-resolved D stays < 0 to ξ(H_i) = 300, and ξ(H_i) = 1000 runs away (Hartree) on every benchmark. **ξ_H = 0.1 lies inside every version, by 5% at +1σ without plasma** | C [E] |

**Bottom line.**
- **During inflation** ξ_H = 0.1 is safe on every benchmark, but less comfortably than with constant H. The end-of-inflation
  drop of R/H² raises ξ_min to 0.067–0.079 at +1σ.
- **After inflation, gravitational only,** the real benchmarks are slightly *safer* than M1r's surrogate in a matter era,
  because ⟨R⟩ = 3H² > 0. They pass M_t + 1σ (D = −22 to −31) and fail above **M_t,crit = 172.92–172.94 GeV (+1.2σ)**.
  With Γ_I = 0.5 they behave like the surrogate (172.875).
- **With the benchmark's actual plasma timing,** a Higgs-coupled hard gas appears within ~1 H_i⁻¹ of the end of inflation,
  before every deadline. ξ_H = 0.1 is then safe to M_t ≥ 174.6 GeV for y ≥ 3×10⁻⁵. This rests on two things:
  - the hard-gas induced mass [E]; and
  - the S1 induced A-term being kept below 12–35 TeV, by small y·M_N², a tuned counterterm, or λ_IH ≥ λ_min.
- **The R-even portal at λ_min removes the danger outright.**

---------------------------------------------------------------------------------------------------------------------

## 1. Code and method

| file | content |
|---|---|
| `code/ih_common.py` | Benchmark loader (bench_*.json of inf_potential; potentials from inf_common.py, imported read-only). Background in H_i units from N_b0 e-folds before ε_H = 1 through the oscillation, with an optional w = 0 fragmentation switch. ⟨R⟩ is the running mean over the local inflaton period. R < 0 window statistics. Cached SM tables (smpot.py) |
| `code/ih_modes.py` | Higgs modes. Bunch–Davies de Sitter modes at τ₀ with the local H(τ₀), a(τ₀) and ξR(τ₀), then the benchmark's own background. Canonical normalisation. Hartree term. Extra masses (λ_IH, A). Occupation numbers. Generalised danger function D, returning the full deadline curve D(τ_th) |
| `code/ih_run.py` | Driver: D at ξ_EW = 0.1 for every history, the M_t grid and the critical M_t (N = 60 and N_obs), criterion and M_h variants, deadlines (`bench_D.json`); lower-edge ξ scans (`xiscan_*.json`); large-ξ scans (`upper_*.json`) |
| `code/ih_inflation.py` | Hawking–Moss with H(N), R(N) (tree RG-improved, and the one-loop de Sitter potential at the instantaneous curvature); ξ_min; the A-term bound; the Gaussian tail from the mode runs (`inflation_*.json`) |
| `code/ih_plasma.py` | inf_reheat's S1 plasma model on the actual oscillation, the λ_IH window and annihilation (`plasma.json`; full curves in the scratchpad) |
| `code/ih_final.py` | D with m_ind²(τ) in the criterion (`final_plasmaD.json`); portal mode runs (`final_portal.json`); numerical controls (`final_controls.json`) |
| `code/ih_summary.py` | → `results/summary.json` |

**Start in the de Sitter phase.** Modes start at N_b0 = 5 e-folds before ε_H = 1, where H = 1.0003 H_i (HT6t). They use
exact Hankel solutions with the local H, a and m² = ξR. On these plateaus H is constant to 3×10⁻⁴ before that point, so
the Bunch–Davies state is accurate. From there to τ = 80 the modes see the benchmark's real H(N) and R(N), the
non-adiabatic end of inflation and the anharmonic oscillation.

**Two convention changes against modes_canon.**
- The Hartree variance subtracts the zeroth-order adiabatic vacuum for every mode, not only for k ≥ 1. Starting deep in
  inflation, the k < aH_end modes are still sub-horizon at τ₀; this is negligible once they are super-horizon.
- D counts only regions whose roll begins after ε_H = 1. The inflationary phase is treated by Hawking–Moss (§2), and
  regions rolling during inflation lie beyond the de Sitter barrier at h ≳ 15 H_i, with x > 10³.

The M1r control keeps modes_canon's conventions exactly (row 1).

**Histories (model choices).**
- **MD:** Γ_I = 0 over τ ≤ 80. inf_reheat's realistic channels give Γ_I/m_I ≲ 10⁻⁶; S1 at the y used in §4 gives
  Γ_I = 10⁻¹⁰–10⁻⁷ H_i.
- **MDf:** MD with the w = 0 fragmentation switch, at τ_end + 0.43 (T1) and τ_end + 0.35 (E1).
- **G05:** Γ_I = 0.5 H_i, inf_potential's convention, which its skeptic found kinematically blocked.
- **N_obs** (e-folds from today's horizon exit to S1) = N_* + 5.40 − 1. For MD it is taken at T_RH = 10¹¹ GeV: 54.7–54.8.
  For G05 it is 57.1–57.3. D(N_obs) = D₆₀ − 1.303 (60 − N_obs).

## 2. During inflation with H(N) [C] (`inflation_HN.json`, `inflation_gaussian_tail.json`, `inflation_Aterm.json`)

- **Criterion.** ⟨N⟩ = (4π/3) ∫₀^{N₀} dN_b e^{3(N₀ − N_b)} (H/H₀)³ e^{−B(N_b)}, with B = 8π²ΔV_b/3H⁴. ΔV_b is the barrier
  at the actual R(N_b) and H(N_b). N₀ = N_obs + 1, or 60 for comparison. For constant H this is higgs_reheat's
  (4π/9)e^{3N₀−B}; reproduced, B = 1460.7 tree and 2402.0 one loop at (172.60, 125.13).
- **The integral is dominated by ε_H → 1.** Mantziris–Markkanen–Rajantie arXiv:2011.03763 [Q] also find decay "most
  likely a few e-foldings before the end of inflation". Their Starobinsky follow-up arXiv:2207.00696 [Q] gets ξ ≳ 0.1 and
  notes sensitivity to the end-of-inflation dynamics.

| M_t | B tree at ε_H = 1 / constant H | log₁₀⟨N⟩ tree / loop (H(N)) | ξ_EW,min tree / loop, H(N) | same, constant H |
|---|---|---|---|---|
| 172.28 | 684 / 2425 | −226 / −734 | 0.060 / 0.040 | 0.043 / 0.034 |
| 172.57 | 430 / 1580 | −115 / −275 | **0.070 / 0.056** | 0.048 / 0.041 |
| 172.86 | 310 / 1160 | −63 / −138 | **0.079 / 0.067** | 0.053 / 0.047 |

- **Robustness of these numbers.** Values are for HT6t at N₀ = N_obs + 1 = 55.8. T1 and E1 agree to ±0.0005. N₀ = 60
  adds +0.002, and H_i = 1.50×10¹¹ adds +0.001.
- **Validity of Hawking–Moss at ε_H → 1** is the literature's assumption [A]. The non-equilibrium alternative is the
  Gaussian count, from the mode runs, of Hubble patches beyond the instantaneous tree de Sitter barrier. It gives
  log₁₀ ≤ −698 (HT6t, +1σ) and −911 to −1969 otherwise. Hawking–Moss is the conservative reading.
- **The A-term** (sign chosen so that Aφ_I < 0). At ξ_EW = 0.1 the bound is |A| ≤ 15.9 / 12.0 TeV (HT6t, 172.57 / 172.86),
  36.5 / 33.6 TeV (T1) and 38.4 / 35.3 TeV (E1) at tree level; the one loop is 10–35% looser. A λ_IH = 10⁻¹⁰ portal
  drives log₁₀⟨N⟩ to −10⁴ to −10⁷.

## 3. After inflation, gravitational only [C] (`bench_D.json`, `xiscan_*.json`, `upper_*.json`)

**D₆₀ at ξ_EW = 0.1** (M_h = 125.20; τ_c is the start of the decisive roll, from S1):

| history | 172.28 | 172.57 | **172.86** | τ_c (172.86) | M_t,crit N = 60 / N_obs | criterion band | M_h 125.09 / 125.31 |
|---|---|---|---|---|---|---|---|
| M1r surrogate (3, 0.5) | −1617 | −310 | −5.1 | 13.0 | 172.874 | 172.82–172.98 | 172.82 / 172.94 |
| HT6t MD | −2024 | −365 | **−23.8** | 18.7 | **172.917 / 172.936** | 172.89–173.00 | 172.86 / 172.98 |
| T1 MD / MDf | −2023 | −371 / −366 | −23.9 / −22.8 | 19.8 | 172.917 / 172.936 | 172.89–173.00 | 172.85 / 172.98 |
| E1 MD / MDf | −2023 | −371 / −362 | −24.0 / −22.0 | 19.8 | 172.917 / 172.936 | 172.89–173.00 | 172.86 / 172.98 |
| HT6t, T1, E1 G05 | −1966 to −1974 | −312 | −6.2 | 14 | 172.875 / 172.884 | 172.82–172.97 | 172.82 / 172.94 |
| HT6t MD, H_i = 1.50×10¹¹ | −1735 | −286 | −4.0 | 22.7 | 172.871 / 172.890 | 172.85–172.95 | 172.82 / 172.94 |
| T1 / E1 MDf, 1.50×10¹¹ | −1733 / −1714 | −283 / −280 | −2.9 / −3.5 | 22.7 | 172.868–172.870 / 172.887 | 172.85–172.95 | — |
| HT6t G05, 1.50×10¹¹ | −1643 | −238 | **+12.1** | 16.6 | 172.831 / 172.840 | 172.78–172.93 | 172.77 / 172.90 |

**Reading.**
- The three potentials are interchangeable for the Higgs at ξ ~ 0.1. What matters is the expansion history: a matter era
  with ⟨R⟩ = 3H², against early radiation with ⟨R⟩ → 0.
- The matter era delays the decisive roll (τ_c ≈ 19 against 13–14) and gains ≈ 0.04 GeV in M_t,crit.
- In σ: 172.92 is +1.20σ, 172.94 is +1.26σ, and 172.875 is +1.05σ.

**Deadlines D(τ_th) = 0** (sharp plasma onset, N = 60, τ from S1):

| M_t | 172.90 | 172.95 | 173.00 | 173.05 | 173.15 | 173.25 | 173.45 |
|---|---|---|---|---|---|---|---|
| MD (all three) | none | 16.0–16.7 | 12.8–13.1 | 10.8–11.0 | 8.5–8.7 | 7.3 | 5.6–5.7 |
| G05 | 12.0–12.5 | 9.6–9.8 | 8.4–8.6 | 7.5–7.7 | 6.5–6.6 | 5.7–5.9 | 4.8–5.0 |
| M1r surrogate | 10.7 | 8.2 | 6.8 | 6.0 | 4.8 | 4.2 | 3.3 |

**Lower edge of ξ_EW** (post-inflation, D₆₀ = 0; band in brackets):

| | 172.28 | 172.57 | 172.86 | 173.15 |
|---|---|---|---|---|
| HT6t / T1 / E1, MD | 0.029 | 0.067 [0.063–0.069] | **0.095** [0.090–0.097] | 0.114 |
| HT6t G05 | 0.030 | 0.068 | 0.098 [0.091–0.103] | 0.119 |
| HT6t MD, 1.50×10¹¹ | 0.033 | 0.071 | 0.099 [0.093–0.101] | 0.116 |

The higgs_reheat skeptic's canonical edges for the surrogate were 0.077 / 0.103 at (172.60 / 172.87, M_h = 125.13).

**Upper edge** (large ξ(H_i); start N_b0 = 2–3; run to τ = 150).
- **Diagnostic.** The produced variance is reconstructed from the final occupation numbers, where the modes are adiabatic.
  The instantaneous subtracted variance is contaminated by the 2ω oscillation of ξR.
- **Mean field, 2σ_prod = h_max, at central (+1σ) M_t.** HT6t ξ(H_i) ≈ 29 (11); T1 ≈ 70 (37); M1r surrogate ≈ 13 (5). The
  surrogate value is consistent with higgs_reheat's 2σ = 0.71 h_max at ξ = 10.
- **Scale-resolved D** stays < 0 through ξ(H_i) = 300.
- **ξ(H_i) = 1000 runs away** (Hartree) at τ ≈ 1.1 on HT6t and T1.
- In ξ_EW these are ×0.845.

**Numerical controls** (HT6t MD; ΔD at 172.86 / 172.92 / 173.00):

| change | ΔD |
|---|---|
| Hartree evaluated at 172.92 instead of 172.57 | +1.2 / +1.1 / +1.0 |
| N_b0 = 3 | +0.1 |
| 400 modes | +0.3 / 0 / 0 |
| K_max = 15 | +1.0 / −0.4 / −0.3 |
| rtol 10⁻¹⁰ | 0 |
| Hartree off | +2.8 |

M_t,crit therefore moves by ≤ 0.004 GeV. The criterion conventions (±0.05–0.08 GeV) and M_h (±0.06) dominate.

## 4. Plasma timing on the benchmarks [E, rates C] (`plasma.json`, `final_plasmaD.json`, `final_portal.json`)

**Model.** This is inf_reheat's plasma_solve, put onto the real oscillation.
- Greene–Kofman production of Majorana N at each zero of m_N = M_N + yφ_I(t), with the real φ̇ at the crossing and a
  Pauli refill rule.
- Perturbative φ_I → NN times the open fraction.
- N → LH at y_ν²|m_N|/8π, with y_ν² = m̃M_N/v² (m̃ = 0.05 eV).
- The hard L, H gas: m_ind² = 0.0228 ρ_SM/E², with E = ρ_SM/n_SM. The coefficient is fixed by κ_T = 0.11 for a thermal SM gas.

y_nat = 2.84×10⁻⁴ (HT6t), 4.30×10⁻⁴ (T1) and 4.40×10⁻⁴ (E1), from inf_reheat's flatness criterion at the benchmark's φ_I*.

| | y | τ(m_ind² ≥ 10⁻³) M_N = m_I/4 / m_I/16 | m_ind²(τ = 8) (H_i²) | A_ind (TeV), m_I/4 / m_I/16 | D < 0 up to (GeV) |
|---|---|---|---|---|---|
| HT6t (A_max 12 TeV) | 10⁻⁵ | 1.66 / 3.36 | 0.005 / 0.002 | 0.4–11 / 0.02–0.8 | 173.39 / 173.13 |
| | 3×10⁻⁵ | 1.10 / 1.24 | 0.04 / 0.013 | 1.2–34 / 0.07–2.3 | > 174.6 |
| | 10⁻⁴ | 1.05 / 1.14 | 0.22 / 0.10 | 3.9–112 / 0.24–7.7 | > 174.6 |
| | y_nat | 1.00 / 1.05 | 0.64 / 0.25 | 11–319 / 0.7–22 | > 174.6 |
| T1 (A_max 34 TeV) | 10⁻⁵ | 1.17 / 1.40 | 0.016 / 0.006 | 3.6–95 / 0.2–6.5 | > 174.6 / 173.44 |
| | 3×10⁻⁵ | 1.01 / 1.08 | 0.15 / 0.055 | 11–284 / 0.7–20 | > 174.6 |
| | y_nat | 0.97 / 0.97 | 21 / 7 | 153–4069 / 9.6–281 | > 174.6 |
| E1 (A_max 35 TeV) | 10⁻⁵ | 1.09 / 1.20 | 0.03 / 0.013 | 14–359 / 0.9–25 | > 174.6 |
| | 3×10⁻⁵ | 0.99 / 1.02 | 0.28 / 0.11 | 43–1077 / 2.7–75 | > 174.6 |
| | y_nat | 0.95 / 0.96 | 60 / 21 | 627–15800 / 39–1096 | > 174.6 |

**Blocking and channel character.**
- **HT6t is blocked for y ≳ 3×10⁻⁵.** yΦ₁/(m_I/2) = 9.5 at y_nat on the steep side (×2 on the hill side), so its early
  plasma is instant preheating, as in inf_reheat.
- **T1 and E1 are perturbatively open for y ≲ 10⁻⁴.** After inflation their amplitude is only 0.03 M_P (E1's steep side
  0.006 M_P), against m_I = 128–256 H_i.
- **Implied widths.** Γ_I(open) = y²m_I/16π = 5×10⁻¹¹–6×10⁻⁷ H_i, so T_RH ~ 10⁹–10¹¹ GeV. A lower T_RH lowers N_obs and
  makes D slightly more negative (−1.3 per e-fold).

**What the plasma result rests on [E/A].**
- **The hard-gas induced mass.** A fresh L + H gas couples to the Higgs condensate through g², g′² and λ (negative at these
  scales). Gauge bosons and tops appear through splitting.
- **The thermalisation time.** inf_reheat's pessimistic estimate puts a fully thermal bath beyond τ = 60 (≈ 2×10³ in
  inf_reheat). Harigaya–Mukaida 1312.3097 [Q] argue for near-instant thermalisation. If only a thermal bath counts, §3
  applies unchanged.

**The R-even portal** (inf_reheat skeptic's convention m_H² = λφ_I²/2).
- λ_min, q and λ_max are as in row 9.
- Mode runs at λ_min: D ≤ −4200 (HT6t) and ≤ −3300 (T1, homogeneous background, which keeps ⟨φ_I²⟩) at every
  M_t ∈ [172.57, 174.6].
- At 0.1 λ_min: D ≤ −820.
- The S1 loop itself induces λ_ind ≈ y²y_ν²/16π² × [1, ln] ~ 10⁻¹²–2×10⁻¹⁰ at y_nat, comparable to λ_min [E].

## 5. What is not established

- **The criterion D** is argued, not simulated, and it is steep (dD/dM_t ≈ 750 per GeV near the edge). Its convention band
  moves M_t,crit by −0.03/+0.08 GeV (matter era) or −0.05/+0.10 GeV (Γ = 0.5). The verdict at +1σ (−22 to −31; band
  [−71, −12]) holds across the band in the matter era, but not for Γ = 0.5 (band reaches +14). The literature's
  end-of-inflation Hubble-patch count (Kearney–Yoo–Zurek 1503.05193; Espinosa et al. 1505.04825 [Q]) is more
  pessimistic, as higgs_reheat noted.
- **The inflation side** uses the de Sitter Hawking–Moss rate up to ε_H = 1 [A]. The one-loop curved-space potential is
  taken at the instantaneous curvature (H_eff² = R/12), which is a choice; the tree/loop spread (0.070 vs 0.056) is the
  scheme error.
- **The homogeneous T1/E1 backgrounds** after fragmentation (0.35–0.43 H_i⁻¹ after the end) are indicative. The w = 0
  switch changes D by ≤ 1.2. The upper edges on T1 use the homogeneous oscillation.
- **The plasma model** (§4) is inf_reheat's [E], with the composition and thermalisation caveats above.
  - S2 (axial) was not redone on these benchmarks; inf_reheat's τ ≈ 0.02–0.03 rests on an EFT-strained estimate.
  - Non-perturbative N production in fragmented T1/E1 is taken from the homogeneous oscillation.
- **A_ind** is a natural-size estimate. Its value is a renormalisation condition, like the hidden-sector induced portals.
- **The calibration loop is not closed.** H_i was not recalibrated for the slow-decay histories. inf_reheat and the
  inf_potential skeptic show that slow decay moves the dark-matter-calibrated H_i (for the T-model stand-ins to
  ~6×10¹² GeV). The Higgs numbers at a different H_i must be redone: h_max/H_i changes.
- **Not included.** M_t and M_h enter only through the SM fit grid (M_t 170.6–174.6). Values above 174.6 are
  extrapolated. α_s was not varied.

## 6. Reproduction (one core; ≈ 25 min in total)

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/inf_higgs/code; export OMP_NUM_THREADS=1
python3 -B ih_run.py bench                       # ~10 min (two chunks; resumable)
python3 -B ih_run.py xiscan HT6t MD 0; python3 -B ih_run.py xiscan T1 MDf 0; python3 -B ih_run.py xiscan E1 MDf 0
python3 -B ih_run.py xiscan HT6t G05 0; python3 -B ih_run.py xiscan HT6t MD 1
python3 -B ih_run.py upper HT6t MD 0; python3 -B ih_run.py upper T1 MD 0 10,30,100,300,1000,3000
python3 -B ih_run.py upper M1r G05 0 1,3,10,30
python3 -B ih_inflation.py main; python3 -B ih_inflation.py tail; python3 -B ih_inflation.py Aterm
python3 -B ih_plasma.py; python3 -B ih_final.py plasmaD; python3 -B ih_final.py portal; python3 -B ih_final.py controls
python3 -B ih_summary.py
```

Mode runs and SM tables are cached in
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/inf_higgs/` (`runs/*.npz`,
`smcache/`), together with the full plasma curves (`plasma_full.json`).

**Credits.**
- SM inputs: Buttazzo et al. 1307.3536; SMDR, Martin–Robertson 1907.02500 (through cf_vacuum and higgs_reheat).
- Curved-space potential: Markkanen–Nurmi–Rajantie–Stopyra 1804.02020.
- Inflationary bound with H(N): Mantziris–Markkanen–Rajantie 2011.03763 and 2207.00696.
- Post-inflation: Herranen et al. 1506.04065; Figueroa–Rajantie–Torrenti 1709.00398; Ema–Mukaida–Nakayama 1602.00483.
- Preheating: Greene–Kofman hep-ph/9807339; Felder–Kofman–Linde hep-ph/9812289; Harigaya–Mukaida 1312.3097.
- Oscillons: Lozanov–Amin 1608.01213 [R].
- Our tracks, read-only: higgs_reheat and its skeptic, inf_potential, inf_reheat and its skeptic.
