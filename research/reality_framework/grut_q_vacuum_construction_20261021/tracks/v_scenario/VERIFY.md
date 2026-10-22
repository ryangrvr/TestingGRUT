# Verifier `v_scenario`: independent check of the scenario track (QV1 solved example)

22 October 2026. Stage 2 of the q-theory construction pass (issue #13). Briefs: `../../BRIEF.md`,
`../../BRIEF_STAGE2.md`. Object of the check: `../scenario/NOTES.md`, its `code/` and `results/` (cited "scenario
§n"). Stage-1 records consulted: `../equations/NOTES.md`, `../literature/NOTES.md`, `../ledger/NOTES.md`.

**Independence.** I read the scenario's code to learn its toy parameters and conventions. Every number below
comes from my own code in `code/`, written from scratch, and I did not import or copy the scenario's code. My
formulations differ on purpose:
- H is solved from the Friedmann constraint at every step. The scenario integrates dH/dt and relies on Bianchi
  propagation.
- Initial data and the quasi-static track come from my own coupled (y, H²) adiabatic series to O(m⁻¹⁰). The
  scenario's track is O(m⁻⁴).
- The light-q relic uses a linearised Einstein-frame oscillator in τ = m t.
- θ* uses a distance ratio at fixed r_s. The scenario uses fitting formulas.
- The turnaround times are integrated from T = 10¹⁰ GeV. The scenario starts them at 1 TeV.
- The public data were downloaded again by me, to the scratchpad only. They are byte-identical to the scenario's
  copies.

**Labels.** [D] derived here (sympy); [C] computed here; [Q] read in the source at the level stated; [R] recalled.
Every comparison with data is retrospective. Neither track claims anything as held out.

---------------------------------------------------------------------------------------------------------------------

## 1. Verdict table

| # | Claim (scenario §) | Scenario | Verifier | Verdict |
|---|---|---|---|---|
| 1 | FRW constraint, space-space and q-equations follow from the action (§1.1, 1.4) | [D] | Minisuperspace Euler-Lagrange: E_N = a³C, E_a = 3a²(space), E_q = −a³(q-eq); residuals 0 [D] | CONFIRMED |
| 2 | Bianchi: dC/dt = −4HC; with dissipation the constraint holds iff Q = Γ_q q̇² (§1.4, 1.6) | residual 0 | dC/dt + 4HC = Γ_q q̇² − Q for the trace-form system; −3HC + (same) if the space-space equation replaces the trace [D] | CONFIRMED |
| 3 | KV forms: q0 = 3M², μ0 = −1/(3χq0), ε − μ0q = 0 at q0, m_q = √(2/3) M = 1.988e18 GeV, b3 = 4/3, b4 = 1/3, β = 1; Γ_q = (3/2)ΓK0K'²/K² damps the energy at Γ (§1.1, 1.4, 1.6) | as stated | identical; m_q = 1.98817e18 GeV; damping coefficient Γ K0/K [D] | CONFIRMED |
| 4 | R1: 3K0H² = ρ + [S²/(12K0) − H²S − HṠ]/m² + O(m⁻⁴) (§1.4) | [D] | the same from the coupled expansion; shape coefficients enter only at O(m⁻⁴) [D] | CONFIRMED |
| 5 | ρ_V = S²/(12K0m²) ≡ KV (9/8)κ²H⁴K'²/ε'' (§2.4) | [C] | ratio = 1 identically [D] | CONFIRMED |
| 6 | eps0 equals the predicted R1 correction; R1 removes it (§2.2), m/H = 101 / 303 / 1010 / 3029 | eps0 2.83e-6 / 3.129e-7 / 2.743e-8 / 3.048e-9; eps1 7.0e-8 / 2.0e-10 / 1.2e-11 / 2e-11; fit 0.0298 (H/m)^2.007 | eps0 2.73e-6 / 2.908e-7 / 2.769e-8 / 3.077e-9 (predicted 2.743e-6 / 2.909e-7 / 2.769e-8 / 3.077e-9); eps1 6.1e-8 / 1.1e-10 / 8.5e-12 / 3.2e-11; fit 0.0262 (H/m)^1.991 [C] | CONFIRMED. The coefficient is the max over windows from N ≥ −0.25, so it depends on window placement. It is not a physical constant. |
| 7 | Deposit (4/3)Δ²/(Km²)(x/sinh x)², x ≈ 0.1, 1, 2, 4 (§2.3) | m/H = 100: 0.99997, 0.99969, 0.99983, 1.0072; 300: 1.0000072, 0.99990, 0.99993, 1.0012 | m/H = 100: 1.00015, 0.99990, 1.00012, 1.0081; 300: 1.000005, 0.99990, 0.99992, 1.0012 (R1 H_c, the scenario's convention) [C] | CONFIRMED to ≤ 1e-3 (H_c convention, §2.3) |
| 8 | Deposit at x = 8, m/H = 100 (§2.3, §0 row 3) | 1.028 | **1.077 ± 0.005** (1.075–1.082 over three initial-data variants; rtol 1e-12/1e-13 identical) [C] | **CORRECTED** |
| 9 | x = 8 at m/H = 300; decay exponent κ (x ≥ 2) | 1.0117; κ = 1.999 (100), 2.003 (300) | 1.0111; κ = 1.993 (100), 2.004 (300) [C] | CONFIRMED |
| 10 | Δ² scaling (§2.3) | 0.9989, 0.9997, 0.9985 | 1.00039, 1.00038, 1.00042 (R0 H_c) [C] | CONFIRMED (better than stated) |
| 11 | KV relaxation, w_M = 1/3: t⟨H⟩, t²⟨ρ_V⟩ (§3.1) | 0.6655, 0.6620 | 0.66549, 0.66195; ⟨ρ_V⟩/(3⟨H⟩²) = 0.498 (scenario 0.496) [C] | CONFIRMED |
| 12 | Γ = 4m³/(192πM²) = 4.42e-3 m_q; GP formula T_reh = 5.03e16 GeV (§1.6, 3.2) | as stated | 4.4210e-3 m_q; 5.030e16 GeV. GP 1009.2448 eqs (T4+), (T5+) read in TeX; M_P reduced there [Q] | CONFIRMED |
| 13 | Decay runs, y_i = −0.1, −0.03, +0.03, +0.1, +0.25: f_inj; T_rh (§3.2) | 1.062, 1.016, 0.986, 0.959, 0.917; 5.05, 4.82, 4.83, 5.06, 5.12 e16 | 1.0620, 1.0160, 0.9859, 0.9585, 0.9173; 5.049, 4.819, 4.827, 5.060, 5.124 e16; E decays at 1.0000005 Γ; late tH = 0.511 [C] | CONFIRMED |
| 14 | Natural-m_q epoch table (§4.1, §0 row 5) | 1 TeV: 2.6e-64, 6.9e-55, −2.0e-64, 4.8e-44; EW: 1.4e-65, 6.8e-58, −6.2e-66, 5.8e-47; BBN: 3.2e-92, 1.4e-111; today: 1.8e-121, 4.7e-169, −1.6e-121, 1.9e-131 /yr | 2.59e-64, 6.94e-55, −1.99e-64, 4.77e-44; 1.39e-65, 6.82e-58, −6.23e-66, 5.79e-47; 3.18e-92, 1.39e-111; 1.84e-121, 4.72e-169, −1.64e-121, 1.90e-131 [C] | CONFIRMED. At 10 MeV the derivative columns (H/H_GR − 1, Ġ/G) differ by about 1.8× because S_N straddles the 16% table/ideal-gas junction. They are not load-bearing. |
| 15 | Light 4 meV relic Ω_q h² (§4.2, §0 row 6) | 378 (adiabatic, 1 TeV); q0 at rest at 1/10/100/1000 TeV: 502 / 4.1e4 / 1.9e5 / 4.1e5; S/ρ = 0 above 300 GeV: 115; +1σ: 2.1e6; 1 eV: 2.2e3 / 6.2e5 | 375 (first-order track, as in the scenario's code); 499 / 4.12e4 / 1.86e5 / 4.04e5; 117; 2.10e6; 2.23e3 / 6.24e5; **O(m⁻⁴) track at 1 TeV: 4.0e3** [C] | CONFIRMED. The "adiabatic start" at m/H = 2.9 is ill-defined to a factor 10.7. The verdict needs only Ω > 0.12. |
| 16 | Decay before BBN needs m_q ≳ 2e5 GeV; Γ/H(1 MeV) = 1.6e-49 | 2e5 [E]; 1.6e-49 | 1.85e5 GeV; 1.62e-49 [C] | CONFIRMED |
| 17 | QV1 late: H0, age, q0, tensions, 100θ*, forced-H0 case (§5.1) | 37.74, 17.27 Gyr, +0.500, −34σ, +42σ, 0.943 (−327σ), ω_m 0.533, 8.92 Gyr (−8.8σ) | 37.741, 17.267, +0.5001, −33.9σ, +42.1σ, 0.9430 (−327.1σ, ratio method), 0.5334, 8.924 Gyr (−8.77σ) [C] | CONFIRMED |
| 18 | QV1 SN, BAO χ² (§5.1) | SN 2101.4 (Δχ² 699); BAO shape 1421.9 (Δχ² 1412); BAO abs 15333 (r_d 147.04) | 2101.44 (698.5); 1421.87 (1411.6); 15297 (r_d = 147.09 [R]) [C] | CONFIRMED |
| 19 | QV1-cal: H0, age, q0, θ*, χ² (§5.2) | 67.36, 13.81, −0.529, 1.0410, SN 1403.88, BAO shape 13.80, abs 33.4 | 67.364, 13.813, −0.5291, 1.04111, 1403.88, 13.80, 32.4 (r_d 147.09) [C] | CONFIRMED |
| 20 | QV1-cal second-order term (§0 row 8, §5.2) | (4/3)ρ_Λ/(Km²) = 1.4e-120; "confirmed" by the fixed-μ toy run, "same algebra" | μ-offset (type ii): ρ_V = c[1 + c/(K0m²) + …] → **1.08e-120** (mpmath). Toy m/H = 100: type ii 3.0018e-6 = c/m²; type i 4.0018e-6 = (4/3)c/m² [D][C] | **CORRECTED** |
| 21 | Initial-equilibrium EW reading (§5.3) | 44.6 GeV, 2.02e-10 s, crunch 4.05e-10 s | 44.60 GeV, 2.024e-10 s, 4.047e-10 s [C] | CONFIRMED |
| 22 | 10 TeV scalar, fixed μ, c = −2.37e13 GeV⁴ (§5.3, §0 row 9) | T_turn 912 GeV, t 4.0e-13 s, crunch ≈ 8e-13 s | 912.4 GeV, **t_turn 4.51e-13 s, crunch 9.01e-13 s** [C] | **CORRECTED** |
| 23 | +c de Sitter; −ρ_Λ turnaround (§5.3) | 1.16e-12 GeV, 912 GeV; a 0.770, 3.54 K, 18.4 / 36.7 Gyr | 1.155e-12, 912.4; 0.7705, 3.537 K, 18.35 / 36.70 Gyr [C] | CONFIRMED |
| 24 | Fixed-μ static roots ρ_V,∞/c − 1 (§5.3) | 1.4e-120, 6.8e-66, 1.35e-60 | 1.4336e-120, 6.7631e-66, 1.3510e-60 (400 digits; checked against a cancellation-free form) [C] | CONFIRMED |
| 25 | P_q re-absorption: (q0'/q0)² = (1 + √(1+8χc))/2; δK/K = χc; m² +5χc; m²/K +4χc; c ≥ −M⁴/8 (§5.3) | δK/K 3.4e-66, 6.8e-61, 7.2e-121 | identical [D]; δK/K = χc − (5/2)χ²c²; 3.38e-66, 6.75e-61 (in magnitude), 7.17e-121 | CONFIRMED |
| 26 | Sign of δG/G for the 10 TeV scalar under P_q (§0 row 9) | −6.8e-61 | **+6.8e-61** (c < 0 and δG/G = −χc; the scenario's JSON has the right sign) | **CORRECTED** |
| 27 | "Relative size of the offset 7.2e-121" (§5.2) | 7.2e-121 | c/M⁴ = 7.17e-121; but δμ/μ0 = 3c/M⁴ = **2.15e-120** for the KV quartic (q0\|μ0\| = M⁴/3) | CLARIFICATION |
| 28 | Toy dust fraction at N = 0 (§1.5) | "6% of ρ" | **2.05%** (ρ_d = 0.06, ρ = 3H_c² = 2.925) | **CORRECTED** (documentation) |
| 29 | Saikawa-Shirai p(T = 0) = 0; GP rate with reduced M_P (§1.3, 1.6) | [Q] | both read in TeX at section level [Q] | CONFIRMED |
| 30 | Deposit floor "always below the smallest deposit" (§2.3) | E_floor/E_sudden = 2.5e-8 (m/H = 100) | True in energy. But the floor enters the ratio through amplitude interference: ±2√(floor/deposit) = **±6% at x = 8**, ±1.1% at x = 6 (m/H = 100) | **CORRECTED** (error estimate) |
| 31 | Credit for §4.2's mechanism | moduli problem, Cembranos | The trace-driven "kick" of a light T-coupled (chameleon/f(R)) scalar is Brax et al., astro-ph/0408415, and Erickcek, Barnaby, Burrage, Huang, arXiv:1304.0009 (PRL 110, 171101) [Q abstract] | MISSING CREDIT |
| 32 | QV1 H0 and age in §4.1 vs §5.1 | 37.72 / 17.28 vs 37.74 / 17.27 | 37.741 / 17.267 | CONFIRMED (§5.1 values; §4.1 differs at 5e-4 from grid quadrature) |
| 33 | Verdict-level conclusions (§0 rows 5–10, §6) | — | every number they rest on is reproduced | CONFIRMED |

**Not checked:**
- the §2.4 numerical window ratios (the identity is checked symbolically only);
- the §5.4 toy runs other than the two offsets (initial-equilibrium turnaround t = 3.51/H_c, binned H with c re-absorbed);
- the §3.1 w_M = 0 cases (the literature track reproduced them);
- the §3.2 runs with KV boundary data and with y_i ≤ −0.5 or ≥ +0.5;
- the §2.2 rows at m/H = 10 and 30. These are dominated by initial data, and my higher-order series start is worse
  there (asymptotic series at m/H_i ≈ 6), so they are not comparable.

**Net.** No verdict of the scenario track changes:
- QV1 fails the exposed late-universe data retrospectively.
- QV1-cal reproduces ΛCDM only because ρ_Λ is calibrated in.
- The light branch overcloses for untuned initial data.
- P_q absorbs every constant exactly.

The corrections are second-order or bookkeeping details.

---------------------------------------------------------------------------------------------------------------------

## 2. Details of the checks

### 2.1 Equations [D] — `code/v1_symbolic.py` → `results/v1_symbolic.json`

- **A (action → FRW).** L = N a³[K(q)R/2 − ε(q) + μq] − N a³ρ(a), with R for lapse N. The Euler-Lagrange equations
  in N, a and q equal +a³C, 3a²·(space-space) and −a³(ε' − K'R/2 − μ), where
  C = 3H²K + 3HK̇ − ρ − (ε − μq). Generic K and ε; residual 0.
- **B (Bianchi, dissipation).** I use the system the scenario integrates:
  - Ḣ from ε' − K'R/2 = μ − Γ_q q̇;
  - q̈ from the trace K̈ + 3HK̇ = (S + 4ρ_V − KR)/3;
  - ρ̇ = −3H(ρ+p) + Q.

  This gives dC/dt + 4HC = Γ_q q̇² − Q exactly. If the space-space equation replaces the trace equation, the result
  is dC/dt + 3HC = Γ_q q̇² − Q. Either way the constraint surface is invariant iff Q = Γ_q q̇².
- **C (KV quartic, M = 1).**
  - ρ_V(Y) = Y²/2 + 2Y³/3 + Y⁴/6 + c(1+Y), i.e. (3/4)m²Y²(1 + 4Y/3 + Y²/3) + c(1+Y) with m² = 2/3.
  - R(Y) = 2Y + 4Y² + 4Y³/3 + 2c.
  - The trace force is S/3 + (c/3)(2+2Y) − m²F(Y), with F = Y + 2Y² + 4Y³/3 + Y⁴/3.
  - The damping coefficient is (2/3)KΓ_q/K'² = ΓK0/K.

  All match scenario §1.4.
- **D (R1).** I expand y = Σeᵏy_k (e = 1/m²) in q-eq + constraint. This gives y1 = S/3 and
  y2 = −2S²/9 − S̈/3 − HṠ, and the O(e) constraint 3H² = ρ + e[S²/12 − H²S − HṠ].
- **G (P_q, KV quartic).** P(q) ≡ qε' − ε = (1/2χ)(u² − u), with u = q²/q0². Then:
  - u = (1 + √(1+8χc))/2;
  - δK/K = χc − (5/2)χ²c²;
  - m² ∝ 1 + 5χc and m²/K ∝ 1 + 4χc.
- **H (static roots).** Type i (ε → ε + c, μ = μ0) gives ρ_V = c + 2c² − (32/3)c³ (M = 1). Type ii
  (μ = μ0 − c/q0, QV1-cal) gives ρ_V = c + (3/2)c², i.e. c[1 + c/(K0m²)], since K0m² = 2/3.

  Generic check against the equations track (§5.2): c* = −c + c²/(2q0²ε''), so
  dc = c − c²/(2q0²ε''). The type-ii coefficient is therefore 2K1²/(K0²ε'') − 1/(2q0²ε'') = 2 − 1/2 = 3/2.

  Full-equation confirmation at m/H = 100 (`code/v7_offsets.py`, c = 0.03 H_c² units):
  - ρ_V,late/c − 1 = 4.0018e-6 for type i;
  - ρ_V,late/c − 1 = 3.0018e-6 for type ii.

### 2.2 Compressed-hierarchy toy [C] — `code/v2_toy_full.py`, `code/v2_qs_series.py` → `results/v2_toy_full.json`

The toy is identical to the scenario's:
- ρ_r = 3(1 − 0.02 − 0.05)e⁴ and ρ_d = 3·0.02·e³ at N_i = −1;
- ΔV = Δ(1 − tanh(N/w))/2, Δ = 0.15, with the latent heat going to radiation;
- runs at N_i = −1 (reduction) and −0.7 (deposit), Nf = 1.5 / max(1, 6w + 0.9);
- windows of max(1, 0.002m) periods, max taken over N ≥ −0.25.

The integration uses DOP853 with rtol 1e-12 (1e-13 identical). At N = 0, ρ_d/ρ = 0.06/2.925 = 2.05%, not the "6%"
of scenario §1.5.

### 2.3 Deposit measurement and the corrections at x = 8 [C]

E a³ is measured about my O(m⁻⁶) quasi-static track (Einstein-frame energy, averaged over one period, five
post-step windows). The spread is ≤ 8e-6. Δ = 0 controls give E/E_sudden = 1.4e-11 (m/H = 100) and 6e-21 (300).
The scenario's floor is 2.5e-8 and 4.4e-12.

- **H_c convention.** x = πm w/(2H_c). At the centre of a sharp step, H_c from R0 (0.99958), R1 (0.99997) and the
  full solution (0.99985) differ by up to 4e-4, through the large S_N term. This moves (x/sinh x)² by up to 8e-4 at
  x = 2. Agreement at 1e-4 (scenario: "1.0000 ± 0.0004") therefore presupposes the R1 convention. Convention-free,
  the formula holds to ~1e-3 for x ≤ 2 at m/H ≥ 100.
- **x = 8, m/H = 100.** I get 1.0765 (N_i = −0.7, O(m⁻¹⁰) start), 1.0824 (O(m⁻⁶) start) and 1.0763 (N_i = −0.5).
  The scenario gets 1.028. The scenario's own floor (2.5e-8 E_sudden) interferes with the deposit amplitude
  (2.9e-5 E_sudden at x = 8), giving ±2√(2.5e-8/2.9e-5) ≈ ±6%. That covers the gap. My floor gives ±0.1–0.4%. The
  corrected value is **1.077 ± 0.005**. The qualitative statement, a finite-H correction vanishing with m/H, is
  unchanged: 1.077 (100), 1.011 (300). The scenario's value at 300 (1.0117) is confirmed.

### 2.4 Realistic scales [C] — `code/v3_realistic.py` → `results/v3_realistic.json`

The background uses:
- the SS tables (g_ρ, g_s from the arXiv ancillary file; w from the authors' page);
- entropy conservation;
- a dust density (ω_b + ω_c)ρ_c,100;
- below 10 MeV, my own Fermi-Dirac trace for e± and μ± at the photon temperature;
- R1 evaluated with mpmath (260 digits).

The quantities follow from R1:
- H/H_GR − 1 = √{(1 + S²/(12K0m²ρ))/(1 + (S+S_N)/(3K0m²))} − 1;
- δG/G = −S/(3K0m²);
- Ġ/G = −H S_N/(3K0m²).

Today (dust): H/H_GR − 1 = (9/24)ρ_d/(K0m²) = 1.84e-121.

**10 TeV threshold.** With c = −2.37e13 GeV⁴, the scenario starts t at 1 TeV with t = 1/(2H_GR) = 2.37e-13 s,
ignoring c, which is 68% of ρ there. The true t(1 TeV) in that universe is 2.80e-13 s (integrated from 10¹⁰ GeV).
Hence t_turn = 4.51e-13 s and the crunch = 9.01e-13 s. The analytic value for pure radiation at constant g*,
(π/4)√(3K/|c|) = 4.47e-13 s, agrees. The EW case is unaffected (c/ρ(1 TeV) = 3e-6): 2.024e-10 s against the
scenario's 2.023e-10 s.

### 2.5 Light q [C] — `code/v4_light_q.py` → `results/v4_light_q.json`

The equation is y'' + 3(H/m)y' + y = (S/ρ)(H/m)², with τ = mt, on the GR + SS background. |y| ≤ 2.6e-3 in every run,
so the dropped KV nonlinearities are below 1%. E = (3/4)K0m²[(y' − σ')² + (y − σ)²]. E a³ is constant to ≤ 4e-4
from m/H = 300 to 2850, and is then carried to today.

All of the scenario's values are reproduced to 1–3%. One sensitivity is new: the "adiabatic start at 1 TeV" (m/H =
2.9) gives 375 with the first-order track, which is what the scenario's code uses, and 4.0e3 with the O(m⁻⁴) track.
The start is not adiabatic there, and the number is defined only by convention. The verdict (Ω_q ≫ 0.12 for untuned
data; the amplitude needed for Ω_q h² < 0.12 is ≲ 5e-7) is unaffected.

### 2.6 Planck era [C] — `code/v5_early.py` → `results/v5_early.json`

The constraint-solved formulation includes the dissipative source. The KV run uses Y_i = −0.5, H_i = 0.5/√3 and
ρ_M = 0.05, with t = √3 τ_KV. GP 1009.2448 was read in TeX (eq. (T4+) Γ = μ³/(192πM_P²) with M_P reduced; eq. (T5+)
T_reh ≈ 4.5e-2 g*^(−1/4)(N_s μ³/M_P)^(1/2)). Note that GP define T_reh at equality of the inflaton and matter
densities, while the scenario reads ρ_r at t = 1/Γ. The two definitions agree to 5% here.

### 2.7 Data [C] — `code/v6_data.py` → `results/v6_data.json`

- **Pantheon+.** 1590 SNe with z_HD > 0.01, STAT+SYS covariance, M profiled, d_L = (1+z_hel)D_M(z_HD). The flat
  ΛCDM best fit is Ω_m = 0.3316, χ² = 1402.92.
- **DESI DR2.** The flat ΛCDM shape best fit is Ω_m = 0.2975, χ² = 10.27. The absolute BAO uses r_d = 147.09 Mpc
  (Planck 2018 VI, [R]), against the scenario's Aubourg-formula value of 147.04 Mpc. That difference explains the
  absolute χ² differences (15297 vs 15333; 32.4 vs 33.4).

---------------------------------------------------------------------------------------------------------------------

## 3. Audit items (conventions, calibration, credit, consistency)

- **Planck mass and 8π.** These are handled consistently: m_q, Γ (GP use the reduced M_P, confirmed in TeX), ρ_crit
  = 3M²H² and the T_rh conversions all use the reduced M = 2.435e18 GeV. The non-reduced alternative (χ⁻¹ = E_P⁴ →
  m_q × 8π) is stated correctly.
- **Calibration.** Only QV1-cal is calibrated (datum ρ_Λ,obs), and it is labelled so everywhere. QV1's H0 = 37.74
  uses ΛCDM-fitted ω_b, ω_c. The scenario says so, and the SN/BAO shape exclusions, which I confirm, do not depend on
  them.
- **Held-out vs retrospective.** No held-out claim is made. Every data comparison is labelled retrospective.
- **Regimes.** The reduction is applied only where m/H ≥ 300 (validated). The light branch is integrated in full. The
  §3 Planck-era numbers carry the stated EFT caveat. No approximation is used outside its stated regime.
- **Missing credit (§4.2).** The mechanism by which the SM trace drives a light, T-coupled f(R)/chameleon scalar off
  its minimum, followed by a relic oscillation, is published:
  - Brax, van de Bruck, Davis, Khoury, Weltman, PRD 70, 123518 (2004), astro-ph/0408415;
  - Erickcek, Barnaby, Burrage, Huang, PRL 110, 171101 (2013), arXiv:1304.0009 ("the chameleon couples to the trace
    of the stress-energy tensor").

  Both were read at abstract level. They should be credited next to the moduli-problem references.
- **Internal consistency with stage 1.** I found no contradiction:
  - equations §5.4's 30% gap at w = 0.1 belongs to a toy with m/H ≈ 40 and a high floor;
  - BRIEF_STAGE2's 5e-169 GeV⁴ matches 4.72e-169;
  - the ledger's 4 meV Eöt-Wash bracket is used as stated.

  One small internal inconsistency inside the scenario: §5.2 writes ρ_V = c + c²/(Km²) for QV1-cal (correct) but
  quotes the 4/3 value 1.4e-120 (row 20).

## 4. Inputs and provenance

| Input | Value | Source |
|---|---|---|
| M, ρ_Λ,obs | 2.435e18 GeV; 2.52e-47 GeV⁴ | brief [P] |
| ε(q), K(q) | KV quartic, K = q/3, χ = 1/M⁴ | scenario §1.1 / literature track [Q]; re-derived consequences here [D] |
| ω_b, ω_c, T0 | 0.0224, 0.120, 2.7255 K | as in the scenario (Planck 2018 VI abstract; Fixsen 2009) [Q], retrospective |
| SS EoS | standardmodel2018.dat (arXiv:1803.01038 ancillary), eos_w.dat (member.ipmu.jp/satoshi.shirai) | downloaded by me; p(T=0) = 0 normalisation read in the TeX, EW section [Q] |
| Γ, T_reh formula | GP 1009.2448 (T4+), (T5+), N_s = 4 | TeX, section level [Q] |
| Δ_EW, heavy c | m_H²v²/8 (125.25, 246.22 GeV); −(3/2)(10 TeV)⁴/(64π²) | as in the scenario [R] |
| DESI DR2, Pantheon+ | CobayaSampler/bao_data desi_bao_dr2; PantheonPlusSH0ES DataRelease | downloaded by me; byte-identical to the scenario's |
| SH0ES H0, q0; Planck θ*; GC age | 73.04 ± 1.04, −0.51 ± 0.024; 1.0411 ± 0.0003; 13.5 (+0.16/−0.14 stat, ±0.5 sys) Gyr | as quoted by the scenario (abstract level) |
| z*, r_d | 1089.92; 147.09 Mpc | Planck 2018 VI table 2 [R] (θ* ratio insensitive to z* at 1e-4) |
| h for the ΛCDM reference in θ* | 0.6736 | Planck 2018 [R] |

Nothing in this verifier is calibrated.

## 5. Reproduction and limits

    python3 -I code/run_all.py      # downloads to the scratchpad (v_scenario_dl), then v1..v7; ~6 min on 4 cores

| File | Content |
|---|---|
| `code/v1_symbolic.py` → `results/v1_symbolic.json` | minisuperspace; Bianchi with dissipation; KV forms; R1; kick identity; form factor; P_q; static roots (types i and ii) |
| `code/v2_qs_series.py` | coupled adiabatic series for y and H² to O(m⁻¹⁰) (sympy) |
| `code/v2_toy_full.py` → `results/v2_toy_full.json` | toy reduction scan, deposit scan, Δ² scaling, controls |
| `code/v3_realistic.py` → `results/v3_realistic.json` | natural-m_q epochs; QV1 / QV1-cal late universe; turnarounds; light-q bookkeeping |
| `code/v4_light_q.py` → `results/v4_light_q.json` | light-q relic |
| `code/v5_early.py` → `results/v5_early.json` | KV relaxation; decay runs |
| `code/v6_data.py` → `results/v6_data.json` | Pantheon+ and DESI DR2 χ² |
| `code/v7_offsets.py` → `results/v7_offsets.json` | type i vs type ii late vacuum at m/H = 100 |
| `results/verdicts.json` | the verdict table in machine-readable form |

**Limits.**
- My toy windows are not placed exactly as the scenario's. Max-over-window quantities therefore agree only to about
  10%, while the predicted-vs-measured comparisons agree to 1e-3 or better.
- The light-q check linearises in y, which holds at |y| ≤ 2.6e-3.
- θ* is checked as a ratio at fixed r_s, not with a CMB code.
- The 10 MeV derivative columns depend on the treatment of the table/ideal-gas junction.
