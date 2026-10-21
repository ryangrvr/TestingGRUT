# Module residual_bias: the systematic bias m0 = 2p − 1 that each restoration route leaves at the Z2 exit

21 October 2026. Claude module of the issue #13 wall-history packet (`BRIEF.md`). It answers the blocking review
finding that the report sorts hot histories into "unbiased network, rejected" against "wall-free" without having
computed the residual bias. The action is **unchanged**. Inflation, H_inf, N_*, reheating, any hidden-bath producer,
ξ = T_h/T_SM, the production time, the start condition of the observable patch and any hidden-charge supply are
**outside-action preparation data [P]**. Nothing here excludes or validates Candidate 1. Every statement is a
conditional history statement.

Labels: **[D]** derived here, **[C]** computed by the code here, **[Q]** quoted (sibling module, BRIEF,
literature), **[E]** estimate, **[P]** preparation datum outside the action.

The thresholds ε_persist ≈ 1e-20 and ε_safe ≈ 1e-17 are **placeholders** for the numbers that `modules/bias_threshold`
is deriving. This module reports m0 values with bands. "Below/between/above" refers only to those placeholders.
Sign convention: + is the previously selected vacuum, so m0 > 0 means the old vacuum is favoured.

## 0. Result in brief

| Route | Preparation [P] | m0 at the exit (central; band) | Status | Side of the placeholders |
|---|---|---|---|---|
| A, H★ | cold +v vacuum; μ_N = 0 bath of the four dof, ξ = 5e-4, present from T_SM = 1e8 GeV; exit T_SM = 3.1e6 GeV | log10\|m0\| = **−1297** (−650 to −2600) | [C] rate; [E] band | far below ε_persist |
| A, same history without plasma damping (Hubble friction only) | as above | **\|m0\| = 1** (coherent amplitude 34× the soft rms; one sign everywhere) | [C, E] | memory kept: no network |
| A, production_map WKB formula evaluated at H★ | as above | 1.7e-5 | [Q formula], **not applicable** here (§3.5) | — |
| A, any μ_N = 0 thermal bath, production at T_h ≥ 10 T_c | ξ free | \|m0\| < 1e-17 / 1e-20 / 1e-30 iff **ξ > 6.2e-5 / 6.7e-5 / 8.1e-5** (f = 10; 5.9e-5 / 6.4e-5 / 7.8e-5 at f = 32) | [C] | below ε_persist for every kinetically thermal bath (ξ ≥ 6.6e-5) produced well above T_c |
| A, bath produced just above T_c | T_h,prod/T_c = 1.1 (1.01) | below 1e-20 iff ξ > 1.4e-4 (3.0e-4) | [C] | — |
| A, free-streaming (non-thermal) restoring gas | exit T_SM,end, A = Π/H | ln-suppression **D = 4.68e-6 (A/GeV)**; below 1e-20 iff A ≳ 1.0e7 GeV, i.e. **T_SM,end ≲ 2.6e7 GeV** | [D] formula, [E] use | heavy-inflaton producer exits at 1e8–6e9 GeV: D = 3.3 to 1e-3, memory essentially kept |
| B, trapping only (no later restoration), fresh +v start, N_* = 50 | H_inf = 1.4 / 2 / 3 / 5 / 10 / 20 TeV | **+0.39 / +0.11 / −0.027 / −2.2e-3 / −3.6e-4 / −4.5e-5**; asymptotically **−0.36 (v/H)³** | [C] | far above ε_safe |
| B, same, relaxed + start | same | +0.43 / +0.10 / −0.095 / −4.5e-3 / −2.1e-3 / −5.2e-4; asymptotically **−0.21 (v/H)²** | [C] | far above ε_safe |
| B then thermal restoration (H★) | H_inf = 1e8 GeV, then the H★ bath | log10\|m0\| ≤ −1297 | [C, E] | far below ε_persist |
| B then a free-streaming restoring gas | H_inf = 1e5–1e8 GeV, exit 1e8–6e9 GeV | log10\|m0\| between −16.9 and −2.2 (§4.6) | [E] | above ε_safe (lowest edge −16.9 at H_inf = 1e8 GeV); an exit at 2.6e7 GeV would instead give < 1e-23 |
| C, Z2-symmetric thermal ensemble, no prior selection | — | **m0 = 0 exactly** | [D] | below anything |
| D, charge forced into the +v vacuum, n_s1 < n_i ≲ 2 n_c | physical slow expansion (Minkowski limit) | **+0.83 ± 0.07 (n_i = 1.5 n_c), −0.44 ± 0.11 (2 n_c)** at formation; minority then disappears | [C] 2D lattice | far above ε_safe |
| D, n_i ≥ 3 n_c | same | coherent memory decays at 0.09–0.13 m_X to the lattice floor (1–2%), then about 1e9 ℓ remain before the n_c exit | [C] to the floor; [E] beyond | below ε_persist if the decay stays exponential [E] |

Answers to the questions put to this module:

1. **The thermal-damping route (A) does erase the memory, by a huge margin, for every kinetically thermal bath
   produced well above T_c.** The k = 0 damping width of the soft φ mode, Γ/T_h = 3.84e-5 at T ≫ T_c (rising as
   1/m_φ near T_c), is 77 times the hard-particle elastic rate that `thermal_potential` lists. Integrated over the
   restored epoch of H★ it gives 2992 amplitude e-folds. The H★ classification "unbiased network" is therefore
   supported (conditional on its [P] inputs). Without plasma damping the same history would keep its sign: |m0| = 1.
2. **The inflation_spectator deviation p_− = 0.4997–0.4999 is not a statistical effect, and its quoted values are a
   grid artifact.** The trapping map is odd, F(−y) = −F(y), so a symmetric end distribution gives m0 = 0 exactly.
   The fresh +v start leaves an odd part of the end distribution, and its overlap with the oscillating trapping map
   gives a **systematic** bias. Converged values are p_− = 0.50109 (5 TeV, N_* = 50) and 0.50018 (10 TeV); with N_* = 60,
   0.50076 and 0.50015. They lie on the opposite side of ½ from the module's numbers, at distances from ½ between
   0.4 and 15 times the module's. The module's cell-centre sign
   assignment on h = 0.02 reproduces its own numbers (0.499935, 0.499720) and does not converge.
3. **Trapping alone leaves |m0| ≥ 3.6e-5 (fresh start) and ≥ 4.6e-4 (relaxed start) for 3 TeV ≤ H_inf ≤ 20 TeV, and
   0.08–0.82 for 0.7–2 TeV.** m0 changes sign once between 2 and 3 TeV. Only an H_inf tuned to that crossing to about
   1e-16 relative would come near the placeholders. Otherwise |m0| sits 12 or more decades above them. A "regime D"
   network made by inflation alone, with no later thermal restoration, is therefore biased. In the naive K10
   extrapolation (2D PRS fit, orientation only) it decays after a comoving coarsening of about 10^5.3–10^5.4,
   not 10^17.
4. **Forced charge (D) is two-sided.** Just above n_s1 the vacuum domains form during the decoherence of the coherent
   φ oscillation, and the bias is O(0.5) with a phase-set sign. From about 3 n_c the decision is postponed to the
   n_c exit, about 1e9 ℓ later, and the coherent remnant is erased [E].

## 1. Reproduce

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/residual_bias/code
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 python3 run_all.py
```

| Script | Wall time | Output |
|---|---|---|
| `a_damping.py` | about 75 s | `results/a_damping.json` |
| `b_inflation.py` | about 11 min | `results/b_inflation.json` |
| `d_charged.py validate` | 3 s | `results/d_validation.json` |
| `d_charged.py worker i 2` (two single-thread workers) | about 12 min for the 36 runs, 25–60 s each (0.5 s when the run files already exist and are reused) | `results/d_runs/*.json` |
| `e_synthesis.py` | 1 s | `results/headline.json` (per-route m0, classification, lattice summary) |
| `run_all.py` | — | `results/runtimes.json` (wall times, sha256 of every script) |

numpy/scipy only, deterministic seeds, every job under 15 min. `PYTHONDONTWRITEBYTECODE=1` keeps Python from writing
cache files. No other module's code is imported. The thermal_potential Kibble–Zurek table is read from its results
file, read only.

## 2. From a coherent remnant to m0 [D]

At a continuous exit the sign in each correlation volume is that of the coarse-grained soft field when the
instability takes over, near KZ freeze-out. Write that field as a homogeneous coherent remnant φ_coh plus Gaussian
soft fluctuations of rms σ on the correlation scale. Then

  p(+) = P(φ > 0) = ½[1 + erf(r/√2)],  **m0 = erf(r/√2) ≈ √(2/π) r**,  r = φ_coh/σ.

Choices made here:
- **σ.** Classical soft thermal fluctuation coarse-grained with a Gaussian window of radius 1/m̂:
  σ² = c_w T_c m̂, with c_w = [√π/2 − (π/2) e erfc(1)]/(2π²) = 0.01087. m̂ = 1/ξ_KZ is taken from thermal_potential's
  Ginzburg-matched z = 1 table [Q]. For H★: ξ_KZ = 2.03 GeV⁻¹ and σ = 2.9 GeV.
- **Coherent amplitude.** From production (φ_0 = v, a sudden bath on the cold +v vacuum, the conservative choice)
  it evolves by linear WKB, ∝ a^{−3/2} m_φ^{−1/2}, down to m̂, times exp(−D), where D = ∫(amplitude damping rate) dt.
  If φ_coh oscillates, |m0| also carries |cos(phase)| ≤ 1; the bound uses 1.
- **Hubble-only case.** Where λφ_coh² > m_φ² the amplitude is capped and then falls ∝ 1/a (quartic, conformal).

Every Z2-odd expectation value of the state must relax for memory to be lost. The slowest odd excitation is the soft
φ quasi-particle at k = 0. No conserved Z2-odd density exists (T^μν and j_X are even), so hydrodynamic long-time
tails do not feed ⟨φ⟩ [D/E argument].

## 3. Route A: thermal restoration of a cold +v vacuum

### 3.1 Width of a φ quantum at rest [D]

For a → b + … with a at rest (m_a = m_φ(T)), the 2↔2 processes a + b ↔ c + d give the net width (Weldon)

  Γ = Σ_ch mult |M|² S /(64π³ m_a) ∫dE_c dE_d F_stat Θ_kin,
  F_stat = f_b(1+f_c)(1+f_d) − (1+f_b)f_c f_d = f_b(1+f_c)(1+f_d)(1 − e^{−m_a/T}),

with E_b = E_c + E_d − m_a and **p**_b = **p**_c + **p**_d. Θ_kin = [|K| ≤ p_c p_d], where
K = E_cE_d − m_a(E_c+E_d) + (m_a² − m_b² + m_c² + m_d²)/2 is p_c p_d cosθ. Steps: the δ³ removes **p**_b; the δ in
cosθ gives E_b/(p_c p_d); and p dp/E = dE. The coherent amplitude decays at Γ/2. Where Γ > 2m_φ (overdamped) the slow
root m²/Γ is used.

Channels from U, with φ_0 the zero-mode quantum:

| Channel | Amplitude | S | Multiplicity |
|---|---|---|---|
| φ_0 φ → φ φ | 6λ_φ | ½ | 1 |
| φ_0 q_i → φ q_i | 2g_X² | 1 | 2 |
| φ_0 φ → q_i q_i | 2g_X² | ½ | 2 |
| φ_0 D → φ D | 2h_D² | 1 | 1 |
| φ_0 φ → D D | 2h_D² | ½ | 1 |

Thermal masses at φ = 0 (leading order, μ_N = 0): m_φ² = Π_φ(T² − T_c²) with T_c = 1.5575 v (inside the resummed
1.5571–1.5641 v [Q]); m_q² = Π_X T²; m_D² = m_D0² + Π_D T².

**1↔2 processes** [C]:
- At linear order in φ_coh there is no cubic vertex (exact Z2).
- With the background, φ → qq stays closed at every T (m_φ/2m_q ≤ 0.946).
- φ → DD is kinematically open for T ≳ 3 T_c. It enters only at second order in φ_coh, through the vertex
  h_D²φ_c φ D D. At H★ production it is 3.4e-6 of the linear width. It is not credited.

**Check against the literature** [C, Q]:
- One field with all masses m: the classical (f = T/E) integral equals **π²/(6m) to 4e-13**. That reproduces exactly
  Parwani's γ = g²m/(64π) for L ⊃ g²φ⁴/4!, m² = g²T²/24, pole Ω = ω − iγ, i.e. γ = g⁴T²/(1536π m) and Γ = 2γ.
  Source: hep-ph/9204216, eq. (3.29); read at section level in scratch, text not copied.
- For the pure φ⁴ part of U at its own thermal mass, the full-Bose width is 0.88 of that leading-order value.

### 3.2 The width [C]

| T_h/T_c − 1 | m_φ/T | Γ/T full Bose | classical | Boltzmann nσv (f_b only) | Γ/m_φ |
|---|---|---|---|---|---|
| 6.9e-4 (t_G) | 0.0024 | 6.4e-4 | 6.5e-4 | 8.7e-4 | 0.27 |
| 0.01 | 0.0091 | 1.86e-4 | 1.94e-4 | 2.49e-4 | 0.021 |
| 0.08 | 0.0246 | 8.3e-5 | 9.0e-5 | 9.3e-5 | 3.4e-3 |
| 0.93 | 0.0553 | 4.4e-5 | 5.0e-5 | 4.0e-5 | 8.0e-4 |
| ≫ 1 | 0.0646 | **3.84e-5** | 4.43e-5 | 3.36e-5 | 6.0e-4 |

- **Channel split at high T:** φφ→φφ 1.97e-5, φq→φq 1.26e-5, φφ→qq 6.1e-6, D channels ~1e-9.
- **The three estimates agree to 15% for ε ≳ 0.1, and to 35% near t_G.** These are the soft-region λ²T²/m form,
  the full Bose integral, and the Boltzmann collision rate of the k = 0 quantum, n⟨σv⟩ with s ≈ 2m_φE_b.
- **The soft-mode rate is not the hard-particle kinetic rate.** thermal_potential's Γ_φ,el = 4.98e-7 T is for
  collisions of hard particles, ⟨s⟩ = 2⟨E⟩². It is 77 times smaller and is not the rate at which the soft mode damps.
- The mode stays underdamped down to t_G.

### 3.3 H★ [C; P inputs]

| Quantity | Value |
|---|---|
| Restored-phase damping A | Π_φ ξ² M_P/√(π²g*/90) = 7.3e8 GeV |
| Exit | T_SM = 3.115e6 GeV, H = 1.36e-5 GeV |
| Γ/H near the exit (ε = 7e-3) | 2.5e4 |
| **D = ∫(Γ/2)dt from T_SM = 1e8 GeV to t_G** | **2992** (classical 3336, nσv 2959) |
| Dilute-gas formula of §3.6 at the same A | 3432 |
| Undamped linear amplitude at m̂ = 0.49 GeV | 446 GeV |
| r undamped | 155 |
| **log10 m0 (damped)** | **−1297** (−1447 classical, −1283 nσv) |
| Band, Γ × [½, 2] | −648 to −2597 [E] |

- The ×½ to ×2 band covers the unresolved O(ε_φ ≈ 0.1–0.3) next-order corrections and the scheme spread.
- **Hubble friction only:** the amplitude turns quartic at T_h = 1565 GeV and is 99 GeV at the decision.
  r = 34.4 and m0 = erf(24) = 1. The whole patch takes one sign, set by the oscillation phase.
- **Floor (not appropriate):** with thermal_potential's hard rate, D = 27.5 and log10 m0 = −9.9. If that rate had
  been used, H★ would have fallen "above ε_safe". This is why the soft-mode width had to be derived.

### 3.4 Scan in ξ and production temperature [C]

`a_damping.json → scan_xi_Tprod`: 17 values of ξ from 1e-6 to 1e-2, and T_h,prod/T_c = 1.001 to 1e3. D ∝ ξ² because
Γ/H ∝ a at fixed ξ. Selected values of log10|m0| (production at f = T_h,prod/T_c):

| ξ | f = 32.1 (H★-like) | f = 1.1 | f = 1.01 |
|---|---|---|---|
| 1e-5 | −0.004 | 0 | 0 |
| 3e-5 | −3.5 | 0 | 0 |
| 4.5e-5 | −9.2 | 0 | 0 |
| 6.6e-5 (ξ_kin) | −21.2 | −2.1 | 0 |
| 1e-4 | −50 | −8.2 | −0.07 |
| 2e-4 | −206 | −41 | −7.4 |
| 5e-4 | −1297 | −272 | −61 |
| 1e-3 | −5195 | −1096 | −254 |

Minimum ξ for log10|m0| below −17 / −20 / −30:

| Production | ξ_min |
|---|---|
| f = 32 (H★-like) | 5.9e-5 / 6.4e-5 / 7.8e-5 |
| f = 10 | 6.2e-5 / 6.7e-5 / 8.1e-5 |
| f = 2 | 7.7e-5 / 8.3e-5 / 9.9e-5 |
| f = 1.1 | 1.35e-4 / 1.44e-4 / 1.73e-4 |
| f = 1.01 | 2.8e-4 / 3.0e-4 / 3.6e-4 |

- The erasure threshold for production well above T_c, ξ ≈ 6.4e-5, coincides with the kinetic-equilibrium threshold
  ξ_kin = 6.6e-5 [Q]. Below ξ_kin the thermal rate is formal, because the bath is not kinetically thermal.
- The window 4.5e-5 ≲ ξ ≲ 6e-5, which would give m0 between ε_persist and ε_safe, therefore has no self-consistent
  thermal bath. Such baths belong to §3.6.

### 3.5 The production_map WKB estimate [Q formula]

- production_map's memory estimate is a coherent amplitude 1.7 v √λ v/A against a vacuum soft level 0.11 √λ v.
  At the H★ value A = 7.3e8 GeV it gives a ratio of 2.1e-5 (m0 ≈ 1.7e-5).
- That formula assumes the bath existed while H ≫ A, so that the zero mode stays frozen until H ~ A. In H★ the bath
  appears at H = 1.4e-2 GeV ≪ A. The zero mode then starts oscillating at production with amplitude about v and is
  far larger at the exit (§3.3).
- The ~1e-3 quoted in the review is not reproduced by any setting here. In any case the WKB number is superseded once
  plasma damping is included.

### 3.6 Free-streaming (non-thermal) restoring gas [D, E]

For a dilute gas with equal sharing among the four dof and no soft Bose enhancement, each collision removes a k = 0
quantum, and the inverse process is incoherent:

- Γ_d = C J/m_φ, with C = (18λ_φ² + 12g_X⁴ + 6h_D⁴)/(128π) = 7.46e-6 and J = Σ∫f/E.
- With Π = c_eff J, c_eff = 0.012525 [Q production_map], and Π = A H in the radiation era, the integral is
  ∫Π/m_φ dt = (A/(2√λ_φ v))·2 arctan√(u_0 − 1), u = AH/λ_φv² [D].
- Hence **D = 4.68e-6 (A/GeV)** for a zero mode that responds from early on.
- Applied to the thermal H★ value of A this gives 3432, against the full 2992.
- **Memory erased below 1e-20 iff A ≳ 1.0e7 GeV ⇔ T_SM,end ≲ 2.6e7 GeV (v/TeV)^{1/2}** [E]. This matches the
  thermal ξ_min (T_c/6.4e-5 = 2.4e7 GeV).
- At production_map's lower B_h edge (A ≈ 200 GeV), D = 9e-4. The heavy-inflaton producer (H★′; exits at
  T_SM = 1e8–6e9 GeV) gives D = 3.3 to 1e-3. **A cold +v start keeps its sign through such a restoration.**

### 3.7 Caveats [E]

- The leading-order width carries O(ε_φ) and O(m/T) corrections, covered by the ×½–×2 band.
- Inside the Ginzburg region (t < t_G) no damping is credited.
- Non-exponential components (threshold singularities of the two-loop spectral function) are smeared by the thermal
  widths of the internal lines. The slowest conceivable relaxation, at the hard-particle rate, still gives D = 27.5
  for H★.
- The amplitude at production is taken as v (sudden). An adiabatic switch-on gives a smaller amplitude, so m0 is
  smaller still.

## 4. Route B: randomized prior from inflation, then trapping

### 4.1 Identity [D]

The RD trapping ODE χ'' + λχ³ − λ(v/H)²τ²χ = 0, with φ frozen at the end of inflation, is odd under χ → −χ. So
F(−y) = −F(y), and

  **m0 = ∫₀^∞ [P_end(a) − P_end(−a)] s(a) da**, with s(a) = ±1 the sign-keeping map.

- A symmetric end distribution, such as equilibrium draws averaged over the ensemble, gives m0 = 0 exactly.
- A fresh +v start, or a relaxed + start, leaves an odd part of P_end. Its overlap with the rapidly alternating s(a)
  is a **systematic bias of the ensemble**. Both the FP evolution and the map are deterministic, so there is no
  statistical component in m0.
- Asymptotic structure [D]: s = +1 on a < a₁ = 2.953 v/H, then alternating with spacing ≈ v/(0.313 H). A smooth odd
  part P_odd ≈ κ₁ a gives m0 ∝ κ₁ (v/H)².
- For the fresh start κ₁ ∝ v/H, so m0 ∝ (v/H)³. For the relaxed start κ₁ = O(1), so m0 ∝ (v/H)².

### 4.2 Method [C]

- The 1D Starobinsky–Yokoyama FP is re-implemented here: X = D = 0, couplings and ξ_φ run to μ = H with the quoted
  one-loop betas. At 10 TeV it gives ξ_φ = −4.889e-4, a2 = −5.97e-3 and y_min = 0.770, the same as inflation_spectator
  [Q].
- Instead of propagating P forward, the **backward** equation m(y, n) = [e^{nQ}F](y) is solved, with F given as exact
  cell averages of the threshold map.
- The thresholds (up to 37 at 20 TeV) are located by bisection of the crossing parity, to 2^-40 of the scan step.
- Grid h = min(0.01, v/(10H)), checked at h × 2 and × ½. The fresh start is evaluated at y₀ = v/H exactly, by
  interpolation of the smooth odd m.

### 4.3 Results [C] (`b_inflation.json`)

| H_inf (GeV) | p_− end of inflation, fresh, N50 | m0 after trapping, fresh N50 (N60) | m0 relaxed N50 (N60) | rms of m over equilibrium draws | sd over one-e-fold sub-patches |
|---|---|---|---|---|---|
| 700 | 0.0904 | +0.819 (+0.778) | +0.684 (+0.651) | 0.73 | 0.042 |
| 1000 | 0.184 | +0.629 (+0.581) | +0.597 (+0.558) | 0.65 | 0.070 |
| 1400 | 0.264 | +0.391 (+0.346) | +0.433 (+0.399) | 0.47 | 0.068 |
| 2000 | 0.330 | +0.109 (+0.080) | +0.096 (+0.082) | 0.11 | 0.027 |
| 3000 | 0.385 | −0.0267 (−0.0297) | −0.095 (−0.090) | 0.11 | 0.013 |
| 5000 | 0.430 | **−2.18e-3 (−1.52e-3)** | −4.49e-3 (−3.90e-3) | 4.8e-3 | 1.7e-3 |
| 7000 | 0.450 | −1.07e-3 (−8.4e-4) | −4.08e-3 (−3.66e-3) | 4.4e-3 | 1.2e-3 |
| 10000 | 0.465 | **−3.63e-4 (−2.90e-4)** | −2.07e-3 (−1.86e-3) | 2.2e-3 | 5.7e-4 |
| 14000 | 0.475 | −1.32e-4 (−1.05e-4) | −1.05e-3 (−9.5e-4) | 1.1e-3 | 2.9e-4 |
| 20000 | 0.482 | **−4.51e-5 (−3.61e-5)** | −5.16e-4 (−4.66e-4) | 5.5e-4 | 1.4e-4 |

- **Asymptotic laws for H ≥ 7 TeV:**
  - fresh start: m0 = −0.36 (v/H)³ (N50) and −0.29 (v/H)³ (N60); the coefficient is flat to 1% from 7 to 20 TeV;
  - relaxed start: m0 = −0.206 (v/H)² (N50) and −0.186 (v/H)² (N60).
- **Grid convergence:** h = 0.005, 0.01 and 0.02 agree to ≤ 0.5%, e.g. −3.6455e-4, −3.6299e-4 and −3.6422e-4 at
  10 TeV (0.15% at 5 TeV, 0.06% at 3 TeV).
- **Cell-centre sign assignment** (the module's method) does not converge: −8.7e-4 (h = 0.005) against +5.6e-4
  (h = 0.02) at 10 TeV.
- **Emulating inflation_spectator** (h = 0.02, cell-centre signs) gives p_− = 0.499935 (5 TeV) and 0.499720 (10 TeV),
  against the module's 0.499930 and 0.499691 [Q].
- **Second artifact in the module:** its start cell sits at y = 0.11 instead of 0.10 at 10 TeV. That inflates the
  end-of-inflation odd part by about 10% (its p_− = 0.4612 against 0.4647 here).
- **Model limit:** like inflation_spectator, the map ignores preheating restoration, which that module flags as
  plausible for H_inf ≳ 15–20 v [Q]. The 14 and 20 TeV rows sit at that edge; above it, §4.5–4.6 apply.

### 4.4 The observable patch [E]

Our patch's bias is a realization: m(y*, N_*) for the field value y* that the patch had N_* e-folds before the end.

- Over equilibrium draws (long prior inflation) the ensemble mean is 0, but the rms is 2.2e-3 at 10 TeV and 5.5e-4
  at 20 TeV.
- For a fresh start, the spread between sub-patches one e-fold later (sd 5.7e-4 at 10 TeV) bounds the realization
  scatter of the patch average from above.
- **Either way the bias our patch sees has the same order as |m0|, and no realization comes near 1e-17.**

### 4.5 B followed by thermal restoration [C, E]

- The ensemble-mean coherent field entering a later restoration is ⟨φ_end⟩ = (0.86–1.25) v for a fresh start at
  H_inf = 1e4–1e13 GeV. Relative to the patch rms it is ⟨φ⟩/φ_rms = (0.85–1.04) v/H (9.3e-6 at 1e8 GeV).
- Through the H★ bath this gives log10|m0| ≤ −1297 + log10(⟨φ⟩/v) ≈ −1297.
- Dephasing of the large-amplitude quartic oscillations only lowers it.
- **H★ as specified (H_inf = 1e8 GeV plus the ξ = 5e-4 bath) stays in the persistent class.**

### 4.6 B followed by a free-streaming restoration (e.g. H★′, or the spectator gas alone) [E]

Weak damping, D = 4.68e-6 A. Band: lower edge 0.36 (v/H)³ e^{−D} (oscillating map), upper edge √(2/π)(⟨φ⟩/φ_rms) e^{−D}
(coherent mean, no dephasing).

| H_inf | T_SM,end = 2.6e7 | 1e8 | 1e9 | 6e9 GeV |
|---|---|---|---|---|
| 1e5 | −27.8 … −23.6 | −7.9 … −3.6 | −6.5 … −2.2 | −6.4 … −2.2 |
| 1e6 | −30.8 … −24.5 | −10.9 … −4.6 | −9.5 … −3.2 | −9.4 … −3.2 |
| 1e8 | −36.8 … −26.5 | −16.9 … −6.6 | −15.5 … −5.1 | −15.4 … −5.1 |

- Entries are log10|m0|.
- **For exits at T_SM,end ≥ 1e8 GeV both edges of the band lie above ε_safe.** At 2.6e7 GeV both lie below
  ε_persist, and in between the band straddles the placeholders.
- H★′ (heavy-inflaton producer, exits 1e8–6e9 GeV) is therefore **not** shown to be an unbiased-network history. Its
  classification needs this estimate replaced by a computation.

## 5. Route C: Z2-symmetric thermal ensemble [D]

- P commutes with H and with Q_X, and the state is P-invariant (thermal, grand canonical, or any P-invariant
  preparation). Then Tr(ρ O_odd) = 0, so p(+) = p(−) and **m0 = 0 identically** (also charge_phase §4 [Q]).
- The realized patch average fluctuates as white noise, ~N_dom^{−1/2}. In the observable patch at the H★ exit that is
  about 1e-33 [E] (about 1e52 Hubble volumes × 1e14 domains each). This is the same finite-size noise as in an
  unbiased lattice, not a systematic bias.

## 6. Route D: charge forced into the selected vacuum above n_s1 [C lattice, E extrapolation]

### 6.1 Setup

- **Code:** `d_charged.py`, an independent 2D classical lattice with the same equations and preparation as
  charged_exit_sim C2 [Q]. Fields φ and complex X, D = 0, 384², dx = 0.5 ℓ (L = 192 ℓ), dt = 0.1, leapfrog.
- **Validation:** energy drift 1.3e-3 over 60 ℓ in Minkowski; charge conserved to 2e-16.
- **Preparation:** S = +1 + noise (1e-2, s_c = 1 ℓ), X on its orbit in the +v vacuum at n_i = κ n_c (n_s1 = 1.352 n_c,
  n_* = 5.62 n_c).
- **Expansion rates:** τ_i = 150 (charged_exit_sim's rate), τ_i = 600 (4× slower H), and Minkowski with n fixed.
  Minkowski is the physical limit, since H/m_X ≈ 5.6e-10 makes the decoherence (about 50 ℓ) instantaneous on the
  Hubble time.
- **p(+) measure:** the mean of f+/(f+ + f−) over [τ_form + 5, τ_form + 25]. A single early output can catch a
  coherent swing.

### 6.2 Results (36 runs, `headline.json → D_summary`)

| Run set | κ = n_i/n_c (n_i/n_s1) | Seeds | windowed m0 (formation where a decision occurs) | Late state |
|---|---|---|---|---|
| expanding τ_i = 150 | 2 (1.48) | 3 | **+0.49 ± 0.12** (p = 0.70, 0.73, 0.81; module 0.71 [Q]) | p → 0.79–0.99 |
| expanding τ_i = 150 | 4 (2.96) | 3 | +0.05 ± 0.07 (module 0.47, 0.50 [Q]) | box-limited |
| expanding τ_i = 600 | 2 | 3 | **−0.40 ± 0.03** (the sign flips with the expansion rate) | p → 0 |
| expanding τ_i = 600 | 4 | 3 | +0.006 ± 0.022 | box-limited |
| Minkowski | 1.5 (1.11) | 4 | **+0.83 ± 0.07** over 3 seeds; the 4th never formed a − region (m0 = +1) | p = 1.00 |
| Minkowski | 2 (1.48) | 4 | **−0.44 ± 0.11** | p = 0.02–0.12 |
| Minkowski | 3 (2.22) | 4 | +0.01 ± 0.07 | f_core = 0.94–0.97; p 0.21–0.72 from a few late regions of random sign |
| Minkowski | 4, 5 (2.96, 3.70) | 4 + 4 | +0.007 ± 0.012, +0.012 ± 0.013 | stays restored, f_core 0.90–0.93; p = 0.500 ± 0.001 |
| Minkowski | 7 (5.18, above n_*) | 4 | +0.039 ± 0.003 (transient coherent remnant) | stays restored, f_core 0.88; p = 0.497–0.500 |

For Minkowski κ ≥ 3 no decision occurs within the run: the window then records only the decaying coherent remnant
and the symmetric fluctuating cells of the restored phase. The late columns are the meaningful ones there.

### 6.3 Decoherence picture [C, E]

- The forced state oscillates coherently at ω = 1.03–1.10 m_X, a period of about 6 ℓ.
- The coherent envelope stays near 1 for 5–7 oscillations (30–40 ℓ) and then collapses. The decoherence time is
  reproducible across seeds to about 1 ℓ (phase spread ≲ 1 rad), so the phase at decoherence is a deterministic
  function of κ, the seed level and the expansion. There is no phase averaging.
- **n_s1 < n_i ≲ 2 n_c:** vacuum regions form by phase separation during the collapse, while n > n_c. The sign is
  that of the oscillation at that moment: |m0| = 0.4–0.8, and the minority then disappears.
- **n_i ≳ 3 n_c:** the charged φ = 0 phase survives the collapse (f_core ≈ 0.9). The coherent remnant decays at
  **Γ_coh = 0.09–0.13 m_X** (amplitude; 16 Minkowski runs; 0.11–0.21 in 5 expanding runs) down to the 1–2%
  statistical floor of the box.
- The physical decision then waits for the expansion-driven n_c crossing, Δt ≈ (κ^{2/3} − 1)/(2H) ≈ 1.4e9 ℓ for
  κ = 4 at T_SM ≈ 200 TeV. Exponential continuation gives exp(−Γ_coh Δt) ≈ exp(−1e8).
- Cross-check: the dilute-gas formula of §3.6 for the decohered gas (ρ ≈ λv⁴/4, quanta of about 100 GeV) gives
  Γ_d/H ≈ 3e5 [E].
- Below ε_persist one needs only an average Γ ≥ 3.6e-8 ℓ⁻¹ over Δt. That is a very weak requirement, provided the
  decay stays exponential. A late power law slower than about t^{−2.4} (starting from the 1% floor at about 50 ℓ)
  would break it; none is expected for a non-conserved odd mode (§2).
- The boundary between the two regimes lies between 2 and 3 n_c in this 2D classical model with seed 1e-2. It is not
  resolved further.

## 7. What this changes in the draft report (conditional, for the synthesis)

- **H★ (thermal, ξ = 5e-4, production at 1e8 GeV): supported.** m0 ≈ 1e-1297 (band 1e-650 to 1e-2600), far below
  ε_persist. Statement 2's "a bath kinetically thermal above T_c erases it" now holds with a computed rate, for
  ξ ≳ 6.4e-5 at production well above T_c and ξ ≳ 1.4e-4 (3e-4) at production 10% (1%) above T_c.
- **"Inflation with H_inf ≫ v erases it (computed for the inflationary case)" is not supported as a p = ½
  statement.** Trapping alone leaves |m0| ≥ 3.6e-5 (fresh) or ≥ 4.6e-4 (relaxed) for 3–20 TeV, apart from the
  tuned sign crossing between 2 and 3 TeV. It falls only as (v/H)³ or (v/H)². A regime-D network from inflation
  alone, with no thermal hidden bath afterwards, is biased and decays early (naive K10 ≈ 10^5.3 at 20 TeV), leaving a
  minority-domain relic instead of a persistent network.
- **H★′ (free-streaming producer, exits 1e8–6e9 GeV)** keeps either the cold +v sign (|m0| ≈ 1) or the
  inflationary odd part (|m0| ~ 1e-17 to 1e-2 [E]). Its "same conclusion" is not established.
- **Forced charge:** the p(+) = 0.71 run is reproduced (0.70–0.81). In the physical (Minkowski) limit the bias is
  O(0.5) for n_i ≲ 2 n_c, and effectively 0 for n_i ≳ 3 n_c [E].

## 8. Free inputs [P]

- ξ, μ_N, the production temperature and how fast production happens (sudden vs adiabatic), the dof content of the
  bath.
- The cold +v start, H_inf, N_*, the start condition (fresh / relaxed / equilibrium draw), instant RD after
  inflation, and whether preheating restoration occurs (not modelled; H_inf ≳ 15–20 v).
- The hidden-charge supply (homogeneous, sudden, n_i/n_c) and the seed level of the lattice.
- The placeholder thresholds.

## 9. What failed or was corrected

1. **The first route-D "formation" statistic was wrong.** It used the first output with both signs ≥ 5%, and with
   0.5 ℓ output spacing it caught coherent swings: p = 0.72–0.77 at 4 n_c, τ_i = 150, against 0.50–0.56 decided.
   A second variant ("first output with f+ + f− ≥ 0.5") returned p = 0 or 1 for the same reason. Both were replaced
   by the windowed mean (§6.1). In the Minkowski κ ≥ 3 runs even the window records only the decaying coherent
   remnant, not a decision; the late values are the physical ones there.
2. `b_inflation` takes about 11 min in one job, inside the 15-min budget; most of it is the trapping-threshold
   bisection (up to 37 thresholds × 40 ODE solves at 20 TeV).
3. The production_map WKB amplitude could not be reproduced at the "~1e-3" level quoted in the task. Its formula gives
   2.1e-5 at the H★ value of A, and it does not apply to H★ (§3.5).
4. `common.py` briefly carried an unused stub (`m0_from_ratio`); it was removed before the production runs.

## 10. Not done

- No real-time quantum or Boltzmann simulation of the zero-mode decay through T_c. Critical dynamics inside t_G are not
  credited.
- No two-loop or next-to-leading width. The band ×½–×2 is an estimate.
- No μ_N < 0 kinetic bath with partial Bose enhancement; it is bracketed by §3.2 and §3.6.
- Route B uses the module's model: 1D FP with X = D = 0, instant RD, no preheating or fragmentation. The free-streaming
  propagation (§4.6) is an estimate, not a computation. The 2D FP correction to the odd part is not included.
- Route D: 2D only, classical, D = 0, one box size, seed level 1e-2, no hidden thermal bath. The boundary between 2
  and 3 n_c is not resolved. The exponential continuation of the remnant decay beyond the lattice floor is an
  estimate.
- The network consequence of each m0 (persistence, late decay, relic) belongs to `modules/bias_threshold` and
  `minority_relic`. The naive K10 numbers are orientation only.

## 11. Literature

| Reference | Read level | Used for |
|---|---|---|
| R. R. Parwani, *Resummation in a hot scalar field theory*, PRD 45 (1992) 4695, arXiv:hep-ph/9204216 | sections 2 and 3, eqs. (2.1)–(2.7) and (3.24)–(3.29) (PDF text-extracted in scratch; not copied) | zero-momentum damping γ = g²m/(64π), m² = g²T²/24, pole convention Ω = ω − iγ; reproduced exactly (§3.1) |
| S. Jeon, *Hydrodynamic transport coefficients in relativistic scalar field theory*, PRD 52 (1995) 3591, arXiv:hep-ph/9409250 | abstract | context: the effective Boltzmann description with thermal masses (no damping formula used) |
| H. A. Weldon, PRD 28 (1983) 2007 | not re-read (background knowledge) | Γ = Γ_d − Γ_i, the net width relation |
| Sibling modules thermal_potential, production_map, inflation_spectator, charged_exit_sim, charge_phase; REPORT.md; BRIEF.md | full or sections | quoted inputs [Q] |

## 12. Files

| File | Contents |
|---|---|
| `code/common.py` | constants, thermal masses, H(T), channels, m0 conversion |
| `code/a_damping.py` | route A: width integrals, checks, D integrals, decision ratio, scans, dilute gas |
| `code/b_inflation.py` | route B: FP backward equation, trapping thresholds, scaling, coherent mean |
| `code/d_charged.py` | route D: 2D lattice (`validate`, `worker i n`) |
| `code/e_synthesis.py` | per-route table, classification, lattice metrics → `results/headline.json` |
| `code/run_all.py` | runs everything → `results/runtimes.json` |
| `results/a_damping.json`, `b_inflation.json`, `d_validation.json`, `d_runs/*.json`, `headline.json`, `runtimes.json` | outputs |
