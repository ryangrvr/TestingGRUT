# Module thermal_potential: the zero-charge hidden-sector transition of unchanged Candidate 1

20 October 2026. Calculation module of the issue #13 wall-history packet (Claude side). This module does not say
that Candidate 1 is excluded or validated. It gives the phase transition of the **unchanged** action in one
specified ensemble: a hidden bath in internal equilibrium, with zero U(1)_X charge density and zero number
chemical potential. It also gives conditional inputs to the network history: the Kibble-Zurek domain scale and
the thermalisation range. Whether such a bath ever existed is preparation data outside the action (see "Free
inputs").

Reproduce (about 41 s on one core; every number below is in `results/*.json`):

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/thermal_potential
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 code/run_all.py
```

Labels: **[D]** derived analytically here, **[C]** computed by the code here, **[Q]** quoted from literature or
from the repository, **[E]** order-of-magnitude estimate.

## 0. Result in one paragraph

Suppose the hidden scalars once formed an internally thermal bath with zero charge, hotter than T_c. The Z2_phi
symmetry is then restored. It breaks again at **T_c/v = 1.5571–1.5641** [C]. This range covers the Arnold-Espinosa
(AE), Parwani and dimensional-reduction (DR) schemes, at mu = mu_0/2, mu_0, 2 mu_0 and v, with RG-evolved couplings.
The temperature is 0.61–1.07% above the leading curvature estimate 1.547646 [D]. The main upward shift comes from
the Debye-screened X loop. In every resummed scheme the perturbative transition is **second order**, with
phi_c/T_c = 0 and no latent heat [C]. The unresummed one-loop potential has a massless X at phi = 0. It shows a
spurious weak first-order signal, phi_c/T_c ≤ 0.0133 and L/rho_h ≤ 5.6e-7. That signal is not even stable under a
change of scale, and it lies inside the Ginzburg region [C]. Beyond perturbation theory, the soft phi mode is a 3D Z2
scalar with a positive effective quartic (λ_eff/λ_φ = 0.980–0.984 in DR) [C]. The expected transition is therefore
continuous, in the 3D Ising class. That is an argument, not a lattice computation.

X and D never condense, because their effective curvatures stay ≥ 2.79e-3 v² and ≥ 1.01e-2 v² at every T [C].
Phi is therefore the only order parameter. V_eff is even in phi to 7e-18, so p(+) = p(-) = 1/2 by symmetry of the
ensemble [D, C]. With p = 1/2 > p_c = 0.3116, both vacua percolate. In the uncorrelated-site model about 92% of the
initial wall area separates the two spanning clusters [C].

For T_h = ξ T_SM with the SM radiation dominated, Kibble-Zurek gives log10 N_H ≈ 3.2–25.8 (z = 1) and 3.6–30.0
(z = 2) domains per Hubble volume over ξ ∈ [1e-6, 1] [E]. The causal bound ξ_KZ ≤ 1/H is never saturated.
Internal kinetic equilibrium at T_c requires **ξ ≳ 6.6e-5** (phi elastic; 1.5e-4 for X) [E]. Number-changing
equilibrium (2→4) requires ξ ≳ 0.068, with a band of 0.008–0.27 [E]. D exchanges energy with phi only for
ξ ≳ 0.015 [E]. These are formation inputs only: persistence, scaling, expansion and backreaction are not computed
here.

## 1. Conventions and inputs

- Units: v = 1 (v = 1000 GeV). Metric (−,+,+,+). Fields: phi real, X = (q1 + i q2)/√2, D real.
- Potential (GRUT_CANDIDATE_1.md §2), in the normalisation of `grut_candidate_1_quantum_20261018/REPORT.md` §1,
  with R² = q1² + q2², c = g_X², d = h_D²:
  V0 = C + mp φ²/2 + mx R²/2 + md D²/2 + λφ φ⁴/4 + λX R⁴/4 + λD D⁴/4 + c φ²R²/2 + d φ²D²/2 + α R²D²/2.
- At mu_0 = g_X v = 0.1 v: (λφ, λX, c, d, λD, α) = (0.01, 0.001, 0.01, 1e-4, 0.001, 0), mp = −λφ v² = −0.01,
  mx = 0, md = m_D0² = 0.01, C = 0.0025 [Q, BRIEF]. Running uses the ten one-loop betas of
  `rg_beta_functions.txt` [Q]. Check: β_α(μ0) = 2cd/(16π²) = 1.2665e-8 [D].
- MSbar Coleman-Weinberg: V_CW = Σ_i n_i m_i⁴ [ln(|m_i²|/μ²) − 3/2]/(64π²) for scalars [Q, standard]. The SM does
  not enter the hidden V_eff, because every portal is zero at mu_0 and multiplicatively renormalised at one loop.
- M_P = 2.435323e18 GeV and g* = 106.75 [Q, BRIEF/task].

## 2. Field-dependent masses (4×4 Hessian) [D]

At a general background (φ, q1, 0, D), the Hessian H_ab = ∂_a∂_b V0 in the basis (φ, q1, q2, D) is

```
H_φφ = mp + 3λφ φ² + c q1² + d D²        H_φq1 = 2c φ q1     H_φD = 2d φ D
H_q1q1 = mx + c φ² + 3λX q1² + α D²      H_q2q2 = mx + c φ² + λX q1² + α D²
H_DD = md + d φ² + 3λD D² + α q1²        H_q1D = 2α q1 D
```

On the background (φ, 0, 0) the Hessian is diagonal. The entries are m_φ² = λφ(3φ² − v²), m_X² = g_X²φ² (+ m_X0² = 0
at μ0) twice, and m_D² = m_D0² + h_D²φ². The full matrix is used in `common.hessian_full` for item 2.

## 3. Debye masses [D]

At leading high T, the one-loop thermal potential is V_T ≈ −π²T⁴N/90 + (T²/24) Tr M²(q) + O(T m³). Hence
Π_ab = (T²/24) ∂_a∂_b Tr M². Summing the diagonal of the Hessian gives

Tr M² = φ²(3λφ + 2c + d) + R²(c + 4λX + α) + D²(d + 3λD + 2α) + const.

(The R² coefficient comes from c R² in H_φφ, (3+1)λX R² from the two q's, and α R² from H_DD. The D² coefficient
comes from d D² in H_φφ, 2α D² from the q's, and 3λD D² from H_DD.) Tr M² is quadratic, so Π is field-independent and
diagonal:

- Π_φ/T² = (3λφ + 2c + d)/12 = λφ/4 + g_X²/6 + h_D²/12 = **0.004175** (check against the brief: exact) [D, C]
- Π_X/T² per real mode = (c + 4λX + α)/12 = g_X²/12 + λX/3 = **0.00116667**, so m_X,D(T_c) = 0.03416 T_c [D, C]
- Π_D/T² = (d + 3λD + 2α)/12 = h_D²/12 + λD/4 = **0.000258333** [D, C]

The leading curvature zero is T/v = √(λφ/0.004175) = **1.547646465** [D]. It reproduces the brief.

## 4. Effective potential and schemes

V_eff(φ, T) = V0 + V_CW + (T⁴/2π²) Σ_i n_i J_B(m_i²/T²) [+ resummation], with n = (1, 2, 1) for (φ, X, D).
J_B(y²) = ∫₀^∞ dx x² ln(1 − e^{−√(x²+y²)}). For y² < 0 the real part is taken: ln|1 − e^{−iθ}| = ln(2|sin θ/2|),
θ = √(|y²| − x²). The same real-part prescription is used for (m²)^{3/2} and ln m².

| scheme | definition |
|---|---|
| oneloop | V0 + V_CW(m²) + V_T(m²), no resummation (X massless at φ = 0, so IR-unreliable) |
| AE | oneloop − (T/12π) Σ n_i [Re(m_i² + Π_i)^{3/2} − Re(m_i²)^{3/2}] for all three species |
| Parwani | V0 + V_CW(m² + Π) + V_T(m² + Π) |
| DR | AE for X and D. For φ, only the hard (non-zero Matsubara) part of its own loop: J_B + (π/6) Re(y²)^{3/2}, with no φ daisy. This is the tree potential of the soft φ zero mode in the 3D effective theory. Its soft loop is left to the nonperturbative 3D treatment. |

AE and Parwani are named for Arnold-Espinosa (hep-ph/9212235) and Parwani (hep-ph/9204216), following common usage.
Only the abstracts of those papers were read. The formulas above are the standard forms, not checked against the papers' text.

### 4.1 The exact J_B (item 1) [C]

- **Table:** direct `scipy.integrate.quad` of the defining integral on 5001 nodes, y² ∈ [−25, 25], step 0.01. The
  known non-analytic terms are subtracted, the analytic remainder R(y²) is cubic-spline interpolated, and the terms are
  added back. For y² > 25 the Bessel sum J_B = −Σ_k y² K₂(ky)/k² is used.
- **Convergent series** (production evaluator), derived here. Expand the n ≠ 0 Matsubara modes in m² using
  ∫d³p/(2π)³ (p² + a²)^{−k} = Γ(k − 3/2) a^{3−2k}/((4π)^{3/2} Γ(k)) and Σ_{n≥1} n^{3−2k} = ζ(2k − 3). This gives
  J_B = −π⁴/45 + π²y²/12 − (π/6)y³ − (y⁴/32) ln(y²/a_b) − 2π^{7/2} Σ_{l≥1} (−1)^l ζ(2l+1) Γ(l+½)/(l+2)! · (y²/4π²)^{l+2},
  with ln a_b = 2 ln 4π − 2γ_E + 3/2 = 5.407617. The series converges for |y| < 2π. In particular the y⁶ term is
  +ζ(3)y⁶/(384π²).
- **Checks** (`results/jb_checks.json`):
  - series vs quadrature: max 6.0e-13 on |y²| ≤ 10
  - Bessel vs quadrature: 2.8e-13 on 10 < y² ≤ 25
  - spline vs series at midpoints: 5.6e-13, derivative 9.1e-11
  - 60 random off-grid quadratures: 2.4e-13
  - negative y²: the analytic continuation of the series matches the real-part quadrature to 1.6e-13
  - Fit of J − (−π⁴/45 + π²y²/12) on y ∈ [0.02, 0.35] gives a cubic coefficient of −0.52359875 (expected −π/6 =
    −0.52359878), a log coefficient of −0.03124991 (expected −1/32) and a y⁴ coefficient of 0.1689883 (expected
    ln a_b/32 = 0.1689880).
  - The truncated high-T expansion through the log term has error 3.2e-7 at y² = 0.1 and 3.1e-4 at y² = 1. This
    equals the ζ(3)y⁶/(384π²) term. Dropping the cubic term gives error 1.7e-2 at y² = 0.1.
- T_c computed with the table and with the series agrees to ≤ 4.4e-11 in every scheme.

### 4.2 Zero-temperature check [C]

The numerical minimum of V0 + V1 at μ0 is φ/v = 1.0001218818; the archived strict one-loop value is 1.00012193740.
The curvature is 0.0200127832 v²; the archived value is 0.020012778428. The differences are O(loop²), as expected for
a numerical minimum against a strict expansion.

## 5. Transition temperature, order, φ_c, spinodals, latent heat (items 1 and 3) [C]

Method (`code/transition.py`):
- T0 is the zero of ∂²V/∂φ²|₀ = 2 dW/ds|₀, with W(s = φ²), computed analytically through J_B′.
- Coexisting minima are searched for T0 < T ≤ T0 + 0.02 v on a log grid φ ∈ [1e-7, 2] v (6000 points), and refined
  with brentq.
- If a coexisting minimum exists: T1 (upper spinodal) is found by bisection, T_c from ΔV = 0, and
  L = T_c ∂_T(V_b − V_s) by the envelope theorem.

| couplings | μ | oneloop T_c/v (order, φ_c/T_c) | AE T_c/v | Parwani T_c/v | DR T_c/v |
|---|---|---|---|---|---|
| running | μ0/2 | 1.549948 (first, 0.0064) | 1.564147 | 1.561916 | 1.559958 |
| running | μ0 | 1.549220 (first, 0.0105) | 1.563466 | 1.561520 | 1.559307 |
| running | 2μ0 | 1.548489 (first, 0.0133) | 1.562785 | 1.561123 | 1.558655 |
| running | v | 1.546962 (second, 0) | 1.561199 | 1.560202 | 1.557140 |
| fixed | μ0/2 | 1.549423 (first, 0.0105) | 1.563702 | 1.561470 | 1.559511 |
| fixed | μ0 | 1.549220 (first, 0.0105) | 1.563466 | 1.561520 | 1.559307 |
| fixed | 2μ0 | 1.549018 (first, 0.0106) | 1.563231 | 1.561570 | 1.559103 |
| fixed | v | 1.548546 (first, 0.0106) | 1.562682 | 1.561687 | 1.558628 |

**Adopted result:** T_c/v ∈ [1.557140, 1.564147] (1557.1–1564.1 GeV at v = 1 TeV). This is the range over AE,
Parwani and DR at four scales with running couplings; the fixed-coupling runs fall inside it. The transition is
perturbatively second order in all resummed schemes, so T0 = T1 = T_c, φ_c = 0 and L = 0. Below T_c the broken
minimum grows continuously as φ_b/T_c ≈ 1.10 √|ε| (AE) or 0.915 √|ε| (DR) for 1e-8 ≤ |ε| ≤ 1e-4 (mean-field
β = 1/2); at ε = −1e-2 the values are 0.103 (AE) and 0.092 (DR).

**Unresummed one-loop** (first-order rows, μ0):
- T0 = 1.5491167, T1 = 1.5492335, T_c = 1.5492205
- φ_c/T_c = 0.0105; barrier 1.9e-12 T_c⁴
- L = 2.67e6 GeV⁴, i.e. L/ρ_h = 3.5e-7 with ρ_h = (π²/30)·4·T_c⁴
- m_φ(broken, T_c) = 7.5e-4 T_c

Across scales φ_c/T_c spans 0.0064–0.0133 and L/ρ_h up to 5.6e-7.

**Why the resummed schemes give second order** [D, C]:
- The would-be cubic of the unresummed scheme is −(T/6π)(g_X²φ²)^{3/2}, strength E_X = g_X³/(6π) = 5.3e-5. That
  gives φ_c/T_c ≈ 2E_X/λφ ≈ 0.011, which matches the oneloop row. In AE and Parwani it becomes
  −(T/6π)(g_X²φ² + Π_X)^{3/2}, which is analytic for g_Xφ < m_X,D = 0.034 T, i.e. φ < 0.34 T ≫ 0.011 T.
- The φ self-daisy −(T/12π)(3λφφ² + M²)^{3/2} would give φ_c/T_c ≈ (3λφ)^{3/2}/(6πλφ) ≈ 0.028. At the curvature zero,
  however, M = √(mp + Π_φT²) = 0.0086–0.0092 T_c stays nonzero, because the screened X loop shifts the curvature
  by −g_X² m_X,D T/(2π). This screens the cubic for φ < M/√(3λφ) ≈ 0.05 T.
- The resulting quartic at φ = 0 stays positive: λ_eff/λφ = 0.579–0.597 (AE), 0.565–0.567 (Parwani), 0.980–0.984 (DR).

**Size of the beyond-leading shift.** The screened X loop shifts the curvature at φ = 0 by −g_X²√(Π_X/T²) T²/(2π) =
−5.4e-5 T², i.e. −1.3% of Π_φ, which raises T_c by about 0.65%. The rest is CW/thermal-log and, in AE/Parwani, the
soft φ daisy. The AE slope dm²(0,T)/dT² equals 0.854 Π_φ (Parwani 0.848; DR 0.987). The deficit in AE/Parwani is the
IR-sensitive soft-φ term −(3λφT/4π)M.

### 5.1 What is and is not controlled (item 3)

- **Expansion parameters at T_c:**
  - X: ε_X = max(g_X², 4λX)T/(4π m_X,D) = 0.0233 (controlled)
  - D: ε_D = max(h_D², 3λD)T/(4π m_D(T)) = 0.0036, with m_D(T_c) = 0.066 T_c (controlled)
  - soft φ at T_c (AE/Parwani): ε_φ = 3λφT/(4πM) = 0.258–0.278. This goes to infinity as T → T_c in the exact
    theory (m_φ → 0).
  - In DR, ε_φ(ε) = 3λφT/(4π m_φ) rises from 0.089 (|ε| = 0.1) to 0.83 (1e-3) and 8.3 (1e-5).
- **Ginzburg region:** ε_φ = 1 with mean-field m² = 2κT_c²|t|, κ = 0.0041196 (DR), gives
  **t_G = (3λφ/4π)²/(2κ) = 6.9e-4**, with a band of 1.9e-5–2.8e-3 for O(1) criteria (λ_3/(8πm) = 1 up to
  ε_φ = 1/2). The unresummed first-order window (T1 − T0)/T_c is 4e-6–7.5e-5, entirely inside t_G, with
  ε_φ(φ_c) = 2.2–8.3. The signal is not under control.
- **Nonperturbative order** (argument, not computed): DR integrates out the hard modes and the soft X and D zero modes
  (controlled at ε_X, ε_D ≲ 0.02). It leaves a 3D real scalar with exact Z2 and positive quartic λ_3 = λ_eff T
  (λ_eff/λφ = 0.98). Such a theory has a continuous transition in the 3D Ising class. Becoming first order would need
  λ_eff ≤ 0. The X/D-induced reduction is about 2% (analytic estimate −g_X⁴T/(16π m_X,D) = −5.8e-5 vs λφ/4 = 2.5e-3).
- **T_c is controlled to roughly ±0.3%.**
  - The scheme spread (0.45%) and the AE−DR difference (0.27%, the soft φ loop) set the scale.
  - A two-loop soft-sunset estimate δm² ≈ 6λφ²T²L/(16π²), with L ∈ [1, ln(T_c/m_G)], gives δT_c/T_c = 0.05%–0.28% [E].
  - No lattice value of the 3D critical mass is used.
- **Not controlled:** critical amplitudes, the exact T_c beyond about 0.3%, dynamics near T_c.

## 6. X and D never condense (item 2) [C]

The curvatures ∂²V/∂q1² and ∂²V/∂D² are evaluated by finite differences (h = 1e-3 v) of the full four-field V_eff,
using the eigenvalues of the 4×4 Hessian with exact J_B. They are taken at the scheme's own φ_min(T), for 31
temperatures in [0, 3]v plus T_c(1 ± 1e-3):
- min m²_X,eff: 2.794e-3 (oneloop), 2.802e-3 (AE), 2.802e-3 (Parwani), in units of v², reached at T_c. This is
  0.998 (oneloop), 0.982 (AE) and 0.985 (Parwani) of the leading value g_X²φ_min² + Π_X T²; the 1.5–1.8% deficit is
  the O(T m) loop correction;
- min m²_D,eff: 1.0100e-2 v², reached at T = 0.

A global grid over (φ, R, D) ∈ [0, 1.6]³ (49 × 49 × 33 points) at T = 0, 0.5, 1.0, 1.4, T_c, 1.7, 2.5 v (AE and
Parwani) puts the minimum at R = D = 0 every time. The lowest off-axis grid point lies 1.6e-6–5.6e-6 v⁴ above it.
Nelder-Mead from five off-axis starts converges to |R|, |D| ≤ 3.0e-7 (refinement tolerance). Analytically, with zero
charge density the X and D directions have tree curvature g_X²φ² ≥ 0 and m_D0² + h_D²φ² > 0, plus positive Debye
masses. Phi is therefore the only order parameter of this ensemble.

## 7. Symmetry: p(+) = p(-) = 1/2 (item 3) [D, C]

V0 depends on φ only through φ², except for the off-diagonal Hessian entries 2g_X²φq1 and 2h_D²φD. Flipping φ is an
orthogonal similarity diag(−1, 1, 1, 1) of the Hessian, so every eigenvalue, and hence every scheme's V_eff, is even.
Numerically, |V(φ, q1, D) − V(−φ, q1, D)| ≤ 6.9e-18 v⁴ at 2000 random points × 5 temperatures in all schemes; the
background V(φ) − V(−φ) is exactly zero. For a Z2-invariant state entering the transition, for example any thermal
state above T_c (ρ ∝ e^{−H/T} with [H, Z2] = 0), every Z2-odd expectation vanishes. Each correlation volume
therefore selects ±v with p = 1/2 exactly, whatever the order of the transition. This is a statement about the
ensemble. It fails only if the incoming state is not Z2 invariant, for example a hidden sector that never restored
the symmetry, or an asymmetric non-thermal preparation. Both are outside-action preparation data.

**Supplement** (`code/percolation.py`) [C]:
- Model: uncorrelated ±1 sites on a simple-cubic lattice (spacing = freeze-out length), p = 1/2, open boundaries.
- At L = 32, 64 and 128 (3 seeds each), the largest + and − clusters both span all three directions. They hold
  ≥ 96.7%, ≥ 97.3% and ≥ 97.6% of their sign respectively.
- The fraction of wall area separating the two spanning clusters is 0.902, 0.910 and 0.916 (rising with L). The rest
  bounds finite closed bubbles.
- This is consistent with p = 1/2 > p_c(site, sc) = 0.3116077 [Q].
- It is initial geometry only, not persistence.

## 8. Kibble-Zurek input (item 4) [E]

**Assumptions:**
- T_h = ξ T_SM with both entropies separately conserved and constant g, so T_h ∝ 1/a and dT_h/dt = −H T_h.
- The SM is radiation dominated, H(T_SM) = √(π²g*/90) T_SM²/M_P, g* = 106.75. Including the hidden energy changes H
  by a factor √(1 + 4ξ⁴/106.75) ≤ 1.0186.
- Then τ_Q = T_c/|dT_h/dt| = 1/H(T_c/ξ).

**Microscopic scales** (DR, μ0): m_φ²(0, T) = 2κT_c²ε with κ = 0.0041196, so ξ₀ = 1/(T_c√(2κ)) = 11.02/T_c.
Convention τ₀ = ξ₀ (c = 1); sensitivity to τ₀ × 10^{±2} is tabulated.

**Freeze-out:** τ(ε̂) = ε̂ τ_Q with τ(ε) = τ₀(ξ(ε)/ξ₀)^z. ξ(ε) = ξ₀ε^{−1/2} for ε ≥ t_G, and ξ₀ t_G^{−1/2}(ε/t_G)^{−ν}
below ("matched"); or ξ₀ε^{−ν} throughout ("naive").
- ν = 1/(3 − Δ_ε) = 0.629971, with Δ_ε = 1.412625 [Q].
- z = 1: underdamped, relativistic order-parameter dynamics. z = 2: relaxational, model-A-like dynamics.
  Model C (z ≈ 2.17) is not computed.
- Naive closed form: N_H = [τ_Q^{1−ν/(1+νz)} τ₀^{ν/(1+νz)}/ξ₀]³ with τ_Q = ξ² M_P/(√(π²g*/90) T_c²). This gives
  log-slopes d log10 N_H / d log10 ξ = 6(1 − ν/(1+νz)) = 3.68 (z = 1) and 4.33 (z = 2).

**m_φ(T) near T_c (DR):** m_φ/T_c = 0.0294, 0.0091, 0.00287, 9.08e-4, 2.87e-4 and 9.08e-5 at ε = 1e-1, 1e-2, 1e-3,
1e-4, 1e-5 and 1e-6. It matches the mean-field √(2κε) within 2.5% at ε = 0.1 and within 1e-6 at 1e-6. The broken side
is √2 larger. These are mean-field masses. Inside t_G the true correlation length follows the Ising exponent.

| ξ = T_h/T_SM | T_SM (GeV) | H (GeV) | Γ_φ,el/H | log10 N_H z=1 (matched / naive) | log10 N_H z=2 (matched / naive) | Kibble-Ginzburg | ξ_KZ (z=1, matched) |
|---|---|---|---|---|---|---|---|
| 1e-6 | 1.56e9 | 3.42 | 2.3e-4 | 3.23 / 2.98 | 3.64 / 3.50 | 0.11 | 0.024 GeV⁻¹ |
| 1e-5 | 1.56e8 | 3.4e-2 | 0.023 | 7.23 / 6.66 | 8.14 / 7.83 | 6.11 | 0.113 GeV⁻¹ |
| 1e-4 | 1.56e7 | 3.4e-4 | 2.27 | 11.09 / 10.34 | 12.64 / 12.15 | 12.11 | 0.586 GeV⁻¹ |
| 1e-3 | 1.56e6 | 3.4e-6 | 227 | 14.78 / 14.02 | 17.03 / 16.48 | 18.11 | 3.48 GeV⁻¹ |
| 1e-2 | 1.56e5 | 3.4e-8 | 2.3e4 | 18.46 / 17.70 | 21.35 / 20.81 | 24.11 | 20.6 GeV⁻¹ |
| 1e-1 | 1.56e4 | 3.4e-10 | 2.3e6 | 22.14 / 21.38 | 25.68 / 25.14 | 30.11 | 122 GeV⁻¹ |
| 1 | 1559 | 3.4e-12 | 2.3e8 | 25.82 / 25.06 | 30.01 / 29.46 | 36.11 | 725 GeV⁻¹ |

Further notes:
- With τ₀ × 10^{±2}, log10 N_H spans 23.5–28.1 (z = 1) and 28.3–31.7 (z = 2) at ξ = 1. At ξ = 1e-4 it spans
  8.8–13.2 and 11.0–14.1.
- The causal bound ξ_KZ ≤ H^{−1} (the radiation-era particle horizon 2t = 1/H) is never saturated: N_H ≥ 10³ on the
  whole grid (≥ 17 with τ₀ × 1e-2 at ξ = 1e-6).
- Rows with Γ_φ,el/H < 1 (ξ ≲ 6.6e-5) have no thermal bath at T_c. There the KZ numbers are formal only.
- The freeze-out ε̂ (z = 1, matched) is 2.5e-9 at ξ = 1, deep inside t_G. It lies on the mean-field side
  (ε̂ > t_G) only for ξ ≤ 3.2e-5 (z = 1) or ξ ≤ 2.0e-4 (z = 2) on the grid.

## 9. Hidden self-thermalisation (item 5) [E]

**Rate formula:** Γ_a = Σ_b n_b Σ_f |M|² S_f/(16π⟨s⟩), with
- n_b = ζ(3)T³/π² per real dof;
- ⟨s⟩ = 2⟨E⟩² = 14.59 T², with ⟨E⟩ = π⁴T/(30ζ(3)) = 2.701 T;
- S_f = 1/2 for identical final particles.

Contact amplitudes from U:
- φφ→φφ: 6λφ
- φφ→q_iq_i and φq_i→φq_i: 2g_X²
- q1q1→q1q1: 6λX; q1q1→q2q2 and q1q2→q1q2: 2λX
- φφ→DD and φD→φD: 2h_D²
- DD→DD: 6λD
- α = 0 at μ0

Summed rates:

| rate | coupling sum | Γ/T_h |
|---|---|---|
| Γ_φ | 18λφ² + 12g_X⁴ + 6h_D⁴ | 4.98e-7 |
| Γ_X | 24λX² + 6g_X⁴ | 1.04e-7 |
| Γ_D,exchange with φ | 6h_D⁴ | 1.0e-11 |
| Γ_D,self | 18λD² | 3.0e-9 |

**Number changing.** 2→4 via σ = |M₄|²Φ₄/(2sS_f), with
- Φ₄ = s²/(24576π⁵), the massless four-body phase space;
- |M₄| = N_d(6λφ)²κ_p/s, N_d = 10 tree diagrams, S_f = 4!;
- central κ_p = 4, N_ch = 1; band κ_p ∈ [1, 16], N_ch ∈ [1, 4].

The 2→4 rate is Γ/T = 4.8e-13 (band 3.0e-14–3.1e-11). The 1↔3 channels are kinematically closed, and no cubic
vertex exists in the symmetric phase (exact Z2_φ).

**Criterion:** Γ(T_h = T_c) > H(T_SM = T_c/ξ). Γ/H ∝ ξ²/T_h for contact interactions, so a bath that is thermal at
T_c stays thermal below it. The resulting thresholds:

| channel | ξ_min |
|---|---|
| φ kinetic | 6.6e-5 |
| X kinetic | 1.45e-4 |
| D self | 8.6e-4 |
| D–φ energy exchange | 0.0148 |
| 2→4 chemical | 0.068 (band 0.0084–0.27) |

With v free, every ξ_min scales as (T_c/1559 GeV)^{1/2} ∝ (v/TeV)^{1/2} at fixed dimensionless couplings.

**Consequences:**
- For ξ below the chemical threshold, a kinetically thermal bath can carry a number chemical potential μ_N ≤ 0. The
  leading Debye coefficient scales as Li₂(e^{μ_N/T})/ζ(2), so T_c, quoted in kinetic temperature, scales as
  √(ζ(2)/Li₂(z)). For example z = 0.9 gives ×1.125, z = 0.5 gives ×1.681 and z = 0.1 gives ×4.00 [D]. The
  "zero chemical potential" premise is therefore an input unless ξ ≳ 0.07.
- If D is not in contact with φ (ξ ≲ 0.015), dropping D's h_D²/12 from Π_φ raises the leading T_c by a factor of
  √(0.004175/0.0041667) = 1.0010 [D].

## 10. Free inputs on which the results depend

1. The hidden sector was internally thermal with T_h,max > T_c, and the state entering the transition is
   Z2-invariant (outside-action preparation).
2. Zero U(1)_X charge density (charge chemical potential μ_Q = 0); finite charge is the separate `charge_phase` module.
3. Zero number chemical potential, or chemical equilibrium (ξ ≳ 0.07).
4. ξ = T_h/T_SM and a visible reheating temperature above T_c/ξ (outside-action preparation data; inflaton and
   reheating are not part of the action).
5. The scale v. T_c/v and every dimensionless output are v-independent; τ_Q T_c, N_H and ξ_min depend on v/M_P.
6. Matching at μ0 (the brief's zeros) and the choice of resummation scheme and scale, which enter as the quoted spread.
7. KZ dynamics: z, the τ₀ = ξ₀ convention, and naive vs Ginzburg-matched exponents.
8. g* = 106.75 (SM above the EW scale; the ν_R are assumed not thermalised).

## 11. Conditional history constraints (not predictions)

- **C1.** If (1)–(3) hold, unchanged Candidate 1 restores Z2 and re-breaks it continuously at T_h = (1.557–1.564)v
  (±0.3% unresolved soft-loop systematic). Domains of both vacua form with p = 1/2 and both percolate. This history
  supplies a formation input, not a demonstrated persistent network.
- **C2.** If the bath is thermal at T_c (ξ ≳ 7e-5 (v/TeV)^{1/2}), the initial domain density is 10^{8.8}–10^{31.7}
  per Hubble volume (τ₀ band) for ξ ∈ [1e-4, 1]. That is far above one, so the initial configuration is
  sub-horizon and wall-crowded. Survival, scaling and domination must come from the network modules, with
  expansion and backreaction.
- **C3.** Avoiding wall formation by this route requires one of: (a) T_h,max < T_c, (b) the hidden sector never
  reached internal equilibrium above T_c, for example ξ ≲ 7e-5 (v/TeV)^{1/2} with a non-thermal state that is not
  Z2-symmetric, or (c) a nonzero charge density that changes the phase structure. Each is a statement about
  preparation, not about the action.

## 12. What failed or was corrected

1. The first J_B series used the denominator (l+1)!, a form recalled from memory. It failed the quadrature check by
   6.3e-7 at y² = 0.1 (wrong y⁶ coefficient ζ(3)/(128π²)). It was re-derived with (l+2)!, which agrees with both
   quadrature and the Bessel sum to 1e-13. The failing variant is recorded in `jb_checks.json`
   (`series_denominator_note`). I have not checked which published sources carry which form.
2. **Robustness check FAILS:** the unresummed one-loop order is not scale-stable. It is first order at seven settings
   and second order at running μ = v. The cause is that the RG-generated tree m_X²(μ) = −2.9e-6 v² (about (1.7 GeV)²)
   removes the unscreened X cubic. The unresummed result is IR-unreliable and is not adopted. 16 of 17 automated
   checks pass; this is the failing one.
3. The first `transition.py` run exceeded the time budget (more than 10 minutes, killed) because the series was
   evaluated with explicit powers. Horner evaluation cut the total to about 41 s. No numbers changed.
4. The AE and Parwani curvature slope is 15% below the leading Π_φ. This is a symptom of the soft φ IR sensitivity,
   not a bug: DR, which removes the soft φ loop, gives 0.987.
5. The thermalisation and 2→4 numbers are order-of-magnitude. Collinear logarithms, Bose enhancement and thermal
   masses are not resolved; the 2→4 band spans a factor of about 1000 in rate.

## 13. Not done

- No two-loop V_eff, no 3D lattice T_c or critical amplitudes. The order of the transition rests on the
  universality argument.
- No finite charge chemical potential, which is the `charge_phase` module's job.
- No real-time simulation of the quench; z and τ₀ are assumed.
- No production mechanism for the hidden bath and no T_h/T_SM prediction. The gravitational (graviton-exchange)
  production that the action itself contains was not computed.
- No wall-network evolution, scaling, domination, backreaction, CMB or expansion observable.
- No scan in couplings; only v-scaling is stated.
- No model-C or model-H dynamic classes.

## 14. Literature and sources

| id | title | read level | used for |
|---|---|---|---|
| arXiv:hep-ph/9212235 (PRD 47, 3546) | Arnold, Espinosa, The effective potential and first-order phase transitions: beyond leading order | abstract | name/attribution of the AE daisy scheme; the formula is the standard one, not checked in the text |
| arXiv:hep-ph/9204216 (PRD 45, 4695) | Parwani, Resummation in a hot scalar field theory | abstract | name of the Parwani scheme (same caveat); context for scalar thermal mass and damping at O(g³) |
| arXiv:1612.00466 | Curtin, Meade, Ramani, Thermal resummation and phase transitions | abstract | context: AE/Parwani-type resummations carry scheme errors, hence the spread is reported |
| arXiv:hep-ph/9901312 | Quiros, Finite temperature field theory and phase transitions | abstract | context for J_B, high-T expansion and IR breakdown; series re-derived here, not compared with the text |
| arXiv:hep-ph/9508379 (NPB 458, 90) | Kajantie, Laine, Rummukainen, Shaposhnikov, Generic rules for high-temperature dimensional reduction | abstract | logic of the DR split (hard and soft modes, 3D effective theory) |
| arXiv:1603.04436 (JHEP 08 (2016) 036) | Kos, Poland, Simmons-Duffin, Vichi, Precision islands in the Ising and O(N) models | abstract | Δ_ε = 1.412625, giving ν = 0.629971 |
| arXiv:1302.0421 (PRE 87, 052107) | Wang, Zhou, Zhang, Garoni, Deng, Bond and site percolation in three dimensions | abstract | p_c(site, sc) = 0.3116077 |
| doi:10.1088/0305-4470/9/8/029 | Kibble 1976, J. Phys. A 9, 1387 | not re-read (background knowledge) | Kibble/Ginzburg domain estimate |
| doi:10.1038/317505a0 | Zurek 1985, Nature 317, 505 | not re-read | KZ freeze-out logic |
| doi:10.1103/RevModPhys.49.435 | Hohenberg, Halperin 1977 | not re-read | model A z ≈ 2 |
| doi:10.1103/PhysRevD.9.3320 | Dolan, Jackiw 1974 | not re-read | leading thermal masses (re-derived in §3) |
| repo: grut_candidate_1_20261018/GRUT_CANDIDATE_1.md §1–3 | action and matching | full (§1–5) | potential, symmetries, matching |
| repo: grut_candidate_1_quantum_20261018/REPORT.md, rg_beta_functions.txt | one-loop potential, betas | sections 1–6; full | normalisation, betas, zero-T check values |
| repo: grut_issue13_cross_review_20261020/REVIEW.md, code/review.py (R6) | cross-review | full; section | leading-estimate provenance |
| git 135d535: grut_issue13_response_20261020/README.md, DERIVATION.md, FEASIBILITY.md; CLAIM_EVIDENCE.md | Codex packet | sections; grep only for CLAIM_EVIDENCE | requirements and claim row "Hidden transition occurs at 1.55 TeV: estimate only" |

Not read: SOURCES.md and RESULTS.json at 135d535; the 479f111 3D solver packet (not needed for this module).

## 15. Files

- `code/common.py`: couplings, RG, J_B (quadrature, table, series, Bessel), Debye masses, V_eff schemes, 4-field V_eff.
- `code/jb_checks.py` → `results/jb_checks.json` (the transient `jb_table.npz` cache is deleted by `run_all.py`).
- `code/transition.py` → `results/transition.json`.
- `code/xd_positivity.py` → `results/xd_positivity.json`.
- `code/kz_thermalization.py` → `results/kibble_zurek.json`, `results/thermalization.json`.
- `code/percolation.py` → `results/percolation.json`.
- `code/run_all.py` → `results/summary.json`: headline numbers, 17 assertions (16 pass, 1 robustness FAIL), and
  sha256 hashes of the scripts.
