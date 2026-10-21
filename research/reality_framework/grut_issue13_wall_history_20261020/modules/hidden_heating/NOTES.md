# Module hidden_heating: the sector coupling, irreducible hidden heating and the outside-action reheating chronology

21 October 2026. Calculation module of the issue #13 wall-history packet (Claude side). This module does not
say that Candidate 1 is excluded or validated. It gives (i) what couples the hidden scalars to the visible
sector inside the **unchanged** action, (ii) the hidden energy and vacuum randomisation that gravity alone
produces in a stated history, and (iii) **conditional history constraints** on preparation data that lie outside
the action: the inflationary Hubble rate, the reheating chronology and a direct inflaton branching.

Reproduce (about 22 s; deterministic; MC seed 20261020; every number below is in `results/*.json`):

```
cd research/reality_framework/grut_issue13_wall_history_20261020/modules/hidden_heating
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 python3 code/run_all.py
```

Labels: **[D]** derived analytically here; **[C]** computed by the code here; **[Q]** quoted (literature or
repository); **[E]** order-of-magnitude estimate; **[P]** outside-action preparation data (inflaton, inflation,
reheating). Headline numbers are in `results/headline.json`.

## 0. Result in one paragraph

Inside the action, the two sectors couple **only through the metric** [D, C]. The three Higgs portals are zero at μ0,
and they stay zero under the non-gravitational RG to all orders, because the theory factorises. A one-loop tensor check
gives zero mixed betas (residue ≤ 7e-20 roundoff). Graviton-induced local portals are ≲ 3e-34 [E]. The running
curvature coupling ξ_φ(μ) is negative for μ > μ0, so ξ R φ² deepens the broken minimum and does not restore it [C].
What remains is graviton exchange and the shared background. I computed the s-channel freeze-in rate from the SM bath
into the four hidden real scalars: dρ_h/dt = **β_E T⁹/M_P⁴, β_E = 0.06815** [C]. Its angular shapes agree exactly with
Clery et al. 2112.15214, and its channel ratios match the Weyl-anomaly central charges (1 : 3 : 12). Clery's total
β₀ is reproduced only if their 45-Weyl channel is doubled (ratio of ours to theirs 0.7242).

If the hidden gas self-thermalised, T_h,max < T_c = 1.5476 v would need **T_RH < 1.35e10 GeV** for instantaneous
reheating (1.10–1.57e10 GeV for T_c ± 30%, ∝ v^{4/7}). For inflaton-dominated reheating the condition is
**T_max < 6.0e9 (H_e/GeV)^{1/9} (v/TeV)^{4/9} GeV**, which is 1.3e10–1.65e11 GeV for H_e = 1e3–1e13 GeV [C]. The gas
does **not** self-thermalise, however: Γ₂₂/H ≈ 2e-26 [E]. The physical, nonthermal criterion
Π_φ > λ_φ v² needs T_RH ≳ 9.1e12 GeV, or 8.3e14 GeV with Hubble damping [C].

The strongest irreducible effect is gravitational, and it comes from inflation [C, P]. In de Sitter space the light
X fluctuations alone restore Z2_φ in the stochastic equilibrium for **H_inf > 436 GeV = 0.436 v**. The origin first
becomes locally stable at 0.379 v. With a finite duration and X starting at zero, the threshold is 0.57 v after 60
e-folds. Below that the φ distribution is bimodal, and wrong-vacuum regions of macroscopic size are avoided only for
**H_inf ≲ 184–275 GeV** (N_obs = 60; domains within 20–5 e-folds of the observable scale; heuristic count).

The inflaton condensate scatters gravitationally into the hidden quartet with **the same energy as into the Higgs
quartet**. Its formal thermal threshold is H_e ≈ 3.1–5.6 TeV, but that gas does not thermalise either. Its
nonthermal Π threshold is H_e ≳ 4e5–9e8 GeV. For a direct branching, T_h,max = (B_h g*/g_h)^{1/4} T_max gives
**B_h < (g_h/g*)(T_c/T_max)⁴ = 2.15e-29 (1e10 GeV/T_max)⁴**. That form is self-consistent only for T_max ≲ 2e7 GeV,
where the hidden bath can be kinetically thermal. Gravitational production is C-symmetric, so every gravity-only
channel gives zero hidden U(1)_X charge and μ = 0.

These are formation and heating conditions only. Persistence of a network is not addressed here.

## 1. Conventions and inputs

- Action and benchmark: `GRUT_CANDIDATE_1.md` §2–3, BRIEF.md [Q]. Metric (−,+,+,+) in the action. The amplitude code
  uses (+,−,−,−) internally for the Dirac algebra; all outputs are metric independent. M_P = 2.435323e18 GeV.
  v = 1000 GeV, g_X = 0.1, λ_φ = 0.01, λ_X = 0.001, h_D = 0.01, λ_D = 0.001, m_D0 = 100 GeV, μ0 = g_X v = 100 GeV.
- Real basis (φ, q1, q2, D), X = (q1 + i q2)/√2, R_X² = q1² + q2². Quartic normalisation follows the quantum packet:
  c = g_X², d = h_D².
- Hidden induced curvature at φ = 0: Π_φ = 3λ_φ⟨δφ²⟩ + g_X²(⟨q1²⟩ + ⟨q2²⟩) + h_D²⟨D²⟩, with
  ⟨χ²⟩ = ∫d³p f(p)/((2π)³E_p) [D]. Thermal check: (3λ_φ + 2g_X² + h_D²)/12 = 0.004175 = BRIEF value [C].
  With equal sharing over the four dof, C_avg = (3λ_φ + 2g_X² + h_D²)/4 = 0.012525.
- Restoration temperature: T_c = 1.547646 v (BRIEF, leading curvature zero), as instructed, with a ±30% band. The
  sibling `thermal_potential` module finds the resummed T_c/v = 1.557–1.564 (+0.6–1.1%) [Q]. Since
  T_RH,max ∝ T_c^{4/7}, the bound would rise by 0.35–0.6%.
- SM bath: g* = 106.75, with ν_R not thermalised (Dirac Yukawa ~1e-12, stated assumption). The bath is massless and
  thermalises instantly, with Bose/Fermi initial statistics. Final-state hidden Bose enhancement is neglected
  (hidden occupation ≪ 1).
- Hidden dof: g_h = 4 (φ, q1, q2, D), all minimally coupled at μ0 (ξ = 0 by matching).

## 2. Item 1: what couples the sectors

**Portals.** Use eight real scalars: four Higgs components and (φ, q1, q2, D). The code builds the exact quartic
tensor λ_abcd by the polarisation identity. It evaluates the one-loop scalar beta
16π²β_abcd = λ_abef λ_efcd + 2 perms.

- With all Higgs portals zero, every mixed Higgs–hidden entry of λ is zero (0 entries above 1e-14), and every mixed
  β entry is at roundoff level (max 7.1e-20) [C]. Gauge and Yukawa terms enter as diagonal anomalous dimensions, or
  carry only Higgs indices. The portal subspace is therefore fixed.
- Beyond one loop the same holds to all orders in non-gravitational perturbation theory: with no shared vertex the
  generating functional factorises, Z = Z_SM Z_hid, so no mixed counterterm can appear [D].
- Literature check [C, Q]: the single-singlet scalar part of 16π²β_λhs equals Lebedev arXiv:1203.0156 eq.(15),
  4λ_hs² + 12λ_hλ_hs + 6λ_sλ_hs, with relative difference 0. Lebedev's gauge/Yukawa part,
  6y_t² − (3/2)(3g² + g′²), multiplies λ_hs. All terms are proportional to the portal.
- The hidden betas from the tensor reproduce `rg_beta_functions.txt`: B_lp, B_lx, B_ld, B_c, B_d, B_alpha [C].
- The linearised portal matrix at μ0 (16π² units, scalar part, λ_h = 0.13) has diagonal 1.620, 1.568, 1.566 and
  off-diagonal entries ≤ 0.04 [C]. It is multiplicative, so zero is preserved.

**Graviton-induced local portals [E].** A local H†Hφ² counterterm at order κ² needs two powers of a mass. Rodigast
and Schuster (0908.2422) state that gravitational running of couplings needs massive fields. This gives
δλ ~ max(m_φ, m_h)²/(16π²M_P²) = 2.1e-35. Eliminating R on shell in ξ_H R H†H gives
λ_eff ~ ξ_H m_φ²/M_P² ≤ 3.4e-34 for |ξ_H| ≤ 0.1. Making such a portal's thermal mass equal to λ_φ v² would need
T ≈ 9.4e18 GeV = 3.9 M_P. Both are negligible.

**Nonminimal couplings.** Use the convention L ⊃ −½ξ_φRφ² − ξ_X R|X|² − ½ξ_D RD², with conformal value ξ = 1/6. The
one-loop running is 16π² dξ_a/dlnμ = λ_aacd(ξ_cd − δ_cd/6) [D]. At ξ = 0 [C] the rates per e-fold are
β_ξφ = −1.0575e-4, β_ξX = −2.9552e-5 and β_ξD = −6.544e-6. The β_ξX value agrees with the cross-review
derivation_scope agent [Q].

Hence ξ(μ) < 0 for μ > μ0. For example ξ_φ(1 TeV) = −2.44e-4 and ξ_φ(1e13 GeV) = −2.68e-3. A negative ξ_φ with
R > 0 adds a tachyonic mass, which deepens the broken minimum and is **anti-restoring**.

- **Radiation era.** R = (ρ − 3p)/M_P² comes from the trace anomaly, using the O(g²) pressure of each gauge factor:
  (ρ − 3p)/T⁴ = Σ_G (2/144)[d_A C_A + (5/2)Σ_Dirac d_R C_R] b_G g_G⁴/(16π²) [D]. Here b_3 = 7, b_2 = 19/6,
  b_1 = −41/6, one-loop running from M_Z, μ = 2πT. This gives 1 − 3w = 2.1e-3 at 100 GeV and 2.6e-4 at 1e10 GeV
  [C]. Davoudiasl et al. hep-ph/0403019 eq.(4) is a factor 2 larger in prefactor [Q]; scalar, Yukawa and quartic
  pieces are omitted [E].
  The magnitude |ξ_φR| reaches λ_φv² only at T ≈ 2.4e11 GeV, or 1.4e11 GeV with a 10× larger trace. It has the
  anti-restoring sign there, and at those temperatures H ≫ m, so the φ zero mode is frozen.
- **Inflation and matter-like reheating [P].** R = 12H² and 3H² respectively. Then ξ_φ12H²/(λ_φv²) = −0.29 at
  H = 1 TeV, so it is included in item 3. X becomes tachyonic at φ = v during de Sitter for H > 2.9 TeV [C], which
  lies above the item-3 thresholds.

**Verdict:** gravity only, through s-channel graviton exchange and the common FRW background. No portal, no
curvature-induced restoration.

## 3. Item 2: gravitational freeze-in from the SM bath

**Amplitude [D, C].** The vertex is −i h_μν T^μν/M_P (g = η + 2h/M_P), with propagator
iP/k², P = ½(ηη + ηη − ηη). The Newtonian limit, V = −Gm1m2/r, is checked in code. Then

M = (T^i·T^f − ½ tr T^i tr T^f)/(M_P² s),  T^f = k1k2 + k2k1 − η k1·k2 − 2ξ(qq − ηq²).

I computed it with explicit Dirac spinors (Weyl P_L projection; completeness to 2e-8) and polarisation vectors.
The fermion and vector T are traceless (0 at machine precision).

- Real scalar initial state: Σ|M|² = t²u²/(s²M_P⁴) exactly (relative error 3e-16).
- Ratios to the Clery et al. eqs.(83)–(85) shapes are constant: 4.000 (scalar), 1.000 (Weyl), 4.000 (vector) at all
  angles. The factor 4 is Clery's 1/2 × 1/2 identical-particle factor.
- Define σ̂ = ½∫dΦ₂Σ|M|² = κ s²/M_P⁴. The results [C] are κ_H = 6.6315e-4 (= 1/(480π) analytically),
  κ_Weyl = 1.6579e-4 and κ_V = 1.3263e-3.
- **Central-charge check [C]:** with symmetry factors S (½ for identical initial particles), S·κ for conformal real
  scalar : Weyl : vector = 1 : 3.000 : 12.000. This equals the Weyl-anomaly c ∝ N_s + 3N_Weyl + 12N_v. That formula
  is standard and was recalled, not fetched.

**Rates [D, C].** Massless initial states with distributions f give
Γ_pairs = Σ_i N_i S_i κ_i J3_i²/(3π⁴) T⁸/M_P⁴ and dρ/dt = Σ_i N_i S_i κ_i 2J3_iJ4_i/(3π⁴) T⁹/M_P⁴.
The thermal integrals are J3 = π⁴/15 or 7π⁴/120, and J4 = 24ζ(5) or (45/2)ζ(5). The species are 4 Higgs (S = ½),
45 Weyl (S = 1) and 12 vectors (S = ½).

| per hidden real dof | Higgs | Weyl | vector | total |
|---|---|---|---|---|
| pair rate / (T⁸/M_P⁴) | 1.914e-4 | 8.243e-4 | 1.148e-3 | **2.1641e-3** |
| energy rate / (T⁹/M_P⁴) | 1.467e-3 | 6.769e-3 | 8.802e-3 | **1.7038e-2** |

For the four hidden dof: **β_E = 0.068150**, pair rate 8.656e-3, particle-number rate 1.731e-2 T⁸/M_P⁴. The mean
pair energy is 7.873 T. A Monte Carlo of the produced quanta (boosted final states, quantum statistics, 2.4e6
quanta) gives ⟨E/T⟩ = 3.9363, against 3.9365 analytically, and **⟨T/E⟩ = 0.3420**. Maxwell–Boltzmann statistics would
lower the pair rate by 4.7%. The ξ dependence (Higgs channel only, since conformal initial states are blind to the
final ξ) is +0.30% at ξ = −0.0034 and −7.2% at ξ = 1/6 [C].

**Literature comparison [Q, C].** Clery et al. eq.(86) gives β₀ = 3997π³/41472000 = 2.98833e-3. Ours plus a second
copy of the Weyl channel reproduces it to 1e-17, i.e. exactly. Their R_0^T is therefore the pair-reaction rate per
real scalar with the 45-Weyl channel counted twice, equivalent to a Dirac spin sum per Weyl field. Our counting is
the one consistent with the central charges. This is recorded as an unresolved factor-2 discrepancy in one channel;
it moves T bounds by (1/0.724)^{1/4} ≈ 8%.

**(a) Instantaneous reheating at T_RH [D, C].** In radiation domination (H = c_H T²/M_P, c_H = √(π²g*/90)),
d(ρ_h/T⁴)/dT = −β_E T²/(c_H M_P³), so ρ_h/T⁴ = β_E(T_RH³ − T³)/(3c_H M_P³). For a formally self-thermalised hidden
gas, T_h⁴ = (30/π²g_h)ρ_h has its maximum at T = (4/7)^{1/3} T_RH:

T_h,max = [10β_E C7/(π² g_h c_H M_P³)]^{1/4} T_RH^{7/4} = 2.9027e-15 GeV^{-3/4} T_RH^{7/4},  C7 = (3/7)(4/7)^{4/3}.

The condition T_h,max < T_c gives **T_RH,max = 1.348e10 GeV**, or 1.099e10 and 1.566e10 GeV for T_c ∓ 30%. It
scales as **T_RH,max = 1.348e10 GeV (v/TeV)^{4/7}** (v = 0.1, 1, 10, 100 TeV give 3.6e9, 1.35e10, 5.0e10, 1.9e11).
The ODE check agrees to 1.0e-6. Asymptotically ξ_T ≡ T_h/T_SM → 1.71e-7 at the bound (row table in
`c3_freezein_history.json`). Instantaneous reheating to 1.35e10 GeV needs H_inf ≈ 255 GeV.

**(b) Inflaton-dominated reheating [C, P].** The inflaton decays perturbatively with width Γ_I into the SM bath. The
background ODE is in ln a; the hidden energy and ⟨χ²⟩ follow by quadrature (backreaction ≤ 1e-29). Check against
Mambrini–Olive eq.(21): T_max/T_RH = 1557 numerically versus 1576 [C, Q].

The production that sets T_h,max happens near T_max: d(ρ_h a⁴)/da ∝ a^{9/8}, while ρ_h itself peaks just after
a_max. Boundary (bath freeze-in only):

| H_e [GeV] | 1e3 | 1e5 | 1e7 | 1e9 | 1e11 | 1e13 |
|---|---|---|---|---|---|---|
| T_RH,max [GeV] | 1.09e10 | 2.93e9 | 8.15e8 | 2.27e8 | 6.32e7 | 1.76e7 |
| T_max,max [GeV] | 1.30e10 | 2.13e10 | 3.56e10 | 5.94e10 | 9.90e10 | 1.65e11 |

The fit is T_max,max ≈ 6.0e9 GeV (H_e/GeV)^{0.111} (v/TeV)^{4/9}. The slope in v is 0.44445, against an expected
4/9. With T_c ± 30%: at H_e = 1e9, T_max,max = 5.07e10–6.67e10; at H_e = 1e13, 1.41e11–1.86e11 GeV.

**(c) The gas does not self-thermalise [E].** At T_RH,max: n_h = 700 GeV³, E ≈ 3.9 T, and the largest contact
|M|² is 36λ_φ². This gives Γ₂₂/H = 1.7e-26, and Γ/H ∝ 1/T afterwards. The formal T_h above is therefore only an
upper characterisation. The physical criterion is item 5.

## 4. Item 3: inflation-era gravity with zero direct hidden branching [P]

**(i) Condensate scattering.** For a quadratic inflaton minimum, a minimally coupled massless χ in the conformal frame
has m_eff² = (ξ − 1/6)R. Here R = (ρ + 3ρcos2m_It)/M_P², so δm² = A cos(2m_I t) with A = −(1−6ξ)ρ/(2M_P²). The
golden rule gives a pair rate A²/(64π) [D]. Energy per real dof is (1−6ξ)² m_I ρ_I²/(128πM_P⁴). This agrees with
Clery et al. eqs.(22),(55),(62) [Q]. Mambrini–Olive eq.(7), as fetched, is a factor 8 lower (unresolved).

A numerical mode-function check (RK4, 1601 modes) gives a ratio of **1.0006** at A = 0.001 [C]. At A = 0.01 the ratio
grows with t_f (1.064 at t_f = 200 and 1.33 at t_f = 400); that growth is parametric (Bose) enhancement
(μ = A/4k), which the golden rule does not contain.

In matter domination without decay, ρ_h(a) = 0.0448 N_h m_I H_e³ x⁴(1 − √x) with x = a_e/a. The peak is
**ρ_h,peak = 7.75e-3 m_I H_e³** at a = 1.266 a_e (N_h = 4) [D]. Numerical ODE/analytic ratio 0.9993 [C].

The same process into the Higgs quartet, thermalised with g* = 106.75, gives 3.07e12 GeV for ρ_end = 1e64 GeV⁴ and
m = 3e13. Clery eq.(61) quotes 3.0e12 [C, Q]. **The hidden quartet receives exactly the Higgs-quartet energy.**

| m_I/H_e | H_e max, formal thermal (T_c −30% / nom / +30%) | H_e max, nonthermal Π, MD | H_e max, nonthermal Π, RD |
|---|---|---|---|
| 1 | 3.91 / 5.59 / 7.26 TeV | 2.7e7 GeV | 4.5e5 GeV |
| 3 | 2.97 / 4.24 / 5.52 TeV | 1.4e8 GeV | 1.3e6 GeV |
| 10 | 2.20 / 3.14 / 4.08 TeV | 8.6e8 GeV | 4.5e6 GeV |

That gas has Γ₂₂/H ≈ 2e-7 (m = H_e) at the peak, independent of H_e [E], so the formal thermal column does not
apply. The Π columns use χ² = (9N_h·0.4/128π)(H_e³/m)(a_e/a)², Π = C_avg χ², and the requirement
Π − λ_φv² > H² at some a.

**(ii) De Sitter fluctuations of the light hidden spectators.** The Starobinsky–Yokoyama equilibrium is
P_eq ∝ exp(−8π²U_eq/3H⁴) for (φ, q1, q2, D). U_eq includes ½ξ_a(μ = H)R χ_a² with R = 12H², and α is set to 0
(1e-8). The code computes:

- the marginal F(φ) = U_φ − a⁻¹ln∫d²q dD e^{−aU}, using an erfcx closed form for q and quadrature for D;
- the barrier S = a[F(0) − F_min];
- the odd-mode FP eigenvalue Λ₁ per e-fold, from a tridiagonal Schrödinger form in x = φ/H;
- the Kramers form of SY eq.(57), Λ₁/H = 1/(4π²N_x I_x);
- ⟨U⟩ from a factorised (φ, r) and (φ, D) quadrature, which agrees with a full 3D grid to 1e-4.

Checks on the pure double well [C, Q]: FP, Kramers and SY eq.(58), √2|m²|/(3πH) e^{−2π²m⁴/3λH⁴}, agree to 0.1% and
3% at S = 14.4, and to 0.4% and 7% at S = 6.3.

Results [C]:

| H/v | 0.10 | 0.16 | 0.20 | 0.24 | 0.28 | 0.32 | 0.36 | 0.40 | ≥ 0.44 |
|---|---|---|---|---|---|---|---|---|---|
| S | 653 | 96.7 | 37.8 | 16.8 | 7.99 | 3.80 | 1.67 | 0.54 | 0 (restored) |
| Λ₁ per e-fold | 4e-285 | 6e-44 | 1.2e-18 | 9.8e-10 | 4.4e-6 | 1.8e-4 | 1.1e-3 | 2.6e-3 | 4e-3 → 0.0175 |

- **Restoration by X fluctuations.** At φ = 0, X is massless in de Sitter. Its equilibrium variance is
  ⟨q1² + q2²⟩ = √(3/(2π³λ_X)) H² = 6.955 H², which induces Π_φ = g_X²·6.955H². The leading estimate is
  H/v = √(λ_φ/(6.955 g_X²)) = **0.3792**. The full marginal puts the origin's local stability at H = 378.99 GeV and
  global restoration at **H = 435.5 GeV**, with a trimodal P_eq in between [C].
- X equilibrates on 1/Λ₁(X) ≈ 356 e-folds [D, Q]. For a finite duration with X starting at zero,
  ⟨r²⟩ = NH²/(2π²), which gives restoration thresholds H/v = 0.81, 0.57 and 0.44 for N = 30, 60 and 100 [D].
- **Bimodal regime** (H < 0.436 v). Sign flips occur at rate Λ₁/2 per Hubble patch per e-fold. The count of flip
  events whose domains exited within Δ e-folds of the observable scale is
  N_dom(Δ) = (Λ₁/2)(e^{3Δ} − 1)/3 [E, heuristic O(1)]. Setting N_dom = 1 gives **H = 274.6, 221.0 and 183.8 GeV**
  for Δ = 5, 10 and 20 (N_obs = 60 [P]). Λ₁N_obs = 0.1 at H = 378 GeV.
- **Unimodal regime.** φ at the end of inflation is a nearly scale-invariant field centred on 0. Under a Gaussian
  level-set estimate, the probability that the sign differs between points 60 e-folds apart is ≈ 0.22–0.39 [E]. Its
  zero set is a wall network on all scales at formation. Persistence is not computed here.
- **Energy.** ⟨U_true⟩ = c H⁴ with c = 0.042–0.094 over H/v = 0.08–10 [C]. The formal thermal criterion
  T_h,max < T_c gives H < 3.67 TeV (2.56–4.77 TeV for ∓30%). It is weaker than the restoration criteria above.

All hidden scales are ∝ v, apart from the ξ(μ = H) logarithms, so these H thresholds scale ∝ v [D].

## 5. Item 4: direct inflaton → hidden branching B_h [P]

Both sectors are sourced by Γ_Iρ_I and redshift as radiation, so ρ_h/ρ_R = B_h/(1 − B_h) and
**T_h,max = (B_h g*/((1−B_h)g_h))^{1/4} T_max** [D]. An ODE check gives ρ_h/ρ_R = 1.0000e-12 and T_h,max/T_max to
3e-6 [C]. Requiring T_h,max < T_c gives **B_h,max = (g_h/g*)(T_c/T_max)⁴/(1+…)** = 2.15e-5, 2.15e-13, 2.15e-21,
2.15e-29 and 2.15e-37 for T_max = 1e4, 1e6, 1e8, 1e10 and 1e12 GeV. The ±30% T_c band multiplies these by
0.24 and 2.86 [C].

This assumes the hidden bath is thermal. Kinetic equilibrium at T_h = T_c needs ξ_T ≳ 4.3e-5 by my estimate, or
6.6e-5 from the thermal_potential module. Number-changing equilibrium needs ξ_T ≳ 0.068 [Q]. The thermal form is
therefore self-consistent only for T_max ≲ 2.3e7–3.6e7 GeV (kinetic) or ≲ 2.3e4 GeV (number-changing). Above that,
use the nonthermal bound below.

## 6. Item 5: nonthermal hidden populations

Π_φ = Σ_i c_i ∫d³p f_i/((2π)³E_p), with c = (3λ_φ, g_X², g_X², h_D²). This is the thermal part of the coincident
propagator: ⟨χ²⟩ = ∫d³p(1 + 2f)/((2π)³2E) [D]. The restoration criterion requested is Π_φ > λ_φv² = 1e4 GeV². I also
give a Hubble-damped variant, Π_φ − λ_φv² > H² at the same time, which the zero mode needs before it can respond [E].

- **SM-bath freeze-in gas [D, C].** Π(T) = 4C_avg γ_χ(T_RH³ − T³)T²/(3c_H M_P³) with γ_χ = 2·2.1641e-3·0.3420 =
  1.4804e-3. The thresholds are T_RH > **9.07e12 GeV** (maximum of Π) or **8.25e14 GeV** (Hubble-damped). On the
  formal thermal boundary of §3, Π/λv² ≤ 8.4e-15.
- **Condensate-scattering gas:** H_e thresholds as in §4(i).
- **Branching gas (E = m_I/2) [C].** With equal sharing, B_h,max(no Hubble) = 3.4e-24, 3.4e-20 and 3.4e-16 for
  (H_e, m_I) = (1e5, 3e5), (1e9, 3e9) and (1e13, 3e13) GeV at T_RH ≈ 6.4e7 GeV. The Hubble-damped values are 8.0e-22,
  7.2e-14 and 7.2e-6. Coupling only to φ lowers these by 2.4; coupling only to D raises them by 125. Full rows are in
  `c5_branching_nonthermal.json`.
- **Thermalisation [E].** Γ₂₂ = n|M|²/(32πs), with |M|² ≤ 36λ_φ² and s ≈ 2E². The freeze-in gas gives Γ/H = 1.7e-26.
  The condensate gas gives 2.0e-7 (m = H_e) and 2.2e-8 (m = 3H_e); Γ/H grows as a^{1/2} in MD, reaching 1 only at
  a/a_e ~ 3e13. The branching gas at its boundary gives 1e-9 to 1e-25. Number-changing processes are slower by a
  further ~λ/16π². None of these gravitational or dilute populations self-thermalises before it redshifts. Its
  effect on φ is the Π above.

## 7. Conditional history constraints (no exclusion or validation claimed)

For the unchanged action, Einstein gravity, T_c = 1.5476 v (±30%) and N_obs = 60, avoiding hidden Z2 restoration or
vacuum randomisation through these channels requires:

1. H_inf ≲ 0.436 v, the equilibrium de Sitter restoration by X fluctuations. With finite durations the threshold is
   0.44–0.81 v.
2. For no macroscopic wrong-vacuum regions in the bimodal regime: H_inf ≲ 184–275 GeV.
3. The pre-inflationary sign of φ fixes our vacuum, and it is preparation data.
4. If the hidden gas could thermalise, which it cannot via gravity: T_RH < 1.35e10 GeV (v/TeV)^{4/7}, or
   T_max < 6.0e9 (H_e/GeV)^{1/9}(v/TeV)^{4/9} GeV.
5. Nonthermal Π: T_RH < 9e12 GeV (bath) and H_e ≲ 4e5–9e8 GeV (condensate). Both are implied by items 1–2, since
   instantaneous reheating gives T_RH ≤ 1.15e10 GeV at H = 184 GeV.
6. Direct branching: B_h < 2.15e-29 (1e10 GeV/T_max)⁴ where the hidden bath is thermal, otherwise the nonthermal
   values of §6.

Conversely, H_inf ≳ 0.44 v restores Z2_φ during inflation through hidden-sector gravity alone. Whether that leads to
a persistent macroscopic network needs the post-inflation formation and persistence calculation; that is not done
here. Every gravity-only channel gives zero hidden U(1)_X charge.

## 8. Literature and sources (read level)

| Source | Read level | Used for |
|---|---|---|
| Clery, Mambrini, Olive, Verner, arXiv:2112.15214 | sections: eqs.(11)–(14), (20)–(22), (51)–(62), (79)–(90) via web-fetch extraction of the HTML (not the full PDF) | amplitude shapes, β₀, condensate rate, T_max^h eq.(61) checks |
| Mambrini, Olive, arXiv:2102.06214 | sections: eqs.(7), (10)–(23) via web-fetch | T_max/T_RH (eq.21) check; condensate rate (factor-8 discrepancy) |
| Starobinsky, Yokoyama, astro-ph/9407016 | full text of §§III–IV (pdftotext, local scratch) | P_eq, ⟨φ²⟩, Λ₁ = 1.36859√(λ/24π²)H, double-well eqs.(57)–(58) |
| Lebedev, arXiv:1203.0156 | §2 (eqs.2, 15), pdftotext | Higgs-portal β convention check |
| Davoudiasl, Kitano, Kribs, Murayama, Steinhardt, hep-ph/0403019 | eq.(4) only, pdftotext | trace-anomaly 1−3w (factor-2 prefactor difference noted) |
| Garny, Sandora, Sloth, arXiv:1511.03278 | abstract | context (gravitational DM) |
| Garny, Palessandro, Sandora, Sloth, arXiv:1709.09688 | abstract | context |
| Bernal et al., arXiv:1803.01866 | abstract | context |
| Giudice, Kolb, Riotto, hep-ph/0005123 | abstract | T ∝ a^{-3/8} scaling (taken via MO and own ODE) |
| Ema, Nakayama, Tang, arXiv:1804.07471 | abstract | context, end-of-inflation production |
| Herring, Boyanovsky, Zentner, arXiv:1912.10859 | abstract | IR 1/k³ spectrum of light minimally coupled spectators |
| Chung, Kolb, Riotto, hep-ph/9802238 | abstract | context |
| Ford, PRD 35, 2955 (1987), doi:10.1103/PhysRevD.35.2955 | not fetched | citation only |
| Clery, Mambrini, Olive, Shkerin, Verner, arXiv:2203.02004 | abstract | nonminimal couplings (our (1−6ξ)² derived independently) |
| Rodigast, Schuster, arXiv:0908.2422 | abstract | gravitational running needs masses |
| Harigaya, Mukaida, arXiv:1312.3097 | abstract | gauge-theory thermalisation (not directly applicable to a contact-quartic hidden sector) |
| Weyl-anomaly c ∝ N_s + 6N_Dirac + 12N_v | recalled, not fetched | central-charge check |
| Repository: GRUT_CANDIDATE_1.md §2–3; quantum REPORT.md, rg_beta_functions.txt; cross-review REVIEW.md and derivation_scope FINDINGS; Codex 135d535 README/FEASIBILITY (git show); sibling thermal_potential NOTES.md | read | inputs, T_c, ξ_T thresholds |

## 9. What failed or disagrees

- **Clery et al. β₀:** our pair rate is 0.7242 × theirs. The difference is exactly a doubled 45-Weyl channel
  (reproduced to 1e-17). Unresolved, and treated as a convention question; it shifts T bounds by about 8%.
- **Mambrini–Olive eq.(7), as fetched:** the condensate rate is 8× below Clery et al. and our derivation. Unresolved.
- **Mode-function check at A = 0.01:** the ratio drifts above 1 (1.064 → 1.33 for t_f = 200 → 400) because of
  parametric resonance. The golden rule therefore omits Bose enhancement. For m_I ~ H_e the gravitational resonance
  parameter is O(0.1–1), which leaves an O(1) uncertainty in the condensate numbers (§4(i)) at m_I ≈ H_e.
- **Trace anomaly:** our prefactor is half that of DKKMS eq.(4). The conclusion (anti-restoring, irrelevant below
  ~1e11 GeV) is insensitive to this.
- **"Λ₁N_obs = 1" threshold:** taken naively it gives 1.69 TeV, which is meaningless in the unimodal regime. It was
  replaced by the regime-aware thresholds of §4(ii).
- A first ODE implementation in log variables with atol = 1e-300 was too slow (> 4 min per scan). It was replaced by
  a background ODE plus quadrature, which is valid because the hidden sector does not backreact (checked ≤ 1e-29).

## 10. Not done

- Post-inflation evolution: φ oscillation and fragmentation, domain formation scale, wall persistence, network
  scaling, backreaction and observables. This module supplies only the formation and heating conditions.
- Nonperturbative gravitational preheating and Bose enhancement in condensate scattering; thermal masses of the SM
  bath; g* variation; hidden masses in the rates.
- Explicit graviton-loop portal diagrams (dimensional estimate only); two-loop ξ running; coupling running inside the
  stochastic potential (frozen at μ0 apart from ξ).
- Stochastic multi-field dynamics beyond the equilibrium marginal. Λ₁ uses a 1D reduction (adiabatic elimination of
  X and D). No lattice of the end of inflation was run.
- A full Boltzmann treatment of hidden thermalisation (only Γ₂₂ estimates); non-quadratic inflaton potentials (φ^k
  harmonics).
- A hidden charge or chemical potential: none is produced by gravity; any yield needs an outside-action source.

## 11. Free inputs on which the results depend

- **[P]** H_inf during the last N_obs e-folds, N_obs (60 used), the total duration of inflation (equilibration of X
  takes ~356 e-folds), and the pre-inflationary hidden field configuration (sign of φ, X amplitude).
- **[P]** H_e, inflaton mass m_I and potential shape near the minimum (quadratic used), Γ_I (sets T_RH, T_max), B_h,
  and which hidden dof the inflaton couples to.
- **[P]** Instant thermalisation of the SM bath at all times during reheating.
- The scale v (results scale as stated), and T_c: 1.5476 v nominal ±30%; resummed 1.557–1.564 v from the
  thermal_potential module.
- Matching ξ = 0 at μ0 and the RG scale choice μ = H (inflation) or μ = 2πT (radiation).
- The heuristic restoration criteria: Π > λv², the Hubble-damped variant, and the N_dom count.

## 12. Files

- `code/common.py`: constants, ξ running.
- `code/c1_couplings.py`: item 1, writes `results/c1_couplings.json`.
- `code/c2_graviton_rates.py`: amplitudes, rates and MC (`c2_graviton_rates.json`).
- `code/c3_freezein_history.py`: freeze-in histories and Π for the bath (`c3_freezein_history.json`).
- `code/c4_inflation_sources.py`: condensate, mode check and stochastic spectators (`c4_inflation_sources.json`).
- `code/c5_branching_nonthermal.py`: B_h, Π and thermalisation (`c5_branching_nonthermal.json`).
- `code/run_all.py`: everything, plus `results/headline.json`.
