# VERIFY: module hidden_heating (verifier key v_heating)

21 October 2026. Adversarial verification of `modules/hidden_heating`, with a cross-check against
`modules/inflation_spectator`. All code here is my own and was written before I read the module code. I read the
module's NOTES.md for the claim statements. After I had my own numbers, I looked at the module code only to locate
discrepancies, at `c3_freezein_history.py:176-182`. Only numpy and scipy are used, with OMP_NUM_THREADS=2. Every job
ran in under 6 min.

Reproduce from `verification/v_heating/code/`:
`python3 vh1_portals.py; python3 vh2_rates.py; python3 vh3_history.py; python3 vh4_desitter.py ABC; python3 vh4_desitter.py D;
python3 vh4_desitter.py E; python3 vh4_desitter.py F; python3 vh5_reheating.py; python3 vh6_summary.py`.
Results are in `results/vh*.json`, with a summary in `results/vh_summary.json`.

Labels: [D] derived by hand here, [C] computed by my code, [Q] quoted from literature, with the read level stated.

## Claim 1. The Higgs portals stay exactly zero under one-loop running if they are zero at mu_0

**Claim.** The three Higgs portals phi^2 H^+H, |X|^2 H^+H and D^2 H^+H vanish at mu_0 and stay zero. Every mixed
one-loop beta is zero. The single-singlet limit agrees with Lebedev eq.(15).

**Method [D, C].** I take V ⊃ λ_h(H†H)² + H†H[λ_Hφ φ² + λ_HX|X|² + λ_HD D²] plus the unchanged U, with eight real
fields. The general one-loop formula is
16π²β_abcd = λ_abef λ_efcd + 2 perms + Σ(γ_ae λ_ebcd + 3 perms) + G_abcd. Here γ and G (the gauge⁴ and Yukawa⁴ boxes)
carry Higgs indices only, because the hidden fields are gauge singlets with no Yukawas. Derived by hand:

- 16π²β_λHφ = λ_Hφ(12λ_h + 6λ_φ + 8λ_Hφ + 2γ_H) + 2g_X²λ_HX + 2h_D²λ_HD
- 16π²β_λHX = λ_HX(12λ_h + 8λ_X + 4λ_HX + 2γ_H) + 4g_X²λ_Hφ + 4αλ_HD
- 16π²β_λHD = λ_HD(12λ_h + 6λ_D + 8λ_HD + 2γ_H) + 2h_D²λ_Hφ + 2αλ_HX

Here 2γ_H = 6y_t² − (9/2)g² − (3/2)g′². Each beta is a homogeneous polynomial with no term free of portals. A single
beta is not proportional to its own portal alone (β_λHφ ∋ 2g_X²λ_HX, for example), but the portal vector is
multiplicatively renormalised. Zero is therefore a fixed point (ODE uniqueness).

**Numbers [C].**
- The explicit tensor matches the polarization tensor to 2e-16.
- In 20 random coupling sets, with random Higgs-only γ and G and zero portals, the mixed entries of λ and of β are
  exactly 0.0.
- The hand formulas match the tensor to 1.9e-13 (relative) with random nonzero portals.
- The single-singlet limit equals the scalar part of Lebedev 1203.0156 eq.(15), 4λ_hs² + 12λ_hλ_hs + 6λ_sλ_hs, with
  relative difference 1e-16. Lebedev's gauge/Yukawa piece, −(3/2)(3g² + g′²) + 6y_t², multiplies λ_hs, the same as my
  2γ_H. Lebedev's potential (eq. 2) has h² = 2H†H, so λ_hs = 2λ_Hφ.
- The linearised portal matrix (λ_h = 0.13) is diag(1.620, 1.568, 1.566) with off-diagonal entries 0.02, 0.04 and
  2e-4. It is identical to the module's.
- RK4 running of zero portals from 100 GeV to M_P gives exactly 0.

**Literature.** Lebedev arXiv:1203.0156, section level: §2.1 eq.(2) and §2.3 eq.(15) via pdftotext.

**Verdict: CONFIRMED.** One wording caveat: each beta is a linear combination of portals, not proportional to its
own portal. The module's matrix statement is correct. The all-orders statement (no shared vertex, so Z factorises)
holds for non-gravitational perturbation theory. Graviton-loop portals are a separate [E] estimate, and I did NOT
CHECK them.

## Claim 2. Gravitational freeze-in energy rate beta_E = 0.06815 T^9/M_P^4 into 4 real minimally coupled dof

**Method [C].**
- Canonical graviton g = η + 2h/M_P with reduced M_P, so L = −h_μνT^μν/M_P, and a de Donder propagator. The Newton
  limit check gives V = −m₁m₂/(8πM_P²r): ratio 1.000000. Reduced M_P is used consistently.
- Amplitudes are computed numerically with my own chiral spinors (completeness to 2e-15), a Weyl projector and
  explicit transverse polarizations. The vector T is traceless, gauge invariant and conserved (all exactly 0).
- Σ|M|² per species: scalar t²u²/s², Weyl (t−u)²tu/(4s²), Dirac twice the Weyl value, vector 2t²u²/s². All are in
  units of 1/M_P⁴ and agree with the analytic forms to 1e-15.
- κ (½ final-state factor included): scalar 1/(480π) = 6.6315e-4, Weyl 1/(1920π) = 1.6579e-4, vector
  1/(240π) = 1.3263e-3, conformal scalar initial state 1/(2880π).
- Central-charge check: S·κ for conformal scalar : Weyl : vector = 1 : 3.000 : 12.000. This is the Weyl-anomaly
  c ∝ N_s + 3N_Weyl + 12N_v. It is an independent normalization check that the Weyl (not Dirac) count is correct.
- Rates: Γ = N S κ J3²/(3π⁴), ρ̇ = N S κ 2J3J4/(3π⁴). My thermal integrals agree with π⁴/15, 7π⁴/120, 24ζ(5) and
  22.5ζ(5).

**Numbers [C].**
- Pair rate per hidden dof: 2.16405e-3 T⁸/M_P⁴ (analytically 5789π³/82944000).
- Energy rate per dof: 1.70376e-2. **β_E (4 dof) = 0.068150**.
- Particle-number rate (4 dof): 1.7312e-2 T⁸/M_P⁴. Mean pair energy 7.873 T.
- My own Monte Carlo (1.8e6 quanta, CM-frame angular distributions, boosts): ⟨E/T⟩ = 3.937 against 3.9365
  analytically, and ⟨T/E⟩ = 0.3422 against the module's 0.3420.
- ξ_f dependence: +0.30% at ξ = −0.0034 and −7.2% at ξ = 1/6 (Higgs channel only). Both agree with the module.

**Literature [Q].**
- Clery et al. arXiv:2112.15214 (section level: eqs.(11)–(14), (26)–(29), App. A eqs.(79)–(89), pdftotext of the
  arXiv PDF). Eq.(86): β₀ = 3997π³/41472000 = 2.98833e-3. Their |M̄^{1/2}|², eq.(84), equals ½ × the Dirac spin sum,
  but it multiplies N_{1/2} = 45, which counts Weyl fields ("6 (anti)quarks × 3 colours, 3 (anti)charged leptons,
  3 neutrinos"). Their fermion channel is therefore 2× ours.
  - Ours/theirs = 0.72417.
  - Ours with the Weyl channel doubled reproduces β₀ to 1e-16, which confirms the module's identification.
  - Additional point the module does not state: Clery use R^T as the particle-number source, dn_X/dt + 3Hn_X = R_X
    (eq. 26), while their eq.(11) with symmetry factors (½ for the final state) is the reaction rate. Per produced
    particle, ours/theirs is 1.448, not 0.724. Their condensate rates, eqs.(20) and (23), carry an explicit ×2 for two
    particles, but the thermal rate does not. This affects only number-density comparisons, not β_E.
- Bernal et al. arXiv:1803.01866 (section level: eqs.(4)–(6), (15), (26)–(31), pdftotext). Their eq.(31) is
  α = 3997π³/663552000 = β₀/16 = 1.868e-4, so ours/Bernal = 11.59. The literature is internally inconsistent by a
  factor 16 between the 2018 and 2022 papers of overlapping authors. Their amplitude prefactors, eqs.(26)–(30), are
  8× (fermion) and 4× (vector) below Clery's. I infer a graviton-field normalization issue (h/M_P coupling with
  vertex ½), but this is not verified.
- Garny–Sandora–Sloth arXiv:1511.03278 (section level, pdftotext). It uses the non-reduced M_p = 1.2e19 with an 8π/M_p²
  vertex (equivalent to 1/M_P² reduced). It gives only a massive-DM ⟨σv⟩ with Bessel functions (its eq.(11) region),
  so there is no direct massless T⁸ coefficient to compare. Its counting is N_0 = 4, N_1/2 = 45, N_1 = 12.
- Mambrini–Olive arXiv:2102.06214: grep of the pdftotext found no thermal graviton coefficient. Used only as
  context.

**Verdict: CONFIRMED** (β_E = 0.068150, reduced M_P, Weyl counting). There is a minor overclaim in the
sensitivity statement. The module says the Clery factor "moves T bounds by (1/0.724)^{1/4} ≈ 8%". That holds for T_h
at fixed T_RH. T_RH,max ∝ β_E^{−1/7} moves by only 4.7% (I computed 0.9549 if β_E is divided by 0.7242), and
T_max,max ∝ β_E^{−1/9} moves by about 3.6%.

## Claim 3. Formal thermal bound T_RH,max = 1.348e10 GeV (v/TeV)^{4/7}

**Method [D, C].** Radiation domination with H = c_H T²/M_P and c_H = √(π²g*/90), g* = 106.75. Then
d(ρ_h/T⁴)/dT = −β_E T²/(c_H M_P³), which gives ρ_h = β_E T⁴(T_RH³ − T³)/(3c_H M_P³). This is maximal at
T = (4/7)^{1/3}T_RH, with T_h⁴ = 30ρ_h/(π²g_h) and g_h = 4. Hence
T_h,max⁴ = 10β_E C7 T_RH⁷/(π²g_h c_H M_P³), C7 = (3/7)(4/7)^{4/3}. I also integrated the ODE in ln a independently.

**Numbers [C].**
- T_RH,max = **1.3478e10 GeV**; 1.0993e10 and 1.5659e10 GeV for T_c ∓ 30%.
- v scaling exponent 0.571429 = 4/7. v = 0.1, 1, 10, 100 TeV give 3.62e9, 1.35e10, 5.02e10 and 1.87e11 GeV.
- ODE check: T_h,max/T_c = 1 − 5e-9, at T/T_RH = 0.82986 against 0.82983 analytically. Late ξ_T = 1.710e-7.
- H for instantaneous reheating: 255.2 GeV.
- Inflaton-dominated reheating (my own ODE, matter + perturbative decay, bisection on Γ): T_max,max = 1.30e10, 2.13e10,
  5.94e10 and 1.65e11 GeV for H_e = 1e3, 1e5, 1e9 and 1e13 GeV. The slope is 0.1105 (1/9), and the v exponent is
  exactly 4/9 at H_e = 1e9. All match the module. My T at H = Γ is 1.18× the module's "T_RH" column, a definitional
  difference that leaves T_max unaffected.

**Verdict: CONFIRMED** for the instantaneous and the inflaton-dominated forms. The bound is formal: the module itself
shows the gas does not thermalise (claim 4).

## Claim 4. Nonthermal criterion Pi_phi > lambda_phi v^2 (T_RH > 9.07e12 GeV) and Gamma_22/H ~ 1e-26

**Method [D, C].**
- ∂²U/∂φ² at φ = 0 is −λv² + 3λ⟨δφ²⟩ + g_X²⟨q1² + q2²⟩ + h_D²⟨D²⟩, with ⟨χ²⟩_f = ∫d³p f/((2π)³E). This is the
  thermal part of the coincident propagator ∫(1 + 2f)/(2E).
- The thermal check gives (3λ + 2g² + h²)/12 = 0.004175.
- For the freeze-in gas, Σ1/E per volume from produced particles redshifting as n ∝ T³, E ∝ T, using my MC ⟨T/E⟩:
  Π(T) = 4C_avg γ_χ T²(T_RH³ − T³)/(3c_H M_P³), with C_avg = 0.012525 and γ_χ = 1.4810e-3. A direct quadrature
  agrees to 2e-16.
- Π is maximal at T = (2/5)^{1/3}T_RH.

**Numbers [C].**
- T_RH > **9.069e12 GeV**.
- Hubble-damped variant (Π − λv² > H²): **8.250e14 GeV**, with the optimum at T ≈ 8.4e9 GeV, deep in the dilute
  tail. My first grid missed it and gave 3.9e16; I corrected the grid, and this is recorded in the code. The closed
  form for T ≪ T_RH is T_RH = (4c_H²λv²/(A²M_P²))^{1/6} = 8.2496e14.
- On the formal thermal boundary, Π/λv² = 7.25e-15. The module's "≤ 8.4e-15" is the maximum over its
  non-instantaneous rows.
- Thermalisation: at the T_h maximum my n_h = 171 GeV³ (all 4 dof) gives Γ₂₂/H = **9.0e-27**. The module quotes
  1.7e-26 with n_h = 700 GeV³. Its code (`c3_freezein_history.py:177`) evaluates the asymptotic ρ_h/T⁴ at T = T_RH,
  which overestimates n_h by about 4×. That direction is conservative, and the conclusion is unchanged.

**Caveat.** At T_RH = 9.07e12 GeV, H(T at max Π) ≈ 6e7 GeV ≫ √(λv²) = 100 GeV. The zero mode is frozen, so
Π > λv² alone cannot move φ. The Hubble-damped 8.25e14 GeV is the physically more relevant threshold. Both are
heuristic, necessary-type criteria for randomization. Using 9e12 as a constraint is the conservative direction.

**Verdict: CONFIRMED** (9.07e12 and 8.25e14 GeV; Γ₂₂/H ~ 1e-26, with the module's n_h overestimated by about 4×).

## Claim 5. De Sitter restoration for H > 0.436 v, and no macroscopic wrong-vacuum domains for H < 184–275 GeV

**Method [D, C].**
- Starobinsky–Yokoyama equilibrium P ∝ exp(−8π²U/3H⁴) for (φ, q1, q2, D), with ½ξ_a R z_a² and R = 12H².
- ξ is run with my own one-loop matrix (16π²β_ξ = Σ λ_aacd(ξ_cd − 1/6)), giving ξ_φ(379 GeV) = −1.41e-4.
- Marginal over X: closed form with erfcx. Marginal over D: closed form with K_{1/4}.
- Λ₁ comes from a detailed-balance chain on the half line with odd-sector killing (eigenvalue), and from the
  continuous MFPT in log space for tiny Λ₁. The two agree to 1e-4 where both work.
- Validation: pure quartic Λ₁/√λ = 0.088924 against 0.088924 from SY (1.36859/√(24π²)). For the double well at
  250 and 300 GeV, eigenvalue/MFPT/SY eq.(58) agree to 0.02%/2–5%.
- Langevin (Euler–Maruyama, dN = 0.01, 60000 paths, seeded) for the fresh start φ = +v, X = D = 0.

**Numbers [C].**
- ⟨q1² + q2²⟩_eq = 6.9554 H². The leading estimate gives H/v = 0.37917.
- Full marginal: origin locally stable at **H = 378.9 GeV**, global restoration at **H = 435.5 GeV** (431.1 GeV with
  ξ = 0).
- S and Λ₁ at H/v = 0.10, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36 and 0.40:
  - S = 653.0, 96.66, 37.81, 16.84, 7.99, 3.80, 1.67, 0.54.
  - Λ₁ = 3.7e-285, 5.6e-44, 1.2e-18, 9.8e-10, 4.4e-6, 1.8e-4, 1.05e-3, 2.55e-3.
  - All match the module's table.
- N_dom(Δ) = (Λ₁/2)(e^{3Δ} − 1)/3 = 1 gives H = **274.7, 221.0 and 183.9 GeV** for Δ = 5, 10 and 20.
- Λ₁N₆₀ = 0.1 at 378.1 GeV.
- Finite duration from X = 0: the module assumes free growth, ⟨r²⟩ = NH²/2π². My Langevin with the λ_X quartic gives
  ⟨r²⟩/H² = 1.49, 2.82, 4.22, 6.05 and 6.67 at N = 30, 60, 100, 200 and 300. The g²⟨r²⟩ = λv² thresholds are
  therefore H/v = 0.82, 0.60 and 0.49 for N = 30, 60 and 100, against the module's 0.81, 0.57 and 0.44.
  ⟨r²⟩ equilibrates on about 200–300 e-folds. The module's 356 is the odd-mode time, which is consistent in order.

**Reconciliation with inflation_spectator.** The definitions differ.
1. **All-sizes count < 1e-3.** inflation_spectator: N_reg = (Λ₁/2)(e^{3N*} − 1)/3 < 1e-3 for the relaxed start, giving
   144.4 GeV (N* = 50) and 138.1 GeV (N* = 60). Mine: 145.4/138.5 GeV with its X = D = 0 path, and 144.2/137.9 GeV with
   hidden_heating's X, D marginal. My 1D Λ₁ reproduces its table within 30% (200, 251, 298, 335, 398 GeV), and my 1D
   Langevin reproduces its fresh-start p_-, for example 0.1837/0.2079 against 0.186/0.210 at 1 TeV.
2. **Late re-entry.** inflation_spectator: re-entry after 1 MeV (n′ = 17.77) < 1e-3 at 186.2 GeV (1D) and
   183.2 GeV (entropic). Mine: 187.1 and 183.7 GeV. This is essentially hidden_heating's Δ = 20, N_dom = 1 value
   (183.9 GeV).

Once the definitions match, the two modules agree to ≤ 1–2%. hidden_heating's "184–275 GeV" is a Δ-band:
- Only the 184 GeV end corresponds to "no domain re-entering after BBN".
- The 275 GeV end (Δ = 5) counts only domains larger than e^{-5} of the observable scale. At H = 275 GeV, O(1) or more
  domains are still expected to re-enter between about 1 MeV and the eV scale, including below T_dom = 8 keV, which
  inflation_spectator flags as wall-dominated or black-hole candidates.
- "No macroscopic wrong-vacuum domains for H < 275 GeV" is therefore definition-dependent and not conservative. The
  conservative value is about 184 GeV for N_obs = 60, relaxed start.

**Restoration (0.436 v) against percolation (1.35–1.78 TeV).** These are different questions:
- hidden_heating's 435.5 GeV is the shape of the equilibrium marginal. It needs X equilibrated (about 200–300
  e-folds) and φ relaxed (1/Λ₁ ≈ 60–250 e-folds just above threshold).
- inflation_spectator's p_- = p_c threshold is for a fresh +v start with X frozen over N* = 50/60.
- With matched definitions (fresh +v start, X and D fluctuating from zero, N* = 50/60), my Langevin gives
  p_- = 0.0051/0.0086 at 400 GeV, 0.029/0.044 at 500 GeV, 0.107/0.135 at 700 GeV, 0.203/0.231 at 1 TeV,
  0.279/0.305 at 1.4 TeV and 0.339/0.359 at 2 TeV.
- p_- reaches p_c = 0.3116 at **H = 1.70 TeV (N* = 50) and 1.47 TeV (N* = 60)**. With X frozen I get 1.80/1.64 TeV,
  reproducing inflation_spectator's 1.78/1.61 TeV.
- X fluctuations therefore lower the percolation threshold by only 6–10%. Above 0.436 v and below about 1.5 TeV, a
  +v patch is not randomised in 50–60 e-folds.

The two modules agree once definitions match. hidden_heating's §7 sentence "Conversely, H_inf ≳ 0.44 v restores Z2_φ
during inflation through hidden-sector gravity alone" holds only in equilibrium (long inflation at that H). For
N_obs-scale durations from a +v start it is an overclaim.

**Verdict: CONFIRMED** for the numbers (378.9, 435.5, 274.7, 221.0 and 183.9 GeV). **CORRECTED** for:
- the finite-duration thresholds: 0.82/0.60/0.49 v, not 0.81/0.57/0.44 v;
- the interpretation: equilibrium-only restoration, and the 275 GeV upper end is not a conservative "no macroscopic
  domain" bound.

The overdamped stochastic picture is marginal at m ~ H (m_φ/H = 0.77 at 184 GeV). Both modules state this.

## Claim 6. Direct-branching bound B_h < (g_h/g*)(T_c/T_max)^4

**Method [D, C].** Both sectors are sourced by Γρ_I and redshift as radiation, so ρ_h/ρ_R = B/(1−B) at all times. My
ODE gives this ratio to 1e-15. Then (T_h/T)⁴ = (g*/g_h)B/(1−B), which gives B < x/(1+x) with
x = (g_h/g*)(T_c/T_max)⁴.

**Numbers [C].**
- B_h,max = 2.1497e-5, 2.1497e-13, 2.1497e-21, **2.1497e-29** and 2.1497e-37 for T_max = 1e4, 1e6, 1e8, 1e10 and
  1e12 GeV.
- The ±30% factors are 0.2401 and 2.8561.
- Kinetic self-consistency (my order-of-magnitude contact-scattering estimate): ξ_T ≳ 4.3e-5, so T_max ≲ 3.6e7 GeV.
  This matches the module's 2.3–3.6e7 GeV. The number-changing threshold ξ_T ≳ 0.068 is quoted from a sibling module
  and NOT CHECKED.

**Verdict: CONFIRMED** (conditional on a thermal hidden bath, as the module states).

## Not checked

- Graviton-loop portals (≲ 3e-34).
- The trace-anomaly R in the radiation era.
- The condensate-scattering rates and the Mambrini–Olive factor-8 note.
- The branching-gas Π rows (c5).
- The number-changing threshold 0.068.
- The Gaussian level-set 0.22–0.39 estimate.
