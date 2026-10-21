# v_production: adversarial verification of `modules/production_map`

21 October 2026. Verifier key `v_production`. I wrote all of the code in this directory from scratch, using numpy and scipy with 2 threads. I did not import or copy any module code. I read the module's code and JSON only after I had my own numbers, and only to find where we disagreed. Everything about inflation, reheating and the inflaton is **outside-action preparation data [P]**. Nothing here validates or excludes Candidate 1.

## Files

| File | Contents | Runtime |
|---|---|---|
| `vp_thermal.py` → `vp_thermal.json` | (1)–(3): Hartree gap equations at φ = 0 with a common fugacity, leading T_c(z), MB exit coefficient, G(z), F(z), T_SM,min | < 1 s |
| `vp_producer.py` → `vp_producer.json` | (4): ODE for perturbative MD-inflaton decay, C, ε, ⟨E⟩⟨1/E⟩, A_RD/A_MD, B_min per cell, B_max | 3 s |
| `vp_spectator.py` → `vp_spectator.json` | (5): equilibrium moments from a 180³ Gauss–Legendre quadrature, a cosmic-time check of the frozen-field term, 2000-sample conformal Verlet runs, J, yields, thresholds | 95 s |
| `vp_sr_noise.py` → `vp_sr_noise.json` | slow-roll ε coefficient from 2e6 samples, and the sampling error for 1500 samples | 5 s |
| `vp_collect.py` → `vp_results.json` | collected key numbers | — |

Conventions:
- M_P = 2.435323e18 GeV (reduced), with H² = ρ/(3M_P²).
- g* = g*s = 106.75 (SM only), s = (2π²g*/45)T³, ρ_SM = (π²g*/30)T⁴.
- Γ_I = H(T_RH) = √(π²g*/90) T_RH²/M_P.
- g_h = 4 real hidden dof (φ, q₁, q₂, D).

## Claims and verdicts

### (1) Kinetic T_c(z) — CONFIRMED

**Derivation.**
- For one massless real dof, ⟨χ²⟩ = ∫d³k/(2π)³ f/k = T²Li₂(z)/(2π²), with f = 1/(z⁻¹e^{k/T} − 1).
- Summing the four dof with a common μ_N gives Π_φ = (3λ_φ + 2g_X² + h_D²) T_h² Li₂(z)/(2π²).
- Setting Π_φ = λ_φv² gives T_c(z) = v √(λ_φ/0.004175) √(ζ(2)/Li₂(z)) = **1.547646 v √(ζ(2)/Li₂(z))**.
- At z = 1 this reduces to Li₂ = ζ(2), so Π/T² = 0.004175.

**Hartree check.** I solved the self-consistent gap equations at φ = 0, with m_φ = 0 on the boundary, m_q² = g²I_φ + 4λ_X I_q and m_D² = m_D0² + h²I_φ + 3λ_D I_D. Thermal parts only.
- T_c(1) = **1.557546 v**, which equals the module's 1.5575.
- The ratios T_c(z)/T_c(1), Hartree to Hartree:
  - z = 0.9: **1.11949**
  - z = 0.5: **1.67062**
  - z = 0.1: **3.97843**
  - These match the claimed 1.119, 1.671 and 3.978.
- The leading-order factors √(ζ(2)/Li₂) are 1.12499, 1.68083 and 4.00371.
- Hartree/leading is at most 1.0064 (at z = 1) and tends to 1 as z → 0. This matches the module's "(1 + ≤ 0.7%)".

**Caveats and positive checks.**
- Mixing the leading prefactor 1.5476 with Hartree ratios is harmless at the ≤ 0.7% level.
- The fixed-N versus fixed-μ curvature question: ∂²F/∂φ²|_n − ∂²Ω/∂φ²|_μ = (∂n/∂φ)²/(∂n/∂μ). This vanishes at φ = 0 by Z2. So the curvature criterion is ensemble-independent at the symmetric point. The same holds for a collisionless gas at fixed f.

### (2) MB-regime exit T_SM,end = 92.33 GeV √ε/Y_N (v/TeV) — CONFIRMED

**Derivation** (MB, massless, equal sharing over 4 dof):
- ⟨χ²⟩ per dof = n₁/(2T_h), so Π = c_eff n_h/(2T_h) with c_eff = (3λ_φ + 2g_X² + h_D²)/4 = 0.012525.
- Energy: ρ_h = 3n_hT_h. Hence Π = (3/2) c_eff n_h²/ρ_h.
- With n_h = Y_N s and ρ_h = ε ρ_SM: Π = K Y_N² T²/ε, where K = (3/2) c_eff (s/T³)²/(ρ_SM/T⁴) = (3/2) c_eff (8π²g*/135) = **1.172986**.
- Π = λ_φv² then gives T_SM,end = v √(λ_φ/K) √ε/Y_N = **92.3323 GeV √ε/Y_N** (v/TeV).
- Π ∝ T² at fixed (Y_N, ε), so the phase is restored above T_SM,end.
- (Y_N, ε) is conserved only while both sectors are relativistic, g* is constant and SM entropy is conserved.

**Accuracy of the MB form.** The exact massless Bose result relative to MB is (Li₂Li₄/Li₃²)^{−1/2}:
- 0.9009 at μ = 0 (the claim says "10% below");
- 0.9955 at μ/T = −2 (the claim says "< 1%").

Both confirmed. My exit-table rows agree with the module's to ≤ 0.3%. The module includes thermal masses in n and ρ; I use massless n and ρ.

**Convention note.** g* = 106.75 omits the ν_R that the action contains. If 3 Majorana ν_R were relativistic and thermalized, g* would be 112.0. The coefficient would then be 2.4% lower, and T_SM,min 1.6% lower. Whether ν_R thermalize is not fixed by the action. This is not an error.

### (3) T_SM,min = 2.078e6 GeV (μ_N = 0, Y_N = 4.373e-12) — CONFIRMED

- G(z) ≡ T_SM,end Y_N^{1/3} = T_c(z)(n_h/(T_h³ s/T³))^{1/3}.
- G(1) = 340.02 GeV here (module 339.9 with thermal masses in n).
- G increases monotonically as z decreases: 340, 364, 434, 593, 877 and 4075 GeV for z = 1 down to 1e-6. So the minimum is at z = 1 and Y_N = Y_max.
- T_SM,min = **2.0793e6 GeV** (Hartree, massless n), against the module's 2.0783e6. With the resummed T_c range: **2.0787–2.0880e6 GeV**, against the module's 2.0777–2.0871e6.
- F(1) = ε_th/Y^{4/3} = 16.50, against the module's 16.53. The 0.2% difference is the mass correction in n.
- ξ at exit = 7.49e-4. Confirmed.
- This is a bound only for kinetically thermal baths with μ_N ≤ 0. The module states this and flags the overpopulated (condensate) case as unbounded.

### (4) Heavy-inflaton producer [P]

**Coefficient C in Y_N = 1.5 C B_h T_RH/m_I.**

Verdict: **CONFIRMED.**

- A dimensionless ODE (ρ in units of 3M_P²Γ²) shows T_SM/T_RH = (ρ_r/3M_P²Γ²)^{1/4} exactly. So Y m_I/(B T_RH) = (3/4)ñ/(B r_r^{3/4}), and **g* cancels**.
- C = **0.939201** for H_e/Γ ≥ 1e8. It is 0.9391 at 1e4, 0.9332 at 100 and **0.9202** at 30.
- So 1.5C = **1.4088**.
- This holds for B ≪ 1. In general Y ∝ B(1 − B)^{−3/4}: C = 1.016 at B = 0.1.

**ε = B/(1 − B).**

Verdict: **CONFIRMED.** I get 0.111111 at B = 0.1.

**ξ_formal = m_I/(6C T_RH) = m_I/(5.6 T_RH).**

Verdict: **CONFIRMED.** I get 6C = 5.635, ⟨E⟩/T_SM = 0.5324 m_I/T_RH and ⟨E⟩⟨1/E⟩ = 2.267, which gives Π_free/Π_th = 1.511.

**A = 24 c_eff B Γ M_P²/m_I², constant through reheating.**

Verdict: **CONFIRMED, and strengthened.** Write J = J̃·3M_P²Γ²/m_I² and h = H/Γ.
- **J̃ = 8Bh is an exact solution** of J̃′ = −2J̃ + 4B r_I/h for *any* mix of ρ_I and ρ_r (B ≪ 1), because h′ = −(3r_I + 4r_r)/(2h).
- So A is exactly constant through MD, reheating and RD once the start-up transient (h/h_e)^{1/3} has died away.
- Numerically A_RD/A_MD = 0.99993 at h_e = 1e12.

**B_max = 3.10e-12 (m_I/T_RH) (relic edge).**

Verdict: **CONFIRMED.** Y_max/(1.5C) = **3.1041e-12**, and 3.168e-12 at C = 0.920.

**B_min = 8.0e-17 (m_I/T_RH)².**

Verdict: **CONFIRMED given the stated criterion.**
- Analytic: 2√λ_φ v m_I²/(24 c_eff √(π²g*/90) M_P T_RH²) = **7.985e-17**.
- Full-ODE max_t[c_eff J − λv² − H²] = 0 in seven cells, covering MD and RD optima, gives 7.99–8.04e-17. The optimum is at H ≈ 100 GeV.
- The hidden_heating c5 cell (Γ = H(1e8 GeV), H_e = 1e9, m_I = 3e9) gives 7.2190e-14, against the module's 7.21996e-14.
- **Window closure** at m_I/T_RH = 3.89e4. Confirmed.
- **Caveat:** the "zero mode responds" condition Π − λv² > H² is a heuristic. With (3H/2)² (critical damping of φ̈ + 3Hφ̇) B_min becomes **1.20e-16**, a factor of 1.5. The lower edge is therefore an O(1)-criterion edge, not a sharp computed threshold.

**μ_N/T_h ≈ −23 to −40.**

Verdict: **CORRECTED (minor).**
- The module's own late-RD formula z = 115.5 Y/ξ³ ≈ 2.9e4 B(T_RH/m_I)⁴ at B = B_max gives:
  - −23.1 at m_I/T_RH = 10;
  - −36.9 at 1e3;
  - **−43.8 at 1e4**, which is the (1e13, 1e9) cell of the module's table;
  - −47.9 at window closure.
- The module's −40 lower end is the formal state at an exit that happens during MD reheating. As a statement about "the windows" the range is −23 to −44 (−48 at closure).

**Internal inconsistency in the JSON.**
- `results/headline.json` and `a2_producer.json` report "k (B_min fit) = 1.258", "B_min = 1e-16 (m_I/T_RH)²" and "window exists iff m_I/T_RH < 3.09e4". These are medians over *both* H_e variants.
- NOTES quotes k = 1.00, 8.0e-17 and 3.9e4 for H_e = m_I/3, and 1.4–1.8 for low H_e. The NOTES numbers are the correct per-variant ones.
- The JSON summary strings should be relabelled.

### (5) Spectator energy fraction after instant reheating [P]

**Frozen-field term a′²φ₀²/2 (the "critic omitted it" claim).**

Verdict: **CONFIRMED.**
- In RD (a″ = 0) with a_e = 1 and χ = aφ, the equation is χ″ + λχ³ = 0, which is autonomous.
- The initial χ′ = a′φ₀ + a²φ̇₀ = Hφ₀ + φ̇₀. So E_c = H²φ₀²/2 + U + O(φ̇) is conserved, and late-time ρ_h a⁴ = E_c.
- An independent **cosmic-time** integration of φ̈ + (3/2t)φ̇ + λφ³ = 0 confirms this:
  - late ρa⁴/E_c = 0.99985 (λ = 0.01, φ₀ = 3H);
  - 1.00040 (λ = 0.001, φ₀ = 6H);
  - 0.99974 (with slow-roll velocity).
  - The potential-only estimate is 23–57× too small.
- Physically, the frozen field keeps V(φ₀) while the radiation redshifts until H ~ √λ φ₀. That gives ε ~ φ₀²/M_P², independent of λ.
- The critic's (c/3)(H/M_P)² is the end-of-inflation ratio, and it misses this.

**Equilibrium variances.**

Verdict: **CONFIRMED.** These come from Starobinsky–Yokoyama P ∝ exp(−8π²U/3H⁴), exact quartic, 180³ quadrature (stable against 120³ to 1e-14).

| Quantity | Value (H units) |
|---|---|
| ⟨φ²⟩ | 0.68130 |
| ⟨q₁² + q₂²⟩ | 5.45678 |
| ⟨D²⟩ | 4.14843 |
| Σ | **10.2865** |
| ⟨U₄⟩ | **0.0379954** |

- ⟨U₄⟩ equals the exact identity ⟨φ·∇U⟩ = 4U ⇒ ⟨U⟩ = 4·3/(32π²) H⁴.
- The single-field result is 0.131765 H²/√λ.
- A 400k-sample Monte Carlo agrees to 0.2%.

**ε coefficient 1.727.**

Verdict: **CONFIRMED for zero initial velocity: 1.72708.** It is **CORRECTED to 1.7103 with the slow-roll velocity.**
- φ̇₀ = −∇U/(3H) adds H φ₀·φ̇₀ = −⟨φ·∇U⟩/3 = −4⟨U⟩/3, plus |∇U|²/18.
- So c_ε = Σ/6 − ⟨U⟩/9 + ⟨|∇U|²⟩/54 = **1.71028** (quadrature). A 2e6-sample Monte Carlo gives 1.7108.
- The module's formula "+ O(slow roll)" is mislabelled. The slow-roll cross term is the **same order** as the retained ⟨U⟩ term.
- The module's "sampled incl. slow-roll 1.7268 agrees with analytic 1.7272" comes from **1500 samples with SEM 0.027**. The agreement is noise and does not test the slow-roll term.
- Numerical effect: the floor threshold moves 2.013e11 → 2.023e11 GeV (+0.5%). The quoted thresholds are unaffected at two significant figures.

**Ratios.**

Verdict: **CONFIRMED.**
- Frozen/potential term = 135.4.
- c_ε/(0.07/3) = 74.0. This depends on the critic's c: for c = 0.04–0.09 it is 58–130.

**Number J and yields.**

Verdict: **CONFIRMED.**
- From 2000 samples, Verlet dt = 0.02, average over τ ∈ [1000, 3000]:
  - J = **20.86 ± 0.26 H_e³**, against the module's 20.9 ± 0.3;
  - shares φ/X/D = 11/43/46%;
  - E_c/J = 0.248 H.
- Energy drift is ≤ 4e-5.
- The estimator calibrates to the exact single-field adiabatic invariant to 0.08%.
- Yields: Y_J = **2.82 (H/M_P)^{3/2}** and J/n_eq = 15.3.
- Floor: ξ_E = 2.6056 (H/M_P)^{1/2} and Y_floor = **0.18404 (H/M_P)^{3/2}**.

**Thresholds 2.0e11 / 3.3e10 GeV.**

Verdict: **CONFIRMED.**
- Y_N = 4.373e-12 is reached at H = **2.013e11 GeV** (floor) and **3.26e10 GeV** (J).
- Y at 1e13 GeV is 1.53e-9 (floor) and 2.35e-8 (J).
- √Σφ² ≲ 6.47e11 GeV. Confirmed.
- The thresholds scale ∝ M_P(Y_max/coef)^{2/3}, using the reduced M_P consistently.

**Critique "number freezes at amplitude ~ v is not right".**

Verdict: **CONFIRMED in substance, overstated in tone.** The comoving adiabatic invariant is fixed from oscillation onset. But the critic's ρ/m evaluated at amplitude v differs from it only by O(2): J/n_crit = (4/3)m_φ/ω_quartic(v) = 2.23 for single-field φ. So the critic's rule is a factor ~2 prescription, not a qualitative error. The module also says the critic's ε and its quantum-energy errors partly cancel, which I find consistent.

**Not checked.** Non-instant (MD) reheating Y ∝ T_RH, the running-ξ thresholds (1.78e11 / 2.92e10), the 60-e-fold variant, the isocurvature estimate [E, R], the §5.5 restoring-bath statements, the memory [E] criterion, the overpopulated 5.8e5 GeV placeholder, and the elastic-rate kinetic-exit region.

## Overclaims and convention flags (summary)

1. The ε formula's "+ O(slow roll)" is mislabelled. The slow-roll velocity enters at the order of the retained ⟨U⟩ term. The "sampled vs analytic within 2e-4" agreement is within a ±0.027 sampling error. Use c_ε = 1.71–1.73 (zero vs slow-roll end velocity).
2. The B_min lower edge depends on the heuristic Π − λv² > κH². κ = 1 gives 8.0e-17 and κ = 9/4 gives 1.2e-16. It should carry an O(1.5) criterion label.
3. The JSON summary (k = 1.258, 1e-16, 3.09e4) mixes the two H_e variants and contradicts NOTES (1.00, 8.0e-17, 3.9e4).
4. The μ_N/T range "−23 to −40" refers to exit-time states. The late-RD formula gives down to −44 in the table's windows.
5. g* = 106.75 excludes the action's ν_R. If they were thermalized, T_SM,end shifts by −2.4% and T_SM,min by −1.6%.
6. The reduced Planck mass is used consistently in H² = ρ/3M_P², Γ = H(T_RH) and ρ_r = 3M_P²H_e². No convention error found.

## What failed in my own runs

- My first producer run used atol = 1e-300 and stalled in the stiff post-decay tail. I killed it and rescaled the tolerances.
- My first B_min for the (1e13, 1e12) cell (9.5e-17, optimum at H = 181 GeV) was a truncation artifact: the integration stopped before H = 100 GeV. After extending the range it is 8.04e-17.
