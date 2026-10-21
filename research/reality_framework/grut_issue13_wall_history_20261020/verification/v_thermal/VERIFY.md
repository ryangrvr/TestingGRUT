# VERIFY: module thermal_potential (verifier key v_thermal)

21 October 2026. Adversarial verification of `modules/thermal_potential` (NOTES.md and results/*.json).
I wrote all code here from scratch and did not import or copy module code. I read only the module's NOTES.md and
results JSON, and the read came after I had my own numbers. Units: v = 1 = 1000 GeV.

Reproduce (about 4.5 min total on 2 threads):

```
cd research/reality_framework/grut_issue13_wall_history_20261020/verification/v_thermal
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 vt_transition.py   # 3.2 min -> results/transition_v.json
python3 vt_xd.py           # 8 s    -> results/xd_v.json
python3 vt_kz_therm.py     # 1 s    -> results/kz_therm_v.json
python3 vt_gap_scheme.py   # 45 s   -> results/gap_scheme_v.json
```

`results/verify_summary_v.json` holds the side-by-side comparison with the module's summary.json.

## Independent method (what differs from the module)

- **Couplings and RG.** I build the quartic tensor λ_ijkl = ∂⁴V by polarization in the basis (φ, q1, q2, D), with
  X = (q1 + i q2)/√2, so g_X²φ²|X|² = cφ²R²/2 and λ_X|X|⁴ = λ_X R⁴/4. I use the *generic* one-loop scalar betas:
  16π² β_λijkl = λ_ijmn λ_mnkl + 2 perms, and 16π² β_m²ij = λ_ijkl m²_kl. These reproduce every line of the
  repo's rg_beta_functions.txt exactly. I get β_α(μ0) = 1.26651e-8 and m_X²(μ = v) = −2.92085e-6 v².
- **J_B.** I derived the series myself from the n ≠ 0 Matsubara sum. The coefficient of (y²)^k is
  (T/2)(−1)^{k+1}/k · Γ(k−3/2)/((4π)^{3/2}Γ(k)) · 2ζ(2k−3)/(2π)^{2k−3} · 2π². The y⁶ coefficient is
  ζ(3)/(384π²) = 3.1717e-4, which equals (l+2)! in the module's form; (l+1)! would be 3× too large. Above
  y² = 4 I use the Bessel sum. Checks against direct quadrature:
  - y² ∈ [1e-4, 25]: ≤ 2e-15
  - negative y², using the real-part quadrature ln(2|sin θ/2|): ≤ 5e-15
  - J_B′: ≤ 2e-16, and ≤ 9e-9 for y² < 0 against a quadrature finite difference
- **T0** comes from root-finding the analytic curvature 2 dV/d(φ²) at φ = 0.
- **Order.** I scan for local minima on 20000 log-spaced φ ∈ [1e-7, 0.5] v at 44 temperatures in
  (T0, T0 + 2e-3]. When a first-order transition shows up, T_c comes from ΔV = 0.
- **X and D curvatures** are analytic. I use the Daleckii-Krein second derivative of Σ F(eig) (divided
  differences across the φ–X and φ–D mixing), not finite differences. A finite-difference cross-check on a 4×4
  eigenvalue V_eff agrees to 5e-8. The global search is a 33³ grid plus Nelder-Mead.
- **Adversarial extra.** I added a self-consistent (gap-equation) soft-φ daisy scheme that the module did not run.

## Claims

| # | Claim (as stated in module) | My method | My number | Verdict |
|---|---|---|---|---|
| 1a | Π_φ/T² = λφ/4 + g²/6 + h²/12 = 0.004175 | (1/24)∂²Tr H by hand, plus a numeric finite-difference Laplacian of Tr H | 0.004175000 (numeric 0.0041750000004) | CONFIRMED |
| 1b | Π_X/T² per real mode = g²/12 + λX/3 = 0.00116667 (complex X = 2 real modes) | same, using q1 and q2 separately; cross-check in complex form: 4λX⟨\|X\|²⟩ + g²⟨φ²⟩ with ⟨\|X\|²⟩ = T²/12 | 0.0011666667 for each of q1 and q2 | CONFIRMED (no factor-2 error) |
| 1c | Π_D/T² = h²/12 + λD/4 = 0.000258333 | same | 0.000258333 | CONFIRMED |
| 1d | Leading curvature zero T/v = 1.547646465 | √(λφ/Π_φ) | 1.5476464651 | CONFIRMED |
| 2a | AE T_c/v = 1.563466 at μ0 | own one-loop MSbar CW + exact J_B + AE daisy for all three species | 1.5634661935 | CONFIRMED (Δ < 3e-10) |
| 2b | Parwani T_c/v = 1.561520 at μ0 | own Parwani (m² → m² + Π in CW and J_B) | 1.5615198418 | CONFIRMED |
| 2c | DR T_c/v = 1.559307 at μ0 | own DR (AE for X and D; hard-only φ: J_B + (π/6)Re y³) | 1.5593071664 | CONFIRMED |
| 2d | Adopted range 1.5571–1.5641 (AE, Parwani and DR at μ0/2, μ0, 2μ0 and v, with running couplings); 0.61–1.07% above leading | own generic RG plus all 32 table cells | [1.5571397, 1.5641468]; +0.613% to +1.066%. All 32 cells (fixed and running, 4 schemes) agree to ≤ 2.8e-10 | CONFIRMED |
| 2e | Second order with φ_c = 0 in AE, Parwani and DR at all 8 (couplings × μ) settings | minima scan above T0 | No coexisting minimum in (T0, T0 + 2e-3] for any of the 24 resummed cells. λ_eff/λφ at T0: AE 0.579–0.597, Parwani 0.565–0.567, DR 0.980–0.984 (matches) | CONFIRMED for these three schemes |
| 2f | "In every resummed scheme the perturbative transition is second order" (generalization) | added a gap-equation soft-φ scheme: V_DR − (T/12π)Re(M² + 3λφφ²)^{3/2}, with M² = c_DR − (3λφT/4π)M, so curvature(0) = M² | **First order (spurious):** T_c/v = 1.560166, φ_c/T_c = 0.0258, window (T1 − T0)/T_c = 6.0e-4, M/T at T_c = 1.25e-3 | CORRECTED: the perturbative order is scheme-dependent. The second-order result holds for AE, Parwani and DR only. A self-consistent soft-φ resummation gives a spurious weak first-order transition with φ_c/T_c ≈ 0.026; the module's own estimate is 0.028. The window is about t_G (6.9e-4), so it is not deep inside the Ginzburg region the way the unresummed signal is. T_c = 1.56017 lies inside the adopted range. |
| 2g | Unresummed one-loop: first order with φ_c/T_c = 0.0105 at μ0 (0.0064–0.0133 across scales); T0, T1, T_c = 1.5491167, 1.5492335, 1.5492205; L/ρ_h = 3.5e-7 (≤ 5.6e-7); second order at running μ = v | own oneloop | T0 = 1.5491167, T1 = 1.5492335, T_c = 1.5492205, φ_c/T_c = 0.010541; scale range 0.00637–0.01329; \|L\| = 2.676e6 GeV⁴, i.e. L/ρ_h = 3.53e-7 (max 5.60e-7); running μ = v gives second order | CONFIRMED. My script has a sign-convention slip that prints L < 0; the magnitude agrees. |
| 2h | Zero-T minimum φ/v = 1.0001218818, curvature 0.0200127832 | own CW at μ0 | 1.0001218847, 0.0200127836 | CONFIRMED; the 3e-9 difference is minimizer tolerance |
| 2i | Nonperturbative order: continuous, 3D Ising, because the soft φ is a 3D Z2 scalar with λ_eff > 0 | not computed (no lattice); only λ_eff^DR/λφ = 0.98 confirmed | — | NOT CHECKED. The argument is standard and the module labels it an argument. NOTES §11 C1 says "re-breaks it continuously" without that caveat (mild overclaim). |
| 2j | T_c controlled to about ±0.3% | partial: the gap scheme gives 1.56017, inside the range; no two-loop | — | NOT CHECKED (no two-loop or lattice computation) |
| 3a | X never condenses: min m²_X,eff = 2.794e-3 (oneloop), 2.802e-3 (AE), 2.802e-3 (Parwani) v², reached at T_c, which is 0.998, 0.982 and 0.985 of g²φ² + Π_X T² | analytic Daleckii-Krein curvature along own φ_min(T), 31 T in [0, 3] plus T_c(1 ± 1e-3) | 2.79415e-3, 2.80146e-3, 2.80233e-3, all at T_c; ratios 0.9979, 0.9823, 0.9851 | CONFIRMED |
| 3b | D never condenses: min m²_D,eff = 1.0100e-2 v², reached at T = 0 | same | 1.00998e-2 at T = 0 | CONFIRMED |
| 3c | Global minimum at R = D = 0 for T ∈ {0.5, 1.0, 1.4, T_c, 1.7, 2.5} (AE and Parwani) | 33³ grid on [0, 1.6]³ plus Nelder-Mead from 5 off-axis starts | grid minimum always at R = D = 0; off-axis gap 3.5e-6 to 1.15e-5 v⁴; NM gives \|R\| ≤ 7.3e-7 and \|D\| ≤ 7.0e-7 (tolerance) | CONFIRMED. At tree level this is exact: every R and D term in U is ≥ 0 at μ0 (m_X0² = 0, m_D0² > 0). This holds for **zero charge only**. |
| 4a | N_H = [τ_Q^{1−ν/(1+νz)} τ₀^{ν/(1+νz)}/ξ₀]³; slopes 3.68 (z = 1) and 4.33 (z = 2) | re-derived: τ₀ε̂^{−νz} = ε̂τ_Q and N_H = (τ_Q/ξ̂)³ | slopes 3.6810 and 4.3275 | CONFIRMED. N_H counts domains per H⁻³, not per (4π/3)H⁻³ (+0.62 dex) — a convention that should be stated. |
| 4b | τ_Q = 1/H(T_c/ξ); H = √(π²g*/90) T²/M_P with reduced M_P = 2.435323e18 GeV and g* = 106.75; H(T_SM = T_c) = 3.4e-12 GeV | own | H = 3.41601e-12 GeV at T_c = 1559.307 GeV; τ_Q = 2.9274e11 GeV⁻¹ | CONFIRMED (reduced-M_P form is consistent) |
| 4c | ξ = 1 values: log10 N_H = 25.82 / 25.06 (z = 1, matched / naive) and 30.01 / 29.46 (z = 2); ξ_KZ = 725 GeV⁻¹; ε̂ = 2.5e-9; Kibble-Ginzburg 36.11 | own matched and naive freeze-out, κ_DR = 0.0041196 (own), t_G = 6.917e-4, ν = 0.629971 | 25.819 / 25.063; 30.010 / 29.464; 724.7 GeV⁻¹; ε̂ = 2.476e-9. All 7 ξ rows reproduce. | CONFIRMED |
| 4d | τ₀ band ×10^{±2} gives log10 N_H 23.5–28.1 (z = 1) and 28.3–31.7 (z = 2) at ξ = 1; C2 range 10^{8.8}–10^{31.7} | own; extra test with a collision-limited τ₀ ≈ 1/Γ_φ = 1.8e5 ξ₀ | with collision-limited τ₀: 31.9 (z = 1, matched) and 34.4 (z = 2, matched) | CONFIRMED for the stated band. GAP: a physically motivated τ₀ lies 5.3 dex outside the band and pushes N_H *higher*, so the upper end of C2 is not conservative. The qualitative N_H ≫ 1 is unaffected. |
| 5a | φ kinetic thermalization ξ_min = 6.6e-5 (Γ_φ/T = 4.98e-7 from 18λφ² + 12g⁴ + 6h⁴) | (i) reproduced the module's kinematics; (ii) my exact massless contact kinematics σv_Møl = \|M\|²S/(32πE1E2), with Bose ⟨1/E⟩, a typical-energy particle, or Boltzmann | (i) 6.63e-5 (coupling sum and channel counting confirmed); (ii) Γ/T = 7.6e-7 to 1.7e-6, giving ξ_min = 3.6e-5 to 5.4e-5 | CONFIRMED at order of magnitude. The module's ⟨s⟩ = 2⟨E⟩² average underestimates the rate by 1.5–3.4×, so 6.6e-5 sits at the conservative (high) end; the better value is about (4–5)e-5. |
| 5b | other thresholds: X 1.45e-4, D–φ exchange 0.0148 | module kinematics | 1.454e-4, 0.01483 | CONFIRMED (same kinematic caveat) |
| 5c | fugacity factors √(ζ(2)/Li₂(z)) = 1.125, 1.681, 4.00 for z = 0.9, 0.5, 0.1 | scipy spence | 1.12499, 1.68083, 4.00371 | CONFIRMED |
| 6 | Percolation: largest + and − clusters both span; ≥ 96.7% (L = 32) and ≥ 97.3% (L = 64) of their sign; wall fraction between spanning clusters 0.902 (L = 32), 0.910 (L = 64) | own scipy.ndimage labeling, open boundaries, 3 seeds | span: all True; fractions ≥ 0.9673 and ≥ 0.9734; wall fraction 0.893–0.903 (mean 0.898) at L = 32 and 0.908–0.909 at L = 64 | CONFIRMED (seed-level scatter; uncorrelated-site model only) |
| 7 | Literature: AE hep-ph/9212235, Parwani hep-ph/9204216, CMR 1612.00466, KLRS hep-ph/9508379, KPSV 1603.04436 (Δ_ε = 1.412625), WZZGD 1302.0421 (p_c = 0.3116077) | fetched the arXiv abstract pages | titles and authors match. Δ_ε = 1.412625(10) and p_c(site) = 0.3116077(2) are both in the abstracts. The AE abstract is about two-loop corrections to the ring-improved potential in gauge theory; calling the daisy scheme "AE" is conventional, and the module says so. | CONFIRMED (abstract level) |

## Overclaims and gaps

1. **"In every resummed scheme the perturbative transition is second order"** (NOTES §0, §5) is too broad. It holds for the
   three schemes implemented. A gap-equation soft-φ resummation gives a spurious first-order transition with
   φ_c/T_c = 0.026. The second-order result in AE and Parwani depends on the leading-order M(T0) = 0.0092 T. The φ-daisy
   quartic only turns the total quartic negative when M < 9λφT/(8π) = 0.0036 T, a margin of 2.6×. The exact soft mass
   vanishes at T_c, so perturbation theory cannot settle the order. Only the 3D-Ising argument can, and it is not computed.
2. **C1 (§11)** says "re-breaks it continuously" as a conditional result without repeating that the order is an
   argument, not a lattice computation. The module never claims the order is established nonperturbatively; §0, §5.1,
   §13 and summary.json say "argument, not a lattice computation". That is fine, but C1 should carry the same caveat.
3. **The KZ τ₀ band** (×10^{±2}) leaves out the collision-limited τ₀ ≈ 1/Γ_φ ≈ 10^{5.3} ξ₀. That value raises log10 N_H by
   about 6 (z = 1) or 4.4 (z = 2). The conclusion N_H ≫ 1 gets stronger, but the quoted upper ends (28.1, 31.7) are not bounds.
   N_H also counts per H⁻³ (convention, 0.62 dex).
4. **The thermalisation kinematics** (⟨s⟩ = 2⟨E⟩²) underestimate the 2→2 rate by 1.5–3.4×, so ξ_min is high by about
   1.2–1.85×. This is still within the order of magnitude, as the module claims.
5. **The positivity of X and D holds only at zero charge density.** The module says so; it is not a general statement.
6. **"A bath that is thermal at T_c stays thermal below it"** holds only while T_h ≳ the masses. Formation is not
   affected.

## Not checked

- Two-loop V_eff, the lattice 3D T_c and critical amplitudes, the ±0.3% T_c systematic
- 2→4 number-changing rates
- KZ dynamics classes (model C)
- the L = 128 percolation figure
- hidden_heating and other modules' inputs
