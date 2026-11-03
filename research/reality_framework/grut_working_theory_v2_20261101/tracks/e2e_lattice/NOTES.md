# e2e_lattice: the post-waterfall history of HFa/HFb on a 3D lattice

Claude, end-to-end pass (GitHub issue #13 context). Track brief: `BRIEF_E2E.md`, section "e2e_lattice". This session resumed
work cut off by a container restart. The binary was still current, finished runs were reused, the 96³ run was resumed
from its restart file, and one more convergence run (`C_N96L63`) was added.

Claim labels: [D] derived, [C] computed on the stated model, [Q] quoted, [E] estimate, [A] assumption or closure,
[R] retrospective comparison.

Every potential term, closure, cutoff, seed and switch time here is a **model choice**. None of them follows from
closure. Nothing here is a validation, an exclusion, a new theory or a GRUT prediction. Data comparisons are
retrospective.

## 0. Headline (design scale H = 1.30×10¹¹ GeV unless stated; times in H_d⁻¹ from σ = σ_c; a = 1 there)

1. **No walls form in any closure [C].**
   - The −v volume fraction peaks at 1.8% (local Pauli), 3.0% (skeptic saturation) and 5.7% (no drain). This is at
     t ≈ 0.35–0.49, far below the ≈ 0.31 site-percolation threshold of the simple cubic lattice (textbook value, not
     re-verified).
   - No −v cluster ever wraps the box. The largest one holds ≤ 9×10⁻⁴ of the volume.
   - The 0.05 H⁻¹ time-averaged χ is never negative anywhere (wall area density from it is 0 in P0). The −v sites are
     sub-period overshoots of χ oscillating about +v. They die off with an e-folding time of about 0.3 H⁻¹, and none
     remain at t = 1.9 (P0) or t ≈ 3 (P0GM).
   - Bias extrapolation [D/E]: ΔV = 2c_b v gives ℓ_bias = σ_w/ΔV = 306 H⁻¹. Had a scaling network formed, it would
     dominate at H = 0.021 (t ≈ 23) before decaying at H = 0.0033 (t ≈ 150). This reproduces mf_critic's 23/150 and
     is why the "no network" result matters.
   - Super-box tail [E]: the bias displacement is x₀ = 2.0×10³ (H/2π) at the start, and χ is light for only
     ΔN ≈ 4×10⁻⁴ e-folds above σ_c. So super-Hubble −v domains need a ≳ 2000σ fluctuation.
   - Sub-Hubble seeds [E]: at the start the box-mean χ is 9.4 times its UV vacuum rms. Rare −v seeds of size ~1/m_χ
     therefore exist somewhere, but tension and pressure both close them (−v is the higher vacuum) and the lattice shows
     none surviving.
   - So the "walls not excluded" item of mf_critic §6 is answered on this model [C, closures A]: no wall network forms.
2. **Fragmentation [C].**
   - ε_H = 1 at t = 0.194.
   - The coherent (box-mean-field) fraction f_coh drops below 0.1 at t = 0.49 (P0) / 0.47 (P0GM).
   - The gradient fraction peaks at 0.148 (P0) / 0.221 (P0GM) at t ≈ 0.5. It falls to 0.019 / 0.028 at t = 1.5 and to
     0.010 / 0.028 at t = 3.
   - The left-over field energy is matter-like: w_cyc → f_rad/3.
3. **Reheating depends on the blocking closure, and the closure is the bracket [C, A].**
   - With *local* Pauli blocking (each site blocked against its own capacity F_sat(E(x))·E(x)): f_rad = 0.61 at
     t = 1.5 and 0.60 at t = 3, with w_cyc = 0.21.
   - With the mf_hybrid skeptic's saturation applied to the *box-mean* condensate energy (the homogeneous F_sat model):
     f_rad = 0.19 / 0.25 and w_cyc = 0.07 / 0.09. This is close to the homogeneous run (0.168 at t = 0.6).
   - Fragmentation itself does not unblock the drain. Only the local blocking closure does.
   - K from the fluid continuation (HB.stage2): **3.6–4.0 (local) and 5.7–6.1 (box-mean)** at Γ_th = 0, both with
     +4.5% from spacing. With Γ_th = 10 H the values are 2.8 / 3.9.
   - For comparison: mf_chain's skeptic-corrected K = 4.3–9.8 [Q], and e2e_kinetic's computed K = 17.7 [Q]. The
     kinetic track's thermalisation is not in this lattice.
4. **⟨R⟩ and its spread [C].**
   - R/H_d² = 12 at the slice.
   - The box rms of R(x) peaks at 3.5 (P0) / 4.6 (P0GM) at fragmentation.
   - Cycle-mean ⟨R⟩ at t = 1.5 is 0.10 (P0) / 0.25 (P0GM), with spatial rms 0.40 / 1.21. At t = 3 these are
     0.03 / 0.07, with rms 0.10 / 0.35.
5. **Spectator [C].**
   - s does not develop structure beyond its own vacuum noise: δs/s_* = 6.2×10⁻⁴ → 1.9×10⁻⁴, redshifting, the same in
     every closure.
   - The mass inhomogeneity sits at k ~ m_χ ≫ m_s, so s responds only to the box mean.
   - Spatial spread of m_s²:
     - HFa (g_s χ)²: 9.2 ± 2.4 at t = 0.5, ± 1.3 at t = 1, ± 0.56 at t = 3 (P0); ±1.8 / ±0.81 at t = 1 / 3 (P0GM).
     - HFb (9 − 3R/4): ± 1.9 at t = 0.5, ± 0.67 at t = 1, ± 0.08 at t = 3 (P0); ±1.6 / ±0.26 (P0GM).
   - The box-mean s history matches across closures to about 1% up to t = 0.5. After that the P0 and P0GM means differ
     (for example −0.117 vs −0.106 at t = 1).
6. **GW [C on the lattice, not converged in spacing].**
   - ρ_GW a⁴ freezes by t ≈ 1 (constant to 1% from t = 1 to 2.5).
   - Ω_GW at emission is 4×10⁻⁶–9×10⁻⁶ (P0) and 2×10⁻⁵ (P0GM).
   - Today: **Ω_GW h² ≈ 1×10⁻¹¹ – 5×10⁻¹¹ with this track's K (3.6–5.7)**, or 0.7×10⁻¹² – 5×10⁻¹² with e2e_kinetic's
     K = 17.7.
   - **Peak frequency 1–3 GHz.** It sits at the lattice cutoff; the χ-scale frequency is 0.4–1.6 GHz.
   - This is ≥ 4 orders below the ΔN_eff ceiling (Ω_GW h² ≲ 5.6×10⁻⁶ ΔN_eff ≈ 1.7×10⁻⁶ for ΔN_eff < 0.3) [R], and far
     above the frequency bands of LIGO, LISA and PTAs.
7. **Strong-holder scale H = 8.40×10¹⁰ GeV [C].**
   - Fragmentation comes later (f_coh < 0.1 at 0.74 / 0.63) and is milder (f_grad max 0.076 / 0.168).
   - −v fraction at most 0.5% / 1.5%, with no walls.
   - f_rad(3) = 0.67 / 0.29, K = 3.4 / 5.4.
8. **HFb vs HFa.** At the design point HFb's waterfall is the same as HFa's to 0.1% in every fraction (B_P0 vs P0). The
   two differ only in the spectator mass mechanism.

## 1. What the code is

`code/lat3d.c` (C, OpenMP, ≤ 2 threads; built into scratch at `.../wt2/e2e_lattice/build/lat3d`) evolves:

- **Fields.** σ (the clock) and the real biased χ, with the mf_hybrid HFa/HFb potential at the design scale and
  mf_chain's rescaling rule, including the one-loop term for σ > σ_c. Leapfrog in π = a³φ̇ with exponential local
  friction; 7-point Laplacian; a(t) self-consistent from the Friedmann constraint with a second-order predictor. The
  energy budget ρa³ + ∫p d(a³) is conserved to |Δ| ≤ 0.6% (P0) and ≤ 1% (all runs).
- **Spectators.** Two lattice spectators in every run:
  - s_a with m² = (g_s χ(x))², the HFa mechanism (g_s v = 3 H);
  - s_b with m² = 9 − (3/4)R(x), the HFb mechanism, with R(x) = −T(x)/M_P² taken locally. This is a closure [A]: no
    metric perturbations.
- **Axial-fermion drain.** Four axial singlet species, as a local transfer term Γ_j(x)φ̇². The ingredients:
  - mf_hybrid's half-period production rates ε(A) and the mf_hybrid skeptic's Pauli saturation F_sat(A), both through
    mf_chain's `hist_build.channels`;
  - evaluated at the local condensate energy E_σ(x), E_χ(x) (local-density closure [A]);
  - limited by a blocking energy B_j against the capacity F_sat·E.
- **Fluids.** N energies and numbers, N → LH decay and SM radiation are box-averaged fluids [A].
- **Tensor perturbations.** Sourced by ∇φ∇φ (the Dufaux et al. method) and TT-projected at read-out with the
  central-difference lattice momentum (the Figueroa, García-Bellido and Rajantie prescription).
- **Initial state.** Homogeneous slow-roll attractor (the hist_build.stage1 recipe), from t₀ = −0.02 (σ just above
  σ_c), plus a Gaussian Wigner vacuum for |k| < k_cut = 1.2 m_χ in σ, χ, s_a and s_b.

**Closures (the brackets) [A].**
| tag suffix | closure |
|---|---|
| `P0` | Pauli, local blocking B_j(x) against the local capacity; Γ_th = 0; perturbative tails free |
| `…GM` | **mf_hybrid skeptic's saturation**: capacity = F_sat(⟨E⟩)·⟨E⟩ on the box-mean condensate energy, blocking box-global (the homogeneous F_sat(A) model) |
| `P10…` | as above with Γ_th = 10 H (band unblocking by thermalisation) |
| `TR` | unblocked tracking (mf_hybrid's original friction): control |
| `HU` | no drain (Hubble friction only): control |

The physical choice between local and box-mean blocking depends on the fermion wavelengths compared with the fragment
size. It is not resolved here. The two are the bracket.

**Reused read-only.** `tracks/mf_chain/code/hist_build.py` is imported unchanged (`model`, `V2`, `grad2`, `hess2`,
`valley`, `M2c`, `channels`, `stage1` as the regression target, and `stage2` for the continuation). Its scratch is
redirected to this track's scratch via `MFC_SCRATCH`. Nothing else is copied apart from the slow-roll start recipe,
re-typed in `latcommon.homogeneous_start`.

**Files.**
- `code/lat3d.c`: the lattice.
- `code/latcommon.py`: parameters, tables, initial state, runner.
- `code/latpost.py`: percolation (periodic union-find), spectra, GW TT projection, redshift.
- `code/jobs.py`: job registry; resumable within 560 s chunks with restart files.
- `code/regress.py`: homogeneous regression.
- `code/analyze.py`: continuation, GW today, walls, convergence, exports.
- Results: `results/regress_homog.json`, `results/runs/<TAG>.json` (thinned series, snapshots, percolation, spectra),
  `results/analysis.json`, `results/hist_HFa.json`, `results/hist_HFb.json`.
- Raw fields and restart files stay in scratch.

**Regression R1 [C].** With fluctuations off (4³, all sites equal) the lattice reproduces mf_chain's homogeneous
stage1 for all three closures. The agreement is to ≤ 3×10⁻⁴ in f_rad and in H at t = 0.3 and 0.6, and to ≤ 1% in
χ/v. The time step is converged: 2×10⁻⁴ / 10⁻⁴ / 5×10⁻⁵ agree to 10⁻³ (`results/regress_homog.json`).

**Cost.** 64³ to t = 3 takes about 8 min on 2 threads. 96³ to t = 1.5 takes about 10–18 min.

## 2. Run table (64³, L = 0.42 H⁻¹ comoving, dx = L/64, dt = dx/25, k_cut = 1.2 m_χ, seed 1)

| run | t(ε=1) | t(f_coh<0.1) | f_grad max (t) | f_−v max | f_rad 1.5 / 3 | w_cyc 3 | ⟨R⟩_cyc 3 | R_rms 1.5 / 3 | m_s² rel. spread 1.5 | K (cont.) |
|---|---|---|---|---|---|---|---|---|---|---|
| P0 | 0.194 | 0.494 | 0.148 (0.52) | 0.0178 | 0.613 / 0.604 | 0.207 | 0.032 | 0.40 / 0.10 | 0.107 | 3.59 |
| P0GM | 0.194 | 0.467 | 0.221 (0.49) | 0.0305 | 0.190 / 0.248 | 0.091 | 0.073 | 1.21 / 0.35 | 0.155 | 5.65 |
| P10 | 0.194 | 0.593 | 0.092 (0.57) | 0.0083 | 0.740 / 0.745 | 0.251 | 0.020 | 0.20 / 0.05 | 0.085 | 2.83 |
| P10GM | 0.194 | 0.515 | 0.140 (0.54) | 0.0160 | 0.446 / 0.491 | 0.171 | 0.044 | 0.86 / 0.19 | 0.121 | 3.87 |
| TR | 0.194 | never | 0.000 | 0 | 0.846 / 0.857 | 0.286 | 0.011 | 0.03 / 0.01 | 0.011 | 2.05 |
| HU | 0.175 | 0.347 | 0.260 (0.37) | 0.0568 | 0 / 0 | 0.009 (t=1.5) | 0.33 (t=1.5) | 1.29 / – | 0.166 | – |
| B_P0 (HFb) | 0.194 | 0.494 | 0.149 (0.52) | 0.0178 | 0.613 / 0.604 | 0.207 | 0.032 | 0.41 / 0.10 | 0.034 (HFb m_s²) | 3.59 |
| B_P0GM (HFb) | 0.194 | 0.467 | 0.221 (0.49) | 0.0303 | 0.190 / 0.248 | 0.091 | 0.072 | 1.22 / 0.36 | 0.102 | 5.65 |
| H84_P0 (8.4×10¹⁰) | 0.207 | 0.738 | 0.076 (0.74) | 0.0047 | 0.685 / 0.671 | 0.228 | 0.026 | 0.31 / 0.08 | 0.097 | 3.43 |
| H84_P0GM | 0.207 | 0.631 | 0.168 (0.63) | 0.0152 | 0.253 / 0.290 | 0.102 | 0.068 | 1.27 / 0.38 | 0.144 | 5.39 |

Notes on the table:
- R is in H_d² units.
- K (cont.) is the HB.stage2 continuation switched at the lattice end, t = 3 [A]. Earlier switches give
  K(P0) = 3.97 / 3.75 / 3.66 / 3.59 and K(P0GM) = 6.12 / 5.87 / 5.77 / 5.65 at t_sw = 1 / 1.5 / 2 / 3. The 8% drift
  is the closure error of matching a fragmented state to a homogeneous fluid.
- TR does not fragment, because the unblocked drain empties χ first. It is a control, not a closure candidate (it
  ignores Pauli blocking).

## 3. Walls in detail [C, closures A]

- **−v fraction and percolation.** The fraction is measured on the instantaneous χ and on the 0.05 H⁻¹ running mean,
  with 6-connected periodic clusters and a union-find with box-shift bookkeeping for wrapping.
  - Maxima: 0.018 (P0), 0.031 (P0GM), 0.057 (HU, no drain), 0.005 (H84_P0).
  - `neg_wraps_ever` = false in every run. The largest cluster is ≤ 9×10⁻⁴ of the box (HU), and ≤ 1×10⁻⁴ in the Pauli
    runs.
  - The time-averaged field never has a −v region (≤ 1×10⁻⁴ in HU).
- **Area density.** Area per comoving volume from link sign changes, (2/3)n_cross/(V dx).
  - For the instantaneous field it falls from ~1 H at t = 1 to 0 by t ≈ 2.5 (GM, HU).
  - For the time-averaged field it is 0 at all times in P0 and P10, and ≤ 0.05 H in HU.
  - Late −v fraction in P0GM: 2.2×10⁻³ (t = 1), 4.5×10⁻⁴ (1.5), 1.0×10⁻⁴ (2), 1.9×10⁻⁵ (2.5), 0 (3).
- **Bias pressure vs tension [D from the potential].**
  - ΔV/v² = 2c_b/v = 0.545 (ε_bias = ΔV/(λv⁴/4) = 7×10⁻⁵); σ_w/v² = (2√2/3)√λ v = 166.7; wall thickness
    √2/(√λ v) = 0.008 H⁻¹.
  - ℓ_bias = σ_w/ΔV = 306 H⁻¹.
  - A scaling network (one wall per Hubble volume, the extrapolation [E]) would reach ρ_wall = 3M_P²H² at
    H = σ_w/(3M_P²) = 0.021 (t ≈ 23), before the bias kills it at H ≈ ΔV/σ_w = 0.0033 (t ≈ 150). Hence
    "dominates first" if a network formed. The lattice shows none forms, so no extrapolation of a network is needed.
    Isolated sub-Hubble −v regions collapse under tension and pressure together, in a time ≲ their size.
- **Beyond the box [E].**
  - The box is 0.42 H⁻¹, so the super-Hubble question is analytic.
  - At t₀ the bias displacement is x₀/v = 8.7×10⁻⁴ = 2.0×10³ H/2π (HFb 2.2×10³, H84 3.7×10³).
  - χ has |M_χ²| < (3H/2)² only for ΔN ≈ 4×10⁻⁴ e-folds above σ_c, so the super-Hubble χ variance is ≪ (H/2π)². A
    Hubble-size −v domain needs a ~2000σ fluctuation and is absent from the observable universe.
  - The UV vacuum gives box-mean/rms = 9.4 at t₀, and this ratio only grows afterwards. That is a per-site −v
    probability of ~10⁻²⁰ at the start: rare tiny seeds exist, and they collapse.
  - Caveat: this depends on k_cut. k_cut = 0.8 / 1.6 m_χ was not run (budget), and the −v overshoot fraction may scale
    with the initial noise.
- **Resolution caveat.** The physical lattice spacing grows to 0.022 H⁻¹ at a = 3.4, larger than the wall thickness.
  Walls formed late would be under-resolved. None form while the spacing is still finer than about 2× the thickness
  (t ≲ 1).

## 4. Convergence [C]

Measured at t = 1.5 against P0 (relative change in brackets).

| check | change | f_rad | f_grad | w_cyc | ⟨R⟩_cyc | R_rms | m_s² spread | f_−v max | K(t_sw = 1.5) | ρ_GW a⁴ (t = 1) |
|---|---|---|---|---|---|---|---|---|---|---|
| time step | dt/2 (`C_dt2`) | −0.1% | +0.5% | +0.3% | +0.7% | +0.2% | +0.2% | +4% | +0.1% | – |
| box | L × 1.5 at fixed dx (96³, `C_N96L63`) | −0.1% | +0.1% | +0.3% | +1.0% | +0.5% | +0.1% | +1% | +0.2% | +0.5% |
| spacing, finer | dx/1.5 (96³, `C_N96`) | +1.5% | +2.3% | +2.7% | −3.6% | −10% | +7% | −13% | +4.5% | −50% |
| spacing, coarser | dx × 1.5 (64³ at L × 1.5, `C_L63`) | −21% | +4% | −21% | +44% | +22% | +23% | +12% | +21% | −84% |

Conclusions from the table:
- Time step and box size are converged. The 1.5× box is also an independent realisation, so realisation scatter of
  the box means is ≲ 1%. Separate seeds (`C_s2`, `C_s3`) were not run.
- Spacing convergence is uneven: 1.5× coarser loses 20% of f_rad, while 1.5× finer changes it by 1.5%. The base
  spacing lies just above the resolution threshold of the fragment momenta (comoving k ≈ 200–500 vs k_Nyq = 479).
- Error estimates at the base spacing: about 2–5% on fractions, w and K; about 10% on R_rms, the m_s² spread and the
  −v fraction.
- The GW amplitude is **not** converged in spacing. The spectrum peaks at the lattice Nyquist in every run, with
  50–90% of the TT power above k_Nyq,1D. Its dx dependence is non-monotone (×0.5 finer, ×0.16 coarser). Treat it as a
  cutoff-dependent estimate, good to a factor of a few to ten.

## 5. Gravitational waves [C lattice, E redshift]

- **Method.** Tensor perturbations sourced by the TT part of ∂φ∂φ, as in Dufaux et al. (arXiv:0707.0875); TT projector
  with the lattice momentum k_eff = sin(k dx)/dx, as discussed by Figueroa, García-Bellido and Rajantie
  (arXiv:1110.0337). The tachyonic-preheating context: Felder et al. (arXiv:hep-ph/0012142), García-Bellido & Figueroa
  (arXiv:astro-ph/0701014), Dufaux, Felder, Kofman & Navros (arXiv:0812.2917). All five arXiv numbers were verified by
  search this session.
- **Redshift.** ρ_GW a⁴ is frozen after t ≈ 1. Then Ω_GW,RD = (dρ_GW/d ln k)a⁴/(3M_P²K²), with K = a²H in the late
  radiation era. Today it is ×Ω_γh²(g_*/2)(g_s0/g_*)^{4/3} = 1.63×10⁻⁵ at g_* = 106.75. The frequency is
  f = k/(2πa₀), with a₀ from T(a = 1) = (90K²M_P²H²/(π²g_*))^{1/4} and entropy conservation [E: standard thermal
  history after reheating; g_* = 106.75 is an assumption].
- **Values (total Ω_GW h² today, peak f).**
  - P0: 3.3×10⁻¹¹ at 64³ and 1.4×10⁻¹¹ at 96³ (K = 3.6–3.9); 1.4×10⁻¹² / 0.7×10⁻¹² at K = 17.7. Peak 2.1 GHz (64³),
    2.9 GHz (96³, again at its Nyquist), 0.95–1.4 GHz at K = 17.7.
  - P0GM: 5.1×10⁻¹¹ (K = 5.65), 5.2×10⁻¹² (K = 17.7).
  - P10: 2.8×10⁻¹¹. H84_P0: 0.9×10⁻¹¹. HU (no drain): 5.8×10⁻¹² at K = 17.7. TR: 10⁻¹⁹.
  - χ-scale frequency (k at the χ power peak): 0.4–1.6 GHz.
- **Retrospective comparison [R].**
  - The ΔN_eff ceiling on a GW background, Ω_GW h² ≲ 5.6×10⁻⁶ ΔN_eff (see the review by Caprini & Figueroa,
    arXiv:1801.04268), is ≈ 1.7×10⁻⁶ for ΔN_eff < 0.3. The model sits ≥ 4.5 orders below it, even allowing a factor
    10 for the unconverged amplitude.
  - The GHz band is outside every operating detector.
  - For the domain-wall GW context (none here), see Saikawa's review (arXiv:1703.02576, verified).

## 6. Exports for e2e_chain

**`results/hist_HFa.json`**, primary P0 (local Pauli), with alternates P0GM (the skeptic-saturation bracket), P10,
P10GM, H84_P0 and H84_P0GM.

`lattice` block, on a time grid of 0.004 for t < 1 and 0.02 after, up to t = 3:
- t, a, H, eps_H, w (raw), w_cyc;
- R (raw box mean), R_cyc, R_rms (spatial);
- f_sig, f_chi, f_grad, f_rad, f_SM, f_s;
- s_mean, s_std (s/s_*), m_s2_mean, m_s2_std;
- fminus, area, K_running.

`fluid` block (primary only): the HB.stage2 continuation to K convergence [A], giving K for this closure. Its R_cyc is
the matter part (ρ_s + ρ_c)/M_P². The alternates carry their own K_cont, and every run carries its `head`, with
headline numbers and GW today.

**`results/hist_HFb.json`**: the same, with primary B_P0 and alternate B_P0GM. Here s and m_s² are the HFb spectator
s_b with 9 − 3R/4.

**Splicing advice [A].** Use the lattice block as the coherent stage up to t_sw ∈ [1.5, 3], where f_coh < 0.01 and
the box means are smooth. Continue with e2e_kinetic's `fluid_HFa.json`, matched by moments, rather than with this
track's stage2 continuation, which carries no SM thermalisation physics.

Retrospective check: e2e_kinetic's f_rad = 0.31 at t = 1 and 0.22 at t = 3 sits inside this track's bracket
[0.17–0.25 (GM), 0.60 (local)], nearer the box-mean closure.

**Concurrency note.** e2e_kinetic's `fluid_HFa.json` was being rewritten during this pass. It gave K = 17.70 when
`analyze.py` read it, and K = 17.20 at 01:33. Every "K = 17.7" GW number here scales as Ω ∝ K⁻² and f ∝ K⁻¹ᐟ²
(+6% and +1.5% at 17.2). The chain should re-evaluate with the final K using `latpost.today`.

## 7. Open, not done (budget)

- **k_cut scan** (`C_kc08`, `C_kc16`) and separate seeds (`C_s2`, `C_s3`). These would firm up the −v overshoot
  fraction and its seed dependence. The 1.5× box run is the only independent realisation.
- **Remaining runs:**
  - Γ_th = 10 at H = 5×10¹⁰ GeV (`H50_P10`, `H50_P10GM`), mf_critic's EFT-edge point;
  - H = 2.4×10¹¹ (`H240_HU`);
  - long runs to t = 8 (`P0_long`, `P0GM_long`).
- **A 96³ GW run at a larger box with a finer spacing**, to separate physical and cutoff GW power.
- **A first-principles criterion for local vs box-mean Pauli blocking.** This needs the fermion momentum (~μ₅)
  compared with the fragment scale (comoving k ≈ 200–500). It is the main uncertainty in K: 3.6 vs 5.7 at Γ_th = 0.
- **Missing physics:** metric perturbations, back-reaction of s on χ (negligible: f_s ≤ 3×10⁻⁵), and spatial
  transport of the N fluids.

## 8. Literature (arXiv numbers verified by web search this session)

- J.-F. Dufaux, A. Bergman, G. Felder, L. Kofman, J.-P. Uzan, *Theory and numerics of gravitational waves from
  preheating after inflation*, PRD 76, 123517 (2007), arXiv:0707.0875.
- J.-F. Dufaux, G. Felder, L. Kofman, O. Navros, *Gravity waves from tachyonic preheating after hybrid inflation*,
  arXiv:0812.2917.
- J. García-Bellido, D. G. Figueroa, *A stochastic background of gravitational waves from hybrid preheating*,
  arXiv:astro-ph/0701014.
- D. G. Figueroa, J. García-Bellido, A. Rajantie, *On the transverse-traceless projection in lattice simulations of
  gravitational wave production*, JCAP 11 (2011) 015, arXiv:1110.0337.
- G. Felder, J. García-Bellido, P. B. Greene, L. Kofman, A. Linde, I. Tkachev, *Dynamics of symmetry breaking and
  tachyonic preheating*, PRL 87, 011601 (2001), arXiv:hep-ph/0012142.
- C. Caprini, D. G. Figueroa, *Cosmological backgrounds of gravitational waves*, CQG 35, 163001 (2018),
  arXiv:1801.04268. The ΔN_eff relation is quoted in its standard form; the review's own section was not opened.
- K. Saikawa, *A review of gravitational waves from cosmic domain walls*, Universe 3, 40 (2017), arXiv:1703.02576.
