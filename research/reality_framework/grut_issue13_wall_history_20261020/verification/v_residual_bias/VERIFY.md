# Verification v_residual_bias: module residual_bias

21 October 2026. Adversarial verifier. I wrote all the code in this directory myself and imported nothing from the
module. I read the module code (`common.py`, `a_damping.py`, `b_inflation.py`, `d_charged.py`, `e_synthesis.py`)
only to learn its definitions. Where possible my numerics differ on purpose.

- numpy/scipy only; at most 3 jobs ran at once, each on 1–2 threads.
- Longest job: `vrb_width.py`, 482 s, almost all of it the classical power-law tail. Every other job took under 7 min.
- No git was used. Nothing outside this directory was written.
- Read-only inputs [Q]:
  - `modules/thermal_potential/results/kibble_zurek.json` (ξ_KZ and ε_freeze against ξ);
  - `grut_candidate_1_quantum_20261018/rg_beta_functions.txt` (betas, re-typed);
  - the module's result files, for comparison only;
  - the derived thresholds of `modules/bias_threshold` (ε_persist = 6e-24, ε_safe = 1.5e-18), for a sensitivity
    check.

| File → output | What it does |
|---|---|
| `vrb_width.py` → `results/width.json` | Zero-momentum width of the soft φ mode. The (E_b, E_c) parametrisation differs from the module's (E_c, E_d). The statistical factor is f_b(1+f_c)(1+f_d) − (1+f_b)f_c f_d, not the KMS-simplified form. Vertices and thermal masses are re-derived from U. The classical tail is integrated to ∞. |
| `vrb_damping.py` → `results/damping.json`, `damping.log` | D = ∫(Γ/2)dt by quad over ln T. c_w by quadrature. ξ_min by brentq. The Hubble-only amplitude from the **exact anharmonic adiabatic invariant**. The dilute gas by direct quadrature. The route C noise. |
| `vrb_trapping.py` → `results/trapping.json`, `coherent_mean.json`, `trapping.log` | Route B. Own RG. **Spectral** Fokker–Planck: Schrödinger transform, eigh_tridiagonal on a node grid with y0 = v/H on a node. Own vectorised RK4 trapping map in conformal time, with dτ-halving and h-doubling checks. Own SqRA emulation of inflation_spectator's h = 0.02 cell-centre method. |
| `vrb_patch.py` → `results/patch.json` | Route B patch spread: equilibrium-draw rms of m and one-e-fold sub-patch sd |
| `vrb_lattice.py` → `results/lattice/*.json` | Route D, Minkowski, 2D, 384², dx = 0.5, velocity Verlet. Own code, own seeds: κ = 1.5 and 2 (3 seeds each), κ = 4 (2 seeds). |
| `vrb_quartic.py` → `results/quartic.json` | Test of the Hubble-only claim. Does a coherent λφ⁴ oscillation (2D, 192², noise 1e-3 and 1e-6) stay coherent for 400 periods? |

## Bottom line

**Every headline number of the module that I recomputed is reproduced independently.**
- The soft-mode width Γ/T = 3.8424e-5 matches to all printed digits. So do the classical 4.428e-5, nσv 3.364e-5, the
  channel split and 6.42e-4 at t_G.
- Parwani's coefficient is reproduced as an exact identity: the classical one-field integral equals π²/6 to 1e-9.
- D(H★) = 2990, against the module's 2992.
- All ξ_min tables agree to 3 digits.
- Route B: every m0 agrees to ≤ 0.4%, using a different (spectral) FP method and my own trapping map. So do the
  asymptotic laws and the patch spreads.
- The inflation_spectator p_− values 0.49993 / 0.49969 / 0.49970 / 0.49967 are reproduced to 1e-12 by an h = 0.02
  cell-centre emulation. The converged values are 0.50109 / 0.50018 (N50) and 0.50076 / 0.50015 (N60).
- The route D lattice reproduces the module at κ = 1.5 and κ = 4. At κ = 2 the sign and the O(0.5) size agree.

**The physics behind route A holds and is robust.**
- The width formula is the correct net (Weldon) rate. The amplitude decays at Γ/2.
- 71% of D comes from the constant high-T rate, and only 5% from ε < 1e-2. The near-critical region is therefore not
  load-bearing.
- The classification against ε_persist needs only D ≈ 51 (1e-20) to 58 (6e-24). That is about 2% of the computed
  D, i.e. the rate could be about 50× below the central one. So the conclusion needs exponential decay over about
  50–60 e-folds, not 3000.
- The value log10|m0| = −1297 is formally right but has no physical content beyond "far below the realization
  floor": that floor is about 1e-33 at formation and about 1e-25 at parity (bias_threshold).

**Corrections** (none changes the H★ classification):
1. **ξ_kin coincidence.** "Any self-consistent (kinetically thermal) bath lands below ε_persist" holds only at the
   central rate and the 1e-20 placeholder, for production at f ≳ 32.
   - At ξ_kin = 6.6e-5: log10|m0| = −21.2 (f = 32), −19.6 (f = 10), −9.9 (f = 32, Γ × ½).
   - With bias_threshold's derived ε_persist = 6e-24, ξ_min = 6.9e-5 (f = 32) and 7.1e-5 (f = 10), both above ξ_kin.
     With Γ × ½, ξ_min = 9.1e-5 to 9.8e-5.
   - A strip 6.6e-5 ≲ ξ ≲ 7e-5 to 1e-4 of kinetically thermal baths is therefore not erased below ε_persist.
   - Separately, "μ_N = 0" requires chemical equilibrium (2→4), which thermal_potential places at ξ ≳ 0.068 (band
     0.0084–0.27). The rate per unit Π is nearly μ_N-insensitive (dilute formula 3432 against thermal 2990 at the same
     A), so this matters less.
2. **Hubble-only |m0| = 1 is an overclaim.** r = 34 (my exact-invariant value is 37) is evaluated at the KZ point.
   - For a 100 GeV quartic oscillation the decision comes only at trapping, about 1040 oscillations after the onset
     of nonlinearity.
   - Amplitude-dependent dephasing (176 rad for δA/A = σ/A) and Lamé self-resonance (about 280 e-folds) act in that
     time.
   - My lattice: the coherent mode fragments after 35–60 periods. After 400 periods it keeps 0.43–0.46 of its
     amplitude, against a fluctuation rms of 0.23–0.29.
   - The expected outcome is |m0| = O(0.5–0.9), with a **phase-set** sign. "Memory not erased" holds; "|m0| = 1, one
     sign everywhere, the old vacuum kept" does not. The same applies to "a cold +v start keeps its sign" under the
     free-streaming restoration (D = 3.3 to 1e-3).
3. **The route B → free-streaming band (§4.6) is not a valid band.**
   - Its lower edge applies the trapping-map law 0.36(v/H)³ at H_inf = 1e5–1e8 GeV. That is 4–7 decades beyond the
     computed range (≤ 2e4 GeV), and inside the regime the module itself says is preheating-restored.
   - Its upper edge multiplies ⟨φ⟩/φ_rms by e^{−D}, although the random long-wavelength part is damped by the same
     factor. The relevant ratio is ⟨φ⟩e^{−D}/max(φ_rms e^{−D}, σ_dec), which is ≥ the module's upper edge.
   - So "exit at 2.6e7 GeV gives < 1e-23" and the lowest edge −16.9 are not supported. The conclusion "H★′ is not
     shown to be unbiased" stands and becomes stronger.
4. Minor points:
   - The asymptotic threshold spacing is 2.6–2.8 v/H, not v/(0.313H) = 3.2 v/H; the latter holds only for the first
     thresholds. It does not enter any m0.
   - In route D, "ω = 1.03–1.10 m_X, period 6 ℓ" is the prepared X orbital frequency. The coherent ⟨S⟩ swings have
     periods of about 10–11 ℓ.
   - At κ = 2 my 3 seeds give m0 = −0.61 ± 0.06, against the module's −0.44 ± 0.11 (4 seeds). The combined 7 seeds
     give −0.51 ± 0.12, so the module's ± understates the scatter.
   - Parwani has an erratum (PRD 48, 5965 (1993)) that the module does not cite. Nothing depends on it, because the
     coefficient is reproduced here from first principles.

## Claims, method, my value, verdict

### Route A: soft-mode damping

| # | Claim (module) | Method (mine) | My value | Verdict | Reason |
|---|---|---|---|---|---|
| A1 | Γ = Σ mult\|M\|²S/(64π³m_a)∫dE dE F Θ; F = net Bose rate; amplitude decays at Γ/2; vertices 6λ, 2g², 2h²; thermal masses Π_φ = λ/4+g²/6+h²/12, Π_q = g²/12+λ_X/3, Π_D = h²/12+λ_D/4 | own derivation from the two-body phase space, eliminating p_d instead of p_b; vertices ∂⁴U; masses from ⟨χ²⟩ = T²/12 | same prefactor 1/(64π³m_a), same vertices and masses | CONFIRMED | Weldon net rate f_b(1+f_c)(1+f_d) − (1+f_b)f_cf_d, Im Π(m) = −mΓ, so the amplitude decays as e^{−Γt/2} |
| A2 | Classical one-field integral = π²/(6m) to 4e-13, i.e. Parwani γ = g⁴T²/(1536πm), Γ = 2γ | my integral ∫F_cl over (E_b, E_c), tail to ∞, m = 0.003–1 | π²/6 to ≤ 6e-9 | CONFIRMED | Read level: arXiv abstract page of hep-ph/9204216 only. It confirms title, journal and an erratum PRD 48, 5965 (1993) that the module does not cite. The coefficient follows here from the identity itself. |
| A3 | Pure φ⁴ full Bose / LO = 0.878 | same | 0.8782 | CONFIRMED | — |
| A4 | Γ/T_h (T ≫ T_c) = 3.84e-5 full, 4.43e-5 classical, 3.36e-5 nσv; channels 1.97e-5 / 1.26e-5 / 6.1e-6 / ~1e-9 | own integrator | 3.8424e-5, 4.4279e-5, 3.3637e-5; 1.9686e-5 / 1.2604e-5 / 6.133e-6 / 8.2e-10, 3.3e-10 | CONFIRMED | Differences < 1e-9. (My first classical run cut the outer integral at E = 60T and was 0.9% low. Fixed: f = T/E has a power-law tail.) |
| A5 | Rises as 1/m near T_c: 6.4e-4 T at t_G; Γ/m = 0.27; underdamped | own table | 6.418e-4, Γ/m = 0.267 | CONFIRMED | — |
| A6 | 77× the hard rate 4.98e-7 T | I re-derived the hard rate as n⟨σv⟩ = (18λ²+12g⁴+6h⁴) ζ(3)T/π² /(16π·14.59) | 4.98e-7 T; ratio 77.1 | CONFIRMED | The soft rate exceeds the hard rate because s ≈ 2m_φE_b ≪ 2E² and the soft Bose occupancy is enhanced |
| A7 | 1↔2: φ→qq closed (m_φ/2m_q ≤ 0.946); φ→DD open for T ≳ 3T_c only via the background vertex, 3.4e-6 of the linear rate | arithmetic | 0.9459; DD threshold between 2T_c and 3T_c; ratio 3.4e-6 | CONFIRMED | Second-order effect, correctly not credited |
| A8 | H★: A = 7.3e8 GeV, T_exit = 3.115e6 GeV, H_exit = 1.36e-5 GeV, Γ/H(ε = 7e-3) = 2.5e4 | arithmetic, own Γ | λv²/H_exit = 7.34e8 (Π_φξ²M_P/√(π²g*/90) = 7.43e8 with the LO T_c); 1.363e-5; 2.47e4 | CONFIRMED | The 1% difference in A comes from T_c = 1557.5 vs the LO 1547.6 |
| A9 | D = ∫(Γ/2)dt from T_SM = 1e8 GeV to t_G = 2992 (classical 3336, nσv 2959); hard-rate floor 27.5 | quad over ln T with own Γ(ε), 47 points | 2990.0 (classical 3357 from a 14-point table, nσv 2957); floor 27.55 | CONFIRMED | Constant high-T rate alone: 2127. Cut at ε = 1e-2: 2842. Robust. |
| A10 | Coherent remnant → m0 = erf(r/√2), σ² = c_w T_c m̂, c_w = 0.01087, σ = 2.9 GeV, A_lin = 446 GeV, r_undamped = 155 | c_w by quadrature; same KZ input [Q] | c_w = 0.0108708; σ = 2.885; A_lin = 445.7 GeV; r = 154.5 | CONFIRMED | The Gaussian-window σ and the WKB a^{−3/2}m^{−1/2} are standard. Order-one choices shift log10 m0 by ≲ 1 only. |
| A11 | log10\|m0\| = −1297 (band −648 to −2597); classical −1447, nσv −1283; hard floor −9.9 | own D | −1296.5 (−647 to −2595); −1456; −1282; −9.87 | CONFIRMED | The 0.8-decade difference = ΔD 1.9 (0.06%). Physically the number means only "≪ realization floor". |
| A12 | Hubble only: amplitude turns quartic at 1565 GeV, 99 GeV at decision, r = 34.4, **\|m0\| = 1, one sign everywhere** | exact anharmonic adiabatic invariant; oscillation count to trapping; 2D quartic lattice | onset 1564 GeV, A = 107 GeV, r = 37; about 1040 oscillations from onset to trapping (ε = 5.6e-3 below T_c); dephasing 176 rad; resonance about 280 e-folds; lattice coherent envelope 0.43–0.46 vs rms 0.23–0.29 after 400 periods | **CORRECTED** | r at the KZ point is not the decision ratio for a large coherent oscillation. Fragmentation and dephasing give \|m0\| = O(0.5–0.9) with a sign set by the oscillation phase [E]. "Memory not erased by Hubble friction" stands. |
| A13 | production_map WKB at H★: ratio 2.1e-5, m0 ≈ 1.7e-5; not applicable (assumes the bath at H > A) | arithmetic; the zero mode is frozen iff m_φ² = AH − λv² < H², i.e. H ≳ A | 2.107e-5; H_prod = 1.4e-2 GeV ≪ A | CONFIRMED | — |
| A14 | ξ_min for log10\|m0\| < −17/−20/−30: f = 32: 5.9/6.4/7.8e-5; f = 10: 6.2/6.7/8.1e-5; f = 2: 7.7/8.3/9.9e-5; f = 1.1: 1.35/1.44/1.73e-4; f = 1.01: 2.8/3.0/3.6e-4 | brentq on my D | 5.95/6.42/7.78e-5; 6.19/6.67/8.05e-5; 7.69/8.26/9.92e-5; 1.347/1.445/1.733e-4; 2.79/2.99/3.59e-4 | CONFIRMED | — |
| A15 | Scan values, e.g. ξ = 1e-4: −50 / −8.2 / −0.07; 2e-4: −206 / −41 / −7.4; 1e-3: −5195 / −1096 / −254 | same | −50.4 / −8.16 / −0.073; −206.0 / −40.9 / −7.42; −5192 / −1095 / −253.7 | CONFIRMED | — |
| A16 | "Erasure threshold ≈ 6.4e-5 coincides with ξ_kin = 6.6e-5 → any self-consistent thermal bath lands below ε_persist" | sensitivity to Γ band, f and the derived ε_persist | log10\|m0\|(ξ_kin) = −21.2 (f = 32), −19.6 (f = 10), −9.9 (Γ×½); ξ_min(6e-24) = 6.9e-5 (f = 32), 7.1e-5 (f = 10); Γ×½: 9.1e-5 to 9.8e-5 | **CORRECTED** | It holds only at the central rate, the 1e-20 placeholder and f ≳ 32. The margin is a factor 1.0–1.5 in ξ, not a safety factor. A kinetically thermal strip 6.6e-5 ≲ ξ ≲ 7e-5 to 1e-4 is not covered. |
| A17 | Dilute gas: C = 7.46e-6, c_eff = 0.012525, D = 4.68e-6 (A/GeV); A ≳ 1.0e7 GeV and T_SM,end ≲ 2.6e7 GeV for 1e-20; D(same A as H★) = 3432 | re-derived Γ_d = Σ mult\|M\|²S J_b/(32π m) and Π = Σ coupling × J; D by direct quadrature, not the arctan form | 7.4605e-6, 0.012525; 4.6782e-6 (numeric = formula to 1e-13); A = 9.8e6 to 1.08e7, T_end = 2.70e7 to 2.57e7; 3431.7 | CONFIRMED | — |
| A18 | Heavy-inflaton producer (exits 1e8–6e9 GeV): D = 3.3 to 1e-3, "memory essentially kept; a cold +v start keeps its sign" | arithmetic plus A12 | D = 3.33 to 9.2e-4 | CONFIRMED for D; **CORRECTED** for the wording | The bias stays O(1), but the sign is phase-set, not "kept" (A12) |

### Route B: inflationary prior and trapping

| # | Claim | Method | My value | Verdict | Reason |
|---|---|---|---|---|---|
| B1 | RG at 10 TeV: ξ_φ = −4.889e-4, a2 = −5.97e-3, y_min = 0.770 | own RG ODE; ξ betas derived as the mirror of the mass betas | −4.8889e-4, −5.967e-3, 0.7700 | CONFIRMED | — |
| B2 | Trapping map odd ⇒ m0 = ∫₀^∞[P(a) − P(−a)]s(a)da | Z2 of χ'' + λχ³ − λu²τ²χ = 0 with odd initial data | identity holds | CONFIRMED | — |
| B3 | First threshold a₁ = 2.953 v/H; spacing ≈ v/(0.313H) | own RK4 map, vectorised bisection; dτ halving moves thresholds ≤ 9e-4 | a₁ = 2.950–2.953 u (5–20 TeV); spacing 3.26u at the first gap, **2.6–2.8u asymptotically** | CONFIRMED (a₁); minor CORRECTION (spacing) | For y ≲ 14, f0 ≈ (2/λ)^{1/4}√y gives N_cross ≈ √2 y/(3.708u) = 0.381 y/u. No m0 depends on the description. |
| B4 | m0 after trapping, fresh N50 (N60): +0.391 (+0.346) / +0.109 (+0.080) / −0.0267 (−0.0297) / −2.18e-3 (−1.52e-3) / −3.63e-4 (−2.90e-4) / −4.51e-5 (−3.61e-5) at 1.4/2/3/5/10/20 TeV | spectral FP (node grid, y0 on a node) × own map; cubic-spline integration between exact thresholds | +0.3910 (+0.3457) / +0.1086 (+0.0797) / −0.02670 (−0.02970) / −2.1758e-3 (−1.5219e-3) / −3.645e-4 (−2.912e-4) / −4.519e-5 (−3.611e-5) | CONFIRMED | ≤ 0.4%. h-doubling and dτ-halving change m0 by ≤ 1e-5 relative. |
| B5 | Relaxed + start: +0.43 / +0.10 / −0.095 / −4.5e-3 / −2.1e-3 / −5.2e-4 | same | +0.4335 / +0.0955 / −0.0952 / −4.494e-3 / −2.080e-3 / −5.144e-4 | CONFIRMED | — |
| B6 | Asymptotic −0.36(v/H)³ (N50), −0.29 (N60); relaxed −0.206(v/H)² and −0.186 | ratio at 10 and 20 TeV | −0.3645, −0.3615; −0.2912, −0.2889; relaxed −0.208, −0.2057; −0.1873, −0.1856 | CONFIRMED | The laws are fitted for 7–20 TeV only. Their use at H_inf = 1e5–1e8 GeV is an extrapolation (see B11). |
| B7 | One sign change between 2 and 3 TeV | B4 | +0.109 → −0.027 | CONFIRMED | — |
| B8 | inflation_spectator p_− = 0.49993 / 0.49969 is a grid artifact (h = 0.02 cell-centre signs; start cells 0.11 instead of 0.10). Converged: 0.50109 / 0.50018 (N60: 0.50076 / 0.50015), on the other side of ½, systematic not statistical. | own SqRA chain with h = 0.02 cell-centre signs and my thresholds; start cells 0.21 / 0.11 | 0.4999299890642 / 0.4996914438763 (N50); 0.4996963137 / 0.4996733492 (N60): inflation_spectator's numbers to 1e-12. Start 0.09 instead of 0.11 gives 0.49975, h = 0.005 gives 0.50045 (not converged). Converged: 0.501088 / 0.500182; N60: 0.500761 / 0.500146. | CONFIRMED | The cell-centre sign assignment dominates. The start-cell offset (also 0.21 instead of 0.20 at 5 TeV) is secondary. The FP is deterministic, so the deviation is a systematic ensemble bias. Distances from ½: 0.45–15.5× the module's ✓. |
| B9 | End-of-inflation p_− (fresh N50) 0.430 / 0.465 at 5 / 10 TeV | spectral | 0.43015 / 0.46473 | CONFIRMED | — |
| B10 | Patch: rms over equilibrium draws 4.8e-3 / 2.2e-3 / 5.5e-4; sub-patch sd 1.7e-3 / 5.7e-4 / 1.4e-4 (5/10/20 TeV) | spectral backward m(y0, N), cell-averaged F | 4.80e-3 / 2.22e-3 / 5.49e-4; 1.67e-3 / 5.83e-4 / 1.45e-4 | CONFIRMED | "No realization comes near 1e-17" should read "probability ~1e-17/rms ≈ 1e-14". Minor. |
| B11 | B → free-streaming: log10\|m0\| in [−16.9, −2.2]; exit at 2.6e7 GeV gives < −23 | arithmetic of the band reproduced (e.g. H = 1e8, T_end = 1e8: −16.9 … −6.6) | band arithmetic ✓; coherent mean ⟨φ⟩/v = 0.858 / 0.894 / 1.015 and ⟨φ⟩/φ_rms = 0.085 / 8.7e-3 / 9.3e-6 at 1e4 / 1e5 / 1e8 GeV ✓ | **CORRECTED** (construction) | Lower edge: the trapping law is used 4–7 decades outside its computed range, inside the preheating-restored regime. Upper edge: e^{−D} is applied to a ratio whose denominator is damped too. The true upper value is ⟨φ⟩e^{−D}/max(φ_rms e^{−D}, σ_dec) ≥ the module's. The "< −23" is unsupported. The H★′ conclusion is unchanged or strengthened. |
| B12 | B then H★ thermal bath: log10\|m0\| ≤ −1297 + log10(⟨φ⟩/v) | A11 and B11 coherent mean | ⟨φ⟩/v ≈ 1.0 at 1e8 GeV | CONFIRMED | With fluctuation–dissipation the thermal σ is regenerated, so here e^{−D} applies to the ratio, unlike B11 |

### Routes C and D

| # | Claim | Method | My value | Verdict | Reason |
|---|---|---|---|---|---|
| C1 | Z2-invariant state ⇒ m0 = 0 exactly | symmetry argument | — | CONFIRMED | — |
| C2 | Realized white noise about 1e-33 (about 1e52 Hubble volumes × 1e14 domains) | own count with T_CMB, g_s0 = 3.909, H0 = 67.4 | 1.34e52 × 4.7e13, noise 1.26e-33 | CONFIRMED | Noise, not systematic |
| D1 | Minkowski κ = 1.5: m0 = +0.83 ± 0.07 (3 seeds; the 4th all +); minority disappears | own 2D lattice, own seeds | +0.911 / +0.834 / +0.799, i.e. +0.85 ± 0.06; late p = 1.000 | CONFIRMED | — |
| D2 | Minkowski κ = 2: m0 = −0.44 ± 0.11; late p = 0.02–0.12 | same | −0.579 / −0.681 / −0.566, i.e. −0.61 ± 0.06; late p = 0.003–0.026 | CONFIRMED (sign and O(0.5)); value differs | The combined 7 seeds give −0.51 ± 0.12. The quoted ± understates the scatter. |
| D3 | κ ≥ 3 (κ = 4): stays restored, f_core 0.90–0.93, late p = 0.500 ± 0.001; coherent remnant decays at 0.09–0.13 m_X to a 1–2% floor | same, κ = 4 | f_core 0.928 / 0.926; late p 0.4992 / 0.5004; decay 0.108 / 0.095 per ℓ; floor 0.010 | CONFIRMED | — |
| D4 | "Oscillates coherently at ω = 1.03–1.10 m_X, period about 6 ℓ" | module's stored ⟨S⟩ series and mine | ⟨S⟩ swing period about 10–11 ℓ; 1.02–1.10 is the prepared X orbital ω | minor CORRECTION | The decoherence times (30–40 ℓ) are right |
| D5 | Expanding τ_i = 150 reproduces 0.71 (0.70–0.81); τ_i = 600 gives −0.40 | only re-read the module's stored runs | p_window 0.697 / 0.730 / 0.809; −0.40 | NOT CHECKED (no independent rerun) | — |
| D6 | n_i ≳ 3n_c: remnant erased by exponential continuation exp(−Γ_coh Δt) ≈ exp(−1e8), Δt ≈ 1.4e9 ℓ | arithmetic | Δt = 1.35e9 ℓ at T_SM = 200 TeV (κ = 4); Γ needed ≥ 3.7e-8 per ℓ | arithmetic CONFIRMED; extrapolation NOT CHECKED | 8 decades beyond the 1–2% lattice floor; the module labels it [E] |

### Out of scope here

The bias_threshold items listed in the task (phase-ordering law, λ values, own fits, K, relic-safe temperature) are
verified in `verification/v_bias_threshold`. Only one bias_threshold number is used by this module: the naive 2D
fit K10 = 0.262(|m0|/2)^−1.255 from network_sim, which is quoted correctly.
- My arithmetic gives log10 K10 = 5.25 / 5.37 at 20 TeV fresh N50 / N60, matching "10^5.3–10^5.4".
- That is an extrapolation of 4 decades below the fitted range 0.05 ≤ |m0|/2 ≤ 0.3. The module labels it
  "orientation only".

## Details

**Width integral.** For a at rest with p_d = p_b − p_c and E_d = m_a + E_b − E_c, the angular δ gives
E_d/(p_b p_c). The allowed E_c interval for fixed E_b is the region between the roots of
(m_b² + 2m_aE_b + m_a²)E_c² + 2(E_b + m_a)β E_c + β² + p_b² m_c² ≤ 0, with
β = (m_d² − m_a² − m_b² − m_c²)/2 − m_aE_b. The leading coefficient is always positive, so the region is a single
interval. This is algebraically distinct from the module's quadratic in E_d, and the results agree to < 1e-9.

**Why D is robust.** Γ/H ∝ ξ/T_SM at a fixed Γ/T_h, so ∫(Γ/2H)dT/T is dominated by the lowest temperatures
of the rate-constant region:
- the constant-rate analytic value is 2127;
- the 1/m_φ rise near T_c adds about 860, of which only about 150 comes from ε < 1e-2;
- reaching ε_persist needs D ≈ 51 (1e-20) to 58 (6e-24).
A rate overestimate of 50× would be needed to change the H★ class.

**Long-time behaviour.** I agree with the module's argument. No Z2-odd density is conserved, so there is no
hydrodynamic long-time tail at linear order in ⟨φ⟩. Composite "φ(k) × hydro(−k)" channels decay at least at the
soft φ rate, because hydrodynamics applies only for k below the hard-particle inverse mean free path. This is
argument level only, not computed. As shown above, it needs to hold over about 50–60 e-folds, not 3000.

**Hubble-only counterfactual (A12).** From the exact adiabatic invariant J(A, m²)a³ = const:
- the onset of nonlinearity is at T_h = 1564 GeV, with A = 92 GeV;
- the amplitude at the KZ point is 107 GeV;
- trapping requires λA² = Π_φ(T_c² − T²), i.e. ε = 5.6e-3 below T_c;
- with the quartic frequency 9.1 GeV, that is about 1040 oscillations from onset to trapping.

In my 2D λφ⁴ lattice (A = 1, √λA = 1) the zero mode stays at 1.00 until about 200 (noise 1e-3) or 370 (noise 1e-6)
time units. It then fragments, and settles at an envelope of 0.62–0.65 (t ≈ 740–1100) and 0.43–0.46 (t ≈ 3000,
400 periods), with an rms of 0.23–0.29. The coherent part survives partially, so |m0| = O(1) is plausible, but not
1 and not with the original sign. The counterfactual has no bearing on the H★ classification.

**Inputs not re-derived [Q].** The following are quoted, not re-derived here:
- the KZ table (ξ_KZ = 2.03 GeV⁻¹, ε̂ = 2.78e-5 at ξ = 5e-4);
- T_c = 1557.5 GeV and t_G;
- ξ_kin;
- the RG betas;
- the bias_threshold thresholds.
