# Consolidated review: Codex 373d37e, the lower-scale supplied-phase point (track sp_critic)

25 October 2026 · Claude (issue #13) · consolidates tracks sp_point, sp_yledger, sp_vector and sp_phaseprep.
The Codex packet (`origin/codex/grut-inflation-scale-scan-20261025`,
`research/reality_framework/grut_inflation_scan_20261025/`) was read only through `git archive 373d37e` into the
scratchpad. Our own code is in `code/` and outputs are in `results/`. No Codex code was run and nothing was committed.
The prior review is REVIEW_GAUGEPOINT.md together with tracks/gp_critic/CONSOLIDATED.md.

**This is one supplied tree-level working model. Nothing here is evidence that X, Y, Σ, φ, the vector or the
inflationary inputs exist. It is not a validation, an exclusion or a GRUT prediction.** All five tracks are
same-session AI work, so this is not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

Stated couplings, unchanged [Q]: M_Y = 1e12, h = 1e11, μ = 5e11 GeV, κ = 1e-3, e = 1e-7, λ_Y = λ_Σ = 1.
Light target: g_X = 0.1, λ_X = 1e-3, λ_φ = 0.03, v = 0.924089637 GeV, Q = 2e5. Required yield Y_Q = 4.8294e-9;
the 1% allowance is 0.0049296 pairs per net charge (f = 0.01).
Point: H* = 1.1211167304e13 GeV, c_H = 1, aligned radii, Σ = −|Σ| e^{iπ/2}, instant release into radiation.

## Bottom line

- **Codex's point reproduces.** Four independent integrators give it to every quoted digit. Its final charge
  fraction is ε = −0.2402, and its yield is 66 Y_Q.
- **It is a narrow, low-sensitivity island, not a regular window in H_inf.** It holds only for a supplied CP-odd
  preparation that the baseline action does not produce.
- **The leftover requirement fails at this point in every limit we computed, with any preparation.** The release's
  own debris is the cause:
  - the Y_L decays;
  - the Σ radial mode decaying to two vectors;
  - the heavy-to-light Y decay plus a vector.
- **Size of the failure.**
  - Floor: about 27 times the allowance (debris pairs alone).
  - Frozen limit without erasure: about 350 times.
  - The equilibrium benchmark the rate estimates point to: about 1500 times.
- **"Clears the X-only 0.19 test" is true only inside a benchmark (maximal equilibrium) that the release's own debris
  energy then fails.**
- **Lump erasure is the only benchmark that tolerates debris.** It holds only while the thermalised debris energy
  stays below a lump threshold of about 20–175 E_X (for ε = 0.07–0.24).
  - Only one regular published release is nominally below it: the 8.03e12 GeV candidate, with a Y-only R of 47.
    Its σ = 0.001 ensemble average is 357, well above the threshold.
  - None of Codex's 31 robustness ensembles is below it.
  - The nearest are Codex's c_H = 10, H ≈ 4.35e11 GeV ensembles, at 115–124 E_X against 72–109 E_X, counting Y alone.
- **The c_S = 1 preparation fails the k_* vector piece.** That is the hypothetical preparation that rescues the
  homogeneous release. Under it, the release-induced longitudinal vector at the GMR peak scale k_* gives 13–24 times
  the allowance.
- **Our own earlier figure of about 1060 inflationary vector pairs per net charge is withdrawn.** It used the wrong
  (constant vacuum) mass; with the evolving mass the figure at 4.7e13 GeV is 1.4e-4.

## 1. Verdicts on Codex's handoff items

| # | Codex claim (373d37e) | Consolidated finding | Verdict |
|---|---|---|---|
| 1 | Point: ε −0.24066124, N_X −238.64868, E_X 741.64988, J_Y 13.26903, Y_X −65.9823 Y_Q | Same to all quoted digits with four codes: sp_point (Yoshida-8 and DP5, conformal), sp_yledger (cosmic-time DP5), sp_vector and sp_phaseprep [C] | CONFIRMED |
| 2 | δ + 1e-9 changes ε by 2.6e-7; rtol effect 6.7e-9 | 2.572e-7 (exact tangent 257.2 per rad); spread across integrator settings 1.9e-9 [C] | CONFIRMED |
| 3 | a = 2000–3000 mean shifts by 0.00046; "late fraction about 0.241" | 0.000460. The final value is ε = −0.240202, frozen to 2e-7 from a = 3000 to 3e4 with two integrators [C] | CONFIRMED (use −0.2402) |
| 4 | σ = 0.03, N = 128, seed 2026102600: −0.24909 ± 0.00333 | Identical draws give −0.249090 ± 0.003333 [C] | CONFIRMED |
| 5 | "A small-noise ensemble preserves it" | True only for σ ≲ 0.03–0.08 rad. σ = 0.1 gives −0.18 ± 0.01. The natural spread gives a net of −0.05 to −0.07, and full de Sitter draws give about 0 [C, E] | CORRECTED (scope) |
| 6 | Local map: structured, narrow (±2% in H, ±0.06 rad) | The productive region is a tongue about 1.5% in H by 0.2 rad (area 0.0034). In Codex's own σ = 0.03 metric only H/H* = 1.00 (within ±1%) and δ0 = π/2 (within ±π/12) pass [C] | CONFIRMED, narrower |
| 7 | The point is regular (passes the 1e-9 test) | Only locally. S = \|dε/dδ\| = 257 sits 12 e-folds below the trend ln S ≈ 2.1 (H/1e12 GeV) − 6.3, and S = 2e5–7e6 at ±5% in H. Of 161 published releases with \|ε\| ≥ 0.19, only 3 pass the 1e-9 test and only this one passes σ = 0.03 [C; Q data] | CORRECTED (island, not regime) |
| 8 | Conditional re-minimisation gives \|ε\| ≈ 0.0516 | 0.058 from the exact conditional state (Σ = h/√2); 0.048 with Codex's radii (Σ off by 3e-6) [C] | CORRECTED (0.05–0.06, not sharp; conclusion holds) |
| 9 | Late forced-solution subtraction: J_Y = 13.269 | Reproduced. The dynamic correction is 1e-5 in energy. The eigenmode-resolved adiabatic invariant is 13.6915, constant to 1e-7 [C, D] | CONFIRMED; count CORRECTED (+3.2%) |
| 10 | 0.0834 Y pairs per net charge = 16.9 times the allowance | Y alone: 0.0862 (17.5 times). Adding Σ-radial → AA (0.0411) and Y_H → Y_L A (0.0036–0.0123), with each A → X X̄, gives 0.131–0.140 (26.5–28.3 times) [C, E] | CORRECTED |
| 11 | Γ_Y = 377.95 GeV; a_dec = 1.72e5; E_Y ≈ 275 E_X | Γ is confirmed; the light eigenstate has 364.3 GeV. The energy deposited by the exponential decay is **349 E_X** (recomputed here: 348.99). H = Γ is low by 1.25; sp_point's t = 1/Γ (389) is not an average [D, C] | CORRECTED |
| 12 | Energy-equivalent fraction 0.00356 (diagnostic) | 0.0030 for Y alone; 0.0021–0.0022 with the vectors. Superseded by the computed relic (row 13) [C] | CORRECTED, superseded |
| 13 | "Do not apply 0.19 after the injection; relic not solved" | Benchmark relic computed: f = 15.2–15.9, about 1500 times, with case L equal to case H because no condensed phase forms. Thresholds are in §3 and §6c [C] | CONFIRMED caveat; result FAIL |
| 14 | Unprocessed Y source ≈ 17% of DM | 0.27–0.28 with the vectors. The frozen total with the cold X remainder is f ≈ 3.5 [C] | CORRECTED |
| 15 | Mean Y charge zero (assumed) | Y_L is real (CP-even), and there is no physical phase at dimension ≤ 4. 3N_Y/N_X = −8e-4 and falling [C, D] | CONFIRMED |
| 16 | 66-fold yield, not diluted | An inflaton matter era with T_R = 4.28e13 GeV cuts ε to 0.025. Over the 132 matter-era runs, none has \|ε\| ≥ 0.19 with Y ≥ Y_Q. A late visible-only entropy release needs a new field with hidden branching ≲ 1e-8 and leaves the leftovers unchanged [C, E] | NOT RESOLVED (every route checked fails or adds physics) |
| 17 | Counterexample to a universal X-only no-window claim, for a supplied mismatch | Holds as stated: a supplied, CP-odd, instantaneous, aligned-radii preparation [C] | CONFIRMED (conditional) |
| 18 | The baseline minimum is CP-even; no phase selector was tested | Inflation gives the required box with probability 0.03–0.1% per Hubble patch, and nil per observable universe [C, E] | CONFIRMED |
| 19 | Instant removal of the X Hubble mass (assumption) | The window needs removal within ≲ 0.003–0.01 H_i⁻¹, from single releases [C] | NOT CHECKED by Codex; new constraint |
| 20 | VECTOR_INITIALIZATION: F_g/H ≈ 34.83, m_A/H ≈ 3.483e-6, √P_α ≈ 0.00457 | 34.8322, 3.4832e-6, 0.004569 [C] | CONFIRMED |
| 21 | Do not substitute a constant-vacuum-mass relic | A constant mass overstates the result by (F_inf/F_vac)² = 4.2e5. The evolving-mass Bunch–Davies piece gives 1.6e-6–7.1e-6 pairs per net charge [C; re-derived here, §6b] | CONFIRMED |
| 22 | Masses at H_i (ENSEMBLE_SUMMARY) | Same to 5 digits [C] | CONFIRMED |

## 2. The lower-scale window

- **Reproduced.** Yes (row 1). The charge is frozen from a = 3000 onward, at ε = −0.240202 and N_X = −238.199.
- **Regular.** Only locally.
  - The exact tangent gives S = 257 at H*. The trend at that H is ln S ≈ 17.6, i.e. S ≈ 4e7.
  - There is no change of regime between 1.1e13 and 4.7e13 GeV. The charge forms when the Y mass crosses 3ω_X
    (the κ resonance, at a ≈ 37). Sensitivity keeps building until the Σ vacuum mass passes ω_X (a ≈ 88). Both
    epochs scale ∝ H/M_Y (sp_point [C, A]).
  - Scan-wide, from Codex's published columns [Q, C]:
    - of 2,203 grid releases, 161 have \|ε\| ≥ 0.19, and all of these have Y ≥ Y_Q;
    - none of the 33 coarse ones passes the 1e-9 test (Δε = 6e-3 to 0.4);
    - 3 refined ones pass, and they are Codex's three candidates;
    - only the leading point passes σ = 0.03.
- **Width.**
  - The tongue is about 1.5% in H by 0.2 rad, with thin streaks.
  - In Codex's own σ = 0.03 metric the full width at half maximum is 1.5–2% in H. Only H/H* = 1.00 clears 0.19:
    0.99 and 1.01 give −0.12.
  - The phase band is δ ≈ 1.46–1.64 rad.
- **Background requirements.**
  - An instant end of inflation, with the X Hubble mass gone within ≲ 0.003–0.01 H_i⁻¹ (sp_phaseprep).
  - No inflaton matter era (sp_yledger: ε → 0.025).
  - H_end within about 1%.

## 3. Combined leftover ledger at the point (per actual net charge, N_X = 238.20; f = 2.0285 × pairs)

| Source | Pairs per net charge | f | × 1% allowance | Basis |
|---|---:|---:|---:|---|
| Cold X remainder, frozen, no erasure | 1.58 | 3.21 | 321 | gp_leftovers number count [C] |
| Y_L → XXX / X̄X̄X̄ (N_Y = 13.69) | 0.0862 | 0.175 | 17.5 | [C] |
| Σ radial (N_σ = 4.89) → A A, each A → X X̄ | 0.0411 | 0.083 | 8.3 | [C, D] |
| Y_H → Y_L + A | 0.0036–0.0123 | 0.007–0.025 | 0.7–2.5 | bracketed [E] |
| Inflationary longitudinal vector (piece A, evolving mass) | 1.6e-6–7.1e-6 | ≤ 1.4e-5 | ≤ 0.0014 | [C], checked here |
| Release-induced Goldstone, Hubble-scale (B2) | 0.001–0.04 at σ = 0.03; 0.03–1.0 at natural noise | — | 0.2–8; 6–210 | [E], coefficients bracketed ×30 |
| Release-induced Goldstone at k_* (B1) | natural noise: unresolved; under the c_S = 1 preparation: 0.062–0.117 (±2 SE: 0.043–0.16) | — | 13–24 (9–32) | §6e [C inputs, E] |
| **Computed debris (rows 2–5)** | **0.131–0.140** | **0.27–0.28** | **26.5–28.3** | |
| **Frozen total (remainder + debris)** | | **3.47–3.49** | **≈ 350** | |
| **Maximal equilibrium, R = 510–542 E_X (case H = case L)** | | **15.2–15.9** | **≈ 1500** | gp_leftovers benchmark [C] |

**No double counting.**
- The pair rows are disjoint particle populations:
  - each Y quantum gives 1.5 pairs;
  - a heavy quantum is counted once in N_Y, and its vector separately;
  - each σ quantum gives two vectors;
  - each vector gives one pair.
- The cold remainder is the symmetric part of the X condensate.
- In equilibrium only energies enter, from disjoint pieces of the Hamiltonian (E_X, E_Y,exc, E_σ). The vector energy
  is taken out of the Y and σ budgets. The onset budget ρ_c + V0 remains the single budget for remainder, binding
  pairs and quench debris (gp_critic).

**Why no limit passes [C, A].**
1. *Frozen.* The cold remainder needs ε ≳ 0.99.
2. *Lumps that erase the remainder (case L).* This needs a condensed phase at the light-scale onset, which survives
   only if the thermalised debris energy is small. Our new threshold scan (§6c) puts the threshold at R_L ≈ 150–200 E_X for ε = 0.24
   (100–150 at 0.19, above 200 at 0.33). The point carries
   R = 510–542 E_X, so no condensate forms and f = 15.
3. *Homogeneous equilibrium (case H).* f = 0.01 is reached at R ≈ 0.69 E_X (ε = 0.24) and 0.014 E_X (ε = 0.19).
4. *The most favourable combination.* Suppose the remainder is erased by lumps while the hard debris never thermalises
   and is never captured. sp_yledger's rate estimates disfavour this. The debris pairs alone are still 27 times the
   allowance.

## 4. Phase-preparation verdict

- **The baseline action does not prepare the state** [C, E].
  - The light mode is the Σ phase (the μ source phase), with mass 0.61 H, f_eff = 0.612 H and a kick of 0.26 rad per
    e-fold.
  - The Fokker–Planck equilibrium has ⟨cos ψ₂⟩ = 0.85 and rms 0.58 rad, and relaxes in 10 e-folds.
  - The probability of the band and the |Σ| box is 0.03–0.1% per Hubble patch, and nil per observable universe (the
    memory factor is 1e-3 over 50–60 e-folds).
  - Natural ensembles give a net ε of 0 to −0.07 (sp_point, sp_yledger, sp_vector and sp_phaseprep agree within their
    coarse-grainings).
- **What a preparation must do** [C]:
  - give a phase spread ≲ 0.03–0.08 rad, which needs a phase mass ≥ 0.95–1.15 H (sp_phaseprep; this corrects
    sp_point's 1.5 H);
  - hold |Σ| within a few percent (+3% already fails) and |Y| within 0.2–0.5%, because |Σ| is itself a light mode
    (m² = 1.11 H²): its noise alone gives ε = −0.048 (sp_vector);
  - stay on until ≲ 0.02 e-fold before the end of inflation;
  - switch off within ≲ 0.03 H_i⁻¹, with the X Hubble mass gone within ≲ 0.003–0.01 H_i⁻¹.
- **Existence by hypothesis** [D, C].
  - A Dine–Randall–Thomas-type term, c_S H²|Σ|² + c_Y H²|Y|² + [μ_P Σ*Y² + c.c.], with c_S = 1, c_Y = 0.1035 and
    |μ_P| = 1.91e12 GeV at arg 1.835, has Codex's state as its exact minimum. It gives a phase mass of 1.153 H and
    ε_net = −0.242 ± 0.004.
  - Its costs:
    - a new CP phase;
    - three tuned coefficients;
    - a sudden end of inflation. Standard Hubble-tracking terms give |ε| ≤ 0.11.
  - It leaves the Y debris at 0.09–0.11 pairs per net charge (18–22 times).
  - Under it, B1 is not suppressed: 0.06–0.12 pairs per net charge, which is 13–24 times the allowance (§6e).
- **Verdict.** A preparation is possible only as a tuned hypothesis. **It is not the binding gate:** the leftover
  failure in §3 holds for any preparation.

## 5. Cross-track consistency and overreach

**Agreements.**
- The point is reproduced by four codes.
- The σ ensembles agree:
  - σ = 0.03: −0.2491, −0.2498 and −0.2462 (different seed);
  - σ = 0.1: −0.179 and −0.189;
  - σ = 0.26: −0.092 and −0.082.
- The natural-spread nets agree: −0.053 (0.52 rad), −0.058 (0.489 rad), −0.072 (0.29 rad, all modes) and −0.048
  (radial noise only).
- The Y pairs agree: 0.0861 and 0.0862.
- The frozen-angle factor agrees: sp_vector 2.357e-6 from the mode equations, sp_yledger about 2.4e-6, and this track
  2.356e-6 (§6b).

**Disputes resolved.**
1. *Y energy at decay* (275 / 389 / 349). The deposit of an exponential decay into a radiation-like sector is
   N m ⟨a⟩ with ⟨a⟩ = Γ(3/2)√(2H_i/Γ). Recomputed: 348.99 E_X. Use 349.
2. *Required phase mass* (1.5 H against 0.95–1.15 H). The Starobinsky–Yokoyama Gaussian overestimates the spread at
   m ~ H. Use sp_phaseprep's exact Bunch–Davies value at x_c = 1/40. B1 adds its own requirement (§6e).
3. *Freshness* (0.05 against 0.02 e-fold). sp_point's own runs give −0.202 at 0.02 e-fold and −0.167 at 0.05, so the
   requirement is ≲ 0.02–0.03 e-fold.
4. *Radial tolerance.* sp_point's "radii within about 5%" means 5% of the path to the conditional minimum, which is
   4.9% in |Σ| and 0.01% in |Y|. The individual tolerances are |Σ| about −3% to −5% and below +3%, |Y| 0.2–0.5%, and
   |X| 0.1%. sp_yledger's 2e-7 radial sensitivity measures how reproducible one trajectory is; it is not a physical
   tolerance.
5. *Natural spread* (0.52 / 0.489 / 0.29 rad). These are different coarse-grainings, and every one gives a net
   ε ≤ 0.07.

**Overreach flags.**
1. *sp_yledger, "≥ 350 times in every limit".* The 350 assumes no lump erasure of the cold remainder. The floor across
   limits is about 27 times (§3, point 4). The verdict, FAIL, is unchanged.
2. *sp_yledger, maximal equilibrium as "the relevant regime".* It rests on order-of-magnitude rates [E]. Read 1500
   times as conditional on the benchmark, and 27–350 times as the floor.
3. *sp_vector, "the vector sector fails the point".* B2 rests on a white-noise model whose coefficients are bracketed
   (×30). B1 is unresolved at natural noise, and at natural noise ε has already failed. Under Codex's σ = 0.03, B2 is
   0.2–8 times, which is inconclusive. Read it as a probable additional failure, not a certified one. Under the
   c_S = 1 preparation see §6e.
4. *sp_point.* Its "≳ 1.5 H", "about 0.05 e-fold" and "radii within 5%" are corrected as above.
5. *sp_phaseprep.* Its switch-off thresholds come from single chaotic releases, and its "nil" for the observable
   universe is an estimate within stochastic equilibrium [E].
6. *Our own earlier review* (REVIEW_GAUGEPOINT, gp_critic, gp_leftovers).
   - "Longitudinal vector ≈ 1060 pairs per net charge, sufficient on its own" is **withdrawn** for the inflationary
     Bunch–Davies source. At 4.7e13 GeV that source gives 1.4e-4 (sp_vector, confirmed in §6b).
   - The 4.7e13 GeV leftover FAIL still stands, on net ε ≤ 0.03 alone.
   - gp_critic's suggested remedy e ≳ 0.029 triggers a Floquet resonance off the X orbit (μ ≈ 0.15 per unit a) and is
     not a free fix.

## 6. New consolidation results (this track)

**a. Y-decay energy [D, C]** (`critic_checks.py` A).
- Γ = 377.95 GeV and Γ_L = 364.34 GeV.
- The decay epoch is a = 1.722e5 at H = Γ and 2.436e5 at t = 1/Γ. The decay average is ⟨a⟩_L = 2.199e5.
- R_Y = 274.85 (H = Γ), 388.7 (t = 1/Γ), 344.5 (decay average with J_Y and M_Y), and **348.99** (decay average with
  N_Y,eig and m_L).

**b. Frozen superhorizon gauge angle [D, C]** (`critic_checks.py` C).
- Derivation. In the relativistic superhorizon regime, (z²B′)′ + a²m²B = 0 reduces to (z²B′)′ ≈ 0. Its non-decaying
  solution is B = const for any history m(τ), so the late energy scales as (m_vac/m_inf)².
- Independent ODE: Bunch–Davies start at k/aH = 1000, de Sitter followed by radiation, out to 30 a_*. Ratio of P1 to
  P0:

  | k [H_i] | P1/P0 | instant step / P0 |
  |---|---:|---:|
  | 3e-6 | 2.3565e-6 | 2.3607e-6 |
  | 1e-5 | 2.3561e-6 | 2.3605e-6 |
  | 3e-5 | 2.3530e-6 | 2.3589e-6 |

  The expected value is (m_vac/m_inf)² = 2.3607e-6. The late number changes by ≤ 0.14% between 20 a_* and 30 a_*.
  **sp_vector's piece A is confirmed, and our earlier 1060 is withdrawn.**

**c. Debris-energy thresholds in the gp_leftovers benchmark [C]** (`rmax_relic.py`, reusing sp_yledger/gp_leftovers).

| ε | Case H (homogeneous): R at f = 0.01 | Case L (lumps + conversion): f = 0.01 crossed between |
|---|---:|---:|
| 0.0677 | — | R = 10 and 30 (the condensate is gone at 30) |
| 0.10 | — | 30 and 50 |
| 0.15 | — | 75 and 100 |
| 0.19 | 0.014 | 100 and 150 |
| 0.2402 | 0.69 | 150 and 200 (sp_yledger: f = 9.3 at 275) |
| 0.3267 | 2.06 | above 200 |
| 0.598 | 7.3 | — |

- **Case H tolerates essentially no debris.**
- **Case L tolerates R_L ≈ 20–175 E_X for ε = 0.07–0.24.** Below R_L the condensed phase survives the hot hidden
  plasma (T_h/m ≈ 8–12), and f is negligible. Just above it, f rises steeply, to 0.05–10.
- **The point fails both.** It carries R = 349 E_X (Y alone) to 510–542 E_X (Y plus vectors). Codex's own σ = 0.001
  and σ = 0.03 ensembles there give a Y-only R of 318 and 421 E_X.

**d. Scan-wide Y debris [Q inputs, C]** (`critic_checks.py` D; Codex J_Y and E_X, decay-averaged).
- Over the 161 releases with |ε| ≥ 0.19, R_Y = 16–2506 E_X (median 454), and the Y pairs per net charge are
  0.0041–2.4 (median 0.15). This is Y only, so it is a lower bound.
- The lowest coarse values (R_Y = 16–50, at c_H = 10, where E_X is about 40 times larger) all fail the 1e-9 test.
- The three regular candidates have R_Y = 47 (8.03e12 GeV, ε = 0.190; fails σ = 0.03), 161 (1.58e12 GeV, c_H = 10;
  fails σ = 0.03) and 344 (the leading point).
- **Against the thresholds [C].**
  - All 161 releases fail case H; the smallest ratio R_Y/R_max,H is 7.8.
  - 16 have a Y-only R_Y below the case-L lower threshold:
    - one passes Codex's 1e-9 test. It is the 8.03e12 GeV candidate, whose σ = 0.001 ensemble R is 357, not 47, so the
      nominal value was a fluctuation;
    - 10 fail the 1e-9 test;
    - 5 are local-map points beside H*. Their σ = 0.03 neighbours give |ε_net| ≈ 0.12 (sp_point).
  - Among the 159 regular, productive releases (\|ε\| ≥ 0.01, Y ≥ Y_Q), the smallest R_Y is 40.7, at ε = 0.068 where
    R_L ≈ 10–30.
- **Codex's robustness ensembles** (31 ensembles; `critic_checks.py` E): **none** has a Y-only R_ens below the case-L
  lower threshold. The closest are candidates 10 and 17 (c_H = 10, H = 4.35e11 GeV, smooth map, ε_net = 0.157 and
  0.147, Y = 5.7 and 5.3 Y_Q): R_ens = 115 and 124, against R_L ≈ 80–109 and 72–97. This is before the Σ and vector
  debris are added; at the leading point they raised R by 46%.

**e. B1 under the best hypothetical preparation [C inputs, E conversion]** (`b1_ppi_slope.py`).

- **Response.** Phase noise of σ_S = 0.026 rad is the residual spread of the c_S = 1 term. Under it, the
  noise-averaged response over ψ0 = π/2 ± 0.04 (5 offsets × 256 common draws) is resolved:
  - d⟨θ_Σ(∞)⟩/dψ0 = 3624 ± 310 rad/rad (paired estimate 4030 ± 359);
  - d⟨ε⟩/dψ0 = 0.73 ± 0.04 per rad, with ⟨ε⟩ = −0.250 ± 0.002 at π/2.

  At natural noise, by contrast, sp_vector found 65 ± 132, consistent with zero.
- **Size at k_*.** The preparation's light phase (1.153 H, ν = 0.960, f_eff = 0.635 H) leaves
  δψ(k_*) = 1.16e-3 rad per ln k at the end of inflation. That gives δθ_G(k_*) = 0.70 rad per ln k, against
  sp_vector's criterion of 0.226.
- **Vectors.** Through sp_vector's P0 transfer, with the preparation's blue spectrum (P_ψ ∝ k^1.08), B1 gives
  0.062 (k ≤ 1/3000) to 0.117 (k ≤ 1/1000) pairs per net charge. That is **13–24 times the allowance**, or 8.7–32 times
  across ±2 SE.
- **Reading.** The c_S = 1 term fixes the homogeneous release but not the long-wavelength Goldstone source.
  - A heavier term (c_S = 2, phase mass 1.41 H, ε_net = −0.239 in sp_phaseprep) would lower the k_* power by about
    3e-4, if the smoothed slope stayed the same [E].
  - Its smaller residual spread (0.002 rad) smooths less, so the slope may rise toward the local value of 1.3e6
    rad/rad. That case is not computed.
  - The preparation's radial modes (1.40 H and 1.45 H) are not included. Neither is the |Σ| response (−3.65e6 per
    unit ln|Σ| in the linear regime).
  - **Scope:** the conversion inherits sp_vector's separate-universe premise and frozen-angle transfer [E].

## 7. Corrections for Codex

1. **Window.** Report the final fraction −0.2402. Describe the window as an island: S = 257 against a trend of about
   4e7, a tongue about 1.5% × 0.2 rad, and only H/H* = 1.00 ± 1% with δ0 = π/2 passing σ = 0.03. Of the 161 releases
   with |ε| ≥ 0.19, 3 are regular and 1 is robust.
2. **Noise.** "A small-noise ensemble preserves it" holds only for σ ≲ 0.03–0.08 rad. The natural inflationary state
   gives ε_net ≈ 0–0.07.
3. **Conditional re-minimisation.** |ε| ≈ 0.05–0.06, which is preparation-sensitive at 1e-6.
4. **Y ledger.**
   - N_Y = 13.69 (Majorana eigenmodes).
   - Energy deposited at decay: 349 E_X, from the decay average.
   - Y pairs: 0.0862 per net charge.
   - Add the Σ radial mode (N_σ = 4.89, σ → AA) and Y_H → Y_L A. The debris then totals 0.131–0.140 pairs per net
     charge (26.5–28.3 times the allowance), 0.27–0.28 of DM, and an injection of 510–542 E_X.
5. **Relic.** Replace the energy-equivalent diagnostic with the computed benchmark. The relic is f = 15.2–15.9 (case L
   equals case H), and the frozen total is f ≈ 3.5. State that the X-only 0.19 pass holds only in a benchmark that the
   release's own debris then fails.
6. **Release history.** State the instant-release requirements: the X Hubble mass gone within ≲ 0.003–0.01 H_i⁻¹, and
   no inflaton matter era (the dilution that removes the factor 66 drives ε to 0.025).
7. **Preparation.** Any preparation must fix |Σ| (a light mode, m² = 1.11 H²) and |Y|, not just the phase. The phase
   mass must be ≥ 0.95–1.15 H for the homogeneous release. At c_S = 1 (1.153 H) the k_* vector piece B1 is still 13–24
   times the allowance (§6e).
8. **Vectors.**
   - The evolving-mass Bunch–Davies piece is negligible (≤ 7e-6 pairs per net charge).
   - List the release-induced Goldstone (patch-to-patch Σ winding, about 2300 rad) as an open source.
   - e ≳ 0.029 is not a free remedy.
   - CORRECTIONS item 8 should cite our 1060 figure as withdrawn: that was our error.
9. **Thresholds.** Quote the debris-inclusive thresholds:
   - case H: R ≤ 0.01–7 E_X;
   - case L: R ≤ 20–175 E_X.

   Report the ensemble-averaged J_Y, N_σ and vector energies for each candidate. The nominal J_Y can be
   unrepresentative: for the 8.03e12 GeV candidate the ensemble value is 7.6 times larger.

## 8. The single most decisive next calculation

**The debris-inclusive case-L ledger on the regular low-scale branch (c_H = 10, H ≈ 3–6e11 GeV).**

*Why this calculation.*
- Every published release fails the frozen and homogeneous-equilibrium limits by large factors (§3, §6c–d).
- Lump erasure (case L) is the only benchmark that tolerates debris.
- The only published ensembles near its threshold are Codex's candidates 10 and 17 (ε_net ≈ 0.15, Y ≈ 5.5 Y_Q). There
  the release map is smooth (local ε std 3e-5), so the chaotic-preparation problem of the 1.12e13 GeV island does not
  arise. A supplied CP phase is still needed.

*What to compute there.*
1. The full decay-averaged injection R_tot: the Y eigenmodes, the Σ radial quanta (N_σ) → AA, and Y_H → Y_L A, using
   sp_yledger's ledger code.
2. The case-L threshold R_L at the actual ε_net and undiluted yield (ymult ≈ 5.5).
3. ε_net under the computed inflationary spread at that H, and under a stated preparation.
4. Then:
   - **If R_tot < R_L,** the leftover gate stays open in the case-L benchmark. The decisive physics then becomes the
     case-L premise: lump formation plus the bb211d1 conversion in a T_h/m ≈ 10 plasma, which needs a
     lattice/Boltzmann check. The preparation and the k_* vector piece follow.
   - **If R_tot > R_L,** no release in Codex's published scan passes any leftover benchmark. The one-mediator
     "creates and protects" link then fails at these couplings for every supplied preparation, before fragmentation,
     strings, walls or the vector piece B are considered.

*Cost.* The existing release codes (Codex's, sp_point, sp_yledger) and rmax_relic.py can do this in a few hours.

## 9. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/sp_critic/code
P=<scratch>/pkt/research/reality_framework/grut_inflation_scan_20261025   # git archive 373d37e extraction
python3 -B critic_checks.py $P ../.. ../results/critic_checks.json ../results/rmax_relic.json ../results/rmax_relic_L.json \
    ../results/rmax_relic_Llow.json ../results/rmax_relic_Hhigh.json ../../sp_yledger/results/thermal_relic.json  # ~1 min
python3 -B rmax_relic.py ../results/rmax_relic.json                                     # ~25 min, 4 cores
python3 -B rmax_relic.py ../results/rmax_relic_L.json 0.240201575085664,0.19,0.3267 30,50,75,100,150,200 L
python3 -B rmax_relic.py ../results/rmax_relic_Llow.json 0.0677,0.1,0.15 3,10,30,50,75,100 L
python3 -B rmax_relic.py ../results/rmax_relic_Hhigh.json 0.3267,0.598 0,1,3,10,30 H
python3 -B b1_ppi_slope.py ../../sp_vector ../results/b1_ppi_slope.json                 # ~10 min, 4 cores
```
rmax_relic.py imports sp_yledger/thermal_relic.py and gp_leftovers (writes no bytecode there). b1_ppi_slope.py
compiles sp_vector/code/vrel.c into the scratchpad.

Credit: Affleck–Dine; Dine–Randall–Thomas (Hubble-induced terms); Graham–Mardon–Rajendran (longitudinal vector from
inflation); Starobinsky–Yokoyama; Bunch–Davies; Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki (Q-ball
formation and fragmentation); Greene–Kofman–Linde–Starobinsky (λφ⁴ resonance); Carlson–Machacek–Hall (cannibal
freeze-out); Goldstone equivalence (Cornwall–Levin–Tiktopoulos; Lee–Quigg–Thacker) [R]. Mixed or decaying dark
matter structure constraints are not used.
