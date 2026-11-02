# Skeptic of track fx_formation (η and the object distribution under M1)

1 November 2026 · Claude, issue #13 repair pass (BRIEF.md of this directory). Adversarial check of `NOTES.md` and
`results/*.json` of this track. Same-session AI work on one supplied model. It is not outside replication, not a
validation and not an exclusion. Labels: [D] derived, [C] computed with our own code, [E] estimate, [A] argued,
[R] recalled. **A repair is a model choice, not a consequence of closure.** Git was read-only and nothing was
committed. No 3D runs were made. Code: `skeptic/code/`, results: `skeptic/results/` (140 kB in total).

## Verdict in one paragraph

The central result survives an independent recomputation. I computed the vapour thermodynamics with a Bessel series
instead of the track's momentum quadrature. At the same μ the two agree to 1.5×10⁻¹⁵. With my own coexistence
condition I get Y_crit = 0.288054 at T_h = 0.3360 m_X (S/Q < 3.4716) [C]. A homogeneous linearisation that I derived
myself gives the spinodal growth law Γ = 0.775 δ, k = 0.958 √δ [C]. My own freeze-out ODE gives 1 − η = 0.026–0.065
(R = 0) and 0.044–0.112 (R = 0.91) [E]. So the physics holds: the condensed phase evaporates when S/Q > 3.47, the
heavy e = 10⁻⁶ release has S/Q = 4.95, and the low-heat branch gives macroscopic Q-matter at E/Q ≈ 0.776.

The track has three concrete defects:
1. The proposed decisive 3D run cannot be done. A 24/m_X box holds 1.8×10³ units of vapour entropy, but S/Q = 5 for
   one Q = 2×10⁵ object needs 10⁶.
2. One of its two Y_h values is an artefact (S/Q = 5.4).
3. The pair leftover p is taken from the Codex reproduction case, not from fx_chain's M1 relaxed history.

It also misreads its own verdict. The e = 10⁻⁶ band fails only the *frozen test*, which was built for the supplied
Q = 2×10⁵ object. Under the same recalibration that the track offers for the low-heat branch, the e = 10⁻⁶ band
becomes a clean free-X cold relic: M = 1.0018 m_X per charge, H_i × 0.957, λ_fs ≈ 0.3–0.5 pc. Repairs (ii) and (iii)
are not alternatives to recalibration, because neither restores 0.9374.

## Checks, result by result

| Track result | Recomputation (own code) | Status |
|---|---|---|
| family: fit, σ, first law | Fit residual 3.24×10⁻⁵ and σ = 98.78 m_X³ reproduce. The first-law match "≤ 2×10⁻⁵" holds only for Q ≥ 5×10⁵; at the fit's edge Q = 3×10⁵ the derivative is off by 2.6×10⁻⁴ [C] | confirmed; digits corrected |
| Y_crit = 0.288 at 0.336 | Bessel series against quadrature: identical to 1.5×10⁻¹⁵ at the track's μ_c. Own μ_c = Ω_b + p_vap/n_B gives **0.288054 at 0.3360**. Finite objects: 0.3415 (2×10⁵), 0.3143 (10⁶), 0.2930 (10⁸), 0.2885 (10¹¹). Thermal masses δm²/T² = 10⁻³ / 10⁻² / 10⁻¹ shift it by −0.006% / −0.06% / −0.6%. dY_crit/dμ_c = +0.14% per 10⁻³ [C] | **confirmed** |
| Y_h and S/Q from fx_chain | **Y_h(spinodal) = 0.1854 is an artefact.** It divides the *condensate* density (233.2) by s. Since T_h ∝ 1/a and n_Q ∝ a⁻³ along fx_chain's history, the conserved value is Y_h = n_Q/s = 0.2021 = 0.1854 / 0.9174 (the condensate fraction). So S/Q = 4.95, not 4.95–5.4. The massless 3-dof entropy is accurate: with masses fitted to the Hartree populations (m_φ = 1.73, m_X = 0.57 m_X) s/s₀ = 0.994–0.999, and the pilot's S/Q is 3.604 [C/E] | corrected (S/Q = 4.95) |
| T_ev / T_s | Root-found (not grid): heavy e1e6 **0.802 / 0.0962**. Pilot **0.452 / 0.238** (track 0.47 / 0.24). e1e7 4.04 / 0.0116. The ends 0.94 and 0.078 come only from the artefact [C] | corrected digits |
| evaporation is complete | The track uses the gross exchange rate ("Q < 2×10³⁷ evaporate within 1/H"). Net Hertz–Knudsen flux with the computed fugacity deficit, integrated over the actual window: the deficit peaks at 8.4% (pilot) and 65% (heavy), over Δln a = 0.58 / 2.08. That evaporates objects up to **Q ≈ 10³⁰ (pilot, s = 1), 10²⁴ (s = 0.01)** and 7×10³⁴ / 7×10²⁸ (heavy) [E]. Still far above ⟨Q⟩ ≤ 7×10¹⁹ | conclusion confirmed; number overstated by 10³–10¹³ |
| spinodal Γ, k, Q_frag | Own linearisation (ψ = 10F, φ = 10S, V = 100U) gives the same cubic. **c₁ = 0.775, c₂ = 0.958.** Q_frag = 3.8×10¹¹–1.15×10¹³, H/m_X = 1.31×10⁻¹³ (1.39×10⁻¹³ if n_Q = 254 rather than n_c is used to set the clock) [C] | confirmed |
| ⟨Q⟩ after ripening | Wagner spot check at T_h = 3: R² = (64/81) k_W t gives Q ≈ 2×10¹⁸, consistent. ⟨E/Q⟩ is at the bulk for any Q ≳ 10¹² (0.7767), so the verdict does not depend on it [E] | confirmed (order of magnitude) |
| η (low heat) | Own two-phase path and relaxation ODE, Q = 8.5×10¹²–10²⁰, s = 1 / 0.01. **R = 0: 1 − η = 0.026–0.065. R = 0.91: 0.044–0.112.** The entropy-floor formula matches to 3–5%. T_f = 0.005–0.015 m_X and T_vis(T_f) = **2 keV–20 MeV** (the track gives 0.005–0.008 and 3 keV–0.3 MeV, which leaves out its own s = 0.01, Q = 10²⁰ corner) [E] | confirmed; range widened |
| warm tail | λ_fs (comoving, v ∝ 1/a after T_f, ΛCDM background) = **0.015–62 kpc**. Fraction and λ_fs are anti-correlated: 2.6–4.7% at 29–62 kpc (late freeze-out, small Q, s = 1), and 6.5–11% at 15 pc (s = 0.01, Q = 10²⁰). "3–11% at 0.5–40 kpc" pairs the extremes wrongly, and 62 kpc exceeds the quoted upper end [E] | corrected |
| verdict, e = 10⁻⁶ band | The track uses p = 0.0045855/2.0285 = 2.26×10⁻³, the Codex-reproduction (old-placement) case. fx_chain's M1 relaxed gas-only value is **9.10×10⁻⁴**, so **M = 1.0018 m_X (+6.87%) and the test reads 0.0687 (×6.9)**. For e = 10⁻⁷ the corresponding numbers are 0.723 and M = 2.446. Neither network was run on the evaporated history, so p remains [E] there | digits corrected; conclusion holds |
| threshold R_crit, e_crit | R_crit = 2.757 / 4.668 reproduce. With the cascade vectors scaled ∝ 1/e: **e_crit = 2.13×10⁻⁶** heavy (track 2.20×10⁻⁶) and 1.077×10⁻⁶ pilot | digits corrected |
| repair options | (iii) gives η → 1 at E/Q → 0.776, i.e. −17.2%, so H_i × 1.134 is needed anyway. (ii) gives −14.5…−16.6%. **Recalibration is required in every branch**, and it applies equally to the e = 10⁻⁶ band: M = 1.0018, H_i × 0.957 = 1.285×10¹¹ GeV [C arithmetic on E inputs] | **refuted as alternatives** |
| decisive 3D test | At T_h = 0.5, μ = 0.8817: s_vap = 0.158 m_X³ and n_net = 0.050 m_X³. The 24³ box (dx = 0.5) minus the object holds S_vap = 1.79×10³, so **S/Q in the box = 0.009, not 5 or 3**. S/Q = 5 (3) needs L ≈ 185 (156)/m_X, i.e. 370³ at dx = 0.5, about 8× the cap in each direction. Even at L = 48, dx = 1, S/Q ≥ 3 allows only Q ≤ 5.7×10³ < Q_min. In a small box at fixed energy the object cannot evaporate into a 1.4×10⁴ m_X⁻³ volume either (n would be 14 m_X³, deep in BEC) [C] | **refuted (infeasible as designed)** |

## Corrected physical picture

Each item below is conditional on M1, the supplied inputs and fx_chain's heat budget.

1. **e = 10⁻⁶ band (S/Q = 4.95 heavy, 3.60 pilot).**
   - The condensed phase evaporates between T_ev ≈ 0.80 (0.45) and T_s ≈ 0.096 (0.24) m_X and never re-forms.
   - The dark matter is a **free, cold, asymmetric X gas**: fugacity 0.10–0.39 at T_s, T_vis(T_s) ≈ 1–3 GeV,
     v ∝ 1/a afterwards, λ_fs ≈ 0.3–0.5 pc (100× colder than Codex's gas).
   - The relic has M = 1.0018 m_X per charge, plus pairs p ≈ 9×10⁻⁴ (borrowed from a network that did not include the
     evaporation history [E]).
   - It fails the frozen test (0.069), but that test assumes Q = 2×10⁵ objects. Recalibrating on M gives
     H_i × 0.957 [C].
2. **Low heat (S/Q < 3.47; e ≳ 2.13×10⁻⁶ heavy, 1.08×10⁻⁶ pilot).**
   - The relic is macroscopic Q-matter, ⟨Q⟩ ~ 10¹²–10²⁰ and E/Q ≈ 0.776, plus a free-X tail of 2.6–11%.
   - M = 0.776 + 0.224(1 − η) = 0.782–0.801, which needs H_i × 1.11–1.13.
   - That spread is ±1.2% from 1 − η alone, wider than the program's 1% tolerance. So the recalibration cannot be
     quoted to 1% until sticking and ⟨Q⟩ are computed [E].
3. **Neither branch contains Q = 2×10⁵ objects.** The supplied object and the frozen-test form are retired in both.
   Calibration must be redone on the computed relic in both, with the D-to-H_i feedback iterated (R(e), Y_h and the
   clocks were all computed at H_i = 1.343×10¹¹ GeV) [A].

## New issues

1. **Y_h(spinodal) artefact.** The track divides n_condensate by s. Remove it: S/Q = 4.95, T_ev = 0.802 and T_s =
   0.096 for the heavy e1e6 band.
2. **Wrong p in the evaporated branch.** p is the Codex-reproduction case (2.26×10⁻³), not fx_chain's M1 relaxed
   9.10×10⁻⁴. Neither was computed on the evaporated history (all charge in gas from T_ev, latent-heat cooling).
3. **The decisive 3D run is infeasible.** Replace it with:
   - a thermostatted, grand-canonical measurement of the saturated n_X(T, Ω(Q)) and the sticking probability at
     T_h = 0.5;
   - or a quantum-kinetic vapour coupled to the object.

   S/Q is a property of the cosmological vapour and enters only through the ideal-gas thermodynamics checked here.
   A classical lattice at T = 0.5 m_X with dx = 0.5 is in the Rayleigh–Jeans regime: the vapour density in such a box
   is not the Bose value.
4. **The e = 10⁻⁶ "failure" is a frozen-test failure, not a model failure.** After recalibration that branch is the
   simpler relic, with no warm tail, no nuggets and negligible pairs.
5. **Repairs (ii) and (iii) do not restore the calibration.** (iii) also adds a light state with its own ΔN_eff and
   coupling ledger. It would open XX̄ → light-pair annihilation, which lowers p, and a light scalar would add a long-range
   force between X and nuggets [A].
6. **The pilot is marginal.** S/Q = 3.604 against 3.472 (3.8%). R_crit = 4.67 is within 6% of R = 4.97, and the band
   versus broad vector source differ by 15% in R. A H_i refit or a better heat budget could move the pilot across the
   threshold. The heavy release is robust (×1.43).
7. **Warm tail.** Fraction and λ_fs are anti-correlated (2.6–4.7% at 29–62 kpc). Mixed warm-plus-cold bounds at the
   few-percent level are usually weak [R], but the 62 kpc corner should be put through a transfer calculation.
8. **Minor.** The first law is off by 2.6×10⁻⁴ at Q = 3×10⁵. e_crit is 2.13×10⁻⁶. The H at the spinodal is set with
   n_c instead of n_Q (a 6% change in H, 2% in Q_frag).

## Reproduction (OMP_NUM_THREADS=1; SP scratch = $SP/wt2/fx_formation_skeptic)
```
T=research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_formation
python3 -I $T/skeptic/code/bessel_vapour.py $T/skeptic/results/bessel_vapour.json $T/results/eq_history.json   # ~2 min
SKEP_SCR=$SP/wt2/fx_formation_skeptic python3 -I $T/skeptic/code/skeptic_checks.py $T/skeptic/results/skeptic_checks.json  # ~6 min (restart file for spinodal/window)
python3 -I $T/skeptic/code/onset_entropy.py $T/skeptic/results/onset_entropy.json
python3 -I $T/skeptic/code/evaporated_fs.py $T/skeptic/results/evaporated_relic_fs.json
```
`skeptic/results/skeptic_summary.json` collects the key numbers, and `evaporated_relic_fs.json` holds the free-X
relic's free streaming. The track's files were read only. The g*(T) table is approximate [R].
