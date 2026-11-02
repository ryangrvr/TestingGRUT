# SKEPTIC of track fx_observables (key fx_observables_skeptic)

1 November 2026 · Claude, issue #13 (BRIEF.md of grut_working_theory_v2_20261101) · same-session adversarial AI review
of one supplied model. It is not outside replication, a validation, or an exclusion of the real world. Git was
read-only, nothing was committed, and no Codex code was run. All codes are new, under `skeptic/code/`. The one track
module that is reused (`code/vapour_transfer.py`) is imported read-only for its background and history functions; its
loop and marching are copied before they are changed. Labels: [D] derived, [C] computed with our code, [E] estimate,
[A] argued, [R] recalled.

## Verdict in one paragraph

Every particle-physics number the track gives reproduces with independent methods:
- the Q = 2×10⁵ object, by collocation instead of Numerov plus Newton;
- the form-factor zeros and σ_n0 per coupling²;
- the LZ calibration and coupling limits, from a Monte-Carlo velocity integral;
- the light-WIMP mimic, from a Poisson Asimov test;
- the h → φφ freeze-in yield;
- the trilinear freeze-in, from my own 2→2 collision estimate.

Two things fail.
1. **The physical picture in §2B is wrong where it matters most.** "While coupled the vapour makes the whole DM a single
   fluid with c_eff ≈ 0.01c" needs the vapour mean free path to be short against the wavelength. But the vapour moves at
   v̄ ≈ 0.07–0.12c, about 10× c_eff. For every ζ ≳ 5×10⁻⁴ case (Q ≤ 10¹⁰ at s = 1), the modes at k = 10–50/Mpc are not
   tightly coupled when their acoustic phase starts: ε_TC = k v̄/(aΓ_m) = 1.7–9.
   - A two-fluid re-integration changes P(k) at k ≳ 1/r_s by large factors. At Q = 2×10⁵, s = 1, P(10) is 0.42 (0.20),
     not 0.05 (0.41). At ζ = 2×10⁻³, P(50) is 0.52 (0.31), not 0.06 (0.02).
   - The regime map survives. ζ ≲ 2×10⁻⁵ and ζ ≈ 10⁻⁴ are unchanged, because those modes enter after decoupling. The
     reference Q = 2×10⁵ with s = 1 stays strongly suppressed. The ζ ≈ 2×10⁻³ exclusion becomes marginal: P(10) = 0.84–0.90
     against 0.94.
2. **The merger mass-conversion figure is low by ×2.** It should be 3.2×10⁻³, not 1.6×10⁻³, in the inner kpc. The
   conclusion (negligible) is unchanged.

New issues:
- X̄ burn-out is not free: it deposits ΔS/Q ≈ 0.4–1.5 into the hidden sector for p = 0.002–0.005.
- The "2–6% vapour" holds only for s ≈ 1.
- The 1%-CDM threshold is ζ ≲ 1×10⁻⁵, not 2×10⁻⁵, in the medium-heat branch.
- The trilinear's survival hinges on the Higgs–φ mixing through the electroweak crossover.

## 1. Checks of the track's results

| id | claim | my recomputation | verdict |
|---|---|---|---|
| O1 | Q = 2×10⁵: M = 17 324.83 GeV, r₅₀/r₉₀ = 3.79/5.94, σ/M = 2.7–6.6×10⁻⁴ cm²/g | `indep_object_dd.py`: solve_bvp collocation on (F, F′, S, S′) at Ω = 0.881720 with Robin tails, from a tanh guess. **Q = 2.00002×10⁵, E/Qm = 0.937399, M = 17 325.0 GeV**, S(0) = 0.1258, r₅₀/r₉₀/r₉₉ = 3.81/5.96/8.44, σ/M = 2.69/6.58×10⁻⁴. The "φ-half radius 4.25" is the midpoint between S(0) and 1 (S = 0.563); S = 0.5 sits at 4.02 | holds [C] |
| O2 | merger releases 3.7% (1.29 TeV); ≤ 1.6×10⁻³ of local DM mass converted | The release reproduces from the E(Q) fit (E/Q(4×10⁵) = 0.90244, giving 3.73%). **The conversion is wrong by ×2.** Per volume, ½n²σv collisions each convert f_rel·2M, so the converted fraction per t_H is nσv·t_H·f_rel. `late_and_portal.py` divides by a further 2. Corrected: 2.4×10⁻⁴ (Sun), **3.2×10⁻³ (inner kpc)**, 2.7×10⁻³ (cluster core), 1.5×10⁻⁴ (dwarf) | sub-number refuted, conclusion holds [C] |
| O3 | Codex's elastic gas is CDM (decoupling 2.64 keV) | Γ/H ∝ T² with the g* steps gives T_dec ≈ 2.3–2.6 keV from 1.2×10⁶ at 2 MeV. Tight coupling is trivially satisfied (c_eff ~ 10⁻⁴c, r_s = 3.7 pc) | holds [D] |
| O4 | equilibrium vapour; single fluid with c_eff ≈ 0.01c until a_m; ζ table | see §2 | regime map holds; single-fluid picture and every number at k ≳ 1/r_s refuted for Q ≤ 10¹⁰, s = 1 [C/E] |
| O5 | X̄ burned out, Γ/H = 1.2×10⁵ at 2 MeV | n_core(2 MeV) = 1.21×10⁻¹⁸ m_X³ × 1.95/m_X² over H = 1.79×10⁻²⁴ GeV gives 1.22×10⁵ | holds; hidden cost in §3 [D] |
| O6 | ΔN_eff = 2×10⁻⁹–2×10⁻⁸ | DM/ν₁ = 5.0×10⁻⁵ at 0.1 MeV × vapour kinetic fraction 3.3×10⁻⁴ gives 1.7×10⁻⁸ | holds [D] |
| O7 | zeros: φ² 27.56/48.80/90.27, trilinear 34.11/60.38/111.7, λ_HX 41.68/73.79/136.5 keV | identical to 4 digits from my own profile. Isotope note: Ge is ⁷⁴Ge (natural Ge gives 49.72/61.52/75.18). Xe is A = 131 (natural Xe gives 27.50/34.03/41.59) | holds [C] |
| O8 | σ_n0 = 1.40×10⁻²⁶ b² cm²; b < 6.5×10⁻¹⁰; XLZD 10 events at 2.2×10⁻¹⁰ | Own profile integrals: ∫(S²−1) = −581.72, ∫(S−1) = −380.20, ∫F² = 2268.3, giving σ_n0/b² = 1.398×10⁻²⁶, σ/A² = 6.99×10⁻²⁷, σ/λ² = 5.31×10⁻²⁶. My own Monte-Carlo SHM rate code (`indep_rates.py`) gives N₉₀ = 1.865, a point limit of 4.23×10⁻⁴⁶ at 17.3 TeV, a composite limit of 5.98×10⁻⁴⁵ (×14.2), **b_LZ = 6.54×10⁻¹⁰ and b_XLZD(10 events) = 2.195×10⁻¹⁰** | holds [C] |
| O9 | detectable φ² portal: ΔS/Q = 3×10³–4×10⁴ vs 1.4–2.2 | Y_φ(b = 5×10⁻¹¹) = 2.230×10⁻⁸ analytically (K₁-weighted, 2 φ per decay), against the track's 2.237×10⁻⁸. The entropy integrand uses the energy rate Γ m_h n_h (K₂), which is right: each decay deposits E_h, not m_h. b_max ∝ ΔS^{−1/2} follows | holds [D/E] |
| O10 | trilinear ΔS/Q ≈ 0.4–4 at A = 2.7×10⁻¹⁰ GeV | `trilinear_freezein.py`: own reduced cross sections, Σ\|M\|²(q q̄ → gφ) = 16g²y²s²/(tu) and Σ\|M\|²(qg → qφ) = 16g²y²t²/(s(−u)), with a thermal-mass cutoff, top Boltzmann factor and T < 160 GeV. Fixed zero-T mixing gives **ΔS/Q = 0.98** (equivalent c = 0.25). Mean-field mixing θ(T) ∝ v(T)/m_h(T)² gives **2.5** (c = 0.64), peaking at T ≈ 150 GeV. Allowed XLZD events: 14–22 (fixed) or 6–9 (mean field) | holds as [E]; decided by EW-crossover physics (§3) |
| O11 | Xe mimics a 28 GeV WIMP; 31 events for 3σ | Poisson-likelihood Asimov (not Pearson): **27.0 GeV, 28 events** | holds [C] |
| O12 | kinetic mixing forbidden by hidden C | B·F_X is odd under hidden C alone and every portal operator is even. With α ≠ 0 the hidden side is C(P)-odd, so a CP-odd SM correlator ⟨J_Y O_even…⟩ is needed | holds [A] |
| O13 | M = 3×10⁻²⁰–1.3×10⁻⁵ g; σ/M(10²⁰) = 6×10⁻⁹ | bulk sphere R = 3.76×10⁵ m_X⁻¹, M = 1.28×10⁻⁵ g, σ/M = 6.3×10⁻⁹ cm²/g | holds [D] |

## 2. The vapour transfer: the single-fluid picture fails exactly where it produces the suppression

**Criterion [D].** The track's c_eff² = γ n_v T/ρ_DM is the tight-coupling (TC) limit of objects plus vapour. It needs the
vapour's momentum relaxation to beat its crossing of the wavelength: ε_TC = k v̄/(aΓ_m) ≪ 1, with v̄ = √(8T_h/π m_X).
- The vapour moves at v̄ = 0.07–0.12c.
- The loaded fluid moves at c_eff ≈ c_v/√(ρ_obj/ρ_v) ≈ 0.007–0.014c.
- So the vapour diffusion length near freeze-out exceeds the sound horizon. This is the R ≫ 1 limit of photon–baryon
  Silk damping, which is worse by √(1+R) ≈ 6.

Since Γ_m/H ∝ T, modes reach their acoustic phase k r_s(a) ≈ 1 late, when coupling is already weak. At acoustic onset
(`two_fluid_check.py`, `tc_diagnostic`):

| case (S/Q = 1.235) | ζ | ε_TC at onset, k = 10 / 20 / 50 /Mpc |
|---|---|---|
| Q = 2×10⁵, s = 1 | 1.7×10⁻² | 3.0 / 1.3 / 0.47 |
| Q = 10⁶, s = 1 | 1×10⁻² | 6.8 / 6.4 / 2.3 |
| Q = 10⁸, s = 1 | 2×10⁻³ | 1.8 / 3.5 / 8.8 |
| Q = 10¹⁰, s = 1 | 4.6×10⁻⁴ | 0.47 / 0.94 / 2.4 |
| Q ≥ 10¹² (s = 1); Q = 2×10⁵, s = 10⁻² | ≤ 10⁻⁴ | modes enter after a_m, or ε < 1 |

**Two-fluid re-integration [C/E].** I integrate pressureless objects plus a vapour fluid before a_m:
- momentum exchange Γ_m for the vapour and Γ_m ρ_v/ρ_o for the objects;
- chemical relaxation at Γ_c toward local two-phase equilibrium δ_v = rδ_o, with r = (3/2+x)/(3/2+x−x/(5/2+x));
- a vapour pressure normalised so the TC limit is exactly the track's c_eff²;
- a variant with a vapour viscosity ν = c_v²/(aΓ_m + kc_v).

After a_m the track's Gilbert marching is reused (copied), started from the separate object and gas states. The CDM
reference runs through the same pipeline.

**Control:** with every rate multiplied by 10⁶, the code reproduces the track's single-fluid P(k):
- 0.9614 vs 0.9649 at k = 50 and 1.407 vs 1.413 at k = 100 for Q = 2×10⁵;
- identical to 10⁻⁶ for Q = 10¹², 10¹⁴;
- half-step change ≤ 3×10⁻⁶.

Results (`skeptic/results/two_fluid_transfer.json`; track / two-fluid / two-fluid + viscosity):

| case | P(10) | P(20) | P(50) |
|---|---|---|---|
| Q = 2×10⁵, s = 1, S/Q 1.235 | 0.052 / **0.415** / 0.435 | 0.700 / **0.097** / 0.121 | 0.965 / **0.154** / 0.109 |
| same, S/Q 2.07 | 0.405 / **0.199** / 0.221 | 0.000 / 0.000 / 0.003 | 0.957 / **0.373** / 0.316 |
| Q = 10⁶, s = 1 | 0.420 / **0.701** / 0.719 | 0.019 / **0.442** / 0.486 | 1.088 / **0.050** / 0.090 |
| Q = 10⁸, s = 1 (ζ = 2×10⁻³), S/Q 1.235 | 0.879 / **0.901** / 0.900 | 0.709 / **0.819** / 0.819 | 0.059 / **0.521** / 0.561 |
| same, S/Q 2.07 | 0.805 / **0.839** / 0.837 | 0.557 / **0.713** / 0.713 | 0.018 / **0.312** / 0.362 |
| Q = 10¹⁰, s = 1 (ζ = 4.6×10⁻⁴) | 0.962 / 0.962 / 0.962 | 0.922 / 0.922 / 0.922 | 0.806 / **0.861** / 0.859 |
| Q = 10¹² (ζ = 10⁻⁴), S/Q 1.235 / 2.07 | 0.997 / 0.994 (same) | 0.984 / 0.973 (same) | **0.934 / 0.891 (same)** |
| Q = 10¹⁴ (ζ = 2×10⁻⁵) | 1.000 | 0.999 | 0.991 (same) |
| Q = 2×10⁵, s = 10⁻² (ζ = 1.7×10⁻⁴) | 0.995 / 0.995 | 0.978 / 0.979 | 0.859 / 0.870 |

**What changes.**
1. Every number the track reports at k ≳ 1/r_s for Q ≤ 10¹⁰ at s = 1 is a TC artefact. That includes the "oscillating"
   pattern and the headline "Q = 2×10⁵, s = 1: P(10) = 0.05". The physical modes leave the coupled fluid one by one as
   ε_TC crosses 1, not all at the global a_m.
2. **The regime conclusions survive.**
   - ζ ≲ 2×10⁻⁵ and ζ ≈ 10⁻⁴ are unchanged. Those modes enter after decoupling. The 7–11% deficit at k = 50 is the
     warm-fraction asymptote (1 − f_v)² = 0.933 / 0.889 for f_v = 3.4 / 5.8%, not acoustic physics.
   - Q = 2×10⁵ with s = 1 stays far below the WDM curves at every k ≥ 10.
   - The ζ ≈ 2×10⁻³ exclusion becomes **marginal**: P(10) = 0.90 / 0.84 against 0.94 (5.3 keV), and P(20) = 0.82 / 0.71
     against 0.74.
3. The two-fluid model is still a fluid. Landau (free-streaming) damping of the vapour is only mimicked by ν, and the
   ν-on/ν-off spread (≤ 10%) is a lower bound on its uncertainty. A BGK or Boltzmann treatment of the vapour, with
   velocity-dependent capture and isotropic re-emission in the object frame, is the calculation that would decide the
   k ≥ 10/Mpc values for ζ ≳ 5×10⁻⁴ [E].

## 3. New issues

1. **TC failure** (§2). For ζ ≳ 5×10⁻⁴ the "single fluid with c_eff ≈ 0.01c" description must be replaced by a
   two-component (or kinetic) one. The size of the k ≳ 10/Mpc suppression changes by up to ×10 in P.
2. **The burn-out has an entropy cost.** Each X̄ absorbed on an object releases 2Ω m_X = 1.76 m_X into the hidden bath,
   so Δ(S/Q) = 2Ωp m_X/T_h:
   - 0.41–0.68 for p = 0.0023 (REVIEW_KINETIC, R = 4.97) at T_h = 0.01–0.006;
   - 0.88–1.47 for p = 0.005;
   - 2.5–4.2 for p = 0.0143 (R = 7).

   In fx_formation's equilibrium branch this heat is already inside S/Q. But any X̄ remnant that reaches the object era
   (the case §3a costs "at face value") raises S/Q, and with it f_v ∝ S/Q and the P(k) deficit. It can push the medium
   branch (budget 1.4) toward the evaporation threshold. O5 and O4 are therefore not independent [D/E].
3. **The vapour fraction is not 2–6% in general.** f_v is 2–6% only for s ≈ 1. Small ζ from poor sticking (s ≤ 10⁻²)
   freezes out early at higher T_h, giving f_v = 7–35% (Q = 10¹⁰–10²⁰, s = 10⁻²–10⁻⁴, track's own table). Some of those
   rows freeze at T_h ≈ 0.046–0.05, the edge of the non-relativistic model. So the "CDM" end carries the largest warm
   fraction: P → (1 − f_v)² ≈ 0.42–0.86 at k ≳ 10³–10⁴/Mpc. That is unobservable today, but the digest's "2–6%" and
   "CDM" should be qualified.
4. **The 1% threshold depends on S/Q.** At ζ = 2.15×10⁻⁵ (Q = 10¹⁴, s = 1), P(50) = 0.991 for S/Q = 1.235 but 0.985 for
   2.07. "CDM within 1% at k ≤ 50" needs ζ ≲ 1×10⁻⁵ in the medium branch.
5. **The trilinear verdict rests on the electroweak crossover.** The φ yield is dominated by T ≈ 90–150 GeV, where θ(T) =
   A v(T)/m_h(T)² is not the zero-T value. Mean field raises Y by ×2.5 (ΔS/Q 1.0 → 2.5, at or above the 1.4–2.2
   budget), so "about 5–50 XLZD events" should read **≈ 6–22** [E]. A lattice-informed v(T), m_h(T) through the crossover
   decides whether the one viable extension survives at 10 events.
6. Minor:
   - ⁷⁴Ge rather than natural Ge in the zero table;
   - the φ-half radius is defined as the centre–vacuum midpoint;
   - the LZ toy's N₉₀ = 1.86 is a calibration, not a background-free 2.3. Its high-mass point limits (2.5×10⁻⁴⁷ at 1 TeV)
     were not checked against the published curve [R].

## Reproduction (scratch $SP/wt2/fx_observables_skeptic; OMP_NUM_THREADS=1)

```
S=research/reality_framework/grut_working_theory_v2_20261101/tracks/fx_observables/skeptic
python3 -I $S/code/indep_object_dd.py $SP/wt2/fx_observables_skeptic/indep_obj.json      # 60 s; copy json + _profile.npz to $S/results
python3 -I $S/code/indep_rates.py $S/results/indep_object_profile.npz $S/results/indep_rates.json   # 40 s
python3 -I $S/code/trilinear_freezein.py $S/results/trilinear_freezein.json                         # 1 s
python3 -I $S/code/two_fluid_check.py $S/results/two_fluid_transfer.json                            # 170 s
python3 -I $S/code/two_fluid_check.py $S/results/two_fluid_control.json                             # 60 s (forced-TC control)
```
