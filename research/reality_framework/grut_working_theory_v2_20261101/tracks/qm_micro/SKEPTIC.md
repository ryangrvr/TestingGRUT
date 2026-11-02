# SKEPTIC of track qm_micro (Q-matter branch microphysics, M1r)

2 November 2026 · Claude (adversarial skeptic, key qm_micro_skeptic), decisive pass (BRIEF_DECISIVE.md; issue #13
context). This is same-session AI work on one supplied model (M1r), not outside replication. Nothing here is a
validation, an exclusion of the real world or a new theory.

Every input is a model choice, not a consequence of closure: M1r, e, the kinetic closures, and the trapping and
vector-heat estimates. Labels: [D] derived, [C] computed with our own code, [E] estimate, [A] argued, [Q] quoted
from a reviewed track, [R] recalled or literature.

Ground rules followed:
- Git was read-only and nothing was committed.
- At most 2 processes ran, in the foreground. One chunk of about 1 minute ran past the tool's foreground limit; it
  finished before any further work.
- No lattice runs.
- Code is in `skeptic/code/`, results in `skeptic/results/` (656 kB in total).
- qm_micro's own modules were imported read-only (`python3 -B`) or copied before modification (`coarsen_heat.py`).

## Verdict in one paragraph

**The microphysics of item 1 stands.**
- I recomputed the planar sticking with a different method: my own collocation profile, then a collocation linear
  BVP with exact Robin radiation conditions instead of Numerov block elimination.
- It reproduces every quoted s(T, θ) and X̄ destruction value to ≤ 10⁻⁹, along with σ_planar = 99.037, the threshold
  law 7.89√T, c_s, the gaps 1.740/2.134, u = 0.170 and the Beliaev mean free paths.
- The thermal-flux averages agree to ≤ 0.1–0.25%.

**The kinetics of items 3–4 omit one physical process and contain one coding inconsistency.**
- **Latent-heat (non-isothermal) limitation.** Q-matter objects have no light neutral channel and negligible heat
  capacity, so every condensing X must export its heat Δ + ε through re-emitted quanta. In linear response the
  exchange flux is multiplied by F_h = 1 − W₁²/(W₀W₂). This is the Feder–Russell–Lothe–Pound factor; it equals
  0.004 at T = 0.01, 0.10 at 0.1, 0.28 at 0.3 and 0.8 at 1.
- **Coding inconsistency.** coarsen.py's exponential-μ exchange law drops the Bose factor B = T(∂Φ/∂μ)/Φ = 1.14–1.45
  at T = 0.3–8. The NOTES say "a few %".

**Corrected results.**
- **⟨Q⟩:** with both corrected it falls by ×0.54–0.70, to 3.6–3.8×10¹⁹, 7.7–8.0×10¹⁹, 1.17–1.27×10²⁰ and
  1.45–1.50×10²⁰ for S/Q = 1.5, 2.5, 3.0, 3.2.
- **f_v:** the post-freeze-out decrease of f_v in both of qm_micro's variants is forbidden by energy conservation on
  zero-heat-capacity objects. My energy-conserving two-moment model gives f_v = 0.056 / 0.097–0.098 / 0.120–0.121 /
  0.129–0.132 (perfect trapping), and 0.083 / 0.141 / 0.175 / 0.189 with Beliaev trapping. Each is above qm_micro's
  quoted upper end by 12–24%.
- **M/m_X and H_i:** M/m_X rises by 0.4–1.0%, and H_i shifts by −0.2 to −0.7%, to 1.48–1.52×10¹¹ GeV.
- **Unchanged:** ζ ≈ 0.9–1.4×10⁻⁷, so the verdict "CDM-like on Lyman-α and satellite scales" is unchanged.
- **α:** the "LSW selection" reading of R_c² = 0.372 K_W t is refuted. Five different initial edges give
  α = 0.371–0.373.

## Verdict table

| Result id | Verdict | Recomputed (this pass) |
|---|---|---|
| item1_s_normal | **confirmed** [C] | s = 0.4322177 (T = 0.005), 0.6799233, 0.8384031, 0.9274273, 0.9779456, 0.9914635, 0.9980790, 0.9993107 (T = 1); identical to scatter1d.json to ≤ 10⁻⁹; threshold s/√T = 7.887 (T = 10⁻⁵) |
| item1_s_oblique | **confirmed** [C] | 60°: 0.45068 / 0.78724 / 0.95523 / 0.99874; 85°: 0.10065 / 0.24738 / 0.44189 / 0.0971 (T = 0.5, partial TIR) / 0 (T = 1) |
| item1_thermal | **confirmed** [C] | independent quadrature on the verified table: s̄ = 0.8646 / 0.9494 / 0.9671 (Boltzmann, T_h = 0.1 / 0.3 / 0.5) vs 0.8643 / 0.9485 / 0.9659; X̄ destroyed 0.9088 / 0.9661 / 0.9796 vs 0.9083 / 0.9649 / 0.9770 |
| item1_lattice | not rerun; consistent | a third, independent method agrees with their frequency-domain values, which the lattice reproduced |
| item1_finite_object | not rerun | — |
| item1_trapping [E] | **inputs confirmed**; uncertainty larger than stated | u = 0.17024, c_s = 0.54669, e+P = 348.606; mfp 2.05×10⁷ (T = 0.005) … 2.32×10⁵ (0.3) reproduced; branch convex (decay allowed); see issues 6–7 |
| item2_GT | **K_W, n_sat confirmed; α and its interpretation corrected** | K_W(0.1) = 2.66×10⁻⁷ m_X by hand; α = 0.372 independent of initial edge (issue 3); fx_formation's 64/81 is 2.0–2.2× (not 2.6×) too fast |
| item3_Q | **corrected ×0.54–0.70** (not overturned) | 3.6–3.8×10¹⁹ / 7.7–8.0×10¹⁹ / 1.17–1.27×10²⁰ / 1.45–1.50×10²⁰ (Bose + heat; s_φ = 0…1) |
| item3_fv | **refuted as stated** | perfect trapping 0.0555–0.0560 / 0.0971–0.0983 / 0.1197–0.1213 / 0.1292–0.1315; Beliaev 0.083 / 0.141 / 0.175 / 0.189 |
| item3_freeze_vapour | numbers reproduced on their criterion; kinetic freeze-out earlier | Γ_a = H in the energy-conserving run at T_vis = 0.54–0.99 MeV (Beliaev 10.9–12.9 MeV), with T_v ≈ 0.5 T_o there |
| item3_zeta | **confirmed in substance** | ζ ≈ 0.9–1.4×10⁻⁷ with the corrected ⟨Q⟩ |
| item3_e [E] | **confirmed** (now [C] via dq_heat) | dq_heat's computed R(e): e = 5.5×10⁻⁶ (S/Q = 2.5), 2.98×10⁻⁶ (3.0), 2.55×10⁻⁶ (3.2), all inside qm_micro's band–broad brackets; S/Q = 1.5 unreachable (floor 1.84) |
| item4_M | **refuted mildly**: the range moves up | M/m_X = 0.7883–0.7944 / 0.7974–0.8076 / 0.8023–0.8152 / 0.8048–0.8182; H_i = 1.509–1.516 / 1.492–1.505 / 1.483–1.499 / 1.479–1.496 ×10¹¹ GeV |

## 1. Sticking recomputed by a different method [C] (`code/sk_scatter.py`, `results/sk_scatter.json`)

**Equations.** I derived the linearisation independently from the rescaled Lagrangian:
- L = (1/g²)[½|∂ψ|² + ½(∂S)² − U];
- ψ = e^{−iΩt}(F + u + iv), with (u, v, s) = (U, iV, S_x)e^{−iwt}.

The resulting real-symmetric matrix M(x; w, k⊥) is identical to scatter1d.py's.

**Numerics (all different from qm_micro):**
- **Bulk.** Nested brentq on P(Ω) = 0, giving Ω_b = 0.7758647074595398 and n_B = 449.3133679.
- **Interface.** scipy collocation on two half-lines glued at the pin F(0) = F_B/2, with Ω free and exact asymptotic
  Robin conditions at ±40. Results: Ω − Ω_b = 2×10⁻¹⁶ and σ = 99.03699 by both the equipartition and the direct
  formula.
- **Scattering.** A complex linear two-point BVP (collocation, adaptive mesh to 3.6×10⁴ nodes) with Robin
  conditions Y′ = ZY built from the exact asymptotic eigenmodes. The outgoing direction comes from the numerically
  differentiated dispersion, not from a norm formula.
- **Unitarity** is ≤ 3×10⁻¹².
- A first attempt by shooting failed. The vacuum tail F ∝ e^{−0.63x} is still 2×10⁻⁴ at x = 14, and the closed φ
  channel swamps the solutions. The BVP formulation avoids both problems.

**Comparison.** Every normal-incidence s and X̄-destruction value matches scatter1d.json to ≤ 10⁻⁹:
- X: 0.9274273121 vs 0.92742731 (T = 0.1); 0.9914635368 vs 0.99146354 (T = 0.3).
- X̄: 0.8130478832 vs 0.81304788 (T = 0.02).
- The oblique values match the item-1 table to all quoted digits.

The bulk dispersion gives:
- c_s(k → 0) = 0.54670 and gaps 1.74048 and 2.13359;
- w/k = 0.567, 0.588 and 0.636 at k = 0.4, 0.6 and 1.0. The phonon branch is convex, so Beliaev 1 → 2 decay is
  kinematically allowed.

**Thermal averages** (`code/sk_thermal_check.py`). I used Gauss–Laguerre in kinetic energy with PCHIP interpolation
on qm_micro's verified s(E, μ) table, plus a 4000-point μ re-integration. Agreement is ≤ 0.1% for X and ≤ 0.25% for
X̄ (table above). The 16-node μ rule is adequate even across the T ≳ 0.4 total-internal-reflection kink (≤ 0.07%).

## 2. NEW: the latent-heat (non-isothermal) limitation of every exchange flux [D/C] (`code/sk_heat.py`)

**Physics.**
- M1r has no light hidden neutral state (m_φ = √6; A, Y, Σ heavy).
- At T ≲ 1 the objects' heat capacity is that of the phonon gas: about 1.3×10⁻⁸ per charge at T = 0.009, against a
  latent heat Δ = 1 − Ω_b = 0.224 per condensed charge.
- So an object cannot stay at the vapour temperature while it grows or shrinks. Its temperature T_o adjusts until
  the net heat flux vanishes.
- The thermal relaxation time of an object is shorter than its ripening time by ~(λ² + 2)/(c_o F_h g) ≫ 1, with
  g = 2σ/(n_B R T) ≲ 10⁻⁴. The heat-balanced steady state therefore holds at every T, even where the liquid's heat
  capacity is O(1) per charge.

**Linear response** [D] (Onsager, detailed balance, Bose weights; per species i with charge q_i):
- weight w_i(E) = s_i f_i(1 + f_i) p²/(8π²);
- heat h_i = E − q_i Ω;
- moments W₀ = Σ q_i² ∫w_i, W₁ = Σ q_i ∫w_i h_i/T, W₂ = Σ ∫w_i h_i²/T².

The heat-balanced charge flux is J = W₀(1 − W₁²/(W₀W₂))Δa, against qm_micro's isothermal J = W₀Δa. So
**F_h = 1 − W₁²/(W₀W₂)**.

For a pure Maxwellian X vapour with s = 1 this is 2T²/((Δ + 2T)² + 2T²). That is exactly the Feder–Russell–Lothe–
Pound non-isothermal factor b²/(b² + q²) with no carrier gas (Adv. Phys. 15, 111 (1966)) [R]. My numerical moments
reproduce the closed form to ≤ 4% at T ≤ 0.1.

| T | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 1 | 3 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| F_h (s_φ = 1 / 0) | 0.0010 | 0.0037 | 0.0122 | 0.046 | 0.101 | 0.186 | 0.288 / 0.280 | 0.54 / 0.50 | 0.84 / 0.79 | 0.98 / 0.97 | 0.997 / 0.995 |

- Above T ≈ 0.5 the X̄, and φ if it is absorbed, carry heat without the X's charge–heat correlation, so F_h → 1.
- The φ absorption probability on Q-matter was not computed by qm_micro. I bracket it with s_φ = 0 and 1; this moves
  ⟨Q⟩ by ≤ 8%.
- **Consequence.** The ripening coefficient is K_W F_h, and every object's growth rate carries F_h. The global
  relaxation of the vapour is not slowed, because condensation heats the vapour; issue 4 covers this.

## 3. NEW: coarsen.py's exchange law misses the Bose kinetic coefficient [C] (`code/sk_bose.py`)

**What the code does.** coarsen.py writes the off-coexistence flux as Φ(μ_c)[e^{d/T} − e^{g}]. That is exponential
μ scaling, exact only for Boltzmann statistics. Its linear coefficient is Φ/T.

**What detailed balance gives.** The coefficient is ∂Φ/∂μ = ∫w s f(1+f)/T. This is the L of qm_micro's own K_W
table in item 2.

**The missing factor.** B_X = T(∂Φ/∂μ)/Φ, and the X̄ equivalent:

| T | 0.1 | 0.2 | 0.3 | 0.5 | 1 | 2 | 3–8 |
|---|---|---|---|---|---|---|---|
| B_X | 1.024 | 1.087 | 1.144 | 1.227 | 1.337 | 1.416 | 1.44–1.45 |
| B_X̄ | 1.000 | 1.000 | 1.001 | 1.006 | 1.035 | 1.09 | 1.13–1.22 |

This is not "a few %" (NOTES, Limits). The ripening that sets the object number happens at T ≈ 0.2–8, and there it
runs 14–45% too slowly in the code. So the K_W of item 2 and the K that coarsen.py integrates differ.

## 4. The self-similar ripening coefficient [C] (`code/sk_wagner.py`)

**Test.** I took qm_micro's own `wagner_test` kinetics (imported unchanged) with five initial distributions:
- lognormal cut at ±4σ (theirs);
- cut at ±2σ;
- extended to ±8σ;
- N ∝ (1 − R/R_max)^p with p = 2 and p = 6.

**Result.** α = R_c²/(K_W t) = 0.3720, 0.3706, 0.3725, 0.3725 and 0.3719 respectively (late-time spread ≤ 0.013), with
R_max/R_c = 1.466–1.470 throughout.

**Reading.**
- The coefficient does **not** depend on the initial edge. The NOTES' explanation ("below 0.5, as expected for
  compactly supported initial data, the LSW selection effect") is refuted.
- At T = 0.5 the coded/true kinetic coefficient is 0.8206, so relative to the coded kinetics α = 0.453. The rest
  (0.453 vs 0.5) is slow convergence or numerics of the re-resolution step. I did not resolve it.
- Against Wagner's asymptote ⟨R⟩² = (32/81)K_W t (Wagner 1961; Lifshitz–Slyozov 1961 [R]), fx_formation's
  (64/81)K_W is 2.0× too fast. Against the coded state it is 2.2×, not 2.6×.

## 5. ⟨Q⟩ with heat-balanced and Bose-correct exchange [C on the stated model] (`code/coarsen_heat.py`)

**Method.**
- coarsen_heat.py is a copy of qm_micro's coarsen.py. Its only physics change is that Φ_X and Φ_X̄ are multiplied by
  F_h(T), optionally also by B_X(T) and B_X̄(T).
- The gross capture rate in the freeze-out test is left unscaled.
- The histories were regenerated with qm_micro's qm_vapour.py into my scratch.
- The control run reproduces qm_micro to all printed digits: S/Q = 3.0 gives T_f = 0.00893, f_v = 0.1161 and
  ⟨Q⟩ = 1.899×10²⁰.

Charge-weighted ⟨Q⟩ at freeze-out:

| S/Q | qm_micro (control) | heat only (s_φ 1 / 0) | Bose only | **Bose + heat** (s_φ 1 / 0) | objects per charge (control → Bose+heat) |
|---|---|---|---|---|---|
| 1.5 | 5.50×10¹⁹ | — | — | **3.82 / 3.61×10¹⁹** (×0.70 / 0.66) | 3.1 → 4.3–4.6×10⁻²⁰ |
| 2.5 | 1.19×10²⁰ | 5.71 / 5.42×10¹⁹ (×0.48 / 0.45) | 1.66×10²⁰ (×1.39) | **7.96 / 7.67×10¹⁹** (×0.67 / 0.64) | 1.2 → 1.7–1.8×10⁻²⁰ |
| 3.0 | 1.90×10²⁰ | 8.29 / 8.26×10¹⁹ (×0.44 / 0.43) | 2.64×10²⁰ (×1.39) | **1.27 / 1.17×10²⁰** (×0.67 / 0.61) | 7.0 → 9.4–10.2×10⁻²¹ |
| 3.2 | 2.70×10²⁰ | — | — | **1.45 / 1.50×10²⁰** (×0.54 / 0.55) | 4.6 → 7.3–7.5×10⁻²¹ |

- Insensitivity to the first lumps survives: Q₀ = 1.2×10¹³ gives 1.17×10²⁰ (S/Q = 3.0) and 8.1×10¹⁹ (2.5).
- The two corrections partly cancel. The ripening-set ⟨Q⟩ ~ 10²⁰ stands at a factor ≲ 2.
- The f_v these runs print (0.15–0.21) is **not** physical. With T imposed, the slowed global condensation cannot be
  compensated by vapour heating. f_v is taken from §6 instead.

## 6. NEW (independent confirmation of qm_transfer §3.1): energy conservation fixes f_v after freeze-out [C] (`code/sk_freeze.py`)

**qm_micro's two post-freeze-out variants both continue net capture:**
- (a) on the equilibrium T(a);
- (b) capture only, with T_v ∝ a⁻².

On zero-heat-capacity objects with no neutral channel, each absorbed X must be balanced by re-emission at T_o. Net
capture stops when absorption does, and it does not continue below the coexistence line. Variant (b) violates energy
conservation outright; variant (a) imposes an equilibrium that no longer holds.

**My model.** A Maxwellian vapour (n_v, T_v) with exact flux moments of the sticking-weighted spectrum. The objects
are zero-heat-capacity bodies whose T_o solves Φ_a(Δ + ⟨ε⟩_a) = Φ_e(T_o)(Δ + ⟨ε⟩_e) at every instant, with fixed
object number and area ∝ η^{2/3} from the coarsening run at T = 0.05. It runs on the same histories and visible
clock as qm_micro.

**Validation.** The same code in qm_micro's mode reproduces qm_micro's variants (b)/(a) to ≤ 3%:
- 0.0392/0.0498 vs 0.0383/0.0496 (S/Q = 1.5);
- 0.0690/0.0868 vs 0.0668/0.0864 (2.5);
- 0.0847/0.1066 vs 0.0821/0.1061 (3.0);
- 0.0912/0.1154 vs 0.0885/0.1148 (3.2);
- Beliaev: 0.1233/0.1451 vs 0.1180/0.1429 (3.0).

| S/Q | qm_micro f_v (b … a) | **energy-conserving f_v** (qm_micro objects / Bose+heat objects) | ratio to qm_micro upper / lower | Beliaev: qm_micro → energy-conserving |
|---|---|---|---|---|
| 1.5 | 0.038–0.050 | **0.0560 / 0.0555** | ×1.13 / ×1.46 | 0.056–0.069 → **0.083** |
| 2.5 | 0.067–0.086 | **0.0983 / 0.0971** | ×1.14 / ×1.47 | 0.097–0.118 → **0.141** |
| 3.0 | 0.082–0.106 | **0.1213 / 0.1197** | ×1.14 / ×1.48 | 0.118–0.143 → **0.175** |
| 3.2 | 0.089–0.115 | **0.1315 / 0.1292** | ×1.15 / ×1.49 | 0.127–0.153 → **0.189** |

**What the model shows:**
- f_v tracks f_v,eq to ≤ 2% while Γ_a/H ≳ 20.
- T_v lags T_o: T_v ≈ 0.51–0.54 T_o at Γ_a = H. qm_transfer finds 0.57 at Γ/H ≈ 2.
- f_v is then constant to the end. The objects stay at T_o ≈ 0.7–1.0 T_eq while the vapour cools.
- Kinetic freeze-out (Γ_a = H) comes at T_vis = 0.54–0.99 MeV (perfect trapping) and 10.9–12.9 MeV (Beliaev), against
  qm_micro's 0.41–0.67 and 6.1–6.6 MeV.

**Comparison with qm_transfer.** This agrees in direction with qm_transfer's q-resolved kinetic background:
f_v = 0.064 / 0.113 / 0.151 for S/Q = 1.5 / 2.5 / 3.2 [Q]. Theirs is 13–15% above mine. They take the charge-weighted
⟨Q⟩ as a single size; the capture-effective Q is 1.2–1.4× smaller, so they have about 10% less area. The two methods
therefore bracket f_v at 0.098–0.113 (2.5) and 0.129–0.151 (3.2).

## 7. Mass per charge and H_i [C on E inputs]

The formula is qm_micro's: M = ⟨E/Q⟩ + Δ f_v with ⟨E/Q⟩ = 0.775866, p ≈ 0, and
H_i = 1.342842735×10¹¹ (0.9374/(0.991023 M))^{2/3} GeV.

| S/Q | qm_micro M (full range) | **this pass** (perfect … Beliaev) | qm_micro H_i | **this pass H_i** |
|---|---|---|---|---|
| 1.5 | 0.7845–0.7913 | **0.7883–0.7944** | 1.513–1.521 | **1.509–1.516×10¹¹** |
| 2.5 | 0.7908–0.8022 | **0.7974–0.8076** | 1.499–1.513 | **1.492–1.505×10¹¹** |
| 3.0 | 0.7943–0.8079 | **0.8023–0.8152** | 1.492–1.509 | **1.483–1.499×10¹¹** |
| 3.2 | 0.7957–0.8102 | **0.8048–0.8182** | 1.489–1.507 | **1.479–1.496×10¹¹** |

The spread is now set by interior trapping alone. The "post-freeze-out cooling" half of qm_micro's spread was the
energy-violating variant. Y/Y_Q = 0.991023 is dq_entropy's updated value; fx_holder's 0.990733 would move H_i by
2×10⁻⁴.

## 8. e, ζ and the observable

- **e(S/Q).** dq_heat's computed R(e) table gives e = 5.5×10⁻⁶, 2.98×10⁻⁶ and 2.55×10⁻⁶ for S/Q = 2.5, 3.0 and 3.2.
  All lie inside qm_micro's reading-B band–broad brackets (4.9–5.8, 2.7–3.2 and 2.27–2.69×10⁻⁶). S/Q = 1.5 stays
  unreachable (e → ∞ floor 1.84, reproduced from the closed form) [Q/C].
- **ζ = s_eff⟨Q⟩^{−1/3}** with the corrected ⟨Q⟩ is 1.4, 1.1, 0.95 and 0.9×10⁻⁷. The Beliaev variant keeps ~6×10⁻⁹.
- **Free streaming.** My hand estimate (v = 0.16c at T_vis = 0.5 MeV, v ∝ 1/a, radiation era, ln(a_eq/a_f) ≈ 13.7)
  gives λ_fs ≈ 0.34 kpc, against qm_micro's 0.34–0.61. Earlier freeze-out and qm_transfer's ×3 re-emission heating of
  σ_p change this by O(1) only.
- **qm_transfer's kinetic Boltzmann code** puts the 20-keV reach at ζ ~ 10⁻⁴ (f_v ~ 3–10%). **The CDM-like verdict
  is robust by 3 decades** [Q/E].

## 9. New issues (beyond qm_micro's own caveats)

1. **Latent-heat-limited exchange omitted** (§2). Every growth and ripening flux carries F_h, from 0.004 at T = 0.01
   to 0.28 at 0.3 and 0.8 at 1. ⟨Q⟩ ×0.43–0.48 from this alone.
2. **Bose factor dropped by the exponential-μ exchange law** (§3). B_X = 1.14–1.45 at T = 0.3–8; ⟨Q⟩ ×1.39 alone.
   Combined with issue 1, ⟨Q⟩ ×0.54–0.70. The K of item 2 and the K actually integrated were inconsistent.
3. **α = 0.372 is not an LSW selection effect** (§4). It is independent of the initial edge, and 18% of the gap to
   0.5 is issue 2. fx_formation's error factor is 2.0–2.2×, not 2.6×.
4. **Post-freeze-out net capture violates energy conservation** (§6). Variants (a) and (b) both understate f_v. The
   energy-conserving f_v is 12–24% above qm_micro's upper ends and 42–49% above variant (b). Kinetic freeze-out comes
   1.3–2× earlier in T_vis.
5. **The base case contradicts the track's own trapping estimate.** qm_micro quotes P_trap ≈ 0.05–0.1 at T_f for
   Q ~ 10²⁰ but takes perfect trapping as the base case. If Beliaev decay is the only trapping channel, the
   Beliaev-trapping branch is the central one: f_v = 0.08–0.19, M up to 0.818, H_i down to 1.479×10¹¹ GeV. Its ⟨Q⟩
   (0.2–1.6×10²¹) was not recomputed with issues 1–2.
6. **Beliaev inputs [E], uncertain by ×2–3:**
   - the phonon energy is taken as Δ + 1.5T, while the flux-weighted absorbed energy is Δ + ~2.5T (mfp ×0.85 at T_f,
     ×0.37 at T = 0.1);
   - the relativistic three-phonon vertex is not fixed; a (1 + u − c_s²) vertex would give mfp ×1.8;
   - k = 0.41–2.3 lies outside the k⁵ regime (w/k exceeds c_s by 4–39%).
7. **φ absorption on Q-matter was not computed.** It is the only carrier gas at T ≳ 0.3. Bracketed s_φ = 0…1, it moves
   ⟨Q⟩ by ≤ 8%. φ → X conversion at the interface is also absent.
8. **σ(T).** The T = 0 surface tension is used at T_h = 1–8, where ~30–40% of the ripening integral accrues. The
   thermal ⟨|ψ|²⟩/F_B² ≈ 0.04 at T = 8 suggests ≲ 10% [E].
9. **Cross-track** (qm_transfer, read 2 Nov): its single charge-weighted ⟨Q⟩ gives about 10% too little capture area,
   so its f_v is high by ~5%. Its ⟨Q⟩ input (qm_micro's) should be divided by ~1.5–1.9 (issues 1–2), which raises ζ
   by ~1.2. Neither changes its verdict.

## 10. Limits of this skeptic pass

- No lattice or finite-object reruns. I trust the frequency-domain agreement plus my third method.
- The two-moment vapour is a Maxwellian closure. qm_transfer's q-resolved solution is the finer one, and the two
  differ by ~15% in f_v.
- The F_h and B corrections enter coarsen.py as rate factors on the T-imposed history. That is exact for the
  differential (ripening) flux in linear response but not for the global condensation near freeze-out, which is why
  f_v comes from §6.
- The Beliaev-variant ⟨Q⟩ was not rerun.
- The early dense stage, the Fuchs–Epstein kernel and dq_heat's R(e) are taken as given.

## Reproduction (OMP_NUM_THREADS=1; SCR = $SP/wt2/qm_micro_skeptic; about 15 min at ≤ 2 processes)

```
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/qm_micro/skeptic/code
python3 -B sk_scatter.py ../results/sk_scatter.json                      # 4 min
python3 -B sk_thermal_check.py ../results/sk_thermal_check.json
python3 -B sk_beliaev.py ../results/sk_beliaev.json
python3 -B sk_heat.py ../results/sk_heat.json
python3 -B -W ignore sk_bose.py ../results/sk_bose.json
(cd ../../code; for sq in 1.5 2.5 3.0 3.2; do python3 -B -W ignore qm_vapour.py $SCR/hist_$sq.json $sq; done)
python3 -B -W ignore coarsen_heat.py $SCR/cheat_3.0.json 3.0          # (and 2.5) control + heat cases, 20 s each
SK_CASES=bose_only,bose_heat_sphi1,bose_heat_sphi0,bose_heat_sphi1_Q0_1.2e13 python3 -B -W ignore coarsen_heat.py $SCR/cheat2_3.0.json 3.0   # (and 2.5)
SK_CASES=control_computed,bose_heat_sphi1,bose_heat_sphi0 python3 -B -W ignore coarsen_heat.py $SCR/cheat_1.5.json 1.5                       # (and 3.2)
python3 -B -W ignore sk_wagner.py ../results/sk_wagner.json trunc4,trunc2,exp,smooth2,smooth6
python3 -B -W ignore sk_freeze.py ../results/sk_freeze.json           # 2 min
python3 -B sk_summary.py                                               # sk_summary.json, sk_coarsen.json
```

Credits [R]:
- Feder, Russell, Lothe & Pound, Adv. Phys. 15, 111 (1966) (non-isothermal growth factor).
- Lifshitz & Slyozov 1961; Wagner 1961; Niethammer & Pego, J. Stat. Phys. 95, 867 (1999).
- Beliaev 1958 and Landau–Khalatnikov (phonon damping).

No numbers were taken from the literature.
